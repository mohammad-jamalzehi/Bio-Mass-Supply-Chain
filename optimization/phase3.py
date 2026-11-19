"""Multi-feedstock, multi-plant MILP and independent dispatch audit."""

from dataclasses import asdict, dataclass
from math import inf, isclose
from typing import Mapping, Sequence

import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_array

from validation.phase1 import (
    InventoryDay, daily_emissions_tco2e, nonnegative, useful_heat_gj, validate_inventory,
)

TOL = 1e-7


@dataclass(frozen=True)
class Feedstock:
    truck_m3: float
    wet_density: float
    moisture: float
    heating_value: float
    cost_per_dry_tonne: float

    def dry_tonnes(self, volume_m3: float) -> float:
        return volume_m3 * self.wet_density * (1 - self.moisture)

    def heat(self, volume_m3: float, efficiency: float) -> float:
        return useful_heat_gj(volume_m3, self.wet_density, self.moisture, self.heating_value, efficiency)


@dataclass(frozen=True)
class Plant:
    efficiency: Mapping[str, float]
    emissions_per_useful_gj: Mapping[str, float]


@dataclass(frozen=True)
class Limits:
    demand_gj: float
    truck_limit: int
    total_storage_m3: float
    plant_heat_capacity_gj: Mapping[str, float]
    pellet_silo_m3: float | None = None
    gas_capacity_gj: float | None = None


@dataclass(frozen=True)
class Scenario:
    name: str
    feedstocks: Mapping[str, Feedstock]
    plants: Mapping[str, Plant]
    days: tuple[Limits, ...]
    initial_m3: Mapping[str, float]
    terminal_m3: Mapping[str, float] | None
    gas_cost_per_useful_gj: float
    gas_emissions_per_useful_gj: float

    def __post_init__(self) -> None:
        fuels, plants = set(self.feedstocks), set(self.plants)
        if not fuels or not plants or not self.days:
            raise ValueError("feedstocks, plants and days are required")
        for fuel in self.feedstocks.values():
            fuel.heat(1, 1)
            nonnegative("truck volume", fuel.truck_m3)
            nonnegative("biomass price", fuel.cost_per_dry_tonne)
            if fuel.truck_m3 == 0:
                raise ValueError("truck volume must be positive")
        for plant in self.plants.values():
            if set(plant.efficiency) != fuels or set(plant.emissions_per_useful_gj) != fuels:
                raise ValueError("plant efficiencies and emissions factors must cover all feedstocks")
            for key, fuel in self.feedstocks.items():
                fuel.heat(1, plant.efficiency[key])
                nonnegative("biomass emissions factor", plant.emissions_per_useful_gj[key])
        for stocks in (self.initial_m3, self.terminal_m3):
            if stocks is not None:
                if set(stocks) != fuels:
                    raise ValueError("inventory targets must cover all feedstocks")
                for value in stocks.values():
                    nonnegative("inventory target", value)
        nonnegative("gas price", self.gas_cost_per_useful_gj)
        nonnegative("gas emissions factor", self.gas_emissions_per_useful_gj)
        for day in self.days:
            if type(day.truck_limit) is not int or day.truck_limit < 0:
                raise ValueError("truck limit must be a nonnegative integer")
            if set(day.plant_heat_capacity_gj) != plants:
                raise ValueError("daily limits must cover all plants")
            for value in (day.demand_gj, day.total_storage_m3, *day.plant_heat_capacity_gj.values()):
                nonnegative("daily limit", value)
            for optional in (day.pellet_silo_m3, day.gas_capacity_gj):
                if optional is not None:
                    nonnegative("optional daily limit", optional)
            if day.pellet_silo_m3 is not None and "wp" not in fuels:
                raise ValueError("a pellet silo requires the wood-pellet feedstock wp")


@dataclass(frozen=True)
class Burn:
    fresh_m3: float
    stored_m3: float


@dataclass(frozen=True)
class DayPlan:
    trucks: Mapping[str, int]
    closing_m3: Mapping[str, float]
    burns: Mapping[str, Mapping[str, Burn]]  # plant -> feedstock -> burn
    gas_gj: float


@dataclass(frozen=True)
class DayAudit:
    inventory_by_fuel: dict[str, InventoryDay]
    biomass_gj_by_plant_and_fuel: dict[str, dict[str, float]]
    cost: float
    emissions_tco2e: float


@dataclass(frozen=True)
class SolveResult:
    status: int
    message: str
    best_bound: float | None
    mip_gap: float | None
    plans: tuple[DayPlan, ...] = ()
    audit: tuple[DayAudit, ...] = ()


def verify_plan(scenario: Scenario, plans: Sequence[DayPlan]) -> tuple[DayAudit, ...]:
    """Recalculate physical and monetary totals independently of the solver matrix."""
    fuels, plants = scenario.feedstocks, scenario.plants
    opening = dict(scenario.initial_m3)
    if len(plans) != len(scenario.days):
        raise ValueError("expected one plan per day")
    audited = []
    for number, (plan, limits) in enumerate(zip(plans, scenario.days, strict=True), start=1):
        if set(plan.trucks) != set(fuels) or set(plan.closing_m3) != set(fuels):
            raise ValueError(f"day {number}: incomplete fuel plan")
        if set(plan.burns) != set(plants) or any(set(burns) != set(fuels) for burns in plan.burns.values()):
            raise ValueError(f"day {number}: incomplete plant allocation")
        if any(type(n) is not int or n < 0 for n in plan.trucks.values()) or sum(plan.trucks.values()) > limits.truck_limit:
            raise ValueError(f"day {number}: truck counts violate integrality or shared limit")

        inventory = {}
        cost = 0.0
        for key, fuel in fuels.items():
            fresh = [plan.burns[plant][key].fresh_m3 for plant in plants]
            stored = [plan.burns[plant][key].stored_m3 for plant in plants]
            for volume in (*fresh, *stored):
                nonnegative(f"day {number} burn", volume)
            delivered = plan.trucks[key] * fuel.truck_m3
            inventory[key] = InventoryDay(opening[key], delivered, sum(fresh), sum(stored), plan.closing_m3[key])
            cost += fuel.dry_tonnes(delivered) * fuel.cost_per_dry_tonne
        if sum(plan.closing_m3.values()) > limits.total_storage_m3 + TOL:
            raise ValueError(f"day {number}: total storage capacity exceeded")
        if limits.pellet_silo_m3 is not None and plan.closing_m3["wp"] > limits.pellet_silo_m3 + TOL:
            raise ValueError(f"day {number}: pellet silo capacity exceeded")

        heat = {}
        for plant_key, plant in plants.items():
            heat[plant_key] = {}
            for key, fuel in fuels.items():
                burn = plan.burns[plant_key][key]
                heat[plant_key][key] = fuel.heat(burn.fresh_m3 + burn.stored_m3, plant.efficiency[key])
            if sum(heat[plant_key].values()) > limits.plant_heat_capacity_gj[plant_key] + TOL:
                raise ValueError(f"day {number}: {plant_key} heat capacity exceeded")
        nonnegative(f"day {number} gas heat", plan.gas_gj)
        if limits.gas_capacity_gj is not None and plan.gas_gj > limits.gas_capacity_gj + TOL:
            raise ValueError(f"day {number}: gas heat capacity exceeded")
        if not isclose(sum(sum(values.values()) for values in heat.values()) + plan.gas_gj,
                       limits.demand_gj, rel_tol=0, abs_tol=TOL):
            raise ValueError(f"day {number}: useful heat does not meet demand")
        cost += plan.gas_gj * scenario.gas_cost_per_useful_gj
        emissions = daily_emissions_tco2e(
            ((heat[p][f], plants[p].emissions_per_useful_gj[f]) for p in plants for f in fuels),
            gas_heat_gj=plan.gas_gj, gas_factor_tco2e_per_gj=scenario.gas_emissions_per_useful_gj,
        )
        audited.append(DayAudit(inventory, heat, cost, emissions))
        opening = dict(plan.closing_m3)

    for key in fuels:
        terminal = None if scenario.terminal_m3 is None else scenario.terminal_m3[key]
        try:
            validate_inventory((day.inventory_by_fuel[key] for day in audited),
                               initial=scenario.initial_m3[key], terminal=terminal, tolerance=TOL)
        except ValueError as error:
            raise ValueError(f"{key}: {error}") from error
    return tuple(audited)


def solve(scenario: Scenario) -> SolveResult:
    fuels, plants = scenario.feedstocks, scenario.plants
    columns = {}
    costs, integrality, variable_upper = [], [], []
    row_indices, column_indices, values, lower, upper = [], [], [], [], []

    def variable(key: tuple, cost: float = 0, integer: bool = False, maximum: float = inf) -> None:
        columns[key] = len(costs)
        costs.append(cost)
        integrality.append(int(integer))
        variable_upper.append(maximum)

    def row(terms: dict[tuple, float], lo: float = -inf, hi: float = inf) -> None:
        for key, coefficient in terms.items():
            if coefficient != 0:
                row_indices.append(len(lower))
                column_indices.append(columns[key])
                values.append(coefficient)
        lower.append(lo)
        upper.append(hi)

    for t, day in enumerate(scenario.days):
        variable(("gas", t), scenario.gas_cost_per_useful_gj,
                 maximum=inf if day.gas_capacity_gj is None else day.gas_capacity_gj)
        for f, fuel in fuels.items():
            variable(("trucks", t, f), fuel.dry_tonnes(fuel.truck_m3) * fuel.cost_per_dry_tonne, integer=True)
            variable(("stock", t, f))
            for p in plants:
                variable(("fresh", t, f, p))
                variable(("stored", t, f, p))

    for t, day in enumerate(scenario.days):
        row({("trucks", t, f): 1 for f in fuels}, hi=day.truck_limit)
        row({("stock", t, f): 1 for f in fuels}, hi=day.total_storage_m3)
        if day.pellet_silo_m3 is not None:
            row({("stock", t, "wp"): 1}, hi=day.pellet_silo_m3)
        for f, fuel in fuels.items():
            fresh = {("fresh", t, f, p): 1 for p in plants}
            stored = {("stored", t, f, p): 1 for p in plants}
            delivery = {("trucks", t, f): -fuel.truck_m3}
            prior = {("stock", t - 1, f): -1} if t else {}
            initial = 0 if t else scenario.initial_m3[f]
            row({**fresh, **delivery}, hi=0)
            row({**stored, **prior}, hi=initial)
            row({("stock", t, f): 1, **fresh, **stored, **delivery, **prior}, lo=initial, hi=initial)

        demand = {("gas", t): 1.0}
        for p, plant in plants.items():
            heat = {(stream, t, f, p): fuel.heat(1, plant.efficiency[f])
                    for f, fuel in fuels.items() for stream in ("fresh", "stored")}
            row(heat, hi=day.plant_heat_capacity_gj[p])
            demand.update(heat)
        row(demand, lo=day.demand_gj, hi=day.demand_gj)

    if scenario.terminal_m3 is not None:
        for f, value in scenario.terminal_m3.items():
            row({("stock", len(scenario.days) - 1, f): 1}, lo=value, hi=value)

    matrix = coo_array((values, (row_indices, column_indices)), shape=(len(lower), len(costs))).tocsc()
    answer = milp(
        costs, integrality=integrality, bounds=Bounds(np.zeros(len(costs)), variable_upper),
        constraints=LinearConstraint(matrix, lower, upper), options={"mip_rel_gap": 0.0},
    )
    bound, gap = getattr(answer, "mip_dual_bound", None), getattr(answer, "mip_gap", None)
    if answer.status != 0:
        return SolveResult(answer.status, answer.message, bound, gap)

    def value(*key: str | int) -> float:
        return float(answer.x[columns[key]])

    plans = []
    for t in range(len(scenario.days)):
        trucks = {}
        for f in fuels:
            count = value("trucks", t, f)
            trucks[f] = round(count)
            if not isclose(count, trucks[f], rel_tol=0, abs_tol=TOL):
                raise ValueError("solver returned a fractional truck count")
        plans.append(DayPlan(
            trucks, {f: value("stock", t, f) for f in fuels},
            {p: {f: Burn(value("fresh", t, f, p), value("stored", t, f, p)) for f in fuels} for p in plants},
            value("gas", t),
        ))
    audit = verify_plan(scenario, plans)
    cost = sum(day.cost for day in audit)
    if (not isclose(cost, answer.fun, rel_tol=0, abs_tol=TOL)
            or bound is None or not isclose(cost, bound, rel_tol=0, abs_tol=TOL)
            or gap is None or not isclose(gap, 0, rel_tol=0, abs_tol=TOL)):
        raise ValueError("solver objective, bound or gap disagrees with the independently checked cost")
    return SolveResult(answer.status, answer.message, bound, gap, tuple(plans), audit)


def synthetic_result(scenario: Scenario, output: SolveResult) -> dict:
    """Package the synthetic example for JSON export without writing a file."""
    return {
        "schema_version": 1,
        "data_kind": "synthetic_test_fixture",
        "scenario": scenario.name,
        "feedstocks": list(scenario.feedstocks),
        "plants": list(scenario.plants),
        "units": {"volume": "m3", "heat": "GJ useful", "cost": "synthetic USD", "emissions": "tCO2e"},
        "solver": {"backend": "HiGHS via SciPy milp", "scipy_version": scipy.__version__,
                   "status": output.status, "message": output.message,
                   "best_cost_bound": output.best_bound, "mip_gap": output.mip_gap},
        "totals": {"cost": sum(day.cost for day in output.audit),
                   "emissions_tco2e": sum(day.emissions_tco2e for day in output.audit)} if output.audit else None,
        "days": [
            {"day": index, "limits": asdict(limits), **asdict(plan), **asdict(audit)}
            for index, (limits, plan, audit) in enumerate(
                zip(scenario.days[:len(output.plans)], output.plans, output.audit, strict=True), start=1
            )
        ],
    }
