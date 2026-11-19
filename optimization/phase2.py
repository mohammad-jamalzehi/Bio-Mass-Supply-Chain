"""Small cost-minimization MILP for one biomass fuel and one plant."""

from dataclasses import dataclass
from enum import IntEnum
from math import isclose
from typing import Sequence

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

from validation.phase1 import (
    InventoryDay,
    daily_emissions_tco2e,
    nonnegative,
    useful_heat_gj,
    validate_inventory,
)

TOLERANCE = 1e-7


@dataclass(frozen=True)
class Fuel:
    truck_volume_m3: float
    wet_density_t_per_m3: float
    moisture_fraction: float
    heating_value_gj_per_dry_t: float
    efficiency: float
    cost_per_dry_t: float
    emissions_tco2e_per_useful_gj: float

    def heat(self, volume_m3: float) -> float:
        return useful_heat_gj(
            volume_m3, self.wet_density_t_per_m3, self.moisture_fraction,
            self.heating_value_gj_per_dry_t, self.efficiency,
        )


@dataclass(frozen=True)
class DayLimits:
    demand_gj: float
    truck_limit: int
    storage_capacity_m3: float
    plant_capacity_gj: float
    gas_capacity_gj: float | None = None


@dataclass(frozen=True)
class Scenario:
    days: tuple[DayLimits, ...]
    fuel: Fuel
    initial_inventory_m3: float
    terminal_inventory_m3: float | None
    gas_cost_per_useful_gj: float
    gas_emissions_tco2e_per_useful_gj: float

    def __post_init__(self) -> None:
        if not self.days:
            raise ValueError("at least one day is required")
        self.fuel.heat(1)
        fields = {
            "truck volume": self.fuel.truck_volume_m3,
            "biomass price": self.fuel.cost_per_dry_t,
            "biomass emissions factor": self.fuel.emissions_tco2e_per_useful_gj,
            "initial inventory": self.initial_inventory_m3,
            "gas price": self.gas_cost_per_useful_gj,
            "gas emissions factor": self.gas_emissions_tco2e_per_useful_gj,
        }
        if self.terminal_inventory_m3 is not None:
            fields["terminal inventory"] = self.terminal_inventory_m3
        for name, value in fields.items():
            nonnegative(name, value)
        if self.fuel.truck_volume_m3 == 0:
            raise ValueError("truck volume must be positive")
        for day in self.days:
            if type(day.truck_limit) is not int or day.truck_limit < 0:
                raise ValueError("truck limit must be a nonnegative integer")
            nonnegative("heat demand", day.demand_gj)
            nonnegative("storage capacity", day.storage_capacity_m3)
            nonnegative("plant capacity", day.plant_capacity_gj)
            if day.gas_capacity_gj is not None:
                nonnegative("gas capacity", day.gas_capacity_gj)


@dataclass(frozen=True)
class DayPlan:
    trucks: int
    inventory: InventoryDay
    gas_heat_gj: float


@dataclass(frozen=True)
class Totals:
    cost: float
    emissions_tco2e: float


@dataclass(frozen=True)
class SolveResult:
    status: int
    message: str
    days: tuple[DayPlan, ...] = ()
    totals: Totals | None = None
    best_bound: float | None = None
    mip_gap: float | None = None


def verify_plan(
    scenario: Scenario, days: Sequence[DayPlan], *, expected_cost: float | None = None,
) -> Totals:
    """Audit an operating plan using physical quantities, outside solver rows."""
    if len(days) != len(scenario.days):
        raise ValueError("plan must contain one entry per day")
    validate_inventory(
        (day.inventory for day in days), initial=scenario.initial_inventory_m3,
        terminal=scenario.terminal_inventory_m3, tolerance=TOLERANCE,
    )
    fuel = scenario.fuel
    cost = emissions = 0.0
    for index, (day, limits) in enumerate(zip(days, scenario.days, strict=True), start=1):
        stock = day.inventory
        nonnegative("gas heat", day.gas_heat_gj)
        if type(day.trucks) is not int or not 0 <= day.trucks <= limits.truck_limit:
            raise ValueError(f"day {index}: truckload count violates integrality or limit")
        if not isclose(stock.delivered, day.trucks * fuel.truck_volume_m3, rel_tol=0, abs_tol=TOLERANCE):
            raise ValueError(f"day {index}: deliveries do not match full truckloads")
        biomass_heat = fuel.heat(stock.fresh_burn) + fuel.heat(stock.stored_burn)
        if biomass_heat > limits.plant_capacity_gj + TOLERANCE:
            raise ValueError(f"day {index}: biomass heat exceeds plant capacity")
        if limits.gas_capacity_gj is not None and day.gas_heat_gj > limits.gas_capacity_gj + TOLERANCE:
            raise ValueError(f"day {index}: gas heat exceeds gas capacity")
        if not isclose(biomass_heat + day.gas_heat_gj, limits.demand_gj, rel_tol=0, abs_tol=TOLERANCE):
            raise ValueError(f"day {index}: useful heat does not match demand")
        if stock.closing > limits.storage_capacity_m3 + TOLERANCE:
            raise ValueError(f"day {index}: closing inventory exceeds storage capacity")
        dry_tonnes = stock.delivered * fuel.wet_density_t_per_m3 * (1 - fuel.moisture_fraction)
        cost += dry_tonnes * fuel.cost_per_dry_t + day.gas_heat_gj * scenario.gas_cost_per_useful_gj
        emissions += daily_emissions_tco2e(
            [(biomass_heat, fuel.emissions_tco2e_per_useful_gj)],
            gas_heat_gj=day.gas_heat_gj,
            gas_factor_tco2e_per_gj=scenario.gas_emissions_tco2e_per_useful_gj,
        )
    if expected_cost is not None and not isclose(cost, expected_cost, rel_tol=0, abs_tol=TOLERANCE):
        raise ValueError("recomputed cost disagrees with solver objective")
    return Totals(cost, emissions)


class Column(IntEnum):
    TRUCKS = 0
    FRESH = 1
    STORED = 2
    CLOSING = 3
    GAS = 4


def solve(scenario: Scenario) -> SolveResult:
    """Build and solve the toy MILP. Return a checked plan only on optimal status."""
    shape = (len(scenario.days), len(Column))
    cost = np.zeros(shape)
    integrality = np.zeros(shape, dtype=int)
    upper_bounds = np.full(shape, np.inf)
    rows: list[np.ndarray] = []
    lower: list[float] = []
    upper: list[float] = []
    fuel = scenario.fuel
    truck_volume = fuel.truck_volume_m3
    heat = fuel.heat(1)
    dry_tonnes_per_truck = truck_volume * fuel.wet_density_t_per_m3 * (1 - fuel.moisture_fraction)
    cost[:, Column.TRUCKS] = dry_tonnes_per_truck * fuel.cost_per_dry_t
    cost[:, Column.GAS] = scenario.gas_cost_per_useful_gj
    integrality[:, Column.TRUCKS] = 1

    def add_row(
        day: int, terms: dict[Column, float], lo: float, hi: float,
        *, previous_stock: float = 0,
    ) -> None:
        row = np.zeros(shape)
        for column, coefficient in terms.items():
            row[day, column] = coefficient
        if day > 0:
            row[day - 1, Column.CLOSING] = previous_stock
        rows.append(row.ravel())
        lower.append(lo)
        upper.append(hi)

    for day, limits in enumerate(scenario.days):
        upper_bounds[day, Column.TRUCKS] = limits.truck_limit
        upper_bounds[day, Column.CLOSING] = limits.storage_capacity_m3
        if limits.gas_capacity_gj is not None:
            upper_bounds[day, Column.GAS] = limits.gas_capacity_gj
        initial = scenario.initial_inventory_m3 if day == 0 else 0

        # Closing + fresh + stored - deliveries - prior closing = initial on day 1, else 0.
        add_row(day, {Column.CLOSING: 1, Column.FRESH: 1, Column.STORED: 1,
                      Column.TRUCKS: -truck_volume}, initial, initial, previous_stock=-1)
        add_row(day, {Column.FRESH: 1, Column.TRUCKS: -truck_volume}, -np.inf, 0)
        add_row(day, {Column.STORED: 1}, -np.inf, initial, previous_stock=-1)
        add_row(day, {Column.FRESH: heat, Column.STORED: heat, Column.GAS: 1},
                limits.demand_gj, limits.demand_gj)
        add_row(day, {Column.FRESH: heat, Column.STORED: heat}, -np.inf, limits.plant_capacity_gj)

    terminal = scenario.terminal_inventory_m3
    if terminal is not None:
        add_row(len(scenario.days) - 1, {Column.CLOSING: 1}, terminal, terminal)

    result = milp(
        cost.ravel(), integrality=integrality.ravel(),
        bounds=Bounds(np.zeros(cost.size), upper_bounds.ravel()),
        constraints=LinearConstraint(np.asarray(rows), lower, upper),
        options={"mip_rel_gap": 0.0},
    )
    best_bound = getattr(result, "mip_dual_bound", None)
    gap = getattr(result, "mip_gap", None)
    if result.status != 0:
        return SolveResult(result.status, result.message, best_bound=best_bound, mip_gap=gap)

    days = []
    opening = scenario.initial_inventory_m3
    for trucks, fresh, stored, closing, gas in result.x.reshape(shape):
        rounded_trucks = round(float(trucks))
        if not isclose(float(trucks), rounded_trucks, rel_tol=0, abs_tol=TOLERANCE):
            raise ValueError("solver returned a fractional truckload count")
        stock = InventoryDay(opening, rounded_trucks * truck_volume, fresh, stored, closing)
        days.append(DayPlan(rounded_trucks, stock, gas))
        opening = closing

    totals = verify_plan(scenario, days, expected_cost=result.fun)
    if best_bound is None or not isclose(totals.cost, best_bound, rel_tol=0, abs_tol=TOLERANCE):
        raise ValueError("optimal status was returned without a matching objective bound")
    return SolveResult(result.status, result.message, tuple(days), totals, best_bound, gap)
