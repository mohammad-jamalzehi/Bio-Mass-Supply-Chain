"""Run audited historical-input experiments and export dashboard data."""

import argparse
import csv
import json
import sys
from dataclasses import asdict, replace
from math import fsum, isclose, isfinite
from pathlib import Path
from time import perf_counter
from unittest.mock import patch

import scipy

from optimization import phase3
from optimization.annual_inputs import DEFAULT_CONFIG, LoadedScenario, load_scenario

COST_SCALE = 1000
COST_TOL = 0.0001
FLOW_TOL = 1e-6


def solver_version() -> dict:
    # The public milp result gives no HiGHS version. Read the pinned wheel's
    # bundled constants if exposed, without depending on them to solve.
    try:
        from scipy.optimize._highspy import _core
        version = ".".join(str(getattr(_core, f"HIGHS_VERSION_{part}")) for part in ("MAJOR", "MINOR", "PATCH"))
    except (ImportError, AttributeError):
        version = None
    return {"backend": "HiGHS via SciPy milp", "scipy_version": scipy.__version__, "highs_version": version}


def annual_audit(scenario: phase3.Scenario, plans: tuple[phase3.DayPlan, ...], objective: float) -> tuple[tuple, dict]:
    """Recheck all daily constraints and reconcile whole-horizon totals."""
    audited = phase3.verify_plan(scenario, plans)
    residuals = {}
    purchase_costs = []
    for f, fuel in scenario.feedstocks.items():
        ledger = [day.inventory_by_fuel[f] for day in audited]
        delivered = fsum(item.delivered for item in ledger)
        burned = fsum(item.fresh_burn + item.stored_burn for item in ledger)
        residuals[f] = fsum((scenario.initial_m3[f], delivered, -burned, -ledger[-1].closing))
        if abs(residuals[f]) > FLOW_TOL:
            raise ValueError(f"{f}: annual material balance residual exceeds tolerance")
        purchase_costs.append(delivered * fuel.wet_density * (1 - fuel.moisture) * fuel.cost_per_dry_tonne)
    gas_gj = fsum(plan.gas_gj for plan in plans)
    biomass_gj = fsum(
        heat for day in audited for plant in day.biomass_gj_by_plant_and_fuel.values() for heat in plant.values()
    )
    demand = fsum(day.demand_gj for day in scenario.days)
    if not isclose(biomass_gj + gas_gj, demand, rel_tol=0, abs_tol=FLOW_TOL):
        raise ValueError("annual useful heat does not equal demand")
    cost = fsum(purchase_costs) + gas_gj * scenario.gas_cost_per_useful_gj
    if not all(isclose(cost, value, rel_tol=0, abs_tol=COST_TOL)
               for value in (objective, fsum(day.cost for day in audited))):
        raise ValueError("annual recomputed cost differs from solver objective")
    return audited, {
        "cost": cost, "demand_gj": demand, "gas_gj": gas_gj, "biomass_gj": biomass_gj,
        "emissions_tco2e": fsum(day.emissions_tco2e for day in audited),
        "fuel_balance_residual_m3": residuals,
    }


def finite(value, scale: float = 1) -> float | None:
    if value is None or not isfinite(value):
        return None
    return float(value) * scale


def decode_incumbent(scenario: phase3.Scenario, vector) -> tuple[phase3.DayPlan, ...]:
    """Read the unchanged Phase 3 column order when it returns no plan."""
    width = 1 + len(scenario.feedstocks) * (2 + 2 * len(scenario.plants))
    plans = []
    for values in vector.reshape(len(scenario.days), width):
        entries = iter(float(value) for value in values)
        gas = next(entries)
        trucks, closing = {}, {}
        burns = {plant: {} for plant in scenario.plants}
        for fuel in scenario.feedstocks:
            count = next(entries)
            trucks[fuel] = round(count)
            if not isclose(count, trucks[fuel], rel_tol=0, abs_tol=phase3.TOL):
                raise ValueError("solver incumbent has a fractional truckload")
            closing[fuel] = next(entries)
            for plant in scenario.plants:
                burns[plant][fuel] = phase3.Burn(next(entries), next(entries))
        plans.append(phase3.DayPlan(trucks, closing, burns, gas))
    return tuple(plans)


def run(loaded: LoadedScenario, time_limit: float) -> dict:
    """Run one experiment synchronously. The CLI uses a dedicated process."""
    if not isfinite(time_limit) or time_limit <= 0:
        raise ValueError("solver time limit must be finite and positive")
    original_milp = phase3.milp
    raw = None
    seconds = None

    def bounded_milp(*args, **kwargs):
        nonlocal raw, seconds
        options = {**kwargs.pop("options", {}), "time_limit": time_limit}
        started = perf_counter()
        raw = original_milp(*args, options=options, **kwargs)
        seconds = perf_counter() - started
        if raw.x is not None:
            # HiGHS can return tiny negative values at zero. Keep material
            # violations visible and canonicalize only within Phase 3 tolerance.
            raw.x = raw.x.copy()
            raw.x[(raw.x < 0) & (raw.x >= -phase3.TOL)] = 0.0
        return raw

    model = loaded.scenario
    scaled = replace(model, gas_cost_per_useful_gj=model.gas_cost_per_useful_gj / COST_SCALE,
                     feedstocks={key: replace(fuel, cost_per_dry_tonne=fuel.cost_per_dry_tonne / COST_SCALE)
                                 for key, fuel in model.feedstocks.items()})
    # Phase 3 remains byte-for-byte unchanged. In the isolated experiment
    # process, add a solver deadline and retain raw termination diagnostics.
    validation_error = None
    with patch.object(phase3, "milp", bounded_milp):
        try:
            result = phase3.solve(scaled)
        except ValueError as error:
            if raw is None or raw.status != 0:
                raise
            validation_error = str(error)
            result = None

    objective = finite(getattr(raw, "fun", None), COST_SCALE)
    bound = finite(getattr(raw, "mip_dual_bound", None), COST_SCALE)
    diagnostics = {
        **solver_version(), "status_code": int(raw.status), "message": raw.message,
        "runtime_seconds": seconds, "time_limit_seconds": time_limit,
        "solver_reported_objective": objective, "best_bound": bound, "mip_gap": finite(getattr(raw, "mip_gap", None)),
    }
    record = {
        "schema_version": 1, "data_kind": "historical_input_experiment", "campus_inputs_validated": False,
        "emissions_accuracy": "unverified_factors_and_heat_basis",
        "scenario_id": model.name, "label": loaded.label, "days_requested": len(model.days),
        "input_sha256": loaded.input_sha256, "assumptions": list(loaded.assumptions),
        "units": {"cost": "source dollars", "heat": "GJ useful", "volume": "m3", "emissions": "tCO2e provisional"},
        "inputs": asdict(model), "solver": diagnostics, "validation": "no_verified_plan", "totals": None, "days": [],
    }
    if raw.status not in (0, 1) or raw.x is None or objective is None:
        return record
    try:
        if validation_error is not None:
            raise ValueError(validation_error)
        decoded = decode_incumbent(model, raw.x)
        if raw.status == 0 and decoded != result.plans:
            raise ValueError("annual decoder disagrees with the Phase 3 column order")
        audited, totals = annual_audit(model, decoded, objective)
        if raw.status == 0 and (bound is None or not isclose(totals["cost"], bound, rel_tol=0, abs_tol=COST_TOL)):
            raise ValueError("annual objective does not match the solver bound")
        if bound is not None and bound > totals["cost"] + COST_TOL:
            raise ValueError("solver bound exceeds independently checked incumbent cost")
    except ValueError as error:
        record.update(validation="audit_failed", validation_error=str(error))
        return record
    record.update(
        validation="verified_against_stated_scenario" if raw.status == 0 else "audited_feasible_incumbent",
        totals=totals,
        days=[{"day": number, "demand_gj": limits.demand_gj, **asdict(plan), **asdict(audit)}
              for number, (plan, audit, limits) in enumerate(zip(decoded, audited, model.days, strict=True), start=1)],
    )
    return record


def export(record: dict, directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "result.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    with (directory / "flows.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("day", "plant", "feedstock", "trucks", "opening_m3", "delivered_m3", "closing_m3",
                         "fresh_burn_m3", "stored_burn_m3", "useful_heat_gj", "day_demand_gj", "day_cost", "day_emissions_tco2e"))
        for day in record["days"]:
            for fuel, inventory in day["inventory_by_fuel"].items():
                writer.writerow((day["day"], "", fuel, day["trucks"][fuel], inventory["opening"],
                                 inventory["delivered"], inventory["closing"], "", "", "", "", "", ""))
            for plant, fuels in day["burns"].items():
                for fuel, burn in fuels.items():
                    writer.writerow((day["day"], plant, fuel, "", "", "", "", burn["fresh_m3"],
                                     burn["stored_m3"], day["biomass_gj_by_plant_and_fuel"][plant][fuel], "", "", ""))
            writer.writerow((day["day"], "gas", "", "", "", "", "", "", "", day["gas_gj"],
                             day["demand_gj"], day["cost"], day["emissions_tco2e"]))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--scenario", default="gams-historical-inputs")
    parser.add_argument("--days", type=int)
    parser.add_argument("--time-limit", type=float, default=30)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = run(load_scenario(args.config, args.scenario, args.days), args.time_limit)
    if args.output:
        export(record, args.output)
    summary_keys = ("scenario_id", "days_requested", "input_sha256", "validation", "solver", "totals", "validation_error")
    print(json.dumps({key: record[key] for key in summary_keys if key in record}, allow_nan=False))
    return 0 if record["validation"] == "verified_against_stated_scenario" else 2


if __name__ == "__main__":
    sys.exit(main())
