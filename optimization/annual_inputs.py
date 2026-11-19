"""Build a Phase 3 scenario from version-controlled annual inputs."""

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from optimization.phase3 import Feedstock, Limits, Plant, Scenario

DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "scenarios" / "annual.json"


@dataclass(frozen=True)
class LoadedScenario:
    scenario: Scenario
    label: str
    assumptions: tuple[str, ...]
    input_sha256: str


def overlay(original: dict[str, Any], changes: dict[str, Any]) -> dict[str, Any]:
    """Merge variant fields and reject misspelled overrides."""
    result = original.copy()
    for key, value in changes.items():
        if key not in original:
            raise ValueError(f"unknown override: {key}")
        if isinstance(original[key], dict) and isinstance(value, dict) and key != "day_overrides":
            value = overlay(original[key], value)
        result[key] = value
    return result


def load_scenario(path: Path, scenario_id: str, days: int | None = None) -> LoadedScenario:
    config = json.loads(path.read_text(encoding="utf-8"))
    horizon = config["horizon_days"]
    if type(horizon) is not int or horizon < 1:
        raise ValueError("configured horizon must be a positive integer")
    if days is None:
        days = horizon
    if type(days) is not int or not 1 <= days <= horizon:
        raise ValueError("requested horizon must be within the configured horizon")

    params = config["parameters"]
    assumptions = tuple(config["assumptions"])
    label = config["label"]
    if scenario_id != config["scenario_id"]:
        if scenario_id not in config["alternatives"]:
            raise ValueError(f"unknown scenario: {scenario_id}")
        variant = config["alternatives"][scenario_id]
        params = overlay(params, variant["overrides"])
        label = variant["label"]
        assumptions += (variant["note"],)
    if days < horizon:
        assumptions += ("Short-horizon check applies the terminal policy at the requested cutoff.",)

    overrides = params["day_overrides"]
    for key, change in overrides.items():
        if not key.isascii() or not key.isdecimal() or str(int(key)) != key or not 1 <= int(key) <= horizon:
            raise ValueError(f"invalid day override: {key}")
        overlay(params["daily"], change)

    tails = []
    for day in range(1, days + 1):
        tails.append(overlay(params["daily"], overrides.get(str(day), {})))
    feeds = {name: Feedstock(*values) for name, values in params["feeds"].items()}
    factors = params["emissions_per_useful_gj"]
    plants = {
        name: Plant(values, {fuel: factors["biomass"] for fuel in feeds})
        for name, values in params["efficiencies"].items()
    }
    scenario = Scenario(
        scenario_id, feeds, plants, tuple(Limits(**day) for day in tails),
        params["initial_m3"], params["terminal_m3"], params["gas_cost_per_useful_gj"], factors["gas"],
    )
    identity = json.dumps(
        {"id": scenario_id, "horizon": days, "inputs": params, "assumptions": assumptions},
        sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode()
    return LoadedScenario(scenario, label, assumptions, hashlib.sha256(identity).hexdigest())
