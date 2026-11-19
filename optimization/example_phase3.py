"""Print a solver-audited synthetic result as dashboard-ready JSON."""

import argparse
import json
from dataclasses import replace

from optimization.phase3 import Feedstock, Limits, Plant, Scenario, solve, synthetic_result


def synthetic_scenario(case: str = "base") -> Scenario:
    feeds = {
        "wc": Feedstock(2, 0.5, 0.2, 10, 3),
        "wp": Feedstock(2.5, 0.5, 0.2, 10, 5),
        "bq": Feedstock(2, 0.5, 0.2, 10, 5),
    }
    factors = {key: 0.1 for key in feeds}
    plants = {
        "existing": Plant({"wc": 0.5, "wp": 0.5, "bq": 0.75}, factors),
        "new": Plant({"wc": 0.75, "wp": 1, "bq": 0.75}, factors),
    }
    first = Limits(6, 1, 1, {"existing": 2, "new": 4}, pellet_silo_m3=0.5)
    second = Limits(2, 0, 1, {"existing": 2, "new": 2}, pellet_silo_m3=0.5)
    if case == "tight-silo":
        first = replace(first, pellet_silo_m3=0.4)
    elif case == "no-storage":
        first = replace(first, total_storage_m3=0)
    elif case != "base":
        raise ValueError(f"unknown synthetic case: {case}")
    zeros = {key: 0.0 for key in feeds}
    return Scenario(case, feeds, plants, (first, second), zeros, zeros.copy(), 3, 0.5)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("base", "tight-silo", "no-storage"), default="base")
    scenario = synthetic_scenario(parser.parse_args().case)
    result = solve(scenario)
    print(json.dumps(synthetic_result(scenario, result), indent=2, allow_nan=False))
    if result.status != 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
