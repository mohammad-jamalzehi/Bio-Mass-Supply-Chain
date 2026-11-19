"""Regenerate the website fixture snapshot from the validated Python solvers."""
import json
import sys
from dataclasses import asdict
from pathlib import Path

WEB = Path(__file__).resolve().parents[1]
ROOT = WEB.parent
sys.path.insert(0, str(ROOT))

from optimization import example, example_phase3, phase2, phase3

def main():
    scenario = example.synthetic_scenario()
    result = phase2.solve(scenario)
    assert result.status == 0 and result.totals
    fixtures = [{
        "phase": 2, "scenario": "one-fuel", "data_kind": "synthetic_test_fixture",
        "solver": {"status": result.status, "best_cost_bound": result.best_bound, "mip_gap": result.mip_gap},
        "totals": asdict(result.totals), "days": [asdict(day) for day in result.days],
    }]
    for case in ("base", "tight-silo", "no-storage"):
        model = example_phase3.synthetic_scenario(case)
        solved = phase3.solve(model)
        assert solved.status == 0
        fixtures.append({"phase": 3, **phase3.synthetic_result(model, solved)})
    annual = []
    for path in sorted((ROOT / "results/phase4").glob("*/result.json")):
        row = json.loads(path.read_text(encoding="utf-8"))
        annual.append({k: row[k] for k in ("scenario_id", "data_kind", "campus_inputs_validated", "days_requested",
                     "emissions_accuracy", "solver", "validation", "input_sha256")})
    (WEB / "src/content/results.json").write_text(json.dumps({"fixtures": fixtures, "annual": annual},
        indent=2, allow_nan=False) + "\n", encoding="utf-8")

if __name__ == "__main__":
    main()
