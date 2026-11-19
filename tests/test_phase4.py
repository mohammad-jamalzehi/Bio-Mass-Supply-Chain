import csv
import io
import json
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import mock_open, patch

from optimization.annual import annual_audit, export, run
from optimization.annual_inputs import DEFAULT_CONFIG, load_scenario, overlay

BASE = "gams-historical-inputs"


class AnnualTests(unittest.TestCase):
    def test_full_horizon_and_alternative_inputs(self):
        base = load_scenario(DEFAULT_CONFIG, BASE)
        self.assertEqual(len(base.scenario.days), 365)
        self.assertEqual(base.scenario.days[364].demand_gj, 12624)
        capacity = load_scenario(DEFAULT_CONFIG, "report-narrative-capacity")
        self.assertEqual(capacity.scenario.days[0].plant_heat_capacity_gj["new"], 1036.8)
        report = load_scenario(DEFAULT_CONFIG, "report-table-emissions")
        self.assertEqual(report.scenario.gas_emissions_per_useful_gj, 0.056)
        self.assertEqual(report.scenario.days, base.scenario.days)
        self.assertNotEqual(report.input_sha256, base.input_sha256)
        self.assertEqual(base.input_sha256, load_scenario(DEFAULT_CONFIG, BASE, 365).input_sha256)
        self.assertNotEqual(base.input_sha256, load_scenario(DEFAULT_CONFIG, BASE, 7).input_sha256)

    def test_config_overrides_apply_on_named_day_only(self):
        config = json.loads(DEFAULT_CONFIG.read_text(encoding="utf-8"))
        config["parameters"]["day_overrides"] = {"2": {"demand_gj": 50}}
        with patch.object(Path, "read_text", return_value=json.dumps(config)):
            model = load_scenario(DEFAULT_CONFIG, BASE, 3).scenario
        self.assertEqual([day.demand_gj for day in model.days], [12624, 50, 12624])
        with self.assertRaisesRegex(ValueError, "unknown override"):
            overlay({"demand_gj": 50}, {"demand": 3})

    def test_bad_horizon_and_scenario_are_rejected(self):
        for days in (0, 366):
            with self.assertRaises(ValueError):
                load_scenario(DEFAULT_CONFIG, BASE, days)
        with self.assertRaisesRegex(ValueError, "unknown scenario"):
            load_scenario(DEFAULT_CONFIG, "missing")

    def test_short_horizon_is_audited_and_exportable(self):
        loaded = load_scenario(DEFAULT_CONFIG, BASE, 3)
        record = run(loaded, time_limit=10)
        self.assertEqual(record["validation"], "verified_against_stated_scenario")
        self.assertEqual(record["solver"]["status_code"], 0)
        self.assertAlmostEqual(record["solver"]["mip_gap"], 0)
        self.assertAlmostEqual(record["solver"]["best_bound"], record["totals"]["cost"])
        self.assertEqual(record["totals"]["demand_gj"], 12624 * 3)
        for residual in record["totals"]["fuel_balance_residual_m3"].values():
            self.assertAlmostEqual(residual, 0)
        self.assertEqual(len(record["days"]), 3)
        json.dumps(record, allow_nan=False)
        with patch.object(Path, "mkdir"), patch.object(Path, "write_text"), patch.object(Path, "open", return_value=mock_open().return_value) as opened:
            export(record, Path("unused"))
        csv_text = "".join(call.args[0] for call in opened.return_value.write.call_args_list)
        rows = list(csv.DictReader(io.StringIO(csv_text)))
        self.assertEqual(len(rows), 30)
        self.assertEqual(sum(row["plant"] == "gas" for row in rows), 3)

    def test_infeasible_result_has_diagnostics_without_a_plan(self):
        loaded = load_scenario(DEFAULT_CONFIG, BASE, 1)
        scenario = loaded.scenario
        loaded = replace(loaded, scenario=replace(scenario, days=(replace(scenario.days[0], gas_capacity_gj=0),)))
        result = run(loaded, time_limit=10)
        self.assertEqual(result["solver"]["status_code"], 2)
        self.assertEqual(result["days"], [])
        self.assertIsNone(result["totals"])
        self.assertEqual(result["validation"], "no_verified_plan")

    def test_time_limit_with_incumbent_is_audited_without_claiming_optimality(self):
        from optimization import phase3

        original = phase3.milp

        def limited(*args, **kwargs):
            answer = original(*args, **kwargs)
            answer.status = 1
            answer.message = "Synthetic termination fixture: time limit"
            answer.mip_dual_bound = answer.fun - 0.05
            answer.mip_gap = 0.05 / answer.fun
            return answer

        with patch.object(phase3, "milp", limited):
            record = run(load_scenario(DEFAULT_CONFIG, BASE, 1), time_limit=10)
        self.assertEqual(record["validation"], "audited_feasible_incumbent")
        self.assertEqual(record["solver"]["status_code"], 1)
        self.assertGreater(record["totals"]["cost"], record["solver"]["best_bound"])
        self.assertEqual(len(record["days"]), 1)

    def test_audit_rejects_incomplete_plan(self):
        model = load_scenario(DEFAULT_CONFIG, BASE, 7).scenario
        with self.assertRaisesRegex(ValueError, "one plan per day"):
            annual_audit(model, (), 0)
