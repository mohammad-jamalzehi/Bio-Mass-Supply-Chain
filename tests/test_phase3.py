import json
import unittest
from dataclasses import replace

from optimization.example_phase3 import synthetic_scenario
from optimization.phase3 import Burn, DayPlan, solve, synthetic_result, verify_plan


class MultiPlantTests(unittest.TestCase):
    def test_three_hand_enumerated_competitions(self):
        # With no day-2 trucks and at most one on day 1, these are all choices.
        # Chips: <=6 GJ over two days, or 4 + (2 - 4/3)*2 GJ without storage.
        # Briquettes: 6 GJ in either case. Pellets need 0.5 m3 closing on day 1.
        candidates = {
            "base": {"gas": 24, "wc": 8.4, "wp": 5, "bq": 10},
            "tight-silo": {"gas": 24, "wc": 8.4, "bq": 10},
            "no-storage": {"gas": 24, "wc": 10.4, "bq": 10},
        }
        for case, costs in candidates.items():
            with self.subTest(case=case):
                scenario = synthetic_scenario(case)
                result = solve(scenario)
                winner = min(costs, key=costs.__getitem__)
                self.assertEqual(result.status, 0)
                self.assertAlmostEqual(result.best_bound, costs[winner])
                self.assertAlmostEqual(result.mip_gap, 0)
                self.assertAlmostEqual(sum(day.cost for day in result.audit), costs[winner])
                self.assertEqual(result.plans[0].trucks, {f: int(f == winner) for f in scenario.feedstocks})
                self.assertEqual(sum(result.plans[1].trucks.values()), 0)

    def test_two_binding_plants_and_pellet_inventory_continuity(self):
        result = solve(synthetic_scenario())
        first, second = result.audit
        self.assertAlmostEqual(first.biomass_gj_by_plant_and_fuel["existing"]["wp"], 2)
        self.assertAlmostEqual(first.biomass_gj_by_plant_and_fuel["new"]["wp"], 4)
        self.assertAlmostEqual(first.inventory_by_fuel["wp"].closing, 0.5)
        self.assertAlmostEqual(second.inventory_by_fuel["wp"].opening, 0.5)
        self.assertAlmostEqual(second.inventory_by_fuel["wp"].stored_burn, 0.5)
        self.assertAlmostEqual(second.inventory_by_fuel["wp"].closing, 0)
        self.assertAlmostEqual(second.biomass_gj_by_plant_and_fuel["new"]["wp"], 2)
        self.assertAlmostEqual(sum(day.emissions_tco2e for day in result.audit), 0.8)

    def test_fresh_and_stored_can_supply_either_plant(self):
        scenario = synthetic_scenario()
        first = replace(scenario.days[0], truck_limit=0, gas_capacity_gj=0)
        initial = {"wc": 0, "wp": 2, "bq": 0}
        scenario = replace(scenario, days=(first,), initial_m3=initial)
        result = solve(scenario)
        self.assertEqual(result.status, 0)
        for plant in scenario.plants:
            self.assertAlmostEqual(result.plans[0].burns[plant]["wp"].stored_m3, 1)
            self.assertAlmostEqual(result.plans[0].burns[plant]["wp"].fresh_m3, 0)

    def test_two_fuels_delivered_and_burned_together(self):
        scenario = synthetic_scenario()
        first = replace(scenario.days[0], truck_limit=2, total_storage_m3=2, gas_capacity_gj=0)
        # Each required closing stock forces its own truck. Both trucks also
        # leave 1 m3 to burn. Chips must supply 3 GJ in the new unit; the
        # briquettes split 2 GJ existing and 1 GJ new to meet demand of 6.
        scenario = replace(
            scenario, days=(first,), terminal_m3={"wc": 1, "wp": 0, "bq": 1},
        )
        result = solve(scenario)
        self.assertEqual(result.status, 0)
        self.assertEqual(result.plans[0].trucks, {"wc": 1, "wp": 0, "bq": 1})
        self.assertAlmostEqual(result.best_bound, 6.4)
        self.assertAlmostEqual(result.mip_gap, 0)
        self.assertAlmostEqual(result.audit[0].biomass_gj_by_plant_and_fuel["new"]["wc"], 3)
        self.assertAlmostEqual(result.audit[0].biomass_gj_by_plant_and_fuel["existing"]["bq"], 2)
        self.assertAlmostEqual(result.audit[0].biomass_gj_by_plant_and_fuel["new"]["bq"], 1)
        self.assertEqual(solve(replace(scenario, days=(replace(first, truck_limit=1),))).status, 2)
        self.assertEqual(solve(replace(scenario, days=(replace(first, total_storage_m3=1.9),))).status, 2)

    def test_mixed_deliveries_cannot_exceed_shared_truck_limit(self):
        scenario = synthetic_scenario()
        first = replace(scenario.days[0], total_storage_m3=5, pellet_silo_m3=3)
        plan = DayPlan({"wc": 1, "wp": 1, "bq": 0}, {"wc": 2, "wp": 2.5, "bq": 0},
                       {p: {f: Burn(0, 0) for f in scenario.feedstocks} for p in scenario.plants}, 6)
        with self.assertRaisesRegex(ValueError, "shared limit"):
            verify_plan(replace(scenario, days=(first,), terminal_m3=None), [plan])

    def test_total_storage_is_a_shared_limit_and_silo_is_pellet_only(self):
        scenario = synthetic_scenario()
        first = replace(scenario.days[0], truck_limit=2, total_storage_m3=4, pellet_silo_m3=3)
        plan = DayPlan({"wc": 1, "wp": 1, "bq": 0}, {"wc": 2, "wp": 2.5, "bq": 0},
                       {p: {f: Burn(0, 0) for f in scenario.feedstocks} for p in scenario.plants}, 6)
        new = replace(scenario, days=(first,), terminal_m3=None)
        with self.assertRaisesRegex(ValueError, "total storage"):
            verify_plan(new, [plan])
        with self.assertRaisesRegex(ValueError, "pellet silo"):
            verify_plan(replace(new, days=(replace(first, total_storage_m3=5, pellet_silo_m3=1),)), [plan])
        nonpellets = replace(plan, trucks={"wc": 1, "wp": 0, "bq": 1}, closing_m3={"wc": 2, "wp": 0, "bq": 2})
        verify_plan(replace(new, days=(replace(first, pellet_silo_m3=0),)), [nonpellets])

    def test_cross_plant_fuel_availability_checked_by_phase_one(self):
        scenario = synthetic_scenario()
        result = solve(scenario)
        two_plants = result.plans[0]
        # Keep the same plant heat and inventory arithmetic but misclassify a
        # fresh draw as stored even though the day started empty.
        burns = {plant: dict(fuels) for plant, fuels in two_plants.burns.items()}
        burns["new"]["wp"] = Burn(0.5, 0.5)
        mutated = replace(two_plants, burns=burns)
        with self.assertRaisesRegex(ValueError, "stored burn exceeds opening"):
            verify_plan(scenario, (mutated, result.plans[1]))

        second = result.plans[1]
        burns = {plant: dict(fuels) for plant, fuels in second.burns.items()}
        burns["new"]["wp"] = Burn(0.5, 0)
        with self.assertRaisesRegex(ValueError, "fresh burn exceeds deliveries"):
            verify_plan(scenario, (two_plants, replace(second, burns=burns)))

    def test_independent_audit_rejects_plant_capacity_violation(self):
        scenario = synthetic_scenario()
        result = solve(scenario)
        first = replace(scenario.days[0], plant_heat_capacity_gj={"existing": 1, "new": 4})
        with self.assertRaisesRegex(ValueError, "existing heat capacity"):
            verify_plan(replace(scenario, days=(first, scenario.days[1])), result.plans)

    def test_infeasible_with_no_pellet_storage_and_no_gas(self):
        scenario = synthetic_scenario("no-storage")
        scenario = replace(scenario, days=tuple(replace(day, gas_capacity_gj=0) for day in scenario.days))
        result = solve(scenario)
        self.assertEqual(result.status, 2)
        exported = synthetic_result(scenario, result)
        self.assertIsNone(exported["totals"])
        self.assertEqual(exported["days"], [])
        json.dumps(exported, allow_nan=False)

    def test_mixed_fuels_and_gas_emissions_added_once(self):
        scenario = synthetic_scenario()
        first = replace(scenario.days[0], truck_limit=2, total_storage_m3=2)
        scenario = replace(scenario, days=(first,), terminal_m3=None)
        burns = {p: {f: Burn(0, 0) for f in scenario.feedstocks} for p in scenario.plants}
        burns["existing"]["wc"] = Burn(1, 0)
        burns["new"]["bq"] = Burn(1, 0)
        plan = DayPlan({"wc": 1, "wp": 0, "bq": 1}, {"wc": 1, "wp": 0, "bq": 1}, burns, 1)
        audit, = verify_plan(scenario, [plan])
        self.assertAlmostEqual(audit.cost, 9.4)
        self.assertAlmostEqual(audit.emissions_tco2e, 1.0)

    def test_json_contract_has_daily_audited_values(self):
        scenario = synthetic_scenario()
        data = json.loads(json.dumps(synthetic_result(scenario, solve(scenario)), allow_nan=False))
        self.assertEqual(data["data_kind"], "synthetic_test_fixture")
        self.assertEqual(data["schema_version"], 1)
        self.assertEqual(data["solver"]["status"], 0)
        self.assertAlmostEqual(data["totals"]["cost"], 5)
        self.assertAlmostEqual(data["days"][0]["inventory_by_fuel"]["wp"]["delivered"], 2.5)
        self.assertAlmostEqual(data["days"][0]["burns"]["existing"]["wp"]["fresh_m3"], 1)
