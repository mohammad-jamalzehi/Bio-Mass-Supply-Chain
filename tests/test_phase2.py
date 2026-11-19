import unittest
from dataclasses import replace
from itertools import product

from optimization.example import synthetic_scenario
from optimization.phase2 import DayLimits, solve, verify_plan
from validation.phase1 import validate_inventory


class SolverValidationTests(unittest.TestCase):
    def setUp(self):
        self.scenario = synthetic_scenario()

    def test_optimum_matches_hand_enumeration_and_solver_bound(self):
        # Separate oracle for this fixture: every day-3 truck delivers 9 GJ
        # against demand of only 3 GJ and zero required closing stock.
        costs = {}
        for first, third in product((0, 1), repeat=2):
            if third == 0:
                costs[first, 0, third] = 6 * first + 3 * (18 - 9 * first)
        self.assertEqual(costs, {(0, 0, 0): 54, (1, 0, 0): 33})
        best_schedule = min(costs, key=costs.__getitem__)

        result = solve(self.scenario)
        self.assertEqual(result.status, 0)
        self.assertEqual(tuple(day.trucks for day in result.days), best_schedule)
        self.assertAlmostEqual(result.totals.cost, costs[best_schedule])
        self.assertAlmostEqual(result.best_bound, costs[best_schedule])
        self.assertAlmostEqual(result.mip_gap, 0)
        self.assertAlmostEqual(result.totals.emissions_tco2e, 5.4)
        self.assertEqual(
            [(d.inventory.opening, d.inventory.delivered, d.inventory.fresh_burn,
              d.inventory.stored_burn, d.inventory.closing, d.gas_heat_gj) for d in result.days],
            [(0, 3, 2, 0, 1, 0), (1, 0, 0, 1, 0, 6), (0, 0, 0, 0, 0, 3)],
        )
        validate_inventory((day.inventory for day in result.days), initial=0, terminal=0)

    def test_tighter_plant_capacity_makes_day_one_truck_unusable(self):
        days = (replace(self.scenario.days[0], plant_capacity_gj=3), *self.scenario.days[1:])
        result = solve(replace(self.scenario, days=days))
        self.assertEqual(result.status, 0)
        self.assertAlmostEqual(result.totals.cost, 54)
        self.assertEqual(tuple(day.trucks for day in result.days), (0, 0, 0))

    def test_tighter_storage_makes_day_one_truck_unusable(self):
        days = (replace(self.scenario.days[0], storage_capacity_m3=0), *self.scenario.days[1:])
        result = solve(replace(self.scenario, days=days))
        self.assertEqual(result.status, 0)
        self.assertAlmostEqual(result.totals.cost, 54)

    def test_insufficient_heat_capacity_is_infeasible(self):
        days = (self.scenario.days[0], replace(self.scenario.days[1], gas_capacity_gj=5),
                self.scenario.days[2])
        result = solve(replace(self.scenario, days=days))
        self.assertEqual(result.status, 2)
        self.assertIsNone(result.totals)
        self.assertEqual(result.days, ())

    def test_opening_stock_available_on_first_day_and_terminal_target_optional(self):
        scenario = replace(
            self.scenario, days=(DayLimits(3, 0, 1, 3, gas_capacity_gj=0),),
            initial_inventory_m3=2, terminal_inventory_m3=None,
        )
        result = solve(scenario)
        self.assertEqual(result.status, 0)
        self.assertEqual(result.days[0].inventory.stored_burn, 1)
        self.assertEqual(result.days[0].inventory.closing, 1)
        self.assertEqual(solve(replace(scenario, terminal_inventory_m3=0)).status, 2)

    def test_terminal_stock_cannot_be_used_as_opening_stock(self):
        scenario = replace(
            self.scenario, days=(DayLimits(9, 1, 3, 9, gas_capacity_gj=0),),
            terminal_inventory_m3=3,
        )
        self.assertEqual(solve(scenario).status, 2)

    def test_audit_rejects_corrupted_heat_and_cost(self):
        result = solve(self.scenario)
        days = (replace(result.days[0], gas_heat_gj=1), *result.days[1:])
        with self.assertRaisesRegex(ValueError, "heat does not match demand"):
            verify_plan(self.scenario, days)
        with self.assertRaisesRegex(ValueError, "disagrees with solver objective"):
            verify_plan(self.scenario, result.days, expected_cost=34)


if __name__ == "__main__":
    unittest.main()
