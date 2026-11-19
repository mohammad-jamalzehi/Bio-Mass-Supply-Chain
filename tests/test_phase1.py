import unittest

from validation.phase1 import (
    InventoryDay,
    daily_emissions_tco2e,
    useful_heat_gj,
    validate_inventory,
)


# Artificial ledger for a single feedstock. These are hand-picked inputs.
THREE_DAYS = (
    InventoryDay(0, 10, 4, 0, 6),
    InventoryDay(6, 0, 0, 5, 1),
    InventoryDay(1, 7, 7, 1, 0),
)


class InventoryTests(unittest.TestCase):
    def test_three_day_ledger(self):
        validate_inventory(THREE_DAYS, initial=0, terminal=0)

    def test_fresh_burn_cannot_use_opening_inventory(self):
        # Total volume balances. The wrong source still makes this impossible.
        with self.assertRaisesRegex(ValueError, "fresh burn exceeds deliveries"):
            validate_inventory([InventoryDay(4, 5, 6, 0, 3)], initial=4)

    def test_stored_burn_cannot_use_todays_deliveries(self):
        with self.assertRaisesRegex(ValueError, "stored burn exceeds opening"):
            validate_inventory([InventoryDay(2, 10, 0, 3, 9)], initial=2)

    def test_first_day_respects_supplied_initial_inventory(self):
        with self.assertRaisesRegex(ValueError, "day 1: opening inventory"):
            validate_inventory(THREE_DAYS, initial=1)

    def test_next_opening_must_equal_prior_closing(self):
        with self.assertRaisesRegex(ValueError, "day 2: opening inventory"):
            validate_inventory(
                [THREE_DAYS[0], InventoryDay(5, 0, 0, 5, 0)], initial=0
            )

    def test_unbalanced_closing_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "does not conserve volume"):
            validate_inventory([InventoryDay(0, 10, 4, 0, 7)], initial=0)

    def test_terminal_target_is_explicit_and_optional(self):
        validate_inventory(THREE_DAYS[:-1], initial=0)
        with self.assertRaisesRegex(ValueError, "terminal target"):
            validate_inventory(THREE_DAYS[:-1], initial=0, terminal=0)

    def test_negative_and_nonfinite_values_are_rejected(self):
        for bad in (-1, float("nan"), float("inf")):
            with self.subTest(bad=bad):
                with self.assertRaisesRegex(ValueError, "delivered must be finite"):
                    validate_inventory([InventoryDay(0, bad, 0, 0, 0)], initial=0)


class HeatAndEmissionsTests(unittest.TestCase):
    def test_fixture_heat(self):
        for day, expected in zip(THREE_DAYS, (12, 15, 24), strict=True):
            with self.subTest(expected=expected):
                self.assertAlmostEqual(
                    useful_heat_gj(day.fresh_burn + day.stored_burn, 0.5, 0.2, 10, 0.75),
                    expected,
                )

    def test_heat_rejects_percent_in_place_of_fraction(self):
        with self.assertRaisesRegex(ValueError, "moisture fraction"):
            useful_heat_gj(4, 0.5, 20, 10, 0.75)

    def test_gas_is_counted_once_regardless_of_biomass_rows(self):
        for rows in ([], [(0, 0.1)] * 3, [(0, 0.1)] * 7):
            with self.subTest(rows=len(rows)):
                self.assertAlmostEqual(
                    daily_emissions_tco2e(rows, gas_heat_gj=8, gas_factor_tco2e_per_gj=0.5),
                    4,
                )

    def test_biomass_and_gas_use_separate_factors(self):
        self.assertAlmostEqual(
            daily_emissions_tco2e(
                [(12, 0.1), (15, 0.2)], gas_heat_gj=8, gas_factor_tco2e_per_gj=0.5
            ),
            8.2,
        )


if __name__ == "__main__":
    unittest.main()
