"""Run the reproducible synthetic fixture with python -B -m optimization.example."""

from optimization.phase2 import DayLimits, Fuel, Scenario, solve


def synthetic_scenario() -> Scenario:
    return Scenario(
        days=(DayLimits(6, 1, 1, 6), DayLimits(9, 0, 1, 6), DayLimits(3, 1, 1, 6)),
        fuel=Fuel(
            truck_volume_m3=3, wet_density_t_per_m3=0.5, moisture_fraction=0.2,
            heating_value_gj_per_dry_t=10, efficiency=0.75, cost_per_dry_t=5,
            emissions_tco2e_per_useful_gj=0.1,
        ),
        initial_inventory_m3=0, terminal_inventory_m3=0, gas_cost_per_useful_gj=3,
        gas_emissions_tco2e_per_useful_gj=0.5,
    )


def main() -> None:
    print("SYNTHETIC TEST FIXTURE — not campus results")
    result = solve(synthetic_scenario())
    print(f"Solver status {result.status}: {result.message}")
    if result.totals is None:
        raise SystemExit(1)
    print("Day  Trucks  Opening  Delivered  Fresh  Stored  Closing  Gas heat")
    print("                 m3         m3     m3      m3       m3       GJ")
    for index, day in enumerate(result.days, start=1):
        inv = day.inventory
        print(f"{index:3} {day.trucks:7} {inv.opening:8g} {inv.delivered:10g} "
              f"{inv.fresh_burn:6g} {inv.stored_burn:7g} {inv.closing:8g} {day.gas_heat_gj:9g}")
    print(f"Cost: ${result.totals.cost:.2f} | Best bound: ${result.best_bound:.2f} | Gap: {result.mip_gap:g}")
    print(f"Emissions accounting: {result.totals.emissions_tco2e:.2f} tCO2e")


if __name__ == "__main__":
    main()
