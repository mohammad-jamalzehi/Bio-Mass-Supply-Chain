"""Solver-free checks for the proposed biomass accounting equations."""

from dataclasses import dataclass
from math import isclose, isfinite
from typing import Iterable


@dataclass(frozen=True)
class InventoryDay:
    """Volumes in m³ for one feedstock and one day."""

    opening: float
    delivered: float
    fresh_burn: float
    stored_burn: float
    closing: float


def nonnegative(name: str, value: float) -> None:
    if not isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative")


def validate_inventory(
    days: Iterable[InventoryDay], *, initial: float,
    terminal: float | None = None, tolerance: float = 1e-9,
) -> None:
    """Check one feedstock's ledger in day order, with tolerance in m³."""
    nonnegative("initial inventory", initial)
    nonnegative("tolerance", tolerance)
    if terminal is not None:
        nonnegative("terminal inventory", terminal)

    previous = initial
    for number, day in enumerate(days, start=1):
        for name, value in vars(day).items():
            nonnegative(f"day {number} {name}", value)
        if not isclose(day.opening, previous, rel_tol=0, abs_tol=tolerance):
            raise ValueError(f"day {number}: opening inventory differs from prior closing")
        if day.fresh_burn > day.delivered + tolerance:
            raise ValueError(f"day {number}: fresh burn exceeds deliveries")
        if day.stored_burn > day.opening + tolerance:
            raise ValueError(f"day {number}: stored burn exceeds opening inventory")
        expected = day.opening + day.delivered - day.fresh_burn - day.stored_burn
        if not isclose(day.closing, expected, rel_tol=0, abs_tol=tolerance):
            raise ValueError(f"day {number}: closing inventory does not conserve volume")
        previous = day.closing

    if terminal is not None and not isclose(previous, terminal, rel_tol=0, abs_tol=tolerance):
        raise ValueError("closing inventory differs from the terminal target")


def useful_heat_gj(
    volume_m3: float, wet_density_t_per_m3: float, moisture_fraction: float,
    heating_value_gj_per_dry_t: float, efficiency: float,
) -> float:
    """Convert delivered or stored volume to useful heat in GJ."""
    nonnegative("volume", volume_m3)
    if not isfinite(wet_density_t_per_m3) or wet_density_t_per_m3 <= 0:
        raise ValueError("wet density must be finite and positive")
    if not isfinite(heating_value_gj_per_dry_t) or heating_value_gj_per_dry_t <= 0:
        raise ValueError("dry heating value must be finite and positive")
    if not 0 <= moisture_fraction < 1:
        raise ValueError("wet-basis moisture fraction must be in [0, 1)")
    if not 0 <= efficiency <= 1:
        raise ValueError("efficiency must be in [0, 1]")
    dry_tonnes = volume_m3 * wet_density_t_per_m3 * (1 - moisture_fraction)
    return dry_tonnes * heating_value_gj_per_dry_t * efficiency


def daily_emissions_tco2e(
    biomass_heat_and_factors: Iterable[tuple[float, float]],
    *, gas_heat_gj: float, gas_factor_tco2e_per_gj: float,
) -> float:
    """Use factors expressed per GJ of useful heat and count gas once."""
    nonnegative("gas heat", gas_heat_gj)
    nonnegative("gas emissions factor", gas_factor_tco2e_per_gj)
    biomass_emissions = 0.0
    for heat, factor in biomass_heat_and_factors:
        nonnegative("biomass heat", heat)
        nonnegative("biomass emissions factor", factor)
        biomass_emissions += heat * factor
    return biomass_emissions + gas_heat_gj * gas_factor_tco2e_per_gj
