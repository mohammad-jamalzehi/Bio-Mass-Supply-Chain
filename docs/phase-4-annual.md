# Phase 4 · Historical input reconciliation

**Authors:** Mohamad Jamalzehi and Kosar Ghazali.

These experiments use the corrected Phase 3 mathematics. Agreement between
the source and the report means transcription agreement, not validation of
campus operations. The GAMS file and report are unchanged.

## Parameter reconciliation

Feedstock triplets are in wood chips, pellets, briquettes order. The PDF
narrative is on page 1, Tables 1–2 on page 2, and Tables 3–5 on page 3.

| Input and unit | GAMS source | Academic report | Experiment treatment |
| --- | --- | --- | --- |
| Horizon | 365 days | Fuel properties held constant for 365 days | 365 daily periods |
| Campus demand, useful GJ/day | 12,624 every day, about 146.1 MW average | No corroborating daily demand in its data tables | All supplied variants retain 12,624 as an unverified source input. |
| Purchase price, source $/dry tonne | 20 / 90 / 136 | Same | Retained. Currency base year and procurement scope remain open. |
| Gas price, source $/useful GJ | 37 | 37 | Retained. Gas supply limit is unspecified and modeled uncapped. |
| Existing / new efficiencies | (0.735, 0.763, 0.761) / (0.768, 0.811, 0.807) | Same | Used as dimensionless useful-heat conversion factors. |
| Truck volume, m³ per full truck | 110 / 110 / 110 | Same | Full truck interpretation retained. |
| Shared truck limit, trucks/day | 6 | 6 | Retained. |
| Total / pellet silo storage, m³ | 412 / 192 | Same volumes | Use a pellet-only cap inside total closing storage. The original equation applied 192 to each fuel. Whether the real capacities are additive remains unresolved. |
| Existing conversion capacity | 518.4 GJ/day | 518.4 GJ and 6 MW | Same if available continuously for 24 hours. |
| New conversion capacity | 1,555.2 GJ/day assigned to new alone | Table 1 gives 1,555.2. Narrative says a 12 MW addition and 18 MW combined. | Base retains table number. Separate narrative alternative sets 1,036.8 GJ/day for new. |
| Heating values, GJ/dry tonne | 12.5 / 18 / 19 | Same | Retained constant. |
| Moisture fractions | 0.30 / 0.10 / 0.12 | Same | Interpreted on a wet-mass basis. |
| Delivered density, wet tonnes/m³ | 0.15 / 0.75 / 1.2 | Same | Retained constant. |
| Biomass / gas emissions factors, declared tCO₂e/GJ heat | 0.96 / 0.5 | 0.094 / 0.056 | Separate factor alternatives. Basis and scope unverified. Provisional per useful GJ and one biomass factor for both plants. |
| Carbon price, allowance, compliance target | 50 $/tCO₂e, 2,000 tCO₂e, 2,000 tCO₂e | Same | Documented only. The original policy settings and four deviation placeholders are inactive in the solved objectives. |
| Opening and terminal stocks | No trustworthy explicit boundary targets | Unresolved | Explicit experimental zero stocks per fuel at both endpoints. |

## Reproduce

Use the pinned Phase 2 Python dependencies. Run from the repository root.

```powershell
python -m pip install -r requirements-phase2.txt
python -B -m unittest discover -s tests -v
python -B -m optimization.annual --days 3 --time-limit 15
python -B -m optimization.annual --days 7 --time-limit 30
python -B -m optimization.annual --scenario gams-historical-inputs --time-limit 30 --output runs/gams
python -B -m optimization.annual --scenario report-narrative-capacity --time-limit 30 --output runs/report-capacity
python -B -m optimization.annual --scenario report-table-emissions --time-limit 30 --output runs/report-factors
```

The [scenario file](../scenarios/annual.json) holds inputs separately from
solver code. Feedstock lists are positional in this order: full truck volume,
wet density, moisture, dry heating value, and price per dry tonne. Daily
overrides use day numbers as strings and replace named fields. The report
alternatives change one assumption at a time. Both retain the unconfirmed
GAMS demand. The emissions alternative has no effect on the objective or
constraints. Alternative optimal plans can still differ if there are ties.
The scenario SHA-256 includes the effective inputs, assumptions and horizon.
A shortened experiment imposes its final stock policy at the shorter cutoff.

The CLI calls the unchanged sparse Phase 3 MILP across the requested days.
In its dedicated process a scoped solver wrapper adds a time limit and retains
the raw SciPy outcome. Prices are temporarily represented in thousands of
source dollars to improve numerical summation and are converted back for audit
and export. Negative solver values within the existing Phase 3 tolerance of
zero are set to zero before audit. Larger violations remain visible. No
physical constraint changes.

For both optimal outcomes and time-limited incumbents, the audit checks each
day's Phase 3 constraints, reconciles whole-horizon material balances by fuel,
and independently recomputes cost, heat, and emissions. The solver objective
must match the recomputed cost. On optimal termination the bound must match it
too. A separate incumbent decoder reads the original Phase 3 variable order,
checked against the Phase 3 plan on optimal runs.

## Exports and statuses

Specify an output directory to produce `result.json` and `flows.csv`.
JSON includes inputs and assumptions, scenario digest, solver version and
elapsed seconds, objective, bound, gap, validation outcome, and audited
plans when available. A null solver field means unavailable. CSV has one
inventory/delivery row per fuel and day, one burn/heat row per plant and fuel
per day, and one gas row per day. The gas row also carries the daily demand,
cost, and emissions totals once. A run without an audited plan exports the
CSV header only. The CLI returns 0 for an audited optimum and 2 for other
outcomes after recording diagnostics. Runtime is machine-dependent.

## Recorded validation runs

The recorded JSON and CSV files are under [results/phase4](../results/phase4).
They use SciPy 1.18.1 with HiGHS 1.12.0. Costs are source dollars.

| Scenario | Days | Termination | Audited cost | Best bound | MIP gap |
| --- | ---: | --- | ---: | ---: | ---: |
| GAMS inputs | 3 | Optimal | 1,206,354.90 | 1,206,354.90 | 0% |
| GAMS inputs | 7 | Optimal | 2,811,863.49 | 2,811,863.49 | 0% |
| GAMS inputs | 365 | 30 s time limit | 146,497,618.57 | 146,455,283.18 | 0.02890% |
| Report narrative plant capacity | 365 | 30 s time limit | 152,217,297.31 | 152,200,046.36 | 0.01133% |
| Report table emissions factors | 365 | 30 s time limit | 146,497,618.57 | 146,455,283.18 | 0.02890% |

Every annual incumbent passed the physical and cost audit. Optimality was
not proven on the 365-day runs. The base and report-factor runs produced
the same cost and heat totals in these runs. Their emissions differ solely
because the factors used for accounting differ. These outputs describe the
specified experiments. The input assumptions and emissions accuracy remain
unverified.
