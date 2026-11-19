# Phase 3 · Multi-feedstock and multi-plant validation

This model tests wood chips (wc), wood pellets (wp), briquettes (bq),
and the existing and new biomass plants. All results in this phase are
**synthetic test fixtures**. The historical model and Phase 1–2 files
are separate. The only objective is cost.

## Reproduce

Reuse the pinned SciPy/HiGHS Python environment from Phase 2. From the
repository root, using its Python interpreter:

```powershell
python -m pip install -r requirements-phase2.txt
python -B -m unittest discover -s tests -v
python -B -m optimization.example_phase3 --case base
python -B -m optimization.example_phase3 --case tight-silo
python -B -m optimization.example_phase3 --case no-storage
```

The examples write JSON to stdout. They do not change model inputs or
generate files automatically.

## Equations and units

For fuel f, plant p and day t, let N[f,t] be a nonnegative integer truck
count, A[f,t] = N[f,t] × truck volume in m³, F[f,p,t] fresh burn and
U[f,p,t] stored burn in m³. O and I are opening and closing inventories
per fuel and day. G[t] is useful gas heat in GJ. The flows are nonnegative.

- **Fuel availability:** sum F over plants ≤ A. Sum U over plants ≤ O.
  Both plants share each fuel's arrivals and opening stock. Each plant
  checking availability alone could count the same stock twice.
- **Inventory:** I = O + A − sum(F + U) over plants. Opening stock on
  day 1 comes from the scenario. Later openings equal the previous day's
  closing stock. An optional final inventory target applies per fuel.
- **Shared limits:** sum N over fuels ≤ daily truck limit. Sum I over
  fuels ≤ total closing storage. An optional silo limit caps wp closing
  stock only. Here the silo is **included in** total storage and is a
  tighter bound on pellets. Its volume is never added to total storage.
  If the historical volumes are confirmed to be additive, revise the
  storage constraints explicitly before using those inputs.
- **Conversion:** useful GJ per m³ equals wet tonnes per m³ × dry mass
  fraction × GJ per dry tonne × efficiency[f,p]. Each plant's useful
  heat from fresh and stored burn across fuels must fit its own capacity.
  The two plants' useful heat plus G equals daily demand. Gas is capped
  only when the scenario supplies a gas capacity.
- **Cost and emissions:** charge dry tonnes delivered at the fuel's
  price and useful gas GJ at the gas price. Initial inventory incurs no
  delivery cost in this window. After solving, multiply each fuel-and-plant
  heat flow by its useful-heat emissions factor, then add gas emissions once
  per day. Emission factors have no effect on the dispatch.

The model assumes full trucks, constant fuel properties, no storage loss or
disposal, and end-of-day capacity accounting. It does not model receiving-time
peak storage. Technically, each decision maps to a named column in a sparse
matrix. The solver sees rows with upper or lower bounds. A separate verifier
reconstructs heat, available material, costs and emissions from the returned
plan and sends each fuel's ledger to the Phase 1 validator. The model only
publishes a plan on status 0 after checking the recomputed cost against the
solver objective and best bound with zero reported gap within tolerance.

## Synthetic cases and independent check

Only one truck may arrive on day 1 and none on day 2. Demand is 6 GJ
then 2 GJ. Existing capacity is 2 GJ both days. New capacity is 4 GJ
then 2 GJ. Initial and terminal stock are zero for every fuel. Both days
have 1 m³ total closing storage in the base case, including a 0.5 m³
pellet silo. Fuel input is 4 GJ/m³ before efficiency. Gas costs $3/GJ.
The toy emissions factors are 0.1 tCO₂e per useful biomass GJ for every
fuel and plant and 0.5 tCO₂e per useful gas GJ.

| Fuel | Truck volume m³ | Cost per truck | Useful GJ/m³ existing | Useful GJ/m³ new |
| --- | ---: | ---: | ---: | ---: |
| Chips | 2 | $2.40 | 2 | 3 |
| Pellets | 2.5 | $5.00 | 2 | 4 |
| Briquettes | 2 | $4.00 | 3 | 3 |

| Day-1 choice | Base cost | Silo tightened to 0.4 m³ | No day-1 storage |
| --- | ---: | ---: | ---: |
| No truck | $24.00 | $24.00 | $24.00 |
| Chips | $8.40 | $8.40 | $10.40 |
| Pellets | $5.00 | infeasible | infeasible |
| Briquettes | $10.00 | $10.00 | $10.00 |

Chips can supply at most 6 GJ if all 2 m³ reach the new plant. This is
possible by burning 4/3 m³ there on day 1 and storing 2/3 m³ for day 2.
With no storage, the remaining 2/3 m³ must burn at the existing plant,
so useful heat drops to 16/3 GJ. Briquettes produce 6 GJ from a full truck
at either plant; the day-1 capacities can absorb all of it. Pellets can
burn at most 2 m³ on day 1 because each plant accepts at most 1 m³.
The remaining 0.5 m³ must be stored and then burned in the new plant on
day 2. This exactly meets demand and binds both day-1 biomass capacities.
Silo space below 0.5 m³ makes this choice infeasible. These four choices
exhaust all integer truck schedules in each fixture.

## Historical inputs requiring evidence

The prior [Phase 1 audit](phase-1-formulation.md) compares the original
[model](../Bio-mass%20Supply%20Chain.gms) with the
[report](../Bio-mass%20Supply%20Chain.pdf). The source assigns 12,624 GJ/day
(about 146 MW average); the report's data tables do not establish a matching
daily figure. Existing capacity of 518.4 GJ/day is equivalent to 6 MW over
24 hours. The report narrative says a 12 MW addition and 18 MW combined,
while its Table 1 and the source label 1,555.2 GJ/day, equivalent to 18 MW,
as new-system capacity. The source's wood/gas factors are 0.96/0.5;
the report lists 0.094/0.056. The factor values, heat basis, and accounting
scope remain unresolved. The source also applies its 192 m³ silo limit
to every fuel despite its pellet description. Confirm its role relative to
the 412 m³ total. Annual demand, gas availability, and initial and terminal
inventory policy remain open. None of these data are used in the fixtures.

## Structured result

Version 1 JSON identifies the scenario, feedstocks, plants, units and
synthetic provenance, then records solver status, backend, cost bound and
MIP gap. Each day has limits, trucks by fuel, inventories, fresh and
stored volumes by plant and fuel, useful heat by plant and fuel, gas
heat, cost, and emissions. Totals sum the audited daily amounts. On
infeasible termination the result contains solver diagnostics, null totals,
and an empty day list. A dashboard can render the audited daily values
without reconstructing the solver or inferring any campus results.
