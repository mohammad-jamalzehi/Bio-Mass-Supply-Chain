# Phase 2 · Small-scale cost MILP validation

**Authors:** Mohamad Jamalzehi and Kosar Ghazali.

This phase uses a **synthetic test fixture** with one fuel, one biomass
plant, and three days. The historical campus model, academic PDF, and
Phase 1 validator are unchanged. No campus optimization or Pareto analysis
is part of this implementation.

## Run the fixture

From the repository root, create a virtual environment and install the
pinned NumPy and SciPy versions. The commands below use PowerShell and
a Python interpreter compatible with both pinned wheels. Python 3.15 was used
for verification.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-phase2.txt
.\.venv\Scripts\python.exe -B -m optimization.example
.\.venv\Scripts\python.exe -B -m unittest discover -s tests -v
```

On macOS or Linux, use `.venv/bin/python` for the last three commands. The
filter `-p test_phase1.py` after `-s tests` runs only the solver-free tests.
The `-B` option avoids adding Python bytecode caches.

SciPy's `milp` wraps open-source HiGHS and takes arrays of costs, bounds,
and linear constraints. It reports a termination code, primal cost, best
bound, and MIP gap. Direct `highspy` is a good option if we later need
lower-level controls and richer solver integration. The compact SciPy API
is enough for this fixture, and its pinned wheel brings its HiGHS backend
without a separate GAMS, CPLEX, or highspy install.
[ SciPy milp API ](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.milp.html)
[ HiGHS Python API ](https://ergo-code.github.io/HiGHS/stable/interfaces/python/)

## Decisions and constraints

There are five variables per day. The truck count N is a nonnegative integer.
Fresh burn F, stored burn U, closing inventory I (all m³), and useful gas heat
G (GJ) are nonnegative continuous variables. Delivered volume A is derived as
truck count times volume per truck. Opening inventory O is given on day 1
and equals the previous day's closing stock thereafter.

| Rule | Equation or bound |
| --- | --- |
| Inventory conservation | I = O + A − F − U |
| Fresh and stored fuel availability | F ≤ A and U ≤ O |
| Heat balance | q × (F + U) + G = daily demand |
| Daily operating limits | N ≤ truck limit, I ≤ closing storage capacity, and q × (F + U) ≤ biomass plant capacity |
| Optional scenario limits | G ≤ gas capacity if supplied, and final I = terminal target if supplied |

The useful heat coefficient q multiplies wet tonnes per m³, dry fraction,
GJ per dry tonne, and plant efficiency. The demand and capacity rows use
the same q. Cost is dry tonnes purchased times biomass price plus useful
gas GJ times gas price, summed over the days. Opening stock has no purchase
charge in this planning window. Emissions are calculated afterward as useful
biomass heat times its factor plus gas heat times its factor, once per day.

The implementation assumes full trucks, no storage loss, constant fuel
properties, and end-of-day storage limits, as proposed in Phase 1. The
optional gas limit and day-specific bounds are explicit scenario inputs.
The annual demand profile, emissions factor values and basis, separate
equipment capacities, gas availability, and annual terminal stock target
remain unresolved. The fixture values are independent of those decisions.

## Independent hand calculation

The synthetic fuel gives 0.5 × (1 − 0.2) × 10 × 0.75 = 3 useful GJ/m³.
A 3 m³ truck costs $6 because it brings 1.2 dry tonnes at $5 per tonne.
Gas costs $3 per useful GJ. Demand is 6, 9, and 3 GJ. The truck limits are
1, 0, and 1 per day. Closing storage is capped at 1 m³, the plant at
6 GJ/day, and both initial and final stock are zero. Gas is uncapped in the
base fixture. The synthetic reporting factors are 0.1 and 0.5 tCO₂e per
useful GJ of biomass and gas heat respectively.

| Truck schedule | Hand calculation |
| --- | --- |
| (0, 0, 0) | 18 GJ of gas costs $54. |
| (1, 0, 0) | Day 1 uses 2 m³ and stores 1 m³. Day 2 burns that 1 m³. Day 2 and day 3 need 6 and 3 GJ of gas. The truck costs $6 and gas costs $27, for $33. |
| (0, 0, 1) | Infeasible. The day-3 truck brings 3 m³, day-3 demand can absorb at most 1 m³, and there is no terminal stock or disposal. |
| (1, 0, 1) | Infeasible for the same day-3 reason, even if stock arrives from day 2. |

This covers every integer truck decision. For (1, 0, 0), 9 useful GJ from
biomass plus 9 GJ from gas meet 18 GJ of demand. The synthetic emissions
total is 9 × 0.1 + 9 × 0.5 = 5.4 tCO₂e. This factor calculation has no effect
on dispatch because cost is the sole objective.

## How the code checks the solve

`optimization/phase2.py` separates scenario inputs, solver construction,
and plan verification. Each constraint row contains coefficients for a
single rule. The previous-day stock coefficient links consecutive rows.
Only truck columns receive an integrality flag. Delivery is derived from
integer trucks, so it cannot drift away from that value. The dense matrix
is deliberate for the tiny validation case.

After optimal termination, the code checks integrality before rounding
truck counts and invokes the Phase 1 ledger validator. It separately
recalculates demand, capacities, delivery quantities, total cost, and
emissions from the returned decisions. It compares the recomputed cost to
the solver's objective and bound within an absolute tolerance of 1e-7.
The fixture independently establishes the optimal cost by hand.
A nonoptimal status returns diagnostics without an endorsed operating plan.
The test suite also exercises first-day opening stock, optional terminal
stock, insufficient heat capacity, and tighter storage and plant limits.
