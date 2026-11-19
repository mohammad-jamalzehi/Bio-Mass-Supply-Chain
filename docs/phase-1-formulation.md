# Phase 1 · Proposed biomass accounting formulation

The original [GAMS model](../Bio-mass%20Supply%20Chain.gms) and
[academic report](../Bio-mass%20Supply%20Chain.pdf) are historical references.
This proposal and its arithmetic checks do not solve the model. The examples
below use artificial inputs. No annual operating plan has been computed.

## Inventory and physical availability

Let A[b,t] = nt[b,t] × VC[b] be delivered volume, F[b,t] = xE[b,t] +
xN[b,t] be burn from same-day deliveries, U[b,t] = xES[b,t] + xNS[b,t]
be burn from opening stock, and I[b,t] = sb[b,t] be closing stock.
All four quantities are measured in m³ per feedstock and day.
A keeps the source's full-truck assumption.

In the original model, co2 says U[t] = I[t], co3 says A[t] − F[t] =
I[t], and co4 says I[t−1] + A[t] − F[t] = I[t]. Subtracting co3 from
co4 gives I[t−1] = 0. Thus days 1–364 cannot retain stock, and co2 can
credit stored burn on day 365 against inventory that is still in closing stock.
GAMS uses a linear lag for t−1, so day 1 has no preceding index.

Replace those three equations with the following proposal for every b,t.

1. **Conservation** · I[b,t] = O[b,t] + A[b,t] − F[b,t] − U[b,t].
2. **Fresh availability** · F[b,t] ≤ A[b,t].
3. **Stored availability** · U[b,t] ≤ O[b,t].
4. **Opening stock** · O[b,1] = initial[b] and O[b,t] = I[b,t−1] for t > 1.
   Require all flows and inventory to be nonnegative.

The balance accounts for every cubic metre. The availability inequalities keep
the bookkeeping categories honest. For example, an overall balance can still
work if fresh burn borrows from opening stock. The fresh constraint rules that
out. A fresh delivery is either burned as fresh fuel or carried into closing
inventory. It cannot be counted as opening stock that same day.

Use an explicit initial inventory input. Zero can be tested as a scenario, but
the case-study opening stock remains unknown. For a future annual solve, choose
and document unrestricted closing inventory, a reserve, or cyclic inventory
I[b,365] = initial[b]. The fixture below imposes zero closing stock only for
that example. The proposal keeps end-of-day storage accounting and assumes no
volume loss or moisture change in storage. A receiving-time capacity limit or
storage loss would require additional data and a model change. Preserve the
daily truck limit and the total storage limit. Confirm whether the pellet silo
is included in total storage and apply its limit only to pellets if its
description is confirmed. The original co7 applies the silo limit to every fuel.

## Useful heat and emissions

For each fuel and each conversion unit, multiply

    m³ × wet tonnes/m³ × (1 − wet-basis moisture fraction)
       × GJ/dry tonne × conversion efficiency = useful GJ

Volume times density gives wet tonnes. The moisture adjustment gives dry
tonnes. Heating value gives fuel energy, and efficiency gives useful heat.
Use the same coefficient for a flow in the demand balance and its plant
capacity constraint. Original co1 omits density from xE, while co9 includes it.
For wood chips with the declared density of 0.15 tonnes/m³, this omission
makes the original demand equation credit about 6.67 times the heat counted
toward capacity for the same fresh existing-system flow.

Define H[b,u,t] as useful heat from fresh and stored fuel for each system u.
Proposed demand balance is sum over b,u of H[b,u,t] + g[t] = HD[t].
The capacity for each system limits the sum over fuels of the same H.
With the source's constant fuel properties, both burn streams use the same
coefficient for a given fuel and system. If properties vary by delivery date,
track inventory quality or receipt cohorts before attributing properties to
withdrawals. Using yesterday's properties for every withdrawal would assume
all stock came from yesterday.

For factors confirmed per GJ of **useful heat**, calculate annual emissions as
sum over b,u,t of H[b,u,t] × biomass_factor[b,u], plus a separate sum
over t of g[t] × gas_factor. The source puts the gas term inside sum((t,b), …),
so the same gas heat is counted once for wc, wp, and bq. At 0.5 tonnes/GJ,
this gives an effective rate of 1.5 tonnes/GJ in the original expression.
The proposal changes the summation only. It does not pick replacement factors.
If evidence shows a factor is per GJ of fuel input or on a lifecycle basis,
document the scope and convert the heat basis before using it. The source
uses the same wood factor for both systems even though its description names
the existing system. Whether that is appropriate remains open.

The source's purchase objective can be retained for a comparable baseline.
Purchased dry tonnes are delivered volume × wet density × dry fraction.
The gas cost is useful gas heat × price per useful GJ as declared. Freight,
losses, lifecycle emissions, carbon policy, and gas capacity have no new
assumptions in this phase.

## Unresolved source discrepancies

| Input | GAMS source | Academic report | What needs confirmation |
| --- | --- | --- | --- |
| Campus heat demand | 12,624 GJ for each day, equivalent to about 146.1 MW average | No matching numerical daily demand found in the report's data tables | Obtain original demand data, units, time aggregation, and whether demand is a flat scenario. |
| Existing capacity | 518.4 GJ/day, equivalent to 6 MW over 24 hours | 6 MW in the case description and 518.4 in Table 1 | Confirm the 24-hour conversion represents available production per modeled day. |
| New capacity | 1,555.2 GJ/day as a separate new-unit limit, equivalent to 18 MW | Narrative says a 12 MW addition and 18 MW combined. Table 1 also lists 1,555.2 for the new system. | Confirm whether 1,555.2 is meant to be combined capacity. A 12 MW increment would be 1,036.8 GJ/day at continuous operation. |
| Biomass emissions factor | 0.96, described as tonnes CO₂e per GJ of heat | Table 1 gives 0.094 | Confirm numerical value, numerator and denominator, system boundary, and whether both units use the same factor. |
| Gas emissions factor | 0.5, described as tonnes CO₂e per GJ of heat | Table 1 gives 0.056 | Confirm numerical value and the heat basis separately from correcting the repeated sum. |

The capacity comparison uses 1 MW × 24 hours = 86.4 GJ. It compares installed
power with daily energy under continuous operation, without assuming that
continuous operation is actually available.

## Three-day hand check

One artificial fuel has a full truck volume of 1 m³, so the delivery volumes
below correspond to 10, 0 and 7 trucks. This fixture does not use the case-study
truck limits. Its sole purpose is to check material and heat accounting.

| Day | Opening m³ | Delivered m³ | Fresh burn m³ | Stored burn m³ | Closing m³ | Useful biomass heat GJ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0 | 10 | 4 | 0 | 6 | 12 |
| 2 | 6 | 0 | 0 | 5 | 1 | 15 |
| 3 | 1 | 7 | 7 | 1 | 0 | 24 |

Toy fuel properties are 0.5 wet tonnes/m³, 20% wet-basis moisture,
10 GJ/dry tonne, and 75% efficiency. Each cubic metre supplies
0.5 × 0.8 × 10 × 0.75 = 3 GJ. On day 2, the five cubic metres must come
from opening stock because no truck arrives. Across the fixture, 17 m³ arrive
and 17 m³ are burned. If day 1 also has 8 GJ of gas heat at a toy factor of
0.5 tonnes CO₂e/GJ, the gas contribution is 4 tonnes CO₂e regardless of how
many feedstock rows are present. These are hand calculations, not solver output.

## Reproducible verification

From the repository root, run:

    python -B -m unittest discover -s tests -v

Use Python 3.10 or later. The standard-library tests require no GAMS, solver,
or package installation. They check continuity, independent availability,
conservation, initial and optional terminal stock, useful heat, and the
gas-emissions regression. The validator accepts one feedstock ledger at a time.
It does not check multi-fuel storage, plant capacity, or annual feasibility.

After the unresolved inputs are reviewed, encode the proposal as a separate
model. If GAMS/CPLEX is unavailable, translate the formulation using a Python
interface to open-source HiGHS. Start with a pinned environment and a small
model with fixture flows fixed, then check both infeasible availability cases,
a binding storage or plant limit, and an enumerated tiny integer instance.
Compare independent recomputation of balances and objectives with the solver
output. For later annual runs, record the input digest, solver version and
settings, model and solve status, objective bounds and gaps, and exported
decisions. An annual solve requires a functioning solver. Pareto analysis is
outside Phase 1.
