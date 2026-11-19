# Biomass Supply Chain Optimization

> A mixed-integer programming model exploring the cost and emissions implications of biomass sourcing, transportation, storage, and heat production.

**Domain:** operations research · energy systems · supply chain analytics  
**Core tools:** GAMS · CPLEX · mixed-integer programming (MIP)

## The decision problem

A heat-supply system can use several biomass feedstocks, existing/new conversion capacity, and natural gas. Choosing delivery quantities, inventory levels, and how much material each system consumes has operational and environmental consequences.

The model represents:
- **365 planning periods** and **three feedstock categories** (`wc`, `wp`, `bq`, as named in the source).
- Integer numbers of daily truck deliveries and continuous biomass/energy flow decisions.
- Heat demand, storage capacity, truck limits, and conversion capacity constraints.
- Separate objective functions for **supply cost** (`z1`) and **emissions** (`z2`).

## Mathematical formulation at a glance

The model chooses transport, allocation, inventory and backup gas variables subject to demand and operational constraints.

| Component | Implementation |
| --- | --- |
| Cost objective | `objectivefunction1` / `z1` |
| Emissions objective | `objectivefunction2` / `z2` |
| Delivery decisions | `nt(b,t)` — integer truckloads |
| Biomass allocation | `xE`, `xN`, `xES`, `xNS` |
| Inventory and backup | `sb(b,t)`, `g(t)` |
| Main constraints | Demand, logistics, inventory and conversion capacity (`co1`–`co9`) |

**Important distinction:** the current GAMS file solves the cost and emissions objectives **independently** in two MIP runs. It does **not** implement or plot a Pareto frontier. A true multi-objective trade-off analysis is a logical next extension.

## Reproduce the model

1. Install a compatible **GAMS** distribution with access to the **CPLEX** MIP solver.
2. Open `Bio-mass Supply Chain.gms` in GAMS Studio.
3. Inspect the parameters and data defined at the top of the model.
4. Run the model using GAMS Studio or from a terminal:

```bash
gams "Bio-mass Supply Chain.gms"
```

5. Inspect the solution listing and the displayed objective values and decision variables. The repository also contains an original report and solver artifacts.

> The model has not yet been revalidated against a fresh solver/environment as part of this documentation update. Treat supplied historical output files as artifacts, not independently verified performance claims.

## Repository structure

```text
Bio-mass Supply Chain.gms   # GAMS model and scenario parameters
Bio-mass Supply Chain.pdf   # Academic report
Bio-mass Supply Chain.log   # Historical solver log
Bio-mass Supply Chain.mpx   # Original generated artifact
README.md                   # Project guide
```

## Engineering notes and next steps

- Add model/data validation, unit checks, and explicit inventory-balance tests.
- Revisit period-boundary handling and verify that units in all heat/emissions equations are consistent.
- Implement an **epsilon-constraint** or weighted-objective experiment to produce and visualize a Pareto frontier.
- Separate model inputs from formulation, then expose scenarios through a reproducible Python or web interface.

**Project status:** Academic model / research prototype; not a deployable decision-support application.
