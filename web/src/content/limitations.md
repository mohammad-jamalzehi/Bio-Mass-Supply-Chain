### Historical formulation findings

The original GAMS model solves cost and emissions separately. It cannot be read as a Pareto frontier. Its inventory equations force prior-day closing stock to zero, its fresh existing-plant term omits density in the demand equation, and its gas emission term sits inside a sum over three fuel types. Its final display can also combine a cost level from the first solve with decisions from the second.

### Research extensions

- Confirm the campus demand profile, capacity availability, and whether the silo is additive or included. Obtain a documented gas limit and defensible initial and terminal inventory targets.
- Establish the emissions measurement boundary and factor basis. Compare provisional assumptions and use primary data before drawing environmental conclusions.
- Run longer annual proof searches and report incumbent, bound and gap. Current Phase 4 exports establish feasibility under stated assumptions, with annual optimality unproven.
- Once the inputs are defensible, build an ε-constraint study by minimizing cost under separately chosen emissions ceilings. Evaluate both outcomes for each saved plan and remove dominated points. A frontier is future work.
- Extend cautiously to seasonal demand, uncertain deliveries, peak receiving storage or fuel-quality tracking.
