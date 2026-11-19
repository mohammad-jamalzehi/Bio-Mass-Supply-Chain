### The validation ladder

1. **Phase 1 · Hand calculations and arithmetic checks.** Track opening stock, receipts, fresh burn, stored burn and closing stock for one fuel. Reject examples that balance in aggregate but borrow from the wrong pool.
2. **Phase 2 · A small MILP.** Enumerate four possible integer truck schedules by hand and compare the winner with the HiGHS result.
3. **Phase 3 · Competing fuels and plants.** Audit each fuel ledger and each plant's useful heat. Compare the winning fuel under changes to storage and verify feasible and infeasible cases.
4. **Phase 4 · Historical-input experiments.** Load named scenarios, preserve unresolved data, and audit both daily and whole-horizon balances.

SciPy's `milp` passes a sparse constraint matrix, bounds and integrality flags to HiGHS. A feasible operating plan satisfies the equations. A proven optimum also has a matching solver bound within tolerance. The MIP gap measures remaining distance from that proof. A time-limited incumbent may still pass the physical audit. Validating the model against real operations also requires reliable demand, capacities, and emissions factors.
