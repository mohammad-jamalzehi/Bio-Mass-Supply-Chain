The equations below describe the **corrected Python model** in Phase 3. Each equation links to its implementation. Phase 2 uses the single-fuel, single-plant special case.

### Sets and parameters

| Symbol | Meaning | Unit |
|---|---|---|
| $f\in\mathcal F$ | Wood chips, wood pellets, briquettes | index |
| $p\in\mathcal P$ | Existing or new conversion plant | index |
| $t\in\{1,\ldots,T\}$ | Ordered planning days | day |
| $V_f,\rho_f,m_f,h_f$ | Full-truck volume, wet density, wet-basis moisture, dry heating value | m³/truck, wet t/m³, fraction, GJ/dry t |
| $\eta_{pf},q_{pf}$ | Conversion efficiency and useful heat per unit volume | fraction, GJ/m³ |
| $c_f,c_g$ | Biomass purchase price and gas cost | source currency/dry t, source currency/GJ |
| $e_{pf},e_g$ | Biomass and gas factors per useful heat | tCO₂e/GJ |
| $D_t,L_t,K_{pt},S_t,S^{wp}_t$ | Heat demand, truck limit, plant capacity, total storage and optional pellet silo limit | GJ, trucks, GJ, m³, m³ |
| $I_f^0,I_f^T,\bar G_t$ | Initial stock, optional terminal target and optional gas cap | m³, m³, GJ |

### Decision variables and derived quantities

| Symbol | Domain | Unit | Meaning |
|---|---|---|---|
| $N_{ft}$ | $\mathbb Z_{\ge 0}$ | trucks | Full trucks received |
| $F_{fpt},U_{fpt}$ | $\mathbb R_{\ge 0}$ | m³ | Burn from fresh receipts and opening stock |
| $I_{ft}$ | $\mathbb R_{\ge 0}$ | m³ | Closing inventory |
| $G_t$ | $\mathbb R_{\ge 0}$ | GJ | Useful gas heat |
| $A_{ft},O_{ft},H_{fpt}$ | Derived | m³, m³, GJ | Deliveries, opening stock and useful biomass heat |

### Objective and constraints

Equations (1)–(4) define heat, cost, and emissions accounting. Equations (5)–(13) set the operating rules. Optional bounds apply only when supplied by the scenario.
