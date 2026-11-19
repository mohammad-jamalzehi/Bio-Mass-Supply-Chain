### Historical inputs · Original GAMS and academic report

Values transcribed from GAMS and compared with the academic report (Tables 1–5). Agreement means agreement between documents. It does not establish operational accuracy.

| Feedstock | Cost per dry t | Wet t/m³ | Moisture | GJ/dry t | Existing efficiency | New efficiency |
|---|---:|---:|---:|---:|---:|---:|
| Wood chips | 20 | 0.15 | 30% | 12.5 | 0.735 | 0.768 |
| Wood pellets | 90 | 0.75 | 10% | 18 | 0.763 | 0.811 |
| Briquettes | 136 | 1.20 | 12% | 19 | 0.761 | 0.807 |

Both sources agree on the values in the feedstock table. Costs use the source’s dollar denomination, with base year and currency treatment unresolved.

| Parameter | Historical value | Evidence and status |
|---|---|---|
| Truck volume | 110 m³ per fuel | GAMS and Table 2 agree. Original model treats deliveries as full trucks. |
| Receiving | 6 trucks/day | GAMS and Table 1 agree. |
| Storage and pellet silo | 412 m³ total, 192 m³ silo | GAMS and Table 1 agree on volumes. Original equation applies 192 m³ to all fuels. Corrected scenarios use a pellet-only cap inside total storage. Site interpretation remains unconfirmed. |
| Existing plant | 518.4 GJ/day | Table 1 agrees. Equivalent to 6 MW operating continuously. Availability unconfirmed. |
| New plant | 1,555.2 GJ/day | Table 1 agrees, but the narrative describes a 12 MW addition and 18 MW combined. The narrative alternative uses 1,036.8 GJ/day for the new plant. |
| Campus heat demand | 12,624 GJ/day over 365 periods | GAMS only. Magnitude and daily profile are unresolved. |
| Gas | 37 source dollars/GJ | GAMS and Table 1 agree on price. Supply cap is unspecified. |
| Emissions | GAMS 0.96 biomass / 0.5 gas; report 0.094 / 0.056 | Conflicting factors. Both pairs are used provisionally per useful GJ in separately labeled experiments. Units and lifecycle scope remain unresolved. |
| Carbon policy inputs | Price 50, allowance 2,000, compliance target 2,000 | GAMS and Table 1 agree. The active source solves do not use them. |
| Opening and terminal stock | No reliable historical targets | Corrected annual experiments assume zero explicitly. |

### Synthetic inputs · Validation fixtures

| Fixture | Explicit test values |
|---|---|
| Phase 2 | One fuel and one plant, three days. Demands 6 / 9 / 3 GJ; truck limits 1 / 0 / 1; 3 m³ per truck; 1 m³ storage; 6 GJ/day plant cap. Fuel: density 0.5, moisture 20%, heating value 10, efficiency 0.75, price 5/dry tonne. Gas price 3/GJ. Zero opening and final stock. |
| Phase 3 | Two days, three fuels and two plants. Demands 6 / 2 GJ and truck limits 1 / 0. Existing capacity 2 / 2 GJ, new 4 / 2 GJ. Truck volumes 2 / 2.5 / 2 m³, prices 3 / 5 / 5 per dry tonne, density 0.5, moisture 20%, heating value 10. Existing efficiencies 0.5 / 0.5 / 0.75, new 0.75 / 1 / 0.75. Gas price 3/GJ. Zero opening and final stock. |
| Phase 3 storage cases | Base total storage 1 m³ and pellet sublimit 0.5 m³. Alternatives reduce the first-day silo to 0.4 m³ or first-day total storage to zero. |
| Synthetic emissions | 0.1 biomass and 0.5 gas tCO₂e/useful GJ. Chosen for accounting tests. |
