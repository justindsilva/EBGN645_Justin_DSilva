## Main results

Extra world-emission change from adding the CBAM on top of the EU ETS (MtCO2):

| Case | China price-taker | Cournot China | CBAM shuts off EU imports from |
|---|---|---|---|
| Base, 2030 (48.5% phase-in) | −0.073 | −0.053 | 2034 / 2032 |
| Base, 2034 (full strength) | −0.133 | −0.072 | |
| Supply elasticity 0.5, 2030 | −3.07 | −1.93 | 2034 / 2032 |
| Supply elasticity 0.5, 2034 | −6.29 | −2.92 | |

- **With a strategic China, the CBAM's extra emission cut is 27–54% smaller** 
- **The CBAM becomes prohibitive two years earlier** (2032 against 2034). Cournot China restricts supply to RoW, so RoW prices are higher, RoW steel reaches the EU at a higher cost, and a smaller charge closes the EU market.
- **EU imports get cleaner (resource shuffling) only when China is strategic.** In the price-taker counterfactual, China does not sell to the EU at all.


Leakage *rates* are also reported, but they are unstable here. Supply is almost perfectly inelastic (elasticity 0.01), so EU emissions barely move and the denominator of the leakage ratio is close to zero. The change in world emissions is the headline measure.

## How to run

From the repository root:


gams model/cbam_cournot.gms
gams model/cbam_phasein.gms



- `cbam_cournot.gms` solves 6 scenarios × 4 sensitivity cases at full CBAM strength and writes `results/cournot_results.csv`.
- `cbam_phasein.gms` solves the same model for each phase-in year from 2026 to 2034 and writes `results/phasein_results.csv`.


## Repository layout

| Path | Contents |
|---|---|
| `data/cbam_data.gms` | All 2025 inputs, with source workbook cell references |
| `data/cbam_steel_gams_inputs.xlsx` | Source workbook: raw data, sources, calibration checks |
| `model/cbam_calib.gms` | Demand and supply calibration (included by every model) |
| `model/cbam_cournot.gms` | Cournot model at full CBAM strength |
| `model/cbam_phasein.gms` | Cournot model over the 2026–2034 CBAM phase-in |
| `model/cbam_qp.gms`, `model/cbam_mcp.gms` | Competitive-calibration benchmark (QP and MCP) |

Scenario names: `c_` = China price-taker, `k_` = Cournot China; `_base` = no policy, `_ets` = EU carbon price, `_cbam` = carbon price + CBAM. 
Sensitivities: `base`, `etaS_mr` (supply elasticity 0.5), `e_alt` (alternative emission intensities), `epsD_lo` (demand elasticity −0.18).

## Main caveats

- **The supply elasticity of 0.5** in `etaS_mr` is an author assumption, not an estimate.
- **`epsD_lo` cannot be calibrated as full Cournot conduct.** The implied markup on China→RoW exceeds the observed price gap, so the case is reported as infeasible rather than solved.
- **Data gaps:** China→EU tonnage is the 2024 figure carried forward, the RoW price is a world export price used as a proxy, and two-way trade is netted.
- **Elasticities and process intensities** come from Clausing, Colmer, Hsiao & Wolfram, *The Global Effects of Carbon Border Adjustment Mechanisms* ([working paper](https://allanhsiao.com/files/CCHW_cbam.pdf)). Full data sources are in the workbook's `Sources` sheet.
