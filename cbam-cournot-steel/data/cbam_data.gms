$onText
===============================================================================
 cbam_data.gms -- Base-year (2025) data for the 3-region steel CBAM model
-------------------------------------------------------------------------------
 Every number below is transcribed by hand from cbam_steel_gams_inputs.xlsx
 (cached values, workbook prepared 7 Oct 2026). The cell each number comes from
 is given in the comment next to it. Source IDs (S1..S13, A1..A4) refer to the
 workbook's Sources sheet.

 Where a workbook value is itself a simple formula of raw inputs (emission
 intensities, the EU carbon price in USD), the raw inputs are defined here and
 the formula is evaluated in GAMS, then asserted against the workbook result.

 Units: Mt, USD/t, USD/tCO2, tCO2/t. Money in million USD.
 Implements CLAUDE.md Section 4 (data) for the toy model of Milestone 2.
===============================================================================
$offText

Set i "regions" / EU "European Union (27)", CN "China", RW "Rest of World" /;
Alias (i, j, k);

*-------------------------------------------------------------------------------
* 1. Prices, production, consumption (base year 2025)
*-------------------------------------------------------------------------------
Parameter
* p0: EU = W. Europe ex-works, CN = China ex-works, RW = world export FOB
*     (proxy: no single RoW price exists)
   p0(i) "base price, 2025 avg of 22 SteelBenchmarker HRB quotes (USD/t)  [Prices!G5:G7, S12]"
         / EU 686.090909090909, CN 395.090909090909, RW 465.181818181818 /
   S0(i) "base crude steel production (Mt)                                 [Trade!E13:E15, S1-S2]"
         / EU 126.0, CN 960.8, RW 762.6 /
   D0(i) "base consumption = column sums of net flow matrix (Mt)           [Trade!B16:D16]"
         / EU 149.2, CN 835.0, RW 865.2 /;

Table x0(i,j) "net bilateral base flows, origin i -> destination j (Mt)  [Trade!B13:D15]"
            EU          CN          RW
   EU    126.000       0.000       0.000
   CN      2.782     835.000     123.018
   RW     20.418       0.000     742.182 ;

*-------------------------------------------------------------------------------
* 2. Transport / trade-barrier costs, t(i,j) = |p0(j) - p0(i)|
*    Calibrated, not observed: includes EU safeguards and anti-dumping duties.
*-------------------------------------------------------------------------------
Table t(i,j) "unit transport + trade-barrier cost (USD/t)  [Calibration!B22:D24]"
                       EU                    CN                    RW
   EU                   0                   291      220.909090909091
   CN                 291                     0      70.0909090909091
   RW    220.909090909091      70.0909090909091                     0 ;

*-------------------------------------------------------------------------------
* 3. Emission intensities (Scope 1), e(i) = BOF share * e_BOF + EAF share * e_EAF
*    RoW shares are the residual of world minus EU minus China.
*-------------------------------------------------------------------------------
Scalars
   Sworld  "world crude steel production (Mt)          [Inputs!B5,  S1]" / 1849.4 /
   eBOF    "BF-BOF avg intensity, worldsteel (tCO2/t)  [Inputs!B16, S2]" / 2.34 /
   eEAF    "scrap-EAF avg intensity, worldsteel        [Inputs!B17, S2]" / 0.69 /
   eBOFalt "BF-BOF intensity, Clausing et al.          [Inputs!B19, S5]" / 2.40 /
   eEAFalt "EAF intensity, Clausing et al.             [Inputs!B20, S5]" / 0.90 /;

Parameter
   shBOF(*) "oxygen (BF-BOF) share of production  [Inputs!B9,B11,B13, S2]"
          / World 0.694, EU 0.542, CN 0.894 /
   shEAF(*) "electric (EAF) share of production   [Inputs!B10,B12,B14, S2]"
          / World 0.303, EU 0.458, CN 0.106 /;

* RoW shares as residual (Emissions!D5:D6)
shBOF('RW') = (Sworld*shBOF('World') - S0('EU')*shBOF('EU') - S0('CN')*shBOF('CN')) / S0('RW');
shEAF('RW') = (Sworld*shEAF('World') - S0('EU')*shEAF('EU') - S0('CN')*shEAF('CN')) / S0('RW');

Parameter
   e(i)      "emission intensity used by the model (tCO2/t)        [Emissions!B7:D7]"
   e_ws(i)   "intensity, worldsteel process averages (base)        [Emissions!B7:D7]"
   e_alt(i)  "intensity, Clausing et al. process averages (sens.)  [Emissions!B8:D8]"
   e_wb(i)   "workbook value of e_ws, for the assertion"
         / EU 1.5843, CN 2.1651, RW 1.45575111723053 /
   e_altwb(i) "workbook value of e_alt, for the assertion"
         / EU 1.713,  CN 2.241,  RW 1.59415318646735 /;

e_ws(i)  = shBOF(i)*eBOF    + shEAF(i)*eEAF;
e_alt(i) = shBOF(i)*eBOFalt + shEAF(i)*eEAFalt;
e(i)     = e_ws(i);

*-------------------------------------------------------------------------------
* 4. Carbon prices, CBAM phase factor, elasticities
*-------------------------------------------------------------------------------
Scalars
   cbamEUR "CBAM certificate price Q1-2026 (EUR/tCO2)       [Inputs!B42, S8]" / 75.36 /
   usdEUR  "EUR/USD 2025 average (USD per EUR)              [Inputs!B40, S7]" / 1.13 /
   tauEU   "EU carbon price (USD/tCO2) = cbamEUR*usdEUR     [Scenarios!B5]"
   tauCN   "China origin carbon price credited (USD/tCO2)   [Scenarios!B7]" / 0 /
   tauRW   "RoW origin carbon price (USD/tCO2)              [Scenarios!B8]" / 0 /
   tau0    "effective carbon cost in 2025 prices (USD/tCO2) [Inputs!B45, A3]" / 0 /
   phi     "CBAM phase-in factor, 2034 end state            [Scenarios!B23, S11]" / 1 /
   epsD    "steel demand price elasticity, base             [Inputs!B35, S5]" / -1.34 /
   epsDlo  "alt: low demand elasticity (Fernandez 2018)     [Inputs!B36, S6]" / -0.18 /
   etaS    "steel supply elasticity, base                   [Inputs!B37, S5]" / 0.01 /
   etaSmr  "alt: medium-run supply elasticity (ASSUMPTION)  [Inputs!B38, A2]" / 0.5 /;

tauEU = cbamEUR*usdEUR;

*-------------------------------------------------------------------------------
* 5. Calibrated parameters as computed in the workbook. NOT used by the model:
*    the model recalibrates in GAMS (model/cbam_calib.gms) and asserts that its
*    values match these.                                   [GAMS!B6:E8]
*-------------------------------------------------------------------------------
Table regwb(i,*) "workbook calibration: a, b (demand) and c, d (marginal cost)"
                        a                     b             c                     d
   EU    1198.09905020353      3.43168995383793        -67923      544.516594516595
   CN     689.93487109905     0.353106541327115        -39114      41.1210355007191
   RW    812.332428765265     0.401237413989189        -46053      60.9994516367451 ;

*-------------------------------------------------------------------------------
* 6. Data assertions
*-------------------------------------------------------------------------------
abort$(abs(tauEU - 85.1568) > 1e-9) "tauEU does not match Scenarios!B5", tauEU;
abort$(smax(i, abs(e_ws(i)  - e_wb(i)))    > 1e-9) "e(i) does not match Emissions!B7:D7", e_ws, e_wb;
abort$(smax(i, abs(e_alt(i) - e_altwb(i))) > 1e-9) "e_alt(i) does not match Emissions!B8:D8", e_alt, e_altwb;
abort$(smax(i, abs(sum(j, x0(i,j)) - S0(i))) > 1e-9) "x0 rows do not sum to S0", x0, S0;
abort$(smax(j, abs(sum(i, x0(i,j)) - D0(j))) > 1e-9) "x0 columns do not sum to D0", x0, D0;
abort$(smax((i,j), abs(t(i,j) - abs(p0(j) - p0(i)))) > 1e-9) "t(i,j) is not |p0(j)-p0(i)|", t, p0;
* Base year must be a competitive equilibrium: p0(i) + t(i,j) >= p0(j), = on active routes
abort$(smin((i,j), p0(i) + t(i,j) - p0(j)) < -1e-9) "base year: a route has negative margin", t;
abort$(smax((i,j)$(x0(i,j) > 0), abs(p0(i) + t(i,j) - p0(j))) > 1e-9) "base year: active route not at zero profit", t;

*-------------------------------------------------------------------------------
* 7. Model parameters, filled by model/cbam_calib.gms (CLAUDE.md Section 5.1)
*-------------------------------------------------------------------------------
Parameter
   a(j)   "demand intercept (USD/t)"
   b(j)   "demand slope ((USD/t)/Mt)"
   c(i)   "marginal-cost intercept (USD/t)"
   d(i)   "marginal-cost slope ((USD/t)/Mt)"
   MC0(i) "base-year marginal cost (USD/t)";
