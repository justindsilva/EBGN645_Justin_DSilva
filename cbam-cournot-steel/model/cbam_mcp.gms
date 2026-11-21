$title CBAM x market power -- competitive equilibrium as an MCP, scenarios, report

$onText
===============================================================================
 cbam_mcp.gms -- Competitive spatial price equilibrium posed directly as an MCP
 CLAUDE.md Sections 5.2, 5.3, 6 (steps 2-4, 7), 7 and 8. Competitive case only:
 every producer, China included, is a price-taker. (The Cournot term
 theta(i,j)*b(j)*x(i,j) of Section 5.3 is deliberately not included.)

 Variables:  x(i,j) >= 0 shipments (Mt),  p(j) free price (USD/t)
             production S(i) = sum_j x(i,j),  consumption D(j) = sum_i x(i,j)

 KKT / complementarity system:
   cost(i,j):  c(i) + d(i)*S(i) + tauO(i)*e(i) + t(i,j) + kappa(i,j) - p(j) >= 0
               _|_  x(i,j) >= 0
               (no-arbitrage: delivered marginal cost >= price, = if route used)
   pdef(j):    p(j) - a(j) + b(j)*D(j) = 0
               _|_  p(j) free
               (inverse demand = market clearing for a homogeneous good)

 Policy wedges (Section 5.2):
   tauO(EU) = sw_ets*tauEU*phi,  tauO(CN) = tauCN,  tauO(RW) = tauRW
   kappa(i,EU) = beta*phi*max(tauEU*sw_ets - tauO(i), 0)*e(i)   for i <> EU

 Outputs: results/results.csv and results/results.gdx (base data),
          results/sensitivity.csv (all sensitivity cases, Section 8).
===============================================================================
$offText

$setGlobal RES %system.fp%..%system.dirSep%results%system.dirSep%
$include "%system.fp%..%system.dirSep%data%system.dirSep%cbam_data.gms"
$include "%system.fp%cbam_calib.gms"

abort$(smax(j, abs(a(j) - regwb(j,'a'))) > 1e-6) "a(j) differs from workbook", a;
abort$(smax(j, abs(b(j) - regwb(j,'b'))) > 1e-6) "b(j) differs from workbook", b;
abort$(smax(i, abs(c(i) - regwb(i,'c'))) > 1e-6) "c(i) differs from workbook", c;
abort$(smax(i, abs(d(i) - regwb(i,'d'))) > 1e-6) "d(i) differs from workbook", d;

*-------------------------------------------------------------------------------
* Scenarios (Section 7, competitive rows) and sensitivity cases (Section 8)
*-------------------------------------------------------------------------------
Set sc "competitive scenarios" / c_base "no policy (calibration point)"
                                 c_ets  "EU carbon price, no CBAM"
                                 c_cbam "EU carbon price + CBAM" /;
Table scen(sc,*) "policy switches"
            sw_ets   beta
   c_base     0       0
   c_ets      1       0
   c_cbam     1       1 ;

Set sens "sensitivity cases (one at a time)"
       / base    "etaS 0.01, epsD -1.34, worldsteel e(i)"
         etaS_mr "supply elasticity 0.5 (author assumption A2)"
         epsD_lo "demand elasticity -0.18 (Fernandez 2018)"
         e_alt   "Clausing et al. emission intensities" /;

Scalars sw_ets "EU ETS switch" / 0 /
        beta   "CBAM switch"   / 0 /;
Parameters
   tauO(i)    "origin carbon price actually paid (USD/tCO2)"
   kappa(i,j) "CBAM charge per tonne on route i->j (USD/t)";

*-------------------------------------------------------------------------------
* MCP
*-------------------------------------------------------------------------------
Positive Variable x(i,j) "shipments origin i -> destination j (Mt)";
Free Variable     p(j)   "price in market j (USD/t)";

Equations
   cost(i,j) "zero-profit / no-arbitrage on route i->j   _|_ x(i,j) >= 0"
   pdef(j)   "inverse demand = market clearing in j     _|_ p(j) free";

* KKT stationarity of the competitive producer i in market j, paired with x(i,j)
cost(i,j).. c(i) + d(i)*sum(k, x(i,k)) + tauO(i)*e(i) + t(i,j) + kappa(i,j)
            - p(j) =g= 0;

* Market clearing: price on the inverse demand curve at D(j), paired with p(j)
pdef(j)..   p(j) - a(j) + b(j)*sum(i, x(i,j)) =e= 0;

Model cbam_mcp / cost.x, pdef.p /;
option mcp = path, solprint = off, limrow = 0, limcol = 0;

*-------------------------------------------------------------------------------
* Storage (indexed by sensitivity case and scenario)
*-------------------------------------------------------------------------------
Parameters
   xr(sens,sc,i,j)  "shipments (Mt)"
   Sr(sens,sc,i)    "production (Mt)"
   Dr(sens,sc,j)    "consumption (Mt)"
   pr(sens,sc,j)    "price (USD/t)"
   Er(sens,sc,i)    "emissions, Scope 1 (MtCO2)"
   Ew(sens,sc)      "world emissions (MtCO2)"
   leak(sens,sc)    "leakage rate vs matching base: [dE(CN)+dE(RW)] / [-dE(EU)]"
   revCBAM(sens,sc) "CBAM revenue (million USD)"
   revETS(sens,sc)  "EU carbon revenue from EU producers (million USD)"
   CS(sens,sc,j)    "consumer surplus 0.5*b*D^2 (million USD)"
   prof(sens,sc,i)  "producer profit LEVEL (meaningless alone: c << 0)"
   dPS(sens,sc,i)   "change in producer surplus vs matching base (million USD)"
   dCS(sens,sc,j)   "change in consumer surplus vs matching base (million USD)"
   dW(sens,sc,*)    "change in welfare CS+PS+revenue vs base (million USD; EU gets revenue)"
   impInt(sens,sc)  "avg emission intensity of EU imports (tCO2/t)"
   resid(i,j)       "residual of cost(i,j) at the solution (USD/t)"
   maxres(sens,sc)  "max complementarity violation over all routes"
   calib(sens,*,i)  "c_base minus base-year data (calibration test)"
   eUsed(sens,i)    "emission intensity used in each sensitivity case"
   abc(sens,*,i)    "calibrated a, b, c, d per sensitivity case";

*-------------------------------------------------------------------------------
* Solve loop
*-------------------------------------------------------------------------------
loop(sens,
*  --- sensitivity settings, then recalibrate (Section 8 -> Section 5.1) ---
   etaS = 0.01;  epsD = -1.34;  e(i) = e_ws(i);
   if(sameas(sens,'etaS_mr'), etaS = etaSmr;);
   if(sameas(sens,'epsD_lo'), epsD = epsDlo;);
   if(sameas(sens,'e_alt'),   e(i) = e_alt(i););
$include "%system.fp%cbam_calib.gms"
   eUsed(sens,i) = e(i);
   abc(sens,'a',i) = a(i);  abc(sens,'b',i) = b(i);
   abc(sens,'c',i) = c(i);  abc(sens,'d',i) = d(i);

   loop(sc,
*     --- policy wedges (Section 5.2) ---
      sw_ets = scen(sc,'sw_ets');  beta = scen(sc,'beta');
      tauO('EU') = sw_ets*tauEU*phi;  tauO('CN') = tauCN;  tauO('RW') = tauRW;
      kappa(i,j) = 0;
      kappa(i,'EU')$(not sameas(i,'EU')) = beta*phi*max(tauEU*sw_ets - tauO(i), 0)*e(i);

*     --- start every scenario at the base point ---
      x.l(i,j) = x0(i,j);  p.l(j) = p0(j);

      solve cbam_mcp using mcp;

*     --- solver status and complementarity residuals (Section 6, step 7) ---
      if(cbam_mcp.modelstat > 2 or cbam_mcp.solvestat <> 1,
         display "PATH FAILED", sens, sc, cbam_mcp.modelstat, cbam_mcp.solvestat;
         abort "MCP solve failed";
      );
      resid(i,j) = c(i) + d(i)*sum(k, x.l(i,k)) + tauO(i)*e(i) + t(i,j)
                   + kappa(i,j) - p.l(j);
      maxres(sens,sc) = max(smax((i,j), -resid(i,j)),
                            smax((i,j)$(x.l(i,j) > 1e-6), abs(resid(i,j))));
      if(maxres(sens,sc) > 1e-6,
         display "COMPLEMENTARITY VIOLATED", sens, sc, resid, x.l;
         abort "cost(i,j) residual check failed (see resid above)";
      );
      abort$(smax(j, abs(p.l(j) - a(j) + b(j)*sum(i, x.l(i,j)))) > 1e-6)
           "pdef(j) residual check failed", p.l;

*     --- store levels ---
      xr(sens,sc,i,j) = x.l(i,j);
      Sr(sens,sc,i)   = sum(j, x.l(i,j));
      Dr(sens,sc,j)   = sum(i, x.l(i,j));
      pr(sens,sc,j)   = p.l(j);
      Er(sens,sc,i)   = e(i)*Sr(sens,sc,i);
      Ew(sens,sc)     = sum(i, Er(sens,sc,i));
      revCBAM(sens,sc) = sum(i, kappa(i,'EU')*x.l(i,'EU'));
      revETS(sens,sc)  = tauO('EU')*e('EU')*Sr(sens,sc,'EU');
      CS(sens,sc,j)    = 0.5*b(j)*sqr(Dr(sens,sc,j));
      prof(sens,sc,i)  = sum(j, (p.l(j) - t(i,j) - kappa(i,j) - tauO(i)*e(i))*x.l(i,j))
                       - (c(i)*Sr(sens,sc,i) + 0.5*d(i)*sqr(Sr(sens,sc,i)));
      impInt(sens,sc)$sum(i$(not sameas(i,'EU')), x.l(i,'EU'))
         = sum(i$(not sameas(i,'EU')), e(i)*x.l(i,'EU'))
         / sum(i$(not sameas(i,'EU')), x.l(i,'EU'));
   );

*  --- calibration test (Section 6, step 3) for every recalibration ---
   calib(sens,'p',i) = pr(sens,'c_base',i) - p0(i);
   calib(sens,'S',i) = Sr(sens,'c_base',i) - S0(i);
   calib(sens,'D',i) = Dr(sens,'c_base',i) - D0(i);
   abort$(smax(i, abs(calib(sens,'p',i)) + abs(calib(sens,'S',i)) + abs(calib(sens,'D',i))) > 1e-4)
        "Calibration test failed", calib;
   abort$(abs(Ew(sens,'c_base') - sum(i, e(i)*S0(i))) > 1e-4)
        "Calibration test failed: emissions", Ew;
);

*-------------------------------------------------------------------------------
* Changes vs the matching (competitive) base, leakage, welfare (Section 7)
*-------------------------------------------------------------------------------
dPS(sens,sc,i) = prof(sens,sc,i) - prof(sens,'c_base',i);
dCS(sens,sc,j) = CS(sens,sc,j)   - CS(sens,'c_base',j);
dW(sens,sc,i)  = dCS(sens,sc,i) + dPS(sens,sc,i);
dW(sens,sc,'EU') = dW(sens,sc,'EU') + revETS(sens,sc) + revCBAM(sens,sc);
dW(sens,sc,'World') = sum(i, dW(sens,sc,i));

leak(sens,sc)$((not sameas(sc,'c_base')) and abs(Er(sens,sc,'EU') - Er(sens,'c_base','EU')) > 1e-9)
   = sum(i$(not sameas(i,'EU')), Er(sens,sc,i) - Er(sens,'c_base',i))
   / (Er(sens,'c_base','EU') - Er(sens,sc,'EU'));

Parameter headline(sens) "leakage reduction from adding the CBAM: L(c_ets) - L(c_cbam)";
headline(sens) = leak(sens,'c_ets') - leak(sens,'c_cbam');

*-------------------------------------------------------------------------------
* QP vs MCP agreement (Section 6, step 4) -- base data, all three scenarios
*-------------------------------------------------------------------------------
$ifThen exist "%RES%qp_results.gdx"
Parameters pQP(sc,j), SQP(sc,i), DQP(sc,j), qpdiff(sc,*,i) "QP minus MCP";
$gdxIn "%RES%qp_results.gdx"
$load pQP SQP DQP
$gdxIn
qpdiff(sc,'p',i) = pQP(sc,i) - pr('base',sc,i);
qpdiff(sc,'S',i) = SQP(sc,i) - Sr('base',sc,i);
qpdiff(sc,'D',i) = DQP(sc,i) - Dr('base',sc,i);
option qpdiff:8;
display "QP minus MCP (must be < 1e-4)", qpdiff;
abort$(smax((sc,i), abs(qpdiff(sc,'p',i)) + abs(qpdiff(sc,'S',i)) + abs(qpdiff(sc,'D',i))) > 1e-4)
     "QP and MCP disagree", qpdiff;
$else
display "NOTE: results/qp_results.gdx not found -- run model/cbam_qp.gms first for the QP-MCP check";
$endIf

*-------------------------------------------------------------------------------
* Display (base data first)
*-------------------------------------------------------------------------------
option headline:4, leak:4, pr:3, Sr:3, Dr:3, xr:3:2:2, Er:3, Ew:3, dW:1,
       dPS:1, dCS:1, revCBAM:1, revETS:1, impInt:4, maxres:8;
display "=== HEADLINE: leakage reduction from the CBAM, L(c_ets) - L(c_cbam) ===", headline;
display leak, pr, Sr, Dr, xr, Er, Ew, impInt, revCBAM, revETS, dCS, dPS, dW, maxres, calib, abc;

*-------------------------------------------------------------------------------
* Write results: results.gdx (everything), results.csv (base data, long format),
* sensitivity.csv (all cases, long format with a sens column)
*-------------------------------------------------------------------------------
execute_unload "%RES%results.gdx", xr, Sr, Dr, pr, Er, Ew, leak, headline, revCBAM,
   revETS, CS, dCS, dPS, dW, impInt, maxres, calib, eUsed, abc, kappa, tauO;

File fcsv / "%RES%results.csv" /;
fcsv.pc = 5;  fcsv.nd = 8;  fcsv.nw = 0;  fcsv.nr = 2;
put fcsv 'scenario', 'variable', 'i', 'j', 'value' /;
loop(sens$sameas(sens,'base'),
loop(sc,
   loop((i,j), put sc.tl, 'x', i.tl, j.tl, xr(sens,sc,i,j) /;);
   loop(i, put sc.tl, 'S', i.tl, '-', Sr(sens,sc,i) /;);
   loop(j, put sc.tl, 'D', '-', j.tl, Dr(sens,sc,j) /;);
   loop(j, put sc.tl, 'p', '-', j.tl, pr(sens,sc,j) /;);
   loop(i, put sc.tl, 'E', i.tl, '-', Er(sens,sc,i) /;);
   put sc.tl, 'E_world', '-', '-', Ew(sens,sc) /;
   if(not sameas(sc,'c_base'), put sc.tl, 'leakage', '-', '-', leak(sens,sc) /;);
   put sc.tl, 'rev_CBAM', '-', 'EU', revCBAM(sens,sc) /;
   put sc.tl, 'rev_ETS', 'EU', '-', revETS(sens,sc) /;
   loop(j, put sc.tl, 'CS', '-', j.tl, CS(sens,sc,j) /;);
   loop(j, put sc.tl, 'dCS', '-', j.tl, dCS(sens,sc,j) /;);
   loop(i, put sc.tl, 'dPS', i.tl, '-', dPS(sens,sc,i) /;);
   loop(i, put sc.tl, 'dW', i.tl, '-', dW(sens,sc,i) /;);
   put sc.tl, 'dW', 'World', '-', dW(sens,sc,'World') /;
   if(impInt(sens,sc) > 0, put sc.tl, 'import_intensity_EU', '-', 'EU', impInt(sens,sc) /;);
);
put 'c_ets-c_cbam', 'leakage_reduction', '-', '-', headline(sens) /;
);
putclose fcsv;

File fsen / "%RES%sensitivity.csv" /;
fsen.pc = 5;  fsen.nd = 8;  fsen.nw = 0;  fsen.nr = 2;
put fsen 'sens', 'scenario', 'variable', 'i', 'j', 'value' /;
loop((sens,sc),
   loop((i,j), put sens.tl, sc.tl, 'x', i.tl, j.tl, xr(sens,sc,i,j) /;);
   loop(i, put sens.tl, sc.tl, 'S', i.tl, '-', Sr(sens,sc,i) /;);
   loop(j, put sens.tl, sc.tl, 'D', '-', j.tl, Dr(sens,sc,j) /;);
   loop(j, put sens.tl, sc.tl, 'p', '-', j.tl, pr(sens,sc,j) /;);
   loop(i, put sens.tl, sc.tl, 'E', i.tl, '-', Er(sens,sc,i) /;);
   put sens.tl, sc.tl, 'E_world', '-', '-', Ew(sens,sc) /;
   if(not sameas(sc,'c_base'), put sens.tl, sc.tl, 'leakage', '-', '-', leak(sens,sc) /;);
   put sens.tl, sc.tl, 'rev_CBAM', '-', 'EU', revCBAM(sens,sc) /;
   put sens.tl, sc.tl, 'rev_ETS', 'EU', '-', revETS(sens,sc) /;
   loop(i, put sens.tl, sc.tl, 'dW', i.tl, '-', dW(sens,sc,i) /;);
   put sens.tl, sc.tl, 'dW', 'World', '-', dW(sens,sc,'World') /;
);
loop(sens, put sens.tl, 'c_ets-c_cbam', 'leakage_reduction', '-', '-', headline(sens) /;);
putclose fsen;
