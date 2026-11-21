$title CBAM x market power -- Cournot China, real 2025 data, Cournot calibration

$onText
===============================================================================
 cbam_cournot.gms -- Stage 2: strategic (Cournot) China
 Real 2025 data (data/cbam_data.gms), producer-based shipments x(i,j).

 Conduct: China is a Cournot player on its export routes (theta = 1 on CN->EU
 and CN->RW) and a price-taker at home (theta(CN,CN) = 0). EU and RW are
 price-takers everywhere.

 COURNOT CALIBRATION
 The 2025 base year is treated as a Cournot equilibrium. On China's export
 routes the observed price gap must equal transport cost PLUS the markup:
     p0(CN) + tc(CN,j) = p0(j) - theta(CN,j)*b(j)*x0(CN,j)
 so  tc(i,j) = t(i,j) - theta(i,j)*b(j)*x0(i,j)      (t = |p0(j) - p0(i)|)
 China's MC0 = p0(CN) is unchanged: it is a price-taker at home and serves it.
 Consequences:
   k_base reproduces p0, S0, D0 AND the base flows x0 exactly.
   c_base is the counterfactual "what if China were a price-taker".
 Requires tc >= 0. If a markup exceeds the observed price gap, the data are
 inconsistent with full Cournot conduct; that case is reported, not solved.

 Variables:  x(i,j) >= 0 shipments from producer i's output to market j (Mt)
             p(j) free price in market j (USD/t)

 Complementarity system (b > 0 in this sign convention):
   cost(i,j): c(i) + d(i)*S(i) + tauO(i)*e(i) + tc(i,j) + kappa(i,j)
              + sw_cournot*theta(i,j)*b(j)*x(i,j) - p(j) >= 0   _|_ x(i,j) >= 0
              (delivered MC >= perceived marginal revenue p(j) - theta*b*x)
   pdef(j):   p(j) - a(j) + b(j)*sum(i, x(i,j)) = 0             _|_ p(j) free

 Scenarios: c_base c_ets c_cbam k_base k_ets k_cbam.
 Leakage and welfare vs the no-policy run with the same conduct.
 Sensitivities: base, etaS_mr, epsD_lo, e_alt; recalibrated each.

 Output: results/cournot_results.csv (sens, scenario, variable, i, j, value)
         results/cournot_results.gdx
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
* Conduct and Cournot-calibrated transport costs
*-------------------------------------------------------------------------------
Table theta(i,j) "conduct: 1 = Cournot on route i->j   [GAMS!B21:D23]"
          EU    CN    RW
   EU      0     0     0
   CN      1     0     1
   RW      0     0     0 ;

Parameter tc(i,j) "Cournot-calibrated transport + barrier cost (USD/t)"
          mkup0(i,j) "base-year Cournot markup theta*b*x0 (USD/t)";

*-------------------------------------------------------------------------------
* Scenarios and sensitivity cases
*-------------------------------------------------------------------------------
Set sc "scenarios" / c_base "no policy, China price-taker"
                     c_ets  "EU carbon price, China price-taker"
                     c_cbam "EU carbon price + CBAM, China price-taker"
                     k_base "no policy, Cournot China (= 2025 base)"
                     k_ets  "EU carbon price, Cournot China"
                     k_cbam "EU carbon price + CBAM, Cournot China" /;
Alias (sc, sc2);
Table scen(sc,*) "policy and conduct switches"
            sw_ets   beta   cournot
   c_base     0       0       0
   c_ets      1       0       0
   c_cbam     1       1       0
   k_base     0       0       1
   k_ets      1       0       1
   k_cbam     1       1       1 ;

Set basemap(sc,sc2) "matching no-policy run (same conduct)"
       / (c_base, c_ets, c_cbam).c_base, (k_base, k_ets, k_cbam).k_base /
    compmap(sc,sc2) "competitive run of the same policy (Cournot warm start)"
       / k_base.c_base, k_ets.c_ets, k_cbam.c_cbam /;

Set sens "sensitivity cases (one at a time)"
       / base    "etaS 0.01, epsD -1.34, worldsteel e(i)"
         etaS_mr "supply elasticity 0.5 (author assumption A2)"
         epsD_lo "demand elasticity -0.18 (Fernandez 2018)"
         e_alt   "Clausing et al. emission intensities" /;

Scalars sw_ets "EU ETS switch" / 0 /, beta "CBAM switch" / 0 /,
        sw_cournot "conduct switch" / 0 /;
Parameters tauO(i)    "origin carbon price actually paid (USD/tCO2)"
           kappa(i,j) "CBAM charge on route i->j (USD/t)";

*-------------------------------------------------------------------------------
* MCP
*-------------------------------------------------------------------------------
Positive Variable x(i,j) "shipments from i's output to market j (Mt)";
Free Variable     p(j)   "price in market j (USD/t)";

Equations
   cost(i,j) "route i->j: delivered MC >= perceived MR   _|_ x(i,j) >= 0"
   pdef(j)   "inverse demand = market clearing           _|_ p(j) free";

* KKT stationarity of producer i in market j (Cournot term only where theta=1)
cost(i,j).. c(i) + d(i)*sum(k, x(i,k)) + tauO(i)*e(i) + tc(i,j) + kappa(i,j)
            + sw_cournot*theta(i,j)*b(j)*x(i,j) - p(j) =g= 0;

* Market clearing on the inverse demand curve, paired with p(j)
pdef(j)..   p(j) - a(j) + b(j)*sum(i, x(i,j)) =e= 0;

Model cbam_cournot / cost.x, pdef.p /;
option mcp = path, solprint = off, limrow = 0, limcol = 0;

*-------------------------------------------------------------------------------
* Storage
*-------------------------------------------------------------------------------
Parameters
   calOK(sens)       "1 = Cournot calibration feasible (all tc >= 0)"
   tcr(sens,i,j)     "Cournot-calibrated transport cost (USD/t)"
   mk0(sens,i,j)     "base-year Cournot markup (USD/t)"
   xr(sens,sc,i,j)   "shipments (Mt)"
   Sr(sens,sc,i)     "production (Mt)"
   Dr(sens,sc,j)     "consumption (Mt)"
   pr(sens,sc,j)     "price (USD/t)"
   mkr(sens,sc,i,j)  "Cournot markup theta*b*x (USD/t)"
   Er(sens,sc,i)     "emissions (MtCO2)"
   Ew(sens,sc)       "world emissions (MtCO2)"
   leak(sens,sc)     "leakage rate vs matching base"
   revCBAM(sens,sc)  "CBAM revenue (million USD)"
   revETS(sens,sc)   "EU carbon revenue (million USD)"
   CS(sens,sc,j)     "consumer surplus 0.5*b*D^2 (million USD)"
   prof(sens,sc,i)   "producer profit level, incl. markup rents (million USD)"
   dCS(sens,sc,j)    "change in CS vs matching base"
   dPS(sens,sc,i)    "change in PS vs matching base"
   dW(sens,sc,*)     "change in welfare CS+PS (+ revenue for EU) vs matching base"
   impInt(sens,sc)   "avg emission intensity of EU imports (tCO2/t)"
   resid(i,j)        "cost(i,j) residual at solution (USD/t)"
   maxres(sens,sc)   "max complementarity violation"
   iters(sens,sc)    "PATH iterations"
   calib(sens,*,i)   "k_base minus base-year data"
   calx(sens)        "max |x(k_base) - x0|"
   hand(sens,sc)     "CN->EU: [p(EU) - b*x] - [MC_CN + tc + kappa] (0 if route used)"
   headline(sens,*)  "leakage reduction from the CBAM, by conduct";

* filled by the competitive runs, read by the Cournot warm start
xr(sens,sc,i,j) = 0;  pr(sens,sc,j) = 0;

*-------------------------------------------------------------------------------
* Solve loop
*-------------------------------------------------------------------------------
loop(sens,
*  --- sensitivity settings, recalibrate demand/supply, then transport costs ---
   etaS = 0.01;  epsD = -1.34;  e(i) = e_ws(i);
   if(sameas(sens,'etaS_mr'), etaS = etaSmr;);
   if(sameas(sens,'epsD_lo'), epsD = epsDlo;);
   if(sameas(sens,'e_alt'),   e(i) = e_alt(i););
$include "%system.fp%cbam_calib.gms"
   mkup0(i,j) = theta(i,j)*b(j)*x0(i,j);
   tc(i,j)    = t(i,j) - mkup0(i,j);
   tcr(sens,i,j) = tc(i,j);
   mk0(sens,i,j) = mkup0(i,j);
   calOK(sens) = 1$(smin((i,j), tc(i,j)) >= 0);

   if(calOK(sens) = 0,
      display "COURNOT CALIBRATION INFEASIBLE: base markup exceeds observed price gap",
              sens, mkup0, t, tc;
   else
      loop(sc,
*        --- policy wedges ---
         sw_ets = scen(sc,'sw_ets');  beta = scen(sc,'beta');
         sw_cournot = scen(sc,'cournot');
         tauO('EU') = sw_ets*tauEU*phi;  tauO('CN') = tauCN;  tauO('RW') = tauRW;
         kappa(i,j) = 0;
         kappa(i,'EU')$(not sameas(i,'EU')) = beta*phi*max(tauEU*sw_ets - tauO(i), 0)*e(i);

*        --- start: competitive runs from the base point, Cournot runs from the
*            competitive solution of the same policy ---
         if(sw_cournot = 0,
            x.l(i,j) = x0(i,j);  p.l(j) = p0(j);
         else
            x.l(i,j) = sum(compmap(sc,sc2), xr(sens,sc2,i,j));
            p.l(j)   = sum(compmap(sc,sc2), pr(sens,sc2,j));
         );

         solve cbam_cournot using mcp;

*        --- status and complementarity residuals ---
         if(cbam_cournot.modelstat > 2 or cbam_cournot.solvestat <> 1,
            display "PATH FAILED", sens, sc, cbam_cournot.modelstat, cbam_cournot.solvestat;
            abort "MCP solve failed";
         );
         resid(i,j) = c(i) + d(i)*sum(k, x.l(i,k)) + tauO(i)*e(i) + tc(i,j) + kappa(i,j)
                      + sw_cournot*theta(i,j)*b(j)*x.l(i,j) - p.l(j);
         maxres(sens,sc) = max(smax((i,j), -resid(i,j)),
                               smax((i,j)$(x.l(i,j) > 1e-6), abs(resid(i,j))));
         if(maxres(sens,sc) > 1e-6,
            display "COMPLEMENTARITY VIOLATED", sens, sc, resid, x.l;
            abort "cost(i,j) residual check failed (see resid above)";
         );
         iters(sens,sc) = cbam_cournot.iterusd;

*        --- store ---
         xr(sens,sc,i,j)  = x.l(i,j);
         Sr(sens,sc,i)    = sum(j, x.l(i,j));
         Dr(sens,sc,j)    = sum(i, x.l(i,j));
         pr(sens,sc,j)    = p.l(j);
         mkr(sens,sc,i,j) = sw_cournot*theta(i,j)*b(j)*x.l(i,j);
         Er(sens,sc,i)    = e(i)*Sr(sens,sc,i);
         Ew(sens,sc)      = sum(i, Er(sens,sc,i));
         revCBAM(sens,sc) = sum(i, kappa(i,'EU')*x.l(i,'EU'));
         revETS(sens,sc)  = tauO('EU')*e('EU')*Sr(sens,sc,'EU');
         CS(sens,sc,j)    = 0.5*b(j)*sqr(Dr(sens,sc,j));
*        profit: includes the Cournot markup rent automatically
         prof(sens,sc,i)  = sum(j, (p.l(j) - tc(i,j) - kappa(i,j) - tauO(i)*e(i))*x.l(i,j))
                          - (c(i)*Sr(sens,sc,i) + 0.5*d(i)*sqr(Sr(sens,sc,i)));
         impInt(sens,sc)$sum(i$(not sameas(i,'EU')), x.l(i,'EU'))
            = sum(i$(not sameas(i,'EU')), e(i)*x.l(i,'EU'))
            / sum(i$(not sameas(i,'EU')), x.l(i,'EU'));

*        --- hand check on CN->EU, Cournot runs ---
*            p(EU) - b(EU)*x(CN,EU) = MC_CN + tc(CN,EU) + kappa(CN,EU)
         hand(sens,sc)$(sw_cournot and x.l('CN','EU') > 1e-6)
            = (p.l('EU') - b('EU')*x.l('CN','EU'))
            - (c('CN') + d('CN')*Sr(sens,sc,'CN') + tauO('CN')*e('CN')
               + tc('CN','EU') + kappa('CN','EU'));
         abort$(abs(hand(sens,sc)) > 1e-6) "hand check on CN->EU failed", hand;
      );

*     --- Cournot calibration test: k_base reproduces p0, S0, D0 and x0 ---
      calib(sens,'p',i) = pr(sens,'k_base',i) - p0(i);
      calib(sens,'S',i) = Sr(sens,'k_base',i) - S0(i);
      calib(sens,'D',i) = Dr(sens,'k_base',i) - D0(i);
      calx(sens) = smax((i,j), abs(xr(sens,'k_base',i,j) - x0(i,j)));
      abort$(smax(i, abs(calib(sens,'p',i)) + abs(calib(sens,'S',i)) + abs(calib(sens,'D',i))) > 1e-4)
           "Cournot calibration test failed (p, S, D)", calib;
      abort$(calx(sens) > 1e-4) "Cournot calibration test failed (flows)", calx;
   );
);

*-------------------------------------------------------------------------------
* Leakage, welfare, headline
*-------------------------------------------------------------------------------
dCS(sens,sc,j)$calOK(sens) = sum(basemap(sc,sc2), CS(sens,sc,j)   - CS(sens,sc2,j));
dPS(sens,sc,i)$calOK(sens) = sum(basemap(sc,sc2), prof(sens,sc,i) - prof(sens,sc2,i));
dW(sens,sc,i)  = dCS(sens,sc,i) + dPS(sens,sc,i);
dW(sens,sc,'EU')    = dW(sens,sc,'EU') + revETS(sens,sc) + revCBAM(sens,sc);
dW(sens,sc,'World') = sum(i, dW(sens,sc,i));

leak(sens,sc)$(calOK(sens) and not (sameas(sc,'c_base') or sameas(sc,'k_base')))
   = sum(basemap(sc,sc2),
         sum(i$(not sameas(i,'EU')), Er(sens,sc,i) - Er(sens,sc2,i))
       / (Er(sens,sc2,'EU') - Er(sens,sc,'EU')));

headline(sens,'competitive')$calOK(sens) = leak(sens,'c_ets') - leak(sens,'c_cbam');
headline(sens,'cournot')$calOK(sens)     = leak(sens,'k_ets') - leak(sens,'k_cbam');
headline(sens,'gap')$calOK(sens)         = headline(sens,'competitive') - headline(sens,'cournot');

*-------------------------------------------------------------------------------
* Display
*-------------------------------------------------------------------------------
option headline:4, leak:4, calOK:0, tcr:2, mk0:2, pr:2, Sr:3, Dr:3, xr:3:2:2,
       mkr:2:2:2, Er:2, Ew:2, impInt:4, revCBAM:1, revETS:1, dCS:1, dPS:1, dW:1,
       maxres:8, iters:0, calib:8, calx:8, hand:8;
display "=== HEADLINE: L(ets) - L(cbam), competitive vs Cournot China ===", headline;
display calOK, tcr, mk0, leak, pr, Sr, Dr, xr, mkr, Er, Ew, impInt, revCBAM, revETS,
        dCS, dPS, dW, maxres, iters, calib, calx, hand;

*-------------------------------------------------------------------------------
* Output
*-------------------------------------------------------------------------------
execute_unload "%RES%cournot_results.gdx", calOK, tcr, mk0, xr, Sr, Dr, pr, mkr, Er, Ew,
   leak, headline, revCBAM, revETS, CS, dCS, dPS, dW, impInt, maxres, calib, calx, hand;

File fout / "%RES%cournot_results.csv" /;
fout.pc = 5;  fout.nd = 8;  fout.nw = 0;  fout.nr = 2;
put fout 'sens', 'scenario', 'variable', 'i', 'j', 'value' /;
loop(sens$calOK(sens),
   loop(sc,
      loop((i,j), put sens.tl, sc.tl, 'x', i.tl, j.tl, xr(sens,sc,i,j) /;);
      loop((i,j)$theta(i,j), put sens.tl, sc.tl, 'markup', i.tl, j.tl, mkr(sens,sc,i,j) /;);
      loop(i, put sens.tl, sc.tl, 'S', i.tl, '-', Sr(sens,sc,i) /;);
      loop(j, put sens.tl, sc.tl, 'D', '-', j.tl, Dr(sens,sc,j) /;);
      loop(j, put sens.tl, sc.tl, 'p', '-', j.tl, pr(sens,sc,j) /;);
      loop(i, put sens.tl, sc.tl, 'E', i.tl, '-', Er(sens,sc,i) /;);
      put sens.tl, sc.tl, 'E_world', '-', '-', Ew(sens,sc) /;
      if(not (sameas(sc,'c_base') or sameas(sc,'k_base')),
         put sens.tl, sc.tl, 'leakage', '-', '-', leak(sens,sc) /;);
      put sens.tl, sc.tl, 'rev_CBAM', '-', 'EU', revCBAM(sens,sc) /;
      put sens.tl, sc.tl, 'rev_ETS', 'EU', '-', revETS(sens,sc) /;
      loop(j, put sens.tl, sc.tl, 'dCS', '-', j.tl, dCS(sens,sc,j) /;);
      loop(i, put sens.tl, sc.tl, 'dPS', i.tl, '-', dPS(sens,sc,i) /;);
      loop(i, put sens.tl, sc.tl, 'dW', i.tl, '-', dW(sens,sc,i) /;);
      put sens.tl, sc.tl, 'dW', 'World', '-', dW(sens,sc,'World') /;
      if(impInt(sens,sc) > 0,
         put sens.tl, sc.tl, 'import_intensity_EU', '-', 'EU', impInt(sens,sc) /;);
   );
   put sens.tl, 'all', 'leak_reduction_competitive', '-', '-', headline(sens,'competitive') /;
   put sens.tl, 'all', 'leak_reduction_cournot',     '-', '-', headline(sens,'cournot') /;
   put sens.tl, 'all', 'leak_reduction_gap',         '-', '-', headline(sens,'gap') /;
);
loop(sens$(calOK(sens) = 0), put sens.tl, 'all', 'cournot_calibration_infeasible', '-', '-', 1 /;);
putclose fout;
