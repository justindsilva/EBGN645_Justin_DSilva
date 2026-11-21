$title CBAM x market power -- non-prohibitive CBAM via the 2026-2034 phase-in

$onText
===============================================================================
 cbam_phasein.gms -- Cournot China, real 2025 data, CBAM phase-in path
-------------------------------------------------------------------------------
 Why: at full strength (phi = 1, the 2034 end state) the CBAM is prohibitive.
 The charge on RoW steel, tauEU*e(RW) = 124 USD/t, lifts RoW's delivered cost
 to the EU (~806 USD/t) above the EU's own self-sufficient price (~766 USD/t),
 so EU imports stop, CBAM revenue is zero and no shuffling can be observed.

 Change: phi, the phase-in factor, is looped over the official schedule
 (2026: 2.5% ... 2034: 100%; Scenarios!B12:B20, S11). As in cbam_mcp.gms and
 cbam_cournot.gms, phi scales BOTH the CBAM charge and the EU producers' own
 carbon cost (the CBAM phases in as EU free allocation phases out):
     tauO(EU)    = sw_ets*tauEU*phi
     kappa(i,EU) = beta*phi*max(tauEU*sw_ets - tauO(i), 0)*e(i),   i <> EU
 For low phi the charge is below the import-parity gap, imports continue, and
 the CBAM reallocates them between origins instead of shutting them off.
 Every year is flagged prohibitive or not from the solution itself.

 Everything else is identical to cbam_cournot.gms: producer-based shipments,
 Cournot China on its export routes (theta = 1 on CN->EU, CN->RW), Cournot
 calibration tc(i,j) = t(i,j) - theta(i,j)*b(j)*x0(i,j), so k_base = 2025.

 Complementarity system (b > 0 in this sign convention):
   cost(i,j): c(i) + d(i)*S(i) + tauO(i)*e(i) + tc(i,j) + kappa(i,j)
              + sw_cournot*theta(i,j)*b(j)*x(i,j) - p(j) >= 0   _|_ x(i,j) >= 0
   pdef(j):   p(j) - a(j) + b(j)*sum(i, x(i,j)) = 0             _|_ p(j) free

 Output: results/phasein_results.csv (sens, year, scenario, variable, i, j, value)
         results/phasein_results.gdx
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
* CBAM phase-in schedule (the new dimension)
*-------------------------------------------------------------------------------
Set yr "CBAM phase-in year" / 2026*2034 /;
Parameter phiyr(yr) "share of embedded emissions charged   [Scenarios!B12:B20, S11]"
   / 2026 0.025, 2027 0.050, 2028 0.100, 2029 0.225, 2030 0.485,
     2031 0.610, 2032 0.735, 2033 0.860, 2034 1.000 /;

*-------------------------------------------------------------------------------
* Conduct and Cournot-calibrated transport costs
*-------------------------------------------------------------------------------
Table theta(i,j) "conduct: 1 = Cournot on route i->j   [GAMS!B21:D23]"
          EU    CN    RW
   EU      0     0     0
   CN      1     0     1
   RW      0     0     0 ;

Parameter tc(i,j)    "Cournot-calibrated transport + barrier cost (USD/t)"
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
       / k_base.c_base, k_ets.c_ets, k_cbam.c_cbam /
    cbamsc(sc)      "CBAM scenarios" / c_cbam, k_cbam /
    etsmap(sc,sc2)  "CBAM run -> ETS-only run, same conduct"
       / c_cbam.c_ets, k_cbam.k_ets /;

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
* MCP (unchanged from cbam_cournot.gms)
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

Model cbam_phasein / cost.x, pdef.p /;
option mcp = path, solprint = off, limrow = 0, limcol = 0;

*-------------------------------------------------------------------------------
* Storage (indexed by sensitivity case, phase-in year, scenario)
*-------------------------------------------------------------------------------
Parameters
   calOK(sens)          "1 = Cournot calibration feasible (all tc >= 0)"
   kap(sens,yr,i)       "CBAM charge on i->EU in the CBAM runs (USD/t)"
   xr(sens,yr,sc,i,j)   "shipments (Mt)"
   Sr(sens,yr,sc,i)     "production (Mt)"
   Dr(sens,yr,sc,j)     "consumption (Mt)"
   pr(sens,yr,sc,j)     "price (USD/t)"
   mkr(sens,yr,sc,i,j)  "Cournot markup theta*b*x (USD/t)"
   Er(sens,yr,sc,i)     "emissions (MtCO2)"
   Ew(sens,yr,sc)       "world emissions (MtCO2)"
   imp(sens,yr,sc)      "EU imports (Mt)"
   leak(sens,yr,sc)     "leakage rate vs matching base"
   revCBAM(sens,yr,sc)  "CBAM revenue (million USD)"
   revETS(sens,yr,sc)   "EU carbon revenue (million USD)"
   CS(sens,yr,sc,j)     "consumer surplus 0.5*b*D^2 (million USD)"
   prof(sens,yr,sc,i)   "producer profit level, incl. markup rents (million USD)"
   dW(sens,yr,sc,*)     "change in welfare CS+PS (+ revenue for EU) vs matching base"
   impInt(sens,yr,sc)   "avg emission intensity of EU imports (tCO2/t)"
   prohib(sens,yr,sc)   "1 = CBAM shuts off all EU imports"
   resid(i,j)           "cost(i,j) residual at solution (USD/t)"
   maxres(sens,yr,sc)   "max complementarity violation"
   calx(sens)           "max |x(k_base) - x0| (Cournot calibration test)"
   hand(sens,yr,sc)     "CN->EU: [p(EU) - b*x] - [MC_CN + tc + kappa] (0 if route used)"
   headline(sens,yr,*)  "CBAM effect, competitive vs Cournot China";

* filled by the competitive runs, read by the Cournot warm start
xr(sens,yr,sc,i,j) = 0;  pr(sens,yr,sc,j) = 0;

*-------------------------------------------------------------------------------
* Solve loop: sensitivity x phase-in year x scenario
*-------------------------------------------------------------------------------
loop(sens,
   etaS = 0.01;  epsD = -1.34;  e(i) = e_ws(i);
   if(sameas(sens,'etaS_mr'), etaS = etaSmr;);
   if(sameas(sens,'epsD_lo'), epsD = epsDlo;);
   if(sameas(sens,'e_alt'),   e(i) = e_alt(i););
$include "%system.fp%cbam_calib.gms"
   mkup0(i,j) = theta(i,j)*b(j)*x0(i,j);
   tc(i,j)    = t(i,j) - mkup0(i,j);
   calOK(sens) = 1$(smin((i,j), tc(i,j)) >= 0);

   if(calOK(sens) = 0,
      display "COURNOT CALIBRATION INFEASIBLE: base markup exceeds observed price gap",
              sens, mkup0, t, tc;
   else
      loop(yr,
         phi = phiyr(yr);
         loop(sc,
*           --- policy wedges (Section 5.2), scaled by this year's phi ---
            sw_ets = scen(sc,'sw_ets');  beta = scen(sc,'beta');
            sw_cournot = scen(sc,'cournot');
            tauO('EU') = sw_ets*tauEU*phi;  tauO('CN') = tauCN;  tauO('RW') = tauRW;
            kappa(i,j) = 0;
            kappa(i,'EU')$(not sameas(i,'EU')) = beta*phi*max(tauEU*sw_ets - tauO(i), 0)*e(i);
            kap(sens,yr,i)$beta = kappa(i,'EU');

*           --- warm start: competitive from base, Cournot from competitive ---
            if(sw_cournot = 0,
               x.l(i,j) = x0(i,j);  p.l(j) = p0(j);
            else
               x.l(i,j) = sum(compmap(sc,sc2), xr(sens,yr,sc2,i,j));
               p.l(j)   = sum(compmap(sc,sc2), pr(sens,yr,sc2,j));
            );

            solve cbam_phasein using mcp;

            if(cbam_phasein.modelstat > 2 or cbam_phasein.solvestat <> 1,
               display "PATH FAILED", sens, yr, sc, cbam_phasein.modelstat, cbam_phasein.solvestat;
               abort "MCP solve failed";
            );
            resid(i,j) = c(i) + d(i)*sum(k, x.l(i,k)) + tauO(i)*e(i) + tc(i,j) + kappa(i,j)
                         + sw_cournot*theta(i,j)*b(j)*x.l(i,j) - p.l(j);
            maxres(sens,yr,sc) = max(smax((i,j), -resid(i,j)),
                                     smax((i,j)$(x.l(i,j) > 1e-6), abs(resid(i,j))));
            if(maxres(sens,yr,sc) > 1e-6,
               display "COMPLEMENTARITY VIOLATED", sens, yr, sc, resid, x.l;
               abort "cost(i,j) residual check failed (see resid above)";
            );

*           --- store ---
            xr(sens,yr,sc,i,j)  = x.l(i,j);
            Sr(sens,yr,sc,i)    = sum(j, x.l(i,j));
            Dr(sens,yr,sc,j)    = sum(i, x.l(i,j));
            pr(sens,yr,sc,j)    = p.l(j);
            mkr(sens,yr,sc,i,j) = sw_cournot*theta(i,j)*b(j)*x.l(i,j);
            Er(sens,yr,sc,i)    = e(i)*Sr(sens,yr,sc,i);
            Ew(sens,yr,sc)      = sum(i, Er(sens,yr,sc,i));
            imp(sens,yr,sc)     = sum(i$(not sameas(i,'EU')), x.l(i,'EU'));
            revCBAM(sens,yr,sc) = sum(i, kappa(i,'EU')*x.l(i,'EU'));
            revETS(sens,yr,sc)  = tauO('EU')*e('EU')*Sr(sens,yr,sc,'EU');
            CS(sens,yr,sc,j)    = 0.5*b(j)*sqr(Dr(sens,yr,sc,j));
            prof(sens,yr,sc,i)  = sum(j, (p.l(j) - tc(i,j) - kappa(i,j) - tauO(i)*e(i))*x.l(i,j))
                                - (c(i)*Sr(sens,yr,sc,i) + 0.5*d(i)*sqr(Sr(sens,yr,sc,i)));
            impInt(sens,yr,sc)$(imp(sens,yr,sc) > 1e-6)
               = sum(i$(not sameas(i,'EU')), e(i)*x.l(i,'EU')) / imp(sens,yr,sc);

*           --- hand check (Section 6, step 8) on CN->EU, Cournot runs ---
            hand(sens,yr,sc)$(sw_cournot and x.l('CN','EU') > 1e-6)
               = (p.l('EU') - b('EU')*x.l('CN','EU'))
               - (c('CN') + d('CN')*Sr(sens,yr,sc,'CN') + tauO('CN')*e('CN')
                  + tc('CN','EU') + kappa('CN','EU'));
            abort$(abs(hand(sens,yr,sc)) > 1e-6) "hand check on CN->EU failed", hand;
         );
      );

*     --- Cournot calibration test (no-policy run does not depend on phi) ---
      calx(sens) = smax((yr,i,j), abs(xr(sens,yr,'k_base',i,j) - x0(i,j)));
      abort$(calx(sens) > 1e-4) "Cournot calibration test failed (flows)", calx;
      abort$(smax((yr,i), abs(pr(sens,yr,'k_base',i) - p0(i))) > 1e-4)
           "Cournot calibration test failed (prices)", pr;
   );
);

*-------------------------------------------------------------------------------
* Leakage, welfare, prohibitive flag, headline
*-------------------------------------------------------------------------------
leak(sens,yr,sc)$(calOK(sens) and not (sameas(sc,'c_base') or sameas(sc,'k_base')))
   = sum(basemap(sc,sc2),
         sum(i$(not sameas(i,'EU')), Er(sens,yr,sc,i) - Er(sens,yr,sc2,i))
       / (Er(sens,yr,sc2,'EU') - Er(sens,yr,sc,'EU')));

dW(sens,yr,sc,i)$calOK(sens) = sum(basemap(sc,sc2), CS(sens,yr,sc,i) + prof(sens,yr,sc,i)
                                                  - CS(sens,yr,sc2,i) - prof(sens,yr,sc2,i));
dW(sens,yr,sc,'EU')    = dW(sens,yr,sc,'EU') + revETS(sens,yr,sc) + revCBAM(sens,yr,sc);
dW(sens,yr,sc,'World') = sum(i, dW(sens,yr,sc,i));

prohib(sens,yr,cbamsc)$calOK(sens) = 1$(imp(sens,yr,cbamsc) < 1e-6);

* leakage reduction from the CBAM, by conduct (Section 7 headline)
headline(sens,yr,'dLeak_c')$calOK(sens) = leak(sens,yr,'c_ets') - leak(sens,yr,'c_cbam');
headline(sens,yr,'dLeak_k')$calOK(sens) = leak(sens,yr,'k_ets') - leak(sens,yr,'k_cbam');
headline(sens,yr,'dLeak_gap')           = headline(sens,yr,'dLeak_c') - headline(sens,yr,'dLeak_k');
* extra world emission cut from adding the CBAM on top of the ETS (MtCO2)
headline(sens,yr,'dEw_c')$calOK(sens)   = Ew(sens,yr,'c_cbam') - Ew(sens,yr,'c_ets');
headline(sens,yr,'dEw_k')$calOK(sens)   = Ew(sens,yr,'k_cbam') - Ew(sens,yr,'k_ets');
* resource shuffling: change in EU import intensity and in China->RoW tonnage
headline(sens,yr,'dImpInt_c')$(calOK(sens) and imp(sens,yr,'c_cbam') > 1e-6)
   = impInt(sens,yr,'c_cbam') - impInt(sens,yr,'c_ets');
headline(sens,yr,'dImpInt_k')$(calOK(sens) and imp(sens,yr,'k_cbam') > 1e-6)
   = impInt(sens,yr,'k_cbam') - impInt(sens,yr,'k_ets');
headline(sens,yr,'dCNRW_c')$calOK(sens) = xr(sens,yr,'c_cbam','CN','RW') - xr(sens,yr,'c_ets','CN','RW');
headline(sens,yr,'dCNRW_k')$calOK(sens) = xr(sens,yr,'k_cbam','CN','RW') - xr(sens,yr,'k_ets','CN','RW');
headline(sens,yr,'prohib_c') = prohib(sens,yr,'c_cbam');
headline(sens,yr,'prohib_k') = prohib(sens,yr,'k_cbam');

*-------------------------------------------------------------------------------
* Display
*-------------------------------------------------------------------------------
option headline:3:2:1, kap:2:2:1, leak:4:2:1, prohib:0:2:1, imp:3:2:1, impInt:4:2:1,
       revCBAM:1:2:1, pr:2:3:1, xr:3:3:2, mkr:2:3:2, Ew:2:2:1, dW:1:3:1,
       maxres:8:2:1, calx:8, hand:8:2:1;
display "=== CBAM by phase-in year: competitive (_c) vs Cournot China (_k) ===", headline;
display calOK, kap, imp, prohib, leak, impInt, revCBAM, pr, xr, mkr, Ew, dW, maxres, calx, hand;

*-------------------------------------------------------------------------------
* Output
*-------------------------------------------------------------------------------
execute_unload "%RES%phasein_results.gdx", phiyr, calOK, kap, xr, Sr, Dr, pr, mkr, Er, Ew,
   imp, leak, revCBAM, revETS, CS, dW, impInt, prohib, headline, maxres, calx, hand;

File fout / "%RES%phasein_results.csv" /;
fout.pc = 5;  fout.nd = 8;  fout.nw = 0;  fout.nr = 2;
put fout 'sens', 'year', 'scenario', 'variable', 'i', 'j', 'value' /;
loop((sens,yr)$calOK(sens),
   put sens.tl, yr.tl, 'all', 'phi', '-', '-', phiyr(yr) /;
   loop(sc,
      loop((i,j), put sens.tl, yr.tl, sc.tl, 'x', i.tl, j.tl, xr(sens,yr,sc,i,j) /;);
      loop((i,j)$theta(i,j), put sens.tl, yr.tl, sc.tl, 'markup', i.tl, j.tl, mkr(sens,yr,sc,i,j) /;);
      loop(i, put sens.tl, yr.tl, sc.tl, 'S', i.tl, '-', Sr(sens,yr,sc,i) /;);
      loop(j, put sens.tl, yr.tl, sc.tl, 'D', '-', j.tl, Dr(sens,yr,sc,j) /;);
      loop(j, put sens.tl, yr.tl, sc.tl, 'p', '-', j.tl, pr(sens,yr,sc,j) /;);
      loop(i, put sens.tl, yr.tl, sc.tl, 'E', i.tl, '-', Er(sens,yr,sc,i) /;);
      put sens.tl, yr.tl, sc.tl, 'E_world', '-', '-', Ew(sens,yr,sc) /;
      put sens.tl, yr.tl, sc.tl, 'imports_EU', '-', 'EU', imp(sens,yr,sc) /;
      if(not (sameas(sc,'c_base') or sameas(sc,'k_base')),
         put sens.tl, yr.tl, sc.tl, 'leakage', '-', '-', leak(sens,yr,sc) /;);
      put sens.tl, yr.tl, sc.tl, 'rev_CBAM', '-', 'EU', revCBAM(sens,yr,sc) /;
      put sens.tl, yr.tl, sc.tl, 'rev_ETS', 'EU', '-', revETS(sens,yr,sc) /;
      loop(i, put sens.tl, yr.tl, sc.tl, 'dW', i.tl, '-', dW(sens,yr,sc,i) /;);
      put sens.tl, yr.tl, sc.tl, 'dW', 'World', '-', dW(sens,yr,sc,'World') /;
      if(impInt(sens,yr,sc) > 0,
         put sens.tl, yr.tl, sc.tl, 'import_intensity_EU', '-', 'EU', impInt(sens,yr,sc) /;);
   );
   loop(i$(not sameas(i,'EU')), put sens.tl, yr.tl, 'cbam', 'kappa', i.tl, 'EU', kap(sens,yr,i) /;);
   put sens.tl, yr.tl, 'all', 'leak_reduction_competitive', '-', '-', headline(sens,yr,'dLeak_c') /;
   put sens.tl, yr.tl, 'all', 'leak_reduction_cournot',     '-', '-', headline(sens,yr,'dLeak_k') /;
   put sens.tl, yr.tl, 'all', 'leak_reduction_gap',         '-', '-', headline(sens,yr,'dLeak_gap') /;
   put sens.tl, yr.tl, 'all', 'cbam_dEworld_competitive',   '-', '-', headline(sens,yr,'dEw_c') /;
   put sens.tl, yr.tl, 'all', 'cbam_dEworld_cournot',       '-', '-', headline(sens,yr,'dEw_k') /;
   put sens.tl, yr.tl, 'all', 'prohibitive_competitive',    '-', '-', headline(sens,yr,'prohib_c') /;
   put sens.tl, yr.tl, 'all', 'prohibitive_cournot',        '-', '-', headline(sens,yr,'prohib_k') /;
);
loop(sens$(calOK(sens) = 0), put sens.tl, 'all', 'all', 'cournot_calibration_infeasible', '-', '-', 1 /;);
putclose fout;
