$title CBAM x market power -- Stage 1: competitive equilibrium as a QP

$onText
===============================================================================
 cbam_qp.gms -- Competitive spatial price equilibrium as net-social-payoff QP
 (Samuelson 1952; Takayama & Judge 1971).

   max NSP = sum_j [a(j)*D(j) - b(j)/2*D(j)^2]                area under demand
           - sum_i [c(i)*S(i) + d(i)/2*S(i)^2]                area under MC
           - sum_ij [t(i,j) + tauO(i)*e(i) + kappa(i,j)]*x(i,j) transport + wedges
   s.t.  D(j) = sum_i x(i,j)      multiplier = p(j)   (market clearing)
         S(i) = sum_j x(i,j)      multiplier = MC(i)  (production balance)
         x, D, S >= 0

 KKT (stationarity in x(i,j), paired with x(i,j) >= 0):
   0 <= x(i,j)  _|_  MC(i) + tauO(i)*e(i) + t(i,j) + kappa(i,j) - p(j) >= 0
 which is exactly the cost(i,j) condition of the MCP in cbam_mcp.gms.

 Carbon and CBAM terms are computational wedges, not resource costs; NSP is
 not welfare. Welfare is computed in cbam_mcp.gms with the revenue added back.

 Output: results/qp_results.gdx (p, S, D, x by scenario), read by cbam_mcp.gms
 to check QP-MCP agreement.
===============================================================================
$offText

$include "%system.fp%..%system.dirSep%data%system.dirSep%cbam_data.gms"
$include "%system.fp%cbam_calib.gms"

* Calibration in GAMS must reproduce the workbook's a, b, c, d
abort$(smax(j, abs(a(j) - regwb(j,'a'))) > 1e-6) "a(j) differs from workbook", a;
abort$(smax(j, abs(b(j) - regwb(j,'b'))) > 1e-6) "b(j) differs from workbook", b;
abort$(smax(i, abs(c(i) - regwb(i,'c'))) > 1e-6) "c(i) differs from workbook", c;
abort$(smax(i, abs(d(i) - regwb(i,'d'))) > 1e-6) "d(i) differs from workbook", d;

*-------------------------------------------------------------------------------
* Scenarios (competitive only)
*-------------------------------------------------------------------------------
Set sc "competitive scenarios" / c_base "no policy (calibration point)"
                                c_ets  "EU carbon price, no CBAM"
                                c_cbam "EU carbon price + CBAM" /;
Table scen(sc,*) "policy switches"
            sw_ets   beta
   c_base     0       0
   c_ets      1       0
   c_cbam     1       1 ;

Scalars sw_ets "EU ETS switch" / 0 /
        beta   "CBAM switch"   / 0 /;
Parameters
   tauO(i)    "origin carbon price actually paid (USD/tCO2)"
   kappa(i,j) "CBAM charge per tonne on route i->j (USD/t)";

*-------------------------------------------------------------------------------
* Model
*-------------------------------------------------------------------------------
Variable          NSP     "net social payoff (million USD; includes wedges)";
* (GAMS is case-insensitive: D/S would clash with parameter d and set s)
Positive Variable x(i,j)  "shipments origin i -> destination j (Mt)"
                  Dem(j)  "consumption D(j) in j (Mt)"
                  Sup(i)  "production S(i) in i (Mt)";

Equations
   obj     "net social payoff"
   dbal(j) "market clearing in j        -- multiplier is p(j)"
   sbal(i) "production balance in i     -- multiplier is MC(i)";

obj..     NSP =e= sum(j, a(j)*Dem(j) - 0.5*b(j)*sqr(Dem(j)))
                - sum(i, c(i)*Sup(i) + 0.5*d(i)*sqr(Sup(i)))
                - sum((i,j), (t(i,j) + tauO(i)*e(i) + kappa(i,j))*x(i,j));

* Written so that relaxing the RHS raises NSP: marginal = +p(j) and +MC(i)
dbal(j).. Dem(j)           =e= sum(i, x(i,j));
sbal(i).. sum(j, x(i,j)) =e= Sup(i);

Model cbam_qp / obj, dbal, sbal /;
option nlp = conopt, solprint = off, limrow = 0, limcol = 0;
cbam_qp.optfile = 1;
$onEcho > conopt.opt
* tight reduced-gradient tolerance so marginals (prices) are accurate to << 1e-4
Tol_Optimality 1.e-12
$offEcho

Parameters
   pQP(sc,j)    "price from dbal marginal (USD/t)"
   MCQP(sc,i)   "marginal cost from sbal marginal (USD/t)"
   SQP(sc,i)    "production (Mt)"
   DQP(sc,j)    "consumption (Mt)"
   xQP(sc,i,j)  "shipments (Mt)"
   pchk(sc,j)   "inverse demand a - b*D at QP solution (USD/t)"
   NSPQP(sc)    "objective value";

loop(sc,
*  --- policy wedges ---
   sw_ets = scen(sc,'sw_ets');  beta = scen(sc,'beta');
   tauO('EU') = sw_ets*tauEU*phi;  tauO('CN') = tauCN;  tauO('RW') = tauRW;
   kappa(i,j) = 0;
   kappa(i,'EU')$(not sameas(i,'EU')) = beta*phi*max(tauEU*sw_ets - tauO(i), 0)*e(i);

*  --- start at the base point ---
   x.l(i,j) = x0(i,j);  Dem.l(j) = D0(j);  Sup.l(i) = S0(i);

   solve cbam_qp using nlp maximizing NSP;
   abort$(cbam_qp.modelstat > 2 or cbam_qp.solvestat <> 1)
        "QP solve failed", cbam_qp.modelstat, cbam_qp.solvestat;

   pQP(sc,j)   = dbal.m(j);
   MCQP(sc,i)  = sbal.m(i);
   SQP(sc,i)   = Sup.l(i);
   DQP(sc,j)   = Dem.l(j);
   xQP(sc,i,j) = x.l(i,j);
   pchk(sc,j)  = a(j) - b(j)*Dem.l(j);
   NSPQP(sc)   = NSP.l;

*  The market-clearing marginal must equal the demand price, MC marginal = c + d*S
   abort$(smax(j, abs(pQP(sc,j) - pchk(sc,j))) > 1e-4)
        "QP: dbal marginal differs from a - b*D", pQP, pchk;
   abort$(smax(i, abs(MCQP(sc,i) - c(i) - d(i)*Sup.l(i))) > 1e-4)
        "QP: sbal marginal differs from c + d*S", MCQP;
);

*-------------------------------------------------------------------------------
* Calibration test: no policy reproduces the base
*-------------------------------------------------------------------------------
Parameter calQP(*,i) "c_base minus base-year data";
calQP('p',i) = pQP('c_base',i) - p0(i);
calQP('S',i) = SQP('c_base',i) - S0(i);
calQP('D',i) = DQP('c_base',i) - D0(i);
calQP('E',i) = e(i)*(SQP('c_base',i) - S0(i));
display calQP;
abort$(smax((i), abs(calQP('p',i))) > 1e-4) "Calibration test failed: p", calQP;
abort$(smax((i), abs(calQP('S',i))) > 1e-4) "Calibration test failed: S", calQP;
abort$(smax((i), abs(calQP('D',i))) > 1e-4) "Calibration test failed: D", calQP;
abort$(abs(sum(i, calQP('E',i))) > 1e-4)    "Calibration test failed: emissions", calQP;

option pQP:4, SQP:4, DQP:4, xQP:4:1:2;
display "QP competitive solutions", pQP, SQP, DQP, xQP, NSPQP;

execute_unload "%system.fp%..%system.dirSep%results%system.dirSep%qp_results.gdx",
   pQP, MCQP, SQP, DQP, xQP, NSPQP;
