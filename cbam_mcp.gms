* CBAM toy model -- 3-region steel market, competitive case


set r /EU, CN, RW/ ;
alias(r,rr) ;

parameter
pbar(r)    "base price (USD/t)"              /EU 600, CN 400, RW 450/,
qbar_s(r)  "base production (Mt)"            /EU 120, CN 1000, RW 780/,
qbar_d(r)  "base consumption (Mt)"           /EU 150, CN 990, RW 760/,
em(r)      "emission intensity (tCO2/t)"     /EU 1.4, CN 2.2, RW 0.8/,
e_s(r)     "supply elasticity"               /EU 1.0, CN 0.5, RW 1.0/,
e_d        "demand elasticity"               /-0.3/ ;

* unit transport + trade-barrier cost, origin r -> destination rr (USD/t)
* chosen so used base routes have zero margin: pbar(r) + t(r,rr) = pbar(rr)
table t(r,rr)
        EU    CN    RW
   EU    0   250   200
   CN  200     0    52
   RW  150    98     0 ;

scalar
tauEU  "EU carbon price (USD/tCO2)"          /50/,
tauCN  "China carbon price credited"         /0/,
tauRW  "RoW carbon price credited"           /0/ ;

*!!! calibration: linear demand P = a + b*Qd (b < 0), linear MC = c + d*Qs
parameter a(r), b(r), c(r), d(r) ;

b(r) = pbar(r) / (e_d * qbar_d(r)) ;
a(r) = pbar(r) - b(r) * qbar_d(r) ;
d(r) = pbar(r) / (e_s(r) * qbar_s(r)) ;
c(r) = pbar(r) - d(r) * qbar_s(r) ;

*!!! policy wedges, set per scenario inside the loop
set scen /c_base "no policy", c_ets "EU carbon price only", c_cbam "EU carbon price + CBAM"/ ;

table swtab(scen,*)
         ets   cbam
   c_base  0      0
   c_ets   1      0
   c_cbam  1      1 ;

parameter
tauO(r)     "origin carbon price actually paid (USD/tCO2)",
kappa(r,rr) "CBAM charge on route r -> rr (USD/t)" ;

*!!!!! begin QCP: maximize net social payoff
* carbon and CBAM terms are wedges, not resource costs: revenue is a transfer

positive variable Qd(r), Qs(r) ;
positive variable X(r,rr) "exports r -> rr (r <> rr)" ;
variable W "net social payoff (incl. wedges)" ;

equation objfn, market_clearing(r) ;

objfn.. W =e=
    sum(r, a(r) * Qd(r) + b(r) * Qd(r) * Qd(r) / 2
         - c(r) * Qs(r) - d(r) * Qs(r) * Qs(r) / 2
         - tauO(r) * em(r) * Qs(r) )
    - sum((r,rr)$[not sameas(r,rr)], (t(r,rr) + kappa(r,rr)) * X(r,rr)) ;

* supply + imports >= demand + exports; multiplier is -P(r)
market_clearing(r).. Qs(r) + sum(rr$[not sameas(rr,r)], X(rr,r))
                     =g=
                     Qd(r) + sum(rr$[not sameas(rr,r)], X(r,rr)) ;

model simple /objfn, market_clearing/ ;

*!!!!! begin MCP: the KKT conditions of the QCP, each paired with its variable

positive variable P(r) "regional price (USD/t)" ;

equation foc_qs(r), foc_qd(r), foc_x(r,rr) ;

* producer: marginal cost + own carbon cost >= price          _|_ Qs(r) >= 0
foc_qs(r)..  c(r) + d(r) * Qs(r) + tauO(r) * em(r) =g= P(r) ;
* consumer: price >= marginal willingness to pay               _|_ Qd(r) >= 0
foc_qd(r)..  P(r) =g= a(r) + b(r) * Qd(r) ;
* arbitrage: origin price + transport + CBAM >= dest. price    _|_ X(r,rr) >= 0
foc_x(r,rr)$[not sameas(r,rr)].. P(r) + t(r,rr) + kappa(r,rr) =g= P(rr) ;

model cbam_mcp
/
market_clearing.P,
foc_qs.Qs,
foc_qd.Qd,
foc_x.X
/;

*!!!!! scenario loop: solve QCP, then confirm the point solves the MCP

parameter rep ;

loop(scen,
   tauO("EU") = swtab(scen,"ets") * tauEU ;
   tauO("CN") = tauCN ;
   tauO("RW") = tauRW ;
   kappa(r,rr) = 0 ;
   kappa(r,"EU")$[not sameas(r,"EU")] =
        swtab(scen,"cbam") * max(swtab(scen,"ets") * tauEU - tauO(r), 0) * em(r) ;

   solve simple using QCP maximizing W ;
   abort$[simple.modelstat > 2] "QCP failed", scen ;

*  warm start the MCP from the QCP; iterlim = 0 checks the QCP point is already
*  an MCP solution (same equilibrium, read off the complementarity system)
   rep(scen,"P_qcp",r) = -market_clearing.m(r) ;
   P.l(r) = -market_clearing.m(r) ;
   cbam_mcp.iterlim = 0 ;
   solve cbam_mcp using mcp ;
   abort$[cbam_mcp.modelstat > 2] "QCP point does not solve the MCP", scen ;

   rep(scen,"Qd",r)        = Qd.l(r) ;
   rep(scen,"Qs",r)        = Qs.l(r) ;
   rep(scen,"P",r)         = P.l(r) ;
   rep(scen,"Emis",r)      = em(r) * Qs.l(r) ;
   rep(scen,"Emis","World") = sum(r, em(r) * Qs.l(r)) ;
   rep(scen,"X_toEU",r)    = X.l(r,"EU") ;
   rep(scen,"X_CN",r)      = X.l("CN",r) ;
   rep(scen,"ImpInt","EU")$sum(r, X.l(r,"EU"))
                           = sum(r, em(r) * X.l(r,"EU")) / sum(r, X.l(r,"EU")) ;
   rep(scen,"Rev","ETS")   = tauO("EU") * em("EU") * Qs.l("EU") ;
   rep(scen,"Rev","CBAM")  = sum(r, kappa(r,"EU") * X.l(r,"EU")) ;
   rep(scen,"CS",r)        = -b(r) * Qd.l(r) * Qd.l(r) / 2 ;
*  producers sell at the home price; exporters earn zero arbitrage rent (foc_x)
   rep(scen,"PS",r)        = (P.l(r) - tauO(r) * em(r)) * Qs.l(r)
                             - c(r) * Qs.l(r) - d(r) * Qs.l(r) * Qs.l(r) / 2 ;
) ;

* calibration check: c_base reproduces the base point
abort$[smax(r, abs(rep("c_base","P",r) - pbar(r)) + abs(rep("c_base","Qs",r) - qbar_s(r))
             + abs(rep("c_base","Qd",r) - qbar_d(r))) > 1e-4] "calibration failed", rep ;

* QCP shadow prices must equal MCP prices
abort$[smax((scen,r), abs(rep(scen,"P",r) - rep(scen,"P_qcp",r))) > 1e-4] "QCP and MCP prices differ", rep ;

*!!! leakage and welfare, measured against c_base
rep(scen,"Leakage","rate")$[not sameas(scen,"c_base")] =
     sum(r$[not sameas(r,"EU")], rep(scen,"Emis",r) - rep("c_base","Emis",r))
   / (rep("c_base","Emis","EU") - rep(scen,"Emis","EU")) ;
rep(scen,"dW",r) = rep(scen,"CS",r) + rep(scen,"PS",r)
                 - rep("c_base","CS",r) - rep("c_base","PS",r) ;
rep(scen,"dW","EU") = rep(scen,"dW","EU") + rep(scen,"Rev","ETS") + rep(scen,"Rev","CBAM") ;
rep(scen,"dW","World") = sum(r, rep(scen,"dW",r)) ;
rep("c_cbam","Headline","dLeak") = rep("c_ets","Leakage","rate") - rep("c_cbam","Leakage","rate") ;

option rep:3:2:1 ;
display rep ;

execute_unload 'rep_cbam_mcp.gdx' ;

$exit
