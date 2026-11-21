$onText
 cbam_calib.gms -- One-point calibration (Milestone 2 Sec. 3)
 Assignments only, no declarations, so it can be $included inside a loop when
 epsD, etaS or e(i) change. Uses p0, D0, S0, e, tau0, epsD, etaS.

   Demand  p(j) = a(j) - b(j)*D(j):   at (D0, p0) with point elasticity epsD
   Supply  MC(i) = c(i) + d(i)*S(i):  at (S0, MC0) with point elasticity etaS
   MC0(i) = p0(i) - tau0*e(i): every domestic route is active in the base year,
   so the producer's own price minus the 2025 carbon cost (zero, A3) is its MC.
$offText

b(j)   = p0(j) / (abs(epsD)*D0(j));
a(j)   = p0(j) + b(j)*D0(j);
MC0(i) = p0(i) - tau0*e(i);
d(i)   = MC0(i) / (etaS*S0(i));
c(i)   = MC0(i) - d(i)*S0(i);
