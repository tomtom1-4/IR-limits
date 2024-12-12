#-
#: IncDir procedures
#: SmallExtension 200M
#: WorkSpace      4G
#: MaxTermSize 4M
Off Statistics;

Format 180;


s s12, s23, s13, p1q, p2q, p3q, ep, nl;
s s123, p123q;
s x1, x2, x3, z1, z2, z3;
v E0, E0conj, p1, p2, p3;
s Pi, muR;
cf Log, PolyLog, den, li2, log;
s logmu, marker, Ioperator;

#ifndef `splitting'
#define splitting "qqpqpbar"
#endif

#include P1`splitting'Exp.h
.sort
drop P1AvgExp, P0AvgExp;
id PolyLog(2, s12?) = li2(s12);

* switch to Sudakov parametrization
.sort
mul replace_(s12, x3);
mul replace_(s13, x2);
mul replace_(s23, x1);
.sort
mul replace_(p1q, z1);
mul replace_(p2q, z2);
mul replace_(p3q, z3);

.sort
mul replace_(x3, 1 - x1 - x2);
mul replace_(z3, 1 - z1 - z2);
.sort
l P1PolM2 = P1PolExp;
l P1PolM1 = P1PolExp;
l P1PolM0 = P1PolExp;
l IopM2 = P0PolExp;
l IopM1 = P0PolExp;
l IopM0 = P0PolExp;
.sort

#do i=0,2
  if(expression(P1PolM`i'));
    if(count(ep, 1)!=-`i') discard;
    id Log(-s12 - s13 - s23) = 0;
    id Log(-1) = 0;
  else if(expression(IopM`i'));
    if(count(Ioperator, 1) != 1) discard;
    id Ioperator = 1;
    if(count(ep, 1) != -`i') discard;
  endif;
#enddo

if(expression(P1PolExp));
  id Log(-1) = Log((s12 + s13 + s23)/muR^2) - i_*Pi;
  id Log(-s12 - s13 - s23) = Log((s12 + s13 + s23)/muR^2) - i_*Pi;
else if(expression(P0PolExp));
  id Log(-1) = Log((s12 + s13 + s23)/muR^2) - i_*Pi;
  id Log(-s12 - s13 - s23) = Log((s12 + s13 + s23)/muR^2) - i_*Pi;
  if(count(Ioperator, 1) != 1) discard;
  id Ioperator = 1;
endif;

.sort
mul replace_(i_, 0);
.sort
l test = P1PolExp - (P1PolM0 - Log((s12 + s13 + s23)/muR^2)*P1PolM1*ep + (Log((s12 + s13 + s23)/muR^2)^2/2 - Pi^2/2 )*P1PolM2*ep^2);
if(expression(test)) if(count(ep, 1) != 0) discard;

bracket ep, logmu, Log, marker;
print[] test;
.sort
drop P1PolExp, P0PolExp;
mul replace_(ep, 1);
Format O3;
id Log(s12?) = log(s12);
.sort

#write <P1_`splitting'.m> "std::vector<std::vector<std::complex<double>>> P1_`splitting'(LV<double> p1, LV<double> p2, LV<double> p3, double muR, int nl) {"
#write <P1_`splitting'.m> "  double s12 = 2.*p1*p2;"
#write <P1_`splitting'.m> "  double s13 = 2.*p1*p3;"
#write <P1_`splitting'.m> "  double s23 = 2.*p2*p3;"
#write <P1_`splitting'.m> "  double s123 = s12 + s13 + s23;"
#write <P1_`splitting'.m> "  LV<double> p = p1 + p2 + p3;"
#write <P1_`splitting'.m> "  LV<double> q(std::vector<double>({p.components[0], -p.components[1], -p.components[2], -p.components[3]}));"
#write <P1_`splitting'.m> "  double p1q = p1*q;"
#write <P1_`splitting'.m> "  double p2q = p2*q;"
#write <P1_`splitting'.m> "  double p3q = p3*q;"
#write <P1_`splitting'.m> "  double p123q = p123*q;"
*
#write <P1_`splitting'.m> "  double x1 = s23/s123;"
#write <P1_`splitting'.m> "  double x2 = s13/s123;"
#write <P1_`splitting'.m> "  double z1 = p1q/p123q;"
#write <P1_`splitting'.m> "  double z2 = p2q/p123q;"

#write <P1_`splitting'.m> "  std::complex<double> log = std::log(std::abs(s123/muR/muR)) + I*M_PI;"
#write <P1_`splitting'.m> "  std::vector<std::complex<double>> prefactor = {1, -log, -M_PI*M_PI/12. + 1./2.*std::pow(log, 2)};"
#write <P1_`splitting'.m> "  std::vector<double> P1`splitting' = {0., 0., 0.};"
#if ((`splitting' == "gqqbar") || (`splitting' == "ggg"))
  #write <P1_`splitting'.m> "  std::vector<std::vector<std::complex<double>>> output = {{0., 0.}, {0., 0.}};"
  #write <P1_`splitting'.m> "  for(int i = 0; i <= 1; i++) for(int j = 0; j <= 1; j++) {"
  #write <P1_`splitting'.m> "  LV<std::complex<double>> E0 = polarization(p, (i?1:-1));"
  #write <P1_`splitting'.m> "  LV<std::complex<double>> E0conj = polarization(p, (j?1:-1)).conj();"
#endif
#do i=0,2
  ExtraSymbols array,Z`i';
  Format C;
  #optimize P1PolM`i'
  #write <P1_`splitting'.m> "  double Z`i'[`optimmaxvar_'+1];"
  #write <P1_`splitting'.m> "%O"
  #write <P1_`splitting'.m> "  P1`splitting'[-`i'+2] = %e", P1PolM`i'
  #clearoptimize
  .sort
#enddo

#write <P1_`splitting'.m> "  std::vector<std::complex<double>> Ioperator = {0., 0., 0.};"
#do i=0,0
  ExtraSymbols array,Y`i';
  Format C;
  #optimize IopM`i'
  #write <P1_`splitting'.m> "  double Y`i'[`optimmaxvar_'+1];"
  #write <P1_`splitting'.m> "%O"
  #write <P1_`splitting'.m> "  Ioperator[-`i'+2] = %e", IopM`i'
  #clearoptimize
  .sort
#enddo

#write <P1_`splitting'.m> "  double P1 = 2.*std::real(prefactor[0]*P1`splitting'[2] + prefactor[1]*P1`splitting'[1] + prefactor[2]*P1`splitting'[0] - Ioperator[2]);"
#if (`splitting' == "qqpqpbar") || (`splitting' == "qqqbar") || (`splitting' == "qgg")
  #write <P1_`splitting'.m> "  std::vector<std::vector<std::complex<double>>> output = {{P1, 0.}, {0., P1}};"
#else if((`splitting' == "gqqbar") || (`splitting' == "ggg"))
  #write <P1_`splitting'.m> "  output[i][j] = P1;
  }"
#endif
#write <P1_`splitting'.m> "  return output;
}"
.sort


.end