std::vector<std::vector<std::complex<double>>> P0qqPqPbar(LV<double> p1, LV<double> p2, LV<double> p3) {
  double s12 = 2.*p1*p2;
  double s13 = 2.*p1*p3;
  double s23 = 2.*p2*p3;
  double s123 = s12 + s13 + s23;
  LV<double> p = p1 + p2 + p3;
  LV<double> q(std::vector<double>({p.components[0], -p.components[1], -p.components[2], -p.components[3]}));
  double p1q = p1*q; 
  double p2q = p2*q; 
  double p3q = p3*q;
  double P0qqPqPbar = (2*(2*std::pow(p1q,2)*(p2q + p3q)*s23*(s12 + s13 + s23) +
(p2q + p3q)*(std::pow(p2q,2)*(-2*std::pow(s13,2) + s12*s23 - s13*s23) +
std::pow(p3q,2)*(-2*std::pow(s12,2) - s12*s23 + s13*s23) + 2*p2q*p3q*(2*s12*s13
+ s12*s23 + s13*s23)) + p1q*(std::pow(p2q,2)*(-2*std::pow(s13,2) + 2*s12*s23 +
std::pow(s23,2)) + std::pow(p3q,2)*(-2*std::pow(s12,2) + 2*s13*s23 +
std::pow(s23,2)) + 2*p2q*p3q*(2*s12*s13 + 3*s12*s23 + 3*s13*s23 +
2*std::pow(s23,2)))))/(3.*std::pow(p2q + p3q,2)*(p1q + p2q +
p3q)*std::pow(s23,2));
  std::vector<std::vector<std::complex<double>>> output = {{P0qqPqPbar, 0.}, {0., P0qqPqPbar}};
  return output;
}

std::vector<std::vector<std::complex<double>>> P0_qgg(LV<double> p1, LV<double> p2, LV<double> p3) {
  double s12 = 2.*p1*p2;
  double s13 = 2.*p1*p3;
  double s23 = 2.*p2*p3;
  double s123 = s12 + s13 + s23;
  LV<double> p = p1 + p2 + p3;
  LV<double> q(std::vector<double>({p.components[0], -p.components[1], -p.components[2], -p.components[3]}));
  double p1q = p1*q; 
  double p2q = p2*q; 
  double p3q = p3*q;
  double P0qgg = (2*(-2*std::pow(p1q,3)*(p2q + p3q)*std::pow(s23,2)*(s12 + s13 +
s23)*(p2q*(2*s12 - 7*s13 + s23) + p3q*(-7*s12 + 2*s13 + s23)) +
2*std::pow(p1q,2)*(p2q + p3q)*s23*(s12 + s13 + s23)*(p2q*p3q*(18*std::pow(s12,2)
- 18*s12*s13 + 18*std::pow(s13,2) + 13*s12*s23 + 13*s13*s23 - 2*std::pow(s23,2))
+ std::pow(p2q,2)*(18*s12*s13 + 9*std::pow(s13,2) - 2*s12*s23 + 24*s13*s23 -
std::pow(s23,2)) + std::pow(p3q,2)*(9*std::pow(s12,2) - s23*(2*s13 + s23) +
6*s12*(3*s13 + 4*s23))) + p1q*(std::pow(p2q,4)*s23*(18*std::pow(s13,3) +
std::pow(s12,2)*(36*s13 - 2*s23) + 59*std::pow(s13,2)*s23 +
40*s13*std::pow(s23,2) - std::pow(s23,3) + s12*(54*std::pow(s13,2) + 75*s13*s23
- 3*std::pow(s23,2))) + std::pow(p3q,4)*s23*(18*std::pow(s12,3) +
std::pow(s12,2)*(54*s13 + 59*s23) - s23*(2*std::pow(s13,2) + 3*s13*s23 +
std::pow(s23,2)) + s12*(36*std::pow(s13,2) + 75*s13*s23 + 40*std::pow(s23,2))) +
p2q*std::pow(p3q,3)*(18*std::pow(s12,3)*(2*s13 + 5*s23) +
2*std::pow(s12,2)*s23*(72*s13 + 85*s23) + 3*s12*s23*(18*std::pow(s13,2) +
45*s13*s23 + 28*std::pow(s23,2)) + s23*(36*std::pow(s13,3) +
39*std::pow(s13,2)*s23 + 7*s13*std::pow(s23,2) - 4*std::pow(s23,3))) +
2*std::pow(p2q,2)*std::pow(p3q,2)*(54*std::pow(s12,3)*s23 +
3*s12*s23*(6*std::pow(s13,2) + 14*s13*s23 + 9*std::pow(s23,2)) +
std::pow(s12,2)*(-36*std::pow(s13,2) + 18*s13*s23 + 76*std::pow(s23,2)) +
s23*(54*std::pow(s13,3) + 76*std::pow(s13,2)*s23 + 27*s13*std::pow(s23,2) -
3*std::pow(s23,3))) + std::pow(p2q,3)*p3q*(36*std::pow(s12,3)*s23 +
3*std::pow(s12,2)*s23*(18*s13 + 13*s23) + 2*s23*(45*std::pow(s13,3) +
85*std::pow(s13,2)*s23 + 42*s13*std::pow(s23,2) - 2*std::pow(s23,3)) +
s12*(36*std::pow(s13,3) + 144*std::pow(s13,2)*s23 + 135*s13*std::pow(s23,2) +
7*std::pow(s23,3)))) + (p2q + p3q)*(std::pow(p2q,4)*s13*s23*(18*std::pow(s12,2)
+ 27*s12*s13 + 9*std::pow(s13,2) + 35*s12*s23 + 26*s13*s23 + 17*std::pow(s23,2))
+ std::pow(p3q,4)*s12*s23*(9*std::pow(s12,2) + 27*s12*s13 + 18*std::pow(s13,2) +
26*s12*s23 + 35*s13*s23 + 17*std::pow(s23,2)) +
p2q*std::pow(p3q,3)*(9*std::pow(s12,3)*(4*s13 + 5*s23) +
9*std::pow(s12,2)*s23*(11*s13 + 8*s23) + s13*s23*(18*std::pow(s13,2) +
19*s13*s23 + 9*std::pow(s23,2)) + s12*s23*(36*std::pow(s13,2) + 80*s13*s23 +
35*std::pow(s23,2))) + std::pow(p2q,2)*std::pow(p3q,2)*(45*std::pow(s12,3)*s23 +
9*s12*s23*(3*std::pow(s13,2) + 8*s13*s23 + 3*std::pow(s23,2)) +
s13*s23*(45*std::pow(s13,2) + 56*s13*s23 + 27*std::pow(s23,2)) +
std::pow(s12,2)*(-72*std::pow(s13,2) + 27*s13*s23 + 56*std::pow(s23,2))) +
std::pow(p2q,3)*p3q*(18*std::pow(s12,3)*s23 + std::pow(s12,2)*s23*(36*s13 +
19*s23) + s13*s23*(45*std::pow(s13,2) + 72*s13*s23 + 35*std::pow(s23,2)) +
s12*(36*std::pow(s13,3) + 99*std::pow(s13,2)*s23 + 80*s13*std::pow(s23,2) +
9*std::pow(s23,3))))))/(9.*p2q*p3q*std::pow(p2q + p3q,2)*(p1q + p2q +
p3q)*s12*s13*std::pow(s23,2));
  std::vector<std::vector<std::complex<double>>> output = {{P0qgg, 0.}, {0., P0qgg}};
  return output;
}

std::vector<std::vector<std::complex<double>>> P0_qqqbar(LV<double> p1, LV<double> p2, LV<double> p3) {
  double s12 = 2.*p1*p2;
  double s13 = 2.*p1*p3;
  double s23 = 2.*p2*p3;
  double s123 = s12 + s13 + s23;
  LV<double> p = p1 + p2 + p3;
  LV<double> q(std::vector<double>({p.components[0], -p.components[1], -p.components[2], -p.components[3]}));
  double p1q = p1*q; 
  double p2q = p2*q; 
  double p3q = p3*q;
  double P0qqqbar = (2*(std::pow(p1q,4)*(p2q + p3q)*s13*(6*s13 - s23)*s23*(s12 +
s13 + s23) - p3q*(p2q + p3q)*(std::pow(p2q,3)*s13*(s13 - 6*s23)*s23*(s12 + s13 +
s23) - p2q*std::pow(p3q,2)*(2*std::pow(s12,2)*(s13 - 6*s23)*s23 +
s12*s13*(12*std::pow(s13,2) + 4*s13*s23 - 5*std::pow(s23,2)) +
3*s13*s23*(2*std::pow(s13,2) + s13*s23 + 3*std::pow(s23,2))) +
std::pow(p3q,3)*(3*s12*s13*s23*(s13 + s23) - 3*s13*s23*(std::pow(s13,2) +
std::pow(s23,2)) + std::pow(s12,2)*(6*std::pow(s13,2) - 2*s13*s23 +
6*std::pow(s23,2))) + std::pow(p2q,2)*p3q*(6*std::pow(s13,4) +
3*std::pow(s13,3)*s23 + 6*std::pow(s12,2)*std::pow(s23,2) -
std::pow(s13,2)*s23*(2*s12 + 9*s23) - s13*s23*(std::pow(s12,2) + 5*s12*s23 +
12*std::pow(s23,2)))) + std::pow(p1q,3)*(std::pow(p2q,2)*(-6*std::pow(s13,4) +
std::pow(s13,3)*s23 + s13*(s12 - 3*s23)*std::pow(s23,2) - 6*std::pow(s23,4) +
std::pow(s13,2)*s23*(5*s12 + 4*s23)) + p2q*p3q*(std::pow(s12,2)*s13*s23 +
2*s12*s13*(6*std::pow(s13,2) + 14*s13*s23 + std::pow(s23,2)) +
s23*(31*std::pow(s13,3) + 24*std::pow(s13,2)*s23 - 7*s13*std::pow(s23,2) -
12*std::pow(s23,3))) + std::pow(p3q,2)*(std::pow(s12,2)*s13*(-6*s13 + s23) +
s12*s13*s23*(11*s13 + s23) + 2*s23*(9*std::pow(s13,3) + 7*std::pow(s13,2)*s23 -
2*s13*std::pow(s23,2) - 3*std::pow(s23,3)))) +
std::pow(p1q,2)*(std::pow(p2q,3)*(-6*std::pow(s13,4) - 3*std::pow(s13,3)*s23 -
6*std::pow(s23,4) + s13*std::pow(s23,2)*(5*s12 + s23) + std::pow(s13,2)*s23*(s12
+ 4*s23)) + 3*std::pow(p3q,3)*(std::pow(s12,2)*s13*(-6*s13 + s23) +
2*s12*std::pow(s23,2)*(s13 + 2*s23) + s23*(7*std::pow(s13,3) +
4*std::pow(s13,2)*s23 + s13*std::pow(s23,2) - 2*std::pow(s23,3))) +
2*std::pow(p2q,2)*p3q*(-9*std::pow(s13,4) + std::pow(s12,2)*s13*s23 +
3*std::pow(s13,3)*s23 + 9*std::pow(s13,2)*std::pow(s23,2) +
3*s13*std::pow(s23,3) - 9*std::pow(s23,4) + s12*(6*std::pow(s13,3) +
8*std::pow(s13,2)*s23 + 8*s13*std::pow(s23,2) + 6*std::pow(s23,3))) +
p2q*std::pow(p3q,2)*(std::pow(s12,2)*s13*(-6*s13 + 5*s23) +
2*s23*(27*std::pow(s13,3) + 19*std::pow(s13,2)*s23 + 4*s13*std::pow(s23,2) -
9*std::pow(s23,3)) + s12*(36*std::pow(s13,3) + 39*std::pow(s13,2)*s23 +
17*s13*std::pow(s23,2) + 24*std::pow(s23,3)))) + p1q*(-(std::pow(p2q,4)*s13*(s13
- 6*s23)*s23*(s12 + s13 + s23)) +
std::pow(p3q,4)*(-2*std::pow(s12,2)*(9*std::pow(s13,2) - 2*s13*s23 +
3*std::pow(s23,2)) + 3*s13*s23*(4*std::pow(s13,2) + s13*s23 + 3*std::pow(s23,2))
+ s12*s23*(-8*std::pow(s13,2) + s13*s23 + 12*std::pow(s23,2))) +
std::pow(p2q,3)*p3q*(-12*std::pow(s13,4) - 7*std::pow(s13,3)*s23 +
12*s12*std::pow(s23,3) + 2*std::pow(s13,2)*s23*(s12 + 12*s23) +
s13*s23*(std::pow(s12,2) + 28*s12*s23 + 31*std::pow(s23,2))) +
2*p2q*std::pow(p3q,3)*(std::pow(s12,2)*(-6*std::pow(s13,2) + 4*s13*s23 -
6*std::pow(s23,2)) + s13*s23*(19*std::pow(s13,2) + 14*s13*s23 +
19*std::pow(s23,2)) + 9*s12*(2*std::pow(s13,3) + std::pow(s13,2)*s23 +
s13*std::pow(s23,2) + 2*std::pow(s23,3))) +
std::pow(p2q,2)*std::pow(p3q,2)*(std::pow(s12,2)*(5*s13 - 6*s23)*s23 +
2*s13*(-9*std::pow(s13,3) + 4*std::pow(s13,2)*s23 + 19*s13*std::pow(s23,2) +
27*std::pow(s23,3)) + s12*(24*std::pow(s13,3) + 17*std::pow(s13,2)*s23 +
39*s13*std::pow(s23,2) + 36*std::pow(s23,3))))))/(9.*std::pow(p1q +
p3q,2)*std::pow(p2q + p3q,2)*(p1q + p2q + p3q)*std::pow(s13,2)*std::pow(s23,2));
  std::vector<std::vector<std::complex<double>>> output = {{P0qqqbar, 0.}, {0., P0qqqbar}};
  return output;
}

std::vector<std::vector<std::complex<double>>> P0_gqqbar(LV<double> p1, LV<double> p2, LV<double> p3) {
  double s12 = 2.*p1*p2;
  double s13 = 2.*p1*p3;
  double s23 = 2.*p2*p3;
  double s123 = s12 + s13 + s23;
  LV<double> p = p1 + p2 + p3;
  LV<double> q(std::vector<double>({p.components[0], -p.components[1], -p.components[2], -p.components[3]}));
  double p1q = p1*q; 
  double p2q = p2*q; 
  double p3q = p3*q;
  std::vector<std::vector<std::complex<double>>> output = {{0., 0.}, {0., 0.}};
  for(int i = 0; i <= 1; i++) for(int j = 0; j <= 1; j++) {
    LV<std::complex<double>> E0 = polarization(p, (i?1:-1));
    LV<std::complex<double>> E0conj = polarization(p, (j?1:-1)).conj();
  double output[i][j] = ((-9*std::pow(p2q + p3q,2)*s23*(s12 + s13 +
s23)*(p2q*s13*(3*s12 + s13 + 2*s23) + p3q*s12*(s12 + 3*s13 + 2*s23)) -
9*std::pow(p1q,2)*(p2q + p3q)*s23*(s13*s23*(s13 + s23) + std::pow(s12,2)*(4*s13
+ s23) + s12*(4*std::pow(s13,2) + 6*s13*s23 + std::pow(s23,2))) +
p1q*(std::pow(p3q,2)*(9*std::pow(s12,3)*(4*s13 + s23) +
2*std::pow(s12,2)*s23*(18*s13 + 5*s23) + std::pow(s23,2)*(std::pow(s13,2) +
11*s13*s23 + 2*std::pow(s23,2)) + s12*s23*(-9*std::pow(s13,2) + 27*s13*s23 +
11*std::pow(s23,2))) + p2q*p3q*(9*std::pow(s12,3)*s23 +
std::pow(s12,2)*(-72*std::pow(s13,2) - 45*s13*s23 + 11*std::pow(s23,2)) +
s12*s23*(-45*std::pow(s13,2) + 18*s13*s23 + 22*std::pow(s23,2)) +
s23*(9*std::pow(s13,3) + 11*std::pow(s13,2)*s23 + 22*s13*std::pow(s23,2) +
4*std::pow(s23,3))) + std::pow(p2q,2)*(std::pow(s12,2)*s23*(-9*s13 + s23) +
s23*(9*std::pow(s13,3) + 10*std::pow(s13,2)*s23 + 11*s13*std::pow(s23,2) +
2*std::pow(s23,3)) + s12*(36*std::pow(s13,3) + 36*std::pow(s13,2)*s23 +
27*s13*std::pow(s23,2) +
11*std::pow(s23,3)))))*Dot(E0,E0conj))/(12.*p1q*std::pow(p2q +
p3q,2)*s12*s13*std::pow(s23,2)) + ((s12 + s13 +
s23)*(-36*std::pow(p2q,2)*s13*(2*s12 + s23) + p1q*p3q*s23*(27*s12 - 9*s13 +
2*s23) + p1q*p2q*s23*(9*s12 + 9*s13 + 2*s23) - 18*p2q*p3q*(s13*s23 + s12*(4*s13
+ s23)))*Dot(p1,E0conj)*Dot(p2,E0))/(6.*p1q*(p2q + p3q)*s12*s13*std::pow(s23,2))
- (2*(s12 + s13 + s23)*(9*std::pow(p2q,2)*s13*(2*s12 + s23) - p1q*p3q*s23*(9*s12
+ s23) + 9*std::pow(p3q,2)*s12*(2*s13 + s23) - p1q*p2q*s23*(9*s13 + s23) +
9*p2q*p3q*(s13*s23 + s12*(4*s13 + s23)))*Dot(p2,E0)*Dot(p2,E0conj))/(3.*p1q*(p2q
+ p3q)*s12*s13*std::pow(s23,2)) + Dot(p1,E0)*(((s12 + s13 +
s23)*(-18*std::pow(p2q,2)*s13*(2*s12 + s23) + p1q*p2q*s23*(9*s12 + s23) +
p1q*p3q*s23*(9*s12 + s23))*Dot(p1,E0conj))/(3.*p1q*(p2q +
p3q)*s12*s13*std::pow(s23,2)) + ((s12 + s13 +
s23)*(-36*std::pow(p2q,2)*s13*(2*s12 + s23) + p1q*p3q*s23*(27*s12 - 9*s13 +
2*s23) + p1q*p2q*s23*(9*s12 + 9*s13 + 2*s23) - 18*p2q*p3q*(s13*s23 + s12*(4*s13
+ s23)))*Dot(p2,E0conj))/(6.*p1q*(p2q + p3q)*s12*s13*std::pow(s23,2)));
  }
  return output;
}

std::vector<std::vector<std::complex<double>>> P0_ggg(LV<double> p1, LV<double> p2, LV<double> p3) {
  double s12 = 2.*p1*p2;
  double s13 = 2.*p1*p3;
  double s23 = 2.*p2*p3;
  double s123 = s12 + s13 + s23;
  LV<double> p = p1 + p2 + p3;
  LV<double> q(std::vector<double>({p.components[0], -p.components[1], -p.components[2], -p.components[3]}));
  double p1q = p1*q; 
  double p2q = p2*q; 
  double p3q = p3*q;
  std::vector<std::vector<std::complex<double>>> output = {{0., 0.}, {0., 0.}};
  for(int i = 0; i <= 1; i++) for(int j = 0; j <= 1; j++) {
    LV<std::complex<double>> E0 = polarization(p, (i?1:-1));
    LV<std::complex<double>> E0conj = polarization(p, (j?1:-1)).conj();
  double output[i][j] = (-9*(std::pow(p1q,7)*(p2q +
p3q)*s12*s13*std::pow(s23,2)*(s12 + s13 + s23)*(p3q*(3*s12 + 2*s13 + s23) +
p2q*(2*s12 + 3*s13 + s23)) + std::pow(p2q,2)*std::pow(p3q,2)*std::pow(p2q +
p3q,2)*(std::pow(p2q,2) + std::pow(p3q,2))*s12*s13*s23*(s12 + s13 +
s23)*(p2q*s13*(3*s12 + s13 + 2*s23) + p3q*s12*(s12 + 3*s13 + 2*s23)) +
std::pow(p1q,6)*(p2q + p3q)*s12*s13*s23*(s12 + s13 +
s23)*(2*p2q*p3q*(std::pow(s12,2) - s12*s13 + std::pow(s13,2) + 6*s12*s23 +
6*s13*s23 + 2*std::pow(s23,2)) + std::pow(p2q,2)*(std::pow(s13,2) + 9*s13*s23 +
2*std::pow(s23,2) + 2*s12*(s13 + 2*s23)) + std::pow(p3q,2)*(std::pow(s12,2) +
2*s23*(2*s13 + s23) + s12*(2*s13 + 9*s23))) +
p1q*p2q*p3q*(std::pow(p2q,6)*s12*std::pow(s13,2)*s23*(5*std::pow(s12,2) +
7*s12*s13 + 2*std::pow(s13,2) + 10*s12*s23 + 7*s13*s23 + 5*std::pow(s23,2)) +
std::pow(p3q,6)*std::pow(s12,2)*s13*s23*(2*std::pow(s12,2) + 7*s12*(s13 + s23) +
5*std::pow(s13 + s23,2)) + p2q*std::pow(p3q,5)*s12*s13*s23*(6*std::pow(s12,3) +
3*std::pow(s13,3) + 3*std::pow(s13,2)*s23 + 2*s13*std::pow(s23,2) +
2*std::pow(s23,3) + std::pow(s12,2)*(27*s13 + 22*s23) + s12*(24*std::pow(s13,2)
+ 37*s13*s23 + 18*std::pow(s23,2))) +
std::pow(p2q,5)*p3q*s12*s13*s23*(3*std::pow(s12,3) + 3*std::pow(s12,2)*(8*s13 +
s23) + s12*(27*std::pow(s13,2) + 37*s13*s23 + 2*std::pow(s23,2)) +
2*(3*std::pow(s13,3) + 11*std::pow(s13,2)*s23 + 9*s13*std::pow(s23,2) +
std::pow(s23,3))) +
std::pow(p2q,2)*std::pow(p3q,4)*(2*std::pow(s13,4)*std::pow(s23,2) +
std::pow(s12,3)*s13*s23*(50*s13 + 33*s23) + std::pow(s12,4)*(2*std::pow(s13,2) +
11*s13*s23 + 2*std::pow(s23,2)) + 2*std::pow(s12,2)*s13*s23*(23*std::pow(s13,2)
+ 33*s13*s23 + 12*std::pow(s23,2)) + s12*s13*s23*(9*std::pow(s13,3) +
14*std::pow(s13,2)*s23 + 7*s13*std::pow(s23,2) + 4*std::pow(s23,3))) +
std::pow(p2q,3)*std::pow(p3q,3)*(4*std::pow(s13,4)*std::pow(s23,2) +
std::pow(s12,4)*s23*(13*s13 + 4*s23) +
std::pow(s12,2)*s13*s23*(53*std::pow(s13,2) + 76*s13*s23 + 16*std::pow(s23,2)) +
std::pow(s12,3)*s13*(-4*std::pow(s13,2) + 53*s13*s23 + 29*std::pow(s23,2)) +
s12*s13*s23*(13*std::pow(s13,3) + 29*std::pow(s13,2)*s23 +
16*s13*std::pow(s23,2) + 4*std::pow(s23,3))) +
std::pow(p2q,4)*std::pow(p3q,2)*(2*std::pow(s13,4)*std::pow(s23,2) +
std::pow(s12,4)*s23*(9*s13 + 2*s23) + 2*std::pow(s12,3)*s13*s23*(23*s13 + 7*s23)
+ s12*s13*s23*(11*std::pow(s13,3) + 33*std::pow(s13,2)*s23 +
24*s13*std::pow(s23,2) + 4*std::pow(s23,3)) +
std::pow(s12,2)*s13*(2*std::pow(s13,3) + 50*std::pow(s13,2)*s23 +
66*s13*std::pow(s23,2) + 7*std::pow(s23,3)))) + std::pow(p1q,2)*(p2q +
p3q)*(std::pow(p2q,6)*s12*std::pow(s13,2)*s23*(2*std::pow(s12,2) + 3*s12*s13 +
std::pow(s13,2) + 5*s12*s23 + 4*s13*s23 + 3*std::pow(s23,2)) +
std::pow(p3q,6)*std::pow(s12,2)*s13*s23*(std::pow(s12,2) + 3*s12*s13 +
2*std::pow(s13,2) + 4*s12*s23 + 5*s13*s23 + 3*std::pow(s23,2)) +
std::pow(p2q,5)*p3q*s12*s13*s23*(2*std::pow(s12,3) + 5*std::pow(s13,3) +
23*std::pow(s13,2)*s23 + 21*s13*std::pow(s23,2) + 3*std::pow(s23,3) +
2*std::pow(s12,2)*(8*s13 + s23) + s12*(19*std::pow(s13,2) + 32*s13*s23 +
3*std::pow(s23,2))) + p2q*std::pow(p3q,5)*s12*s13*s23*(5*std::pow(s12,3) +
2*std::pow(s13,3) + 2*std::pow(s13,2)*s23 + 3*s13*std::pow(s23,2) +
3*std::pow(s23,3) + std::pow(s12,2)*(19*s13 + 23*s23) + s12*(16*std::pow(s13,2)
+ 32*s13*s23 + 21*std::pow(s23,2))) +
2*std::pow(p2q,3)*std::pow(p3q,3)*(2*std::pow(s13,3)*(s13 - s23)*std::pow(s23,2)
+ 2*std::pow(s12,4)*s23*(4*s13 + s23) +
std::pow(s12,2)*s13*s23*(17*std::pow(s13,2) + 23*s13*s23 + 6*std::pow(s23,2)) +
std::pow(s12,3)*(-4*std::pow(s13,3) + 17*std::pow(s13,2)*s23 +
11*s13*std::pow(s23,2) - 2*std::pow(s23,3)) + s12*s13*s23*(8*std::pow(s13,3) +
11*std::pow(s13,2)*s23 + 6*s13*std::pow(s23,2) + 5*std::pow(s23,3))) +
std::pow(p2q,2)*std::pow(p3q,4)*(-4*std::pow(s13,3)*std::pow(s23,3) +
4*std::pow(s12,3)*s13*s23*(11*s13 + 10*s23) + std::pow(s12,4)*(4*std::pow(s13,2)
+ 13*s13*s23 + 4*std::pow(s23,2)) + std::pow(s12,2)*s13*s23*(36*std::pow(s13,2)
+ 58*s13*s23 + 31*std::pow(s23,2)) + s12*s13*s23*(9*std::pow(s13,3) +
5*std::pow(s13,2)*s23 + 4*s13*std::pow(s23,2) + 8*std::pow(s23,3))) +
std::pow(p2q,4)*std::pow(p3q,2)*(9*std::pow(s12,4)*s13*s23 +
4*std::pow(s13,4)*std::pow(s23,2) + std::pow(s12,3)*s23*(36*std::pow(s13,2) +
5*s13*s23 - 4*std::pow(s23,2)) + 2*std::pow(s12,2)*s13*(2*std::pow(s13,3) +
22*std::pow(s13,2)*s23 + 29*s13*std::pow(s23,2) + 2*std::pow(s23,3)) +
s12*s13*s23*(13*std::pow(s13,3) + 40*std::pow(s13,2)*s23 +
31*s13*std::pow(s23,2) + 8*std::pow(s23,3)))) +
std::pow(p1q,5)*(2*std::pow(p2q,4)*s12*s13*s23*(std::pow(s13,3) +
7*std::pow(s13,2)*s23 + 7*s13*std::pow(s23,2) + std::pow(s23,3) +
2*std::pow(s12,2)*(s13 + s23) + s12*(3*std::pow(s13,2) + 10*s13*s23 +
3*std::pow(s23,2))) + 2*std::pow(p3q,4)*s12*s13*s23*(std::pow(s12,3) +
std::pow(s12,2)*(3*s13 + 7*s23) + s23*(2*std::pow(s13,2) + 3*s13*s23 +
std::pow(s23,2)) + s12*(2*std::pow(s13,2) + 10*s13*s23 + 7*std::pow(s23,2))) +
p2q*std::pow(p3q,3)*(2*std::pow(s13,2)*std::pow(s23,4) +
std::pow(s12,4)*s13*(2*s13 + 9*s23) + 2*std::pow(s12,3)*s13*s23*(7*s13 + 23*s23)
+ std::pow(s12,2)*s23*(7*std::pow(s13,3) + 66*std::pow(s13,2)*s23 +
50*s13*std::pow(s23,2) + 2*std::pow(s23,3)) + s12*s13*s23*(4*std::pow(s13,3) +
24*std::pow(s13,2)*s23 + 33*s13*std::pow(s23,2) + 11*std::pow(s23,3))) +
std::pow(p2q,2)*std::pow(p3q,2)*(11*std::pow(s12,4)*s13*s23 +
4*std::pow(s13,2)*std::pow(s23,4) + std::pow(s12,3)*s13*(-4*std::pow(s13,2) +
7*s13*s23 + 52*std::pow(s23,2)) + std::pow(s12,2)*s23*(7*std::pow(s13,3) +
90*std::pow(s13,2)*s23 + 63*s13*std::pow(s23,2) + 4*std::pow(s23,3)) +
s12*s13*s23*(11*std::pow(s13,3) + 52*std::pow(s13,2)*s23 +
63*s13*std::pow(s23,2) + 18*std::pow(s23,3))) +
std::pow(p2q,3)*p3q*(4*std::pow(s12,4)*s13*s23 +
2*std::pow(s13,2)*std::pow(s23,4) + std::pow(s12,3)*s13*s23*(7*s13 + 24*s23) +
s12*s13*s23*(9*std::pow(s13,3) + 46*std::pow(s13,2)*s23 + 50*s13*std::pow(s23,2)
+ 11*std::pow(s23,3)) + std::pow(s12,2)*(2*std::pow(s13,4) +
14*std::pow(s13,3)*s23 + 66*std::pow(s13,2)*std::pow(s23,2) +
33*s13*std::pow(s23,3) + 2*std::pow(s23,4)))) + std::pow(p1q,4)*(p2q +
p3q)*(2*std::pow(p2q,4)*s12*s13*s23*(std::pow(s13,3) + 7*std::pow(s13,2)*s23 +
7*s13*std::pow(s23,2) + std::pow(s23,3) + 2*std::pow(s12,2)*(s13 + s23) +
s12*(3*std::pow(s13,2) + 10*s13*s23 + 3*std::pow(s23,2))) +
2*std::pow(p3q,4)*s12*s13*s23*(std::pow(s12,3) + std::pow(s12,2)*(3*s13 + 7*s23)
+ s23*(2*std::pow(s13,2) + 3*s13*s23 + std::pow(s23,2)) + s12*(2*std::pow(s13,2)
+ 10*s13*s23 + 7*std::pow(s23,2))) +
p2q*std::pow(p3q,3)*(4*std::pow(s13,2)*std::pow(s23,4) +
std::pow(s12,4)*s13*(4*s13 + 11*s23) + std::pow(s12,3)*s23*(23*std::pow(s13,2) +
39*s13*s23 - 4*std::pow(s23,2)) + std::pow(s12,2)*s13*s23*(12*std::pow(s13,2) +
56*s13*s23 + 39*std::pow(s23,2)) + s12*s13*s23*(4*std::pow(s13,3) +
12*std::pow(s13,2)*s23 + 23*s13*std::pow(s23,2) + 11*std::pow(s23,3))) +
std::pow(p2q,2)*std::pow(p3q,2)*(14*std::pow(s12,4)*s13*s23 +
4*std::pow(s13,2)*std::pow(s23,3)*(-s13 + s23) +
std::pow(s12,3)*(-8*std::pow(s13,3) + 4*std::pow(s13,2)*s23 +
31*s13*std::pow(s23,2) - 4*std::pow(s23,3)) +
std::pow(s12,2)*s23*(4*std::pow(s13,3) + 48*std::pow(s13,2)*s23 +
39*s13*std::pow(s23,2) + 4*std::pow(s23,3)) + s12*s13*s23*(14*std::pow(s13,3) +
31*std::pow(s13,2)*s23 + 39*s13*std::pow(s23,2) + 18*std::pow(s23,3))) +
std::pow(p2q,3)*p3q*(4*std::pow(s12,4)*s13*s23 -
4*std::pow(s13,3)*std::pow(s23,3) + 12*std::pow(s12,3)*s13*s23*(s13 + s23) +
s12*s13*s23*(11*std::pow(s13,3) + 39*std::pow(s13,2)*s23 +
39*s13*std::pow(s23,2) + 11*std::pow(s23,3)) +
std::pow(s12,2)*(4*std::pow(s13,4) + 23*std::pow(s13,3)*s23 +
56*std::pow(s13,2)*std::pow(s23,2) + 23*s13*std::pow(s23,3) +
4*std::pow(s23,4)))) +
std::pow(p1q,3)*(std::pow(p2q,6)*s12*s13*s23*(2*std::pow(s13,3) +
11*std::pow(s13,2)*s23 + 10*s13*std::pow(s23,2) + std::pow(s23,3) +
2*std::pow(s12,2)*(2*s13 + s23) + 3*s12*(2*std::pow(s13,2) + 5*s13*s23 +
std::pow(s23,2))) + std::pow(p3q,6)*s12*s13*s23*(2*std::pow(s12,3) +
std::pow(s12,2)*(6*s13 + 11*s23) + s23*(2*std::pow(s13,2) + 3*s13*s23 +
std::pow(s23,2)) + s12*(4*std::pow(s13,2) + 15*s13*s23 + 10*std::pow(s23,2))) +
p2q*std::pow(p3q,5)*(2*std::pow(s13,2)*std::pow(s23,4) +
std::pow(s12,3)*s13*s23*(33*s13 + 50*s23) + std::pow(s12,4)*(2*std::pow(s13,2) +
11*s13*s23 + 2*std::pow(s23,2)) + 2*std::pow(s12,2)*s13*s23*(12*std::pow(s13,2)
+ 33*s13*s23 + 23*std::pow(s23,2)) + s12*s13*s23*(4*std::pow(s13,3) +
7*std::pow(s13,2)*s23 + 14*s13*std::pow(s23,2) + 9*std::pow(s23,3))) +
std::pow(p2q,2)*std::pow(p3q,4)*(4*std::pow(s13,2)*std::pow(s23,3)*(-2*s13 +
s23) + std::pow(s12,4)*(8*std::pow(s13,2) + 29*s13*s23 + 4*std::pow(s23,2)) +
std::pow(s12,2)*s13*s23*(43*std::pow(s13,2) + 104*s13*s23 + 70*std::pow(s23,2))
+ std::pow(s12,3)*(-4*std::pow(s13,3) + 62*std::pow(s13,2)*s23 +
78*s13*std::pow(s23,2) - 8*std::pow(s23,3)) + s12*s13*s23*(18*std::pow(s13,3) +
16*std::pow(s13,2)*s23 + 27*s13*std::pow(s23,2) + 25*std::pow(s23,3))) +
2*std::pow(p2q,3)*std::pow(p3q,3)*(std::pow(s13,2)*std::pow(s23,2)*(std::pow(s13
,2) - 8*s13*s23 + std::pow(s23,2)) + std::pow(s12,4)*(std::pow(s13,2) +
17*s13*s23 + std::pow(s23,2)) + std::pow(s12,3)*(-8*std::pow(s13,3) +
25*std::pow(s13,2)*s23 + 25*s13*std::pow(s23,2) - 8*std::pow(s23,3)) +
s12*s13*s23*(17*std::pow(s13,3) + 25*std::pow(s13,2)*s23 +
25*s13*std::pow(s23,2) + 17*std::pow(s23,3)) + std::pow(s12,2)*(std::pow(s13,4)
+ 25*std::pow(s13,3)*s23 + 51*std::pow(s13,2)*std::pow(s23,2) +
25*s13*std::pow(s23,3) + std::pow(s23,4))) +
std::pow(p2q,5)*p3q*(4*std::pow(s12,4)*s13*s23 +
2*std::pow(s13,4)*std::pow(s23,2) + std::pow(s12,3)*s13*s23*(24*s13 + 7*s23) +
s12*s13*s23*(11*std::pow(s13,3) + 50*std::pow(s13,2)*s23 +
46*s13*std::pow(s23,2) + 9*std::pow(s23,3)) + std::pow(s12,2)*(2*std::pow(s13,4)
+ 33*std::pow(s13,3)*s23 + 66*std::pow(s13,2)*std::pow(s23,2) +
14*s13*std::pow(s23,3) + 2*std::pow(s23,4))) +
std::pow(p2q,4)*std::pow(p3q,2)*(18*std::pow(s12,4)*s13*s23 +
4*std::pow(s13,3)*(s13 - 2*s23)*std::pow(s23,2) +
std::pow(s12,3)*(-4*std::pow(s13,3) + 43*std::pow(s13,2)*s23 +
16*s13*std::pow(s23,2) - 8*std::pow(s23,3)) + s12*s13*s23*(29*std::pow(s13,3) +
78*std::pow(s13,2)*s23 + 70*s13*std::pow(s23,2) + 25*std::pow(s23,3)) +
std::pow(s12,2)*(8*std::pow(s13,4) + 62*std::pow(s13,3)*s23 +
104*std::pow(s13,2)*std::pow(s23,2) + 27*s13*std::pow(s23,3) +
4*std::pow(s23,4)))))*Dot(E0,E0conj))/(p1q*p2q*std::pow(p1q +
p2q,2)*p3q*std::pow(p1q + p3q,2)*std::pow(p2q +
p3q,2)*std::pow(s12,2)*std::pow(s13,2)*std::pow(s23,2)) + (18*(s12 + s13 +
s23)*(2*std::pow(p1q,4)*p3q*(p2q + p3q)*s12*(2*s12 + s13)*std::pow(s23,2) +
std::pow(p2q,3)*std::pow(p3q,2)*s12*s13*(p3q*s13*s23 + 2*p2q*s13*(2*s12 + s23) +
p3q*s12*(4*s13 + s23)) +
p1q*std::pow(p2q,2)*p3q*s12*s13*(2*std::pow(p2q,2)*s13*(2*s12 + s23) +
std::pow(p3q,2)*(4*s12*s13 - 2*s12*s23 + 3*s13*s23 + std::pow(s23,2)) +
p2q*p3q*(8*s12*s13 - s12*s23 + 4*s13*s23 + std::pow(s23,2))) +
std::pow(p1q,3)*s23*(p2q*std::pow(p3q,2)*s12*(-(s12*s13) + std::pow(s13,2) +
8*s12*s23 + 4*s13*s23) - std::pow(p2q,2)*p3q*(std::pow(s12,2)*(s13 - 4*s23) +
s12*s13*(s13 - 2*s23) + 4*std::pow(s13,2)*s23) - std::pow(p2q,3)*s13*(4*s13*s23
+ s12*(s13 + s23)) + std::pow(p3q,3)*s12*(s13*s23 + s12*(s13 + 4*s23))) +
std::pow(p1q,2)*p2q*p3q*(std::pow(p3q,2)*s12*s23*(-2*s12*s13 + std::pow(s13,2) +
4*s12*s23 + 3*s13*s23) + std::pow(p2q,2)*s13*(std::pow(s12,2)*(4*s13 - s23) +
s12*(2*s13 - s23)*s23 - 4*s13*std::pow(s23,2)) +
2*p2q*p3q*(-2*std::pow(s13,2)*std::pow(s23,2) + s12*s13*s23*(s13 + s23) +
std::pow(s12,2)*(2*std::pow(s13,2) - s13*s23 +
2*std::pow(s23,2)))))*Dot(p1,E0conj)*Dot(p2,E0))/(p1q*p2q*(p1q + p2q)*p3q*(p1q +
p3q)*(p2q + p3q)*std::pow(s12,2)*std::pow(s13,2)*std::pow(s23,2)) +
(36*(std::pow(p1q,4)*(std::pow(p3q,2)*s12*(2*s12 + s13) +
std::pow(p2q,2)*s13*(s12 + 2*s13) + 2*p2q*p3q*(std::pow(s12,2) + s12*s13 +
std::pow(s13,2)))*std::pow(s23,2)*(s12 + s13 + s23) +
std::pow(p2q,2)*std::pow(p3q,2)*(p2q + p3q)*s12*s13*(s12 + s13 +
s23)*(p2q*s13*(2*s12 + s23) + p3q*s12*(2*s13 + s23)) -
std::pow(p1q,3)*p2q*p3q*s23*(s12 + s13 + s23)*(p2q*(s12*s13*(s13 - 3*s23) -
2*std::pow(s12,2)*s23 - 2*std::pow(s13,2)*s23) + p3q*(std::pow(s12,2)*(s13 -
2*s23) - 3*s12*s13*s23 - 2*std::pow(s13,2)*s23)) +
std::pow(p1q,2)*p2q*p3q*std::pow(p2q + p3q,2)*s12*s13*(2*std::pow(s12,2)*s13 +
std::pow(s23,2)*(s13 + s23) + s12*(2*std::pow(s13,2) + 2*s13*s23 +
std::pow(s23,2))) + p1q*p2q*p3q*s12*s13*(std::pow(p2q,3)*s13*(2*s12 + s23)*(s12
+ s13 + s23) + std::pow(p3q,3)*s12*(s12 + s13 + s23)*(2*s13 + s23) +
std::pow(p2q,2)*p3q*(s23*std::pow(s13 + s23,2) + std::pow(s12,2)*(6*s13 + s23) +
2*s12*(3*std::pow(s13,2) + 4*s13*s23 + std::pow(s23,2))) +
p2q*std::pow(p3q,2)*(s23*std::pow(s13 + s23,2) + std::pow(s12,2)*(6*s13 + s23) +
2*s12*(3*std::pow(s13,2) + 4*s13*s23 +
std::pow(s23,2)))))*Dot(p2,E0)*Dot(p2,E0conj))/(p1q*p2q*(p1q + p2q)*p3q*(p1q +
p3q)*(p2q + p3q)*std::pow(s12,2)*std::pow(s13,2)*std::pow(s23,2)) +
Dot(p1,E0)*((36*(std::pow(p1q,4)*p3q*(p2q + p3q)*s12*(2*s12 +
s13)*std::pow(s23,2)*(s12 + s13 + s23) +
std::pow(p2q,4)*std::pow(p3q,2)*s12*std::pow(s13,2)*(2*s12 + s23)*(s12 + s13 +
s23) + std::pow(p1q,3)*p3q*(p2q + p3q)*s12*s23*(s12 + s13 + s23)*(p3q*s13*s23 +
p3q*s12*(s13 + 4*s23) + p2q*(std::pow(s13,2) + 2*s12*s23)) + p1q*p2q*p3q*(p2q +
p3q)*(s12 + s13 + s23)*(p2q*p3q*s12*s13*(-s12 + s13)*s23 +
std::pow(p3q,2)*std::pow(s12,2)*s23*(s13 + 2*s23) +
2*std::pow(p2q,2)*std::pow(s13,2)*(std::pow(s12,2) + s12*s23 + std::pow(s23,2)))
+ std::pow(p1q,2)*(std::pow(p2q,4)*std::pow(s13,2)*s23*(s12 + s13 + s23)*(s12 +
2*s23) + std::pow(p3q,4)*std::pow(s12,2)*s23*(s12 + s13 + s23)*(s13 + 2*s23) +
std::pow(p2q,3)*p3q*s13*(2*std::pow(s12,3)*s13 + 2*s13*std::pow(s23,2)*(s13 +
s23) + s12*s23*(3*std::pow(s13,2) + 4*s13*s23 - std::pow(s23,2)) +
std::pow(s12,2)*(2*std::pow(s13,2) + 5*s13*s23 - std::pow(s23,2))) +
2*std::pow(p2q,2)*std::pow(p3q,2)*s12*s23*(2*std::pow(s12,2)*s23 +
std::pow(s13,2)*(s13 + s23) + s12*(std::pow(s13,2) + 2*s13*s23 +
2*std::pow(s23,2))) + p2q*std::pow(p3q,3)*s12*s23*(s13*std::pow(s13 + s23,2) +
std::pow(s12,2)*(s13 + 6*s23) + 2*s12*(std::pow(s13,2) + 4*s13*s23 +
3*std::pow(s23,2)))))*Dot(p1,E0conj))/(p1q*p2q*(p1q + p2q)*p3q*(p1q + p3q)*(p2q
+ p3q)*std::pow(s12,2)*std::pow(s13,2)*std::pow(s23,2)) + (18*(s12 + s13 +
s23)*(2*std::pow(p1q,4)*p3q*(p2q + p3q)*s12*(2*s12 + s13)*std::pow(s23,2) +
std::pow(p2q,3)*std::pow(p3q,2)*s12*s13*(p3q*s13*s23 + 2*p2q*s13*(2*s12 + s23) +
p3q*s12*(4*s13 + s23)) +
p1q*std::pow(p2q,2)*p3q*s12*s13*(2*std::pow(p2q,2)*s13*(2*s12 + s23) +
std::pow(p3q,2)*(4*s12*s13 - 2*s12*s23 + 3*s13*s23 + std::pow(s23,2)) +
p2q*p3q*(8*s12*s13 - s12*s23 + 4*s13*s23 + std::pow(s23,2))) +
std::pow(p1q,3)*s23*(p2q*std::pow(p3q,2)*s12*(-(s12*s13) + std::pow(s13,2) +
8*s12*s23 + 4*s13*s23) - std::pow(p2q,2)*p3q*(std::pow(s12,2)*(s13 - 4*s23) +
s12*s13*(s13 - 2*s23) + 4*std::pow(s13,2)*s23) - std::pow(p2q,3)*s13*(4*s13*s23
+ s12*(s13 + s23)) + std::pow(p3q,3)*s12*(s13*s23 + s12*(s13 + 4*s23))) +
std::pow(p1q,2)*p2q*p3q*(std::pow(p3q,2)*s12*s23*(-2*s12*s13 + std::pow(s13,2) +
4*s12*s23 + 3*s13*s23) + std::pow(p2q,2)*s13*(std::pow(s12,2)*(4*s13 - s23) +
s12*(2*s13 - s23)*s23 - 4*s13*std::pow(s23,2)) +
2*p2q*p3q*(-2*std::pow(s13,2)*std::pow(s23,2) + s12*s13*s23*(s13 + s23) +
std::pow(s12,2)*(2*std::pow(s13,2) - s13*s23 +
2*std::pow(s23,2)))))*Dot(p2,E0conj))/(p1q*p2q*(p1q + p2q)*p3q*(p1q + p3q)*(p2q
+ p3q)*std::pow(s12,2)*std::pow(s13,2)*std::pow(s23,2)));
  }
  return output;
}

