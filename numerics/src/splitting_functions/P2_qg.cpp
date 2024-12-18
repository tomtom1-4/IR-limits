#include "../splitting_functions.hpp"

std::vector<std::vector<std::complex<double>>> P2_qg(double z, LV<double> p, LV<double> kperp, double s12, double muR, int nl) {
  double n = 3.;
  double nf = nl;
  double r2NS = li2((-1 + z)/z)*(-1 - pow(n,-2)) + li2(1 - z)*(1 + pow(n,-2)) + log(1/z)*(2 +
  2*pow(n,-2)) + li2(z/(-1 + z))*(1 + pow(n,2)) - 2*log(1/(1 - z))*(1 + pow(n,2))
  + (zeta2*pow(n,-2)*(1 + pow(n,2))*(-2 - z + 2*(1 + 7*z)*pow(n,2)))/(2.*z) +
  log((muR*muR)/s12)*(log(1/z)*(2 + 2*pow(n,-2)) - 2*log(1/(1 - z))*(1 + pow(n,2)) -
  (pow(n,-2)*(1 + pow(n,2))*(11 + 4*n*nf + 13*pow(n,2)))/6.) -
  (li2(z)*pow(n,-2)*(-1 + pow(n,4)))/(2.*z) + log(z)*((zeta2*(3 - 6/(-1 + z) - z +
  (3 - z)*pow(n,-2)))/2. - (3*pow(n,-2)*(1 + pow(n,2)))/2. + (li2(z)*pow(n,-2)*(-3
  + 4*z + pow(n,2)*(3 + 4*z - pow(z,2)) - pow(z,2)))/(2.*(-1 + z))) +
  (3*li3(z)*pow(n,-2)*(3 - 4*z + pow(z,2) + pow(n,2)*(-1 - 4*z +
  pow(z,2))))/(2.*(-1 + z)) + (li3(1 - z)*pow(n,-2)*pow(z,-2)*(-2 + 6*z + (2 +
  6*z)*pow(n,4) - 9*pow(z,2) + 3*(-3 + z)*pow(n,2)*pow(z,2) + 3*pow(z,3)))/2. +
  log(1 - z)*((log(z)*(1 + pow(n,-2)))/2. + (pow(n,-2)*(9*(-1 + z) + 9*(-1 +
  z)*pow(n,2) - 2*nf*pow(n,3) + 2*pow(n,4)))/(6.*z) - (li2(1 -
  z)*pow(n,-2)*pow(z,-2)*((1 + 3*z)*pow(n,4) + pow(-1 + z,3) + (-3 +
  z)*pow(n,2)*pow(z,2)))/2. - (zeta2*pow(n,-2)*pow(z,-2)*(1 - 3*z - (1 +
  3*z)*pow(n,4) - 3*pow(z,2) + (-3 + z)*pow(n,2)*pow(z,2) + pow(z,3)))/2.) - 2*(1
  + pow(n,2))*pow(log((muR*muR)/s12),2) + (pow(n,-2)*pow(z,-2)*(72*zeta3*(-1 + pow(n,4))
  + 144*z*zeta3*(2 + pow(n,4)) - 2*(-62*n*nf + 6*(-19 + 18*zeta3) + (11 -
  216*zeta3)*pow(n,2) - 98*nf*pow(n,3) + (161 + 108*zeta3)*pow(n,4))*pow(z,2) +
  (-219 - 124*n*nf + 40*pow(n,2) - 196*nf*pow(n,3) + 331*pow(n,4))*pow(z,3) -
  9*pow(z,4)*pow(1 + pow(n,2),2)))/(72.*(-1 + z));

  double r2S = li2((-1 + z)/z)*(-7.444444444444445 + (10*nf)/(9.*n) - 12*zeta2 - 2*li2(z/(-1 +
  z))) + log(1/z)*((28*zeta3)/3. + 2*li3(z/(-1 + z)) - 2*li3((-1 + z)/z)*pow(n,-2)
  + (25*zeta2*pow(n,-2)*(11 + 4*n*nf - 11*pow(n,2)))/24.) + (li3((-1 +
  z)/z)*pow(n,-2)*(11 + 4*n*nf - 11*pow(n,2)))/3. + li4(1 - z)*(-8 - 2*pow(n,-2) -
  6*pow(n,2)) - 2*li4(z/(-1 + z))*pow(n,2) + li4((-1 + z)/z)*(6 - 2*pow(n,-2) +
  6*pow(n,2)) + (li3(z/(-1 + z))*(-11 - 4*n*nf + 11*pow(n,2)))/3. + (li3(z)*(18 -
  4*n*nf + 9*pow(n,-2) + 13*pow(n,2)))/6. + li2(z/(-1 + z))*((n*(67*n - 10*nf))/9.
  + 12*zeta2*pow(n,2)) + log(1/(1 - z))*(2*li3((-1 + z)/z) + 14*zeta2*log(1/z) -
  (28*zeta3*pow(n,2))/3. - 2*li3(z/(-1 + z))*pow(n,2) + (25*zeta2*(-11 - 4*n*nf +
  11*pow(n,2)))/24.) + log(-((-1 + z)*z))*((zeta2*(22 + 10*n*nf -
  33*pow(n,2)))/12. + (-737*n + 110*nf + 34*nf*pow(n,2) + (35 -
  108*zeta3)*pow(n,3))/(108.*n)) + log((muR*muR)/s12)*(4*li3((-1 + z)/z) +
  log(1/z)*(7.444444444444445 - (10*nf)/(9.*n) + 26*zeta2 + 4*li2(z/(-1 + z)) -
  4*li2((-1 + z)/z)*pow(n,-2)) + (li2((-1 + z)/z)*pow(n,-2)*(11 + 4*n*nf -
  11*pow(n,2)))/3. - 4*li3(z/(-1 + z))*pow(n,2) + (li2(z/(-1 + z))*(-11 - 4*n*nf +
  11*pow(n,2)))/3. + zeta2*(-13.291666666666666 - 5*n*nf + (341*pow(n,2))/24.) +
  log(1/(1 - z))*((n*(-67*n + 10*nf))/9. + 4*li2((-1 + z)/z) - 26*zeta2*pow(n,2) -
  4*li2(z/(-1 + z))*pow(n,2)) - (-737*n + 110*nf + 34*nf*pow(n,2) + 5*(7 +
  180*zeta3)*pow(n,3))/(108.*n)) + li2(1 - z)*(3*zeta2*(1 + pow(n,2)) -
  (pow(n,-2)*(-1 + pow(n,4)))/2.) - ((-1 + z)*li2(z)*pow(n,-2)*(-1 +
  pow(n,4)))/(2.*z) + (zeta2*pow(n,-2)*(-12*(-1 + z) - 50*nf*z*pow(n,3) + (-12 +
  347*z)*pow(n,4)))/(12.*z) + log(z)*(li3(1 - z)*(1 + pow(n,-2)) + li3(z)*(2 +
  2*pow(n,-2)) - (zeta2*pow(n,-2)*(36 + 2*n*nf + 25*pow(n,2) + 10*nf*pow(n,3) -
  55*pow(n,4)))/12. + (pow(n,-2)*(56*n*nf + 27*(3 + 4*zeta3) + (-323 +
  162*zeta3)*pow(n,2) + 38*nf*pow(n,3) + (-386 + 54*zeta3)*pow(n,4)))/54.) + li3(1
  - z)*(6 - (4*n*nf)/3. + pow(n,2)*(4.333333333333333 - 2/z - pow(z,-2)) +
  pow(n,-2)*(3 - 4/z + pow(z,-2))) + (pow(n,-2)*pow(z,-2)*(3240*zeta3*(-1 +
  pow(n,4)) + 6480*z*zeta3*(2 + pow(n,4)) + (-100*nf*(-13 + 90*zeta3)*pow(n,3) +
  36*(-45*(4 + 9*zeta3) + 2*pow(M_PI,4)) + 9*pow(n,2)*(-5*(504 + 293*zeta3) +
  28*pow(M_PI,4)) + pow(n,4)*(-39040 + 41175*zeta3 +
  999*pow(M_PI,4)))*pow(z,2)))/3240. + ((-1 - 2*pow(n,-2) + pow(n,2))*pow(li2(1 -
  z),2))/2. + pow(n,-2)*pow(li2((-1 + z)/z),2) + pow(n,2)*pow(li2(z/(-1 + z)),2) +
  ((-8*log(1/z))/3. + (11 + 4*n*nf - 11*pow(n,2))/18. + (8*log(1/(1 -
  z))*pow(n,2))/3.)*pow(log((muR*muR)/s12),3) + (2*pow(n,2)*pow(log((muR*muR)/s12),4))/3. -
  7*zeta2*pow(n,2)*pow(log(1/(1 - z)),2) - 7*zeta2*pow(n,-2)*pow(log(1/z),2) +
  pow(log((muR*muR)/s12),2)*((n*(-67*n + 10*nf))/18. + 4*li2((-1 + z)/z) + log(1/(1 -
  z))*(-4*log(1/z) + (11 + 4*n*nf - 11*pow(n,2))/6.) - (log(1/z)*pow(n,-2)*(11 +
  4*n*nf - 11*pow(n,2)))/6. - 13*zeta2*pow(n,2) - 4*li2(z/(-1 + z))*pow(n,2) +
  2*pow(n,2)*pow(log(1/(1 - z)),2) + 2*pow(n,-2)*pow(log(1/z),2)) + zeta2*(2 -
  pow(n,-2) + 3*pow(n,2))*pow(log(z),2) + (1 - pow(n,-2) + 2*pow(n,2))*pow(log(1 -
  z),2)*pow(log(z),2) + log(1 - z)*(4*li3(1 - z)*(1 + pow(n,2)) + li3(z)*(2 -
  2*pow(n,-2) + 4*pow(n,2)) + log(z)*(zeta2*(-1 + 2*pow(n,-2) - 3*pow(n,2)) +
  li2(1 - z)*(-1 - 2*pow(n,-2) + pow(n,2)) - (pow(n,-2)*(-1 + pow(n,4)))/2.) +
  (pow(n,-2)*(9 + 3*z*(-3 + 4*zeta3) - 3*(-3 + z*(3 + 4*zeta3))*pow(n,2) -
  2*nf*(-1 + z)*pow(n,3) + (-2 + z*(2 - 24*zeta3))*pow(n,4)))/(6.*z) + ((-1 +
  z)*zeta2*pow(n,-2)*(-1 + 3*z + (1 + 3*z)*pow(n,4))*pow(z,-2))/2. + (li2(1 -
  z)*pow(n,-2)*pow(z,-2)*(-3 + 12*z + pow(n,4)*(3 + 6*z - 13*pow(z,2)) -
  9*pow(z,2) - 18*pow(n,2)*pow(z,2) + 4*nf*pow(n,3)*pow(z,2)))/6. + ((-2 +
  pow(n,-2) - 3*pow(n,2))*pow(log(z),3))/3.) + ((2 - pow(n,-2) +
  3*pow(n,2))*pow(log(z),4))/12.;

  std::vector<std::vector<std::complex<double>>> output = {{0.,0.},{0.,0.}};
  std::vector<std::vector<std::complex<double>>> P0qg = P0_qg(z, p, kperp);
  for(int i = 0; i <= 1; i++) for(int j = 0; j <= 1; j++) {
    output[i][j] = (r2S*P0qg[i][j] + C_F*r2NS*((i==j)?1.:0.));
  }
  return output;
}