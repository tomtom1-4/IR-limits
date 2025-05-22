#include "../splitting_functions.hpp"

// Based on Badger and Glover, 2004
// hep-ph/0405236

std::vector<std::vector<std::complex<double>>> P2_qg(double z, LV<double> p, LV<double> kperp, double s12, double muR, int nl) {
  double n = 3.;
  double nf = nl;
  double mu2 = muR*muR;
  double Pi = M_PI;
  std::cout << "s12 = " << s12 << ", mu2 = " << mu2 << ", z= " << z << std::endl;
  double r2NS = li2((-1 + z)/z)*(-1 - pow(n,-2)) + li2(1 - z)*(1 + pow(n,-2)) + log(1/z)*(1 +
  pow(n,-2)) + log(1/(1 - z))*(-1 - pow(n,2)) + li2(z/(-1 + z))*(1 + pow(n,2)) +
  (zeta2*pow(n,-2)*(1 + pow(n,2))*(-2*(2 + z) + (4 + 27*z)*pow(n,2)))/(4.*z) +
  log(mu2/s12)*((3*log(1/z)*(1 + pow(n,-2)))/2. + ((5*n - 2*nf)*(1 +
  pow(n,2)))/(3.*n) - (3*log(1/(1 - z))*(1 + pow(n,2)))/2.) -
  (li2(z)*pow(n,-2)*(-1 + pow(n,4)))/(2.*z) + (3*li3(z)*pow(n,-2)*(3 - 4*z +
  pow(z,2) + pow(n,2)*(-1 - 4*z + pow(z,2))))/(2.*(-1 + z)) +
  log(z)*((-3*pow(n,-2)*(1 + pow(n,2)))/2. + (li2(z)*pow(n,-2)*(-3 + 4*z +
  pow(n,2)*(3 + 4*z - pow(z,2)) - pow(z,2)))/(2.*(-1 + z)) - (zeta2*pow(n,-2)*(3 -
  4*z + pow(z,2) + pow(n,2)*(9 - 4*z + pow(z,2))))/(2.*(-1 + z))) + (li3(1 -
  z)*pow(n,-2)*pow(z,-2)*(-2 + 6*z + (2 + 6*z)*pow(n,4) - 9*pow(z,2) + 3*(-3 +
  z)*pow(n,2)*pow(z,2) + 3*pow(z,3)))/2. + log(1 - z)*((log(z)*(1 + pow(n,-2)))/2.
  + (pow(n,-2)*(9*(-1 + z) + 9*(-1 + z)*pow(n,2) - 2*nf*pow(n,3) +
  2*pow(n,4)))/(6.*z) - (li2(1 - z)*pow(n,-2)*pow(z,-2)*((1 + 3*z)*pow(n,4) +
  pow(-1 + z,3) + (-3 + z)*pow(n,2)*pow(z,2)))/2. - (zeta2*pow(n,-2)*pow(z,-2)*(1
  - 3*z - (1 + 3*z)*pow(n,4) - 3*pow(z,2) + (-3 + z)*pow(n,2)*pow(z,2) +
  pow(z,3)))/2.) - (5*(1 + pow(n,2))*pow(log(mu2/s12),2))/4. +
  (pow(n,-2)*pow(z,-2)*(72*zeta3*(-1 + pow(n,4)) + 144*z*zeta3*(2 + pow(n,4)) -
  2*(-62*n*nf + 18*(1 + 6*zeta3) + (347 - 216*zeta3)*pow(n,2) - 98*nf*pow(n,3) +
  (365 + 108*zeta3)*pow(n,4))*pow(z,2) + (45 - 124*n*nf + 712*pow(n,2) -
  196*nf*pow(n,3) + 739*pow(n,4))*pow(z,3) - 9*pow(z,4)*pow(1 +
  pow(n,2),2)))/(72.*(-1 + z));

  double r2S = -3*zeta3 - (5*n*nf*(-13 + 90*zeta3))/162. + li2((-1 + z)/z)*(-7.444444444444445
  + (10*nf)/(9.*n) - (23*zeta2)/2. - 2*li2(z/(-1 + z))) + (-7.333333333333333 +
  (4*nf)/(3.*n))*li3((-1 + z)/z) + (2*n*(11*n - 2*nf)*li3(z/(-1 + z)))/3. +
  ((5*n*(-11*n + 2*nf)*zeta2)/12. - (n*(19*nf + n*(-193 +
  27*zeta3)))/27.)*log(-((-1 + z)*z)) + log(1/z)*((25*(-11 + (2*nf)/n)*zeta2)/12.
  + (14*zeta3)/3. + li3(z/(-1 + z)) - li3((-1 + z)/z)*pow(n,-2)) + li4(1 - z)*(-8
  - 2*pow(n,-2) - 6*pow(n,2)) + li4(z)*(-7 - 2*pow(n,-2) - 5*pow(n,2)) - li4(z/(-1
  + z))*pow(n,2) + li4((-1 + z)/z)*(5 - 2*pow(n,-2) + 6*pow(n,2)) + (li3(z)*(18 -
  4*n*nf + 9*pow(n,-2) + 13*pow(n,2)))/6. + li2(z/(-1 + z))*((n*(67*n - 10*nf))/9.
  + (23*zeta2*pow(n,2))/2.) + log(1/(1 - z))*((25*n*(11*n - 2*nf)*zeta2)/12. +
  li3((-1 + z)/z) + 13*zeta2*log(1/z) - (14*zeta3*pow(n,2))/3. - li3(z/(-1 +
  z))*pow(n,2)) + log(mu2/s12)*((5*n*(11*n - 2*nf)*zeta2)/2. - (n*(-19*nf + n*(193
  + 99*zeta3)))/27. + (-7.333333333333333 + (4*nf)/(3.*n))*li2((-1 + z)/z) +
  (2*n*(11*n - 2*nf)*li2(z/(-1 + z)))/3. + 2*li3((-1 + z)/z) +
  log(1/z)*(7.444444444444445 - (10*nf)/(9.*n) + (49*zeta2)/2. + 3*li2(z/(-1 + z))
  - 3*li2((-1 + z)/z)*pow(n,-2)) - 2*li3(z/(-1 + z))*pow(n,2) + log(1/(1 -
  z))*((n*(-67*n + 10*nf))/9. + 3*li2((-1 + z)/z) - (49*zeta2*pow(n,2))/2. -
  3*li2(z/(-1 + z))*pow(n,2))) + li2(1 - z)*(3*zeta2*(1 + pow(n,2)) -
  (pow(n,-2)*(-1 + pow(n,4)))/2.) - ((-1 + z)*li2(z)*pow(n,-2)*(-1 +
  pow(n,4)))/(2.*z) + (zeta2*pow(n,-2)*(-12*(-1 + z) - 50*nf*z*pow(n,3) + (-12 +
  347*z)*pow(n,4)))/(12.*z) + log(z)*(li3(1 - z)*(1 + pow(n,-2)) + li3(z)*(2 +
  2*pow(n,-2)) - (zeta2*pow(n,-2)*(36 + 2*n*nf + 25*pow(n,2) + 10*nf*pow(n,3) -
  55*pow(n,4)))/12. + (pow(n,-2)*(56*n*nf + 27*(3 + 4*zeta3) + (-323 +
  162*zeta3)*pow(n,2) + 38*nf*pow(n,3) + (-386 + 54*zeta3)*pow(n,4)))/54.) +
  (7*pow(Pi,4))/90. + li3(1 - z)*(6 - (4*n*nf)/3. + pow(n,2)*(4.333333333333333 -
  2/z - pow(z,-2)) + pow(n,-2)*(3 - 4/z + pow(z,-2))) +
  pow(n,2)*(-7.049382716049383 + (27*pow(Pi,4))/80. + zeta3*(13.777777777777779 +
  2/z + pow(z,-2))) + pow(n,-2)*(pow(Pi,4)/45. - (zeta3*pow(z,-2)*(2 - 8*z +
  9*pow(z,2)))/2.) + (pow(n,2)*pow(zeta2,2))/8. + ((-1 - 2*pow(n,-2) +
  pow(n,2))*pow(li2(1 - z),2))/2. + pow(n,-2)*pow(li2((-1 + z)/z),2) +
  pow(n,2)*pow(li2(z/(-1 + z)),2) + ((n*(-11*n + 2*nf))/9. - (11*log(1/z))/6. +
  (11*log(1/(1 - z))*pow(n,2))/6.)*pow(log(mu2/s12),3) +
  (11*pow(n,2)*pow(log(mu2/s12),4))/24. - (13*zeta2*pow(n,2)*pow(log(1/(1 -
  z)),2))/2. - (13*zeta2*pow(n,-2)*pow(log(1/z),2))/2. +
  pow(log(mu2/s12),2)*((n*(-67*n + 10*nf))/18. + (5*li2((-1 + z)/z))/2. + log(1/(1
  - z))*((n*(-11*n + 2*nf))/3. - 3*log(1/z)) + (3.6666666666666665 -
  (2*nf)/(3.*n))*log(1/z) - (49*zeta2*pow(n,2))/4. - (5*li2(z/(-1 +
  z))*pow(n,2))/2. + (3*pow(n,2)*pow(log(1/(1 - z)),2))/2. +
  (3*pow(n,-2)*pow(log(1/z),2))/2.) + zeta2*(2 - pow(n,-2) +
  3*pow(n,2))*pow(log(z),2) + (1 - pow(n,-2) + 2*pow(n,2))*pow(log(1 -
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

  std::cout << "r2S = " << r2S << std::endl;
  std::cout << "r2NS = " << r2NS << std::endl;
  std::vector<std::vector<std::complex<double>>> output = {{0.,0.},{0.,0.}};
  std::vector<std::vector<std::complex<double>>> P0qg = P0_qg(z, p, kperp);
  for(int i = 0; i <= 1; i++) for(int j = 0; j <= 1; j++) {
    output[i][j] = (r2S*P0qg[i][j] + C_F*r2NS*((i==j)?1.:0.));
  }
  return output;
}