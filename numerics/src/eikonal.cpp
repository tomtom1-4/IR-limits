#include "eikonal.hpp"

const unsigned nLi2dat = 7;
const double li2dat[] = {
  +1.0000000000000000e0,
  +2.7777777777777418e-2,
  -2.7777777776967998e-4,
  +4.7241117865795084e-6,
  -9.1857320746997598e-8,
  +1.8967653373154463e-9,
  -3.9078560491306871e-11,
};

const double zeta2 = M_PI*M_PI/6.;

double li2(double x)
{
  if (x > 1.) std::cout << "out-of-range argument" << std::endl;
  if (x == 1.) return zeta2;
  if (x > 0.5) return -li2(1.-x)-log(x)*log(1.-x)+zeta2;
  if (x < -1.) {
    const double ln = log(-x);
    return -li2(1./x)-0.5*ln*ln-zeta2;
  }

  double result = 0.0;
  const double y = -log(1.-x);
  const double y2 = y*y;
  for (int i = nLi2dat-1; i > 0; --i)
    {
      result += li2dat[i];
      result *= y2;
    }
  result += li2dat[0]-0.25*y;
  result *= y;

  return result;
}

template <typename Field> struct LV {
  std::vector<Field> components = std::vector<Field>(4);
  LV(){
    components = std::vector<Field>(4);
  }
  LV(std::vector<Field> lv):components(lv){};
  LV(Field* lv) {
    for(int i = 0; i < 4; i++) components[i] = lv[i];
  }
  void print() const {
    std::cout << "{" << components[0] << ", "
                     << components[1] << ", "
                     << components[2] << ", "
                     << components[3] << "}" << std::endl;
  }

  // Define the multiplication operator for Field1 * LV<Field2>
  template <typename Field1>
  LV<typename std::common_type<Field1, Field>::type> operator*(const Field1& scalar) const {
    LV<typename std::common_type<Field1, Field>::type> result;
    for (int i = 0; i < 4; ++i) {
      result.components[i] = scalar*components[i];
    }
    return result;
  }
};

template <typename Field> struct LM {
  std::vector<std::vector<Field>> components = std::vector<std::vector<Field>>(4, std::vector<Field>(4));
  LM(){
    components = std::vector<std::vector<Field>>(4, std::vector<Field>(4));
  }
  LM(std::vector<std::vector<Field>> lm):components(lm){};
  void print() const {
    std::cout << "{";
    for(int i = 0; i < 4; i++) {
      std::cout << "{";
      for(int j = 0; j < 4; j++) {
        std::cout << components[i][j] << ",";
      }
      std::cout << "}," << std::endl;
    }
    std::cout << std::endl;
  }

  // Define the multiplication operator for Field1 * LM<Field2>
  template <typename Field1>
  LM<typename std::common_type<Field1, Field>::type> operator*(const Field1& scalar) const {
    LM<typename std::common_type<Field1, Field>::type> result;
    for (int i = 0; i < 4; ++i) for(int j = 0; j < 4; j++) {
      result.components[i][j] *= scalar;
    }
    return result;
  }

  // Define the multiplication operator for LM<Field1> * LM<Field2>
  template <typename Field1>
  LM<typename std::common_type<Field1, Field>::type> operator*(const LM<Field1>& m) const {
    LM<typename std::common_type<Field1, Field>::type> result;
    for (int i = 0; i < 4; ++i) for(int j = 0; j < 4; j++) for(int k = 0; k < 4; k++) {
      result.components[i][j] += m.components[i][k]*this->components[k][j]*(k==0?1.:-1.);
    }
    return result;
  }

  LM<Field> transpose() {
    LM<Field> output;
    for(int mu1 = 0; mu1 < 4; mu1++) for(int mu2 = 0; mu2 < 4; mu2++) {
      output.components[mu1][mu2] = this->components[mu2][mu1];
    }
    return output;
  }

  Field trace() {
    Field output = 0;
    for(int i = 0; i < 4; i++) {
      output += (i==0?1.:-1.)*components[i][i];
    }
    return output;
  }
};

template <typename Field1, typename Field2>
LV<typename std::common_type<Field1, Field2>::type> operator*(const Field1& scalar, const LV<Field2>& v) {
  LV<typename std::common_type<Field1, Field2>::type> output = v*scalar;
  return output;
}

template <typename Field1, typename Field2>
LV<typename std::common_type<Field1, Field2>::type> operator/(const LV<Field1>& v, const Field2& scalar) {
  LV<typename std::common_type<Field1, Field2>::type> result;
  for (size_t i = 0; i < 4; ++i) {
    result.components[i] = v.components[i]/scalar;
  }
  return result;
}

template <class Field1, class Field2>
LV<typename std::common_type<Field1, Field2>::type> operator+(const LV<Field1>& v1, const LV<Field2>& v2) {
  LV<typename std::common_type<Field1, Field2>::type> output;
  for(int i = 0; i < 4; i++) {
    output.components[i] = v1.components[i] + v2.components[i];
  }
  return output;
}

template <class Field1, class Field2>
LV<typename std::common_type<Field1, Field2>::type> operator-(const LV<Field1>& v1, const LV<Field2>& v2) {
  LV<typename std::common_type<Field1, Field2>::type> output;
  for(int i = 0; i < 4; i++) {
    output.components[i] = v1.components[i] - v2.components[i];
  }
  return output;
}

template <class Field1, class Field2>
typename std::common_type<Field1, Field2>::type operator*(const LV<Field1>& v1, const LV<Field2>& v2) {
  typename std::common_type<Field1, Field2>::type output;
  for(int i = 0; i < 4; i++) {
    output += (i==0?1.:-1.)*v1.components[i]*v2.components[i];
  }
  return output;
}

template <class Field1, class Field2>
LV<typename std::common_type<Field1, Field2>::type> operator*(const LM<Field1>& m, const LV<Field2>& v) {
  LV<typename std::common_type<Field1, Field2>::type> output;
  for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) {
    output.components[i] += m.components[i][j]*v.components[j]*(j==0?1.:-1.);
  }
  return output;
}

LV<double> j1(const LV<double>& p, const LV<double>& q) {
  return p/(p*q);
}

LV<std::complex<double>> gamma11(const LV<double>& p1, const LV<double>& p2, const LV<double>& q) {
  if(p1*p2<1.e-8) return LV<std::complex<double>>({0.,0.,0.,0.});
  std::complex<double> log = std::log((p1*p2)*mu*mu/(p1*q)/(p2*q)/2.);
  if((p1.components[0] < 0) and (p2.components[0] < 0)) log = log + I*M_PI;
  else log = log - I*M_PI;
  std::complex<double> prefactor = -1./12.*(std::pow(M_PI, 2) + 6.*std::pow(log, 2))/std::pow(4.*M_PI, 2);
  LV<std::complex<double>> output = prefactor*(p1/(p1*q) - p2/(p2*q));

  return output;
}

LM<double> gamma20(const LV<double>& p, const LV<double>& q1, const LV<double>& q2) {
  LM<double> output;
  for(int mu1 = 0; mu1 < 4; mu1++) for(int mu2 = 0; mu2 < 4; mu2++) {
    output.components[mu1][mu2] = 1./(p*(q1 + q2))*(p.components[mu1]*p.components[mu2]/2.*(1./(p*q1) - 1./(p*q2))
        + 1./(q1*q2)*(p.components[mu1]*q1.components[mu2] - p.components[mu2]*q2.components[mu1] + 1./2.*(p*q2 - p*q1)*metric[mu1][mu2]));
  }
  return output;
}

LM<std::complex<double>> gamma21(const LV<double>& pi, const LV<double>& pj, const LV<double>& q1, const LV<double>& q2, int nl=5) {
  LV<double> pi_perp = pi - (pi*q2)/(q1*q2)*q1 - (pi*q1)/(q1*q2)*q2;
  LV<double> pj_perp = pj - (pj*q2)/(q1*q2)*q1 - (pj*q1)/(q1*q2)*q2;
  LM<double> g_perp;
  for(int mu1 = 0; mu1 < 4; mu1++) for(int mu2 = 0; mu2 < 4; mu2++) {
    g_perp.components[mu1][mu2] = (metric[mu1][mu2] - (q1.components[mu1]*q2.components[mu2] + q1.components[mu2]*q2.components[mu1])/(q1*q2));
  }

  LM<std::complex<double>> output;
  bool analytic_continuation = false;
  if((pi.components[0] < 0.) and (pj.components[0] < 0.)) analytic_continuation = true;
  std::complex<double> logc = std::log(2.*q1*q2/mu/mu) - I*M_PI;
  double l1 = std::log((pi*q1)/(pi*q1 + pi*q2));
  double l2 = std::log((pj*q2)/(pj*q1 + pj*q2));
  double L1 = li2(1. - pj*q2/(pj*q1 + pj*q2));
  double L2 = li2(1. - pi*q1/(pi*q1 + pi*q2));
  std::complex<double> L3 = li2(1. - (pi*pj)*(q1*q2)/(pi*q1 + pi*q2)/(pj*q1 + pj*q2));
  if(analytic_continuation) L3 = L3 + 2.*M_PI*I*std::log(1. - (pi*pj)*(q1*q2)/(pi*q1 + pi*q2)/(pj*q1 + pj*q2));

  double s1j = 2.*(pj*q1);
  double s1i = 2.*(pi*q1);
  double s2j = 2.*(pj*q2);
  double s2i = 2.*(pi*q2);
  double sij = 2.*(pi*pj);
  double s12 = 2.*(q1*q2);

  std::complex<double> l3 = std::log(sij*s12/s1i/s2j);
  if(analytic_continuation) l3 = l3 - 2.*M_PI*I;
  std::vector<std::complex<double>> JHat(3);
  JHat[0] = (-2*(-(s1j*s2i) + s1i*s2j + s12*sij))/(s1i*s2j);
  JHat[1] = -((l3*(-(s1j*s2i) + s1i*s2j + s12*sij))/(s1i*s2j));
  JHat[2] = (-3.*(2*std::pow(l2,2) + 2*l2*l3 + std::pow(l3,2) + 2.*(L2 + L3)) - (6*std::pow(l1,2)*(-(s1j*s2i) + s1i*s2j + s12*sij))/(s1i*s2j) - (6*l1*(l2 + l3)*(-(s1j*s2i) + s1i*s2j + s12*sij))/(s1i*s2j) +
     (-4*C_A*s1i*s1j*s2i*s2j + 4*C_A*std::pow(s1i,2)*std::pow(s2j,2) + 6*C_A*L1*s1j*s2i*(s1i + s2i)*(s1j + s2j) + 6*C_A*std::pow(l2,2)*s1j*s2i*(s1i + s2i)*(s1j + s2j) + 6*C_A*L2*s1j*s2i*(s1i + s2i)*(s1j + s2j) +
        6*C_A*l2*l3*s1j*s2i*(s1i + s2i)*(s1j + s2j) + 3*C_A*std::pow(l3,2)*s1j*s2i*(s1i + s2i)*(s1j + s2j) + 6*C_A*L3*s1j*s2i*(s1i + s2i)*(s1j + s2j) - 6*C_A*L1*s1i*(s1i + s2i)*s2j*(s1j + s2j) - 6*C_A*L1*s12*(s1i + s2i)*(s1j + s2j)*sij -
        6*C_A*std::pow(l2,2)*s12*(s1i + s2i)*(s1j + s2j)*sij - 6*C_A*L2*s12*(s1i + s2i)*(s1j + s2j)*sij - 6*C_A*l2*l3*s12*(s1i + s2i)*(s1j + s2j)*sij - 3*C_A*std::pow(l3,2)*s12*(s1i + s2i)*(s1j + s2j)*sij -
        6*C_A*L3*s12*(s1i + s2i)*(s1j + s2j)*sij + 8*nl*s1i*s2j*(s1j*s2i - s1i*s2j)*T_F)/(C_A*s1i*(s1i + s2i)*s2j*(s1j + s2j)))/6.;

  std::vector<std::complex<double>> Jpm(3);
  Jpm[0] = (-4*(s1i*s1j*s2i - s1i*s2i*s2j + s12*s2i*sij))/((s1i + s2i)*(s1j*s2i + s1i*s2j - s12*sij));
  Jpm[1] = (-2.*(l3*s1i*s1j*s2i - l3*s1i*s2i*s2j + l3*s12*s2i*sij))/((s1i + s2i)*(s1j*s2i + s1i*s2j - s12*sij));
  Jpm[2] = (-((6*std::pow(l1,2) + 6*l1*l3 + 3.*std::pow(l3,2) + 2*std::pow(M_PI,2))*s12*s2i*sij*(s2i*(s1j + s2j) - s12*sij)) -
     std::pow(s1i,2)*(s1j + s2j)*((6*std::pow(l1,2) + 6*(l1 + l2)*l3 + 3.*std::pow(l3,2) + 2*(3*std::pow(l2,2) + std::pow(M_PI,2)))*s2i*(s1j - s2j) - 6*l1*(l1 + l3)*s12*sij) -
     s1i*((6*std::pow(l1,2) + 6*(l1 + l2)*l3 + 3.*std::pow(l3,2) + 2*(3*std::pow(l2,2) + std::pow(M_PI,2)))*std::pow(s2i,2)*(s1j - s2j)*(s1j + s2j) +
        2*s12*s2i*(-6*l1*(l1 + l3)*s1j + (6*std::pow(l1,2) + 6*l1*l3 + 3.*std::pow(l3,2) + 2*std::pow(M_PI,2))*s2j)*sij + 6*l1*(l1 + l3)*std::pow(s12,2)*std::pow(sij,2)))/
   (3.*(s1i + s2i)*(s1j*s2i + s1i*s2j - s12*sij)*((s1i + s2i)*(s1j + s2j) - s12*sij));

  std::vector<std::complex<double>> JpmSwapped(3);
  JpmSwapped[0] = (-4*(-(s1i*s1j*s2j) + s1j*s2i*s2j + s12*s1j*sij))/((s1j + s2j)*(s1j*s2i + s1i*s2j - s12*sij));
  JpmSwapped[1] = (-2.*(-(l3*s1i*s1j*s2j) + l3*s1j*s2i*s2j + l3*s12*s1j*sij))/((s1j + s2j)*(s1j*s2i + s1i*s2j - s12*sij));
  JpmSwapped[2] = (-((6*std::pow(l2,2) + 6*l2*l3 + 3.*std::pow(l3,2) + 2*std::pow(M_PI,2))*s12*s1j*sij*(s1j*(s1i + s2i) - s12*sij)) -
     (s1i + s2i)*std::pow(s2j,2)*((6*std::pow(l2,2) + 6*(l1 + l2)*l3 + 3.*std::pow(l3,2) + 2*(3*std::pow(l1,2) + std::pow(M_PI,2)))*s1j*(-s1i + s2i) - 6*l2*(l2 + l3)*s12*sij) -
     s2j*((6*std::pow(l2,2) + 6*(l1 + l2)*l3 + 3.*std::pow(l3,2) + 2*(3*std::pow(l1,2) + std::pow(M_PI,2)))*std::pow(s1j,2)*(-s1i + s2i)*(s1i + s2i) +
        2*s12*s1j*((6*std::pow(l2,2) + 6*l2*l3 + 3.*std::pow(l3,2) + 2*std::pow(M_PI,2))*s1i - 6*l2*(l2 + l3)*s2i)*sij + 6*l2*(l2 + l3)*std::pow(s12,2)*std::pow(sij,2)))/
   (3.*(s1j + s2j)*(s1j*s2i + s1i*s2j - s12*sij)*((s1i + s2i)*(s1j + s2j) - s12*sij));

  std::vector<std::complex<double>> Jpp(3);
  Jpp[0] = -2.;
  Jpp[1] = -l3;
  Jpp[2] = (-2*(std::pow(l1,2) + l1*l2 + std::pow(l2,2)) - 2*(l1 + l2)*l3 - std::pow(l3,2) - 2.*(L1 + L2 + L3))/2.;

  std::vector<std::complex<double>> prefactor(3);
  prefactor[0] = 1.;
  prefactor[1] = -logc;
  prefactor[2] = (6.*std::pow(logc, 2) - std::pow(M_PI, 2))/12.;

  std::complex<double> JHat_ep0 = JHat[0]*prefactor[2] + JHat[1]*prefactor[1] + JHat[2]*prefactor[0];
  std::complex<double> Jpm_ep0 = Jpm[0]*prefactor[2] + Jpm[1]*prefactor[1] + Jpm[2]*prefactor[0];
  std::complex<double> JpmSwapped_ep0 = JpmSwapped[0]*prefactor[2] + JpmSwapped[1]*prefactor[1] + JpmSwapped[2]*prefactor[0];
  std::complex<double> Jpp_ep0 = Jpp[0]*prefactor[2] + Jpp[1]*prefactor[1] + Jpp[2]*prefactor[0];

  for(int mu1 = 0; mu1 < 4; mu1++) for(int mu2 = 0; mu2 < 4; mu2++) {
    output.components[mu1][mu2] = (JHat_ep0*g_perp.components[mu1][mu2]
                                + Jpp_ep0*(q1*q2)/(pi*q1)/(pj*q2)*(pi_perp.components[mu1]*pj_perp.components[mu2] - pj_perp.components[mu1]*pi_perp.components[mu2])
                                + Jpm_ep0*(g_perp.components[mu1][mu2] - 2.*pi_perp.components[mu1]*pi_perp.components[mu2]/(pi_perp*pi_perp))
                                + JpmSwapped_ep0*(g_perp.components[mu1][mu2] - 2.*pj_perp.components[mu1]*pj_perp.components[mu2]/(pj_perp*pj_perp)))*2./s12/std::pow(4.*M_PI, 2);

  }

  return output;
}

std::unordered_map<std::string, std::complex<double>> J_g_eikonal(double *pp, double *q, amplitude& A) {
  std::unordered_map<std::string, std::complex<double>> J;
  for (int i = 0; i < A.process.size(); i++) { // particle index
    double pi[4];
    part(pp, pi, 4 * i, 4 * i + 4);
    for (int mu = 0; mu < 4; mu++) {  // Lorentz index
      for (int a = 0; a < 8; a++) {  // Color index
        for (int b = 0; b < A.process[i]; b++) {
          for (int c = 0; c < A.process[i]; c++) {
            std::string key = std::to_string(i) + std::to_string(mu) + std::to_string(a) + std::to_string(b) + std::to_string(c);
            if (A.particle_type[i] == 1) {
              J[key] = - gs * lam[a][3 * b + c]/2. * pi[mu]/minkovski(pi, q);
            }
            else if (A.particle_type[i] == -1) {
              J[key] = gs * lam[a][3 * c + b]/2. * pi[mu]/minkovski(pi, q);
            }
            else if (A.particle_type[i] == 2) {
              J[key] = -gs * fabc[a][8 * c + b] * I * pi[mu]/minkovski(pi, q);
            }
          }
        }

      }
    }
  }
  return J;
}

std::unordered_map<std::string, std::complex<double>> J_gg_eikonal(double *pp, double *q1, double *q2, amplitude& A) {
  // Irreducible component of the double soft current
  std::unordered_map<std::string, std::complex<double>> J;
  for (int i = 0; i < A.process.size(); i++) { // particle index
    double pi[4];
    part(pp, pi, 4 * i, 4 * i + 4);
    for (int mu1 = 0; mu1 < 4; mu1++) for (int mu2 = 0; mu2 < 4; mu2++) {  // Lorentz index
      // momentum dependence
      double gamma = (pi[mu1]*q1[mu2] - pi[mu2]*q2[mu1])/minkovski(q1,q2)/(minkovski(pi, q1) + minkovski(pi, q2))
                    - (minkovski(pi, q1) - minkovski(pi, q2))/2./(minkovski(pi, q1) + minkovski(pi, q2))*(pi[mu1]*pi[mu2]/minkovski(pi, q1)/minkovski(pi, q2) + metric[mu1][mu2]/minkovski(q1,q2));
      for (int a1 = 0; a1 < 8; a1++) for(int a2 = 0; a2 < 8; a2++) {  // Color index
        for (int b = 0; b < A.process[i]; b++) {
          for (int c = 0; c < A.process[i]; c++) {
            std::string key = std::to_string(i) + std::to_string(mu1) + std::to_string(mu2) + std::to_string(a1) + std::to_string(a2) + std::to_string(b) + std::to_string(c);
            J[key] = 0;
            for(int d = 0; d < 8; d++) {
              if (A.particle_type[i] == 1) {
                J[key] += gs*gs*I*fabc[a1][8*a2+d] * lam[d][3 * b + c]/2.*gamma;
              }
              else if (A.particle_type[i] == -1) {
                J[key] += -gs*gs*I*fabc[a1][8*a2+d] * lam[d][3 * c + b]/2.*gamma;
              }
              else if (A.particle_type[i] == 2) {
                J[key] += gs*gs*I*fabc[a1][8*a2+d] * I*fabc[d][8 * c + b]*gamma;
              }
            }
          }
        }

      }
    }
  }
  return J;
}

std::unordered_map<std::string, std::complex<double>> J_qq_eikonal(double *pp, double *q1, double *q2, amplitude& A) {
  // Irreducible component of the double soft current
  std::unordered_map<std::string, std::complex<double>> J;

  for(int s = -1; s <=1; s+=2) {
    std::vector<std::complex<double>> u = spinor(q1, s);
    std::vector<std::complex<double>> v = (-1.)*spinor(q2, s);
    for (int i = 0; i < A.process.size(); i++) { // particle index
      if(A.particle_type[i] == 0) continue;
      double pi[4];
      part(pp, pi, 4 * i, 4 * i + 4);
      // momentum dependence
      std::complex<double> gamma = -(u * (gamma0 * p_slash(pi) * v))/2./minkovski(q1, q2)/(minkovski(pi, q1) + minkovski(pi, q2));
      for (int a1 = 0; a1 < 3; a1++) for(int a2 = 0; a2 < 3; a2++) {  // Color index
        for (int b = 0; b < A.process[i]; b++) {
          for (int c = 0; c < A.process[i]; c++) {
            std::string key = std::to_string(i) + std::to_string(s) + std::to_string(a1) + std::to_string(a2) + std::to_string(b) + std::to_string(c);
            J[key] = 0;
            for(int d = 0; d < 8; d++) {
              if (A.particle_type[i] == 1) {
                J[key] += gs*gs*lam[d][3*a1+a2]/2. * lam[d][3 * b + c]/2.*gamma;
              }
              else if (A.particle_type[i] == -1) {
                J[key] += -gs*gs*lam[d][3*a1+a2]/2. * lam[d][3 * c + b]/2.*gamma;
              }
              else if (A.particle_type[i] == 2) {
                J[key] += gs*gs*lam[d][3*a1+a2]/2. * I*fabc[d][8 * c + b]*gamma;
              }
            }
          }
        }
      }
    }
  }
  return J;
}

std::unordered_map<std::string, std::complex<double>> J_qqg_eikonal(double *pp, double *q1, double *q2, double *q3, amplitude& A) {
  // Irreducible component of the quark-anti-quark-gluon soft current
  std::unordered_map<std::string, std::complex<double>> J;
  double q123[4] = {q1[0]+q2[0]+q3[0],q1[1]+q2[1]+q3[1],q1[2]+q2[2]+q3[2],q1[3]+q2[3]+q3[3]};
  double q13[4] = {q1[0]+q3[0],q1[1]+q3[1],q1[2]+q3[2],q1[3]+q3[3]};
  double q23[4] = {q2[0]+q3[0],q2[1]+q3[1],q2[2]+q3[2],q2[3]+q3[3]};
  double q12[4] = {q1[0]+q2[0],q1[1]+q2[1],q1[2]+q2[2],q1[3]+q2[3]};

  for(int s = -1; s <=1; s+=2) {
    std::vector<std::complex<double>> u = spinor(q1, s);
    std::vector<std::complex<double>> v = spinor(q2, s);
    for(int mu = 0; mu < 4; mu++) {
      for (int i = 0; i < A.process.size(); i++) { // particle index
        if(A.particle_type[i] == 0) continue;
        double pi[4];
        part(pp, pi, 4 * i, 4 * i + 4);
        // momentum dependence
        std::complex<double> gamma_ab = (u*((gamma0*p_slash(pi)*p_slash(q23)*gammas[mu]*(1./minkovski(q23, q23)) - gamma0*gammas[mu]*p_slash(q13)*p_slash(pi)*(1./minkovski(q13, q13)))*v))/minkovski(q123, q123)/minkovski(pi, q123);
        std::complex<double> gamma_nab = (u*(((gamma0*p_slash(pi))*(pi[mu]/minkovski(q12, q12)*(1./minkovski(pi, q3) - 1./minkovski(pi, q12)))
                  + ((gamma0*gammas[mu]*2.*(minkovski(pi, q12) - minkovski(pi, q3)) - gamma0*p_slash(pi)*4.*q12[mu] + gamma0*p_slash(q3)*4.*pi[mu])*(1./minkovski(q12, q12))
                                              - gamma0*gammas[mu]*p_slash(q13)*p_slash(pi)*(1./minkovski(q13,q13)) - gamma0*p_slash(pi)*p_slash(q23)*gammas[mu]*(1./minkovski(q23, q23)))*(1./minkovski(q123, q123)))*v))/minkovski(pi, q123);
        for (int a1 = 0; a1 < 3; a1++) for(int a2 = 0; a2 < 3; a2++) for(int a3 = 0; a3 < 8; a3++) {  // Color index
          for (int ci1 = 0; ci1 < A.process[i]; ci1++) for(int ci2 = 0; ci2 < A.process[i]; ci2++) {
            std::string key = std::to_string(i) + std::to_string(s) + std::to_string(mu) + std::to_string(a1) + std::to_string(a2) + std::to_string(a3) + std::to_string(ci1) + std::to_string(ci2);
            J[key] = 0;
            for(int c = 0; c < 8; c++) for(int d = 0; d < 3; d++) {
              if (A.particle_type[i] == 1) {
                J[key] += gs*gs*gs*(lam[a3][3*a1+d]*lam[c][3*d+a2] - lam[c][3*a1+d]*lam[a3][3*d+a2])/8.*lam[c][3*ci1 + ci2]/2.*gamma_nab;
              }
              else if (A.particle_type[i] == -1) {
                J[key] += -gs*gs*gs*(lam[a3][3*a1+d]*lam[c][3*d+a2] - lam[c][3*a1+d]*lam[a3][3*d+a2])/8.*lam[c][3*ci2 + ci1]/2.*gamma_nab;
              }
              else if (A.particle_type[i] == 2) {
                J[key] += gs*gs*gs*(lam[a3][3*a1+d]*lam[c][3*d+a2] - lam[c][3*a1+d]*lam[a3][3*d+a2])/8.*I*fabc[c][8*ci2 + ci1]*gamma_nab;
              }
            }
            for (int c = 0; c < 8; c++) for(int d = 0; d < 3; d++) {
              if (A.particle_type[i] == 1) {
                J[key] += gs*gs*gs*(lam[a3][3*a1+d]*lam[c][3*d+a2] + lam[c][3*a1+d]*lam[a3][3*d+a2])/8.*lam[c][3*ci1 + ci2]/2.*gamma_ab;
              }
              else if (A.particle_type[i] == -1) {
                J[key] += -gs*gs*gs*(lam[a3][3*a1+d]*lam[c][3*d+a2] + lam[c][3*a1+d]*lam[a3][3*d+a2])/8.*lam[c][3*ci2 + ci1]/2.*gamma_ab;
              }
              else if (A.particle_type[i] == 2) {
                J[key] += gs*gs*gs*(lam[a3][3*a1+d]*lam[c][3*d+a2] + lam[c][3*a1+d]*lam[a3][3*d+a2])/8.*I*fabc[c][8*ci2 + ci1]*gamma_ab;
              }
            }
          }
        }
      }
    }
  }
  return J;
}

std::complex<double> soft_g(double *pp_full, std::unordered_map<std::string, std::complex<double>> J,
    std::unordered_map<std::string, std::complex<double>> M,
    std::vector<int> helicities_full, std::vector<int> colors_full, amplitude& A) {

  double q[4];
  part(pp_full, q, A.process.size() * 4, A.process.size() * 4 + 4);
  std::vector<std::complex<double>> eps = polarization(q, helicities_full[A.process.size()]);

  std::complex<double> M_eikonal = 0;

  for (int i = 0; i < A.process.size(); i++) { // Particle sum
    for(int c = 0; c < A.process[i]; c++) { // Color summataton
      std::string key = "";
      for (int l = 0; l < A.process.size(); l++) {
        key += std::to_string(helicities_full[l]);
      }
      for (int l = 0; l < A.process.size(); l++) {
        key += std::to_string(colors_full[l]);
      }
      key.replace(key.size() - A.process.size() + i, 1, std::to_string(c));
      for(int mu = 0; mu < 4; mu++) {
        std::string key1 = std::to_string(i) + std::to_string(mu) + std::to_string(colors_full[A.process.size()]) + std::to_string(colors_full[i]) + std::to_string(c);
        M_eikonal += J[key1]*M[key]*eps[mu]*metric[mu][mu];
      }
    }
  }
  return M_eikonal;
}

std::complex<double> soft_gg(double *pp_full, std::unordered_map<std::string, std::complex<double>> J1,
    std::unordered_map<std::string, std::complex<double>> J2, std::unordered_map<std::string, std::complex<double>> J12,
    std::unordered_map<std::string, std::complex<double>> M,
    std::vector<int> helicities_full, std::vector<int> colors_full, amplitude& A) {

  double q1[4], q2[4];
  part(pp_full, q1, A.process.size() * 4, A.process.size() * 4 + 4);
  part(pp_full, q2, A.process.size() * 4 + 4, A.process.size() * 4 + 8);
  std::vector<std::complex<double>> eps1 = polarization(q1, helicities_full[A.process.size()]);
  std::vector<std::complex<double>> eps2 = polarization(q2, helicities_full[A.process.size() + 1]);

  std::complex<double> M_eikonal = 0;

  for (int i1 = 0; i1 < A.process.size(); i1++) {
    double pi1[4];
    part(pp_full, pi1, i1 * 4, i1 * 4 + 4);
    for (int i2 = 0; i2 < A.process.size(); i2++) {
      double pi2[4];
      part(pp_full, pi2, i2 * 4, i2 * 4 + 4);
      for (int c1 = 0; c1 < A.process[i1]; c1++) {
        for (int c2 = 0; c2 < A.process[i2]; c2++) {
          std::string key = "";
          for (int l = 0; l < A.process.size(); l++) {
            key += std::to_string(helicities_full[l]);
          }

          for (int l = 0; l < A.process.size(); l++) {
            key += std::to_string(colors_full[l]);
          }

          if (i1 == i2) {
            key.replace(key.size() - A.process.size() + i1, 1, std::to_string(c2));
          }
          else {
            key.replace(key.size() - A.process.size() + i1, 1, std::to_string(c1));
            key.replace(key.size() - A.process.size() + i2, 1, std::to_string(c2));
          }

          for (int mu1 = 0; mu1 < 4; mu1++) {
            for (int mu2 = 0; mu2 < 4; mu2++) {
              if (i1 == i2) {
                std::string key1 = std::to_string(i1) + std::to_string(mu1) + std::to_string(colors_full[A.process.size()]) + std::to_string(colors_full[i1]) + std::to_string(c1);
                std::string key2 = std::to_string(i2) + std::to_string(mu2) + std::to_string(colors_full[A.process.size() + 1]) +  std::to_string(c1) + std::to_string(c2);
                M_eikonal += 1./2. * J1[key1] * J2[key2] * M[key] * eps1[mu1] * eps2[mu2] * metric[mu1][mu1] * metric[mu2][mu2];
                key1 = std::to_string(i1) + std::to_string(mu1) + std::to_string(colors_full[A.process.size()]) +  std::to_string(c1) + std::to_string(c2) ;
                key2 = std::to_string(i2) + std::to_string(mu2) + std::to_string(colors_full[A.process.size() + 1]) + std::to_string(colors_full[i2]) + std::to_string(c1);
                M_eikonal += 1./2. * J2[key2] * J1[key1] * M[key] * eps1[mu1] * eps2[mu2] * metric[mu1][mu1] * metric[mu2][mu2];
              }
              else  {
                std::string key1 = std::to_string(i1) + std::to_string(mu1) + std::to_string(colors_full[A.process.size()]) +  std::to_string(colors_full[i1]) + std::to_string(c1) ;
                std::string key2 = std::to_string(i2) + std::to_string(mu2) + std::to_string(colors_full[A.process.size() + 1]) + std::to_string(colors_full[i2]) + std::to_string(c2);
                M_eikonal += 2 * 1./2. * J1[key1] * J2[key2] * M[key] * eps1[mu1] * eps2[mu2] * metric[mu1][mu1] * metric[mu2][mu2];
              }
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < A.process.size(); i++) {
    double pi[4];
    part(pp_full, pi, i * 4, i * 4 + 4);
    for (int c = 0; c < A.process[i]; c++) {
      std::string key = "";
      for (int l = 0; l < A.process.size(); l++) {
        key += std::to_string(helicities_full[l]);
      }

      for (int l = 0; l < A.process.size(); l++) {
        key += std::to_string(colors_full[l]);
      }
      key.replace(key.size() - A.process.size() + i, 1, std::to_string(c));

      for (int mu1 = 0; mu1 < 4; mu1++) for (int mu2 = 0; mu2 < 4; mu2++) {
        std::string key12 = std::to_string(i) + std::to_string(mu1) + std::to_string(mu2) + std::to_string(colors_full[A.process.size()]) + std::to_string(colors_full[A.process.size() + 1]) + std::to_string(colors_full[i]) + std::to_string(c);
        M_eikonal += J12[key12] * M[key] * eps1[mu1] * eps2[mu2] * metric[mu1][mu1] * metric[mu2][mu2];
      }
    }
  }
  return M_eikonal;
}

std::complex<double> soft_qq(double *pp_full, std::unordered_map<std::string, std::complex<double>> J12,
    std::unordered_map<std::string, std::complex<double>> M,
    std::vector<int> helicities_full, std::vector<int> colors_full, amplitude& A) {

  if(helicities_full[A.process.size()] != helicities_full[A.process.size()])
    return 0.;

  std::complex<double> M_eikonal = 0;

  for (int i = 0; i < A.process.size(); i++) { // Particle sum
    for(int c = 0; c < A.process[i]; c++) { // Color summataton
      std::string key = "";
      for (int l = 0; l < A.process.size(); l++) {
        key += std::to_string(helicities_full[l]);
      }
      for (int l = 0; l < A.process.size(); l++) {
        key += std::to_string(colors_full[l]);
      }
      key.replace(key.size() - A.process.size() + i, 1, std::to_string(c));
      std::string key1 = std::to_string(i) + std::to_string(helicities_full[A.process.size()]) + std::to_string(colors_full[A.process.size()]) + std::to_string(colors_full[A.process.size() + 1]) + std::to_string(colors_full[i]) + std::to_string(c);
      M_eikonal += J12[key1]*M[key];
    }
  }
  return M_eikonal;
}

std::complex<double> soft_qqg(double *pp_full, std::unordered_map<std::string, std::complex<double>> J12, std::unordered_map<std::string, std::complex<double>> J3,
    std::unordered_map<std::string, std::complex<double>> J123, std::unordered_map<std::string, std::complex<double>> M,
    std::vector<int> helicities_full, std::vector<int> colors_full, amplitude& A) {

  double q1[4], q2[4], q3[4];
  part(pp_full, q1, A.process.size() * 4, A.process.size() * 4 + 4);
  part(pp_full, q2, A.process.size() * 4 + 4, A.process.size() * 4 + 8);
  part(pp_full, q3, A.process.size() * 4 + 8, A.process.size() * 4 + 12);
  std::vector<std::complex<double>> eps = polarization(q3, helicities_full[A.process.size() + 2]);
  std::complex<double> M_eikonal = 0;

  for(int i = 0; i < A.process.size(); i++) for(int j = 0; j < A.process.size(); j++) {
    double pi[4], pj[4];
    part(pp_full, pi, i*4, i*4 + 4);
    part(pp_full, pj, j*4, j*4 + 4);
    for(int ci = 0; ci < A.process[i]; ci++) for(int cj = 0; cj < A.process[j]; cj++) {
      std::string key = "";
      for (int l = 0; l < A.process.size(); l++) {
        key += std::to_string(helicities_full[l]);
      }

      for (int l = 0; l < A.process.size(); l++) {
        key += std::to_string(colors_full[l]);
      }

      if (i == j) {
        key.replace(key.size() - A.process.size() + i, 1, std::to_string(ci));
      }
      else {
        key.replace(key.size() - A.process.size() + i, 1, std::to_string(ci));
        key.replace(key.size() - A.process.size() + j, 1, std::to_string(cj));
      }
      for (int mu = 0; mu < 4; mu++) {
        if (i == j) {
          std::string key1 = std::to_string(i) + std::to_string(mu) + std::to_string(colors_full[A.process.size()+2]) + std::to_string(colors_full[i]) + std::to_string(cj);
          std::string key2 = std::to_string(i) + std::to_string(helicities_full[A.process.size()]) + std::to_string(colors_full[A.process.size()]) + std::to_string(colors_full[A.process.size() + 1]) +  std::to_string(cj) + std::to_string(ci);
          M_eikonal += 1./2.*J3[key1]*J12[key2]*M[key]*eps[mu]*metric[mu][mu];
          key1 = std::to_string(i) + std::to_string(mu) + std::to_string(colors_full[A.process.size()+2]) + std::to_string(cj) + std::to_string(ci);
          key2 = std::to_string(i) + std::to_string(helicities_full[A.process.size()]) + std::to_string(colors_full[A.process.size()]) + std::to_string(colors_full[A.process.size() + 1]) +  std::to_string(colors_full[i]) + std::to_string(cj);
          M_eikonal += 1./2.*J12[key2]*J3[key1]*M[key]*eps[mu]*metric[mu][mu];
        }
        else  {
          std::string key1 = std::to_string(i) + std::to_string(mu) + std::to_string(colors_full[A.process.size()+2]) + std::to_string(colors_full[i]) + std::to_string(ci);
          std::string key2 = std::to_string(j) + std::to_string(helicities_full[A.process.size()]) + std::to_string(colors_full[A.process.size()]) + std::to_string(colors_full[A.process.size() + 1]) +  std::to_string(colors_full[j]) + std::to_string(cj);
          M_eikonal += 2.*1./2.*J3[key1]*J12[key2]*M[key]*eps[mu]*metric[mu][mu];
        }
      }
    }
  }
  for(int i = 0; i < A.process.size(); i++) {
    double pi[4], pj[4];
    part(pp_full, pi, i*4, i*4 + 4);
    for(int ci = 0; ci < A.process[i]; ci++) {
      std::string key = "";
      for (int l = 0; l < A.process.size(); l++) {
        key += std::to_string(helicities_full[l]);
      }

      for (int l = 0; l < A.process.size(); l++) {
        key += std::to_string(colors_full[l]);
      }
      key.replace(key.size() - A.process.size() + i, 1, std::to_string(ci));
      for (int mu = 0; mu < 4; mu++) {
        std::string key1 = std::to_string(i) + std::to_string(helicities_full[A.process.size()]) + std::to_string(mu) + std::to_string(colors_full[A.process.size()]) + std::to_string(colors_full[A.process.size() + 1]) + std::to_string(colors_full[A.process.size() + 2]) +  std::to_string(colors_full[i]) + std::to_string(ci);
        M_eikonal += J123[key1]*M[key]*eps[mu]*metric[mu][mu];
      }
    }
  }
  return M_eikonal;
}


std::complex<double> double_soft_tree(double *pp_full, std::unordered_map<std::string, std::complex<double>> J1,
    std::unordered_map<std::string, std::complex<double>> J2, std::unordered_map<std::string, std::complex<double>> M,
    std::vector<int> helicities_full, std::vector<int> colors_full, amplitude& A) {

  double q1[4], q2[4];
  part(pp_full, q1, A.process.size() * 4, A.process.size() * 4 + 4);
  part(pp_full, q2, A.process.size() * 4 + 4, A.process.size() * 4 + 8);
  std::vector<std::complex<double>> eps1 = polarization(q1, helicities_full[A.process.size()]);
  std::vector<std::complex<double>> eps2 = polarization(q2, helicities_full[A.process.size() + 1]);

  std::complex<double> M_eikonal = 0;

  for (int i1 = 0; i1 < A.process.size(); i1++) {
    double pi1[4];
    part(pp_full, pi1, i1 * 4, i1 * 4 + 4);
    for (int i2 = 0; i2 < A.process.size(); i2++) {
      double pi2[4];
      part(pp_full, pi2, i2 * 4, i2 * 4 + 4);
      for (int c1 = 0; c1 < A.process[i1]; c1++) {
        for (int c2 = 0; c2 < A.process[i2]; c2++) {
          std::string key = "";
          for (int l = 0; l < A.process.size(); l++) {
            key += std::to_string(helicities_full[l]);
          }

          for (int l = 0; l < A.process.size(); l++) {
            key += std::to_string(colors_full[l]);
          }

          if (i1 == i2) {
            key.replace(key.size() - A.process.size() + i1, 1, std::to_string(c2));
          }
          else {
            key.replace(key.size() - A.process.size() + i1, 1, std::to_string(c1));
            key.replace(key.size() - A.process.size() + i2, 1, std::to_string(c2));
          }

          for (int mu1 = 0; mu1 < 4; mu1++) {
            for (int mu2 = 0; mu2 < 4; mu2++) {
              if (i1 == i2) {
                std::string key1 = std::to_string(i1) + std::to_string(mu1) + std::to_string(colors_full[A.process.size()]) + std::to_string(colors_full[i1]) + std::to_string(c1);
                std::string key2 = std::to_string(i2) + std::to_string(mu2) + std::to_string(colors_full[A.process.size() + 1]) +  std::to_string(c1) + std::to_string(c2);
                M_eikonal += 1./2. * J1[key1] * J2[key2] * M[key] * eps1[mu1] * eps2[mu2] * metric[mu1][mu1] * metric[mu2][mu2];
                key1 = std::to_string(i1) + std::to_string(mu1) + std::to_string(colors_full[A.process.size()]) +  std::to_string(c1) + std::to_string(c2) ;
                key2 = std::to_string(i2) + std::to_string(mu2) + std::to_string(colors_full[A.process.size() + 1]) + std::to_string(colors_full[i2]) + std::to_string(c1);
                M_eikonal += 1./2. * J2[key2] * J1[key1] * M[key] * eps1[mu1] * eps2[mu2] * metric[mu1][mu1] * metric[mu2][mu2];
              }
              else  {
                std::string key1 = std::to_string(i1) + std::to_string(mu1) + std::to_string(colors_full[A.process.size()]) +  std::to_string(colors_full[i1]) + std::to_string(c1) ;
                std::string key2 = std::to_string(i2) + std::to_string(mu2) + std::to_string(colors_full[A.process.size() + 1]) + std::to_string(colors_full[i2]) + std::to_string(c2);
                M_eikonal += 2 * 1./2. * J1[key1] * J2[key2] * M[key] * eps1[mu1] * eps2[mu2] * metric[mu1][mu1] * metric[mu2][mu2];
                //key1 = std::to_string(i1) + std::to_string(mu1) + std::to_string(colors_full[A.process.size()])+  std::to_string(colors_full[i1]) + std::to_string(c1) ;
                //key2 = std::to_string(i2) + std::to_string(mu2) + std::to_string(colors_full[A.process.size() + 1]) + std::to_string(colors_full[i2]) + std::to_string(c2);
                //M_eikonal += 1./2. * J2[key2] * J1[key1] * M[key] * eps1[mu1] * eps2[mu2] * metric[mu1][mu1] * metric[mu2][mu2];
              }
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < A.process.size(); i++) {
    double pi[4];
    part(pp_full, pi, 4 * i, 4 * i + 4);
    std::string key = "";
    for (int l = 0; l < A.process.size(); l++) {
      key += std::to_string(helicities_full[l]);
    }

    for (int l = 0; l < A.process.size(); l++) {
      key += std::to_string(colors_full[l]);
    }
    std::complex<double> nonabel = 0;
    for (int mu1 = 0; mu1 < 4; mu1++) {
      for (int mu2 = 0; mu2 < 4; mu2++) {
        nonabel += gs * gs * ((pi[mu1] * q1[mu2] - pi[mu2] * q2[mu1])/minkovski(q1, q2)/(minkovski(pi, q1) + minkovski(pi, q2))
        - (minkovski(pi, q1) - minkovski(pi, q2))/2./(minkovski(pi, q1) + minkovski(pi, q2))
            * (pi[mu1] * pi[mu2]/minkovski(pi, q1)/minkovski(pi, q2) + metric[mu1][mu2]/minkovski(q1, q2))) * eps1[mu1] * eps2[mu2] * metric[mu1][mu1] * metric[mu2][mu2];
      }
    }
    std::complex<double> colfac = 0;
    for (int a = 0; a < 8; a++) {
      for (int ci = 0; ci < A.process[i]; ci++) {
        std::complex<double> colfac_con = 0;
        if (A.particle_type[i] == 1) {
          colfac_con = lam[a][3 * colors_full[i] + ci]/2.;
        }
        else if (A.particle_type[i] == -1) {
          colfac_con = -lam[a][3 * ci + colors_full[i]]/2.;
        }
        else if (A.particle_type[i] == 2) {
          colfac_con = I * fabc[a][8 * ci + colors_full[i]];
        }
        key.replace(key.size() - A.process.size() + i, 1, std::to_string(ci));
        colfac_con *= M[key] * I * fabc[a][8 * colors_full[A.process.size()] + colors_full[A.process.size() + 1]];
        colfac += colfac_con;
      }
    }
    M_eikonal += colfac * nonabel;
  }
  return M_eikonal;
}

std::complex<double> double_soft_tree_qq(double *pp_full, std::unordered_map<std::string, std::complex<double>> M,
                     std::vector<int> helicities_full, std::vector<int> colors_full, amplitude& A) {
  double q1[4], q2[4];
  part(pp_full, q1, A.process.size() * 4, A.process.size() * 4 + 4);
  part(pp_full, q2, A.process.size() * 4 + 4, A.process.size() * 4 + 8);
  std::vector<std::complex<double>> eps1 = polarization(q1, helicities_full[A.process.size()]);
  std::vector<std::complex<double>> eps2 = polarization(q2, helicities_full[A.process.size() + 1]);
  std::vector<std::complex<double>> u = spinor(q1, helicities_full[A.process.size()]);
  std::vector<std::complex<double>> v = spinor(q2, -helicities_full[A.process.size() + 1]);

  std::complex<double> M_eikonal = 0;

  for (int i = 0; i < A.process.size(); i++) {
    double pi[4];
    part(pp_full, pi, 4 * i, 4 * i + 4);
    std::string key = "";
    for (int l = 0; l < A.process.size(); l++) {
      key += std::to_string(helicities_full[l]);
    }

    for (int l = 0; l < A.process.size(); l++) {
      key += std::to_string(colors_full[l]);
    }
    std::complex<double> nonabel = 0;

    nonabel += gs * gs * (u * (gamma0 * p_slash(pi) * v))/2./minkovski(q1, q2)/(minkovski(pi, q1) + minkovski(pi, q2));
    std::complex<double> colfac = 0;
    for (int a = 0; a < 8; a++) {
      for (int ci = 0; ci < A.process[i]; ci++) {
        std::complex<double> colfac_con = 0;
        if (A.particle_type[i] == 1) {
          colfac_con = lam[a][3 * colors_full[i] + ci]/2.;
        }
        else if (A.particle_type[i] == -1) {
          colfac_con = -lam[a][3 * ci + colors_full[i]]/2.;
        }
        else if (A.particle_type[i] == 2) {
          colfac_con = I * fabc[a][8 * ci + colors_full[i]];
        }
        key.replace(key.size() - A.process.size() + i, 1, std::to_string(ci));
        colfac_con *= M[key] * lam[a][3 * colors_full[A.process.size()] + colors_full[A.process.size() + 1]]/2.;
        colfac += colfac_con;
      }
    }
    M_eikonal += colfac * nonabel;
  }
  return M_eikonal;
}

double soft_g_squared(double *pp_full, std::unordered_map<std::string, double> M_ij, amplitude& A) {
  double q[4];
  part(pp_full, q, A.process.size()*4, A.process.size()*4 + 4);
  double approx = 0;
  for(int i = 0; i < A.process.size(); i++) {
    if(A.process[i] == 1) continue;
    double pi[4];
    part(pp_full, pi, i*4, i*4 + 4);
    for(int j = 0; j < A.process.size(); j++) {
      if(A.process[j] == 1) continue;
      if(i == j) continue;
      double pj[4];
      part(pp_full, pj, j*4, j*4 + 4);
      approx += -gs*gs*minkovski(pi, pj)/minkovski(pi, q)/minkovski(pj, q)*std::real(M_ij[std::to_string(i) + std::to_string(j)]);
    }
  }
  return approx;
}

double soft_qq_squared(double *pp_full, std::unordered_map<std::string, double> M_ij, amplitude& A) {
  double q1[4], q2[4];
  part(pp_full, q1, A.process.size()*4, A.process.size()*4 + 4);
  part(pp_full, q2, A.process.size()*4 + 4, A.process.size()*4 + 8);
  double q12[4];
  add_arr(q1, q2, q12, 4);
  double approx = 0;
  for(int i = 0; i < A.process.size(); i++) {
    if(A.process[i] == 1) continue;
    double pi[4];
    part(pp_full, pi, 4*i, 4*i + 4);
    for(int j = 0; j < A.process.size(); j++) {
      if(A.process[j] == 1) continue;
      double pj[4];
      part(pp_full, pj, 4*j, 4*j + 4);
      approx += T_F*std::pow(gs, 4)*(minkovski(pi, q1)*minkovski(pj, q2) + minkovski(pi, q2)*minkovski(pj, q1) - minkovski(pi, pj)*minkovski(q1, q2))/std::pow(minkovski(q1, q2), 2)/minkovski(pi, q12)/minkovski(pj, q12)
        *std::real(M_ij[std::to_string(i) + std::to_string(j)]);
    }
  }
  return approx;
}

double soft_gg_squared(double*pp_full, std::unordered_map<std::string, double> Mij, std::unordered_map<std::string, double> Mijkl, amplitude& A) {
  double q1_arr[4], q2_arr[4];
  part(pp_full, q1_arr, A.process.size()*4, A.process.size()*4 + 4);
  part(pp_full, q2_arr, A.process.size()*4 + 4, A.process.size()*4 + 8);
  LV<double> q1(q1_arr);
  LV<double> q2(q2_arr);
  double approx = 0;
  for(int i = 0; i < A.process.size(); i++) {
    if(A.process[i] == 1) continue;
    double pi_arr[4];
    part(pp_full, pi_arr, i*4, i*4 + 4);
    LV<double> pi(pi_arr);
    for(int j = 0; j < A.process.size(); j++) {
      if(A.process[j] == 1) continue;
      double pj_arr[4];
      part(pp_full, pj_arr, 4*j, 4*j + 4);
      LV<double> pj(pj_arr);

      // Reducible part
      for(int k = 0; k < A.process.size(); k++) {
        if(A.process[k] == 1) continue;
        double pk_arr[4];
        part(pp_full, pk_arr, 4*k, 4*k + 4);
        LV<double> pk(pk_arr);
        for(int l = 0; l < A.process.size(); l++) {
          if(A.process[l] == 1) continue;
          double pl_arr[4];
          part(pp_full, pl_arr, 4*l, 4*l + 4);
          LV<double> pl(pl_arr);
          approx += std::pow(gs, 4)*(j1(pi, q1)*j1(pj, q1))*(j1(pk, q2)*j1(pl, q2))*(std::real(Mijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)])
                                                                                   + std::real(Mijkl[std::to_string(k) + std::to_string(l) + std::to_string(i) + std::to_string(j)]))/2. ;
        }
      }
      // Irreducible part
      approx += std::pow(gs, 4)*C_A*std::real(Mij[std::to_string(i) + std::to_string(j)])
          *((gamma20(pi, q1, q2)*gamma20(pj, q1, q2).transpose()).trace()
            +j1(pi, q1)*(gamma20(pi, q1, q2)*j1(pj, q2))
            -j1(pi, q1)*(gamma20(pj, q1, q2)*j1(pj, q2))
            -1./2.*(j1(pi, q1)*j1(pi, q1))*(j1(pi, q2)*j1(pj, q2))
            -1./2.*(j1(pj, q2)*j1(pj, q2))*(j1(pi, q1)*j1(pj, q1))
            +3./4.*(j1(pi, q1)*j1(pj, q1))*(j1(pi, q2)*j1(pj, q2)));
    }
  }
  return approx;
}

double den(double arg) {
  return 1./arg;
}

double soft_gqq_squared(double* pp_full, std::unordered_map<std::string, double> M_ij,
                                         std::unordered_map<std::string, double> M_ijkl,
                                         std::unordered_map<std::string, double> dM_ijk, amplitude& A) {
  double q1[4], q2[4], q3[4], q12[4], q13[4], q23[4], q123[4];
  part(pp_full, q1, A.process.size()*4, A.process.size()*4 + 4);
  part(pp_full, q2, A.process.size()*4 + 4, A.process.size()*4 + 8);
  part(pp_full, q3, A.process.size()*4 + 8, A.process.size()*4 + 12);
  add_arr(q1, q2, q12, 4);
  add_arr(q1, q3, q13, 4);
  add_arr(q2, q3, q23, 4);
  add_arr(q12, q3, q123, 4);
  double q1q2 = minkovski(q1, q2);
  double q1q3 = minkovski(q1, q3);
  double q2q3 = minkovski(q2, q3);

  double approx = 0.;
  for(int i = 0; i < A.process.size(); i++) {
    if(A.process[i] == 1) continue;
    double pi[4];
    part(pp_full, pi, i*4, i*4 + 4);
    for(int j = 0; j < A.process.size(); j++) {
      if(A.process[j] == 1) continue;
      double pj[4];
      part(pp_full, pj, 4*j, 4*j + 4);
      double piq1 = minkovski(pi, q1);
      double piq2 = minkovski(pi, q2);
      double piq3 = minkovski(pi, q3);
      double piq12 = minkovski(pi, q12);
      double piq13 = minkovski(pi, q13);
      double pjq1 = minkovski(pj, q1);
      double pjq2 = minkovski(pj, q2);
      double pjq12 = minkovski(pj, q12);
      double pipj = minkovski(pi, pj);
      double piq23 = minkovski(pi, q23);
      double mi2 = minkovski(pi, pi);
      double pjq3 = minkovski(pj, q3);
      double pjq23 = minkovski(pj, q23);
      double piq123 = minkovski(pi, q123);
      double q123q123 = minkovski(q123, q123);

      double Sij1_1 = pipj/piq1/pjq1;
      for(int k = 0; k < A.process.size(); k++) {
        if(A.process[k] == 1) continue;
        double pk[4];
        part(pp_full, pk, 4*k, 4*k + 4);
        double pkq1 = minkovski(pk, q1);
        double pkq2 = minkovski(pk, q2);
        double pkq3 = minkovski(pk, q3);
        double pipk = minkovski(pi, pk);
        double pjpk = minkovski(pj, pk);

        double Sijk =
          + den(q1q3*q2q3)*den(2*piq1*pjq3*q2q3*pkq3 + 2*piq1*pjq3*q2q3*pkq2 + 2*
          piq1*pjq3*q2q3*pkq1 + 2*piq1*pjq3*q1q3*pkq3 + 2*piq1*pjq3*q1q3*pkq2 + 2*
          piq1*pjq3*q1q3*pkq1 + 2*piq1*pjq3*q1q2*pkq3 + 2*piq1*pjq3*q1q2*pkq2 + 2*
          piq1*pjq3*q1q2*pkq1 + 2*piq1*pjq2*q2q3*pkq3 + 2*piq1*pjq2*q2q3*pkq2 + 2*
          piq1*pjq2*q2q3*pkq1 + 2*piq1*pjq2*q1q3*pkq3 + 2*piq1*pjq2*q1q3*pkq2 + 2*
          piq1*pjq2*q1q3*pkq1 + 2*piq1*pjq2*q1q2*pkq3 + 2*piq1*pjq2*q1q2*pkq2 + 2*
          piq1*pjq2*q1q2*pkq1) * (  - q1q2*pipj*pkq3 + pjq3*q1q2*pipk + piq3*q1q2*
            pjpk - 2*piq3*pjq3*pkq2 - 2*piq3*pjq2*pkq3 - piq3*pjq2*pkq1 - piq3*
            pjq1*pkq2 - piq2*pjq3*pkq1 + piq2*pjq1*pkq3 - piq1*pjq3*pkq2 - piq1*
            pjq2*pkq3 );

          Sijk +=  + den(q1q2*q2q3)*den(2*piq1*pjq3*q2q3*pkq3 + 2*piq1*pjq3*q2q3*
          pkq2 + 2*piq1*pjq3*q2q3*pkq1 + 2*piq1*pjq3*q1q3*pkq3 + 2*piq1*pjq3*q1q3*
          pkq2 + 2*piq1*pjq3*q1q3*pkq1 + 2*piq1*pjq3*q1q2*pkq3 + 2*piq1*pjq3*q1q2*
          pkq2 + 2*piq1*pjq3*q1q2*pkq1 + 2*piq1*pjq2*q2q3*pkq3 + 2*piq1*pjq2*q2q3*
          pkq2 + 2*piq1*pjq2*q2q3*pkq1 + 2*piq1*pjq2*q1q3*pkq3 + 2*piq1*pjq2*q1q3*
          pkq2 + 2*piq1*pjq2*q1q3*pkq1 + 2*piq1*pjq2*q1q2*pkq3 + 2*piq1*pjq2*q1q2*
          pkq2 + 2*piq1*pjq2*q1q2*pkq1) * ( q1q3*pipj*pkq2 - pjq2*q1q3*pipk + piq3
            *pjq2*pkq1 - piq3*pjq1*pkq2 - piq2*q1q3*pjpk + 2*piq2*pjq3*pkq2 +
            piq2*pjq3*pkq1 + 2*piq2*pjq2*pkq3 + piq2*pjq1*pkq3 + piq1*pjq3*pkq2
              + piq1*pjq2*pkq3 );

          Sijk +=  + den(q1q2)*den(2*piq1*pjq3*q2q3*pkq3 + 2*piq1*pjq3*q2q3*pkq2
          + 2*piq1*pjq3*q2q3*pkq1 + 2*piq1*pjq3*q1q3*pkq3 + 2*piq1*pjq3*q1q3*pkq2
          + 2*piq1*pjq3*q1q3*pkq1 + 2*piq1*pjq3*q1q2*pkq3 + 2*piq1*pjq3*q1q2*pkq2
          + 2*piq1*pjq3*q1q2*pkq1 + 2*piq1*pjq2*q2q3*pkq3 + 2*piq1*pjq2*q2q3*pkq2
          + 2*piq1*pjq2*q2q3*pkq1 + 2*piq1*pjq2*q1q3*pkq3 + 2*piq1*pjq2*q1q3*pkq2
          + 2*piq1*pjq2*q1q3*pkq1 + 2*piq1*pjq2*q1q2*pkq3 + 2*piq1*pjq2*q1q2*pkq2
          + 2*piq1*pjq2*q1q2*pkq1) * (  - pipj*pkq1 + pjq1*pipk - 2*piq2*pjpk -
            piq1*pjpk );

          Sijk +=  + den(q1q3)*den(2*piq1*pjq3*q2q3*pkq3 + 2*piq1*pjq3*q2q3*pkq2
          + 2*piq1*pjq3*q2q3*pkq1 + 2*piq1*pjq3*q1q3*pkq3 + 2*piq1*pjq3*q1q3*pkq2
          + 2*piq1*pjq3*q1q3*pkq1 + 2*piq1*pjq3*q1q2*pkq3 + 2*piq1*pjq3*q1q2*pkq2
          + 2*piq1*pjq3*q1q2*pkq1 + 2*piq1*pjq2*q2q3*pkq3 + 2*piq1*pjq2*q2q3*pkq2
          + 2*piq1*pjq2*q2q3*pkq1 + 2*piq1*pjq2*q1q3*pkq3 + 2*piq1*pjq2*q1q3*pkq2
          + 2*piq1*pjq2*q1q3*pkq1 + 2*piq1*pjq2*q1q2*pkq3 + 2*piq1*pjq2*q1q2*pkq2
          + 2*piq1*pjq2*q1q2*pkq1) * ( pipj*pkq1 - pjq1*pipk + 2*piq3*pjpk + piq1
            *pjpk );

          Sijk +=  + den(q2q3)*den(2*piq1*pjq3*q2q3*pkq3 + 2*piq1*pjq3*q2q3*pkq2
          + 2*piq1*pjq3*q2q3*pkq1 + 2*piq1*pjq3*q1q3*pkq3 + 2*piq1*pjq3*q1q3*pkq2
          + 2*piq1*pjq3*q1q3*pkq1 + 2*piq1*pjq3*q1q2*pkq3 + 2*piq1*pjq3*q1q2*pkq2
          + 2*piq1*pjq3*q1q2*pkq1 + 2*piq1*pjq2*q2q3*pkq3 + 2*piq1*pjq2*q2q3*pkq2
          + 2*piq1*pjq2*q2q3*pkq1 + 2*piq1*pjq2*q1q3*pkq3 + 2*piq1*pjq2*q1q3*pkq2
          + 2*piq1*pjq2*q1q3*pkq1 + 2*piq1*pjq2*q1q2*pkq3 + 2*piq1*pjq2*q1q2*pkq2
          + 2*piq1*pjq2*q1q2*pkq1) * (  - pipj*pkq3 + pipj*pkq2 - pjq3*pipk +
            pjq2*pipk + piq3*pjpk - piq2*pjpk );

        approx += - std::pow(gs, 6)*T_F*Sijk*std::real(dM_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)]);
        for(int l = 0; l < A.process.size(); l++) {
          if(A.process[l] == 1) continue;
          double pl[4];
          part(pp_full, pl, 4*l, 4*l + 4);
          double pkq2 = minkovski(pk, q2);
          double pkq3 = minkovski(pk, q3);
          double pkq23 = minkovski(pk, q23);
          double plq2 = minkovski(pl, q2);
          double plq3 = minkovski(pl, q3);
          double plq23 = minkovski(pl, q23);
          double plpk = minkovski(pl, pk);

          double Ikl = (pkq2*plq3 + pkq3*plq2 - plpk*q2q3)/q2q3/q2q3/pkq23/plq23;
          // reducible part
          approx += -std::pow(gs, 6)*T_F*Sij1_1*Ikl*(std::real(M_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)])
                                                 + std::real(M_ijkl[std::to_string(k) + std::to_string(l) + std::to_string(i) + std::to_string(j)]))/2.;

        }
      }
      double Sij1_3 = piq2/4./q2q3/q2q3/piq1/piq23*(pipj/pjq1*(3*pjq3/pjq23 - 2.*piq3/piq23) - 2.*mi2*pjq3/piq1/pjq23)
        + pipj/8./q2q3/piq1/piq23*(-3.*pipj/pjq1/pjq23 + 2.*mi2*(1./piq1/pjq23 + 1./pjq1/piq23))
        + piq3/4./q2q3/q2q3/piq1/piq23*(pipj/pjq1*(3*pjq2/pjq23 - 2.*piq2/piq23) - 2.*mi2*pjq2/piq1/pjq23)
        + pipj/8./q2q3/piq1/piq23*(-3.*pipj/pjq1/pjq23 + 2.*mi2*(1./piq1/pjq23 + 1./pjq1/piq23));

      approx += -std::pow(gs, 6)*T_F*C_A*Sij1_3*std::real(M_ij[std::to_string(i) + std::to_string(j)]);

      double Sij2_3 =
        +1./2./q2q3/piq123*(2./q123q123/q2q3*(pjq2*(2.*mi2*q1q3 - piq123*piq3)/piq1/pjq23
        - piq2*(2.*pipj*q1q3 - (piq1 - piq2 + 3.*piq3)*pjq3)/pjq1/piq23)
        + 1./q123q123/q1q2*((2.*piq12*(pipj*q2q3 - piq2*pjq3 - piq3*pjq2) + mi2*(pjq2*q1q3 - pjq1*q2q3))/piq1/pjq23
        - (mi2*((pjq1 + 2.*pjq2)*q2q3 + pjq2*q1q3) - 2.*(piq2*pjq1 + (piq1 + 2.*piq2)*pjq2)*piq3)/pjq1/piq23)
        + 1./q123q123*((pipj*piq123 + mi2*(pjq2 - 2.*pjq1))/piq1/pjq23 - (3.*mi2*pjq2 - pipj*piq1)/pjq1/piq23)
        + 1./q2q3*(piq2*(1./piq1 - 1./piq23)*(mi2*pjq3/piq1/pjq23 - pipj*piq3/pjq1/piq23))
        + 1./2.*mi2*pipj*(1./piq23 - 1./piq1)*(1./piq1/pjq23 - 1./pjq1/piq23))

        + 1./2./q2q3/piq123*(2./q123q123/q2q3*(pjq3*(2.*mi2*q1q2 - piq123*piq2)/piq1/pjq23
        - piq3*(2.*pipj*q1q2 - (piq1 - piq3 + 3.*piq2)*pjq2)/pjq1/piq23)
        + 1./q123q123/q1q3*((2.*piq13*(pipj*q2q3 - piq3*pjq2 - piq2*pjq3) + mi2*(pjq3*q1q2 - pjq1*q2q3))/piq1/pjq23
        - (mi2*((pjq1 + 2.*pjq3)*q2q3 + pjq3*q1q2) - 2.*(piq3*pjq1 + (piq1 + 2.*piq3)*pjq3)*piq2)/pjq1/piq23)
        + 1./q123q123*((pipj*piq123 + mi2*(pjq3 - 2.*pjq1))/piq1/pjq23 - (3.*mi2*pjq3 - pipj*piq1)/pjq1/piq23)
        + 1./q2q3*(piq3*(1./piq1 - 1./piq23)*(mi2*pjq2/piq1/pjq23 - pipj*piq2/pjq1/piq23))
        + 1./2.*mi2*pipj*(1./piq23 - 1./piq1)*(1./piq1/pjq23 - 1./pjq1/piq23));


      approx += - std::pow(gs, 6)*T_F*C_A*Sij2_3*std::real(M_ij[std::to_string(i) + std::to_string(j)]);

      double d = 4.;
      double Sij3_3_ab =
       + den(4*pow(q2q3,2) + 8*q1q3*q2q3 + 4*pow(q1q3,2) + 8*q1q2*q2q3 + 8*
      q1q2*q1q3 + 4*pow(q1q2,2))*den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*
      pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * (
          - 16*pipj + 4*d*pipj );

      Sij3_3_ab +=  + den(q1q2*q1q3)*den(4*pow(q2q3,2) + 8*q1q3*q2q3 + 4*pow(
      q1q3,2) + 8*q1q2*q2q3 + 8*q1q2*q1q3 + 4*pow(q1q2,2))*den(piq3*pjq3 +
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 +
      piq1*pjq2 + piq1*pjq1) * ( 8*pow(q2q3,2)*pipj - 8*piq3*pjq2*q2q3 - 4*
         piq3*pjq1*q2q3 - 8*piq2*pjq3*q2q3 - 4*piq2*pjq1*q2q3 - 4*piq1*pjq3*
         q2q3 - 4*piq1*pjq2*q2q3 - 16*piq1*pjq1*q2q3 + 4*d*piq1*pjq1*q2q3 );

      Sij3_3_ab +=  + den(q1q2)*den(4*pow(q2q3,2) + 8*q1q3*q2q3 + 4*pow(q1q3,2)
       + 8*q1q2*q2q3 + 8*q1q2*q1q3 + 4*pow(q1q2,2))*den(piq3*pjq3 + piq3*pjq2
       + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2
       + piq1*pjq1) * ( 8*q2q3*pipj - 4*q1q3*pipj - 4*piq3*pjq2 + 4*piq3*pjq1
          - 4*piq2*pjq3 + 8*piq2*pjq2 + 8*piq2*pjq1 + 4*piq1*pjq3 + 8*piq1*
         pjq2 + 2*d*q1q3*pipj - 2*d*piq3*pjq1 - 2*d*piq2*pjq1 - 2*d*piq1*pjq3
          - 2*d*piq1*pjq2 );

      Sij3_3_ab +=  + den(q1q3)*den(4*pow(q2q3,2) + 8*q1q3*q2q3 + 4*pow(q1q3,2)
       + 8*q1q2*q2q3 + 8*q1q2*q1q3 + 4*pow(q1q2,2))*den(piq3*pjq3 + piq3*pjq2
       + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2
       + piq1*pjq1) * ( 8*q2q3*pipj - 4*q1q2*pipj + 8*piq3*pjq3 - 4*piq3*pjq2
          + 8*piq3*pjq1 - 4*piq2*pjq3 + 4*piq2*pjq1 + 8*piq1*pjq3 + 4*piq1*
         pjq2 + 2*d*q1q2*pipj - 2*d*piq3*pjq1 - 2*d*piq2*pjq1 - 2*d*piq1*pjq3
          - 2*d*piq1*pjq2 );

      double dF = 3.;
      approx += - std::pow(gs, 6)*T_F/2*(C_F - T_F/dF)*Sij3_3_ab*std::real(M_ij[std::to_string(i) + std::to_string(j)]);

      double Sij3_3_nab =
       + den(4*q2q3)*den(piq3 + piq2)*den(pjq3 + pjq2)*den(piq3*pjq3 + piq3*
      pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*
      pjq2 + piq1*pjq1) * (  - 4*pow(pipj,2) );

      Sij3_3_nab +=  + den(pow(q2q3,2))*den(4*pow(q2q3,2) + 8*q1q3*q2q3 + 4*pow(
      q1q3,2) + 8*q1q2*q2q3 + 8*q1q2*q1q3 + 4*pow(q1q2,2))*den(piq3*pjq3 +
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 +
      piq1*pjq2 + piq1*pjq1) * ( 32*q1q2*q1q3*pipj );

      Sij3_3_nab +=  + den(pow(q2q3,2))*den(pjq3 + pjq2)*den(2*q2q3 + 2*q1q3 + 2*
      q1q2)*den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 +
      piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * (  - 4*pjq3*q1q2*pipj
          - 4*pjq2*q1q3*pipj + 4*piq3*pow(pjq2,2) + 4*piq2*pow(pjq3,2) + 4*
         piq1*pjq2*pjq3 );

      Sij3_3_nab +=  + den(pow(q2q3,2))*den(piq3 + piq2)*den(2*q2q3 + 2*q1q3 + 2*
      q1q2)*den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 +
      piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * (  - 4*piq3*q1q2*pipj
          + 4*pow(piq3,2)*pjq2 - 4*piq2*q1q3*pipj + 4*piq2*piq3*pjq1 + 4*pow(
         piq2,2)*pjq3 );

      Sij3_3_nab +=  + den(2*pow(q2q3,2))*den(piq3 + piq2)*den(pjq3 + pjq2)*den(
      piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 +
      piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * ( 2*piq3*pjq2*pipj + 2*piq2*pjq3*
         pipj );

      Sij3_3_nab +=  + den(4*pow(q2q3,2) + 8*q1q3*q2q3 + 4*pow(q1q3,2) + 8*q1q2*
      q2q3 + 8*q1q2*q1q3 + 4*pow(q1q2,2))*den(piq3*pjq3 + piq3*pjq2 + piq3*
      pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*
      pjq1) * ( 16*pipj - 4*d*pipj );

      Sij3_3_nab +=  + den(q1q3*q2q3)*den(pjq3 + pjq2)*den(2*q2q3 + 2*q1q3 + 2*
      q1q2)*den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 +
      piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * (  - 2*pjq3*q1q2*pipj
          + piq3*q1q2*mi2 + 2*piq3*pjq2*pjq3 + 2*piq2*pow(pjq3,2) + 2*piq2*
         pjq1*pjq3 + 2*piq1*pjq2*pjq3 );

      Sij3_3_nab +=  + den(q1q3*q2q3)*den(piq3 + piq2)*den(2*q2q3 + 2*q1q3 + 2*
      q1q2)*den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 +
      piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * ( pjq3*q1q2*mi2 - 2*
         piq3*q1q2*pipj + 2*pow(piq3,2)*pjq2 + 2*piq2*piq3*pjq3 + 2*piq2*piq3*
         pjq1 + 2*piq1*piq3*pjq2 );

      Sij3_3_nab +=  + den(q1q2*q2q3)*den(pjq3 + pjq2)*den(2*q2q3 + 2*q1q3 + 2*
      q1q2)*den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 +
      piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * (  - 2*pjq2*q1q3*pipj
          + 2*piq3*pow(pjq2,2) + 2*piq3*pjq1*pjq2 + piq2*q1q3*mi2 + 2*piq2*
         pjq2*pjq3 + 2*piq1*pjq2*pjq3 );

      Sij3_3_nab +=  + den(q1q2*q2q3)*den(piq3 + piq2)*den(2*q2q3 + 2*q1q3 + 2*
      q1q2)*den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 +
      piq2*pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * ( pjq2*q1q3*mi2 - 2*
         piq2*q1q3*pipj + 2*piq2*piq3*pjq2 + 2*piq2*piq3*pjq1 + 2*pow(piq2,2)*
         pjq3 + 2*piq1*piq2*pjq3 );

      Sij3_3_nab +=  + den(q1q2*q1q3)*den(4*pow(q2q3,2) + 8*q1q3*q2q3 + 4*pow(
      q1q3,2) + 8*q1q2*q2q3 + 8*q1q2*q1q3 + 4*pow(q1q2,2))*den(piq3*pjq3 +
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 +
      piq1*pjq2 + piq1*pjq1) * (  - 8*pow(q2q3,2)*pipj + 8*piq3*pjq2*q2q3 + 4*
         piq3*pjq1*q2q3 + 8*piq2*pjq3*q2q3 + 4*piq2*pjq1*q2q3 + 4*piq1*pjq3*
         q2q3 + 4*piq1*pjq2*q2q3 + 16*piq1*pjq1*q2q3 - 4*d*piq1*pjq1*q2q3 );

      Sij3_3_nab +=  + den(piq1)*den(4*q2q3)*den(pjq3 + pjq2)*den(piq3*pjq3 +
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 +
      piq1*pjq2 + piq1*pjq1) * ( 4*pow(pipj,2) );

      Sij3_3_nab +=  + den(piq1)*den(pow(q2q3,2))*den(2*q2q3 + 2*q1q3 + 2*q1q2)*
      den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*
      pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * ( 4*piq3*q1q2*pipj - 4*pow(
         piq3,2)*pjq2 + 4*piq2*q1q3*pipj - 4*piq2*piq3*pjq1 - 4*pow(piq2,2)*
         pjq3 );

      Sij3_3_nab +=  + den(piq1)*den(2*pow(q2q3,2))*den(pjq3 + pjq2)*den(piq3*
      pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*
      pjq3 + piq1*pjq2 + piq1*pjq1) * (  - 2*piq3*pjq2*pipj - 2*piq2*pjq3*pipj
          );

      Sij3_3_nab +=  + den(piq1)*den(q1q3*q2q3)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*
      den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*
      pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * (  - pjq3*q1q2*mi2 + 2*piq3*
         q1q2*pipj - 2*pow(piq3,2)*pjq2 - 2*piq2*piq3*pjq3 - 2*piq2*piq3*pjq1
          - 2*piq1*piq3*pjq2 );

      Sij3_3_nab +=  + den(piq1)*den(q1q2*q2q3)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*
      den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*
      pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * (  - pjq2*q1q3*mi2 + 2*piq2*
         q1q3*pipj - 2*piq2*piq3*pjq2 - 2*piq2*piq3*pjq1 - 2*pow(piq2,2)*pjq3
          - 2*piq1*piq2*pjq3 );

      Sij3_3_nab +=  + den(piq1)*den(pjq1)*den(4*q2q3)*den(piq3*pjq3 + piq3*pjq2
       + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2
       + piq1*pjq1) * (  - 4*pow(pipj,2) );

      Sij3_3_nab +=  + den(piq1)*den(pjq1)*den(2*pow(q2q3,2))*den(piq3*pjq3 +
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 +
      piq1*pjq2 + piq1*pjq1) * ( 2*piq3*pjq2*pipj + 2*piq2*pjq3*pipj );

      Sij3_3_nab +=  + den(piq1)*den(q1q2)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*den(piq3
      *pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1
      *pjq3 + piq1*pjq2 + piq1*pjq1) * ( pjq1*mi2 + 2*piq2*pipj );

      Sij3_3_nab +=  + den(piq1)*den(q1q3)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*den(piq3
      *pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1
      *pjq3 + piq1*pjq2 + piq1*pjq1) * ( pjq1*mi2 + 2*piq3*pipj );

      Sij3_3_nab +=  + den(piq1)*den(q2q3)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*den(piq3
      *pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1
      *pjq3 + piq1*pjq2 + piq1*pjq1) * (  - pjq3*mi2 - pjq2*mi2 + 2*pjq1*mi2
          + 4*piq3*pipj + 4*piq2*pipj - 4*piq1*pipj );

      Sij3_3_nab +=  + den(pjq1)*den(4*q2q3)*den(piq3 + piq2)*den(piq3*pjq3 +
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 +
      piq1*pjq2 + piq1*pjq1) * ( 4*pow(pipj,2) );

      Sij3_3_nab +=  + den(pjq1)*den(pow(q2q3,2))*den(2*q2q3 + 2*q1q3 + 2*q1q2)*
      den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*
      pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * ( 4*pjq3*q1q2*pipj + 4*pjq2*
         q1q3*pipj - 4*piq3*pow(pjq2,2) - 4*piq2*pow(pjq3,2) - 4*piq1*pjq2*
         pjq3 );

      Sij3_3_nab +=  + den(pjq1)*den(2*pow(q2q3,2))*den(piq3 + piq2)*den(piq3*
      pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*
      pjq3 + piq1*pjq2 + piq1*pjq1) * (  - 2*piq3*pjq2*pipj - 2*piq2*pjq3*pipj
          );

      Sij3_3_nab +=  + den(pjq1)*den(q1q3*q2q3)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*
      den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*
      pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * ( 2*pjq3*q1q2*pipj - piq3*
         q1q2*mi2 - 2*piq3*pjq2*pjq3 - 2*piq2*pow(pjq3,2) - 2*piq2*pjq1*pjq3
          - 2*piq1*pjq2*pjq3 );

      Sij3_3_nab +=  + den(pjq1)*den(q1q2*q2q3)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*
      den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*
      pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * ( 2*pjq2*q1q3*pipj - 2*piq3*
         pow(pjq2,2) - 2*piq3*pjq1*pjq2 - piq2*q1q3*mi2 - 2*piq2*pjq2*pjq3 - 2
         *piq1*pjq2*pjq3 );

      Sij3_3_nab +=  + den(pjq1)*den(q1q2)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*den(piq3
      *pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1
      *pjq3 + piq1*pjq2 + piq1*pjq1) * ( 2*pjq2*pipj + piq1*mi2 );

      Sij3_3_nab +=  + den(pjq1)*den(q1q3)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*den(piq3
      *pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1
      *pjq3 + piq1*pjq2 + piq1*pjq1) * ( 2*pjq3*pipj + piq1*mi2 );

      Sij3_3_nab +=  + den(pjq1)*den(q2q3)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*den(piq3
      *pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1
      *pjq3 + piq1*pjq2 + piq1*pjq1) * ( 4*pjq3*pipj + 4*pjq2*pipj - 4*pjq1*
         pipj - piq3*mi2 - piq2*mi2 + 2*piq1*mi2 );

      Sij3_3_nab +=  + den(q1q2)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*den(piq3*pjq3 +
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 +
      piq1*pjq2 + piq1*pjq1) * (  - 8*pipj );

      Sij3_3_nab +=  + den(q1q2)*den(4*pow(q2q3,2) + 8*q1q3*q2q3 + 4*pow(q1q3,2)
       + 8*q1q2*q2q3 + 8*q1q2*q1q3 + 4*pow(q1q2,2))*den(piq3*pjq3 + piq3*pjq2
       + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2
       + piq1*pjq1) * (  - 8*q2q3*pipj - 4*q1q3*pipj + 20*piq3*pjq2 + 8*piq3*
         pjq1 + 20*piq2*pjq3 - 8*piq2*pjq2 - 4*piq2*pjq1 + 8*piq1*pjq3 - 4*
         piq1*pjq2 + 24*piq1*pjq1 + 2*d*q1q3*pipj + 4*d*piq2*pjq1 + 4*d*piq1*
         pjq2 - 4*d*piq1*pjq1 );

      Sij3_3_nab +=  + den(q1q2)*den(pjq3 + pjq2)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*
      den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*
      pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * (  - 2*pjq2*pipj - piq1*mi2
          );

      Sij3_3_nab +=  + den(q1q2)*den(piq3 + piq2)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*
      den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*
      pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * (  - pjq1*mi2 - 2*piq2*pipj
          );

      Sij3_3_nab +=  + den(q1q2)*den(q2q3)*den(4*pow(q2q3,2) + 8*q1q3*q2q3 + 4*
      pow(q1q3,2) + 8*q1q2*q2q3 + 8*q1q2*q1q3 + 4*pow(q1q2,2))*den(piq3*pjq3
       + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3
       + piq1*pjq2 + piq1*pjq1) * ( 12*piq3*pjq2*q1q3 + 12*piq2*pjq3*q1q3 - 8*
         piq2*pjq2*q1q3 - 12*piq2*pjq1*q1q3 - 12*piq1*pjq2*q1q3 - 2*d*piq3*
         pjq2*q1q3 - 2*d*piq2*pjq3*q1q3 - 4*d*piq2*pjq2*q1q3 + 2*d*piq2*pjq1*
         q1q3 + 2*d*piq1*pjq2*q1q3 );

      Sij3_3_nab +=  + den(q1q3)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*den(piq3*pjq3 +
      piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 +
      piq1*pjq2 + piq1*pjq1) * (  - 8*pipj );

      Sij3_3_nab +=  + den(q1q3)*den(4*pow(q2q3,2) + 8*q1q3*q2q3 + 4*pow(q1q3,2)
       + 8*q1q2*q2q3 + 8*q1q2*q1q3 + 4*pow(q1q2,2))*den(piq3*pjq3 + piq3*pjq2
       + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2
       + piq1*pjq1) * (  - 8*q2q3*pipj - 4*q1q2*pipj - 8*piq3*pjq3 + 20*piq3*
         pjq2 - 4*piq3*pjq1 + 20*piq2*pjq3 + 8*piq2*pjq1 - 4*piq1*pjq3 + 8*
         piq1*pjq2 + 24*piq1*pjq1 + 2*d*q1q2*pipj + 4*d*piq3*pjq1 + 4*d*piq1*
         pjq3 - 4*d*piq1*pjq1 );

      Sij3_3_nab +=  + den(q1q3)*den(pjq3 + pjq2)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*
      den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*
      pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * (  - 2*pjq3*pipj - piq1*mi2
          );

      Sij3_3_nab +=  + den(q1q3)*den(piq3 + piq2)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*
      den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*
      pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * (  - pjq1*mi2 - 2*piq3*pipj
          );

      Sij3_3_nab +=  + den(q1q3)*den(q2q3)*den(4*pow(q2q3,2) + 8*q1q3*q2q3 + 4*
      pow(q1q3,2) + 8*q1q2*q2q3 + 8*q1q2*q1q3 + 4*pow(q1q2,2))*den(piq3*pjq3
       + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3
       + piq1*pjq2 + piq1*pjq1) * (  - 8*piq3*pjq3*q1q2 + 12*piq3*pjq2*q1q2 -
         12*piq3*pjq1*q1q2 + 12*piq2*pjq3*q1q2 - 12*piq1*pjq3*q1q2 - 4*d*piq3*
         pjq3*q1q2 - 2*d*piq3*pjq2*q1q2 + 2*d*piq3*pjq1*q1q2 - 2*d*piq2*pjq3*
         q1q2 + 2*d*piq1*pjq3*q1q2 );

      Sij3_3_nab +=  + den(q2q3)*den(4*pow(q2q3,2) + 8*q1q3*q2q3 + 4*pow(q1q3,2)
       + 8*q1q2*q2q3 + 8*q1q2*q1q3 + 4*pow(q1q2,2))*den(piq3*pjq3 + piq3*pjq2
       + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*pjq1 + piq1*pjq3 + piq1*pjq2
       + piq1*pjq1) * ( 16*q1q3*pipj + 16*q1q2*pipj - 16*piq3*pjq3 + 16*piq3*
         pjq2 - 12*piq3*pjq1 + 16*piq2*pjq3 - 16*piq2*pjq2 - 12*piq2*pjq1 - 12
         *piq1*pjq3 - 12*piq1*pjq2 + 8*piq1*pjq1 + 2*d*piq3*pjq1 + 2*d*piq2*
         pjq1 + 2*d*piq1*pjq3 + 2*d*piq1*pjq2 - 4*d*piq1*pjq1 );

      Sij3_3_nab +=  + den(q2q3)*den(pjq3 + pjq2)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*
      den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*
      pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * (  - 4*pjq3*pipj - 4*pjq2*
         pipj + 4*pjq1*pipj + piq3*mi2 + piq2*mi2 - 2*piq1*mi2 );

      Sij3_3_nab +=  + den(q2q3)*den(piq3 + piq2)*den(2*q2q3 + 2*q1q3 + 2*q1q2)*
      den(piq3*pjq3 + piq3*pjq2 + piq3*pjq1 + piq2*pjq3 + piq2*pjq2 + piq2*
      pjq1 + piq1*pjq3 + piq1*pjq2 + piq1*pjq1) * ( pjq3*mi2 + pjq2*mi2 - 2*
         pjq1*mi2 - 4*piq3*pipj - 4*piq2*pipj + 4*piq1*pipj );

      approx += -std::pow(gs, 6)*T_F*C_A/4.*Sij3_3_nab*std::real(M_ij[std::to_string(i) + std::to_string(j)]);
    }
  }
  return approx;
}

double Sij_2(double *pi, double *pj, double *q1, double *q2) {
  double q1q2 = minkovski(q1, q2);
  double q12[4];
  add_arr(q1, q2, q12, 4);
  double piq1 = minkovski(pi, q1);
  double piq2 = minkovski(pi, q2);
  double piq12 = piq1 + piq2;
  double pjq1 = minkovski(pj, q1);
  double pjq2 = minkovski(pj, q2);
  double pjq12 = pjq1 + pjq2;
  double pipj = minkovski(pi, pj);
  double Sij12_2 = 1./q1q2/q1q2*(piq1*pjq2 + piq2*pjq1)/piq12/pjq12
        - pipj*pipj/2./piq1/pjq2/piq2/pjq1*(2. - (piq1*pjq2 + piq2*pjq1)/piq12/pjq12)
        + pipj/2./q1q2*(2./piq1/pjq2 + 2./pjq1/piq2 - 1./piq12/pjq12*(4. + std::pow(piq1*pjq2 + piq2*pjq1, 2)/piq1/pjq2/piq2/pjq1));
  return Sij12_2;
}

double soft_ggg_squared(double* pp_full, std::unordered_map<std::string, double> M_ij,
                                         std::unordered_map<std::string, double> M_ijkl,
                                         std::unordered_map<std::string, double> dM_ijk,
                                         std::unordered_map<std::string, double> M_ijklab, amplitude& A) {
  double q1[4], q2[4], q3[4], q12[4], q13[4], q23[4], q123[4];
  part(pp_full, q1, A.process.size()*4, A.process.size()*4 + 4);
  part(pp_full, q2, A.process.size()*4 + 4, A.process.size()*4 + 8);
  part(pp_full, q3, A.process.size()*4 + 8, A.process.size()*4 + 12);
  add_arr(q1, q2, q12, 4);
  add_arr(q1, q3, q13, 4);
  add_arr(q2, q3, q23, 4);
  add_arr(q12, q3, q123, 4);
  double q1q2 = minkovski(q1, q2);
  double q1q3 = minkovski(q1, q3);
  double q2q3 = minkovski(q2, q3);

  double approx = 0.;
  for(int i = 0; i < A.process.size(); i++) {
    if(A.process[i] == 1) continue;
    double pi[4];
    double piq1 = minkovski(pi, q1);
    double piq2 = minkovski(pi, q2);
    double piq3 = minkovski(pi, q3);
    double piq12 = minkovski(pi, q12);
    double piq13 = minkovski(pi, q13);
    double piq23 = minkovski(pi, q23);
    double piq123 = minkovski(pi, q123);
    part(pp_full, pi, i*4, i*4 + 4);
    for(int j = 0; j < A.process.size(); j++) {
      if(A.process[j] == 1) continue;
      double pj[4];
      part(pp_full, pj, 4*j, 4*j + 4);
      double pjq1 = minkovski(pj, q1);
      double pjq2 = minkovski(pj, q2);
      double pjq12 = minkovski(pj, q12);
      double pjq13 = minkovski(pj, q13);
      double pjq23 = minkovski(pj, q23);
      double pipj = minkovski(pi, pj);
      double mi2 = minkovski(pi, pi);
      double pjq3 = minkovski(pj, q3);
      double q123q123 = minkovski(q123, q123);

      double wij1_1 = 2.*pipj/piq1/pjq1;
      double wij2_1 = 2.*pipj/piq2/pjq2;
      double wij3_1 = 2.*pipj/piq3/pjq3;

      double wij12_2 = Sij_2(pi, pj, q1, q2) + Sij_2(pj, pi, q1, q2) - Sij_2(pi, pi, q1, q2) - Sij_2(pj, pj, q1, q2);
      double wij13_2 = Sij_2(pi, pj, q1, q3) + Sij_2(pj, pi, q1, q3) - Sij_2(pi, pi, q1, q3) - Sij_2(pj, pj, q1, q3);
      double wij23_2 = Sij_2(pi, pj, q2, q3) + Sij_2(pj, pi, q2, q3) - Sij_2(pi, pi, q2, q3) - Sij_2(pj, pj, q2, q3);

      for(int k = 0; k < A.process.size(); k++) {
        if(A.process[k] == 1) continue;
        double pk[4];
        part(pp_full, pk, 4*k, 4*k + 4);
        double pkq1 = minkovski(pk, q1);
        double pkq2 = minkovski(pk, q2);
        double pkq3 = minkovski(pk, q3);
        double pipk = minkovski(pi, pk);
        double pjpk = minkovski(pj, pk);
        for(int l = 0; l < A.process.size(); l++) {
          if(A.process[l] == 1) continue;
          double pl[4];
          part(pp_full, pl, 4*l, 4*l + 4);
          double plq1 = minkovski(pl, q1);
          double plq2 = minkovski(pl, q2);
          double plq3 = minkovski(pl, q3);
          double pkpl = minkovski(pk, pl);
          double wkl1_1 = 2.*pkpl/pkq1/plq1;
          double wkl2_1 = 2.*pkpl/pkq2/plq2;
          double wkl3_1 = 2.*pkpl/pkq3/plq3;

          for(int a = 0; a < A.process.size(); a++) {
            if(A.process[a] == 1) continue;
            double pa[4];
            part(pp_full, pa, 4*a, 4*a + 4);
            double paq1 = minkovski(pa, q1);
            double paq2 = minkovski(pa, q2);
            double paq3 = minkovski(pa, q3);
            for(int b = 0; b < A.process.size(); b++) {
              if(A.process[b] == 1) continue;
              double pb[4];
              part(pp_full, pb, 4*b, 4*b + 4);
              double pbq3 = minkovski(pb, q3);
              double papb = minkovski(pa, pb);
              double wab3_1 = 3.*papb/paq3/pbq3;
              approx += -1./6.*wij1_1*wkl2_1*wab3_1*(std::real(M_ijklab[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l) + std::to_string(a) + std::to_string(b)])
                                                   + std::real(M_ijklab[std::to_string(i) + std::to_string(j) + std::to_string(a) + std::to_string(b) + std::to_string(k) + std::to_string(l)])
                                                   + std::real(M_ijklab[std::to_string(k) + std::to_string(l) + std::to_string(i) + std::to_string(j) + std::to_string(a) + std::to_string(b)])
                                                   + std::real(M_ijklab[std::to_string(k) + std::to_string(l) + std::to_string(a) + std::to_string(b) + std::to_string(i) + std::to_string(j)])
                                                   + std::real(M_ijklab[std::to_string(a) + std::to_string(b) + std::to_string(k) + std::to_string(l) + std::to_string(i) + std::to_string(j)])
                                                   + std::real(M_ijklab[std::to_string(a) + std::to_string(b) + std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]));

            }
          }

          approx += 1./2.*C_A/8.*(wij12_2*wkl3_1 + wij13_2*wkl2_1 + wij23_2*wkl1_1)*(std::real(M_ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)])
                                                                            + std::real(M_ijkl[std::to_string(k) + std::to_string(l) + std::to_string(i) + std::to_string(j)]));

        }
      }
    }
  }
  return approx*std::pow(gs, 6);
}

double soft_g_squared_1l(double *pp_full, std::unordered_map<std::string, double> M0_ij, std::unordered_map<std::string, double> fM_ijk,
  std::unordered_map<std::string, double> M1_ij, amplitude& A) {
  double q_arr[4];
  part(pp_full, q_arr, A.process.size()*4, A.process.size()*4 + 4);
  LV<double> q(q_arr);
  double approx = 0;
  for(int i = 0; i < A.process.size(); i++) {
    if(A.process[i] == 1) continue;
    double pi_arr[4];
    part(pp_full, pi_arr, i*4, i*4 + 4);
    LV<double> pi(pi_arr);
    if(i < 2) pi = (-1.)*pi;
    for(int j = 0; j < A.process.size(); j++) {
      if(A.process[j] == 1) continue;
      if(i == j) continue;
      double pj_arr[4];
      part(pp_full, pj_arr, j*4, j*4 + 4);
      LV<double> pj(pj_arr);
      if(j < 2) pj = (-1.)*pj;
      approx += -gs*gs*(j1(pi, q)*j1(pj, q))*std::real(M1_ij[std::to_string(i) + std::to_string(j)]);
      approx += gs*gs*gs*gs*C_A*2.*std::real(gamma11(pi, pj, q)*j1(pi, q))*std::real(M0_ij[std::to_string(i) + std::to_string(j)]);
      for(int k = 0; k < A.process.size(); k++) {
        if(A.process[k] == 1) continue;
        if((k == i) or (k == j)) continue;
        double pk_arr[4];
        part(pp_full, pk_arr, k*4, 4*k + 4);
        LV<double> pk(pk_arr);
        if(k < 2) pk = (-1.)*pk;
        approx += gs*gs*gs*gs*2.*std::imag(gamma11(pj, pk, q)*j1(pi, q))*std::real(fM_ijk[std::to_string(i) + std::to_string(j) + std::to_string(k)]);
      }
    }
  }
  return approx;
}

double soft_gg_squared_1l(double*pp_full, std::unordered_map<std::string, double> M0ij,
    std::unordered_map<std::string, double> M1ij, std::unordered_map<std::string, double> M0ijkl,
    std::unordered_map<std::string, double> M1ijkl, std::unordered_map<std::string, double> M0ijkla,
    std::unordered_map<std::string, double> Q0ijkl, amplitude& A, int nl) {
  double q1_arr[4], q2_arr[4];
  part(pp_full, q1_arr, A.process.size()*4, A.process.size()*4 + 4);
  part(pp_full, q2_arr, A.process.size()*4 + 4, A.process.size()*4 + 8);
  LV<double> q1(q1_arr);
  LV<double> q2(q2_arr);
  double approx = 0;
  for(int i = 0; i < A.process.size(); i++) {
    if(A.process[i] == 1) continue;
    double pi_arr[4];
    part(pp_full, pi_arr, i*4, i*4 + 4);
    LV<double> pi(pi_arr);
    for(int j = 0; j < A.process.size(); j++) {
      if(A.process[j] == 1) continue;
      double pj_arr[4];
      part(pp_full, pj_arr, 4*j, 4*j + 4);
      LV<double> pj(pj_arr);

      for(int k = 0; k < A.process.size(); k++) {
        if(A.process[k] == 1) continue;
        double pk_arr[4];
        part(pp_full, pk_arr, 4*k, 4*k + 4);
        LV<double> pk(pk_arr);
        for(int l = 0; l < A.process.size(); l++) {
          if(A.process[l] == 1) continue;
          double pl_arr[4];
          part(pp_full, pl_arr, 4*l, 4*l + 4);
          LV<double> pl(pl_arr);

          // Reducible part (tree level current)
          if((i != j) and (k != l))
            approx += std::pow(gs, 4)*(j1(pi, q1)*j1(pj, q1))*(j1(pk, q2)*j1(pl, q2))*(std::real(M1ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)])
                                                                                     + std::real(M1ijkl[std::to_string(k) + std::to_string(l) + std::to_string(i) + std::to_string(j)]))/2. ;

          // Independent emmision
          if((i != j) and (k != l)) {
            approx += std::pow(gs, 6)*(((-1.)*j1(pk, q1)*j1(pl, q1))*C_A*2.*std::real(gamma11(pi, pj, q2)*j1(pi, q2))
                                      +((-1.)*j1(pk, q2)*j1(pl, q2))*C_A*2.*std::real(gamma11(pi, pj, q1)*j1(pi, q1)))
                        *(std::real(M0ijkl[std::to_string(k) + std::to_string(l) + std::to_string(i) + std::to_string(j)])
                        + std::real(M0ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]))/2.;
            for(int a = 0; a < A.process.size(); a++) {
              if(A.process[a] == 1) continue;
              double pa_arr[4];
              part(pp_full, pa_arr, 4*a, 4*a + a);
              LV<double> pa(pa_arr);
              if((a != k) and (a != l)) {
                //approx += std::pow(gs, 6)*std::real(M0ijkla[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l) + std::to_string(a)])/2.
                //  *(((-1.)*j1(pi, q1)*j1(pj, q1))*std::imag(gamma11(pl, pa, q2)*j1(pk, q2))
                //  + ((-1.)*j1(pi, q2)*j1(pj, q2))*std::imag(gamma11(pl, pa, q1)*j1(pk, q1)));
              }
            }
          }

          // Quadupole operator contribution
          approx += std::pow(gs, 6)*(-1.)*2.*std::real((gamma11(pk, pl, q1)*j1(pi, q1))*(j1(pj, q2)*j1(pi, q2))/4.
                                                      +(gamma11(pk, pl, q2)*j1(pi, q2))*(j1(pj, q1)*j1(pi, q1))/4.
                                                      +(gamma11(pj, pk, q1)*(gamma20(pi, q1, q2)*j1(pl, q2)))/4.
                                                      +(gamma11(pj, pk, q2)*(gamma20(pi, q2, q1)*j1(pl, q1)))/4.)
                    *std::real(Q0ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]);
          if(i != j)
            approx += std::pow(gs, 6)*2.*std::real(j1(pl, q1)*(gamma21(pi, pj, q1, q2, nl)*j1(pk, q2))
                                                        +j1(pl, q2)*(gamma21(pi, pj, q2, q1, nl)*j1(pk, q1)))/8.
                    *std::real(Q0ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]);



        }
      }
      // Irreducible part (tree level current)
      approx += std::pow(gs, 4)*C_A*std::real(M1ij[std::to_string(i) + std::to_string(j)])
          *((gamma20(pi, q1, q2)*gamma20(pj, q1, q2).transpose()).trace()
            +j1(pi, q1)*(gamma20(pi, q1, q2)*j1(pj, q2))
            -j1(pi, q1)*(gamma20(pj, q1, q2)*j1(pj, q2))
            -1./2.*(j1(pi, q1)*j1(pi, q1))*(j1(pi, q2)*j1(pj, q2))
            -1./2.*(j1(pj, q2)*j1(pj, q2))*(j1(pi, q1)*j1(pj, q1))
            +3./4.*(j1(pi, q1)*j1(pj, q1))*(j1(pi, q2)*j1(pj, q2)));

      if(i != j) {
        approx += std::pow(gs, 6)*std::pow(C_A, 2)*2.*std::real(j1(pi, q1)*(gamma21(pi, pj, q1, q2)*j1(pi, q2))/8.
                                                               +j1(pi, q2)*(gamma21(pi, pj, q2, q1)*j1(pi, q1))/8.
                                                               -j1(pi, q1)*(gamma21(pi, pj, q1, q2)*j1(pj, q2))/8.
                                                               -j1(pi, q2)*(gamma21(pi, pj, q2, q1)*j1(pj, q1))/8.
                                                               -(gamma11(pi, pj, q1)*j1(pi, q1))*(j1(pi, q2)*j1(pj, q2))
                                                               -(gamma11(pi, pj, q2)*j1(pi, q2))*(j1(pi, q1)*j1(pj, q1))
                                                               +(gamma20(pi, q1, q2)*gamma21(pi, pj, q1, q2).transpose()).trace()/4.
                                                               +(gamma20(pi, q2, q1)*gamma21(pi, pj, q2, q1).transpose()).trace()/4.)
                *std::real(M0ij[std::to_string(i) + std::to_string(j)]);
      }
    }
  }
  return approx;
}

double soft_gg_squared_1l_quad(double*pp_full, std::unordered_map<std::string, double> M0ij,
    std::unordered_map<std::string, double> M1ij, std::unordered_map<std::string, double> M0ijkl,
    std::unordered_map<std::string, double> M1ijkl, std::unordered_map<std::string, double> M0ijkla,
    std::unordered_map<std::string, double> Q0ijkl, amplitude& A, int nl) {
  double q1_arr[4], q2_arr[4];
  part(pp_full, q1_arr, A.process.size()*4, A.process.size()*4 + 4);
  part(pp_full, q2_arr, A.process.size()*4 + 4, A.process.size()*4 + 8);
  LV<double> q1(q1_arr);
  LV<double> q2(q2_arr);
  double approx = 0;
  for(int i = 0; i < A.process.size(); i++) {
    if(A.process[i] == 1) continue;
    double pi_arr[4];
    part(pp_full, pi_arr, i*4, i*4 + 4);
    LV<double> pi(pi_arr);
    for(int j = 0; j < A.process.size(); j++) {
      if(A.process[j] == 1) continue;
      double pj_arr[4];
      part(pp_full, pj_arr, 4*j, 4*j + 4);
      LV<double> pj(pj_arr);
      for(int k = 0; k < A.process.size(); k++) {
        if(A.process[k] == 1) continue;
        double pk_arr[4];
        part(pp_full, pk_arr, 4*k, 4*k + 4);
        LV<double> pk(pk_arr);
        for(int l = 0; l < A.process.size(); l++) {
          if(A.process[l] == 1) continue;
          double pl_arr[4];
          part(pp_full, pl_arr, 4*l, 4*l + 4);
          LV<double> pl(pl_arr);
          // Quadupole operator contribution
          approx += std::pow(gs, 6)*(-1.)*2.*std::real((gamma11(pk, pl, q1)*j1(pi, q1))*(j1(pj, q2)*j1(pi, q2))/4.
                                                      +(gamma11(pk, pl, q2)*j1(pi, q2))*(j1(pj, q1)*j1(pi, q1))/4.
                                                      +(gamma11(pj, pk, q1)*(gamma20(pi, q1, q2)*j1(pl, q2)))/4.
                                                      +(gamma11(pj, pk, q2)*(gamma20(pi, q2, q1)*j1(pl, q1)))/4.)
                    *std::real(Q0ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]);
          if(i != j)
            approx += std::pow(gs, 6)*2.*std::real(j1(pl, q1)*(gamma21(pi, pj, q1, q2, nl)*j1(pk, q2))
                                                        +j1(pl, q2)*(gamma21(pi, pj, q2, q1, nl)*j1(pk, q1)))/8.
                    *std::real(Q0ijkl[std::to_string(i) + std::to_string(j) + std::to_string(k) + std::to_string(l)]);



        }
      }
    }
  }
  return approx;
}

