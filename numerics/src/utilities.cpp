#include "utilities.hpp"

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

double li2(double x) {
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
