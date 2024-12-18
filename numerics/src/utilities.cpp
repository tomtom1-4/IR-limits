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

const unsigned nLi3dat = 28;
const double li3dat[] = {
  +9.9999999999999999e-1,
  -3.7499999999999999e-1,
  +7.8703703703703703e-2,
  -8.6805555555555555e-3,
  +1.2962962962962962e-4,
  +8.1018518518518518e-5,
  -3.4193571608537594e-6,
  -1.3286564625850340e-6,
  +8.6608717561098513e-8,
  +2.5260875955320399e-8,
  -2.1446944683640647e-9,
  -5.1401106220129789e-10,
  +5.2495821146008234e-11,
  +1.0887754406636337e-11,
  -1.2779396094489561e-12,
  -2.3698241773096296e-13,
  +3.1043578877652041e-14,
  +5.2617586301809793e-15,
  -7.5384794808586713e-16,
  -1.1862322625901778e-16,
  +1.8316963091716708e-17,
  +2.7068173862683523e-18,
  -4.4551487031372274e-19,
  -6.2374702753511882e-20,
  +1.0819745927876939e-20,
  +1.4470579164468938e-21,
  -2.4357087950089783e-22,
  -3.1858705075010159e-23};

// Coefficients for the quadrilogarithm for double arguments in the range
// [-1,1/2]. A Chebychev economization has been performed.

const unsigned nLi4dat = 13;
const double li4dat[] =
  {
    +1.0000000000000000e0,
    -4.3749999999999811e-1,
    +1.1651234567901286e-1,
    -1.9820601851914722e-2,
    +1.9279320987537474e-3,
    -3.1057097978981336e-5,
    -1.5624008998892365e-5,
    +8.4850766609360565e-7,
    +2.2909557483519552e-7,
    -2.1818241515288884e-8,
    -3.8812086294383839e-9,
    +5.2269024271386895e-10,
    +6.7350275938466138e-11,
  };

const unsigned nS12dat = 14;
static const double s12dat[] =
  {
    +2.50000000000000000000000000000000017e-1,
    +1.0416666666666666666666666666663283e-2,
    -1.15740740740740740740740740506743412e-4,
    +2.06679894179894179894179078508253495e-6,
    -4.13359788359788359786678514412975012e-8,
    +8.69864874494504121896219505207926127e-10,
    -1.88721076381695982470401421658486232e-11,
    +4.18204266583767160404184289236370136e-13,
    -9.41577860030863104679826142453384509e-15,
    +2.1465154948596402604813220890299258e-16,
    -4.94287892802404593212536273410514027e-18,
    +1.14763628191847339520489420246104811e-19,
    -2.6757467386214950316359096425484803e-21,
    +5.81053106716134930692373192778044681e-23
  };

// Coefficients for s13 for real arguments in the range [0,1/2].
// A Chebychev economization has been performed. Only the coefficients of
// odd powers, y^(2n+1), are specified, while the y^4 term is added
// separately. y = -Log[1-x]

const unsigned nS13dat = 6;
const double s13dat[] =
  {
    +5.5555555555557924e-2,
    +2.7777777777244863e-3,
    -3.3068782539740705e-5,
    +6.1238217344237218e-7,
    -1.2518672557852010e-8,
    +2.5732481317698610e-10
  };

/**************************************************************************
 *                                                                        *
 * s12                                                                    *
 *                                                                        *
 **************************************************************************/

// helper function only valid in [0,1/2]

double s12(double x)
{
  double result = 0.0;
  const double y = -log(1.-x);
  const double y2 = y*y;
  for (int i = nS12dat-1; i > 0; --i)
    {
      result += s12dat[i];
      result *= y2;
    }
  result += s12dat[0]-y/double(12);
  result *= y2;

  return result;
}

/**************************************************************************
 *                                                                        *
 * s13                                                                    *
 *                                                                        *
 **************************************************************************/

// helper function only valid in [0,1/2]

double s13(double x)
{
  double result = 0.0;
  const double y = -log(1.-x);
  const double y2 = y*y;
  const double y3 = y2*y;
  for (int i = nS13dat-1; i > 0; --i)
    {
      result += s13dat[i];
      result *= y2;
    }
  result += s13dat[0]-y/48.;
  result *= y3;

  return result;
}

/****************************************************************************
 *                                                                          *
 * li2                                                                      *
 *                                                                          *
 ****************************************************************************/

double li2(double x)
{
  if (x > 1.) std::cout << "out-of-range argument" << std::endl;
  if (x == 1.) return zeta2;
  if (x > 0.5) return -li2(1.-x)-log(x)*log(1.-x)+zeta2;
  if (x < -1.)
    {
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

/****************************************************************************
 *                                                                          *
 * li3                                                                      *
 *                                                                          *
 ****************************************************************************/

double li3(double x)
{
  if (x > 1.) std::cout << "out-of-range argument" << std::endl;
  if (x == 1.) return zeta3;
  if (x > 0.5)
    {
      const double ln = log(x);
      return -s12(1.-x)+ln*(li2(x)+0.5*log(1.-x)*ln)+zeta3;
    }
  if (x < -1.)
    {
      const double ln = log(-x);
      return li3(1/x)-ln*(zeta2+ln*ln/double(6));
    }

  double result = 0.0;
  const double y = -log(1.-x);
  for (int i = nLi3dat-1; i >= 0; --i)
    {
      result += li3dat[i];
      result *= y;
    }

  return result;
}

/****************************************************************************
 *                                                                          *
 * li4                                                                      *
 *                                                                          *
 ****************************************************************************/

double li4(double x)
{
  if (x > 1.) std::cout << "out-of-range argument" << std::endl;
  if (x == 1.) return zeta4;
  if (x > 0.5)
    {
      const double ln = log(x);
      return zeta4-s13(1.-x)
        -ln*(s12(1.-x)-0.5*ln*(zeta2-li2(1.-x)-log(1.-x)*ln/3.)-zeta3);
    }
  if (x < -1.)
    {
      const double ln = log(-x);
      const double ln2 = ln*ln;
      return -li4(1./x)-0.5*ln2*(zeta2+ln2/12.)-1.75*zeta4;
    }

  double result = 0.0;
  const double y = -log(1.-x);
  for (int i = nLi4dat-1; i >= 0; --i)
    {
      result += li4dat[i];
      result *= y;
    }

  return result;
}
