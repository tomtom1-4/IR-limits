#ifndef PARAMS_HEADER
#define PARAMS_HEADER

#include <complex>

// Kinematic parameters
static double COM = 1000;
static double Emin = 1.e-1;
static double Emax = 10;
static int events = 1;
static double exp_growth = pow(Emax/Emin, 1/double(events));

// Physical parameters
static double gs = 0.1;
static double n_f = 5.; //????
static double mu = 100;
static double e = 1.;
static double costhetaW = 8.0357973609878002E+01/9.1153480619183000E+01;
static double C_F = 4./3.;
static double C_A = 3.;
static double T_F = 1./2.;
static double NC = 3;

// Unphysical parameters
static std::complex<double> I(0, 1);
#endif