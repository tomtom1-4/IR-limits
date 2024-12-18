#ifndef UTILS_HEADER
#define UTILS_HEADER

#include <cmath>
#include <iostream>

static double zeta2 = M_PI*M_PI/6.;
static double zeta3 = 1.202056903159594;
static double zeta4 = 1.082323233711138;

/// Dilogarithm for double argument in (-Infty,1]
double li2(double x);

/// Trilogarithm for double argument in (-Infty,1]
double li3(double x);

/// Quadrilogarithm for double argument in (-Infty,1]
double li4(double x);

#endif