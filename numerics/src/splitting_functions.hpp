#ifndef SPLITTING_FUNCTIONS_HEADER
#define SPLITTING_FUNCTIONS_HEADER

#include <complex>
#include <unordered_map>
#include <vector>
#include <type_traits>
#include "amplitude.hpp"
#include "color.hpp"
#include "dirac.hpp"
#include "Lorentz.hpp"
#include "utilities.hpp"

std::vector<std::vector<std::complex<double>>> P0_qg(double z, LV<double> p, LV<double> kperp);

std::vector<std::vector<std::complex<double>>> P0_gg(double z, LV<double> p, LV<double> kperp);

std::vector<std::vector<std::complex<double>>> P0_qq(double z, LV<double> p, LV<double> kperp);

std::vector<std::vector<std::complex<double>>> P1_qg(double z, LV<double> p, LV<double> kperp, double s12, double muR);

std::vector<std::vector<std::complex<double>>> P1_gg(double z, LV<double> p, LV<double> kperp, double s12, double muR);

std::vector<std::vector<std::complex<double>>> P1_qq(double z, LV<double> p, LV<double> kperp, double s12, double muR);

std::vector<std::vector<std::complex<double>>> P0_qqPqPbar(LV<double> p1, LV<double> p2, LV<double> p3);

std::vector<std::vector<std::complex<double>>> P0_qgg(LV<double> p1, LV<double> p2, LV<double> p3);

std::vector<std::vector<std::complex<double>>> P0_qqqbar(LV<double> p1, LV<double> p2, LV<double> p3);

std::vector<std::vector<std::complex<double>>> P0_gqqbar(LV<double> p1, LV<double> p2, LV<double> p3);

std::vector<std::vector<std::complex<double>>> P0_ggg(LV<double> p1, LV<double> p2, LV<double> p3);

#endif