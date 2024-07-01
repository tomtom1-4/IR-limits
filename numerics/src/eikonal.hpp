#ifndef EIKONAL_HEADER
#define EIKONAL_HEADER

#include <complex>
#include <unordered_map>
#include "amplitude.hpp"
#include "color.hpp"
#include "dirac.hpp"

std::unordered_map<std::string, std::complex<double>> J_g_eikonal(double *pp, double *q, amplitude& A);

std::unordered_map<std::string, std::complex<double>> J_gg_eikonal(double *pp, double *q1, double *q2, amplitude& A);

std::unordered_map<std::string, std::complex<double>> J_qq_eikonal(double *pp, double *q1, double *q2, amplitude& A);

std::unordered_map<std::string, std::complex<double>> J_qqg_eikonal(double *pp, double *q1, double *q2, double *q3, amplitude& A);

std::complex<double> soft_g(double *pp_full, std::unordered_map<std::string, std::complex<double>> J,
    std::unordered_map<std::string, std::complex<double>> M,
    std::vector<int> helicities_full, std::vector<int> colors_full, amplitude& A);

std::complex<double> soft_gg(double *pp_full, std::unordered_map<std::string, std::complex<double>> J1,
  std::unordered_map<std::string, std::complex<double>> J2, std::unordered_map<std::string, std::complex<double>> J12,
    std::unordered_map<std::string, std::complex<double>> M,
    std::vector<int> helicities_full, std::vector<int> colors_full, amplitude& A);

std::complex<double> soft_qq(double *pp_full, std::unordered_map<std::string, std::complex<double>> J12,
    std::unordered_map<std::string, std::complex<double>> M,
    std::vector<int> helicities_full, std::vector<int> colors_full, amplitude& A);


std::complex<double> double_soft_tree_qq(double *pp_full, std::unordered_map<std::string, std::complex<double>> M,
                     std::vector<int> helicities_full, std::vector<int> colors_full, amplitude& A);

std::complex<double> soft_qqg(double *pp_full, std::unordered_map<std::string, std::complex<double>> J12, std::unordered_map<std::string, std::complex<double>> J3,
    std::unordered_map<std::string, std::complex<double>> J123, std::unordered_map<std::string, std::complex<double>> M,
    std::vector<int> helicities_full, std::vector<int> colors_full, amplitude& A);

double soft_g_squared(double *pp_full, std::unordered_map<std::string, std::complex<double>> M_ij, amplitude& A);

double soft_qq_squared(double *pp_full, std::unordered_map<std::string, std::complex<double>> M_ij, amplitude& A);

double soft_gg_squared(double*pp_full, std::unordered_map<std::string, std::complex<double>> Mij, std::unordered_map<std::string, std::complex<double>> Mijkl, amplitude& A);

double soft_gqq_squared(double* pp_full, std::unordered_map<std::string, std::complex<double>> M_ij,
                                         std::unordered_map<std::string, std::complex<double>> M_ijkl,
                                         std::unordered_map<std::string, std::complex<double>> dM_ijk, amplitude& A);

#endif