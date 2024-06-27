#ifndef ARRAY_MANIPULATION_HEADER
#define ARRAY_MANIPULATION_HEADER

#include <vector>
#include <iostream>
#include <complex>
#include "params.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for generall array manipulations
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class A>
static A sum_of_arr(A *arr, int length) {
    A res = 0;
    for (int i = 0; i < length; i++) {
        res += arr[i];
    }
    return res;
}


static bool arrayeq(int* arr1, int*arr2, int length) {
    for (int i = 0; i < length; i++) {
        if (arr1[i] != arr2[i]) {
            return false;
        }
    }
    return true;
}

template <class A>
static void vector2arr(std::vector<A> v, A *a) {
    for (int i = 0; i < v.size(); i++) {
        a[i] = v[i];
    }
}

static void print_vec_int(std::vector<int> a) {
    for (int i = 0; i < a.size(); i++) {
        std::cout << a[i] << ", ";
    }
}

static void print_vec_double(std::vector<double> a) {
    for (int i = 0; i < a.size(); i++) {
        std::cout << a[i] << ", ";
    }
}

static void print_int_arr(int *arr, int length) {
    for (int i = 0; i < length; i++ ) {
        std::cout << arr[i] << ", ";
    }
}

static int delta(int a, int b) {
    if (a == b) {
        return 1;
    }
    else {
        return 0;
    }
}

// Get a part of an other array
static void part(double *arr, double *part_of_arr, int index1, int index2) {
    // Updates part_of_arr to be arr[index1:index2]
    for (int i = 0; i < index2 - index1; i++) {
        part_of_arr[i] = arr[i + index1] ;
    }
    return ;
}

// Add two arrays
static void add_arr(double *arr1, double *arr2, double *sum_arr, int length) {
    // updates sum_arr to be arr1 + arr2
    for (int i = 0; i < length; i++) {
        sum_arr[i] = arr1[i] + arr2[i];
    }
}

// Subtract two arrays
static void sub_arr(double *arr1, double *arr2, double *sum_arr, int length) {
    // updates sum_arr to be arr1 - arr2
    for (int i = 0; i < length; i++) {
        sum_arr[i] = arr1[i] - arr2[i];
    }
}

// Multiply array by scalar
static void mult_arr(double *arr, double *arr_out, double scalar, int length) {
    // updates arr_out to be scalar * arr
    for (int i = 0; i < length; i++) {
        arr_out[i] = arr[i] * scalar;
    }
}

// Make a matrix multiplication of two square matrices
static void matrix_mult(double *A, double *B, double *Out, int length) {
    // updates Out to be A * B in a matrix sense in arbitrary dimensions
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < length; j++) {
            Out[length * i + j] = 0;
            for (int dummy = 0; dummy < length; dummy++) {
                Out[length * i + j] = Out[length * i + j] + A[length * i + dummy] * B[length * dummy + j];
            }
        }
    }
}

// Matrix multiplication with complex entries
static void complex_mult(std::complex<double> *matrix1, std::complex<double> *matrix2, std::complex<double> *matrixOUT, int length) {
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < length; j++) {
            matrixOUT[length * i + j] = 0;
            for (int dummy = 0; dummy < 3; dummy++) {
                matrixOUT[length * i + j] = matrixOUT[length * i + j] + matrix1[length * i + dummy] * matrix2[length * dummy + j];
            }
        }
    }
}
// Print a square matrix A of length length
static void print_m(double *A, int length) {
    // prints the matrix A in a nice form
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < length; j++) {
            std::cout << A[length * i + j] << "\t" ;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Print an entire event
static void print_event(double *P, int length) {
    // prints a hole event
    for (int i = 0; i < (length); i++) {
        std::cout << "P" << "[" << i + 1 <<"] = \t";
        for (int component = 0; component < 4; component++) {
            std::cout << (P[4 * i + component]) << ", ";
        }
        std::cout << "\t" << P[4 * i + 0] * P[4 * i + 0] - P[4 * i + 1] * P[4 * i + 1] - P[4 * i + 2] * P[4 * i + 2] - P[4 * i + 3] * P[4 * i + 3]  << std::endl;
    }
    std::cout << std::endl;
}

static double absolut(double value) {
    return sqrt(pow(value,2));
}

static void reset_arr(double *arr, int length) {
    for (int i = 0; i < length; i++) {
        arr[i] = 0.;
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for 4-vectors
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static void rcl_momenta(double *pp, double pp_rcl[][4], int Number_of_particles) {
    for (int i = 0; i < Number_of_particles; i++) {
        for (int j = 0; j < 4; j++) {
            pp_rcl[i][j] = pp[4 * i + j];
        }
    }
}


// print a 4-vector
static void print_p(double *p) {
    for (int i = 0; i < 4; i++) {
        std::cout << "P" << i << "\t" << p[i] << std::endl;
    }
}

// Returns the Minkovski scalar product
template <class R, class C>
static C minkovski(R *v1, C *v2) {
    // The function accepts any type as long as it is a pointer
    // For two arrays v1 = arr1, v2 = arr2
    // For vectors use vi = &vi[0]
    // Notice that the return value must be of the same type as the second argument
    // returns the minkovski scalar product of two Lorentz vectors
    C ans = v1[0] * v2[0] - v1[1] * v2[1] - v1[2] * v2[2]  - v1[3] * v2[3];
    return ans;
}

static std::complex<double> minkovski_arr_vec(double *a, std::vector<std::complex<double> > v) {
    std::complex<double> ans = a[0] * v[0] - a[1] * v[1] - a[2] * v[2] - a[3] * v[3];
    return ans;
}

// Multiplies a 4x4 matrix by a 4-vector
static void matrix_mult4(double *Lambda, double *v, double *v_out) {
    // updates v_out to be Lambda * v in a matrix sense.
    // The vectors must have dim 4
    for (int i = 0; i < 4; i++) {
        v_out[i] = 0;
        for (int j = 0; j < 4; j++) {
            v_out[i] = v_out[i] + Lambda[4 * i + j] * v[j];
        }
    }
    return ;
}

// Lam becomes rotation matrix
static void rotation(double x, double y, double z, double* Lam) {
  Lam[0] = 1; Lam[1] = 0; Lam[2] = 0; Lam[3] = 0;
  Lam[4] = 0; Lam[5] = std::cos(z)*std::cos(y);
  Lam[6] = std::cos(z)*std::sin(y)*std::sin(x)-std::sin(z)*std::cos(x);
  Lam[7] = std::cos(z)*std::sin(y)*std::cos(x)+std::sin(z)*std::sin(x);
  Lam[8] = 0; Lam[9] = std::sin(z)*std::cos(y);
  Lam[10] = std::sin(z)*std::sin(y)*std::sin(x)+std::cos(z)*std::cos(x);
  Lam[11] = std::sin(z)*std::sin(y)*std::cos(x)-std::cos(z)*std::sin(x);
  Lam[12] = 0; Lam[13] = -std::sin(y);
  Lam[14] = std::cos(y)*std::sin(x); Lam[15] = std::cos(y)*std::cos(x);
}

static void print_vec(std::vector<std::complex<double> > v) {
    for (int i = 0; i < v.size(); i++) {
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}

static void print_matrix(std::vector<std::vector<std::complex<double> > > M, int length) {
    for (int i = 0; i < length; i++) {
        for (int j = 0; j< length; j++) {
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// Utilities
static int factorial(int a) {
    int output = 1;
    for(int i = 1; i <= a; i++) output*=i;
    return output;
}

#endif