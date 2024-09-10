// /**
//  * Function definitions of the FFT tool
//  * For more information, see fft.c
//  */
//
// #ifndef FFT_DEFINED
// #define FFT_DEFINED
// /**
//  * Compute FFT of two real inputs and store into two complex outputs of sizes N
//  * Assumes arrays are allocated and of size N
//  */
// void FFT_real(double* input_real_1, double* input_real_2, double* output_real_1, double* output_imag_1, double* output_real_2, double* output_imag_2, int N);
// /**
//  * Compute FFT of two real inputs and store into two complex outputs of sizes N
//  * Assumes arrays are allocated and of size N
//  * Only gives up nonredundant part for real FFT which have c_(N-n)=c_n
//  */
// void FFT_real_short(double* input_real_1, double* input_real_2, double* output_real_1, double* output_imag_1, double* output_real_2, double* output_imag_2, int N);
// /**
//  * Compute FFT of single complex input and stores into single complex output of size N
//  * Assumes arrays are allocated and of size N
//  */
// void FFT(double* input_real, double* input_imag, double* output_real, double* output_imag, int N, int stepsize);
// void DCT_real(double* input_real,double* input_imag,double* output_real,double* output_imag,int N);
// #endif
#include "common.h"
#include "class_sz.h"
# include <fftw3.h>

#ifndef FFT_H
#define FFT_H


/****************************************************************

This is the famous FFTLog.

First imlplemented by the living legend Andrew Hamilton:

http://casa.colorado.edu/~ajsh/FFTLog/

This version is a C version that was adapted from the C++ version found
in Copter JWG Carlson, another big loss for the cosmology community.

https://github.com/jwgcarlson/Copter

I've transformed this from C++ to C99 as the lowest common denominator
and provided bindings for C++ and python.

These are the C++ bindings

*****************************************************************/

// /* Compute the correlation function xi(r) from a power spectrum P(k), sampled
//  * at logarithmically spaced points k[j]. */
// void pk2xi(int N,  const double k[],  const double pk[], double r[], double xi[]);
//
// /* Compute the power spectrum P(k) from a correlation function xi(r), sampled
//  * at logarithmically spaced points r[i]. */
// void xi2pk(int N,  const double r[],  const double xi[], double k[], double pk[]);
//
// /* Compute the function
//  *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
//  * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
//  * in this notation.  The input k-values must be logarithmically spaced.  The
//  * resulting xi_l^m(r) will be evaluated at the dual r-values
//  *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
// void fftlog_ComputeXiLM(int l, int m, int N, const double k[],  const double pk[], double r[], double xi[]);
//
//
// #include <complex.h>
//
// /* Compute the discrete Hankel transform of the function a(r).  See the FFTLog
//  * documentation (or the Fortran routine of the same name in the FFTLog
//  * sources) for a description of exactly what this function computes.
//  * If u is NULL, the transform coefficients will be computed anew and discarded
//  * afterwards.  If you plan on performing many consecutive transforms, it is
//  * more efficient to pre-compute the u coefficients. */
// void fht(int N,  const double r[],  const double complex a[], double k[], double complex b[], double mu,
//          double q, double kcrc, int noring, double complex* u);
// //         double q = 0, double kcrc = 1, bool noring = true, double complex* u = NULL);
//
// /* Pre-compute the coefficients that appear in the FFTLog implementation of
//  * the discrete Hankel transform.  The parameters N, mu, and q here are the
//  * same as for the function fht().  The parameter L is defined (for whatever
//  * reason) to be N times the logarithmic spacing of the input array, i.e.
//  *   L = N * log(r[N-1]/r[0])/(N-1) */
// void compute_u_coefficients(int N, double mu, double q, double L, double kcrc, double complex u[]);
//


// int pkfft() ;
// void xi2pk(int N, const double r[], const double xi[], double k[], double pk[]);
// void pk2xi(int N, const double k[], const double pk[], double r[], double xi[]);
// void fftlog_ComputeXiLM(int l, int m, int N, const double k[], const double pk[], double r[], double xi[]);
// void fht(int N, const double r[], const double complex a[], double k[], double complex b[], double mu,
//          double q, double kcrc, int noring, double complex* u);
// void compute_u_coefficients(int N, double mu, double q, double L, double kcrc, double complex u[]);
// static double goodkr(int N, double mu, double q, double L, double kr);

int test01 ( );
int test02 ( );
int test03 ( );
int test04 ( );



/****************************************************************

This is the famous FFTLog.

First imlplemented by the living legend Andrew Hamilton:

http://casa.colorado.edu/~ajsh/FFTLog/

This version is a C version that was adapted from the C++ version found
in Copter JWG Carlson, another big loss for the cosmology community.

https://github.com/jwgcarlson/Copter

I've transformed this from C++ to C99 as the lowest common denominator
and provided bindings for C++ and python.

These are the C++ bindings

*****************************************************************/

/* Compute the correlation function xi(r) from a power spectrum P(k), sampled
 * at logarithmically spaced points k[j]. */
void pk2xi(int N,  const double k[],  const double pk[], double r[], double xi[], struct class_sz_structure * pclass_sz);
void cl2gamma(int N,  const double k[],  const double pk[], double r[], double xi[], struct class_sz_structure * pclass_sz);

/* Compute the power spectrum P(k) from a correlation function xi(r), sampled
 * at logarithmically spaced points r[i]. */
void xi2pk(int N,  const double r[],  const double xi[], double k[], double pk[], struct class_sz_structure * pclass_sz);
void gamma2cl(int N,  const double r[],  const double xi[], double k[], double pk[], struct class_sz_structure * pclass_sz);

/* Compute the function
 *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
 * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
 * in this notation.  The input k-values must be logarithmically spaced.  The
 * resulting xi_l^m(r) will be evaluated at the dual r-values
 *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
void fftlog_ComputeXiLM(int l, int m, int N, const double k[],  const double pk[], double r[], double xi[], struct class_sz_structure * pclass_sz);
void fftlog_ComputeXiLM_cl2gamma(int l, int m, int N, const double k[],  const double pk[], double r[], double xi[], struct class_sz_structure * pclass_sz);
void fftlog_ComputeXiLMsloz(int l, int m, int N, const double k[], const double pk[], double r[], double xi[],struct class_sz_structure * pclass_sz);

#include <complex.h>

/* Compute the discrete Hankel transform of the function a(r).  See the FFTLog
 * documentation (or the Fortran routine of the same name in the FFTLog
 * sources) for a description of exactly what this function computes.
 * If u is NULL, the transform coefficients will be computed anew and discarded
 * afterwards.  If you plan on performing many consecutive transforms, it is
 * more efficient to pre-compute the u coefficients. */
 void fht(int N, const double r[], const double complex a[], double k[], double complex b[], double mu,
          double q, double kcrc, int noring, double complex* u, struct class_sz_structure * pclass_sz);
//         double q = 0, double kcrc = 1, bool noring = true, double complex* u = NULL);

/* Pre-compute the coefficients that appear in the FFTLog implementation of
 * the discrete Hankel transform.  The parameters N, mu, and q here are the
 * same as for the function fht().  The parameter L is defined (for whatever
 * reason) to be N times the logarithmic spacing of the input array, i.e.
 *   L = N * log(r[N-1]/r[0])/(N-1) */
void compute_u_coefficients(int N, double mu, double q, double L, double kcrc, double complex u[]);



#endif // FFTLOG_H


/******************************************************************************/

// int all_tests ( )
