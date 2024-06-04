# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include  <math.h>
# include  "fft.h"
// #include <complex.h>
#include <gsl/gsl_sf_gamma.h>
# include <fftw3.h>
#include "omp.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

static double complex lngamma(double complex z)
{
  gsl_sf_result lnr, phi;
  gsl_sf_lngamma_complex_e(creal(z), cimag(z), &lnr, &phi);
  return lnr.val + I*phi.val;
}


static double complex polar (double r, double phi) {
  return (r*cos(phi) +I*(r*sin(phi)));
}

// static double complex lngamma(double complex z) {
//     return clog(gamma(z));
// }


static void lngamma_4(double x, double y, double* lnr, double* arg) {
    double complex w = lngamma(x+y*I);
    if(lnr) *lnr = creal(w);
    if(arg) *arg = cimag(w);
}

static double goodkr(int N, double mu, double q, double L, double kr) {
    double xp = (mu+1+q)/2;
    double xm = (mu+1-q)/2;
    double y = M_PI*N/(2*L);
    double lnr, argm, argp;
    lngamma_4(xp, y, &lnr, &argp);
    lngamma_4(xm, y, &lnr, &argm);
    double arg = log(2/kr) * N/L + (argp + argm)/M_PI;
    double iarg = round(arg);
    if(arg != iarg)
        kr *= exp((arg - iarg)*L/N);
    return kr;
}

void compute_u_coefficients(int N, double mu, double q, double L, double kcrc, double complex u[]) {
    double y = M_PI/L;
    double k0r0 = kcrc * exp(-L);
    double t = -2*y*log(k0r0/2);
    int m;
    if(q == 0) {
        double x = (mu+1)/2;
        double lnr, phi;
        int m;
        for(m = 0; m <= N/2; m++) {
            lngamma_4(x, m*y, &lnr, &phi);
            u[m] = polar(1.0,m*t + 2*phi);
        }
    }
    else {
        double xp = (mu+1+q)/2;
        double xm = (mu+1-q)/2;
        double lnrp, phip, lnrm, phim;
        for(m = 0; m <= N/2; m++) {
            lngamma_4(xp, m*y, &lnrp, &phip);
            lngamma_4(xm, m*y, &lnrm, &phim);
            u[m] = polar(exp(q*log(2) + lnrp - lnrm), m*t + phip - phim);
        }
    }

    for(m = N/2+1; m < N; m++)
        u[m] = conj(u[N-m]);
    if((N % 2) == 0)
      u[N/2] = (creal(u[N/2]) + I*0.0);
}

void fht(int N, const double r[], const double complex a[], double k[], double complex b[], double mu,
         double q, double kcrc, int noring, double complex* u, struct class_sz_structure * pclass_sz)
{
    double L = log(r[N-1]/r[0]) * N/(N-1.);
    double complex* ulocal = NULL;
    if(u == NULL) {
        if(noring)
            kcrc = goodkr(N, mu, q, L, kcrc);
        ulocal = malloc (sizeof(complex double)*N);
        compute_u_coefficients(N, mu, q, L, kcrc, ulocal);
        u = ulocal;
    }

// printf("fftlog coefficients computed.\n");

    int id = omp_get_thread_num();
    // omp_set_lock(lock); //Only a single thread writes
    // printf("My Thread num in fht is: %d\n", id);
    /* Compute the convolution b = a*u using FFTs */
    // fftw_init_threads();
    // fftw_make_planner_thread_safe();
    // lock();
    // fftw_plan_with_nthreads(1);
    // fftw_make_planner_thread_safe();
    // omp_lock_t lock;
    // #pragma omp single
    // omp_set_lock(lock); //Only a single thread writes
    // fftw_plan forward_plan = fftw_plan_dft_1d(N, (fftw_complex*) a, (fftw_complex*) b,  -1, FFTW_ESTIMATE);
    // printf("My Thread num in fht 1a is: %d\n", id);
    // fftw_plan reverse_plan = fftw_plan_dft_1d(N, (fftw_complex*) b, (fftw_complex*) b, +1, FFTW_ESTIMATE);
    // printf("My Thread num in fht 1b is: %d\n", id);
    // fftw_execute(forward_plan);

    /* Compute the convolution b = a*u using FFTs */
    fftw_execute_dft(pclass_sz->forward_plan, (fftw_complex*) a, (fftw_complex*) b);
// printf("first plan executed.\n");
    // printf("My Thread num in fht 2 is: %d\n", id);
    int m;
    for(m = 0; m < N; m++)
      b[m] *= u[m] / (double)(N);       // divide by N since FFTW doesn't normalize the inverse FFT
    // fftw_execute(reverse_plan);
    fftw_execute_dft(pclass_sz->reverse_plan, (fftw_complex*) b, (fftw_complex*) b);
    // fftw_destroy_plan(forward_plan);
    // fftw_destroy_plan(reverse_plan);
    // omp_unset_lock(lock);
// unlock();
    /* Reverse b array */
    double complex tmp;
    int n;
    for(n = 0; n < N/2; n++) {
        tmp = b[n];
        b[n] = b[N-n-1];
        b[N-n-1] = tmp;
    }

    /* Compute k's corresponding to input r's */
    double k0r0 = kcrc * exp(-L);
    k[0] = k0r0/r[0];

    for(n = 1; n < N; n++)
        k[n] = k[0] * exp(n*L/N);

    free(ulocal);
}

void fftlog_ComputeXiLM(int l, int m, int N, const double k[], const double pk[], double r[], double xi[], struct class_sz_structure * pclass_sz) {
  double complex* a = malloc(sizeof(complex double)*N);
  double complex* b = malloc(sizeof(complex double)*N);
      // fftw_complex* a=NULL;
      // fftw_complex* b=NULL;
      // a = fftw_alloc_complex(N);
      // b = fftw_alloc_complex(N);
    int i;
    for(i = 0; i < N; i++)
        // a[i] = pow(k[i], m - 0.5) * pk[i];
        a[i] = pow(k[i], 1 ) * pk[i]; // m = 1 in our case
    fht(N, k, a, r, b, 0, 0, 1, 1, NULL, pclass_sz);
    for(i = 0; i < N; i++)
        xi[i] = creal(pow(2*M_PI*r[i], -1) * b[i]);
        // xi[i] = creal(pow(2*M_PI*r[i], -1.5) * b[i]);
    //
    free(a);
    free(b);

      // fftw_free(a);
      // fftw_free(b);
}

void fftlog_ComputeXiLMsloz(int l, int m, int N, const double k[], const double pk[], double r[], double xi[], struct class_sz_structure * pclass_sz) {
  double complex* a = malloc(sizeof(complex double)*N);
  double complex* b = malloc(sizeof(complex double)*N);
  int i;
    for(i = 0; i < N; i++)
        a[i] = pow(k[i], m - 0.5) * pk[i];

    fht(N, k, a, r, b, l + 0.5, 0, 1, 1, NULL,pclass_sz);
    for(i = 0; i < N; i++)
        xi[i] = creal(pow(2*M_PI*r[i], -1.5) * b[i]);

    free(a);
    free(b);
}



void fftlog_ComputeXiLM_cl2gamma(int l, int m, int N, const double k[], const double pk[], double r[], double xi[], struct class_sz_structure * pclass_sz) {
  double complex* a = malloc(sizeof(complex double)*N);
  double complex* b = malloc(sizeof(complex double)*N);
      // fftw_complex* a=NULL;
      // fftw_complex* b=NULL;
      // a = fftw_alloc_complex(N);
      // b = fftw_alloc_complex(N);
    int i;
    for(i = 0; i < N; i++)
        // a[i] = pow(k[i], m - 0.5) * pk[i];
        a[i] = pow(k[i], 1 ) * pk[i]; //dim = 2 in our case
    fht(N, k, a, r, b, 2, 0, 1, 1, NULL, pclass_sz);
    for(i = 0; i < N; i++)
        xi[i] = creal(pow(2*M_PI*r[i], -1) * b[i]);
        // xi[i] = creal(pow(2*M_PI*r[i], -1.5) * b[i]);
    //
    free(a);
    free(b);

      // fftw_free(a);
      // fftw_free(b);
}


void pk2xi(int N, const double k[], const double pk[], double r[], double xi[], struct class_sz_structure * pclass_sz) {
    // fftlog_ComputeXiLM(0, 2, N, k, pk, r, xi, pclass_sz);
    fftlog_ComputeXiLM(0, 1, N, k, pk, r, xi, pclass_sz);
}

void xi2pk(int N, const double r[], const double xi[], double k[], double pk[], struct class_sz_structure * pclass_sz) {
    // static const double TwoPiCubed = 8*M_PI*M_PI*M_PI;
    double TwoPiCubed = pow(2.*M_PI,2.);///2./pow(2.*M_PI,-3)/4./pow(M_PI,3);
    // fftlog_ComputeXiLM(0, 2, N, r, xi, k, pk, pclass_sz);
    fftlog_ComputeXiLM(0, 1, N, r, xi, k, pk, pclass_sz);
    int j;
    for(j = 0; j < N; j++)
        pk[j] *= TwoPiCubed;
}

void cl2gamma(int N, const double k[], const double pk[], double r[], double xi[], struct class_sz_structure * pclass_sz) {
    // fftlog_ComputeXiLM(0, 2, N, k, pk, r, xi, pclass_sz);
    fftlog_ComputeXiLM_cl2gamma(0, 1, N, k, pk, r, xi, pclass_sz);
}

void gamma2cl(int N, const double r[], const double xi[], double k[], double pk[], struct class_sz_structure * pclass_sz) {
    // static const double TwoPiCubed = 8*M_PI*M_PI*M_PI;
    double TwoPiCubed = pow(2.*M_PI,2.);///2./pow(2.*M_PI,-3)/4./pow(M_PI,3);
    // fftlog_ComputeXiLM(0, 2, N, r, xi, k, pk, pclass_sz);
    fftlog_ComputeXiLM_cl2gamma(0, 1, N, r, xi, k, pk, pclass_sz);
    int j;
    for(j = 0; j < N; j++)
        pk[j] *= TwoPiCubed;
}
