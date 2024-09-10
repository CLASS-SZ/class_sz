//==================================================================================================
// Author Jens Chluba Jan 2010
// Last modification: April 2011
//==================================================================================================
#ifndef PATTERSON_H
#define PATTERSON_H


#include "common.h"
#include "lensing.h"
//==================================================================================================
// patterson formulae & integration
// Integral_-1^1 f(x) dx = f(0)+Sum _1^order w[i]*(f(0-x[i])+f(0+x[i]))
//==================================================================================================
// int Integrate_using_Patterson(double a, double b, double epsrel, double epsabs,
//                               double (*fptr)(double), int *neval, double *r);


int Integrate_Patterson_refine(int reclev, double a, double b, double epsrel, double epsabs,
                               double (*fptr)(double, void *p), double *r, void *p, double * Patterson_fvals,  int Patterson_fvals_size);


int compute_integral_function_Patterson(int order, double xc, double Dx,
                                        double (*f)(double, void *p), int *neval, double *r,
                                        double * Patterson_fvals, int *Patterson_fvals_size, void *p);

int Integrate_using_Patterson(double a, double b, double epsrel, double epsabs,
                              double (*fptr)(double, void *p), int *neval, double *r, void *p, double * Patterson_fvals, int Patterson_fvals_size);

// double Integrate_using_Patterson_adaptive(double a, double b, double epsrel, double epsabs,
//                                           double (*fptr)(double));

double Integrate_using_Patterson_adaptive(double a, double b, double epsrel, double epsabs,
                                          double (*fptr)(double, void *p), void *p, int show_neval);



double DC_sumprod(const double *yarr, const double *zarr, int M);
#endif

//==================================================================================================
//==================================================================================================
