# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include  <math.h>
# include  "fft.h"
// #include <complex.h>
#include <gsl/gsl_sf_gamma.h>
# include <fftw3.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


// #include "fft.h"
// static const int N = 1000;
// #define PK_FNAME "/Users/boris/Work/CLASS-SZ/SO-SZ/FFTLog/test/pklogspaced.txt"
// #define OUTNAME "/Users/boris/Work/CLASS-SZ/SO-SZ/FFTLog/test/test_result.txt"

// // int main ( );
// void test01 ( );
// void test02 ( );
// void test03 ( );
// void test04 ( );
// void timestamp ( );
// int all_tests ( );


// void test01 ( );
// void test02 ( );
// void test03 ( );
// void test04 ( );
// // void timestamp ( );
// int all_tests ( );
//
// /******************************************************************************/
//
// int all_tests ( )
//
// /******************************************************************************/
// /*
//   Purpose:
//
//     MAIN is the main program for FFTW_TEST.
//
//   Discussion:
//
//     FFTW_TEST tests the FFTW library.
//
//   Licensing:
//
//     This code is distributed under the GNU LGPL license.
//
//   Modified:
//
//     05 November 2007
//
//   Author:
//
//     John Burkardt
// */
// {
//   // timestamp ( );
//   printf ( "\n" );
//   printf ( "FFTW_TEST\n" );
//   printf ( "  C version\n" );
//   printf ( "  Test the FFTW library.\n" );
//
//   test01 ( );
//   test02 ( );
//   test03 ( );
//   test04 ( );
// /*
//   Terminate.
// */
//   printf ( "\n" );
//   printf ( "FFTW_TEST\n" );
//   printf ( "  Normal end of execution.\n" );
//   printf ( "\n" );
//   // timestamp ( );
//
//   return 0;
// }
// /******************************************************************************/
//
// int test01 ( )
//
// /******************************************************************************/
// /*
//   Purpose:
//
//     TEST01: apply FFT to complex 1D data.
//
//   Discussion:
//
//     In this example, we generate N=100 random complex values stored as
//     a vector of type FFTW_COMPLEX named "IN".
//
//     We have FFTW compute the Fourier transform of this data named "OUT".
//
//     We have FFTW compute the inverse Fourier transform of "OUT" to get
//     "IN2", which should be the original input data, scaled by N.
//
//   Licensing:
//
//     This code is distributed under the GNU LGPL license.
//
//   Modified:
//
//     13 March 2017
//
//   Author:
//
//     John Burkardt
// */
// {
//   int i;
//   fftw_complex *in;
//   fftw_complex *in2;
//   int n = 100;
//   fftw_complex *out;
//   fftw_plan plan_backward;
//   fftw_plan plan_forward;
//   unsigned int seed = 123456789;
//
//   printf ( "\n" );
//   printf ( "TEST01\n" );
//   printf ( "  Demonstrate FFTW on a single vector of complex data.\n" );
//   printf ( "\n" );
//   printf ( "  Transform data to FFT coefficients.\n" );
//   printf ( "  Backtransform FFT coefficients to recover data.\n" );
//   printf ( "  Compare recovered data to original data.\n" );
// /*
//   Create the input array.
// */
//   in = fftw_malloc ( sizeof ( fftw_complex ) * n );
//
//   srand ( seed );
//
//   for ( i = 0; i < n; i++ )
//   {
//     in[i][0] = ( double ) rand ( ) / ( double ) RAND_MAX;
//     in[i][1] = ( double ) rand ( ) / ( double ) RAND_MAX;
//   }
//
//   printf ( "\n" );
//   printf ( "  Input Data:\n" );
//   printf ( "\n" );
//
//   for ( i = 0; i < n; i++ )
//   {
//     if ( i < 10 || n - 10 <= i )
//     {
//       printf ( "  %3d  %12f  %12f\n", i, in[i][0], in[i][1] );
//     }
//     if ( i == 10 )
//     {
//       printf ( "  ...  ............  ............\n" );
//     }
//   }
// /*
//   Create the output array.
// */
//   out = fftw_malloc ( sizeof ( fftw_complex ) * n );
//
//   plan_forward = fftw_plan_dft_1d ( n, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
//
//   fftw_execute ( plan_forward );
//
//   printf ( "\n" );
//   printf ( "  Output FFT Coefficients:\n" );
//   printf ( "\n" );
//
//   for ( i = 0; i < n; i++ )
//   {
//     if ( i < 10 || n - 10 <= i )
//     {
//       printf ( "  %3d  %12f  %12f\n", i, out[i][0], out[i][1] );
//     }
//     if ( i == 10 )
//     {
//       printf ( "  ...  ............  ............\n" );
//     }
//   }
// /*
//   Recreate the input array.
// */
//   in2 = fftw_malloc ( sizeof ( fftw_complex ) * n );
//
//   plan_backward = fftw_plan_dft_1d ( n, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE );
//
//   fftw_execute ( plan_backward );
//
//   printf ( "\n" );
//   printf ( "  Recovered input data divided by N:\n" );
//   printf ( "\n" );
//
//   for ( i = 0; i < n; i++ )
//   {
//     if ( i < 10 || n - 10 <= i )
//     {
//       printf ( "  %3d  %12f  %12f\n", i, in2[i][0] / ( double ) n, in2[i][1] / ( double ) n );
//     }
//     if ( i == 10 )
//     {
//       printf ( "  ...  ............  ............\n" );
//     }
//   }
// /*
//   Free up the allocated memory.
// */
//   fftw_destroy_plan ( plan_forward );
//   fftw_destroy_plan ( plan_backward );
//
//   fftw_free ( in );
//   fftw_free ( in2 );
//   fftw_free ( out );
//
//   return 1;
// }
// /******************************************************************************/
//
// int test02 ( )
//
// /******************************************************************************/
// /*
//   Purpose:
//
//     TEST02: apply FFT to real 1D data.
//
//   Licensing:
//
//     This code is distributed under the GNU LGPL license.
//
//   Modified:
//
//     13 March 2017
//
//   Author:
//
//     John Burkardt
// */
// {
//   int i;
//   double *in;
//   double *in2;
//   int n = 100;
//   int nc;
//   fftw_complex *out;
//   fftw_plan plan_backward;
//   fftw_plan plan_forward;
//   unsigned int seed = 123456789;
//
//   printf ( "\n" );
//   printf ( "TEST02\n" );
//   printf ( "  Demonstrate FFTW on a single vector of real data.\n" );
//   printf ( "\n" );
//   printf ( "  Transform data to FFT coefficients.\n" );
//   printf ( "  Backtransform FFT coefficients to recover data.\n" );
//   printf ( "  Compare recovered data to original data.\n" );
// /*
//   Set up an array to hold the data, and assign the data.
// */
//   in = fftw_malloc ( sizeof ( double ) * n );
//
//   srand ( seed );
//
//   for ( i = 0; i < n; i++ )
//   {
//     in[i] = ( double ) rand ( ) / ( double ) RAND_MAX;
//   }
//
//   printf ( "\n" );
//   printf ( "  Input Data:\n" );
//   printf ( "\n" );
//
//   for ( i = 0; i < n; i++ )
//   {
//     if ( i < 10 || n - 10 <= i )
//     {
//       printf ( "  %3d  %12f\n", i, in[i] );
//     }
//     if ( i == 10 )
//     {
//       printf ( "  ...  ............\n" );
//     }
//   }
// /*
//   Set up an array to hold the transformed data,
//   get a "plan", and execute the plan to transform the IN data to
//   the OUT FFT coefficients.
// */
//   nc = ( n / 2 ) + 1;
//
//   out = fftw_malloc ( sizeof ( fftw_complex ) * nc );
//
//   plan_forward = fftw_plan_dft_r2c_1d ( n, in, out, FFTW_ESTIMATE );
//
//   fftw_execute ( plan_forward );
//
//   printf ( "\n" );
//   printf ( "  Output FFT Coefficients:\n" );
//   printf ( "\n" );
//
//   for ( i = 0; i < nc; i++ )
//   {
//     if ( i < 10 || nc - 10 <= i )
//     {
//       printf ( "  %3d  %12f  %12f\n", i, out[i][0], out[i][1] );
//     }
//     if ( i == 10 )
//     {
//       printf ( "  ...  ............  ............\n" );
//     }
//   }
// /*
//   Set up an arrray to hold the backtransformed data IN2,
//   get a "plan", and execute the plan to backtransform the OUT
//   FFT coefficients to IN2.
// */
//   in2 = fftw_malloc ( sizeof ( double ) * n );
//
//   plan_backward = fftw_plan_dft_c2r_1d ( n, out, in2, FFTW_ESTIMATE );
//
//   fftw_execute ( plan_backward );
//
//   printf ( "\n" );
//   printf ( "  Recovered input data divided by N:\n" );
//   printf ( "\n" );
//
//   for ( i = 0; i < n; i++ )
//   {
//     if ( i < 10 || n - 10 <= i )
//     {
//       printf ( "  %3d  %12f\n", i, in2[i] / ( double ) ( n ) );
//     }
//     if ( i == 10 )
//     {
//       printf ( "  ...  ............\n" );
//     }
//   }
// /*
//   Release the memory associated with the plans.
// */
//   fftw_destroy_plan ( plan_forward );
//   fftw_destroy_plan ( plan_backward );
//
//   fftw_free ( in );
//   fftw_free ( in2 );
//   fftw_free ( out );
//
//   return 1;
// }
// /******************************************************************************/
//
// int test03 ( )
//
// /******************************************************************************/
// /*
//   Purpose:
//
//     TEST03: apply FFT to complex 2D data.
//
//   Discussion:
//
//     In this example, we generate NX=8 by NY=10 random complex values
//     stored as an NX by NY array of type FFTW_COMPLEX named "IN".
//
//     We have FFTW compute the Fourier transform of this data named "OUT".
//
//     We have FFTW compute the inverse Fourier transform of "OUT" to get
//     "IN2", which should be the original input data, scaled by NX * NY.
//
//     For a 2D complex NX by NY array used by FFTW, we need to access elements
//     as follows:
//
//       a[i*ny+j][0] is the real      part of A(I,J).
//       a[i*ny+j][1] is the imaginary part of A(I,J)..
//
//   Licensing:
//
//     This code is distributed under the GNU LGPL license.
//
//   Modified:
//
//     13 March 2017
//
//   Author:
//
//     John Burkardt
// */
// {
//   int i;
//   fftw_complex *in;
//   fftw_complex *in2;
//   int j;
//   int nx = 8;
//   int ny = 10;
//   fftw_complex *out;
//   fftw_plan plan_backward;
//   fftw_plan plan_forward;
//   unsigned int seed = 123456789;
//
//   printf ( "\n" );
//   printf ( "TEST03\n" );
//   printf ( "  Demonstrate FFTW on a %d by %d array of complex data.\n",
//     nx, ny );
//   printf ( "\n" );
//   printf ( "  Transform data to FFT coefficients.\n" );
//   printf ( "  Backtransform FFT coefficients to recover data.\n" );
//   printf ( "  Compare recovered data to original data.\n" );
// /*
//   Create the input array.
// */
//   in = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );
//
//   srand ( seed );
//
//   for ( i = 0; i < nx; i++ )
//   {
//     for ( j = 0; j < ny; j++ )
//     {
//       in[i*ny+j][0] = ( double ) rand ( ) / ( double ) RAND_MAX;
//       in[i*ny+j][1] = ( double ) rand ( ) / ( double ) RAND_MAX;
//     }
//   }
//
//   printf ( "\n" );
//   printf ( "  Input Data:\n" );
//   printf ( "\n" );
//
//   for ( i = 0; i < nx; i++ )
//   {
//     for ( j = 0; j < ny; j++ )
//     {
//       if ( ( i < 3 && j < 3 ) || ( nx - 3 <= i && ny - 3 <= j ) )
//       {
//         printf ( "  %3d  %3d  %12f  %12f\n", i, j, in[i*ny+j][0], in[i*ny+j][1] );
//       }
//       if ( i == 3 && j == 3 )
//       {
//         printf ( "  ...  ...  ............  ............\n" );
//       }
//     }
//   }
// /*
//   Create the output array.
// */
//   out = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );
//
//   plan_forward = fftw_plan_dft_2d ( nx, ny, in, out, FFTW_FORWARD,
//     FFTW_ESTIMATE );
//
//   fftw_execute ( plan_forward );
//
//   printf ( "\n" );
//   printf ( "  Output FFT Coefficients:\n" );
//   printf ( "\n" );
//
//   for ( i = 0; i < nx; i++ )
//   {
//     for ( j = 0; j < ny; j++ )
//     {
//       if ( ( i < 3 && j < 3 ) || ( nx - 3 <= i && ny - 3 <= j ) )
//       {
//         printf ( "  %3d  %3d  %12f  %12f\n", i, j, out[i*ny+j][0], out[i*ny+j][1] );
//       }
//       if ( i == 3 && j == 3 )
//       {
//         printf ( "  ...  ...  ............  ............\n" );
//       }
//     }
//   }
// /*
//   Recreate the input array.
// */
//   in2 = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );
//
//   plan_backward = fftw_plan_dft_2d ( nx, ny, out, in2, FFTW_BACKWARD,
//     FFTW_ESTIMATE );
//
//   fftw_execute ( plan_backward );
//
//   printf ( "\n" );
//   printf ( "  Recovered input data divided by NX * NY:\n" );
//   printf ( "\n" );
//
//   for ( i = 0; i < nx; i++ )
//   {
//     for ( j = 0; j < ny; j++ )
//     {
//       if ( ( i < 3 && j < 3 ) || ( nx - 3 <= i && ny - 3 <= j ) )
//       {
//         printf ( "  %3d  %3d  %12f  %12f\n", i, j,
//           in2[i*ny+j][0] / ( double ) ( nx * ny ),
//           in2[i*ny+j][1] / ( double ) ( nx * ny ) );
//       }
//       if ( i == 3 && j == 3 )
//       {
//         printf ( "  ...  ...  ............  ............\n" );
//       }
//     }
//   }
// /*
//   Free up the allocated memory.
// */
//   fftw_destroy_plan ( plan_forward );
//   fftw_destroy_plan ( plan_backward );
//
//   fftw_free ( in );
//   fftw_free ( in2 );
//   fftw_free ( out );
//
//   return 1;
// }
// /******************************************************************************/
//
// int test04 ( )
//
// /******************************************************************************/
// /*
//   Purpose:
//
//     TEST04: apply FFT to real 2D data.
//
//   Discussion:
//
//     In this example, we generate NX=8 by NY=10 random real values
//     stored as an NX by NY array of type DOUBLE named "IN".
//
//     We have FFTW compute the Fourier transform of this data named "OUT".
//
//     We have FFTW compute the inverse Fourier transform of "OUT" to get
//     "IN2", which should be the original input data, scaled by NX * NY.
//
//     The Fourier coefficients are stored in an NX by NYH array where
//     NYH = (NY/2) + 1.  We only compute about half the data because
//     of real data implies symmetric FFT coefficients.
//
//       a[i*nyh+j][0] is the real      part of A(I,J).
//       a[i*nyh+j][1] is the imaginary part of A(I,J)..
//
//   Licensing:
//
//     This code is distributed under the GNU LGPL license.
//
//   Modified:
//
//     13 March 2017
//
//   Author:
//
//     John Burkardt
// */
// {
//   int i;
//   double *in;
//   double *in2;
//   int j;
//   int nx = 100;
//   int ny = 100;
//   int nyh;
//   fftw_complex *out;
//   fftw_plan plan_backward;
//   fftw_plan plan_forward;
//   unsigned int seed = 123456789;
//
//   printf ( "\n" );
//   printf ( "TEST04\n" );
//   printf ( "  Demonstrate FFTW on a %d by %d array of real data.\n",
//     nx, ny );
//   printf ( "\n" );
//   printf ( "  Transform data to FFT coefficients.\n" );
//   printf ( "  Backtransform FFT coefficients to recover data.\n" );
//   printf ( "  Compare recovered data to original data.\n" );
// /*
//   Create the input array, an NX by NY array of doubles.
// */
//   in = ( double * ) malloc ( sizeof ( double ) * nx * ny );
//
//   srand ( seed );
//
//   for ( i = 0; i < nx; i++ )
//   {
//     for ( j = 0; j < ny; j++ )
//     {
//       in[i*ny+j] = ( double ) rand ( ) / ( double ) RAND_MAX;
//     }
//   }
//
//   printf ( "\n" );
//   printf ( "  Input Data:\n" );
//   printf ( "\n" );
//
//   for ( i = 0; i < nx; i++ )
//   {
//     for ( j = 0; j < ny; j++ )
//     {
//       if ( ( i < 3 && j < 3 ) || ( nx - 3 <= i && ny - 3 <= j ) )
//       {
//         printf ( "  %3d  %3d  %12f\n", i, j, in[i*ny+j] );
//       }
//       if ( i == 3 && j == 3 )
//       {
//         printf ( "  ...  ...  ............\n" );
//       }
//     }
//   }
// /*
//   Create the output array OUT, which is of type FFTW_COMPLEX,
//   and of a size NX * NYH that is roughly half the dimension of the input data
//   (ignoring the fact that the input data is real, and the FFT
//   coefficients are complex).
// */
//   nyh = ( ny / 2 ) + 1;
//
//   out = fftw_malloc ( sizeof ( fftw_complex ) * nx * nyh );
//
//   plan_forward = fftw_plan_dft_r2c_2d ( nx, ny, in, out, FFTW_ESTIMATE );
//
//   fftw_execute ( plan_forward );
//
//   printf ( "\n" );
//   printf ( "  Output FFT Coefficients:\n" );
//   printf ( "\n" );
//
//   for ( i = 0; i < nx; i++ )
//   {
//     for ( j = 0; j < nyh; j++ )
//     {
//
//       if ( ( i < 3 && j < 3 ) || ( nx - 3 <= i && nyh - 3 <= j ) )
//       {
//         printf ( "  %3d  %3d  %12f  %12f\n", i, j, out[i*nyh+j][0], out[i*nyh+j][1] );
//       }
//       if ( i == 3 && j == 3 )
//       {
//         printf ( "  ...  ...  ............  ............\n" );
//       }
//     }
//   }
// /*
//   Recreate the input array.
// */
//   in2 = ( double * ) malloc ( sizeof ( double ) * nx * ny );
//
//   plan_backward = fftw_plan_dft_c2r_2d ( nx, ny, out, in2, FFTW_ESTIMATE );
//
//   fftw_execute ( plan_backward );
//
//   printf ( "\n" );
//   printf ( "  Recovered input data divided by NX * NY:\n" );
//   printf ( "\n" );
//
//   for ( i = 0; i < nx; i++ )
//   {
//     for ( j = 0; j < ny; j++ )
//     {
//       if ( ( i < 3 && j < 3 ) || ( nx - 3 <= i && ny - 3 <= j ) )
//       {
//         printf ( "  %3d  %3d  %12f\n", i, j, in2[i*ny+j] / ( double ) ( nx * ny ) );
//       }
//       if ( i == 3 && j == 3 )
//       {
//         printf ( "  ...  ...  ............\n" );
//       }
//     }
//   }
// /*
//   Free up the allocated memory.
// */
//   fftw_destroy_plan ( plan_forward );
//   fftw_destroy_plan ( plan_backward );
//
//   free ( in );
//   free ( in2 );
//   fftw_free ( out );
//
//   return 1;
// }



// /* Computes the Gamma function using the Lanczos approximation */
// static double complex gamma(double complex z) {
//     /* Lanczos coefficients for g = 7 */
//     static double p[] = {
//         0.99999999999980993227684700473478,
//         676.520368121885098567009190444019,
//        -1259.13921672240287047156078755283,
//         771.3234287776530788486528258894,
//        -176.61502916214059906584551354,
//         12.507343278686904814458936853,
//        -0.13857109526572011689554707,
//         9.984369578019570859563e-6,
//         1.50563273514931155834e-7
//     };
//
//     if(creal(z) < 0.5)
//         return M_PI / (sin(M_PI*z)*gamma(1. - z));
//     z -= 1;
//     double complex x = p[0];
//     for(int n = 1; n < 9; n++)
//       x += p[n] / (z + (double)(n));
//     double complex t = z + 7.5;
//     return sqrt(2*M_PI) * cpow(t, z+0.5) * cexp(-t) * x;
// }

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
         double q, double kcrc, int noring, double complex* u, struct tszspectrum * ptsz)
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
    fftw_execute_dft(ptsz->forward_plan, (fftw_complex*) a, (fftw_complex*) b);

    // printf("My Thread num in fht 2 is: %d\n", id);
    int m;
    for(m = 0; m < N; m++)
      b[m] *= u[m] / (double)(N);       // divide by N since FFTW doesn't normalize the inverse FFT
    // fftw_execute(reverse_plan);
    fftw_execute_dft(ptsz->reverse_plan, (fftw_complex*) b, (fftw_complex*) b);
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

void fftlog_ComputeXiLM(int l, int m, int N, const double k[], const double pk[], double r[], double xi[], struct tszspectrum * ptsz) {
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
    fht(N, k, a, r, b, 0, 0, 1, 1, NULL, ptsz);
    for(i = 0; i < N; i++)
        xi[i] = creal(pow(2*M_PI*r[i], -1) * b[i]);
        // xi[i] = creal(pow(2*M_PI*r[i], -1.5) * b[i]);
    //
    free(a);
    free(b);

      // fftw_free(a);
      // fftw_free(b);
}

void fftlog_ComputeXiLM_cl2gamma(int l, int m, int N, const double k[], const double pk[], double r[], double xi[], struct tszspectrum * ptsz) {
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
    fht(N, k, a, r, b, 2, 0, 1, 1, NULL, ptsz);
    for(i = 0; i < N; i++)
        xi[i] = creal(pow(2*M_PI*r[i], -1) * b[i]);
        // xi[i] = creal(pow(2*M_PI*r[i], -1.5) * b[i]);
    //
    free(a);
    free(b);

      // fftw_free(a);
      // fftw_free(b);
}


void pk2xi(int N, const double k[], const double pk[], double r[], double xi[], struct tszspectrum * ptsz) {
    // fftlog_ComputeXiLM(0, 2, N, k, pk, r, xi, ptsz);
    fftlog_ComputeXiLM(0, 1, N, k, pk, r, xi, ptsz);
}

void xi2pk(int N, const double r[], const double xi[], double k[], double pk[], struct tszspectrum * ptsz) {
    // static const double TwoPiCubed = 8*M_PI*M_PI*M_PI;
    double TwoPiCubed = pow(2.*M_PI,2.);///2./pow(2.*M_PI,-3)/4./pow(M_PI,3);
    // fftlog_ComputeXiLM(0, 2, N, r, xi, k, pk, ptsz);
    fftlog_ComputeXiLM(0, 1, N, r, xi, k, pk, ptsz);
    int j;
    for(j = 0; j < N; j++)
        pk[j] *= TwoPiCubed;
}

void cl2gamma(int N, const double k[], const double pk[], double r[], double xi[], struct tszspectrum * ptsz) {
    // fftlog_ComputeXiLM(0, 2, N, k, pk, r, xi, ptsz);
    fftlog_ComputeXiLM_cl2gamma(0, 1, N, k, pk, r, xi, ptsz);
}

void gamma2cl(int N, const double r[], const double xi[], double k[], double pk[], struct tszspectrum * ptsz) {
    // static const double TwoPiCubed = 8*M_PI*M_PI*M_PI;
    double TwoPiCubed = pow(2.*M_PI,2.);///2./pow(2.*M_PI,-3)/4./pow(M_PI,3);
    // fftlog_ComputeXiLM(0, 2, N, r, xi, k, pk, ptsz);
    fftlog_ComputeXiLM_cl2gamma(0, 1, N, r, xi, k, pk, ptsz);
    int j;
    for(j = 0; j < N; j++)
        pk[j] *= TwoPiCubed;
}
