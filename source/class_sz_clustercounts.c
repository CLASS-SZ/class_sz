/** @file szpowerspectrum.c Documented SZ module. 2017-2021
 *
 * Boris Bolliet with inputs from Florian Ruppin, Thejs Brinckmann, Eunseong Lee++
 * based on the original Planck code szcounts.f90 in cosmomc
 *
 *This module is dedicated to the computation of
 *the number counts from Halo Mass Functions (HMF)
 */

#include "class_sz_clustercounts.h"
#include "class_sz_tools.h"
#include "Patterson.h"
#include "r8lib.h"
# include  "fft.h"
# include <fftw3.h>

int szcount_init(struct background * pba,
                 struct nonlinear * pnl,
                 struct primordial * ppm,
                 struct tszspectrum * ptsz,
                 struct szcount * pcsz)
{
  // ptsz->sz_verbose = ptsz->sz_verbose;
  // ptsz->has_sz_counts = _FALSE_;
  pcsz->has_sz_counts = ptsz->has_sz_counts;
  if (ptsz->has_sz_counts == _FALSE_
   || ptsz->has_sz_rates == _TRUE_)
  {
    if (ptsz->sz_verbose > 0)
      printf("->No SZ cluster counts requested. SZ cluster counts module skipped.\n");
      return _SUCCESS_;
  }

  else
  {

    if (ptsz->sz_verbose > 0)
      printf("->Computing SZ cluster counts.\n");

      if (ptsz->sigmaM_ym == 0.){
        if (ptsz->sz_verbose>0){
          printf("--> No scatter in ym relation.\n");
        }}

   // // if ((ptsz->experiment == 0 && ptsz->has_completeness_for_ps_SZ == 1)
   // //  || (ptsz->experiment == 0 && ptsz->has_sz_counts  == 1))
   //    read_Planck_noise_map(ptsz);

    initialise_and_allocate_memory_cc(ptsz,pcsz);


    if (ptsz->has_sz_counts_fft){
    class_call(compute_counts_sz_fft(pba,
                                    pnl,
                                    ppm,
                                    ptsz,
                                    pcsz),
                   pcsz->error_message,
                   pcsz->error_message);
    }
    else{
    class_call(compute_count_sz(pba,
                                pnl,
                                ppm,
                                ptsz,
                                pcsz),
               pcsz->error_message,
               pcsz->error_message);
    }

  }

  return _SUCCESS_;

}


int szcounts_free(struct szcount * pcsz,struct tszspectrum * ptsz)
{
  if (pcsz->has_sz_counts == _TRUE_){
  free(pcsz->redshift);
  int index_z;
  for (index_z=0;index_z<pcsz->Nbins_z;index_z++){
    free(pcsz->dNdzdy_theoretical[index_z]);
  }
  free(pcsz->dNdzdy_theoretical);
  free(pcsz->z_center);
  free(pcsz->steps_m);
  free(pcsz->steps_z);
  free(pcsz->logy);
  }

  if (ptsz->has_sz_counts_fft == _TRUE_){
    free(ptsz->array_y_to_m_redshift);
    free(ptsz->array_y_to_m_y);
    free(ptsz->array_y_to_m_at_z_y);
  }

  return _SUCCESS_;
}




int compute_counts_sz_fft(struct background * pba,
                          struct nonlinear * pnl,
                          struct primordial * ppm,
                          struct tszspectrum * ptsz,
                          struct szcount * pcsz)
{
  //clock_t begin = clock();
  if (ptsz->sz_verbose > 0)
    printf("->SZ_counts starting computations using ffts.\n");


// test m - >y at z;
double z = 0.3;
double m = 5e14;
double y = get_y_at_m_and_z(m,z,ptsz,pba);
printf("z = %.5e m = %.5e y = %.5e\n",z,m,y);
tabulate_y_to_m(pba,pnl,ppm,ptsz);
double mrec;
mrec = get_y_to_m_at_z_and_y(z,y,ptsz);
printf("z = %.5e m = %.5e y = %.5e mrec = %.6e\n",z,m,y,mrec);
double der = get_dlnm_dlny(log(y),z,ptsz);
printf("z = %.5e m = %.5e y = %.5e mrec = %.6e der = %.5e\n",z,m,y,mrec,der);
double dNdlny = get_dNdlny_at_z_and_y(z,y,pba,ptsz);
printf("z = %.5e m = %.5e y = %.5e mrec = %.6e der = %.5e dNdlny = %.5e\n",z,m,y,mrec,der,dNdlny);



    if(ptsz->sz_verbose>1) printf("constructing fftw plan\n");
    fftw_plan forward_plan, reverse_plan;
    // ptsz->N_samp_fftw = 100;
    fftw_complex* a_tmp;
    fftw_complex* b_tmp;
    a_tmp = fftw_alloc_complex(2*ptsz->N_samp_fftw);
    b_tmp = fftw_alloc_complex(2*ptsz->N_samp_fftw);
    forward_plan = fftw_plan_dft_1d(2*ptsz->N_samp_fftw,
                                    (fftw_complex*) a_tmp,
                                    (fftw_complex*) b_tmp,
                                    -1, FFTW_ESTIMATE);
    reverse_plan = fftw_plan_dft_1d(2*ptsz->N_samp_fftw,
                                    (fftw_complex*) b_tmp,
                                    (fftw_complex*) b_tmp,
                                    +1, FFTW_ESTIMATE);

    // forward_plan = fftw_plan_dft_1d(ptsz->N_samp_fftw,
    //                                 (fftw_complex*) a_tmp,
    //                                 (fftw_complex*) b_tmp,
    //                                 -1, FFTW_ESTIMATE);
    // reverse_plan = fftw_plan_dft_1d(ptsz->N_samp_fftw,
    //                                 (fftw_complex*) b_tmp,
    //                                 (fftw_complex*) b_tmp,
    //                                 +1, FFTW_ESTIMATE);

    fftw_free(a_tmp);
    fftw_free(b_tmp);

int index_parallel = 0;

int abort;

#ifdef _OPENMP
double tstart, tstop;
#endif

abort = _FALSE_;
/* number of threads (always one if no openmp) */
int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
    //omp_set_num_threads(number_of_threads);
  }
#endif

int id;
omp_lock_t lock;


#pragma omp parallel \
   shared(abort,pba,ptsz,ppm,pnl,lock,forward_plan,reverse_plan)\
   private(id,index_parallel)\
   num_threads(number_of_threads)
	 {

#ifdef _OPENMP
	   tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
for (index_parallel=0;index_parallel<1;index_parallel++)
	     {
#pragma omp flush(abort)
// double *in,*out;
// fftw_plan myplan;
int N = ptsz->N_samp_fftw;
// double z_fft = 0.3;
//

double sig2 = pow(0.5,2.);
double fac =1./sqrt(2.*_PI_*sig2);

double lny_fft[N];
double klny_fft[N];
double lnymin_fft =  ptsz->lnymin;//ptsz->lnymin;
double lnymax_fft =  ptsz->lnymax;//ptsz->lnymax;
// double dNdlny_fft[N];
double lognormal_fft[N];
// double test_func[N];
// double rout[N];
// double rin2[N];
// // double fac =1./sqrt(2.*_PI_*pow(ptsz->sigmaM_ym,2));
// double fac =1./sqrt(2.*_PI_);
// // in  = fftw_malloc (N*sizeof(double));
// // out = fftw_malloc (N*sizeof(double));
// //
// // myplan = fftw_plan_r2r_1d(N,in,out,FFTW_R2HC,FFTW_FORWARD);
// // exit(0);
  fftw_complex in[2*N], out[2*N], in2[2*N],test[2*N],out_test[2*N]; /* double [2] */
  fftw_complex product[2*N],product_out[2*N];
//   // double complex* in = malloc(sizeof(complex double)*N);
//   // double *in = malloc(sizeof(complex double)*N);
//   // double complex* out = malloc(sizeof(complex double)*N);
//
//   double complex* out_lognormal_fft = malloc(sizeof(complex double)*N);
//   double complex* out_dNdlny_fft = malloc(sizeof(complex double)*N);
//
//   double complex* in_lognormal_fft = malloc(sizeof(complex double)*N);
//   double complex* in_dNdlny_fft = malloc(sizeof(complex double)*N);
//
//   // double complex* in2 = malloc(sizeof(complex double)*N);
//
//   // fftw_plan p, q;
  int i;
//
//   // double * y
//
//   /* prepare a cosine wave */

printf("preparing arrays\n");
double dx = (lnymax_fft-lnymin_fft)/(N-1.);
  for (i = 0; i < N; i++) {
    double x = lnymin_fft+i*(lnymax_fft-lnymin_fft)
                 /(N-1.);
    lny_fft[i] = x;
//
//     dNdlny_fft[i] = get_dNdlny_at_z_and_y(z_fft,exp(lny_fft[i]),pba,ptsz);
//     in_dNdlny_fft[i] = dNdlny_fft[i];
//
//
//     // double arg0 = lny_fft[i]/(sqrt(2.)*ptsz->sigmaM_ym);
    double arg0 = x/sqrt(2.*sig2);///(sqrt(2.));
    lognormal_fft[i] = fac*exp(-arg0*arg0);
//     in_lognormal_fft[i] = lognormal_fft[i];
//
//     test_func[i] = pow(x,3)+pow(x,2);
//
    in[i][0]= lognormal_fft[i];//fac*exp(-arg0*arg0);//cos(3 * 2*M_PI*i/N);
    // in[i][0]= cos(3 * 2*M_PI*i/N);//fac*exp(-arg0*arg0);//cos(3 * 2*M_PI*i/N);
    in[i][1]= 0.;

    in[N+i][0]= 0.;
    in[N+i][1]= 0.;
//
    // test[i][0] = pow(x,3)+pow(x,2);
    double w = 10.;
    test[i][0] = sin(2 * M_PI * x/100.)*0.5*(1.-tanh((x-w/2.)))*0.5*(1.+tanh((x+w/2.)));
    test[i][1] = 0.;

    test[N+i][0] = 0.;
    test[N+i][1] = 0.;
//     // in[i] = dNdlny_fft[i];//fac*exp(-arg0*arg0);//cos(3 * 2*M_PI*i/N);
//
//     // printf("i = %d %.5e %.5e %.5e %.5e\n",i,lny_fft[i],arg0,lognormal_fft[i],in_lognormal_fft[i]);
//     printf("i = %d %.5e\n",i,in[i][0],in[i][1]);
//
//     // printf("N = %d y[i] = %.3e dN[i] = %.3e  arg0 = %.8e gauss = %.5e\n",
//     // N,
//     // // lnymin_fft,
//     // // lnymax_fft,
//     // lny_fft[i],
//     // in_dNdlny_fft[i],
//     // lny_fft[i],
//     // lny_fft[i]);
//     // in[i] =
//     // in[i][1] = 0;
  }
//
printf("array filled.\n");

  id = omp_get_thread_num();
//   /* forward Fourier transform, save the result in 'out' */
  // p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  // fftw_execute(p);
//   // fftw_execute_dft(forward_plan, (fftw_complex*) in_dNdlny_fft, (fftw_complex*) out_dNdlny_fft);
//   //

  printf("doing ffts\n");
  fftw_execute_dft(forward_plan, (fftw_complex*) in, (fftw_complex*) out);
  fftw_execute_dft(forward_plan, (fftw_complex*) test, (fftw_complex*) out_test);
printf("ffts done\n");

  // exit(0);
//   // fftw_execute_dft(reverse_plan, (fftw_complex*) out, (fftw_complex*) out);
//
// product[0][0]= out[][0]*out_test[i][0]
  for (i = 0; i < 2*N; i++){
 product[i][0]= (out[i][0]*out_test[i][0]-out[i][1]*out_test[i][1])/(2.*N);
 product[i][1]= (out[i][0]*out_test[i][1]+out[i][1]*out_test[i][0])/(2.*N);
  }
//
//
  fftw_execute_dft(reverse_plan, (fftw_complex*) product, (fftw_complex*) product_out);
//
//   // fftw_execute_dft(forward_plan, (fftw_complex*) in2, (fftw_complex*) in2);
//   //
//   // for (i = 0; i < N; i++){
//   //   // printf("thread  %d freq: %3d %+9.5f %+9.5f %+9.5f %+9.5f %+9.5f I\n",id, i, out[i], out[i],out_dNdlny_fft[i],out_dNdlny_fft[i],in_dNdlny_fft[i]);
//   //   // printf("thread  %d freq: %3d out example = %.6e  in = %.6e I\n",id, i, out[i], in[i]);
//   //   // printf("thread  %d freq: %3d out dN = %.6e  in = %.6e I\n",id, i, out_dNdlny_fft[i],in_dNdlny_fft[i]);
//   //   // printf("thread  %d freq: %3d out gaussian = %.6e\t in = %.6e I\n",id, i, out_lognormal_fft[i],in_lognormal_fft[i]);
//   //   printf("i = %d %.5e %.5e %.5e\n",i,lny_fft[i],lognormal_fft[i],in_lognormal_fft[i]);
//   //
//   // }
//     /* Reverse b array */
//     double complex tmp;
//     int n;
//     // for(n = 0; n < N/2; n++) {
//     //     tmp = out[n];
//     //     out[n] = out[N-n-1];
//     //     out[N-n-1] = tmp;
//     // }
    for(i = 0; i < N; i++){
        // rout[i] = creal(out[i])*1./N;
        // product_out[i][0] = product_out[i][0]/(double)(N);
        // product_out[i][1] = product_out[i][1]/(double)(N);

        printf("i = %d r = %.5e\n",i,product_out[i][0]);

        // rin2[i] = creal(in2[i])*1./N;;
      }

      // int n;
      // double complex tmp;
      // for(n = 0; n < N; n++) {
      //     tmp = product_out[n];
      //     product_out[n] = product_out[2*N-n-1];
      //     product_out[2*N-n-1] = tmp;
      // }

//
//   for (i = 0; i < N; i++){
//     // printf("thread  %d freq: %3d %+9.5f %+9.5f %+9.5f %+9.5f %+9.5f I\n",id, i, out[i], out[i],out_dNdlny_fft[i],out_dNdlny_fft[i],in_dNdlny_fft[i]);
//     printf("thread  %d freq: %3d in = %.6e  test = %.6e out = %.6e + %.6e I\n",id, i, lognormal_fft[i], test_func[i], product_out[i][0], product_out[i][1]);
//     // printf("thread  %d freq: %3d out dN = %.6e  in = %.6e I\n",id, i, out_dNdlny_fft[i],in_dNdlny_fft[i]);
//     // printf("thread  %d freq: %3d out gaussian = %.6e\t in = %.6e I\n",id, i, out_lognormal_fft[i],in_lognormal_fft[i]);
//       // printf("i = %d %.5e %.5e %.5e %.5e I\n",i,lny_fft[i],lognormal_fft[i],in_lognormal_fft[i],out_lognormal_fft[i]);
//       // printf("i = %3d %.5e %.5e %.5e %.5e I\n",i,lny_fft[i],lognormal_fft[i],in[i],out[i]);
//   // fprintf(fp,"%.5e\t%.5e\t%.5e\t%.5e\n",lny_fft[i],lognormal_fft[i],klny_fft[i],rout[i],rin2[i]);
//   }



        /* Compute k's corresponding to input r's */
    // double k0r0 = 1.;//kcrc * exp(-L);
    // klny_fft[0] = k0r0/lny_fft[0];

  //   for(n = 0; n < N; n++)
  //       klny_fft[n] = 1./klny_fft[n];// * exp(n*L/N);
  //
  //   for(i = 0; i < N; i++){
  //       rout[i] = creal(out[i])*1./N;
  //       rin2[i] = creal(in2[i])*1./N;;
  //     }
  //
  FILE *fp;
  char Filepath[_ARGUMENT_LENGTH_MAX_];
  sprintf(Filepath,"%s%s%s",ptsz->root,"test_ffts",".txt");
  fp=fopen(Filepath, "w");
  // fftw_execute_dft(forward_plan, (fftw_complex*) in_lognormal_fft, (fftw_complex*) out_lognormal_fft);
  for (i = 0; i < N; i++){
    // printf("thread  %d freq: %3d %+9.5f %+9.5f %+9.5f %+9.5f %+9.5f I\n",id, i, out[i], out[i],out_dNdlny_fft[i],out_dNdlny_fft[i],in_dNdlny_fft[i]);
    // printf("thread  %d freq: %3d in = %.6e  out = %.6e in2 = %.5e \n",id, i, lognormal_fft[i], rout[i],rin2[i]);
    // printf("thread  %d freq: %3d out dN = %.6e  in = %.6e I\n",id, i, out_dNdlny_fft[i],in_dNdlny_fft[i]);
    // printf("thread  %d freq: %3d out gaussian = %.6e\t in = %.6e I\n",id, i, out_lognormal_fft[i],in_lognormal_fft[i]);
      // printf("i = %d %.5e %.5e %.5e %.5e I\n",i,lny_fft[i],lognormal_fft[i],in_lognormal_fft[i],out_lognormal_fft[i]);
      // printf("i = %3d %.5e %.5e %.5e %.5e I\n",i,lny_fft[i],lognormal_fft[i],in[i],out[i]);
  // fprintf(fp,"%.5e\t%.5e\n",lny_fft[i],lognormal_fft[i],klny_fft[i],rout[i],rin2[i]);
  fprintf(fp,"%.5e\t%.5e\t%.5e\t%.5e\n",lny_fft[i],product_out[i][0],test[i][0],in[i][0]);
  }
  fclose(fp);
  //
  exit(0);

  // fftw_destroy_plan(p);

  // /* backward Fourier transform, save the result in 'in2' */
  // printf("\nInverse transform:\n");
 //  fftw_execute_dft(reverse_plan, (fftw_complex*) out, (fftw_complex*) in2);
 //  // q = fftw_plan_dft_1d(N, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE);
 //  // fftw_execute(q);
 //  // /* normalize */
 //  for (i = 0; i < N; i++) {
 //    in2[i] *= 1./N;
 //    // in2[i][1] *= 1./N;
 //  }
 //  for (i = 0; i < N; i++){
 //    printf("thread  %d recover: %d %.5e %.5e I vs. %.5e %.5e I\n",
 //        id,i, in[i], in[i], in2[i], in2[i]);
 //
 //   printf("%i %.3e %.3e\n",i,in[i],in[i]);
 // }
    //
    // printf("thread  %d recover: %3d %+9.5e %+9.5e I vs. %+9.5e %+9.5e I\n",
    //     id,i, in[i], creal(in[i]), in2[i], in2[i]);
  // fftw_destroy_plan(q);

  // fftw_cleanup();

// test chatgpt:


    n = ptsz->N_samp_fftw; // input array size
    double *in1r, *in2r, *outr, *xin;
    fftw_plan p1, p2;

    // allocate input and output arrays
    xin = (double*) fftw_malloc(sizeof(double) * n);
    in1r = (double*) fftw_malloc(sizeof(double) * n);
    in2r = (double*) fftw_malloc(sizeof(double) * n);
    outr = (double*) fftw_malloc(sizeof(double) * n);

    // initialize input arrays
    for (i = 0; i < n; i++) {

    double x = lnymin_fft+i*(lnymax_fft-lnymin_fft)
                 /(n-1.);

    // double arg0 = x;///(sqrt(2.));

       // lognormal_fft[i] = fac*exp(-arg0*arg0);
        xin[i] = x;
        double arg0 = x/sqrt(2.*sig2);///(sqrt(2.));
        in1r[i] = fac*exp(-arg0*arg0);

        // in1[i] = sin(2 * M_PI * i / n);
        // in2[i] = cos(4 * M_PI * i / n);
        // in2[i] = pow(x,3)+pow(x,2);
        double w = 10.;
        in2r[i] = sin(2 * M_PI * x/100.)*0.5*(1.-tanh((x-w/2.)))*0.5*(1.+tanh((x+w/2.)));;
    }

    // create FFTW plans for forward and inverse transforms
    p1 = fftw_plan_r2r_1d(n, in1r, in1r, FFTW_R2HC, FFTW_ESTIMATE);
    p2 = fftw_plan_r2r_1d(n, in2r, in2r, FFTW_R2HC, FFTW_ESTIMATE);

    // perform forward transforms
    fftw_execute(p1);
    fftw_execute(p2);

    // perform convolution in frequency domain
    outr[0] = in1r[0] * in2r[0];
    for (i = 1; i < n; i++) {
        outr[i] = in1r[i] * in2r[i];// - in1r[n - i] * in2r[n - i];
    }

    // create FFTW plan for inverse transform
    p1 = fftw_plan_r2r_1d(n, outr, outr, FFTW_HC2R, FFTW_ESTIMATE);

    // perform inverse transform
    fftw_execute(p1);

    // normalize output array
    for (i = 0; i < n; i++) {
        outr[i] /= n;
    }

    // print output array


    // FILE *fp;
    // char Filepath[_ARGUMENT_LENGTH_MAX_];
    sprintf(Filepath,"%s%s%s",ptsz->root,"test_ffts_real",".txt");
    fp=fopen(Filepath, "w");
    for (i = 0; i < n; i++) {
        printf("%.5e\t%.5e\n",xin[i],outr[i]);

        fprintf(fp,"%.5e\t%.5e\n",xin[i],outr[i]);

    }
    fclose(fp);

    // free memory and destroy plans
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(in1r);
    fftw_free(in2r);
    fftw_free(outr);

    // end test chatgpt

          }
#ifdef _OPENMP
      tstop = omp_get_wtime();
      if (ptsz->sz_verbose > 0)
         printf("In %s: time spent in parallel region (loop over cluster counts ffts) = %e s for thread %d\n",
                   __func__,tstop-tstart,omp_get_thread_num());


#endif
   // free(Pvecback);
   // free(Pvectsz);
   // free(b_l1_l2_l_1d);
	} //end of parallel region

   if (abort == _TRUE_) return _FAILURE_;

   fftw_destroy_plan(forward_plan);
   fftw_destroy_plan(reverse_plan);
  exit(0);

  return _SUCCESS_;
}

int compute_count_sz(struct background * pba,
                     struct nonlinear * pnl,
                     struct primordial * ppm,
                     struct tszspectrum * ptsz,
                     struct szcount * pcsz)
{
  //clock_t begin = clock();

// return _SUCCESS_;


  if (ptsz->sz_verbose > 0)
    printf("->SZ_counts starting grid computation.\n");
  ///PARALLEL
  double * pvecsz;
  int abort;
#ifdef _OPENMP
  //instrumentation times
  double tstart, tstop;
#endif
  abort = _FALSE_;
  // beginning of parallel region

int index_y;
#pragma omp parallel \
shared(abort,pba,ppm,pnl,ptsz,pcsz)\
private(tstart, tstop,pvecsz,index_y)
  {


#ifdef _OPENMP
    tstart = omp_get_wtime();
#endif
    class_alloc_parallel(pvecsz,
                         pcsz->pvecsz_size*sizeof(double),
                         pcsz->error_message);
    int i;
    for(i = 0; i<pcsz->pvecsz_size;i++) pvecsz[i] = 0.;

#pragma omp for schedule (dynamic)
    for (index_y=0; index_y<pcsz->Nbins_y; index_y ++){
#pragma omp flush(abort)

      pvecsz[pcsz->index_y] = index_y;
      class_call_parallel(grid_C_2d(pvecsz,
                                    pba,
                                    ppm,
                                    pnl,
                                    ptsz,
                                    pcsz),
                          pcsz->error_message,
                          pcsz->error_message);

    }

#ifdef _OPENMP
    tstop = omp_get_wtime();
    if (ptsz->sz_verbose > 0)
      printf("In %s: time spent in parallel region (loop over s/n's) = %e s for thread %d\n",
             __func__,tstop-tstart,omp_get_thread_num());
#endif
    free(pvecsz);
  }
  if (abort == _TRUE_) return _FAILURE_;


  //end bin in mass
  if (ptsz->sz_verbose > 0)
    printf("->SZ_counts computations done.\n");

  write_output_cluster_counts(pcsz,ptsz);




  return _SUCCESS_;
}




int grid_C_2d(
              double * pvecsz,
              struct background *pba,
              struct primordial * ppm,
              struct nonlinear * pnl,
              struct tszspectrum * ptsz,
              struct szcount * pcsz
              ){
  //int i;

  const int dim_1 = ptsz->Ny;
  const int dim_2 = ptsz->nthetas;




  int l_array[3];
  double theta_array[3];

  double * completeness_2d_to_1d = NULL;
  class_alloc(completeness_2d_to_1d,pcsz->nsteps_m*pcsz->nsteps_z*sizeof(double *),pcsz->error_message);
  int index_z, index_m, index_y;
double ** completeness_2d = NULL;
class_alloc(completeness_2d,pcsz->nsteps_m*sizeof(double *),pcsz->error_message);
for (index_m=0;index_m<pcsz->nsteps_m;index_m++)
{
  class_alloc(completeness_2d[index_m],pcsz->nsteps_z*sizeof(double *),pcsz->error_message);
}




  double fsky = 0.;
  int index_patches;
  for (index_patches=0;index_patches<ptsz->nskyfracs;index_patches++)
    fsky += ptsz->skyfracs[index_patches];

int index_m_z = 0;

  for (index_m=0;index_m<pcsz->nsteps_m;index_m++)
  {

    for (index_z=0;index_z<pcsz->nsteps_z;index_z++){

      completeness_2d_to_1d[index_m_z]=1e-300;
      completeness_2d[index_m][index_z] = 0.;
      index_m_z += 1;

    }
  }



  index_y = (int) pvecsz[pcsz->index_y];

  double y_min,y_max;

  y_min = pow(10., pcsz->logy[index_y] - pcsz->dlogy/2.);
  y_max = pow(10., pcsz->logy[index_y] + pcsz->dlogy/2.);

//   if (index_y != pcsz->Nbins_y){
//   y_min = pow(10., pcsz->logy[index_y] - pcsz->dlogy/2.);
//   y_max = pow(10., pcsz->logy[index_y] + pcsz->dlogy/2.);
//   }
//   else{
//
//   y_min = pow(10., pcsz->logy[index_y] - ptsz->bin_dlog10_snr_last_bin/2.);
//   y_max = 1e100;//pow(10., pcsz->logy[index_y] + ptsz->bin_dlog10_snr_last_bin/2.);
// }

  if (ptsz->sz_verbose > 3){
    printf("->SZ_counts grid_C_2d.\n");
    printf("->In signal-to-noise bin:\n");
    printf("->bin id = %d y_min = %.3e y_max = %.3e\n",index_y,y_min,y_max);
    }

if (pcsz->has_completeness == 1){
  if (ptsz->sigmaM_ym == 0.){
    if (ptsz->sz_verbose>0){
      printf("--> No scatter in ym relation.\n");
    }
    int index_m_z = 0;

// for (index_m=0;index_m<pcsz->nsteps_m;index_m++){

    for (index_z=0;index_z<pcsz->nsteps_z;index_z++){
      double zp = pcsz->steps_z[index_z];

      for (index_m=0;index_m<pcsz->nsteps_m;index_m++){

        //compute_theta_and_y_at_z_and_m
        double mp= exp(pcsz->steps_m[index_m]);

        double m500c = 0.;
        double m_ym = mp;



        if (ptsz->integrate_wrt_m200m == 1){
          m500c = get_m200m_to_m500c_at_z_and_M(zp,mp,ptsz);
        }
        if (ptsz->integrate_wrt_m200c == 1){
          m500c = get_m200c_to_m500c_at_z_and_M(zp,mp,ptsz);
        }
        //
        if (ptsz->integrate_wrt_m500c == 1){
          m500c = mp;
        }

        if (ptsz->use_m500c_in_ym_relation == 1){
        m_ym = m500c;
        }

        // printf("m_ym = %.5e\n",m_ym);
        // exit(0);
        // zp = 1.00000e-05;
        double yp = get_y_at_m_and_z(m_ym,zp,ptsz,pba);

        // if (index_y == 2  && index_z == 0 && index_m == 0){
        // printf("\n");
        // printf("\n");
        // printf("m_ym = %.5e zp = %.5e yp = %.5e\n",m_ym,zp,yp);
        // printf("\n");
        // printf("\n");
        // exit(0);
        // }


        double thp = get_theta_at_m_and_z(m500c,zp,ptsz,pba);
        //Planck

        // if not planck, apply the mismatch function with C correction
        if (ptsz->experiment == 1){
          double m_pivot = ptsz->m_pivot_ym*pba->h;// 1. convert to msun/h //not that it used to be *0.7 as in hasselfield paper.
          double m_over_m_pivot_500c = mp/m_pivot;
          thp = thp*pow(m_over_m_pivot_500c,ptsz->C_ym);
        }


        // find_theta_bin(ptsz,thp,l_array,theta_array);
        // int l1 = l_array[1];
        // int l2 = l_array[2];
        // double th1 = theta_array[1];
        // double th2 = theta_array[2];


        int index_patches;
        double comp_sky_tot = 0.;
        for (index_patches =0;index_patches<ptsz->nskyfracs;index_patches++){

          // double y1 = ptsz->ylims[index_patches][l1];
          // double y2 = ptsz->ylims[index_patches][l2];
          // double y = y1 + (y2-y1)/(th2-th1)*(thp-th1);

          double y_interp = pwl_value_1d(ptsz->nthetas,
                                         ptsz->thetas,
                                         ptsz->ylims[index_patches],
                                         thp);
          // double y_interp = pwl_value_1d(ptsz->nthetas,
          //                                ptsz->thetas,
          //                                ptsz->sky_averaged_ylims,
          //                                thp); // ~5% difference

          double y = y_interp;
          // printf("y = %.5e y_interp = %.5e r = %.5e\n",y,y_interp,y/y_interp);

          double c2;

          if (ptsz->use_planck_binned_proba == 1){
          c2 = erf_compl(yp,y,ptsz->sn_cutoff);
          c2 *= erf_compl(yp,y,y_min);
          c2 *= (1.-erf_compl(yp,y,y_max));

          // if (index_y == 2 && index_patches == 0 && index_z == 0 && index_m == 0){
          // printf("\n");
          // printf("mp = %.5e zp = %.5e yp = %.5e y = %.5e c2 = %.5e\n",mp,zp,yp,y,c2);
          // printf("\n");
          // exit(0);
          // }


          if (index_y == 0){
            c2 = erf_compl(yp,y,ptsz->sn_cutoff);
            c2 *= (1.-erf_compl(yp,y,y_max));


          }

          if (index_y == pcsz->Nbins_y-1){
            c2 = erf_compl(yp,y,ptsz->sn_cutoff) ;
            c2 *= erf_compl(yp,y,y_min);


          }}

          else {
            c2 = erf_compl_nicola(yp,y,ptsz->sn_cutoff,y_min,y_max);
          }

          completeness_2d[index_m][index_z] += c2*ptsz->skyfracs[index_patches];///fsky;
          // printf("%d\n",index_patches);

          // completeness_2d_to_1d[index_m_z] += c2*ptsz->skyfracs[index_patches];
          // completeness_2d_to_1d[index_m_z] += 1.*ptsz->skyfracs[index_patches];

        } // end loop patches
        // completeness_2d_to_1d[index_m_z] = log(completeness_2d_to_1d[index_m_z]);
        index_m_z += 1;

        if (completeness_2d[index_m][index_z]>1.) completeness_2d[index_m][index_z] = 1.;
        if (completeness_2d[index_m][index_z]<0.) completeness_2d[index_m][index_z] = 0.;

        // double comp_sky_tot = 0.;
        // for (index_patches =0;index_patches<ptsz->nskyfracs;index_patches++){
        // }


      }//end m loop
    }//end z loop
  }//end if sigmaM=0

  else {
    double fac =1./sqrt(2.*_PI_*pow(ptsz->sigmaM_ym,2));

    int index1,index2;

    double ** erfs = NULL;
    double * erfs_2d_to_1d = NULL;



    class_alloc(erfs_2d_to_1d,
                ptsz->Ny*ptsz->nthetas*sizeof(double *),
                ptsz->error_message);

    class_alloc(erfs,
                dim_1*sizeof(double *),
                pcsz->error_message);

    int index_y_th = 0;
    for (index1=0;index1<dim_1;index1++)
    {
      class_alloc(erfs[index1],dim_2*sizeof(double*),pcsz->error_message);


      for (index2=0;index2<dim_2;index2++){

          erfs[index1][index2]=0.;
          erfs_2d_to_1d[index_y_th] = 1.e-300;
          index_y_th += 1;


      }
    }



    if (ptsz->sz_verbose > 3)
      printf("->SZ_counts grid_C_2d debug 1.\n");
    //double fsky = 0.;
    //int index_patches;
    // int index2,index1;

    ////// tabulate erfs as a function of theta and y in each s/n bin
    int index_th_y = 0;
    for (index2=0;index2<ptsz->Nth;index2++){
      // double lny = pcsz->lnymin;
      //double th1 = exp(ptsz->erfs_2d_to_1d_y_array[index2])

      for (index1=0;index1<ptsz->Ny;index1++)
      {
        double y0=exp(ptsz->erfs_2d_to_1d_y_array[index1]);
        for (index_patches=0;index_patches<ptsz->nskyfracs;index_patches++){
          //fsky += ptsz->skyfracs[index_patches];

          double y1 = ptsz->ylims[index_patches][index2];
          // double y1 = get_ylim_of_theta(th1,ptsz->ylims[index_patches][index2];
          // int k = index_y;
          //
          // double qmin=pcsz->logy[k]-pcsz->dlogy/2.;
          // double qmax=pcsz->logy[k]+pcsz->dlogy/2.;
          // double q_min=log10(y_min);
          // double q_max=log10(y_max);

          double c2;

          if (ptsz->use_planck_binned_proba == 1){
          if (index_y==0)  {c2=erf_compl(y0,y1,ptsz->sn_cutoff)*(1.-erf_compl(y0,y1,y_max));}
          else if (index_y==pcsz->Nbins_y-1) {c2=erf_compl(y0,y1,y_min)*erf_compl(y0,y1,ptsz->sn_cutoff);}
          else {c2=erf_compl(y0,y1,ptsz->sn_cutoff)*erf_compl(y0,y1,y_min)*(1.-erf_compl(y0,y1,y_max));}
          }
          else{
            c2 = erf_compl_nicola(y0,y1,ptsz->sn_cutoff,y_min,y_max);
          }
          erfs[index1][index2]=erfs[index1][index2]+c2*ptsz->skyfracs[index_patches];

          // erfs_2d_to_1d[index_th_y] += c2*ptsz->skyfracs[index_patches]/fsky;


        } //end loop patches
        // erfs_2d_to_1d[index_th_y] = log(erfs_2d_to_1d[index_th_y]);
        index_th_y += 1;
      } //end loop y
    } //end loop thetas

    ////// end tabulate erfs as a function of theta and y


    // tabulate completeness as a function of z and m
    // integrate erfs wrt y at all (z,M) to get completeness in a z,m grid
  if (ptsz->sz_verbose > 3)
      printf("->SZ_counts grid_C_2d debug 2.\n");

    int index_m_z = 0;
    for (index_z=0;index_z<pcsz->nsteps_z;index_z++){

      double zp = pcsz->steps_z[index_z];
      // if (ptsz->sz_verbose > 3)
      //   printf("->SZ_counts grid_C_2d debug 3, z = %.4e.\n",zp);

      for (index_m=0;index_m<pcsz->nsteps_m;index_m++){

        double mp= exp(pcsz->steps_m[index_m]);
        if (ptsz->integrate_wrt_m200m == 1){
          mp = get_m200m_to_m500c_at_z_and_M(zp,mp,ptsz);
        }
        if (ptsz->integrate_wrt_m200c == 1){
          mp = get_m200c_to_m500c_at_z_and_M(zp,mp,ptsz);
        }

        double yp = get_y_at_m_and_z(mp,zp,ptsz,pba);
        double thp = get_theta_at_m_and_z(mp,zp,ptsz,pba);

        // if not planck, apply the mismatch function with C correction
        if (ptsz->experiment == 1){
          double m_pivot = ptsz->m_pivot_ym*pba->h;
          double m_over_m_pivot_500c = mp/m_pivot;
          thp = thp*pow(m_over_m_pivot_500c,ptsz->C_ym);
        }

        find_theta_bin(ptsz,thp,l_array,theta_array);
        int l1 = l_array[1];
        int l2 = l_array[2];
        double th1 = theta_array[1];
        double th2 = theta_array[2];

        double mu = log(yp);
        double int_comp =1.e-300;

        // double mu_high = mu + 5.*(sqrt(2.)*ptsz->sigmaM_ym);
        // int l1y_high, l2y_high;
        // int l1y_low, l2y_low;
        // find_y_bin(ptsz,mu_high,l_array,theta_array);
        // l1y_high = l_array[1];
        // l2y_high = l_array[2];
        // // double y1 = theta_array[1];
        // // double y2 = theta_array[2];
        //
        // double mu_low = mu - 5.*(sqrt(2.)*ptsz->sigmaM_ym);
        // find_y_bin(ptsz,mu_low,l_array,theta_array);
        // l1y_low = l_array[1];
        // l2y_low = l_array[2];
        // // double y1 = theta_array[1];
        // // double y2 = theta_array[2];

        //
        // if (l1y_low != l2y_low && l1y_high != l2y_high){



        // double y = yp;




        // double lny=pcsz->lnymin;
        int k;

        // if (ptsz->sz_verbose > 3)
        //   printf("->SZ_counts grid_C_2d debug 5, z = %.4e.\n",zp);

        // at a fixed theta(z,m)
        // integrate over y, erf(theta,y)*fac/y*exp(-arg(y))

        for (k=0;k<ptsz->Ny-1;k++){
        // for (k=l1y_low;k<l2y_high;k++){
          // printf("k = %d int_comp1 = %e\n",k,int_comp);
          // double y0=exp(lny);
          double y0 = exp(ptsz->erfs_2d_to_1d_y_array[k]);

          // printf("k = %d int_comp2 = %e\n",k,int_comp);
          // double y=exp(lny+pcsz->dlny);
          double y = exp(ptsz->erfs_2d_to_1d_y_array[k+1]);
          // y = exp(ptsz->erfs_2d_to_1d_y_array[k]);
          // printf("k = %d int_comp3 = %e\n",k,int_comp);
          double dy=y-y0;
          // printf("k = %d int_comp4 = %e\n",k,int_comp);
          double arg0=((ptsz->erfs_2d_to_1d_y_array[k]-mu)/(sqrt(2.)*ptsz->sigmaM_ym));

          double win0=erfs[k][l1]+(erfs[k][l2]-erfs[k][l1])/(th2-th1)*(thp-th1);
          double win=erfs[k+1][l1]+(erfs[k+1][l2]-erfs[k+1][l1])/(th2-th1)*(thp-th1);

          // double ekl1 = get_detection_proba_at_y_and_theta(y0,th1,erfs_2d_to_1d,ptsz);
          // double ekl2 = get_detection_proba_at_y_and_theta(y0,th2,erfs_2d_to_1d,ptsz);
          // double win0 = ekl1+(ekl2-ekl1)/(th2-th1)*(thp-th1);
          //
          // ekl1 = get_detection_proba_at_y_and_theta(y,th1,erfs_2d_to_1d,ptsz);
          // ekl2 = get_detection_proba_at_y_and_theta(y,th2,erfs_2d_to_1d,ptsz);
          //
          // double win = ekl1+(ekl2-ekl1)/(th2-th1)*(thp-th1);
          //

          //double arg=((lny+pcsz->dlny-mu)/(sqrt(2.)*ptsz->sigmaM_ym));
          double arg=((ptsz->erfs_2d_to_1d_y_array[k+1]-mu)/(sqrt(2.)*ptsz->sigmaM_ym));
          double py=(win0*fac/y0*exp(-arg0*arg0)+win*fac/y*exp(-arg*arg))*0.5;
          // if (fabs(py)<1e-100)
          //   py = 0.;
          //lny=lny+pcsz->dlny;

          int_comp=int_comp+py*dy;
          // if (py*dy<0)
          // printf("k = %d int_comp15 = %.5e py = %.5e dy = %.5e win0 = %.5e win = %.5e th1 = %.5e th2 = %.5e thp = %.5e\n",k,int_comp,py,dy,win0,win,th1,th2,thp);
          //
          }
        // }
//
// // printf("int_compe = %.3e\n",int_comp);
// struct Parameters_for_integrand_cluster_counts_completeness X;
//   X.ptsz = ptsz;
//   X.erfs_2d_to_1d = erfs_2d_to_1d;
//   X.theta = thp;
//   X.theta1 = th1;
//   X.theta2 = th2;
//   X.y = yp;
//
//   void * params = &X;
//
//
//   double epsrel = 1e-2;
//   double epsabs = 1e-80;
//
//   double lny_min = ptsz->erfs_2d_to_1d_y_array[0];
//   double lny_max = ptsz->erfs_2d_to_1d_y_array[ptsz->Ny-1];
//
//   int_comp=Integrate_using_Patterson_adaptive(lny_min,
//                                               lny_max,
//                                               epsrel, epsabs,
//                                               integrand_cluster_counts_completeness,
//                                               params,0);
// printf("int_compi = %.3e\n",int_comp);
// exit(0);


        if (int_comp > fsky) {
          printf("int_comp larger than fsky.\n");
          int_comp=fsky;
        }
        if (int_comp <= 0. || isinf(int_comp) || isnan(int_comp)) {
          if (ptsz->sz_verbose>3){
        if (int_comp <= 0.) printf("int_comp<0 thp = %.5e th2 = %.5e mp = %.5e zp = %.5e\n",thp,th2,mp,zp); // This is not problematic. I forgot why. theta just too large i think.
        if (isinf(int_comp)) printf("int_comp=infty.\n");
        if (isnan(int_comp)) printf("comp=nan.\n");
      }
        int_comp=1.e-300;
}

        completeness_2d[index_m][index_z] = int_comp;
        completeness_2d_to_1d[index_m_z] = log(int_comp);
        // printf("c = %.4e %.4e\n",completeness_2d_to_1d[index_m_z],int_comp);
        index_m_z += 1;

      }//end m loop

    }//end z loop

// exit(0);
      // end tabulate completeness as a function of z and m
  //  exit(0);
      if (ptsz->sz_verbose > 3)
      printf("->SZ_counts grid_C_2d debug 3.\n");

  // freeing memory
  for (index1=0;index1<ptsz->Ny;index1++)
  {
    free(erfs[index1]);
  }
  free(erfs);
  free(erfs_2d_to_1d);

  }// end else sigmaM != 0


}


    for (index_z=0;index_z<pcsz->Nbins_z;index_z++){
      double z_bin_min = pcsz->z_center[index_z]-0.5*pcsz->dz;
      double z_bin_max = pcsz->z_center[index_z]+0.5*pcsz->dz;

// if (index_z == 0){
//   z_bin_min = pcsz->z_0;
// }
// else{
//   z_bin_min = pcsz->z_center[j];
// }


// struct Parameters_for_integrand_cluster_counts_redshift V;
//   V.ptsz = ptsz;
//   V.pba = pba;
//   V.completeness_2d_to_1d = completeness_2d_to_1d;
//
//   void * params = &V;
  double r; //result of the integral
//
//   double epsrel = ptsz->redshift_epsrel_cluster_counts;
//   double epsabs = ptsz->redshift_epsabs_cluster_counts;
//   //int show_neval = ptsz->patterson_show_neval;
//
//   double z_min = z_bin_min;
//   double z_max = z_bin_max;


  // r=Integrate_using_Patterson_adaptive(z_min,
  //                                      z_max,
  //                                      epsrel, epsabs,
  //                                      integrand_cluster_counts_redshift,
  //                                      params,0);
  //
int index_z_steps_z;
int index_z_steps_z_min,index_z_steps_z_max;
index_z_steps_z_min = 0;
double test;
test = 1e100;
  for(index_z_steps_z = 0;index_z_steps_z<pcsz->nsteps_z;index_z_steps_z++){
    if (fabs(pcsz->steps_z[index_z_steps_z]-(pcsz->z_center[index_z]-0.5*pcsz->dz))<test){
      test = fabs(pcsz->steps_z[index_z_steps_z]-(pcsz->z_center[index_z]-0.5*pcsz->dz));
      index_z_steps_z_min = index_z_steps_z;
    }
}
test = 1e100;
  for(index_z_steps_z = 0;index_z_steps_z<pcsz->nsteps_z;index_z_steps_z++){
    if (fabs(pcsz->steps_z[index_z_steps_z]-(pcsz->z_center[index_z]+0.5*pcsz->dz))<test){
      test = fabs(pcsz->steps_z[index_z_steps_z]-(pcsz->z_center[index_z]+0.5*pcsz->dz));
      index_z_steps_z_max = index_z_steps_z-1;
    }
}

if (ptsz->sz_verbose>3){
if (index_y == 0)
printf("index_y = %d index_z = %d z_min = %.3e index_z_steps_z_min = %d stepz_min = %.6e\n",index_y,index_z,pcsz->z_center[index_z]-0.5*pcsz->dz,index_z_steps_z_min,pcsz->steps_z[index_z_steps_z_min]);

if (index_y == 0)
printf("index_y = %d index_z = %d z_max = %.3e index_z_steps_z_max = %d stepz_max = %.6e\n",index_y,index_z,pcsz->z_center[index_z]+0.5*pcsz->dz,index_z_steps_z_max,pcsz->steps_z[index_z_steps_z_max]);
}
r = 0.;
for (index_z_steps_z = index_z_steps_z_min;index_z_steps_z<index_z_steps_z_max+1;index_z_steps_z++){
  for (index_m = 0;index_m<pcsz->nsteps_m;index_m++){
    double zp = pcsz->steps_z[index_z_steps_z];
    double zpp = pcsz->steps_z[index_z_steps_z+1];
    double mp = exp(pcsz->steps_m[index_m]);
    double fp = get_volume_at_z(zp,pba)*get_dndlnM_at_z_and_M(zp,mp,ptsz);
    double cp = completeness_2d[index_m][index_z_steps_z];
    double fpp = get_volume_at_z(zpp,pba)*get_dndlnM_at_z_and_M(zpp,mp,ptsz);
    double cpp = completeness_2d[index_m][index_z_steps_z+1];
if (ptsz->sz_verbose>3)
    printf("index_y =  %d cp = %.5e fp = %.5e zp = %.5e  zpp = %.5e  mp = %.5e\n",index_y,cp, fp, zp, zpp, mp);
 r+= 0.5*(fp*cp+fpp*cpp)*pcsz->dlnM*(zpp-zp);
  }
}


      pcsz->dNdzdy_theoretical[index_z][index_y]=4.*_PI_*r;
    //   if (ptsz->has_completeness == 0){
    //
    //   fsky = ptsz->sky_area_deg2/41253.;
    //   pcsz->dNdzdy_theoretical[index_z][index_y]=4.*_PI_*fsky*r;
    //   }
    //   else{
    //   pcsz->dNdzdy_theoretical[index_z][index_y]=4.*_PI_*r;
    // }
    //
     }//end loop z bins for lkl

  free(completeness_2d);
  free(completeness_2d_to_1d);




  return _SUCCESS_;
}



// double integrand_cluster_counts_mass(double lnm, void *p){
//   double result = 0.;
//   struct Parameters_for_integrand_cluster_counts_mass *V = ((struct Parameters_for_integrand_cluster_counts_mass *) p);
//
//           double m_asked = exp(lnm);
//           double f1 = get_volume_at_z(V->z,V->pba)*get_dndlnM_at_z_and_M(V->z,m_asked,V->ptsz);
//           double c1 = get_completeness_at_z_and_M(V->z,m_asked,V->completeness_2d_to_1d,V->ptsz);
//           // if (pcsz->has_completeness == 0){
//           //   c1 = 1.;
//           // }
//           if (isnan(get_dndlnM_at_z_and_M(V->z,m_asked,V->ptsz))){
//             printf("z = %.3e volume = %.3e dn = %.3e c = %.3e\n",V->z,get_volume_at_z(V->z,V->pba),get_dndlnM_at_z_and_M(V->z,m_asked,V->ptsz),c1);
//             exit(0);
//             }
//   result = f1*c1;
//   return result;
// }



// double integrand_cluster_counts_completeness(double lny, void *p){
//   double result = 0.;
//   struct Parameters_for_integrand_cluster_counts_completeness *V = ((struct Parameters_for_integrand_cluster_counts_completeness *) p);
//
//           double y_asked = exp(lny);
//           double win = get_detection_proba_at_y_and_theta(y_asked,V->theta,V->erfs_2d_to_1d,V->ptsz);
//
//           // double ekl1 = get_detection_proba_at_y_and_theta(y_asked,V->theta1,V->erfs_2d_to_1d,V->ptsz);
//           // double ekl2 = get_detection_proba_at_y_and_theta(y_asked,V->theta2,V->erfs_2d_to_1d,V->ptsz);
//           // double win = ekl1+(ekl2-ekl1)/(V->theta2-V->theta1)*(V->theta-V->theta1);
//
//           double mu = log(V->y);
//           double arg=((lny-mu)/(sqrt(2.)*V->ptsz->sigmaM_ym));
//           double fac =1./sqrt(2.*_PI_*pow(V->ptsz->sigmaM_ym,2));
//
//
//   result = win*fac*exp(-arg*arg);
//   return result;
// }



//
// double integrand_cluster_counts_redshift(double z, void *p){
//   double result = 0.;
//   struct Parameters_for_integrand_cluster_counts_redshift *W = ((struct Parameters_for_integrand_cluster_counts_redshift *) p);
//
// struct Parameters_for_integrand_cluster_counts_mass V;
//   V.ptsz = W->ptsz;
//   V.pba = W->pba;
//   V.z = z;
//   V.completeness_2d_to_1d = W->completeness_2d_to_1d;
//
//   void * params = &V;
//   double r; //result of the integral
//
//   double epsrel = W->ptsz->mass_epsrel_cluster_counts;
//   double epsabs = W->ptsz->mass_epsabs_cluster_counts;
//   //int show_neval = ptsz->patterson_show_neval;
//   //
//   // double m_min = W->ptsz->M1SZ;
//   // double m_max = W->ptsz->M2SZ;
//   //
//   double m_min = exp(W->ptsz->steps_m[0]);
//   double m_max = exp(W->ptsz->steps_m[W->ptsz->nsteps_m-1]);
//   // printf("m_min = %.5e, m_max = %.5e epsrel = %.5e epsabs = %.5e\n",m_min,m_max,epsrel,epsabs);
//   //
//   r=Integrate_using_Patterson_adaptive(log(m_min),
//                                        log(m_max),
//                                        epsrel, epsabs,
//                                        integrand_cluster_counts_mass,
//                                        params,0); // 0 is show neval
//
//
//  // gsl_function F;
//  // double result_gsl, error;
//  // F.function = &integrand_cluster_counts_mass;
//  // F.params = params;
//  // int n_subintervals_gsl = 0;
//  // gsl_integration_romberg_workspace * w = gsl_integration_romberg_alloc (n_subintervals_gsl);
//  //
//  // size_t neval;
//  // double xin = log(m_min);
//  // double xout = log(m_max);
//  // gsl_integration_romberg(&F,xin,xout,epsabs,epsrel,&result_gsl,&neval,w);
//  // gsl_integration_romberg_free(w);
// // printf("result_gsl = %.5e r = %.5e\n",result_gsl,r);
//
//
//
//   // result = result_gsl;
//   result = r;
//   // result = result_gsl;
//   return result;
// }



int write_output_cluster_counts(struct szcount * pcsz, struct tszspectrum * ptsz){
int i,j;

if (ptsz->sz_verbose > 0){
  double total_counts = 0.;
  for (j=0;j<pcsz->Nbins_z;j++){
    double N_of_z = 0.;
    if (j== 0) {
        printf("log10(snr)\t",pcsz->Nbins_y);
        for (i=0;i<pcsz->Nbins_y;i++){
          printf("%.3e\t",pcsz->logy[i]);
        }
        printf(" ------ \n");
    }
      printf("z=%.3e\t",pcsz->z_center[j]);
    for (i=0;i<pcsz->Nbins_y;i++){

      N_of_z += pcsz->dNdzdy_theoretical[j][i];

      total_counts += pcsz->dNdzdy_theoretical[j][i];
      printf("%e\t",pcsz->dNdzdy_theoretical[j][i]);
    }
    printf(" ------ N_of_z = %.5e\n",N_of_z);
  }
printf("------------------------------------------------------------\n");

printf("N_of_q = \t");
for (i=0;i<pcsz->Nbins_y;i++){
double N_of_q = 0.;
for (j=0;j<pcsz->Nbins_z;j++){
  N_of_q  += pcsz->dNdzdy_theoretical[j][i];
}
printf("%.5e\t",N_of_q);

}
printf("\n");


  if (pcsz->has_completeness == 0)
    total_counts = total_counts/(pcsz->Nbins_y);

    printf("total counts = %e\n", total_counts);

  }

if (ptsz->write_sz > 0)
{
  char Filepath[_ARGUMENT_LENGTH_MAX_];
  int i,index_m,index_z;
  FILE *fp;
  sprintf(Filepath,"%s%s%s",pcsz->root,"dndzdy",".txt");
  fp=fopen(Filepath, "w");

  if(fp == NULL)
    exit(-1);
    int j;
  for (j=0;j<pcsz->Nbins_y;j++){
    for (i=0;i<pcsz->Nbins_z;i++){

      fprintf(fp,"%e\n",pcsz->dNdzdy_theoretical[i][j]);
    }}




  sprintf(Filepath,"%s%s%s",pcsz->root,"dndzdy_bins_z_center",".txt");
  fp=fopen(Filepath, "w");
  if(fp == NULL)
    exit(-1);

    for (i=0;i<pcsz->Nbins_z;i++){
      fprintf(fp,"%e\n",pcsz->z_center[i]);
    }


  sprintf(Filepath,"%s%s%s",pcsz->root,"dndzdy_bins_y_center",".txt");
  fp=fopen(Filepath, "w");
  if(fp == NULL)
    exit(-1);

    for (i=0;i<pcsz->Nbins_y;i++){
      fprintf(fp,"%e\n",pow(10.,pcsz->logy[i]));
    }


  sprintf(Filepath,"%s%s%s",pcsz->root,"dndzdm_bins_m_center",".txt");
  fp=fopen(Filepath, "w");
  if(fp == NULL)
    exit(-1);
    fprintf(fp,"# Column 1: ln(M/(Msun/h)) at the center of the mass bin.\n");
    fprintf(fp,"# Column 2: M [Msun/h]\n");
    fprintf(fp,"# Note that we use logarithmic binning between M_min = %.4e Msun/h and M_max = %.4e Msun/h.\n\n",exp(pcsz->lnM_min),exp(pcsz->lnM_max));

    for (i=0;i<pcsz->nsteps_m;i++){

      fprintf(fp,"%e\t%e\n",pcsz->steps_m[i],exp(pcsz->steps_m[i]));
    }



  printf("Output written in files\n");
  fclose(fp);

}
  return _SUCCESS_;
}



int initialise_and_allocate_memory_cc(struct tszspectrum * ptsz,struct szcount * pcsz){

  pcsz->nzSZ = ptsz->n_arraySZ_for_integral;

  ptsz->has_completeness = pcsz->has_completeness;

  //pcsz->nzSZ = 20.; //cosmomc settings
  //pcsz->size_logM = 105; //cosmomc settings

  //pcsz->rho_m_at_z = ptsz->Omega_m_0*ptsz->Rho_crit_0*pow((1.+pcsz->redshift_for_dndm),3);


  class_alloc(pcsz->redshift,sizeof(double)*pcsz->nzSZ,pcsz->error_message);




//Planck cut_off = 6.;
//SO cut_off = 5.;
// if(ptsz->experiment == 0) ptsz->sn_cutoff = 6.;
// if(ptsz->experiment == 1) ptsz->sn_cutoff = 5.;
  // ptsz->sn_cutoff = 5.;
  // pcsz->alpha;
  // //pcsz->ystar = pow(10.,pcsz->ystar)/pow(2., pcsz->alpha)*0.00472724;//8.9138435358806980e-004;
  // pcsz->beta = 0.66;
  // pcsz->thetastar = 6.997;
  // pcsz->alpha_theta = 1./3.;


  //grid for mass
  // if (pcsz->mass_range == 0){
    pcsz->lnM_max = log(ptsz->M2SZ);
    pcsz->lnM_min = log(ptsz->M1SZ);
    // printf("lnmmin = %.4e lnmmax = %.4e\n",pcsz->lnM_min,pcsz->lnM_max);
  // }
  //
  // else {
    // pcsz->lnM_max = 37.; //cosmomc/szount.f90 range
    // pcsz->lnM_min = 31.54;
  //
  // }

  pcsz->dlnM = ptsz->dlnM_cluster_count_completeness_grid; //0.05 ref value in szcounts.f90


  pcsz->nsteps_m = floor((pcsz->lnM_max - pcsz->lnM_min) /pcsz->dlnM);
  ptsz->nsteps_m = pcsz->nsteps_m;

  double lnM = pcsz->lnM_min;
  int index_m;

  class_alloc(pcsz->steps_m,
              pcsz->nsteps_m*sizeof(double *),
              pcsz->error_message
              );

  for (index_m=0; index_m<pcsz->nsteps_m; index_m++){
    pcsz->steps_m[index_m] = lnM + pcsz->dlnM/2.;
    lnM += pcsz->dlnM;
  }

  class_alloc(ptsz->steps_m,
                //ptsz->steps_m,
                ptsz->nsteps_m*sizeof(double *),
                ptsz->error_message
                );

  for (index_m=0; index_m<ptsz->nsteps_m; index_m++){
    ptsz->steps_m[index_m] = pcsz->steps_m[index_m];
  }
if (ptsz->sz_verbose>3){
    printf("steps_m [%d]:\n",pcsz->nsteps_m);
    for (index_m=0;index_m<pcsz->nsteps_m;index_m++){

      printf("%e\n",exp(pcsz->steps_m[index_m]));
    }
  }

// printf("nsteps_z=%d\n", 1);
  //grid for redshift
  //# Redshift bin parameters
  pcsz->z_0 = ptsz->bin_z_min_cluster_counts;
  pcsz->z_max = ptsz->bin_z_max_cluster_counts;
  pcsz->dz = ptsz->bin_dz_cluster_counts;

  pcsz->Nbins_z =floor((pcsz->z_max - pcsz->z_0)/pcsz->dz); // BB commented to match planck cc: "-1";

// printf("%d\n",pcsz->Nbins_z);
// exit(0);
  class_alloc(pcsz->z_center,pcsz->Nbins_z*sizeof(double),pcsz->error_message);
  int index_z;
  for (index_z = 0; index_z<pcsz->Nbins_z; index_z ++){
    pcsz->z_center[index_z] = pcsz->z_0 + 0.5*pcsz->dz + index_z*pcsz->dz;
    //printf("index_z=%d, z_center=%e\n",index_z,z_center[index_z]);
  }

  // if(pcsz->z_0==0.) pcsz->z_center[0] += 1.e-5;

  double binz=pcsz->z_center[1]-pcsz->z_center[0];


  double z_max = pcsz->z_max;//pcsz->z_center[pcsz->Nbins_z-1] + 0.5*pcsz->dz;
  double z_i = pcsz->z_0;//commented to match planck cc: + 0.5*pcsz->dz;
  pcsz->nsteps_z = 0;

// if (ptsz->sz_verbose>3)
//   printf("zmax = %.5e\n",z_max);

  while (z_i < z_max) {
    z_i = next_z(z_i,binz,ptsz);
    pcsz->nsteps_z += 1;
  }
  //commented to match planck cc: pcsz->nsteps_z += -2;

  ptsz->nsteps_z = pcsz->nsteps_z;

// exit(0);

  class_alloc(pcsz->steps_z,
              pcsz->nsteps_z*sizeof(double *),
              pcsz->error_message);

  class_alloc(ptsz->steps_z,
              //ptsz->steps_z,
              ptsz->nsteps_z*sizeof(double *),
              ptsz->error_message);

  z_i = pcsz->z_0;// commented to match planck cc: + 0.5*pcsz->dz;;

  for(index_z = 0; index_z<pcsz->nsteps_z; index_z++){
    pcsz->steps_z[index_z] = z_i;
    z_i = next_z(z_i,binz,ptsz);
  }

  if (pcsz->steps_z[0]==0) pcsz->steps_z[0] = 1.e-5;

for(index_z = 0; index_z<pcsz->nsteps_z; index_z++)
   ptsz->steps_z[index_z] = pcsz->steps_z[index_z];
  //grid for s/n

if (ptsz->sz_verbose>3){
    printf("steps_z [%d]:\n",pcsz->nsteps_z);
    for (index_z=0;index_z<pcsz->nsteps_z;index_z++){

      printf("%e\n",pcsz->steps_z[index_z]);
    }
}
// exit(0);
  //grid for y
  //# y bin parameters
  //# Logymin corresponds to log10 of S/N_min (5 or 6)
  //# Logymax corresponds to log10 of S/N_max  (~32 for planck) (but the s/n is higher in the next to last bin)
  if (ptsz->experiment==0){
  pcsz->logy_min = 0.7;
  pcsz->logy_max = 1.5;
  pcsz->dlogy = ptsz->bin_dlog10_snr;
}
else if (ptsz->experiment==1){
  pcsz->logy_min = ptsz->log10_snr_min; //log10(ptsz->sn_cutoff);
  pcsz->logy_max = ptsz->log10_snr_max;//1.8124259665302023;
  pcsz->dlogy =ptsz->bin_dlog10_snr;
}
  pcsz->Nbins_y = floor((pcsz->logy_max - pcsz->logy_min)/pcsz->dlogy)+1;
  // printf("%d\n",pcsz->Nbins_y);
  //exit(0);
  double * logy;
  class_alloc(logy,(pcsz->Nbins_y)*sizeof(double),pcsz->error_message);
  int index_y;
  double y_i = pcsz->logy_min + pcsz->dlogy/2.;
  for (index_y = 0; index_y<pcsz->Nbins_y; index_y ++){
    logy[index_y] = y_i;
    y_i += pcsz->dlogy;
    // if (y_i >= 1.5){
    //   break;
    // }
    // // printf("index_y=%d, logy=%e\n",index_y,logy[index_y]);
  }
  // pcsz->Nbins_y = index_y+1;
  // pcsz->Nbins_y = index_y;
  class_alloc(pcsz->logy,(pcsz->Nbins_y)*sizeof(double),pcsz->error_message);
  for (index_y = 0; index_y<pcsz->Nbins_y; index_y ++){
    pcsz->logy[index_y] = logy[index_y];
if (ptsz->sz_verbose>3){
    printf("index_y=%d, logy=%e\n",index_y,logy[index_y]);
  }
}
  // if (pcsz->logy_max <= pcsz->logy[pcsz->Nbins_y-1]){
  //   pcsz->logy_max = pcsz->logy[pcsz->Nbins_y-1] + pcsz->dlogy;
  // }
  ptsz->bin_dlog10_snr_last_bin = pcsz->dlogy;
  // ptsz->bin_dlog10_snr_last_bin = (pcsz->logy_max-pcsz->logy[pcsz->Nbins_y-1]);
  // pcsz->logy[pcsz->Nbins_y] = pcsz->logy[pcsz->Nbins_y-1]+pcsz->dlogy/2.+ptsz->bin_dlog10_snr_last_bin/2.;

// for (index_y = 0; index_y<pcsz->Nbins_y+1; index_y ++){
//   printf("index_y=%d, logy=%e dl=%e dllb=%e\n",index_y,pcsz->logy[index_y],pcsz->dlogy,ptsz->bin_dlog10_snr_last_bin);
// }
// exit(0);
free(logy);
  //y_500 grid
  // pcsz->lnymin = -11.5;
  // pcsz->lnymax = 10.;

  pcsz->lnymin = ptsz->lnymin;
  pcsz->lnymax = ptsz->lnymax;
  pcsz->dlny = ptsz->dlny; // 0.05 in planck

  ptsz->Ny = floor((pcsz->lnymax-pcsz->lnymin)/pcsz->dlny);

  class_alloc(ptsz->erfs_2d_to_1d_y_array,
              //ptsz->erfs_2d_to_1d_y_array,
              ptsz->Ny*sizeof(double *),
              ptsz->error_message);

  int iy;
  for (iy = 0; iy < ptsz->Ny; iy++){
    ptsz->erfs_2d_to_1d_y_array[iy] = pcsz->lnymin + iy*pcsz->dlny;
  }

  // printf("index_y=%e, logy=%e\n",
  // ptsz->erfs_2d_to_1d_y_array[0],
  // ptsz->erfs_2d_to_1d_y_array[ptsz->Ny-1]);
  // exit(0);

  class_alloc(pcsz->dNdzdy_theoretical,
              pcsz->Nbins_z*sizeof(double *),
              pcsz->error_message);


  for (index_z=0;
       index_z<pcsz->Nbins_z;
       index_z++)
  {
    class_alloc(pcsz->dNdzdy_theoretical[index_z],
                (pcsz->Nbins_y)*sizeof(double *),
                pcsz->error_message);
  }


  pcsz->index_y = 0;
  pcsz->pvecsz_size = pcsz->index_y + 1;

    return _SUCCESS_;
}


// double get_ylim_at_theta(double thp,)


int find_theta_bin(struct tszspectrum * ptsz, double thp, int * l_array, double * theta_array){
  int l1,l2;
  double th1,th2;

  if (thp > ptsz->theta_bin_max){
    l1 = ptsz->nthetas - 2;
    l2 = ptsz->nthetas - 1;
    th1 = ptsz->thetas[l1];
    th2 = ptsz->thetas[l2];
    //printf("above\n");

  }

  else if (thp < ptsz->theta_bin_min){
    l1 = 0;
    l2 = 1;
    th1 = ptsz->thetas[l1];
    th2 = ptsz->thetas[l2];
    //printf("bellow\n");
  }

  else{
    //find index where thp is closest to theta
    double dif_theta_min = fabs(ptsz->thetas[0]-thp);
    int P=0;
    int c;
    for (c = 1; c < ptsz->nthetas; c++)
    {
      if (fabs(ptsz->thetas[c] -thp)< dif_theta_min)
      {
        dif_theta_min = fabs(ptsz->thetas[c]-thp);
        P = c;
      }
    }

    l1 = P;
    th1 = ptsz->thetas[l1];
    l2 = l1 +1;
    if (thp<th1){
      l2 = l1;
      l1 = l2 -1;
    }
    th1 = ptsz->thetas[l1];
    th2 = ptsz->thetas[l2];
  }

  l_array[0] = 0;
  l_array[1] = l1;
  l_array[2] = l2;

  theta_array[0] = thp;
  theta_array[1] = th1;
  theta_array[2] = th2;

  return _SUCCESS_;
}

int find_y_bin(struct tszspectrum * ptsz, double thp, int * l_array, double * theta_array){
  int l1,l2;
  double th1,th2;

  if (thp > ptsz->erfs_2d_to_1d_y_array[ptsz->Ny-1]){

    printf("yp above y_max\n");
    l_array[0] = 0;
    l_array[1] = 0;
    l_array[2] = 0;
    exit(0);
    //printf("above\n");

  }

  else if (thp < ptsz->erfs_2d_to_1d_y_array[0]){

    // printf("yp below y_min\n");
    l_array[0] = 0;
    l_array[1] = 0;
    l_array[2] = 0;
    // exit(0);
  }

  else{
    //find index where thp is closest to theta
    double dif_theta_min = fabs(ptsz->erfs_2d_to_1d_y_array[0]-thp);
    int P=0;
    int c;
    for (c = 1; c < ptsz->Ny; c++)
    {
      if (fabs(ptsz->erfs_2d_to_1d_y_array[c] -thp)< dif_theta_min)
      {
        dif_theta_min = fabs(ptsz->erfs_2d_to_1d_y_array[c]-thp);
        P = c;
      }
    }

    l1 = P;
    th1 = ptsz->erfs_2d_to_1d_y_array[l1];
    l2 = l1 +1;
    if (thp<th1){
      l2 = l1;
      l1 = l2 -1;
    }
    th1 = ptsz->erfs_2d_to_1d_y_array[l1];
    th2 = ptsz->erfs_2d_to_1d_y_array[l2];


  l_array[0] = 0;
  l_array[1] = l1;
  l_array[2] = l2;

  theta_array[0] = thp;
  theta_array[1] = th1;
  theta_array[2] = th2;
}

  return _SUCCESS_;
}
