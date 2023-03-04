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
#include <gsl/gsl_rng.h>

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

    fftw_destroy_plan(ptsz->forward_plan_counts_fft);
    fftw_destroy_plan(ptsz->reverse_plan_counts_fft);

    free(ptsz->array_y_to_m_redshift);
    free(ptsz->array_y_to_m_y);
    free(ptsz->array_y_to_m_at_z_y);
    free(ptsz->szcounts_fft_z);
    free(ptsz->szcounts_fft_qobs);
    // free(ptsz->szcounts_fft_sigmayobs);
    // free(ptsz->szcounts_fft_nexpected_qobs);
    free(ptsz->szcounts_fft_dndzdq);
    int izsig;
    for (izsig=0;izsig<ptsz->szcounts_fft_nz;izsig++){
      // free(ptsz->szcounts_fft_index_zsig[izsig]);
      // free(ptsz->szcounts_fft_index_zq[izsig]);
      free(ptsz->szcounts_fft_index_zq_final[izsig]);
    }
    free(ptsz->szcounts_fft_index_zq_final);

    int ipatches;
    for (ipatches = 0;ipatches<ptsz->nskyfracs;ipatches++){
      free(ptsz->szcounts_fft_qmconv_all_patches[ipatches]);
    }
    free(ptsz->szcounts_fft_qmconv_all_patches);
    // free(ptsz->szcounts_fft_index_zsig);
    // free(ptsz->szcounts_fft_index_zq);

    // free(ptsz->szcounts_fft_nexpected_dndzdqgt);
//     int index_qobs;
//     for (index_qobs = 0; index_qobs<ptsz->N_samp_fftw; index_qobs++){
//     free(ptsz->szcounts_fft_rates_at_z_sigy_qobs[index_qobs]);
// }
// free(ptsz->szcounts_fft_rates_at_z_sigy_qobs);

// printf("freeing dndzdq\n");
// // free(ptsz->szcounts_fft_dndzdq);
// printf("dndzdq freed\n");

  }

  return _SUCCESS_;
}

double  get_szcounts_dndzdq_at_z_q(double z_asked, double qobs_asked, struct tszspectrum * ptsz){
  double z = z_asked;
  double nu = qobs_asked;

  // double z = log(1.+z_asked);
  // double m = log(m_asked);
   if (z<ptsz->szcounts_fft_z[0]){
      // z = ptsz->array_z_L_sat[0];
      // printf("redshift min out of range in Lsat asked %.3e bound %.3e.\n",z,ptsz->szcounts_fft_z[0]);
      // exit(0);
      return 0.;
    }
        // printf("dealing with mass conversion in hmf\n");
   if (z>ptsz->szcounts_fft_z[ptsz->szcounts_fft_nz-1]){
      // z =  ptsz->array_z_L_sat[ptsz->n_z_L_sat-1];

      // printf("redshift max out of range in Lsat asked %.3e bound %.3e.\n",z,ptsz->szcounts_fft_z[ptsz->szcounts_fft_nz-1]);
      // exit(0);
      return 0.;
    }


   if (nu<ptsz->szcounts_fft_qobs[0]){
    // nu = ptsz->array_nu_L_sat[0];
      // printf("qobs  min out of range in Lsat asked %.8e bound %.8e.\n",exp(nu),exp(ptsz->szcounts_fft_qobs[0]));
      // exit(0);
      return 0.;
  }
      // printf("dealing with mass conversion in hmf\n");
   if (nu>ptsz->szcounts_fft_qobs[ptsz->szcounts_fft_nqobs-1]){
      // nu =  ptsz->array_nu_L_sat[ptsz->n_nu_L_sat-1];
      // printf("qobs max out of range in Lsat asked %.3e bound %.3e.\n",exp(nu),exp(ptsz->szcounts_fft_qobs[ptsz->szcounts_fft_nqobs-1]));
      // exit(0);
      return 0.;
    }

 double result  = pwl_interp_2d(ptsz->szcounts_fft_nqobs,
                                ptsz->szcounts_fft_nz,
                                ptsz->szcounts_fft_qobs,
                                ptsz->szcounts_fft_z,
                                ptsz->szcounts_fft_dndzdq,
                                1,
                                &nu,
                                &z);

// double result  = pwl_interp_2d(
//                                ptsz->szcounts_fft_nz,
//                                ptsz->szcounts_fft_nqobs,
//                                ptsz->szcounts_fft_z,
//                                ptsz->szcounts_fft_qobs,
//                                ptsz->szcounts_fft_dndzdq,
//                                1,
//                                &z,
//                                &nu);

 return result;

}

//
// double  get_szcounts_dndzdqgt_at_z_q(double z_asked, double qobs_asked, struct tszspectrum * ptsz){
//   double z = z_asked;
//   double nu = qobs_asked;
//
//   // double z = log(1.+z_asked);
//   // double m = log(m_asked);
//    if (z<ptsz->szcounts_fft_z[0]){
//       // z = ptsz->array_z_L_sat[0];
//       printf("redshift min out of range in Lsat asked %.3e bound %.3e.\n",z,ptsz->szcounts_fft_z[0]);
//       exit(0);
//     }
//         // printf("dealing with mass conversion in hmf\n");
//    if (z>ptsz->szcounts_fft_z[ptsz->szcounts_fft_nz-1]){
//       // z =  ptsz->array_z_L_sat[ptsz->n_z_L_sat-1];
//
//       printf("redshift max out of range in Lsat asked %.3e bound %.3e.\n",z,ptsz->szcounts_fft_z[ptsz->szcounts_fft_nz-1]);
//       exit(0);
//     }
//
//
//    if (nu<ptsz->szcounts_fft_nexpected_qobs[0]){
//     // nu = ptsz->array_nu_L_sat[0];
//       printf("qobs  min out of range in Lsat asked %.8e bound %.8e.\n",exp(nu),exp(ptsz->szcounts_fft_nexpected_qobs[0]));
//       exit(0);
//   }
//       // printf("dealing with mass conversion in hmf\n");
//    if (nu>ptsz->szcounts_fft_nexpected_qobs[ptsz->szcounts_fft_nexpected_qobs_n-1]){
//       // nu =  ptsz->array_nu_L_sat[ptsz->n_nu_L_sat-1];
//       printf("qobs max out of range in Lsat asked %.3e bound %.3e.\n",exp(nu),exp(ptsz->szcounts_fft_nexpected_qobs[ptsz->szcounts_fft_nexpected_qobs_n-1]));
//       exit(0);
//     }
//
//  double result  = pwl_interp_2d(ptsz->szcounts_fft_nexpected_qobs_n,
//                                 ptsz->szcounts_fft_nz,
//                                 ptsz->szcounts_fft_nexpected_qobs,
//                                 ptsz->szcounts_fft_z,
//                                 ptsz->szcounts_fft_nexpected_dndzdqgt,
//                                 1,
//                                 &nu,
//                                 &z);
//
//  // double result  = pwl_interp_2d(ptsz->szcounts_fft_nexpected_qobs_n,
//  //                                ptsz->szcounts_fft_nz,
//  //                                ptsz->szcounts_fft_nexpected_qobs,
//  //                                ptsz->szcounts_fft_z,
//  //                                ptsz->szcounts_fft_nexpected_dndzdqgt,
//  //                                1,
//  //                                &nu,
//  //                                &z);
//
//  return result;
//
// }

//
// double  get_szcounts_rates_at_z_sigobs_qobs(double z_asked, double sig_asked, double qobs_asked, struct tszspectrum * ptsz){
//   double z = z_asked;
//   double m = log(sig_asked);
//   double nu = qobs_asked;
//
//   // printf("nu asked = %.3e\n",nu_asked);
//   // exit(0);
//
//
//   // double z = log(1.+z_asked);
//   // double m = log(m_asked);
//    if (z<ptsz->szcounts_fft_z[0]){
//       // z = ptsz->array_z_L_sat[0];
//       printf("redshift min out of range in Lsat asked %.3e bound %.3e.\n",z,ptsz->szcounts_fft_z[0]);
//       exit(0);
//     }
//         // printf("dealing with mass conversion in hmf\n");
//    if (z>ptsz->szcounts_fft_z[ptsz->szcounts_fft_nz-1]){
//       // z =  ptsz->array_z_L_sat[ptsz->n_z_L_sat-1];
//
//       printf("redshift max out of range in Lsat asked %.3e bound %.3e.\n",z,ptsz->szcounts_fft_z[ptsz->szcounts_fft_nz-1]);
//       exit(0);
//     }
//
//    if (m<ptsz->szcounts_fft_sigmayobs[0]){
//     // m = ptsz->array_m_L_sat[0];
//       printf("mass min out of range in Lsat asked %.3e bound %.3e.\n",m,ptsz->szcounts_fft_sigmayobs[0]);
//       exit(0);
//   }
//       // printf("dealing with mass conversion in hmf\n");
//    if (m>ptsz->szcounts_fft_sigmayobs[ptsz->szcounts_fft_nsigmayobs-1]){
//       // m =  ptsz->array_m_L_sat[ptsz->n_m_L_sat-1];
//       printf("mass max out of range in Lsat asked %.3e bound %.3e.\n",m,ptsz->szcounts_fft_sigmayobs[ptsz->szcounts_fft_nsigmayobs-1]);
//       exit(0);
//     }
//
//    if (nu<ptsz->szcounts_fft_qobs[0]){
//     // nu = ptsz->array_nu_L_sat[0];
//       printf("freq min out of range in Lsat asked %.8e bound %.8e.\n",exp(nu),exp(ptsz->szcounts_fft_qobs[0]));
//       exit(0);
//   }
//       // printf("dealing with mass conversion in hmf\n");
//    if (nu>ptsz->szcounts_fft_qobs[ptsz->szcounts_fft_nqobs-1]){
//       // nu =  ptsz->array_nu_L_sat[ptsz->n_nu_L_sat-1];
//       printf("freq max out of range in Lsat asked %.3e bound %.3e.\n",exp(nu),exp(ptsz->szcounts_fft_qobs[ptsz->szcounts_fft_nqobs-1]));
//       exit(0);
//     }
//
//   // if (ptsz->tau_profile == 1){
//   // find the closest l's in the grid:
//   int id_l_low;
//   int id_l_up;
//   int n_nu = ptsz->szcounts_fft_nqobs;
//   int n_m = ptsz->szcounts_fft_nsigmayobs;
//   int n_z = ptsz->szcounts_fft_nz;
//   r8vec_bracket(n_nu,ptsz->szcounts_fft_qobs,nu,&id_l_low,&id_l_up);
//
//   // interpolate 2d at l_low:
//
//  double ln_rho_low = pwl_interp_2d(n_m,
//                                 n_z,
//                                 ptsz->szcounts_fft_sigmayobs,
//                                 ptsz->szcounts_fft_z,
//                                 ptsz->szcounts_fft_rates_at_z_sigy_qobs[id_l_low-1],
//                                 1,
//                                 &m,
//                                 &z);
//
//  double ln_rho_up = pwl_interp_2d(n_m,
//                                 n_z,
//                                 ptsz->szcounts_fft_sigmayobs,
//                                 ptsz->szcounts_fft_z,
//                                 ptsz->szcounts_fft_rates_at_z_sigy_qobs[id_l_up-1],
//                                 1,
//                                 &m,
//                                 &z);
//  double ln_l_low = ptsz->szcounts_fft_qobs[id_l_low-1];
//  double ln_l_up = ptsz->szcounts_fft_qobs[id_l_up-1];
//
//  // printf("lnrh %.5e %.5e %d %d\n",ln_rho_low,ln_rho_up,id_l_low,id_l_up);
//
//  return ln_rho_low + ((nu - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low);
//  // return ln_rho_low + ((l - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low);
//
//
// }


// thanks chatgpt
void shiftArray(double arr[], int size, int s) {
    double temp[s];
    int i;
    for (i = 0; i < s; i++) {
        temp[i] = arr[i];
    }
    for (i = s; i < size; i++) {
        arr[i - s] = arr[i];
    }
    for (i = 0; i < s; i++) {
        arr[size - s + i] = temp[i];
    }
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
// double z = 0.3;
// double m = 5e14;
// double y = get_y_at_m_and_z(m,z,ptsz,pba);
// printf("z = %.5e m = %.5e y = %.5e\n",z,m,y);
// tabulate_y_to_m(pba,pnl,ppm,ptsz);
// double mrec;
// mrec = get_y_to_m_at_z_and_y(z,y,ptsz);
// printf("z = %.5e m = %.5e y = %.5e mrec = %.6e\n",z,m,y,mrec);
// double der = get_dlnm_dlny(log(y),z,ptsz);
// printf("z = %.5e m = %.5e y = %.5e mrec = %.6e der = %.5e\n",z,m,y,mrec,der);
// double dNdlny = get_dNdlny_at_z_and_y(z,y,pba,ptsz);
// printf("z = %.5e m = %.5e y = %.5e mrec = %.6e der = %.5e dNdlny = %.5e\n",z,m,y,mrec,der,dNdlny);

// exit(0);

    if(ptsz->sz_verbose>1) printf("constructing fftw plan\n");
    fftw_plan forward_plan, reverse_plan;
    // ptsz->N_samp_fftw = 100;
    fftw_complex* a_tmp;
    fftw_complex* b_tmp;
    a_tmp = fftw_alloc_complex(2*ptsz->N_samp_fftw);
    b_tmp = fftw_alloc_complex(2*ptsz->N_samp_fftw);
    ptsz->forward_plan_counts_fft = fftw_plan_dft_1d(2*ptsz->N_samp_fftw,
                                    (fftw_complex*) a_tmp,
                                    (fftw_complex*) b_tmp,
                                    -1, FFTW_ESTIMATE);
    ptsz->reverse_plan_counts_fft = fftw_plan_dft_1d(2*ptsz->N_samp_fftw,
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
//
// int index_parallel_zp = 0;
// int index_parallel_sigmayobsp = 0;
//
// int abort;
//
// #ifdef _OPENMP
// double tstart, tstop;
// #endif
//
// abort = _FALSE_;
// /* number of threads (always one if no openmp) */
// int number_of_threads= 1;
// #ifdef _OPENMP
// #pragma omp parallel
//   {
//     number_of_threads = omp_get_num_threads();
//     //omp_set_num_threads(number_of_threads);
//   }
// #endif
//
// int id;
// omp_lock_t lock;
//
//
// #pragma omp parallel \
//    shared(abort,pba,ptsz,ppm,pnl,lock,forward_plan,reverse_plan)\
//    private(id,index_parallel_zp,index_parallel_sigmayobsp)\
//    num_threads(number_of_threads)
// 	 {
//
// #ifdef _OPENMP
// 	   tstart = omp_get_wtime();
// #endif
//
// // #pragma omp for schedule (dynamic)
// // for (index_parallel_zp=0;index_parallel_zp<1;index_parallel_zp++)
// // 	     {
// // #pragma omp flush(abort)
//
// #pragma omp for collapse(2)
// // for (index_parallel_zp=0; index_parallel_zp<1; index_parallel_zp++)
// for (index_parallel_zp=0; index_parallel_zp<ptsz->szcounts_fft_nz; index_parallel_zp++)
// {
// // for (index_parallel_sigmayobsp=0; index_parallel_sigmayobsp<1; index_parallel_sigmayobsp++)
// for (index_parallel_sigmayobsp=0; index_parallel_sigmayobsp<ptsz->szcounts_fft_nsigmayobs; index_parallel_sigmayobsp++)
//   {
//
// // double zp = 0.3;//ptsz->szcounts_fft_z[index_parallel_zp];
// // double sigmay_obsp = 5e-3;//exp(ptsz->szcounts_fft_sigmayobs[index_parallel_sigmayobsp]);
// double zp =ptsz->szcounts_fft_z[index_parallel_zp];
// double sigmay_obsp = exp(ptsz->szcounts_fft_sigmayobs[index_parallel_sigmayobsp]);
//
// if (ptsz->sz_verbose>10)
// printf("#########  zp = %.5e sigp = %.5e\n",zp,sigmay_obsp);
//
// int index_z_sigma = ptsz->szcounts_fft_index_zsig[index_parallel_zp][index_parallel_sigmayobsp];
// // double *in,*out;
// // fftw_plan myplan;
// int N = ptsz->N_samp_fftw;
// double a = 2 * M_PI/5.;
//
// double z_fft = zp; // pick a redshift just for testing
// //
//
// // double sigma = 1.;
// double sigma = ptsz->sigmaM_ym;
//
// double sig2 = pow(sigma,2.);
// double fac =1./sqrt(2.*_PI_*sig2);
//
// // double lny_fft[N];
// // double klny_fft[N];
// double lnymin_fft =  ptsz->lnymin;//ptsz->lnymin;
// double lnymax_fft =  ptsz->lnymax;//ptsz->lnymax;
// // double dNdlny_fft[N];
// // double lognormal_fft[N];
// double xarr[2*N];
//
// fftw_complex in[2*N], out[2*N], in2[2*N],test[2*N],out_test[2*N]; /* double [2] */
// fftw_complex product[2*N],product_out[2*N];
//
//
// // printf("preparing arrays\n");
//
// // the kernel is centered at 0, we need to pad the arrays.
// // so we preserve symmetry around 0.
//
// double lnymin_fft_padded = ptsz->lnymin-3.;
// double lnymax_fft_padded = -lnymin_fft_padded; //lnymin is negative
//
//
// double L = (lnymax_fft_padded-lnymin_fft_padded);
// double dx = L/(N);
// int i;
//   for (i = 0; i < N; i++) {
//     double x = lnymin_fft_padded+i*dx;
//     // lny_fft[i] = x;
//     xarr[i] = x;
//     xarr[i+N] = lnymax_fft_padded+i*dx;
// //
// //     dNdlny_fft[i] = get_dNdlny_at_z_and_y(z_fft,exp(lny_fft[i]),pba,ptsz);
// //     in_dNdlny_fft[i] = dNdlny_fft[i];
// //
// //
// //     // double arg0 = lny_fft[i]/(sqrt(2.)*ptsz->sigmaM_ym);
//     double arg0 = (x)/sqrt(2.*sig2);///(sqrt(2.));
//
//     // in[i][0]= lognormal_fft[i];
//     // if (x==0)
//     //   in[i][0]=1.;
//     // else
//     //   in[i][0]=sin(M_PI*x)/(M_PI*x);
//     // double arg0 = x/sigma;
//     in[i][0] = fac*exp(-arg0*arg0)/exp(x); // divided by exp(x) for lognormal
//     //fac*exp(-0.5*arg0*arg0);//cos(3 * 2*M_PI*i/N);
//     // in[i][0]= cos(3 * 2*M_PI*i/N);//fac*exp(-arg0*arg0);//cos(3 * 2*M_PI*i/N);
//     in[i][1]= 0.;
//
//     in[N+i][0]= 0.;
//     in[N+i][1]= 0.;
// //
//     // test[i][0] = pow(x,3)+pow(x,2);
//     // double w = 5.;
//     // test[i][0] = pow(1.+sin(a * (x+19.)),2.)*0.5*(1.-tanh((x+19.-w/0.3)))*0.5*(1.+tanh((x+19.+w/2.)));
//     // test[i][0] = sin(a * x);//*0.5*(1.-tanh((x-w/2.)))*0.5*(1.+tanh((x+w/2.)));
//     // if (x==0)
//     //   test[i][0] = 2.;
//     // else
//     //   test[i][0] = 2.*sin(2.*M_PI*x)/(2.*M_PI*x);
//
//
//     test[i][0] = get_dNdlny_at_z_and_y(z_fft,exp(x),pba,ptsz);
//     // printf("lny %.5e test %.5e\n",x,test[i][0]);
//     if (isinf(test[i][0]))
//       exit(0);
//
//     test[i][1] = 0.;
//     //
//     test[N+i][0] = 0.;
//     test[N+i][1] = 0.;
//
//   }
//
//   // exit(0);
//
//
//   id = omp_get_thread_num();
//
//   fftw_execute_dft(forward_plan, (fftw_complex*) in, (fftw_complex*) out);
//   fftw_execute_dft(forward_plan, (fftw_complex*) test, (fftw_complex*) out_test);
// // printf("ffts done\n");
//
//   // exit(0);
// //   // fftw_execute_dft(reverse_plan, (fftw_complex*) out, (fftw_complex*) out);
// //
// // product[0][0]= out[][0]*out_test[i][0]
// for (i = 0; i < 2*N; i++){
//    product[i][0]= (out[i][0]*out_test[i][0]-out[i][1]*out_test[i][1])/(double)(2*N);
//    product[i][1]= (out[i][0]*out_test[i][1]+out[i][1]*out_test[i][0])/(double)(2*N);
//     }
// //
// //
//   fftw_execute_dft(reverse_plan, (fftw_complex*) product, (fftw_complex*) product_out);
//
//     for(i = 0; i < 2*N; i++){
//         // rout[i] = creal(out[i])*1./N;
//         // product_out[i][0] = product_out[i][0]/(double)(N);
//         // product_out[i][1] = product_out[i][1]/(double)(N);
//         product_out[i][0] *= dx;
// if(ptsz->sz_verbose>10)
//         printf("fft thread %d i = %d r = %.5e\n",id,i,product_out[i][0]);
//
//         // rin2[i] = creal(in2[i])*1./N;;
//       }
//
//       // double tmp;
//
//       double result_conv[2*N];
//
//       for(i = 0; i < 2*N; i++) {
//         // int shift = int(N/2);
//         //   tmp = product_out[i][0];
//           result_conv[i] = product_out[i][0];// =product_out[shift+i][0];
//
//       }
//
//       int shift = (int) (N/2);
//       shiftArray(result_conv, 2*N, shift);
// if (ptsz->sz_verbose>10){
//   FILE *fp;
//   char Filepath[_ARGUMENT_LENGTH_MAX_];
//   sprintf(Filepath,"%s%s%s",ptsz->root,"test_ffts_complex",".txt");
//   fp=fopen(Filepath, "w");
//   // fftw_execute_dft(forward_plan, (fftw_complex*) in_lognormal_fft, (fftw_complex*) out_lognormal_fft);
//   for (i = 0; i < 2*N; i++){
//   fprintf(fp,"%.5e\t%.5e\t%.5e\t%.5e\n",xarr[i],result_conv[i],test[i][0],in[i][0]);
//   }
//   fclose(fp);
// }
//   //
//   // exit(0);
//
//
//   // at this stage we have tilde F(lny,sigma_lny):
//   // This function/data is stored in the array result_conv.
//   // Now we need to convolve it with the obs kernel.
// if (ptsz->sz_verbose>10)
// printf("done with the scatter convolution.\n");
// // exit(0);
//   // example:
//   //z = 3.00000e-01 m = 5.00000e+14 y = 1.37830e-03 mrec = 4.999999e+14 der = 5.61798e-01
//   // typically for planck, sigmyobs varies between: 5e-5 and 1e-2
// double sigmay_obs =  sigmay_obsp; // one of the sigma values
//
//
// double qp[2*N];
// // double function1[2*N];
// // double function2[2*N];
// // double out1[2*N];
// // double out2[2*N];
// // double prod[2*N];
// // double prod_out[2*N];
//
// double qmin_fft_padded = -50;
// double qmax_fft_padded = 50;
//
// double L_q = (qmax_fft_padded-qmin_fft_padded);
// double dq = L_q/(N);
// fac = 1./sqrt(2.*_PI_);
// // int i;
// for (i = 0; i < N; i++) {
//  double x = qmin_fft_padded+i*dq;
//  qp[i] = x;
//  qp[i+N] = qmax_fft_padded+i*dq;
//  double arg0 = x/sqrt(2.);
//  in[i][0] = fac*exp(-arg0*arg0);
//  in[i][1]= 0.;
//  in[N+i][0]= 0.;
//  in[N+i][1]= 0.;
//
//    double lnyp =  log(x*sigmay_obs);
//    double conv1;
//    if (lnyp<xarr[0])
//     conv1 = 0.;
//    else if (lnyp>xarr[2*N-1])
//     conv1 = 0.;
//    else
//     conv1 =  pwl_value_1d(2*N,
//                           xarr,
//                           result_conv,
//                           lnyp);
//  test[i][0] = conv1;
//  test[i][1] = 0.;
//  test[N+i][0] = 0.;
//  test[N+i][1] = 0.;
// }
//
// fftw_execute_dft(forward_plan, (fftw_complex*) in, (fftw_complex*) out);
// fftw_execute_dft(forward_plan, (fftw_complex*) test, (fftw_complex*) out_test);
//
// for (i = 0; i < 2*N; i++){
//    product[i][0]= (out[i][0]*out_test[i][0]-out[i][1]*out_test[i][1])/(double)(2*N);
//    product[i][1]= (out[i][0]*out_test[i][1]+out[i][1]*out_test[i][0])/(double)(2*N);
//     }
//
// fftw_execute_dft(reverse_plan, (fftw_complex*) product, (fftw_complex*) product_out);
//
//     for(i = 0; i < 2*N; i++){
//         product_out[i][0] *= dq;
//       if (ptsz->sz_verbose>10)
//         printf("fft thread %d i = %d r = %.5e\n",id,i,product_out[i][0]);
//       }
//
// double result_qconv[2*N];
//   for(i = 0; i < 2*N; i++) {
//           result_qconv[i] = product_out[i][0];
//
//
//       }
// shiftArray(result_qconv, 2*N, shift);
//
// if (ptsz->sz_verbose>=1){
// for (i = 0; i < N; i++){
//   if ((qp[i]>0.7) && (qp[i]<0.9) && (fabs(zp-0.3)<0.05))
//     printf("fft thread %d i = %d z = %.3e sigma = %.3e q = %.5e r = %.5e\n",id,i,zp,sigmay_obs,qp[i],sigmay_obs*result_qconv[i]);
// }
// }
//
//
// if (ptsz->sz_verbose>10){
//   FILE *fp;
//   char Filepath[_ARGUMENT_LENGTH_MAX_];
//   sprintf(Filepath,"%s%s%s",ptsz->root,"test_ffts_q_complex",".txt");
//   fp=fopen(Filepath, "w");
//   // fftw_execute_dft(forward_plan, (fftw_complex*) in_lognormal_fft, (fftw_complex*) out_lognormal_fft);
//   for (i = 0; i < 2*N; i++){
//   fprintf(fp,"%.5e\t%.5e\t%.5e\t%.5e\n",qp[i],result_qconv[i],test[i][0],in[i][0]);
//   }
//   fclose(fp);
// }
// // exit(0);
//
// /* start fftlog stuff
// // Here we use the FFTLog algorithm
// // because the y values are in a log grid
// double rp[N];//the log-spaced y_values
// double function1[N]; //the 1st function
// double function2[N]; // the second function (kernel)
// double out1[N]; // the fourier transform of the second function
// double out2[N]; // the fourier transform of the first function
// double kp[N]; // the frequency grid
// double prod[N];
// double prod_out[N];
//
//
// printf("allocating arrays for fftlogin'\n");
// double lnymin_fftlog = ptsz->lnymin-1.;
// double lnymax_fftlog = ptsz->lnymax+3.;
//
// // double lnqmin_fftlog = log(1e-2);
// // double lnqmax_fftlog = log()
// double Lfftlog = (lnymax_fftlog-lnymin_fftlog);
// double dxfftlog = Lfftlog/(N);
// for (i=0;i<N;i++){
//   double lnyp = lnymin_fftlog+i*dxfftlog;
//   rp[i] = exp(lnyp);
//
//   // printf("rp = %.5e\n",rp[i]);
//
//   double conv1 =  pwl_value_1d(2*N,
//                           xarr,
//                           result_conv,
//                           lnyp);
//
//   function1[i] = conv1;
//   function2[i] = 1./sqrt(2.*_PI_*sigmay_obs*sigmay_obs)
//                  *exp(-0.5*pow(rp[i]/sigmay_obs,2.));
// }
// // exit(0);
// printf("done allocating arrays for fftlogin'\n");
//   // Compute the function
//   // *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
//   // * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
//   // * in this notation.  The input k-values must be logarithmically spaced.  The
//   // * resulting xi_l^m(r) will be evaluated at the dual r-values
//   // *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0].
//   // void fftlog_ComputeXiLM(int l, int m, int N, const double k[],  const double pk[], double r[], double xi[]);
//   // here is what we need:
//   // int l = -1
//   // int m = 0
//
//   fftlog_ComputeXiLMsloz(-1, 0, N, rp,  function1, kp, out1,ptsz);
//
// printf("computed fftlog1.\n");
//
//   fftlog_ComputeXiLMsloz(-1, 0, N, rp,  function2, kp, out2,ptsz);
// for (i=0;i<N;i++){
//   // prod[i] = (2.*M_PI*kp[i]*out1[i])*(2.*M_PI*kp[i]*out2[i]);
//   // prod[i] = (2.*M_PI*kp[i]*out1[i]);//*(2.*M_PI*kp[i]*out2[i]);
//   prod[i] = (out1[i]);
// }
//
//  fftlog_ComputeXiLMsloz(-1, 0, N, kp, prod, rp, prod_out,ptsz);
//
//  sprintf(Filepath,"%s%s%s",ptsz->root,"test_fftlog",".txt");
//  fp=fopen(Filepath, "w");
//  // fftw_execute_dft(forward_plan, (fftw_complex*) in_lognormal_fft, (fftw_complex*) out_lognormal_fft);
//  for (i = 0; i < N; i++){
//
//  // prod_out[i] = (2.*M_PI*rp[i]*prod_out[i]);
//  prod_out[i] = (prod_out[i]);
//
//  printf("fftlog thread %d i = %d r = %.5e\n",id,i,prod_out[i]);
//
//  fprintf(fp,"%.5e\t%.5e\t%.5e\t%.5e\n",rp[i],prod_out[i],function1[i],function2[i]);
//  }
//  fclose(fp);
//
//  */ // end fftlog stuff
//
//
//  int index_qobs;
//  for (index_qobs = 0; index_qobs<N; index_qobs++){
//    ptsz->szcounts_fft_rates_at_z_sigy_qobs[index_qobs][index_z_sigma] = sigmay_obs*result_qconv[index_qobs];
//    ptsz->szcounts_fft_qobs[index_qobs] = qp[index_qobs];
//  }
//           }
//         }
// #ifdef _OPENMP
//       tstop = omp_get_wtime();
//       if (ptsz->sz_verbose > 0)
//          printf("In %s: time spent in parallel region (loop over cluster counts ffts) = %e s for thread %d\n",
//                    __func__,tstop-tstart,omp_get_thread_num());
//
//
// #endif
//    // free(Pvecback);
//    // free(Pvectsz);
//    // free(b_l1_l2_l_1d);
// 	} //end of parallel region
//
//    if (abort == _TRUE_) return _FAILURE_;
//
//
//   // exit(0);
//
//
// if (ptsz->sz_verbose >=1){
// double test_counts = 0.;
// double ztest = 0.3;
// double sigtest = 5e-3;
// double qtest;// = 1.2;
// int index_qtest;
// int Nqtest = 50;
// for (index_qtest = 0;index_qtest<Nqtest;index_qtest++){
// qtest = index_qtest*20./Nqtest;
// test_counts = get_szcounts_rates_at_z_sigobs_qobs(ztest,sigtest,qtest,ptsz);
// printf("ztest = %.5e sigtest = %.5e qtest = %.5e rate = %.5e\n",
//         ztest,sigtest,qtest,test_counts);
// }
// }
//
// if (ptsz->sz_verbose >=1)
//   printf("we now computed the quantities necessary for nexpected.\n");
//
//
// // first we create a parallel region.
//   int index_parallel_qp = 0;
//
//   // int abort;
//
//   // #ifdef _OPENMP
//   // // double tstart, tstop;
//   // #endif
//
//   abort = _FALSE_;
//   /* number of threads (always one if no openmp) */
//   number_of_threads= 1;
//   #ifdef _OPENMP
//   #pragma omp parallel
//     {
//       number_of_threads = omp_get_num_threads();
//       //omp_set_num_threads(number_of_threads);
//     }
//   #endif
//
//   // int id;
//   // omp_lock_t lock;
//
//
//   #pragma omp parallel \
//      shared(abort,pba,ptsz,ppm,pnl,lock,forward_plan,reverse_plan)\
//      private(id,index_parallel_zp,index_parallel_qp)\
//      num_threads(number_of_threads)
//   	 {
//
//   #ifdef _OPENMP
//   	   tstart = omp_get_wtime();
//   #endif
//
//
// #pragma omp for collapse(2)
// // for (index_parallel_zp=0; index_parallel_zp<1; index_parallel_zp++)
// for (index_parallel_zp=0;
//      index_parallel_zp<ptsz->szcounts_fft_nz;
//      index_parallel_zp++)
// {
// // for (index_parallel_sigmayobsp=0; index_parallel_sigmayobsp<1; index_parallel_sigmayobsp++)
// for (index_parallel_qp=0;
//      index_parallel_qp<ptsz->szcounts_fft_nexpected_qobs_n;
//      index_parallel_qp++)
//   {
//
//
//     double zp =ptsz->szcounts_fft_z[index_parallel_zp];
//     double qp = ptsz->szcounts_fft_nexpected_qobs[index_parallel_qp];
//
//     // zp = 0.3; // test
//     // qp = 5.2;
//
//
//     id = omp_get_thread_num();
//     int N = ptsz->N_samp_fftw;
//     double result_ymconv_all_patches[2*N];
//     int i;
//     for (i = 0; i < N; i++) {
//         result_ymconv_all_patches[i] = 0.;
//       }
//
//   double lnymin_fft =  ptsz->lnymin;//ptsz->lnymin;
//   double lnymax_fft =  ptsz->lnymax;//ptsz->lnymax;
//   // the kernel is centered at 0, we need to pad the arrays.
//   // so we preserve symmetry around 0.
//   double lnymin_fft_padded = ptsz->lnymin-3.;
//   double lnymax_fft_padded = -lnymin_fft_padded; //lnymin is negative
//
//
//   double L = (lnymax_fft_padded-lnymin_fft_padded);
//   double dx = L/(N);
//
//
//     int index_patches = 0;
//     for (index_patches =0;
//          index_patches<1;//ptsz->nskyfracs;
//          index_patches++){
//
//     // printf("thread %d index_patches = %d zp = %.5e qp %.5e\n",id,index_patches,zp,qp);
//     index_patches = 5; // text
//
//
//
//     double sigma = ptsz->sigmaM_ym;
//
//     double sig2 = pow(sigma,2.);
//     double fac =1./sqrt(2.*_PI_*sig2);
//
//
//     double xarr[2*N];
//
//     fftw_complex in[2*N], out[2*N], in2[2*N],test[2*N],out_test[2*N]; /* double [2] */
//     fftw_complex product[2*N],product_out[2*N];
//
//
//
//
//
//
//   for (i = 0; i < N; i++) {
//     double x = lnymin_fft_padded+i*dx;
//     xarr[i] = x;
//     xarr[i+N] = lnymax_fft_padded+i*dx;
//
//     double arg0 = (x)/sqrt(2.*sig2);
//     in[i][0] = fac*exp(-arg0*arg0)/exp(x); // divided by "y" as it is a lognormal
//     in[i][1]= 0.;
//
//     in[N+i][0]= 0.;
//     in[N+i][1]= 0.;
//
//     test[i][0] = get_dNdlny_at_z_and_y(zp,exp(x),pba,ptsz);
//     // get the Cumulative:
//     double thetap = get_theta_at_y_and_z(exp(x),zp,ptsz,pba);
//     double sn = get_szcountsz_sigma_at_theta_in_patch(thetap,index_patches,ptsz);
//     double arg = (y - qp * sn)/(sqrt(2.) * sn);
//     double phi = 0.5*(1.+erf(arg));
//
//     test[i][0] *= phi;
//
//     if (isinf(test[i][0]))
//       exit(0);
//     test[i][1] = 0.;
//     //
//     test[N+i][0] = 0.;
//     test[N+i][1] = 0.;
//
//   }
//
//
//
//     fftw_execute_dft(forward_plan, (fftw_complex*) in, (fftw_complex*) out);
//     fftw_execute_dft(forward_plan, (fftw_complex*) test, (fftw_complex*) out_test);
//
//   for (i = 0; i < 2*N; i++){
//      product[i][0]= (out[i][0]*out_test[i][0]-out[i][1]*out_test[i][1])/(double)(2*N);
//      product[i][1]= (out[i][0]*out_test[i][1]+out[i][1]*out_test[i][0])/(double)(2*N);
//       }
//
//
//     fftw_execute_dft(reverse_plan, (fftw_complex*) product, (fftw_complex*) product_out);
//
//     for(i = 0; i < 2*N; i++){
//         product_out[i][0] *= dx;
//       // if (ptsz->sz_verbose>=1)
//       //   printf("fft nexpected thread %d i = %d r = %.5e\n",id,i,product_out[i][0]);
//       }
//
//   double result_ymconv[2*N];
//   for(i = 0; i < 2*N; i++) {
//           result_ymconv[i] = product_out[i][0];
//           result_ymconv[i] *= 4.*_PI_*ptsz->skyfracs[index_patches];
//       }
//   int shift = (int) (N/2);
//   shiftArray(result_ymconv, 2*N, shift);
//
//   for(i = 0; i < 2*N; i++) {
//   result_ymconv_all_patches[i] += result_ymconv[i];
//
//   }
//
//
// } // end loop over patches
//
//
//
// // now we integrate over lny:
//
//
// // double lnymin_fft =  ptsz->lnymin;//ptsz->lnymin;
// // double lnymax_fft =  ptsz->lnymax;//ptsz->lnymax;
// double xarr[2*N];
// // int i;
// double int_conv_lny = 0.;
// for (i = 0; i < N; i++) {
//   double x = lnymin_fft_padded+i*dx;
//   xarr[i] = x;
//   if ((x>ptsz->lnymin) || (x<ptsz->lnymax))
//     int_conv_lny += result_ymconv_all_patches[i]*dx;
// }
//
//
// // printf("thread %d zp = %.5e (%d) qp = %.5e (%d) integral_conv_lny = %.5e\n",
// //        id,zp,index_parallel_zp,
// //        qp,index_parallel_qp,
// //        int_conv_lny);
//
//  int index_zq = ptsz->szcounts_fft_index_zq[index_parallel_zp][index_parallel_qp];
//  ptsz->szcounts_fft_nexpected_dndzdqgt[index_zq] = int_conv_lny;
//
//
//
// }// end loop2
// } // end loop 1
// #ifdef _OPENMP
// tstop = omp_get_wtime();
// if (ptsz->sz_verbose > 0)
// printf("In %s: time spent in parallel region (loop over cluster counts ffts) = %e s for thread %d\n",
//          __func__,tstop-tstart,omp_get_thread_num());
//
//
// #endif
//
// } //end of parallel region
//
// if (abort == _TRUE_) return _FAILURE_;
//
//
//
// // double qtest = 6.83333e+00;
// // double ztest = 5.75000e-01; //integral_conv_lny = 2.02701e-09
// // double rtest = get_szcounts_dndzdqgt_at_z_q(ztest,qtest,ptsz);
// // printf("got %.5e expected %.5e\n",rtest,2.02701e-09);

// printf("starting the real stuff:\n");

// FILE *fp;
// char Filepath[_ARGUMENT_LENGTH_MAX_];

// set up arrays used for the convolutions:


int npatches;
if (ptsz->use_skyaveraged_noise)
    npatches = 1;
else
    npatches = ptsz->nskyfracs;
int N = ptsz->N_samp_fftw;
double lnqmin_fft = ptsz->szcounts_lnqmin_fft;
double lnqmax_fft = ptsz->szcounts_lnqmax_fft;
double L = (lnqmax_fft-lnqmin_fft);
double dx = L/(double) N;
double xarr[2*N],kernel_scatter[N];
int i;
double sigma = ptsz->sigmaM_ym;
double sig2 = pow(sigma,2.);
double fac_scatter =1./sqrt(2.*_PI_*sig2);
for (i = 0; i < N; i++) {
  double x = lnqmin_fft+i*dx;
  xarr[i] = x;
  xarr[i+N] = lnqmax_fft+i*dx;
  double arg0 = (x)/sqrt(2.*sig2);
  kernel_scatter[i] = fac_scatter*exp(-arg0*arg0);
}
double qmin_fft_padded = ptsz->szcounts_qmin_fft_padded;
double qmax_fft_padded = ptsz->szcounts_qmax_fft_padded;

double L_q = (qmax_fft_padded-qmin_fft_padded);
double dq = L_q/(double) (N);
// double fac = 1./sqrt(2.*_PI_);
double qp[2*N],kernel_qobs[N];
for (i = 0; i < N; i++) {
  double x = qmin_fft_padded+i*dq;
  qp[i] = x;
  qp[i+N] = qmax_fft_padded+i*dq;
  double arg0 = x/sqrt(2.);
  double fac = 1./sqrt(2.*_PI_);
  kernel_qobs[i] = fac*exp(-arg0*arg0);
}

int index_zloop;
for (index_zloop = 0;
     index_zloop<ptsz->szcounts_fft_nz;
     index_zloop++){

double z = ptsz->szcounts_fft_z[index_zloop];
// printf("In zloop computing z = %.5e\n",z);
// at each z we parallelize the loop over patches

int index_patchesloop;

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

// number_of_threads = 1; // just pretend we have one thread
int id;
omp_lock_t lock;
// int npatches = 10;//ptsz->nskyfracs;

// double result_qmconv_all_patches;//[ptsz->nskyfracs][2*N];

#pragma omp parallel \
   shared(N,xarr,kernel_scatter,npatches,index_zloop,abort,pba,ptsz,ppm,pnl,lock)\
   private(id,index_patchesloop)\
   num_threads(number_of_threads)
	 {

#ifdef _OPENMP
	   tstart = omp_get_wtime();
#endif


#pragma omp for schedule (dynamic)
for (index_patchesloop=0;index_patchesloop<npatches;index_patchesloop++) // number of patches: ptsz->nskyfracs;
	     {
#pragma omp flush(abort)

id = omp_get_thread_num();

// printf("thread %d computing patch %d.\n",id,index_patchesloop);



// first we tabulate q as a function of m.
int idpatch = index_patchesloop;
int ntab = ptsz->ntab_dlnm_dlnq;
double lnq_tab[ntab];
double lnm_tab[ntab];

int itab;
double lnm_tab_mmin = log(ptsz->M1SZ); // add some padding
double lnm_tab_mmax = log(ptsz->M2SZ); // add some padding
double dlnm_tab = (lnm_tab_mmax-lnm_tab_mmin)/(ntab-1.);


// printf("saving lnq - lnm array for iz = %d\n",index_zloop);
  // FILE *fp;
  // char Filepath[_ARGUMENT_LENGTH_MAX_];
  // sprintf(Filepath,"%s%s%d%s",ptsz->root,"test_lnq_lnm_iz_",index_zloop,".txt");
  // fp=fopen(Filepath, "w");

for (itab = 0;itab<ntab;itab++){
    lnm_tab[itab] = lnm_tab_mmin+itab*dlnm_tab;
    double mtab = exp(lnm_tab[itab]);
    double ytab = get_y_at_m_and_z(mtab,z,ptsz,pba);
    double thetatab = get_theta_at_m_and_z(mtab,z,ptsz,pba);
    double sigtab = get_szcountsz_sigma_at_theta_in_patch(thetatab,idpatch,ptsz);
    // if (log(sigtab)>10)
    //   lnq_tab[itab] = -100.;
    // else

    lnq_tab[itab] = log(ytab/sigtab);

    if (isinf(lnq_tab[itab])||isnan(lnq_tab[itab])){
      printf("lnq_tab[itab] = %.5e ytab =%.5e sigtab = %.5e\n",lnq_tab[itab],ytab,sigtab);
      exit(0);
    }
    // if (itab % 300 == 0)
    //   printf("idpatch %d z = %.5e lnq = %.5e lnm = %.5e\n",idpatch,z,lnq_tab[itab],lnm_tab[itab]);


      // fprintf(fp,"%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",lnq_tab[itab],lnm_tab[itab],thetatab,sigtab,ytab);

    }
// fclose(fp);
// printf("lnq -  lnm array saved");

// exit(0);
// set up the convolutions

// setup the bounds for the convolutions
// double lnqmin_fft = ptsz->szcounts_lnqmin_fft;
// double lnqmax_fft = ptsz->szcounts_lnqmax_fft;
// double L = (lnqmax_fft-lnqmin_fft);
// double dx = L/(double) N;

// convolution arrays
// double complex* a = malloc(sizeof(complex double)*N);
// double complex* b = malloc(sizeof(complex double)*N);

fftw_complex* in;
fftw_complex* out;

fftw_complex* test;
fftw_complex* out_test; /* double [2] */

fftw_complex* product;
fftw_complex* product_out;


in = fftw_alloc_complex(2*N);
out = fftw_alloc_complex(2*N);

test = fftw_alloc_complex(2*N);
out_test = fftw_alloc_complex(2*N);

product = fftw_alloc_complex(2*N);
product_out = fftw_alloc_complex(2*N);
// fftw_free(a_tmp);
// fftw_free(b_tmp);

    // the probability distribution of qt|qm
    // double sigma = ptsz->sigmaM_ym;
    // double sig2 = pow(sigma,2.);
    // double fac =1./sqrt(2.*_PI_*sig2);


// if (index_zloop == 1){
// printf("saving dlnqdlnm etc arrays for iz = %d\n",index_zloop);
//   // FILE *fp;
//   // char Filepath[_ARGUMENT_LENGTH_MAX_];
//   sprintf(Filepath,"%s%s%d%s",ptsz->root,"test_lnq_dn_iz_",index_zloop,".txt");
//   fp=fopen(Filepath, "w");
// }



// fill the arrays
  int i;
  for (i = 0; i < N; i++) {
    double x = xarr[i];
    // xarr[i] = x;
    // xarr[i+N] = lnqmax_fft+i*dx;
    // double arg0 = (x)/sqrt(2.*sig2);
    in[i][0] = kernel_scatter[i];//fac*exp(-arg0*arg0);// /exp(x); // divided by "q" as it is a lognormal, but canceled since we integrate in logq
    in[i][1]= 0.;
    in[N+i][0]= 0.;
    in[N+i][1]= 0.;

    double lnm;
    // if (x<lnq_tab[0]){
    //   printf("x=%.5e lnqtab[0]= %.5e\n",x,lnq_tab[0]);
    // }
    // else if (x>lnq_tab[ntab-1]){
    //   printf("x=%.5e lnqtab[ntab-1]= %.5e\n",x,lnq_tab[ntab-1]);
    // }
    lnm  = pwl_value_1d(ntab,
                        lnq_tab,
                        lnm_tab,
                        x);
    double m = exp(lnm);
    double dNdlnm = get_dndlnM_at_z_and_M(z,m,ptsz);

    // if (i % 100 == 0)
    //  printf("thread %d z = %.5e lnq = %.5e m = %.5e dNdlnm = %.5e\n",id,z,x,m,dNdlnm);
if (ptsz->use_skyaveraged_noise){
    // printf("using sky averaged sigma.\n");
    dNdlnm *= get_volume_at_z(z,pba)*4.*M_PI*ptsz->fsky_from_skyfracs;
  }
else{
    dNdlnm *= get_volume_at_z(z,pba)*4.*M_PI*ptsz->skyfracs[index_patchesloop];
}

    // compute derivative dlnmdlnq
    double tol = ptsz->tol_dlnm_dlnq;
    double lnqp = x+tol;
    double lnqm = x-tol;
    double lnmp = pwl_value_1d(ntab,
                               lnq_tab,
                               lnm_tab,
                               lnqp);
    double lnmm = pwl_value_1d(ntab,
                              lnq_tab,
                              lnm_tab,
                              lnqm);
    double dlnmdlnq = (lnmp-lnmm)/2./tol;
    if (isnan(dlnmdlnq)||(isinf(dlnmdlnq))){
      // printf("dlnmdlnq = %.5e lnmp = %.5e lnmm %.5e\n",dlnmdlnq,lnmp,lnmm);
      // exit(0);
      dlnmdlnq  = 0.; // kill this case
    }

    // check at th

    // if (i % 100 == 0)
    //  printf("thread %d z = %.5e lnq = %.5e dlnmdlnq = %.5e dNdlnm = %.5e\n",id,z,x,dlnmdlnq,dNdlnm);


    double dNdlnq = dNdlnm*dlnmdlnq;

    test[i][0] = dNdlnq;
    test[i][1] = 0.;

// if(index_zloop == 1){
    // fprintf(fp,"%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",x,dlnmdlnq,dNdlnm,dNdlnq,in[i][0]);
// }

    test[N+i][0] = 0.;
    test[N+i][1] = 0.;



  }// end fill arrays to convolve

// if (index_zloop == 1){

  // fclose(fp);
  // printf("array saved");
// }
// exit(0);
  fftw_execute_dft(ptsz->forward_plan_counts_fft, (fftw_complex*) in, (fftw_complex*) out);
  fftw_execute_dft(ptsz->forward_plan_counts_fft, (fftw_complex*) test, (fftw_complex*) out_test);

for (i = 0; i < 2*N; i++){
   product[i][0]= (out[i][0]*out_test[i][0]-out[i][1]*out_test[i][1])/(double)(2*N);
   product[i][1]= (out[i][0]*out_test[i][1]+out[i][1]*out_test[i][0])/(double)(2*N);
    }


    fftw_execute_dft(ptsz->reverse_plan_counts_fft, (fftw_complex*) product, (fftw_complex*) product_out);


// printf("saving result of qm convolution iz = %d\n",index_zloop);


    double result_qmconv[2*N];
    for(i = 0; i < 2*N; i++){
        product_out[i][0] *= dx;
        result_qmconv[i] = product_out[i][0];
        if (result_qmconv[i]<0.)
          result_qmconv[i] = 0.;
      // if (ptsz->sz_verbose>=1)
      //   printf("fft nexpected thread %d i = %d r = %.5e\n",id,i,product_out[i][0]);
      }
  int shift = (int) (N/2);
  shiftArray(result_qmconv, 2*N, shift);


  // FILE *fp;
  // char Filepath[_ARGUMENT_LENGTH_MAX_];
  // sprintf(Filepath,"%s%s%d%s",ptsz->root,"test_qmconv_iz_",index_zloop,".txt");
  // fp=fopen(Filepath, "w");

  for(i = 0; i < 2*N; i++){
    ptsz->szcounts_fft_qmconv_all_patches[index_patchesloop][i] = result_qmconv[i];
    }

// fclose(fp);


fftw_free(in);
fftw_free(out);
fftw_free(test);
fftw_free(out_test);
fftw_free(product);
fftw_free(product_out);


} // end patches loop
#ifdef _OPENMP
tstop = omp_get_wtime();
if (ptsz->sz_verbose > 0)
printf("In %s: time spent in parallel region (loop over cluster counts final ffts) = %e s for thread %d\n",
         __func__,tstop-tstart,omp_get_thread_num());


#endif

} //end of parallel region

if (abort == _TRUE_) return _FAILURE_;


double result_qmconv_all[2*N];
int i;
for(i = 0; i < 2*N; i++){
  result_qmconv_all[i] = 0.;
for (index_patchesloop=0;index_patchesloop<npatches;index_patchesloop++) // number of patches: ptsz->nskyfracs;
	     {
         result_qmconv_all[i] += ptsz->szcounts_fft_qmconv_all_patches[index_patchesloop][i];
       }


       // if (i % 100 == 0)
       //  printf("z = %.5e i = %d convqm = %.5e\n",z,i,result_qmconv_all[i]);

     }

  // return _SUCCESS_;


// printf("first convolution done.\n");
//   return _SUCCESS_;
// now we do the second convolution
// double qmin_fft_padded = ptsz->szcounts_qmin_fft_padded;
// double qmax_fft_padded = ptsz->szcounts_qmax_fft_padded;
//
// double L_q = (qmax_fft_padded-qmin_fft_padded);
// double dq = L_q/(double) (N);
// double fac = 1./sqrt(2.*_PI_);
// double qp[2*N];

// convolution arrays
// double xarr[2*N];
// setup the bounds for the convolutions
// double lnqmin_fft = ptsz->szcounts_lnqmin_fft;
// double lnqmax_fft = ptsz->szcounts_lnqmax_fft;
// double L = (lnqmax_fft-lnqmin_fft);
// double dx = L/(double) N;


// fftw_complex in[2*N], out[2*N], in2[2*N],test[2*N],out_test[2*N]; /* double [2] */
// fftw_complex product[2*N],product_out[2*N];


fftw_complex* in;
fftw_complex* out;

fftw_complex* test;
fftw_complex* out_test; /* double [2] */

fftw_complex* product;
fftw_complex* product_out;

in = fftw_alloc_complex(2*N);
out = fftw_alloc_complex(2*N);

test = fftw_alloc_complex(2*N);
out_test = fftw_alloc_complex(2*N);

product = fftw_alloc_complex(2*N);
product_out = fftw_alloc_complex(2*N);
// printf("saving final conv arrays etc arrays for iz = %d\n",index_zloop);
//   // FILE *fp;
//   // char Filepath[_ARGUMENT_LENGTH_MAX_];
// FILE *fp;
// char Filepath[_ARGUMENT_LENGTH_MAX_];
// sprintf(Filepath,"%s%s%d%s",ptsz->root,"test_final_array_to_conv_iz_",index_zloop,".txt");
// fp=fopen(Filepath, "w");

// // restore the xarr array
// for (i = 0; i < N; i++) {
//   double xlnq = lnqmin_fft+i*dx;
//   xarr[i] = xlnq;
//   xarr[i+N] = lnqmax_fft+i*dx;
// }


for (i = 0; i < N; i++) {

 // double x = qmin_fft_padded+i*dq;
 // qp[i] = x;
 // // ptsz->szcounts_fft_qobs[i] = qp[i];
 //
 // qp[i+N] = qmax_fft_padded+i*dq;
 double x = qp[i];
 // double arg0 = x/sqrt(2.);
 in[i][0] = kernel_qobs[i];//fac*exp(-arg0*arg0);
 in[i][1]= 0.;
 in[N+i][0]= 0.;
 in[N+i][1]= 0.;

   double lnqp =  log(x);
   double conv1;

   if (x<=0.){
     conv1=0.;
   }
   else{
  lnqp =  log(x);
   if (lnqp<xarr[0])
    conv1 = 0.;
   else if (lnqp>xarr[2*N-1])
    conv1 = 0.;
   else
    conv1 =  pwl_value_1d(2*N,
                          xarr,
                          result_qmconv_all,
                          lnqp)/x;
   }
 test[i][0] = conv1;
 test[i][1] = 0.;
 test[N+i][0] = 0.;
 test[N+i][1] = 0.;

// fprintf(fp,"%.5e\t%.5e\t%.5e\n",qp[i],test[i][0],in[i][0]);


}

// fclose(fp);

fftw_execute_dft(ptsz->forward_plan_counts_fft, (fftw_complex*) in, (fftw_complex*) out);
fftw_execute_dft(ptsz->forward_plan_counts_fft, (fftw_complex*) test, (fftw_complex*) out_test);

    for (i = 0; i < 2*N; i++){
       product[i][0]= (out[i][0]*out_test[i][0]-out[i][1]*out_test[i][1])/(double)(2*N);
       product[i][1]= (out[i][0]*out_test[i][1]+out[i][1]*out_test[i][0])/(double)(2*N);
        }

fftw_execute_dft(ptsz->reverse_plan_counts_fft, (fftw_complex*) product, (fftw_complex*) product_out);

    // for(i = 0; i < 2*N; i++){
    //     product_out[i][0] *= dq;
    //   if (ptsz->sz_verbose>10)
    //     printf("fft thread %d i = %d r = %.5e\n",id,i,product_out[i][0]);
    //   }

  double result_qconv[2*N];
    for(i = 0; i < 2*N; i++) {
            result_qconv[i] = product_out[i][0]*dq;
        }
  int shift = (int) (N/2);
  shiftArray(result_qconv, 2*N, shift);



fftw_free(in);
fftw_free(out);
fftw_free(test);
fftw_free(out_test);
fftw_free(product);
fftw_free(product_out);

  // FILE *fp;
  // char Filepath[_ARGUMENT_LENGTH_MAX_];
  // sprintf(Filepath,"%s%s%d%s",ptsz->root,"test_allpatches_final_qconv_iz_",index_zloop,".txt");
  // fp=fopen(Filepath, "w");


for(i = 0; i < N; i++) {
  int index_zq = ptsz->szcounts_fft_index_zq_final[index_zloop][i];
  ptsz->szcounts_fft_dndzdq[index_zq] = result_qconv[i];


  // fprintf(fp,"%.5e\t%.5e\n",qp[i],result_qconv[i]);
  // if ((index_zloop%3==0)&&(i%100==0)){
  //   printf("z = %.5e q = %.5e dndzdq = %.5e %d\n",z,qp[i],ptsz->szcounts_fft_dndzdq[index_zq],index_zq);
  // }

  // if ((z>0.08) && (z<0.12) && (qp[i]>2.8) && (qp[i]<3.2)){
  //   printf("z = %.5e q = %.5e  conv = %.5e\n",z,qp[i],result_qconv[i]);
  // }

  }
  // fclose(fp);

if (ptsz->sz_verbose>=1)
  printf("second convolution done for z = %.5e.\n",z);

     }// end zloop


// destroy the fft plans.
// fftw_destroy_plan(forward_plan);
// fftw_destroy_plan(reverse_plan);


// test
// double ztest = 8.45000e-02;
// double qtest = 3.02734e+00;
// double rtest = get_szcounts_dndzdq_at_z_q(ztest,qtest,ptsz);
//
// printf("ztest = %.5e qtest = %.5e rtest = %.5e\n",ztest,qtest,rtest);
//

// now we compute the likelihoods over the catalogue.

      // if (ptsz->sz_verbose>=1){

        int icat = 0;
        int imissz = 0;
        double max_z = 0.;
        double min_z = 100.;
        int imissz_arr[ptsz->szcat_size];
        if (ptsz->sz_verbose>=1) printf("szcat: got %d lines.\n",ptsz->szcat_size);
        for (icat=0;icat<ptsz->szcat_size;icat++){
            if (ptsz->sz_verbose>=10) printf("szcat z = %.3e \t snr = %.3e\n",ptsz->szcat_z[icat],ptsz->szcat_snr[icat]);
        if (ptsz->szcat_z[icat]<= 0) {
          imissz_arr[imissz] = icat;
          imissz += 1;
        }
        else{
          double zc = ptsz->szcat_z[icat];
          if (zc<min_z) min_z = zc;
          if (zc>max_z) max_z = zc;
        }
        }
        if (ptsz->sz_verbose>=1) {
          printf("cat min z = %.5e cat max z = %.5e\n",min_z,max_z);
        }
        if (min_z<ptsz->szcounts_fft_z_min){
          printf("You need to decrease szcounts_fft_z_min, min z from cat is %.5e\n",min_z);
          printf("Your current min z is %.5e\n",ptsz->szcounts_fft_z_min);
        }
        if (max_z>ptsz->szcounts_fft_z_max){
          printf("You need to increase szcounts_fft_z_max, max z from cat is %.5e\n",max_z);
          printf("Your current max z is %.5e\n",ptsz->szcounts_fft_z_max);
        }
// randomly assign zs to missing:
       const gsl_rng_type * T;
       gsl_rng * rmissz;
       gsl_rng_env_setup();
       T = gsl_rng_default;
       rmissz = gsl_rng_alloc (T);
       int iimissz;

       for (iimissz = 0; iimissz < imissz; iimissz++)
         {
           double u = gsl_rng_uniform (rmissz);
           ptsz->szcat_z[imissz_arr[iimissz]] = min_z + (max_z-min_z)*u;
           if (ptsz->sz_verbose>=1){
           printf ("cluster %d now has z = %.5f\n",imissz_arr[iimissz], u);
            }
         }

       gsl_rng_free (rmissz);

      if (ptsz->sz_verbose>=1){
        printf("szcat: got %d objects with missing redshift.\n",imissz);
        printf("random z's have been assigned to them, see above.\n");
      }

      // }


   // class_alloc(ptsz->szrate,ptsz->szcat_size*sizeof(double *),ptsz->error_message);
   int irate;
   for (irate = 0; irate < ptsz->szcat_size ; irate++){
     // ptsz->szrate[irate] = 0.;
     double qcat = ptsz->szcat_snr[irate];
     double zcat = ptsz->szcat_z[irate];
     ptsz->szrate[irate] = get_szcounts_dndzdq_at_z_q(zcat,qcat,ptsz);
      if (ptsz->sz_verbose>=1){
        printf("rate of cluster[%d] = %.5e\n",irate,ptsz->szrate[irate]);
      }
   }

//    // now we just need to compute the total number of clusters above the threshold.
//    if (ptsz->sz_verbose>=1)
//     printf("computing nexpected above qcut = %.5e\n",ptsz->sn_cutoff);
//    // ptsz->sn_cutoff
//    int index_z;
//    int index_q;
//    double dz_ntot = (ptsz->szcounts_fft_z_max-ptsz->szcounts_fft_z_min)/(ptsz->szcounts_fft_nz-1.);
//    int nq = 200;
//    double dq_ntot = (ptsz->szcounts_qmax_fft_padded-ptsz->sn_cutoff)/(nq-1.);
//    double ntot = 0.;
//    double arr_dndzdq[ptsz->szcounts_fft_nz][nq];
//    for (index_z = 0; index_z<ptsz->szcounts_fft_nz; index_z++){
//      for (index_q = 0; index_q<nq; index_q++){
//         double zp = ptsz->szcounts_fft_z_min+index_z*dz_ntot;
//         double qp = ptsz->sn_cutoff+index_q*dq_ntot;
//         double dn = get_szcounts_dndzdq_at_z_q(zp,qp,ptsz);
//         arr_dndzdq[index_z][index_q] = dn;//ntot += 0.5*()*dz*dq;
//         // if (index_z == 0 || index_z == ptsz->szcounts_fft_nz-1 || index_q == 0 || index_q == nq-1) {
//         //                 ntot += dn*dz_ntot*dq_ntot;
//         //             } else {
//                         ntot += dn*dz_ntot*dq_ntot;
//                     // }
//
//      }
//    }
//
//
//         // Calculating the integral value
//         // wrt y at each point for x
//         double ax[ptsz->szcounts_fft_nz];
//         for (index_z = 0; index_z < ptsz->szcounts_fft_nz; ++index_z)
//         {
//             ax[index_z] = 0;
//             for (index_q = 0; index_q < nq; ++index_q)
//             {
//                 // if (index_q == 0 || index_q == nq - 1)
//                     ax[index_z] += arr_dndzdq[index_z][index_q];
//                 // else if (index_q % 2 == 0)
//                 //     ax[index_z] += 2 * arr_dndzdq[index_z][index_q];
//                 // else
//                 //     ax[index_z] += 4 * arr_dndzdq[index_z][index_q];
//             }
//             ax[index_z] *= (dq_ntot);// / 3);
//         }
//
//         ntot = 0.;
//
//         // Calculating the final integral value
//         // using the integral obtained in the above step
//         for (index_z = 0; index_z < ptsz->szcounts_fft_nz; ++index_z)
//         {
//             // if (index_z == 0 || index_z == ptsz->szcounts_fft_nz - 1)
//                 ntot += ax[index_z];
//             // else if (i % 2 == 0)
//             //     ntot += 2 * ax[index_z];
//             // else
//             //     ntot += 4 * ax[index_z];
//         }
//         ntot *= (dz_ntot);// / 3);
//
//
//   ptsz->szcounts_ntot = ntot;
//
//
//
// // if (ptsz->has_sz_rates){
//   ptsz->szunbinned_loglike = - ntot;
//   int index_rate;
//   for (index_rate=0;index_rate<ptsz->szcat_size;index_rate++){
//     ptsz->szunbinned_loglike += log(ptsz->szrate[index_rate]);
//   // printf("cluster id = %d\trate = %.4e \n",index_rate,ptsz->szrate[index_rate]);
//   }
//   if (ptsz->sz_verbose >= 1 ) printf("loglike unbinned cc = %.3e\n",ptsz->szunbinned_loglike);
//   if (ptsz->sz_verbose >= 1 ) printf("Ntot above qcut = %.3e\n",ptsz->szcounts_ntot);
// // }



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


  int npatches;
  if (ptsz->use_skyaveraged_noise)
      npatches = 1;
  else
      npatches = ptsz->nskyfracs;

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
        for (index_patches =0;index_patches<npatches;index_patches++){

          // double y1 = ptsz->ylims[index_patches][l1];
          // double y2 = ptsz->ylims[index_patches][l2];
          // double y = y1 + (y2-y1)/(th2-th1)*(thp-th1);

          double y_interp = get_szcountsz_sigma_at_theta_in_patch(thp,index_patches,ptsz);;

          // double y_interp = pwl_value_1d(ptsz->nthetas,
          //                                ptsz->thetas,
          //                                ptsz->ylims[index_patches],
          //                                thp);
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

if (ptsz->use_skyaveraged_noise){
          completeness_2d[index_m][index_z] += c2*ptsz->fsky_from_skyfracs;///fsky;

        }
else{
     completeness_2d[index_m][index_z] += c2*ptsz->skyfracs[index_patches];///fsky;
}
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
        for (index_patches=0;index_patches<npatches;index_patches++){
          //fsky += ptsz->skyfracs[index_patches];

          double y1;
          if (ptsz->use_skyaveraged_noise == 0){
            y1 = ptsz->ylims[index_patches][index2];
          }
          else{
            y1 =  ptsz->sky_averaged_ylims[index2];
          }
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
          if (ptsz->use_skyaveraged_noise == 0){
          erfs[index1][index2]=erfs[index1][index2]+c2*ptsz->skyfracs[index_patches];
        }
        else{
              erfs[index1][index2]=erfs[index1][index2]+c2*ptsz->fsky_from_skyfracs;
        }
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
          if (ptsz->sz_verbose>3) printf("int_comp larger than fsky.\n");
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


     for (index_m=0;index_m<pcsz->nsteps_m;index_m++){
       free(completeness_2d[index_m]);
     }


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
// read_sz_catalog(ptsz);
  //pcsz->nzSZ = 20.; //cosmomc settings
  //pcsz->size_logM = 105; //cosmomc settings

  //pcsz->rho_m_at_z = ptsz->Omega_m_0*ptsz->Rho_crit_0*pow((1.+pcsz->redshift_for_dndm),3);

// printf("allocating memory for szcounts ffts.\n");
// exit(0);
  if (ptsz->has_sz_counts_fft){
if (ptsz->sz_verbose>1)
printf("allocating memory for szcounts ffts.\n");
ptsz->szcounts_fft_nqobs = ptsz->N_samp_fftw;
    class_alloc(ptsz->szcounts_fft_qobs,sizeof(double)*ptsz->szcounts_fft_nqobs,ptsz->error_message);
    // checked its well freed.

    class_alloc(ptsz->szcounts_fft_z,sizeof(double)*ptsz->szcounts_fft_nz,ptsz->error_message);
    // checked it's well freed.

    // class_alloc(ptsz->szcounts_fft_sigmayobs,sizeof(double)*ptsz->szcounts_fft_nsigmayobs,ptsz->error_message);
    // class_alloc(ptsz->szcounts_fft_index_zsig,sizeof(int *)*ptsz->szcounts_fft_nz,ptsz->error_message);
    // class_alloc(ptsz->szcounts_fft_index_zq,sizeof(int *)*ptsz->szcounts_fft_nz,ptsz->error_message);
    class_alloc(ptsz->szcounts_fft_index_zq_final,sizeof(int *)*ptsz->szcounts_fft_nz,ptsz->error_message);
    // checked.

    // class_alloc(ptsz->szcounts_fft_rates_at_z_sigy_qobs,sizeof(double *)*ptsz->szcounts_fft_nqobs,ptsz->error_message);
    class_alloc(ptsz->szcounts_fft_qmconv_all_patches,sizeof(double *)*ptsz->nskyfracs,ptsz->error_message);
    //checked

ptsz->szcounts_lnqmax_fft = log(ptsz->szcounts_qmax_fft_padded);
ptsz->szcounts_qmin_fft_padded = -ptsz->szcounts_qmax_fft_padded;
ptsz->szcounts_lnqmin_fft = -ptsz->szcounts_lnqmax_fft;

class_alloc(ptsz->array_y_to_m_redshift,sizeof(double *)*ptsz->n_z_y_to_m,ptsz->error_message);
class_alloc(ptsz->array_y_to_m_y,sizeof(double *)*ptsz->n_y_y_to_m,ptsz->error_message);
class_alloc(ptsz->array_y_to_m_at_z_y,sizeof(double *)*ptsz->n_z_y_to_m*ptsz->n_y_y_to_m,ptsz->error_message);
// all three checked


int i ;
double L_q = (ptsz->szcounts_qmax_fft_padded-ptsz->szcounts_qmin_fft_padded);
double dq = L_q/(double) (ptsz->szcounts_fft_nqobs);
for (i = 0; i < ptsz->szcounts_fft_nqobs; i++) {
   ptsz->szcounts_fft_qobs[i] = ptsz->szcounts_qmin_fft_padded+i*dq;
}

    int ipatches;
    for (ipatches = 0;ipatches<ptsz->nskyfracs;ipatches++){
      class_alloc(ptsz->szcounts_fft_qmconv_all_patches[ipatches],sizeof(double *)*2*ptsz->szcounts_fft_nqobs,ptsz->error_message);
    }
    // checked.

    int index_parallel_zp;
    int index_parallel_sigmayobsp;
    int index_qobs;

    // for (index_qobs = 0; index_qobs<ptsz->szcounts_fft_nqobs; index_qobs++){
    //   class_alloc(ptsz->szcounts_fft_rates_at_z_sigy_qobs[index_qobs],sizeof(double *)*ptsz->szcounts_fft_nz*ptsz->szcounts_fft_nsigmayobs,ptsz->error_message);
    // }



    // printf("allocating dndzdq\n");
    class_alloc(ptsz->szcounts_fft_dndzdq,
                sizeof(double *)*ptsz->szcounts_fft_nz*ptsz->szcounts_fft_nqobs,
                ptsz->error_message);
    // checked

    // printf("dndzdq allocated\n");
    // exit(0);
    // class_alloc(ptsz->szcounts_fft_nexpected_dndzdqgt,
    //             sizeof(double *)*ptsz->szcounts_fft_nz*ptsz->szcounts_fft_nexpected_qobs_n,
    //             ptsz->error_message);
    //
    // for (index_qobs = 0; index_qobs<ptsz->szcounts_fft_nexpected_qobs_n; index_qobs++){
    //   ptsz->szcounts_fft_nexpected_qobs[index_qobs] = ptsz->szcounts_fft_nexpected_qobs_min + index_qobs*(ptsz->szcounts_fft_nexpected_qobs_max-ptsz->szcounts_fft_nexpected_qobs_min)/(double) ptsz->szcounts_fft_nexpected_qobs_n;
    //
    //       if (ptsz->sz_verbose>=1){
    //         printf("%d qobs = %.5e\n",index_qobs,ptsz->szcounts_fft_nexpected_qobs[index_qobs]);
    //       }
    //   }

    int index_zq = 0;
    // int index_zqfinal = 0;
    for (index_parallel_zp=0; index_parallel_zp<ptsz->szcounts_fft_nz; index_parallel_zp++)
    {
      // class_alloc(ptsz->szcounts_fft_index_zq[index_parallel_zp],sizeof(int)*ptsz->szcounts_fft_nexpected_qobs_n,ptsz->error_message);
      class_alloc(ptsz->szcounts_fft_index_zq_final[index_parallel_zp],sizeof(int)*ptsz->szcounts_fft_nqobs,ptsz->error_message);
      //checked


      // for (index_qobs = 0; index_qobs<ptsz->szcounts_fft_nexpected_qobs_n; index_qobs++){
      //     ptsz->szcounts_fft_index_zq[index_parallel_zp][index_qobs] = index_zq;
      //     index_zq += 1;
      // }

      // for (index_qobs = 0; index_qobs<ptsz->szcounts_fft_nqobs; index_qobs++){
      //     ptsz->szcounts_fft_index_zq_final[index_parallel_zp][index_qobs] = index_zqfinal;
      //     index_zqfinal += 1;
      // }
    }
      // exit(0);
    int index_zqfinal = 0;
    for (index_parallel_zp=0; index_parallel_zp<ptsz->szcounts_fft_nz; index_parallel_zp++){
      for (index_qobs = 0; index_qobs<ptsz->szcounts_fft_nqobs; index_qobs++){

          ptsz->szcounts_fft_index_zq_final[index_parallel_zp][index_qobs] = index_zqfinal;
          index_zqfinal += 1;
      }
    }

    for (index_parallel_zp=0; index_parallel_zp<ptsz->szcounts_fft_nz; index_parallel_zp++)
    {
      ptsz->szcounts_fft_z[index_parallel_zp] = ptsz->szcounts_fft_z_min + index_parallel_zp*(ptsz->szcounts_fft_z_max-ptsz->szcounts_fft_z_min)/(double) (ptsz->szcounts_fft_nz-1.);

      // class_alloc(ptsz->szcounts_fft_index_zsig[index_parallel_zp],sizeof(int)*ptsz->szcounts_fft_nsigmayobs,ptsz->error_message);

    }
    // for (index_parallel_sigmayobsp=0; index_parallel_sigmayobsp<ptsz->szcounts_fft_nsigmayobs; index_parallel_sigmayobsp++)
    // {
    //   ptsz->szcounts_fft_sigmayobs[index_parallel_sigmayobsp] = log(ptsz->szcounts_fft_sigmayobs_min) + index_parallel_sigmayobsp*(log(ptsz->szcounts_fft_sigmayobs_max)-log(ptsz->szcounts_fft_sigmayobs_min))/ptsz->szcounts_fft_nsigmayobs;
    // }
    // int index_zsig = 0;
    // for (index_parallel_zp=0; index_parallel_zp<ptsz->szcounts_fft_nz; index_parallel_zp++){
    //   for (index_parallel_sigmayobsp=0; index_parallel_sigmayobsp<ptsz->szcounts_fft_nsigmayobs; index_parallel_sigmayobsp++){
    //     ptsz->szcounts_fft_index_zsig[index_parallel_zp][index_parallel_sigmayobsp] =index_zsig;
    //     index_zsig+=1;
    //   }
    // }
if (ptsz->sz_verbose>1)
printf("allocated memory for szcounts ffts.\n");
  }

  class_alloc(pcsz->redshift,sizeof(double)*pcsz->nzSZ,pcsz->error_message);
  // checked



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
