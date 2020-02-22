/** @file szpowerspectrum.h Documented includes for sz module. */

#ifndef __SZ__
#define __SZ__

#include "common.h"
#include "lensing.h"


#define _tSZ_power_spectrum_ ((ptsz->has_sz_ps == _TRUE_) && (index_md == ptsz->index_md_sz_ps))
#define _tSZ_trispectrum_ ((ptsz->has_sz_trispec == _TRUE_))
#define _tSZ_2halo_ ((ptsz->has_sz_2halo == _TRUE_))
#define _tSZ_te_y_y_ ((ptsz->has_sz_te_y_y == _TRUE_))
#define _cov_N_Cl_ ((ptsz->has_sz_cov_N_Cl == _TRUE_))
#define _mean_y_ ((ptsz->has_mean_y == _TRUE_) && (index_md == ptsz->index_md_mean_y))
#define _hmf_ ((ptsz->has_hmf == _TRUE_) && (index_md == ptsz->index_md_hmf))


struct tszspectrum {


  double hmf_int;
  double y_monopole;
  double * cl_sz;
  double * cl_te_y_y;
  double ** tllprime_sz;
  double * cl_sz_2h;
  double ** cov_N_cl;
  double ** r_N_cl;
  double ** r_cl_clp;
  double ** trispectrum_ref;
  double * cov_cl_cl;

  int index_md;

  int create_ref_trispectrum_for_cobaya;

  int has_sz_2halo;
  int has_sz_te_y_y;
  int has_sz_cov_N_Cl;

  int has_hmf;
  int index_md_hmf;
  int index_integrand_id_hmf;


  int has_mean_y;
  int index_md_mean_y;
  int index_integrand_id_mean_y;

  int has_sz_ps;
  int index_md_sz_ps;
  int index_integrand_id_sz_ps_first;
  int index_integrand_id_sz_ps_last;


  int has_sz_trispec;
  int index_md_sz_trispec;

  int index_integrand_id_sz_trispec_first;
  int index_integrand_id_sz_trispec_last;


  int number_of_integrands;
  int index_integrand;
  int index_integrand_te_y_y;
  int index_integrand_2halo_term;

  int index_integrand_2h_first; //for trispectrum
  int index_integrand_2h_last;  //for trispectrum

  int index_integrand_cov_N_cl_first;
  int index_integrand_cov_N_cl_last;


  int index_integrand_N_for_cov_N_cl_first;
  int index_integrand_N_for_cov_N_cl_last;


  int index_integrand_id;

  int number_of_integrals_per_thread;

  int index_integrands_first;
  int index_integrands_last;





  //double  pk;

  FileName root; /**< root for all file names */
  FileName path_to_class; /**< root for all file names */


 /* vector of all SZ quantities function of redshift*/

  int  tsz_size;

  int  index_flag_cov_N_cl;
  int  index_Rho_crit;
  int  index_Delta_c;
  int  index_rVIR;
  int  index_cVIR;
  int  index_m500;
  int  index_r500;
  int  index_l500;
  int  index_ls;
  int  index_rs;
  int  index_m200;
  int  index_Rh;
  int  index_mf;
  int  index_dlognudlogRh;
  int  index_lognu;
  int  index_dlogSigma2dlogRh;
  int  index_dndlogRh;
  int  index_logSigma2;
  int  index_z;
  int  index_m200c;
  int  index_l200c;
  int  index_r200c;
  int  index_multipole;
  int  index_multipole_prime;
  int  index_multipole_for_pressure_profile;
  int  index_pressure_profile;
  int  index_completeness;
  int  index_te_of_m;
  int  index_volume;
  int  index_pk_for_halo_bias;
  int  index_dlnMdeltadlnM;

  int  index_mean_y;
  int  index_hmf;


  int  index_halo_bias;
  int  index_k_value_for_halo_bias;

  int index_integral;
  int index_integral_te_y_y;
  int index_integral_2halo_term;

  int index_integral_2h_first;
  int index_integral_2h_last;

  int index_integral_cov_N_cl_first;
  int index_integral_cov_N_cl_last;

  int index_integral_N_for_cov_N_cl_first;
  int index_integral_N_for_cov_N_cl_last;


  int  index_integrals_over_m_first;
  int  index_integrals_over_m_last;

  int  index_integrals_over_z_first;
  int  index_integrals_over_z_last;



  int  index_integral_over_m;
  int  index_integral_te_y_y_over_m;
  int  index_integral_2halo_term_over_m;
  int  index_integral_2h_first_over_m;
  int  index_integral_2h_last_over_m;
  int  index_integral_cov_N_cl_first_over_m;
  int  index_integral_cov_N_cl_last_over_m;
  int  index_integral_N_for_cov_N_cl_first_over_m;
  int  index_integral_N_for_cov_N_cl_last_over_m;




  //mass bins for covariance between cluster counts and power spectrum
  int nbins_M;
  double * M_bins;
  double dlogM;




  //units for tSZ spectrum
  double exponent_unit;

  //completeness
  double theta_bin_min;
  double theta_bin_max;
  int nthetas;
  double * thetas;

  double *skyfracs;
  int nskyfracs;


  double ** ylims;


  int experiment;
  //SO completeness
  double * SO_thetas;
  double * SO_Qfit;
  int  SO_Q_size;

  double * SO_RMS;
  double * SO_skyfrac;
  int  SO_RMS_size;

  //INPUT PARAMETERS
  int nlSZ;
  int n_ell_independent_integrals;



  /*Redshift limits for the integration*/
  double z1SZ;
  double z2SZ;

  /*Array size*/
  int n_arraySZ;//number of z in the interpolation
  int n_arraySZ_for_integral;//number of z in the integration

  //mass limits: h^-1 Msun
  double M1SZ;
  double M2SZ;

  double delta_alpha;
  double alpha_p;

  double alpha_b;
  double Ap;
  int mass_dependent_bias;

  //Planck pressure profile
  double P0GNFW;
  double c500;
  double gammaGNFW;
  double alphaGNFW;
  double betaGNFW;

  double ln_x_size_for_pp;
  double * ln_x_for_pp;


  //Battaglia pressure profile
  double P0_B12;
  double xc_B12;
  double beta_B12;

  double alpha_m_P0_B12;
  double alpha_m_xc_B12;
  double alpha_m_beta_B12;

  double alpha_z_P0_B12;
  double alpha_z_xc_B12;
  double alpha_z_beta_B12;


  /*Pressure profile is considered between x_in and x_out*/
  double x_inSZ;
  double x_outSZ;

  double HSEbias;

  /*For the computation of sigma2*/
  int  ndimSZ;
  double logR1SZ; // 0.0034Mpc/h, 1.8e4  solar mass
  double logR2SZ; // 54.9Mpc/h, 7.5e16 solar mass
  double delta_cSZ;

  /*Multplicity function Tinker 2010*/

  double alphaSZ;
  double beta0SZ;
  double gamma0SZ;

  double phi0SZ;
  double eta0SZ;


  /*Multplicity function Bocquet 2015*/

  double Ap0;
  double a0;
  double b0;
  double c0;

  int MF;
  //1:Tinker 2010 (T10)
  //2:Bocquet 2015 (B15)
  //3:Jenkins 2001 (J01)
  //4:Tinker 2008 (T08)
  //5:Tinker 2008 interpolated @ M500 (T08@M500)

  //Precision Parameters For qromb_sz_integrand
  int K;
  double EPS;
  double JMAX;


  //Precision Parameters For qromb_sz_sigma
  int K_sigma;
  double EPS_sigma;
  double JMAX_sigma;


  //Foreground parameters
  double A_cib, A_rs, A_ir, A_cn;

  //Cl spectrum
  double * ell;
  double * ell_plc;
  double * ell_plc_low;
  double * ell_mock;
  double * ell_trispectrum;
  double * x_gauss;
  double * w_gauss;



  double dlogell;
  double ell_min_mock;
  double ell_max_mock;




  double Tcmb_gNU;

  double Rho_crit_0;
  double D_z1SZ;
  double Omega_m_0;
  double Omega_ncdm_0;

  double Sigma8OmegaM_SZ;
  double sigma8_Pcb;

  short has_tszspectrum;  //do we need tSZ spectrum? */
  short sz_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */
  short write_sz;  //do we need SZ quatitiies vs redshift? */

  short has_completeness_for_ps_SZ;
  short which_ps_sz;
  double H0_in_class_units;
  int ell_sz;
  // Figure 7 of KS02 -> KS02
  // Planck 2015 effective multipoles -> P15
  // SZFASTDKS -> DKS

  int concentration_parameter;
  //Duffy 2008: D08
  //Seljak 2000: S00

  int pressure_profile;
  //Planck 2013 (P13)
  //Arnaud et al 2010 (A10)
  //Custom. GNFW

  int HMF_prescription_NCDM;
  int effective_temperature;
  int temperature_mass_relation;
  int mean_y;

  double * PP_lnx;
  double * PP_lnI;
  double * PP_d2lnI;

  int PP_lnx_size;


  double * CM_redshift;
  double * CM_logM;

  int CM_redshift_size;
  int CM_logM_size;
  double * CM_logC;

  double * array_redshift;
  double * array_radius;
  double * array_sigma_at_z_and_R;
  double * array_dsigma2dR_at_z_and_R;



  ErrorMsg error_message; /**< zone for writing error messages */


};

/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int szpowerspectrum_init(struct background * pba,
                           struct nonlinear * pnl,
                           struct primordial * ppm,
                           struct tszspectrum * ptsz);


  int szpowerspectrum_free(struct tszspectrum *ptsz);


  int compute_sz(struct background * pba,
                 struct nonlinear * pnl,
                 struct primordial * ppm,
                 struct tszspectrum * ptsz,
                 double * pvecback,
                 double * Pvectsz);



  //This evaluates the integrand which will be integrated
  //over M and then over z
  int integrand_at_m_and_z(double logM ,
                           double * pvecback,
                           double * pvectsz,
                           struct background * pba,
                           struct primordial * ppm,
                           struct nonlinear * pnl,
                           struct tszspectrum * ptsz);



  int evaluate_HMF(double logM ,
                   double * pvecback,
                   double * pvectsz,
                   struct background * pba,
                   struct nonlinear * pnl,
                   struct tszspectrum * ptsz);

  int evaluate_completeness(double * pvecback,
                            double * pvectsz,
                            struct background * pba,
                            struct tszspectrum * ptsz);

  int evaluate_pressure_profile(double * pvecback,
                                double * pvectsz,
                                struct background * pba,
                                struct tszspectrum * ptsz);

  int write_output_to_files_ell_indep_ints(struct nonlinear * pnl,
                                           struct background * pba,
                                           struct tszspectrum * ptsz);

  int write_output_to_files_cl(struct nonlinear * pnl,
                               struct background * pba,
                               struct tszspectrum * ptsz);


  int show_preamble_messages(struct background * pba,
                             struct nonlinear * pnl,
                             struct primordial * ppm,
                             struct tszspectrum * ptsz);

  int show_results(struct background * pba,
                   struct nonlinear * pnl,
                   struct primordial * ppm,
                   struct tszspectrum * ptsz);

  int select_multipole_array(struct tszspectrum * ptsz);

  int evaluate_halo_bias(double * pvecback,
                         double * pvectsz,
                         struct background * pba,
                         struct primordial * ppm,
                         struct nonlinear * pnl,
                         struct tszspectrum * ptsz);

  int initialise_and_allocate_memory(struct tszspectrum * ptsz);


  int evaluate_temperature_mass_relation(double * pvecback,
                                         double * pvectsz,
                                         struct background * pba,
                                         struct tszspectrum * ptsz);

double evaluate_dlnMdeltadlnM(double logM,
                             double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct nonlinear * pnl,
                             struct tszspectrum * ptsz);

#ifdef __cplusplus
}
#endif

#endif
