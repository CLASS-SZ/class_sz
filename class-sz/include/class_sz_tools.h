
#ifndef __SZ_TOOLS__
#define __SZ_TOOLS__

#include "common.h"
#include "r8lib.h"
#include <time.h>
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf_bessel.h"

#define _M_sun_  1.98855e30 /*solar mass in kg*/



/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif


double get_delta_mean_from_delta_crit_at_z(double delta_crit,
                                           double z,
                                           struct class_sz_structure * pclass_sz);

double get_Omega_m_nonu_at_z(double z,
                             struct class_sz_structure * pclass_sz);


double get_z_of_chi(
                    double chi,
                    struct class_sz_structure * pclass_sz,
                    struct background * pba);


int y_to_m (
            double xout,
            double * mRES,
            double z,
            double m,
            // double rd,
            struct class_sz_structure * pclass_sz,
            struct background * pba,
            struct nonlinear * pnl,
            struct primordial * ppm
          );

int solve_y_to_m(
                    double * result,
                    double z,
                    double m,
                    struct class_sz_structure * pclass_sz,
                    struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm
                  );

int zbrent_y_to_m(
              double x1,
              double x2,
              double tol,
              double fa,
              double fb,
              double * knl,
              double z,
              double m,
              // double rd,
              struct class_sz_structure * pclass_sz,
              struct background * pba,
              struct nonlinear * pnl,
              struct primordial * ppm
            );


int chi_to_z (
            double xout,
            double * zRES,
            double chi,
            struct class_sz_structure * pclass_sz,
            struct background * pba
          );

int solve_chi_to_z(
                    double * result,
                    double chi,
                    struct class_sz_structure * pclass_sz,
                    struct background * pba
                  );



int zbrent_chi_to_z(
              double x1,
              double x2,
              double tol,
              double fa,
              double fb,
              double * knl,
              double chi,
              struct class_sz_structure * pclass_sz,
              struct background * pba
            );          


int m_to_xout (
            double xout,
            double * mRES,
            double z,
            double m,
            double rd,
            struct class_sz_structure * pclass_sz,
            struct background * pba,
            struct nonlinear * pnl,
            struct primordial * ppm
          );

int solve_m_to_xout(
                    double * result,
                    double z,
                    double m,
                    struct class_sz_structure * pclass_sz,
                    struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm
                  );

int zbrent_m_to_xout(
              double x1,
              double x2,
              double tol,
              double fa,
              double fb,
              double * knl,
              double z,
              double m,
              double rd,
              struct class_sz_structure * pclass_sz,
              struct background * pba,
              struct nonlinear * pnl,
              struct primordial * ppm
            );

int tabulate_m_to_xout(struct background * pba,
                       struct nonlinear * pnl,
                       struct primordial * ppm,
                       struct class_sz_structure * pclass_sz);

double get_dlnq_dlnm(double lnq,
                     double z,
                     int idpatch,
                     struct class_sz_structure * pclass_sz,
                     struct background * pba);


double get_theta_at_y_and_z(double y,
                            double z,
                            struct class_sz_structure * pclass_sz,
                            struct background * pba);
double get_szcountsz_sigma_at_theta_in_patch(double theta,int id_patch,struct class_sz_structure *pclass_sz);
int tabulate_y_to_m(struct background * pba,
                       struct nonlinear * pnl,
                       struct primordial * ppm,
                       struct class_sz_structure * pclass_sz);
double get_y_to_m_at_z_and_y(double z_asked, double y_asked, struct class_sz_structure * pclass_sz);
double get_dlnm_dlny(double lny,double z,struct class_sz_structure * pclass_sz);
double get_dNdlny_at_z_and_y(double z_asked, double y_asked, struct background * pba, struct class_sz_structure * pclass_sz);
double integrate_over_m_at_z(double * pvecback,
                           double * pvectsz,
                           struct background * pba,
                           struct nonlinear * pnl,
                           struct primordial * ppm,
                           struct perturbs * ppt,
                           struct class_sz_structure * pclass_sz);


  int zbrent_sz_delta_to_delta_prime_nfw(double x1,
                                          double x2,
                                          double tol,
                                          double cvir,
                                          double cvir_prime,
                                          double delta,
                                          double fa,
                                          double fb,
                                          double * delta_prime,
                                        struct class_sz_structure * pclass_sz);


  int zbrent_sz(double x1,
                double x2,
                double tol,
                double mVIR,
                double rvir,
                double c,
                double delrho,
                double fa,
                double fb,
                double * logMDEL,
                struct class_sz_structure * pclass_sz);

  int zbrent_D_to_V_sz(double x1,
                        double x2,
                        double tol,
                        double mDEL,
                        double delrho,
                        double fa,
                        double fb,
                        double z,
                        double delc,
                        double rhoc,
                        double * logMVIR,
                        struct class_sz_structure * pclass_sz,
                        struct background * pba);






  int mVIR_to_mDEL(double mVIR ,
                double rvir ,
                double c ,
                double delrho,
                double * result,
                struct class_sz_structure * pclass_sz);


 int mDEL_to_mVIR(
               double mDEL ,
               double delrho,
               double delc,
               double rhoc,
               double z,
               double * result,
               struct class_sz_structure * pclass_sz,
               struct background * pba
             );

double m_nfw(double x);

double delta_to_delta_prime_nfw(
  double delta,
  double cvir,
  double cvir_prime,
  struct class_sz_structure * pclass_sz
);

 int mDEL_to_mDELprime(
               double mDEL ,
               double delrho,
               double delrho_prime,
               double delc,
               double rhoc,
               double z,
               double * mDELprime,
               struct class_sz_structure * pclass_sz,
               struct background * pba
             );

 int mDtomV (
             double logMVIR ,
             double mD,
             double rvir,
             double c,
             double delrho,
             double * mRES,
             struct class_sz_structure * pclass_sz
           );

 int mVtomD (
             double logMD ,
             double mVIR,
             double rvir,
             double c,
             double delrho,
             double * mRES,
             struct class_sz_structure * pclass_sz
           );

 int dtod_prime_nfw( double delta_prime,
                     double delta,
                     double cvir,
                     double cvir_prime,
                     double * dRES
                   );

 double get_planck_sigma_at_theta500(double theta500, struct class_sz_structure * pclass_sz);

 double get_m_to_xout_at_z_and_m(double z_asked, double m_asked, struct class_sz_structure * pclass_sz);
 double evaluate_cvir_of_mvir(double mvir,
                             double z,
                             struct class_sz_structure * pclass_sz,
                             struct background * pba);

 double evaluate_rvir_of_mvir(double mvir,
                             double delc,
                             double rhoc,
                             struct class_sz_structure * pclass_sz);


  int CvirMvirKLYPIN(double * result,
                     double logM ,
                     double z,
                     struct class_sz_structure * pclass_sz);

  int read_Zhao_CM_init(struct class_sz_structure * pclass_sz);

  int load_rho_nfw_profile(struct class_sz_structure * pclass_sz);
  int load_T10_alpha_norm(struct class_sz_structure * pclass_sz);
  int load_normalized_dndz(struct class_sz_structure * pclass_sz);
  int load_normalized_dndz_ngal(struct class_sz_structure * pclass_sz);
  int load_normalized_source_dndz(struct class_sz_structure * pclass_sz);
  int load_normalized_fdndz(struct class_sz_structure * pclass_sz);
  int load_normalized_cosmos_dndz(struct class_sz_structure * pclass_sz);
  int load_unbinned_nl_yy(struct class_sz_structure * pclass_sz);
  int load_unbinned_nl_tt(struct class_sz_structure * pclass_sz);
  int load_ksz_filter(struct class_sz_structure * pclass_sz);
  int load_M_min_of_z(struct class_sz_structure * pclass_sz);
  double get_T10_alpha_at_z(double z_asked, struct class_sz_structure * pclass_sz);
  double get_knl_at_z(double z_asked, struct class_sz_structure * pclass_sz);
  double get_hmf_counter_term_nmin_at_z(double z_asked, struct class_sz_structure * pclass_sz);
  double get_hmf_counter_term_b1min_at_z(double z_asked, struct class_sz_structure * pclass_sz);
  double get_hmf_counter_term_b2min_at_z(double z_asked, struct class_sz_structure * pclass_sz);
  double get_m200m_to_m200c_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz);
  double get_m200c_to_mvir_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz);
  double get_m200m_to_mvir_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz);
  double get_m200c_to_m200m_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz);
  double get_m200m_to_m500c_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz);
  double get_m200c_to_m500c_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz);
  double get_m500c_to_m200c_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz);
  double get_m500c_to_m200m_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz);
  double get_nl_index_at_z_and_k(double z_asked, double k_asked, struct class_sz_structure * pclass_sz, struct nonlinear * pnl);
  double get_nl_index_at_z_and_k_no_wiggles(double z_asked, double k_asked, struct class_sz_structure * pclass_sz, struct nonlinear * pnl);

  double get_normalization_gas_density_profile(double z_asked, double m_asked, struct class_sz_structure * pclass_sz);
  double get_normalization_matter_density_profile(double z_asked, double m_asked, struct class_sz_structure * pclass_sz);


  double get_completeness_at_z_and_M(double z_asked, double m_asked, double * completeness_2d_to_1d, struct class_sz_structure * pclass_sz);
  double get_detection_proba_at_y_and_theta(double y_asked, double th_asked, double * erfs_2d_to_1d, struct class_sz_structure * pclass_sz);

  double get_y_at_m_and_z(double m, double z, struct class_sz_structure * pclass_sz, struct background * pba);
  double get_theta_at_m_and_z(double m, double z, struct class_sz_structure * pclass_sz, struct background * pba);
  double get_volume_at_z(double z, struct background * pba);

  int CvirMvirZHAO(double * result,
                   double logM ,
                   double z,
                   struct class_sz_structure * pclass_sz);

  int  C200M200SC14(double * result,
                    double logM ,
                    double z,
                    struct class_sz_structure * pclass_sz);

  int external_pressure_profile_init(struct precision * ppr, struct class_sz_structure * pclass_sz);

  int evaluate_redshift_int_lensmag(double * pvectsz,
                                    struct class_sz_structure * pclass_sz);

  int evaluate_redshift_int_nlensmag(double * pvectsz,
                                    struct class_sz_structure * pclass_sz);

  int tabulate_redshift_int_lensmag(struct class_sz_structure * pclass_sz,
                                    struct background * pba);
  int tabulate_redshift_int_nlensmag(struct class_sz_structure * pclass_sz,
                                    struct background * pba);
                                    
  double evaluate_redshift_int_gallens_sources(double z,
                                    struct class_sz_structure * pclass_sz);

  int tabulate_redshift_int_gallens_sources(struct class_sz_structure * pclass_sz,
                                    struct background * pba);
int tabulate_L_sat_at_nu_and_nu_prime(struct background * pba,
                                      struct class_sz_structure * pclass_sz);
int tabulate_L_sat_at_z_m_nu(struct background * pba,
                             struct class_sz_structure * pclass_sz);
int tabulate_m200c_to_m200m(struct background * pba,
                            struct class_sz_structure * pclass_sz);

int tabulate_m200m_to_m200c(struct background * pba,
                            struct class_sz_structure * pclass_sz);

int tabulate_m200c_to_mvir(struct background * pba,
                           struct class_sz_structure * pclass_sz);

int tabulate_m200m_to_mvir(struct background * pba,
                           struct class_sz_structure * pclass_sz);

int tabulate_m200m_to_m500c(struct background * pba,
                            struct class_sz_structure * pclass_sz);

int tabulate_m200c_to_m500c(struct background * pba,
                            struct class_sz_structure * pclass_sz);

int tabulate_m500c_to_m200c(struct background * pba,
                            struct class_sz_structure * pclass_sz);

int tabulate_m500c_to_m200m(struct background * pba,
                            struct class_sz_structure * pclass_sz);
  int redshift_int_lensmag(
                    struct class_sz_structure * pclass_sz,
                    struct background * pba,
                    double * pvectsz,
                    double * result
                  );

  int redshift_int_nlensmag(
                    int index_g,
                    struct class_sz_structure * pclass_sz,
                    struct background * pba,
                    double * pvectsz,
                    double * result
                  );

  int redshift_int_gallens_sources(
                    struct class_sz_structure * pclass_sz,
                    struct background * pba,
                    double * pvectsz,
                    double * result
                  );

  int spectra_sigma_prime(struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct class_sz_structure * pclass_sz,
                          double R,
                          double z,
                          double * sigma_prime);

  int spectra_sigma_for_tSZ(
                            struct background * pba,
                            struct primordial * ppm,
                            struct nonlinear * pnl,
                            struct class_sz_structure * pclass_sz,
                            double R,
                            double z,
                            double *sigma
                            );

  int spectra_sigma_ncdm(struct background * pba,
                         struct primordial * ppm,
                         struct nonlinear * pnl,
                         struct class_sz_structure * pclass_sz,
                         double R,
                         double z,
                         double *sigma);

  int spectra_sigma_ncdm_prime(struct background * pba,
                               struct primordial * ppm,
                               struct nonlinear * pnl,
                               struct class_sz_structure * pclass_sz,
                               double R,
                               double z,
                               double *sigma);

  int spectra_vrms2(
                   struct background * pba,
                   struct primordial * ppm,
                   struct nonlinear *pnl,
                   struct class_sz_structure * pclass_sz,
                   double z,
                   double * vrms2
                 );


 int tabulate_vrms2_from_pk(struct background * pba,
                            struct nonlinear * pnl,
                            struct primordial * ppm,
                            struct class_sz_structure * pclass_sz);

 int tabulate_knl(struct background * pba,
                  struct nonlinear * pnl,
                  struct primordial * ppm,
                  struct class_sz_structure * pclass_sz);

 int tabulate_nl_index(struct background * pba,
                       struct nonlinear * pnl,
                       struct primordial * ppm,
                       struct class_sz_structure * pclass_sz);


  int splint(double xa[],
             double ya[],
             double y2a[],
             int npoints,
             double x,
             double *y);


int solve_pkl_to_knl(
              double * result,
              double z,
              struct class_sz_structure * pclass_sz,
              struct background * pba,
              struct nonlinear * pnl,
              struct primordial * ppm
            );

int zbrent_pkl_to_knl(
              double x1,
              double x2,
              double tol,
              double fa,
              double fb,
              double * knl,
              double z,
              struct class_sz_structure * pclass_sz,
              struct background * pba,
              struct nonlinear * pnl,
              struct primordial * ppm
            );

int pkl_to_knl (
            double knl,
            double * mRES,
            double z,
            struct class_sz_structure * pclass_sz,
            struct background * pba,
            struct nonlinear * pnl,
            struct primordial * ppm
          );



  int p_gnfw (double * p_gnfw_x,
                double x ,
                double kl,

                double * pvectsz,
                struct background * pba,
                struct class_sz_structure * pclass_sz);

  int rho_nfw (double * rho_nfw_x,
                double x ,
                double * pvectsz,
                struct background * pba,
                struct class_sz_structure * pclass_sz);
  int rho_gnfw (double * rho_nfw_x,
                double x ,
                double * pvectsz,
                struct background * pba,
                struct class_sz_structure * pclass_sz);

  int MF_T10 (double * result,
              double * lognu ,
              double z ,
              struct class_sz_structure * pclass_sz);

  int MF_B15 (double * result,
              double * lognu ,
              double z ,
              struct class_sz_structure * pclass_sz);

  int MF_B15_M500c(double * result,
                   double * lognu ,
                   double z ,
                   struct class_sz_structure * pclass_sz);

double get_f_tinker10_at_nu_and_z(double nu, double z, struct class_sz_structure * pclass_sz);
double get_T10_alpha_at_z(double z,struct class_sz_structure * pclass_sz);
double get_f_tinker08_at_nu_and_z(double nu, double z, struct class_sz_structure * pclass_sz);

  int  MF_J01(double * result,
              double * lognu ,
              struct class_sz_structure * pclass_sz);

  int MF_T08(double * result,
             double * lognu ,
             double z ,
             struct class_sz_structure * pclass_sz);


  int MF_T08_m500(double * result,
                  double * lognu ,
                  double z ,
                  double delta_crit,
                  struct class_sz_structure * pclass_sz);

  //HMF Tinker et al 2008
  //@ M1600m
  int MF_T08_M1600m(double * result,
                    double * lognu ,
                    double z ,
                    struct class_sz_structure * pclass_sz);

  double erf_compl_ps(double y,
                      double sn,
                      double q);

double erf_compl_nicola(double y,
                        double sn,
                        double q,
                        double ymin,
                        double ymax,
                        double dof);

  //coding number count likelihood
  double erf_compl(double y,
                   double sn,
                   double q,
                   double dof);

  double d_erf_compl_dq(double y,
                        double sn,
                        double q);

  double next_z(double z_i, double dz, struct class_sz_structure * pclass_sz);

  int integrate_over_redshift(struct background * pba,
                              struct nonlinear * pnl,
                              struct primordial * ppm,
                              struct perturbs * ppt,
                              struct class_sz_structure * pclass_sz,
                              double * Pvecback,
                              double * Pvectsz);

  int read_Planck_noise_map(struct class_sz_structure * pclass_sz);

  int read_sz_catalog(struct class_sz_structure * pclass_sz);

  int tabulate_sigma_and_dsigma_from_pk(struct background * pba,
                                        struct nonlinear * pnl,
                                        struct primordial * ppm,
                                        struct class_sz_structure * pclass_sz);


  int spectra_sigma2_hsv(
                   struct background * pba,
                   struct primordial * ppm,
                   struct nonlinear *pnl,
                   struct class_sz_structure * pclass_sz,
                   double z,
                   double * sigma2_hsv
                 );

  int tabulate_sigma2_hsv_from_pk(struct background * pba,
                                  struct nonlinear * pnl,
                                  struct primordial * ppm,
                                  struct class_sz_structure * pclass_sz);




  int read_SO_Qfit(struct class_sz_structure * pclass_sz);

  int read_SO_noise(struct class_sz_structure * pclass_sz);


  int two_dim_ft_pressure_profile(double kl,
                                  struct class_sz_structure * pclass_sz,
                                  struct background * pba,
                                  double * pvectsz,
                                  double * result) ;

  int two_dim_ft_nfw_profile(struct class_sz_structure * pclass_sz,
                              struct background * pba,
                              double * pvectsz,
                              double * result) ;

double integrand_mass(double logM, void *p);
double integrand_sigma2_hsv(double lnk, void *p);
double nl_fitting_function(double lnk,void *p);
double Delta_c_of_Omega_m(double Omega_m);

double get_dndlnM_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz);


double get_L_sat_at_z_and_M_at_nu(double z_asked,
                                  double m_asked,
                                  int index_nu,
                                  struct background * pba,
                                  struct class_sz_structure * pclass_sz);

double  get_L_sat_at_z_M_nu(double z_asked, double m_asked, double nu_asked, struct class_sz_structure * pclass_sz);

// double get_L_sat_at_z_and_M_at_nu_prime(double z_asked,
//                                   double m_asked,
//                                   struct background * pba,
//                                   struct class_sz_structure * pclass_sz);

double get_dydz_at_z(double z_asked, struct class_sz_structure * pclass_sz);


int tabulate_dydz(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct class_sz_structure * pclass_sz);

double get_dcib0dz_at_z_and_nu(double z_asked, double nu_asked, struct class_sz_structure * pclass_sz);


int tabulate_dcib0dz(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct class_sz_structure * pclass_sz);

int tabulate_dndlnM(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct class_sz_structure * pclass_sz);

int bispectrum_condition(double ell_1, double ell_2, double ell_3);

int tabulate_mean_galaxy_bias(struct background * pba,
                              struct nonlinear * pnl,
                              struct primordial * ppm,
                              struct perturbs * ppt,
                              struct class_sz_structure * pclass_sz);


int tabulate_mean_galaxy_number_density(struct background * pba,
                                        struct nonlinear * pnl,
                                        struct primordial * ppm,
                                        struct class_sz_structure * pclass_sz);

int tabulate_mean_galaxy_number_density_ngal(struct background * pba,
                                             struct nonlinear * pnl,
                                             struct primordial * ppm,
                                             struct class_sz_structure * pclass_sz);



double get_HI_density_profile_at_k_M_z(double k_asked, double m_asked, double z_asked, struct class_sz_structure * pclass_sz);

double get_rho_2h_at_k_and_z(double k_asked,
                             double z_asked,
                             struct class_sz_structure * pclass_sz);

double get_rho_2h_at_r_and_m_and_z(double r_asked,
                                   double m_asked,
                                   double z_asked,
                                   struct class_sz_structure * pclass_sz,
                                   struct background * pba);

double get_gas_pressure_2h_at_k_and_z(double k_asked,
                                      double z_asked,
                                      struct class_sz_structure * pclass_sz);

double get_gas_pressure_2h_at_r_and_m_and_z(double r_asked,
                                            double m_asked,
                                            double z_asked,
                                            struct class_sz_structure * pclass_sz,
                                            struct background * pba);


double get_gas_density_profile_at_k_M_z(double l_asked,
                                    double m_asked,
                                    double z_asked,
                                    int normalize,
                                    struct class_sz_structure * pclass_sz);

double get_gas_pressure_profile_at_k_m_z(double l_asked,
                                    double m_asked,
                                    double z_asked,
                                    struct class_sz_structure * pclass_sz);

double get_gas_pressure_profile_at_l_m_z(double l_asked,
                                          double m_asked,
                                          double z_asked,
                                          struct class_sz_structure * pclass_sz);


double get_gas_pressure_profile_at_k(double k_asked,
                                    struct class_sz_structure * pclass_sz);


int tabulate_hmf_counter_terms_nmin(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct class_sz_structure * pclass_sz);

int tabulate_ng_bias_contribution_at_z_and_k(struct background * pba,
                                             struct perturbs * ppt,
                                             struct class_sz_structure * pclass_sz);


int tabulate_hmf_counter_terms_b1min(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct perturbs * ppt,
                                    struct class_sz_structure * pclass_sz);


int tabulate_hmf_counter_terms_b2min(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct class_sz_structure * pclass_sz);


int tabulate_psi_b1g(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct perturbs * ppt,
                    struct class_sz_structure * pclass_sz);

int tabulate_custom1_profile_fft(struct background * pba,
                                 struct class_sz_structure * pclass_sz);

int tabulate_gas_density_profile(struct background * pba,
                                 struct class_sz_structure * pclass_sz);

int tabulate_gas_density_profile_fft(struct background * pba,
                                     struct class_sz_structure * pclass_sz);

int tabulate_matter_nfw_with_power_law_profile_fft(struct background * pba,
                                                   struct class_sz_structure * pclass_sz);

int tabulate_gas_density_profile_2h(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct perturbs * ppt,
                                    struct class_sz_structure * pclass_sz);

int tabulate_gas_density_profile_2h_fft_at_z_and_r(struct background * pba,
                                                   struct nonlinear * pnl,
                                                   struct primordial * ppm,
                                                   struct class_sz_structure * pclass_sz);

int tabulate_gas_pressure_profile_2h(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct perturbs * ppt,
                                    struct class_sz_structure * pclass_sz);

int tabulate_gas_pressure_profile_2h_fft_at_z_and_r(struct background * pba,
                                                   struct nonlinear * pnl,
                                                   struct primordial * ppm,
                                                   struct class_sz_structure * pclass_sz);


double get_psi_b1g_at_k_and_z(double l_asked, double z_asked, struct class_sz_structure * pclass_sz);


int tabulate_psi_b2g(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct class_sz_structure * pclass_sz);

double get_psi_b2g_at_k_and_z(double l_asked, double z_asked, struct class_sz_structure * pclass_sz);


int tabulate_psi_b1kg(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct perturbs * ppt,
                    struct class_sz_structure * pclass_sz);


double get_psi_b1kg_at_k_and_z(double l_asked, double z_asked, struct class_sz_structure * pclass_sz);


int tabulate_psi_b2kg(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct class_sz_structure * pclass_sz);

double get_psi_b2kg_at_k_and_z(double l_asked, double z_asked, struct class_sz_structure * pclass_sz);



int tabulate_psi_b2t(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct class_sz_structure * pclass_sz);

double get_psi_b2t_at_k_and_z(double l_asked, double z_asked, struct class_sz_structure * pclass_sz);


int tabulate_n5k_F1(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct class_sz_structure * pclass_sz);


double get_dkappacmbdz_at_l_and_z(double l,
                                  double z,
                                  struct background * pba,
                                  struct primordial * ppm,
                                  struct nonlinear * pnl,
                                  struct class_sz_structure * pclass_sz);

double get_dkappacmbdz_pklin_at_l_and_z(double l,
                                  double z,
                                  struct background * pba,
                                  struct primordial * ppm,
                                  struct nonlinear * pnl,
                                  struct class_sz_structure * pclass_sz);


double get_n5k_pk_at_z_and_k(double z_asked, double k_asked, struct class_sz_structure * pclass_sz);
int load_n5k_pk_zk(struct class_sz_structure * pclass_sz);
int load_n5k_cl_K1(struct class_sz_structure * pclass_sz);
double get_n5k_cl_K1_at_chi(double chi,struct class_sz_structure * pclass_sz);


int load_cib_Snu(struct class_sz_structure * pclass_sz);
double get_cib_Snu_z_and_nu(double z_asked, double nu_asked, struct class_sz_structure * pclass_sz);

int load_n5k_z_of_chi(struct class_sz_structure * pclass_sz);
double get_n5k_z_of_chi(double chi,struct class_sz_structure * pclass_sz);

int tabulate_psi_b1t(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct perturbs * ppt,
                    struct class_sz_structure * pclass_sz);

double get_psi_b1t_at_k_and_z(double l_asked, double z_asked, struct class_sz_structure * pclass_sz);


int tabulate_psi_b1gt(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct perturbs * ppt,
                    struct class_sz_structure * pclass_sz);

double get_psi_b1gt_at_k1_k2_and_z(double l1_asked, double l2_asked, double z_asked, struct class_sz_structure * pclass_sz);

int tabulate_psi_b1kgg(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct perturbs * ppt,
                    struct class_sz_structure * pclass_sz);

double get_psi_b1kgg_at_k1_k2_and_z(double l1_asked, double l2_asked, double z_asked, struct class_sz_structure * pclass_sz);



int tabulate_psi_b1kgt(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct perturbs * ppt,
                    struct class_sz_structure * pclass_sz);

double get_psi_b1kgt_at_k1_k2_and_z(double l1_asked, double l2_asked, double z_asked, struct class_sz_structure * pclass_sz);

double get_lensing_noise_at_ell(double l,
                                struct class_sz_structure * pclass_sz);

#ifdef __cplusplus
}
#endif
#endif
