/** @file class_sz_clustercounts.h Documented includes for sz module. */

#ifndef __SZC__
#define __SZC__

#include "common.h"
#include "class_sz.h"
#include "class_sz_tools.h"
#include "r8lib.h"



//#define _mean_y_ ((pclass_sz->has_mean_y == _TRUE_) && (index_md == pclass_sz->index_md_mean_y))
//#define _hmf_ ((pclass_sz->has_hmf == _TRUE_) && (index_md == pclass_sz->index_md_hmf))


struct szcount {
  FileName root; /**< root for all file names */

  int nzSZ;
  //double * dndz;
  int has_sz_counts;
  double * redshift;
  double * z_center;
  double * steps_m;
  double * steps_z;
  double * logy;
  double ** dNdzdy_theoretical;
  int nsteps_m;
  int Nbins_z;
  int Nbins_y;
  int Ny;

  double redshift_for_dndm;

  double ystar;
  double alpha;
  double sigmaM;
  double thetastar;
  double sn_cutoff;
  double beta;
  double alpha_theta;

  int size_logM;

  double rho_m_at_z;
  int has_completeness;
  int mass_range;


  double lnM_min;
  double lnM_max;
  double dlnM;

  double z_0;
  double z_max;
  double dz;

  int nsteps_z;

  double logy_min;
  double logy_max;
  double dlogy;


  double lnymin;
  double lnymax;
  double dlny;

  int pvecsz_size;
  int index_y;


    ErrorMsg error_message; /**< zone for writing error messages */
};


/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int szcount_init(struct background * pba,
                   struct nonlinear * pnl,
                   struct primordial * ppm,
                   struct class_sz_structure * pclass_sz,
                   struct szcount * pcsz);


int szcounts_free(struct szcount * pcsz,struct class_sz_structure * pclass_sz);


int compute_count_sz(struct background * pba,
                     struct nonlinear * pnl,
                     struct primordial * ppm,
                     struct class_sz_structure * pclass_sz,
                     struct szcount * pcsz);

int compute_counts_sz_fft(struct background * pba,
                     struct nonlinear * pnl,
                     struct primordial * ppm,
                     struct class_sz_structure * pclass_sz,
                     struct szcount * pcsz);


  int grid_C_2d(double * pvecsz,
                struct background *pba,
                struct primordial * ppm,
                struct nonlinear * pnl,
                struct class_sz_structure * pclass_sz,
                struct szcount * pcsz);

  int write_output_cluster_counts(struct szcount * pcsz, struct class_sz_structure * pclass_sz);
  int initialise_and_allocate_memory_cc(struct class_sz_structure * pclass_sz,struct szcount * pcsz);
  int find_theta_bin(struct class_sz_structure * pclass_sz, double thp, int * l_array, double * theta_array);
  int find_y_bin(struct class_sz_structure * pclass_sz, double thp, int * l_array, double * theta_array);
double integrand_cluster_counts_redshift(double z, void *p);
double integrand_cluster_counts_mass(double lnm, void *p);
double integrand_cluster_counts_completeness(double lny, void *p);
struct Parameters_for_integrand_cluster_counts_redshift{
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * completeness_2d_to_1d;
};

struct Parameters_for_integrand_cluster_counts_mass{
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * completeness_2d_to_1d;
  double z;
};

struct Parameters_for_integrand_cluster_counts_completeness{
  struct class_sz_structure * pclass_sz;
  double * erfs_2d_to_1d;
  double theta;
  double theta1;
  double theta2;
  double y;
};

double  get_szcounts_rates_at_z_sigobs_qobs(double z_asked,
                                            double sig_asked,
                                            double qobs_asked,
                                            struct class_sz_structure * pclass_sz);

double get_szcounts_dndzdqgt_at_z_q(double z_asked, double qobs_asked, struct class_sz_structure * pclass_sz);
double  get_szcounts_dndzdq_at_z_q(double z_asked, double qobs_asked, struct class_sz_structure * pclass_sz);

#ifdef __cplusplus
}
#endif

#endif
