/** @file szclustercount.h Documented includes for sz module. */

#ifndef __SZC__
#define __SZC__

#include "common.h"
#include "szpowerspectrum.h"
#include "sz_tools.h"
#include "r8lib.h"



//#define _mean_y_ ((ptsz->has_mean_y == _TRUE_) && (index_md == ptsz->index_md_mean_y))
//#define _hmf_ ((ptsz->has_hmf == _TRUE_) && (index_md == ptsz->index_md_hmf))


struct szcount {
  FileName root; /**< root for all file names */

  int nzSZ;
  //double * dndz;
  double ** dndlnM;

  double ** dNdzdy_theoretical;
  double ** dNdzdm_theoretical;

  double ** temp_0_theoretical;
  double ** temp_1_theoretical;

  double * dvdz;
  double * redshift;
  double * z_center;
  double * steps_m;
  double * steps_z;
  double * logy;
  int nsteps_m;
  int Nbins_z;
  int Nbins_y;
  int Ny;
  short sz_verbose;
  double redshift_for_dndm;
  double * logM_at_z;

  double ystar;
  double alpha;
  double sigmaM;
  double thetastar;
  double sn_cutoff;
  double beta;
  double alpha_theta;

  int size_logM;

  double rho_m_at_z;
  int has_sz_counts;
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
                   struct tszspectrum * ptsz,
                   struct szcount * pcsz);


int szcount_free(struct szcount * pcsz);


int compute_count_sz(struct background * pba,
                     struct nonlinear * pnl,
                     struct primordial * ppm,
                     struct tszspectrum * ptsz,
                     struct szcount * pcsz,
                     double * pvecback,
                     double * Pvectsz);

  int grid_C_2d(double *** d_completeness_2d_dq,
                double *** completeness_2d,
                double *** d_erfs_dq,
                double *** erfs,
                double * pvecsz,
                struct background *pba,
                struct primordial * ppm,
                struct nonlinear * pnl,
                struct tszspectrum * ptsz,
                struct szcount * pcsz);

  int write_output_cluster_counts(struct szcount * pcsz);
  int initialise_and_allocate_memory_cc(struct tszspectrum * ptsz,struct szcount * pcsz);
  int find_theta_bin(struct tszspectrum * ptsz, double thp, int * l_array, double * theta_array);


#ifdef __cplusplus
}
#endif

#endif
