/** @file szpowerspectrum.c SZ+halo model module.
 *
 * Boris Bolliet, 2020
 */


#define _DEBUG


#include "class_sz.h"
#include "class_sz_tools.h"
#include "Patterson.h"
#include "r8lib.h"
#include "fft.h"


int class_sz_cosmo_init(  struct background * pba,
                          struct thermo * pth,
                          struct perturbs * ppt,
                          struct nonlinear * pnl,
                          struct primordial * ppm,
                          struct spectra * psp,
                          struct lensing * ple,
                          struct tszspectrum * ptsz,
                          struct precision * ppr
                        ){

// ptsz->has_sz_counts = _TRUE_;

  int all_comps = ptsz->has_sz_ps
      + ptsz->has_hmf
      + ptsz->has_n5k
      // + ptsz->has_pk_at_z_1h
      + ptsz->has_pk_at_z_1h
      + ptsz->has_pk_at_z_2h
      + ptsz->has_pk_gg_at_z_1h
      + ptsz->has_pk_gg_at_z_2h
      + ptsz->has_pk_bb_at_z_1h
      + ptsz->has_pk_bb_at_z_2h
      + ptsz->has_pk_b_at_z_2h
      + ptsz->has_gas_pressure_profile_2h
      + ptsz->has_gas_density_profile_2h
      + ptsz->has_pk_em_at_z_1h
      + ptsz->has_pk_em_at_z_2h
      + ptsz->has_pk_HI_at_z_1h
      + ptsz->has_pk_HI_at_z_2h
      + ptsz->has_bk_at_z_1h
      + ptsz->has_bk_at_z_2h
      + ptsz->has_bk_at_z_3h
      + ptsz->has_bk_ttg_at_z_1h
      + ptsz->has_bk_ttg_at_z_2h
      + ptsz->has_bk_ttg_at_z_3h
      + ptsz->has_bk_at_z_hf
      + ptsz->has_mean_y
      + ptsz->has_cib_monopole
      + ptsz->has_cib_shotnoise
      + ptsz->has_dcib0dz
      + ptsz->has_dydz
      + ptsz->has_sz_2halo
      + ptsz->has_sz_trispec
      + ptsz->has_sz_m_y_y_1h
      + ptsz->has_sz_m_y_y_2h
      + ptsz->has_sz_te_y_y
      + ptsz->has_sz_cov_N_N
      + ptsz->has_tSZ_tSZ_tSZ_1halo
      + ptsz->has_tSZ_tSZ_tSZ_2h
      + ptsz->has_tSZ_tSZ_tSZ_3h
      + ptsz->has_kSZ_kSZ_1h
      + ptsz->has_kSZ_kSZ_2h
      + ptsz->has_kSZ_kSZ_tSZ_1h
      + ptsz->has_kSZ_kSZ_tSZ_2h
      + ptsz->has_kSZ_kSZ_tSZ_3h
      + ptsz->has_kSZ_kSZ_gal_1h
      + ptsz->has_kSZ_kSZ_gal_1h_fft
      + ptsz->has_kSZ_kSZ_gal_2h_fft
      + ptsz->has_kSZ_kSZ_gal_3h_fft
      + ptsz->has_kSZ_kSZ_gal_2h
      + ptsz->has_kSZ_kSZ_gal_3h
      + ptsz->has_kSZ_kSZ_gal_hf
      + ptsz->has_kSZ_kSZ_lensmag_1halo
      + ptsz->has_kSZ_kSZ_gallens_1h_fft
      + ptsz->has_kSZ_kSZ_gallens_2h_fft
      + ptsz->has_kSZ_kSZ_gallens_3h_fft
      + ptsz->has_kSZ_kSZ_gallens_hf
      + ptsz->has_kSZ_kSZ_lens_1h_fft
      + ptsz->has_kSZ_kSZ_lens_2h_fft
      + ptsz->has_kSZ_kSZ_lens_3h_fft
      + ptsz->has_kSZ_kSZ_lens_hf
      + ptsz->has_gallens_gallens_1h
      + ptsz->has_gallens_gallens_2h
      + ptsz->has_gallens_lens_1h
      + ptsz->has_gallens_lens_2h
      + ptsz->has_tSZ_gal_1h
      + ptsz->has_tSZ_gal_2h
      + ptsz->has_tSZ_lensmag_1h
      + ptsz->has_tSZ_lensmag_2h
      + ptsz->has_tSZ_cib_1h
      + ptsz->has_tSZ_cib_2h
      + ptsz->has_lens_cib_1h
      + ptsz->has_lens_cib_2h
      + ptsz->has_gal_cib_1h
      + ptsz->has_gal_cib_2h
      + ptsz->has_cib_cib_1h
      + ptsz->has_cib_cib_2h
      + ptsz->has_ngal_ngal_1h
      + ptsz->has_ngal_ngal_2h
      + ptsz->has_ngal_ngal_hf
      + ptsz->has_ngal_lens_1h
      + ptsz->has_ngal_lens_2h
      + ptsz->has_ngal_lens_hf
      + ptsz->has_gal_gal_1h
      + ptsz->has_gal_gal_2h
      + ptsz->has_gal_gal_hf
      + ptsz->has_gal_lens_1h
      + ptsz->has_gal_lens_2h
      + ptsz->has_gal_lens_hf
      + ptsz->has_gal_lensmag_1h
      + ptsz->has_gal_lensmag_2h
      + ptsz->has_gal_gallens_1h
      + ptsz->has_gal_gallens_2h
      + ptsz->has_gal_lensmag_hf
      + ptsz->has_lensmag_lensmag_1h
      + ptsz->has_lensmag_lensmag_2h
      + ptsz->has_lensmag_lensmag_hf
      + ptsz->has_lens_lensmag_hf
      + ptsz->has_lens_lensmag_1h
      + ptsz->has_lens_lensmag_2h
      + ptsz->has_lens_lens_1h
      + ptsz->has_lens_lens_2h
      + ptsz->has_lens_lens_hf
      + ptsz->has_tSZ_lens_1h
      + ptsz->has_tSZ_lens_2h
      + ptsz->has_isw_lens
      + ptsz->has_isw_tsz
      + ptsz->has_isw_auto
      + ptsz->has_dndlnM
      + ptsz->has_vrms2
      + ptsz->has_sz_counts
      + ptsz->has_sz_rates
      + ptsz->need_m200c_to_m500c
      + ptsz->need_m500c_to_m200c
      + ptsz->need_m200m_to_m500c
      + ptsz->need_m200m_to_m200c
      + ptsz->need_m200c_to_m200m
      + ptsz->tabulate_rhob_xout_at_m_and_z
      + ptsz->need_ng_bias;
  int electron_pressure_comps = ptsz->has_sz_ps
      + ptsz->has_mean_y
      + ptsz->has_gas_pressure_profile_2h
      + ptsz->has_dydz
      + ptsz->has_sz_2halo
      + ptsz->has_sz_trispec
      + ptsz->has_sz_m_y_y_1h
      + ptsz->has_sz_m_y_y_2h
      + ptsz->has_sz_te_y_y
      + ptsz->has_tSZ_tSZ_tSZ_1halo
      + ptsz->has_tSZ_tSZ_tSZ_2h
      + ptsz->has_tSZ_tSZ_tSZ_3h
      + ptsz->has_kSZ_kSZ_tSZ_1h
      + ptsz->has_kSZ_kSZ_tSZ_2h
      + ptsz->has_kSZ_kSZ_tSZ_3h
      + ptsz->has_tSZ_gal_1h
      + ptsz->has_tSZ_gal_2h
      + ptsz->has_tSZ_lensmag_1h
      + ptsz->has_tSZ_lensmag_2h
      + ptsz->has_tSZ_cib_1h
      + ptsz->has_tSZ_cib_2h
      + ptsz->has_tSZ_lens_1h
      + ptsz->has_tSZ_lens_2h
      + ptsz->has_isw_tsz;


   // Skip the module if no SZ/halo-model computations are requested:
    if (all_comps == _FALSE_)
   {
      if (ptsz->sz_verbose > 0)
         printf("->No class_sz quantities requested - modules skipped.\n");
         return _SUCCESS_;
   }

   else
   {
     if (ptsz->sz_verbose > 0)
        printf("->Class_sz computations. Initialization.\n");


// printf("entering szp module");
    ptsz->ln_k_size_for_tSZ = (int)(log(ptsz->k_max_for_pk_in_tSZ
                                     /ptsz->k_min_for_pk_in_tSZ)
                                 /log(10.)*ptsz->k_per_decade_for_tSZ) + 2;

  class_alloc(ptsz->ln_k_for_tSZ,ptsz->ln_k_size_for_tSZ*sizeof(double),ptsz->error_message);
  int i;
  for (i=0; i<ptsz->ln_k_size_for_tSZ; i++)
      ptsz->ln_k_for_tSZ[i]=log(ptsz->k_min_for_pk_in_tSZ)+i*log(10.)/ptsz->k_per_decade_for_tSZ;


   // printf("need_hmf = %d\n",ptsz->need_hmf);
   select_multipole_array(ptsz);

   ptsz->chi_star = pth->ra_star*pba->h;
   double tau;
   int first_index_back = 0;
   double * pvecback;
   double OmegaM;

   class_alloc(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);


      class_call(background_tau_of_z(pba,
                                     0.0, //TBC: z1SZ?
                                     &tau
                                     ),
                      ptsz->error_message,
                      ptsz->error_message
                      );

     class_call(background_at_tau(pba,
                                  tau,
                                  pba->long_info,
                                  pba->inter_normal,
                                  &first_index_back,
                                  pvecback),
                      ptsz->error_message,
                      ptsz->error_message
                      );

      ptsz->Rho_crit_0 =
      (3./(8.*_PI_*_G_*_M_sun_))
      *pow(_Mpc_over_m_,1)
      *pow(_c_,2)
      *pvecback[pba->index_bg_rho_crit]
      /pow(pba->h,2);


      ptsz->Omega_m_0 = pvecback[pba->index_bg_Omega_m];
      ptsz->Omega_r_0 = pvecback[pba->index_bg_Omega_r];
      free(pvecback);

      ptsz->Omega_ncdm_0 = ptsz->Omega_m_0
      -pba->Omega0_b
      -pba->Omega0_cdm;

      ptsz->Omega0_b = pba->Omega0_b;

      if (ptsz->f_b_gas == -1.){
        ptsz->f_b_gas = pba->Omega0_b/ptsz->Omega_m_0;
      }

if (ptsz->use_class_sz_fast_mode == 0)
   show_preamble_messages(pba,pth,pnl,ppm,ptsz);
   // this routine also prints sigma8 and stuff like that
   // so it would crash if running by the emulators.


   //compute T_cmb*gNU at 150GHz
   double frequency_in_Hz = 150.0e9;

   ptsz->Tcmb_gNU_at_150GHz =
   pba->T_cmb
   *((_h_P_*frequency_in_Hz
       /(_k_B_*pba->T_cmb))
      *(1./tanh((_h_P_*frequency_in_Hz
                      /(_k_B_*pba->T_cmb))
                     /2.))
      -4.);

    frequency_in_Hz = ptsz->nu_y_dist_GHz*1e9;
    ptsz->Tcmb_gNU =
    pba->T_cmb
    *((_h_P_*frequency_in_Hz
        /(_k_B_*pba->T_cmb))
       *(1./tanh((_h_P_*frequency_in_Hz
                      /(_k_B_*pba->T_cmb))
                      /2.))-4.);


if (ptsz->need_sigma == 1 || ptsz->has_vrms2){


      double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
      // z_min = r8_min(z_min,ptsz->z_for_pk_hm);
      double z_max = 1.0001*r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
      int index_z;

      class_alloc(ptsz->array_redshift,sizeof(double *)*ptsz->n_arraySZ,ptsz->error_message);

      for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
      {
        ptsz->array_redshift[index_z] =
                                        log(1.+z_min)
                                        +index_z*(log(1.+z_max)-log(1.+z_min))
                                        /(ptsz->n_arraySZ-1.); // log(1+z)

                                      }
                                    }



   tabulate_sigma_and_dsigma_from_pk(pba,pnl,ppm,ptsz);

// if (ptsz->use_class_sz_fast_mode){
// for the class_szfast mode
  class_alloc(ptsz->array_pkl_at_z_and_k,
              sizeof(double *)*ptsz->n_arraySZ*ptsz->ndimSZ,
              ptsz->error_message);

  class_alloc(ptsz->array_pknl_at_z_and_k,
              sizeof(double *)*ptsz->n_arraySZ*ptsz->ndimSZ,
              ptsz->error_message);

  class_alloc(ptsz->array_lnk,
              sizeof(double *)*ptsz->ndimSZ,
              ptsz->error_message);

// }

   return _SUCCESS_;
   }
}


int szpowerspectrum_init(
                          struct background * pba,
                          struct thermo * pth,
                          struct perturbs * ppt,
                          struct nonlinear * pnl,
                          struct primordial * ppm,
                          struct spectra * psp,
                          struct lensing * ple,
                          struct tszspectrum * ptsz,
                          struct precision * ppr
			                    )
{

// ptsz->has_sz_counts = _TRUE_;

  int all_comps = ptsz->has_sz_ps
      + ptsz->has_hmf
      + ptsz->has_n5k
      // + ptsz->has_pk_at_z_1h
      + ptsz->has_pk_at_z_1h
      + ptsz->has_pk_at_z_2h
      + ptsz->has_pk_gg_at_z_1h
      + ptsz->has_pk_gg_at_z_2h
      + ptsz->has_pk_bb_at_z_1h
      + ptsz->has_pk_bb_at_z_2h
      + ptsz->has_pk_b_at_z_2h
      + ptsz->has_gas_pressure_profile_2h
      + ptsz->has_gas_density_profile_2h
      + ptsz->has_pk_em_at_z_1h
      + ptsz->has_pk_em_at_z_2h
      + ptsz->has_pk_HI_at_z_1h
      + ptsz->has_pk_HI_at_z_2h
      + ptsz->has_bk_at_z_1h
      + ptsz->has_bk_at_z_2h
      + ptsz->has_bk_at_z_3h
      + ptsz->has_bk_ttg_at_z_1h
      + ptsz->has_bk_ttg_at_z_2h
      + ptsz->has_bk_ttg_at_z_3h
      + ptsz->has_bk_at_z_hf
      + ptsz->has_mean_y
      + ptsz->has_cib_monopole
      + ptsz->has_cib_shotnoise
      + ptsz->has_dcib0dz
      + ptsz->has_dydz
      + ptsz->has_sz_2halo
      + ptsz->has_sz_trispec
      + ptsz->has_sz_m_y_y_1h
      + ptsz->has_sz_m_y_y_2h
      + ptsz->has_sz_te_y_y
      + ptsz->has_sz_cov_N_N
      + ptsz->has_tSZ_tSZ_tSZ_1halo
      + ptsz->has_tSZ_tSZ_tSZ_2h
      + ptsz->has_tSZ_tSZ_tSZ_3h
      + ptsz->has_kSZ_kSZ_1h
      + ptsz->has_kSZ_kSZ_2h
      + ptsz->has_kSZ_kSZ_tSZ_1h
      + ptsz->has_kSZ_kSZ_tSZ_2h
      + ptsz->has_kSZ_kSZ_tSZ_3h
      + ptsz->has_kSZ_kSZ_gal_1h
      + ptsz->has_kSZ_kSZ_gal_1h_fft
      + ptsz->has_kSZ_kSZ_gal_2h_fft
      + ptsz->has_kSZ_kSZ_gal_3h_fft
      + ptsz->has_kSZ_kSZ_gal_2h
      + ptsz->has_kSZ_kSZ_gal_3h
      + ptsz->has_kSZ_kSZ_gal_hf
      + ptsz->has_kSZ_kSZ_lensmag_1halo
      + ptsz->has_kSZ_kSZ_gallens_1h_fft
      + ptsz->has_kSZ_kSZ_gallens_2h_fft
      + ptsz->has_kSZ_kSZ_gallens_3h_fft
      + ptsz->has_kSZ_kSZ_gallens_hf
      + ptsz->has_kSZ_kSZ_lens_1h_fft
      + ptsz->has_kSZ_kSZ_lens_2h_fft
      + ptsz->has_kSZ_kSZ_lens_3h_fft
      + ptsz->has_kSZ_kSZ_lens_hf
      + ptsz->has_gallens_gallens_1h
      + ptsz->has_gallens_gallens_2h
      + ptsz->has_gallens_lens_1h
      + ptsz->has_gallens_lens_2h
      + ptsz->has_tSZ_gal_1h
      + ptsz->has_tSZ_gal_2h
      + ptsz->has_tSZ_lensmag_1h
      + ptsz->has_tSZ_lensmag_2h
      + ptsz->has_tSZ_cib_1h
      + ptsz->has_tSZ_cib_2h
      + ptsz->has_lens_cib_1h
      + ptsz->has_lens_cib_2h
      + ptsz->has_gal_cib_1h
      + ptsz->has_gal_cib_2h
      + ptsz->has_cib_cib_1h
      + ptsz->has_cib_cib_2h
      + ptsz->has_ngal_ngal_1h
      + ptsz->has_ngal_ngal_2h
      + ptsz->has_ngal_ngal_hf
      + ptsz->has_ngal_lens_1h
      + ptsz->has_ngal_lens_2h
      + ptsz->has_ngal_lens_hf
      + ptsz->has_gal_gal_1h
      + ptsz->has_gal_gal_2h
      + ptsz->has_gal_gal_hf
      + ptsz->has_gal_lens_1h
      + ptsz->has_gal_lens_2h
      + ptsz->has_gal_lens_hf
      + ptsz->has_gal_lensmag_1h
      + ptsz->has_gal_lensmag_2h
      + ptsz->has_gal_gallens_1h
      + ptsz->has_gal_gallens_2h
      + ptsz->has_gal_lensmag_hf
      + ptsz->has_lensmag_lensmag_1h
      + ptsz->has_lensmag_lensmag_2h
      + ptsz->has_lensmag_lensmag_hf
      + ptsz->has_lens_lensmag_hf
      + ptsz->has_lens_lensmag_1h
      + ptsz->has_lens_lensmag_2h
      + ptsz->has_lens_lens_1h
      + ptsz->has_lens_lens_2h
      + ptsz->has_lens_lens_hf
      + ptsz->has_tSZ_lens_1h
      + ptsz->has_tSZ_lens_2h
      + ptsz->has_isw_lens
      + ptsz->has_isw_tsz
      + ptsz->has_isw_auto
      + ptsz->has_dndlnM
      + ptsz->has_vrms2
      + ptsz->has_sz_counts
      + ptsz->has_sz_rates
      + ptsz->need_m200c_to_m500c
      + ptsz->need_m500c_to_m200c
      + ptsz->need_m200m_to_m500c
      + ptsz->need_m200m_to_m200c
      + ptsz->need_m200c_to_m200m
      + ptsz->tabulate_rhob_xout_at_m_and_z
      + ptsz->need_ng_bias;
  int electron_pressure_comps = ptsz->has_sz_ps
      + ptsz->has_mean_y
      + ptsz->has_gas_pressure_profile_2h
      + ptsz->has_dydz
      + ptsz->has_sz_2halo
      + ptsz->has_sz_trispec
      + ptsz->has_sz_m_y_y_1h
      + ptsz->has_sz_m_y_y_2h
      + ptsz->has_sz_te_y_y
      + ptsz->has_tSZ_tSZ_tSZ_1halo
      + ptsz->has_tSZ_tSZ_tSZ_2h
      + ptsz->has_tSZ_tSZ_tSZ_3h
      + ptsz->has_kSZ_kSZ_tSZ_1h
      + ptsz->has_kSZ_kSZ_tSZ_2h
      + ptsz->has_kSZ_kSZ_tSZ_3h
      + ptsz->has_tSZ_gal_1h
      + ptsz->has_tSZ_gal_2h
      + ptsz->has_tSZ_lensmag_1h
      + ptsz->has_tSZ_lensmag_2h
      + ptsz->has_tSZ_cib_1h
      + ptsz->has_tSZ_cib_2h
      + ptsz->has_tSZ_lens_1h
      + ptsz->has_tSZ_lens_2h
      + ptsz->has_isw_tsz;


   // Skip the module if no SZ/halo-model computations are requested:
    if (all_comps == _FALSE_)
   {
      if (ptsz->sz_verbose > 0)
         printf("->No class_sz quantities requested - modules skipped.\n");
         return _SUCCESS_;
   }

   else
   {
     if (ptsz->sz_verbose > 0)
        printf("->Class_sz computations. Initialization.\n");


// // printf("entering szp module");
//     ptsz->ln_k_size_for_tSZ = (int)(log(ptsz->k_max_for_pk_in_tSZ
//                                      /ptsz->k_min_for_pk_in_tSZ)
//                                  /log(10.)*ptsz->k_per_decade_for_tSZ) + 2;
//
//   class_alloc(ptsz->ln_k_for_tSZ,ptsz->ln_k_size_for_tSZ*sizeof(double),ptsz->error_message);
//   int i;
//   for (i=0; i<ptsz->ln_k_size_for_tSZ; i++)
//       ptsz->ln_k_for_tSZ[i]=log(ptsz->k_min_for_pk_in_tSZ)+i*log(10.)/ptsz->k_per_decade_for_tSZ;
//
//
//    // printf("need_hmf = %d\n",ptsz->need_hmf);
//    select_multipole_array(ptsz);
//
//
//
//    show_preamble_messages(pba,pth,pnl,ppm,ptsz);
//    if (ptsz->need_sigma == 1
//     || ptsz->has_vrms2){
//
//
//       double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
//       // z_min = r8_min(z_min,ptsz->z_for_pk_hm);
//       double z_max = 1.0001*r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
//       int index_z;
//
//       class_alloc(ptsz->array_redshift,sizeof(double *)*ptsz->n_arraySZ,ptsz->error_message);
//
//       for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
//       {
//         ptsz->array_redshift[index_z] =
//                                         log(1.+z_min)
//                                         +index_z*(log(1.+z_max)-log(1.+z_min))
//                                         /(ptsz->n_arraySZ-1.); // log(1+z)
//
//                                       }
//                                     }
//
//
//
//    tabulate_sigma_and_dsigma_from_pk(pba,pnl,ppm,ptsz);





// // if (ptsz->use_class_sz_fast_mode){
// // for the class_szfast mode
//   class_alloc(ptsz->array_pkl_at_z_and_k,
//               sizeof(double *)*ptsz->n_arraySZ*ptsz->ndimSZ,
//               ptsz->error_message);
//
//   class_alloc(ptsz->array_pknl_at_z_and_k,
//               sizeof(double *)*ptsz->n_arraySZ*ptsz->ndimSZ,
//               ptsz->error_message);
//
//   class_alloc(ptsz->array_lnk,
//               sizeof(double *)*ptsz->ndimSZ,
//               ptsz->error_message);
//
// // }

// begin tk stuff
// start collecting transfer functions
char titles[_MAXTITLESTRINGLENGTH_]={0};
int size_data, number_of_titles;
int index_md;
int index_title;
int index_d_tot;
int index_phi;
int index_psi;

class_call(perturb_output_titles(pba,ppt,class_format,titles),
           ptsz->error_message,
           ptsz->error_message);

// printf("ok titles diones\n");

ptsz->number_of_titles = get_number_of_titles(titles);
// printf("number_of_titles  = %d %s\n",
// number_of_titles,
// titles);
char *pch;
pch = strtok(titles,"\t");
int idp = 0;

while( pch != NULL ) {
      // printf( "%s\n",pch );
      // strcpy(string1,pch);
      // printf( "%s\n",string1 );
      if (strstr(pch,"d_m") != 0)
        ptsz->index_d_tot = idp;
      if (strstr(pch,"phi") != 0)
        ptsz->index_phi = idp;
      if (strstr(pch,"psi") != 0)
        ptsz->index_psi = idp;
      pch = strtok(NULL,"\t");
      idp+=1;
   }

////// end tk stuff



   if (ptsz->has_sz_rates || ptsz->has_sz_counts_fft){
      read_sz_catalog(ptsz);

      if (ptsz->sz_verbose>1){
        int icat = 0;
        int imissz = 0;
        printf("szcat: got %d lines.\n",ptsz->szcat_size);
        for (icat=0;icat<ptsz->szcat_size;icat++){
            if (ptsz->sz_verbose>3) printf("szcat z = %.3e \t snr = %.3e\n",ptsz->szcat_z[icat],ptsz->szcat_snr[icat]);
          if (ptsz->szcat_z[icat]<= 0) imissz += 1;
        }
        printf("szcat: got %d objects with missing redshift.\n",imissz);
      }
    }

   // exit(0);
   if (ptsz->sz_verbose>=1)
    printf("-> allocating class_sz memory...\n");
   initialise_and_allocate_memory(ptsz);
   if (ptsz->sz_verbose>=1)
    printf("-> memory allocated.\n");



   if ((ptsz->has_completeness_for_ps_SZ == 1)  || (ptsz->has_sz_counts  == 1)){
      read_Planck_noise_map(ptsz);
    }
    // exit(0);

   if (ptsz->concentration_parameter == 4)
      read_Zhao_CM_init(ptsz);




      if (ptsz->use_maniyar_cib_model && ptsz->has_cib){
        load_cib_Snu(ptsz);
      //   exit(0);
      double K1_interp = get_cib_Snu_z_and_nu(0.23,95,ptsz);
      // printf("%.5e\n",K1_interp);
      // exit(0);
      }

   tabulate_L_sat_at_z_m_nu(pba,ptsz);
   // printf("%.8e\n",get_L_sat_at_z_M_nu(5.04783496e-01,2.91975583e+16,3.62095835e+02,ptsz)); // should fund  1.78986876e+07
   // exit(0);

      // tabulate_L_sat_at_nu_and_nu_prime(pba,ptsz);



  // printf("tabulating dndlnM quantities %d\n",ptsz->has_sigma2_hsv);
   if (ptsz->has_sigma2_hsv)
   tabulate_sigma2_hsv_from_pk(pba,pnl,ppm,ptsz);


   if (ptsz->need_m200c_to_m200m == 1){
    if (ptsz->sz_verbose>1)
       printf("-> tabulating m200c to m200m...\n");

      tabulate_m200c_to_m200m(pba,ptsz);

    if (ptsz->sz_verbose>1)
      printf("-> m200c to m200m tabulated.\n");
    }
// exit(0);
   if (ptsz->need_m200m_to_m200c == 1){
     if (ptsz->sz_verbose>1)
        printf("-> tabulating m200m to m200c...\n");
      tabulate_m200m_to_m200c(pba,ptsz);
    if (ptsz->sz_verbose>1)
      printf("-> m200m to m200c tabulated.\n");
    }


   if (ptsz->need_m200m_to_m500c == 1){
     if (ptsz->sz_verbose>1)
     printf("-> tabulating m200m to m500c...\n");
      tabulate_m200m_to_m500c(pba,ptsz);
     if (ptsz->sz_verbose>1)
     printf("-> m200m to m500c tabulated.\n");
    }

   if (ptsz->need_m200c_to_m500c == 1){
     if (ptsz->sz_verbose>1)
     printf("-> tabulating m200c to m500c...\n");
      tabulate_m200c_to_m500c(pba,ptsz);
     if (ptsz->sz_verbose>1)
     printf("-> m200c to m500c tabulated.\n");
    }

   if (ptsz->need_m500c_to_m200c == 1){
     if (ptsz->sz_verbose>1)
     printf("-> tabulating m500c to m200c...\n");
      tabulate_m500c_to_m200c(pba,ptsz);
     if (ptsz->sz_verbose>1)
     printf("-> m500c to m200c tabulated.\n");
    }
   //exit(0);
   external_pressure_profile_init(ppr,ptsz);

if (ptsz->MF==1){
   // load alpha(z) normalisation for Tinker et al 2010 HMF
   if (ptsz->T10_alpha_fixed==0){
   if (ptsz->sz_verbose>1)
    printf("-> loading alpha(z) for Tinker al 2010 HMF...\n");
   load_T10_alpha_norm(ptsz);
   if (ptsz->sz_verbose>1)
    printf("-> alpha(z) for Tinker al HMF loaded.\n");
   }
}

   if (ptsz->has_dndlnM == 1
    || ptsz->has_sz_counts
    || ptsz->has_sz_rates
    || ptsz->has_kSZ_kSZ_gal_1h
    || ptsz->has_kSZ_kSZ_gal_1h_fft
    || ptsz->has_kSZ_kSZ_gal_2h_fft
    || ptsz->has_kSZ_kSZ_gal_3h_fft){

if (ptsz->sz_verbose>1)
   printf("-> Tabulating dndlnM HMF in mass and redshift...\n");

  tabulate_dndlnM(pba,pnl,ppm,ptsz);

if (ptsz->sz_verbose>1)
   printf("-> dndlnM HMF tabulated.\n");



  }
// exit(0);
  // printf("tabulating dndlnM quantities -1\n");


  if (ptsz->need_ng_bias){
    if (ptsz->sz_verbose>1)
       printf("-> Tabulating scale dependent bias...\n");
    tabulate_ng_bias_contribution_at_z_and_k(pba,ppt,ptsz);
  }
  // exit(0);

  // printf("ok for now...\n");

if ((ptsz->need_hmf != 0) && (ptsz->hm_consistency==1)){
  if (ptsz->sz_verbose>10)
      printf("counter terms nmin\n");
   tabulate_hmf_counter_terms_nmin(pba,pnl,ppm,ptsz);
  if (ptsz->sz_verbose>10)
  printf("counter terms b1 min\n");
   tabulate_hmf_counter_terms_b1min(pba,pnl,ppm,ppt,ptsz);
   if (ptsz->sz_verbose>10)
   printf("counter terms b1 min done\n");
   // if (ptsz->hm_consistency==1){
   ptsz->hm_consistency_counter_terms_done = 0;
   // tabulate_hmf_counter_terms_nmin(pba,pnl,ppm,ptsz);
   // tabulate_hmf_counter_terms_b1min(pba,pnl,ppm,ptsz);
   tabulate_hmf_counter_terms_b2min(pba,pnl,ppm,ptsz);
   // printf("tabulating dndlnM quantities -1\n");
   ptsz->hm_consistency_counter_terms_done = 1;
   if (ptsz->sz_verbose>10){
   printf("counter terms\n");
   int index_z;
   for (index_z=0; index_z<ptsz->n_z_hmf_counter_terms; index_z++){
   double z =  exp(ptsz->array_redshift_hmf_counter_terms[index_z])-1.;
   double n_min = get_hmf_counter_term_nmin_at_z(z,ptsz);
   double b1_min = get_hmf_counter_term_b1min_at_z(z,ptsz);
   double b2_min = get_hmf_counter_term_b2min_at_z(z,ptsz);
   printf("z = %.3e n_min = %.8e n_min_interp = %.8e b1_min = %.8e b1_min_interp = %.8e b2_min = %.8e b2_min_interp = %.8e\n",
   z, ptsz->array_hmf_counter_terms_nmin[index_z],n_min,
      ptsz->array_hmf_counter_terms_b1min[index_z],b1_min,
      ptsz->array_hmf_counter_terms_b2min[index_z],b2_min);
    }
  }
  // }
}


   if (ptsz->has_vrms2){
if (ptsz->sz_verbose>1)
    printf("-> Tabulating velocity dispersion...\n");
   tabulate_vrms2_from_pk(pba,pnl,ppm,ptsz);
if (ptsz->sz_verbose>1)
   printf("-> Velocity dispersion tabulated.\n");
 }

//  printf("get_vrms2_at_z = %.5e\n",get_vrms2_at_z(0.3,ptsz));
// exit(0);
   // printf("tabulating dndlnM quantities 0\n");

   if (ptsz->has_knl){
if (ptsz->sz_verbose>1)
   printf("-> Tabulating knl...\n");
   tabulate_knl(pba,pnl,ppm,ptsz);
if (ptsz->sz_verbose>1)
  printf("-> knl tabulated.\n");
 }


   // printf("tabulating dndlnM quantities 1\n");

   if (ptsz->has_nl_index){
if (ptsz->sz_verbose>1)
    printf("-> Tabulating nl index...\n");
   tabulate_nl_index(pba,pnl,ppm,ptsz);
if (ptsz->sz_verbose>1)
    printf("-> nl index tabulated.\n");
 }


// printf("-> tabulating xout for Battaglia density profile %d.\n",ptsz->tabulate_rhob_xout_at_m_and_z);
if (ptsz->has_electron_density == 1 || ptsz->tabulate_rhob_xout_at_m_and_z == 1){
      if (ptsz->use_xout_in_density_profile_from_enclosed_mass || ptsz->tabulate_rhob_xout_at_m_and_z){
      if (ptsz->sz_verbose>1)
        printf("-> tabulating xout for Battaglia density profile.\n");
      tabulate_m_to_xout(pba,pnl,ppm,ptsz);
      if (ptsz->sz_verbose>1)
        printf("-> xout for Battaglia density profile tabulated.\n");

      // test:
      // double xout_test = get_m_to_xout_at_z_and_m(5.22863,6.12609e11,ptsz);
      // printf("%.5e\n",xout_test);
      }
  }


   // exit(0);
if (
     ptsz->has_kSZ_kSZ_gallens_1h_fft
  || ptsz->has_kSZ_kSZ_gallens_2h_fft
  || ptsz->has_kSZ_kSZ_gallens_3h_fft
  || ptsz->has_kSZ_kSZ_gallens_hf
  || ptsz->has_gal_gallens_1h
  || ptsz->has_gal_gallens_2h
  || ptsz->has_gallens_gallens_1h
  || ptsz->has_gallens_gallens_2h
  || ptsz->has_gallens_lens_1h
  || ptsz->has_gallens_lens_2h

){

  load_normalized_source_dndz(ptsz);
}

   if (ptsz->has_tSZ_gal_1h
    || ptsz->has_tSZ_gal_2h
    || ptsz->has_kSZ_kSZ_gal_1h
    || ptsz->has_kSZ_kSZ_gal_1h_fft
    || ptsz->has_kSZ_kSZ_gal_2h_fft
    || ptsz->has_kSZ_kSZ_gal_3h_fft
    || ptsz->has_kSZ_kSZ_gal_2h
    || ptsz->has_kSZ_kSZ_gal_3h
    || ptsz->has_kSZ_kSZ_gal_hf
    || ptsz->has_bk_ttg_at_z_1h
    || ptsz->has_bk_ttg_at_z_2h
    || ptsz->has_bk_ttg_at_z_3h
    || ptsz->has_kSZ_kSZ_lensmag_1halo
    || ptsz->has_gal_gal_1h
    || ptsz->has_gal_gal_2h
    || ptsz->has_gal_gal_hf
    || ptsz->has_gal_lens_hf
    || ptsz->has_gal_lens_1h
    || ptsz->has_gal_lens_2h
    || ptsz->has_gal_cib_1h
    || ptsz->has_gal_cib_2h
    || ptsz->has_gal_lensmag_1h
    || ptsz->has_gal_lensmag_2h
    || ptsz->has_gal_gallens_1h
    || ptsz->has_gal_gallens_2h
    || ptsz->has_gal_lensmag_hf
    || ptsz->has_tSZ_lensmag_1h
    || ptsz->has_tSZ_lensmag_2h
    || ptsz->has_lensmag_lensmag_1h
    || ptsz->has_lensmag_lensmag_2h
    || ptsz->has_lensmag_lensmag_hf
    || ptsz->has_lens_lensmag_1h
    || ptsz->has_lens_lensmag_2h
    || ptsz->has_lens_lensmag_hf
  ){

// only performed if requested:
load_normalized_dndz(ptsz);
// if (  ptsz->has_gal_gallens_1h
//    || ptsz->has_gal_gallens_2h
//    || ptsz->has_gallens_gallens_1h
//    || ptsz->has_gallens_gallens_2h
//    || ptsz->has_gallens_lens_1h
//    || ptsz->has_gallens_lens_2h){
// load_normalized_source_dndz(ptsz);
//     }
//unwise
if(ptsz->galaxy_sample==1){
      load_normalized_fdndz(ptsz);
      load_normalized_cosmos_dndz(ptsz);
    }
}

if (ptsz->has_ngal_ngal_1h
  + ptsz->has_ngal_ngal_2h
  + ptsz->has_ngal_ngal_hf
  + ptsz->has_ngal_lens_1h
  + ptsz->has_ngal_lens_2h
  + ptsz->has_ngal_lens_hf
){
 load_normalized_dndz_ngal(ptsz);
}

if (ptsz->has_kSZ_kSZ_gal_1h
 || ptsz->has_kSZ_kSZ_gal_1h_fft
 || ptsz->has_kSZ_kSZ_gal_2h_fft
 || ptsz->has_kSZ_kSZ_gal_3h_fft
 || ptsz->has_kSZ_kSZ_gal_covmat
 || ptsz->has_kSZ_kSZ_gal_lensing_term
 || ptsz->has_kSZ_kSZ_gallens_1h_fft
 || ptsz->has_kSZ_kSZ_gallens_2h_fft
 || ptsz->has_kSZ_kSZ_gallens_3h_fft
 || ptsz->has_kSZ_kSZ_gallens_covmat
 || ptsz->has_kSZ_kSZ_gallens_lensing_term
 || ptsz->has_kSZ_kSZ_gallens_hf
 || ptsz->has_kSZ_kSZ_lens_1h_fft
 || ptsz->has_kSZ_kSZ_lens_2h_fft
 || ptsz->has_kSZ_kSZ_lens_3h_fft
 || ptsz->has_kSZ_kSZ_lens_covmat
 || ptsz->has_kSZ_kSZ_lens_lensing_term
 || ptsz->has_kSZ_kSZ_lens_hf
 || ptsz->has_kSZ_kSZ_lensmag_1halo
 || ptsz->has_kSZ_kSZ_gal_2h
 || ptsz->has_kSZ_kSZ_gal_3h
 || ptsz->has_kSZ_kSZ_gal_hf)
load_ksz_filter(ptsz);


if (ptsz->has_tSZ_gal_1h
 || ptsz->has_tSZ_gal_2h

 || ptsz->has_kSZ_kSZ_gal_1h_fft
 || ptsz->has_kSZ_kSZ_gal_2h_fft
 || ptsz->has_kSZ_kSZ_gal_3h_fft
 || ptsz->has_kSZ_kSZ_gal_1h
 || ptsz->has_kSZ_kSZ_gal_2h
 || ptsz->has_kSZ_kSZ_gal_3h
 || ptsz->has_bk_ttg_at_z_1h
 || ptsz->has_bk_ttg_at_z_2h
 || ptsz->has_bk_ttg_at_z_3h
 || ptsz->has_kSZ_kSZ_gal_hf
 || ptsz->has_kSZ_kSZ_lensmag_1halo //not needed??
 || ptsz->has_gal_gal_1h
 || ptsz->has_gal_gal_2h
 || ptsz->has_gal_cib_1h
 || ptsz->has_gal_cib_2h
 || ptsz->has_pk_gg_at_z_1h
 || ptsz->has_pk_gg_at_z_2h
 || ptsz->has_gal_lens_1h
 || ptsz->has_gal_lens_2h
 || ptsz->has_gal_lensmag_1h
 || ptsz->has_gal_lensmag_2h
 || ptsz->has_gal_gallens_1h
 || ptsz->has_gal_gallens_2h
 // || ptsz->has_tSZ_lensmag_1h //not needed??
 // || ptsz->has_tSZ_lensmag_2h //not needed??
 // || ptsz->has_lensmag_lensmag_1h //not needed??
 // || ptsz->has_lensmag_lensmag_2h //not needed??
 // || ptsz->has_lens_lensmag_1h //not needed??
 // || ptsz->has_lens_lensmag_2h //not needed??
){

tabulate_mean_galaxy_number_density(pba,pnl,ppm,ptsz);
}

if (ptsz->has_ngal_ngal_1h
   +ptsz->has_ngal_ngal_2h
   +ptsz->has_ngal_lens_1h
   +ptsz->has_ngal_lens_2h
 ){
  tabulate_mean_galaxy_number_density_ngal(pba,pnl,ppm,ptsz);
}

if (ptsz->has_mean_galaxy_bias)
{
  tabulate_mean_galaxy_bias(pba,pnl,ppm,ppt,ptsz);
}


if (ptsz->use_fft_for_profiles_transform){

  // printf("-> start tabulation of gas pressure profile  444.\n");
  tabulate_gas_density_profile_fft(pba,ptsz);

  // printf("-> start tabulation of gas pressure profile  444 done.\n");
//
//   int n_k = ptsz->n_k_density_profile;
//   int n_m = ptsz->n_m_density_profile;
//   int n_z = ptsz->n_z_density_profile;
// //
// int index_ztest, index_mtest, index_ktest;
// int ik = 0;
// for (index_ztest = n_z - 5;index_ztest<n_z;index_ztest++){
//   for (index_mtest = n_m - 5;index_mtest<n_m;index_mtest++){
//     for (index_ktest = n_k - 5;index_ktest<n_k;index_ktest++){
//       double m_asked = exp(ptsz->array_profile_ln_m[index_mtest]);
//       double z_asked = exp(ptsz->array_profile_ln_1pz[index_ztest])-1.;;
//       double k_asked = exp(ptsz->array_profile_ln_k[index_ktest]);
//       double result = get_gas_density_profile_at_k_M_z(k_asked,m_asked,z_asked,ptsz);
//       printf("res fft routine i = %d k = %.3e m = %.3e z = %.3e = %.10e\n",ik,k_asked, m_asked, z_asked, result/m_asked/ptsz->f_b_gas);
//       ik ++;
//     }
//   }
// }
 // double m_asked = 5e14;
 // double z_asked = 0.1;
 // double k_asked = 0.6;
 // double result = get_gas_density_profile_at_k_M_z(k_asked,m_asked,z_asked,ptsz);
 //  printf("res = %.10e\n",result/m_asked/ptsz->f_b_gas);
 // exit(0);
     // printf("l_asked %.8e")
 // l = 5.00000000e+04 m = 1.00000000e+18 z = 5.66598637e+00 lnrho = -9.90530714e+00
 // l = 4.31674192e+04 m = 4.34701316e+17 z = 1.25153177e+00 lnrho = -9.22605976e+00
  //    double result = get_gas_density_profile_at_k_M_z(5.00000000e+04,1.00000000e+18,5.66598637e+00,ptsz);
  //    printf("%.8e\n",log(result));
  //
  //
  // printf("##################\n");
  // exit(0);


}
else{


// printf("-> start tabulation of gas pressure profile  444.\n");
 // tabulate density, only when requested (e.g., kSZ)

 tabulate_gas_density_profile(pba,ptsz);

// ik = 0;
// for (index_ztest = n_z - 5;index_ztest<n_z;index_ztest++){
//   for (index_mtest = n_m - 5;index_mtest<n_m;index_mtest++){
//     for (index_ktest = n_k - 5;index_ktest<n_k;index_ktest++){
//       double m_asked = exp(ptsz->array_profile_ln_m[index_mtest]);
//       double z_asked = exp(ptsz->array_profile_ln_1pz[index_ztest])-1.;;
//       double k_asked = exp(ptsz->array_profile_ln_k[index_ktest]);
//       double result = get_gas_density_profile_at_k_M_z(k_asked,m_asked,z_asked,ptsz);
//       printf("res original routine i = %d k = %.3e m = %.3e z = %.3e = %.10e\n",ik, k_asked, m_asked, z_asked, result/m_asked/ptsz->f_b_gas);
//       ik ++;
//     }
//   }
// }

//  if (ptsz->check_consistency_conditions == 1){
 // printf("checking normalization of profile\n");
//  // the normalization of the profile should be m_delta in the limit k->0.
//  result = get_gas_density_profile_at_k_M_z(k_asked,m_asked,z_asked,ptsz);
// //
//
// printf("res = %.10e\n",result/m_asked/ptsz->f_b_gas);
 // exit(0);
// }
// double m_asked = 5e14;
// double z_asked = 0.1;
// double k_asked = 0.6;
// double result = get_gas_density_profile_at_k_M_z(k_asked,m_asked,z_asked,ptsz);
// printf("res = %.10e\n",result/m_asked/ptsz->f_b_gas);
//  exit(0);
}



 // tabulate pressure profile for gnFW
  // printf("tab \n");
  // printf("-> start tabulation of gas pressure profile.\n");
  // printf("electron_pressure_comps = %d\n",electron_pressure_comps);
if (electron_pressure_comps != _FALSE_){
  // printf("-> start tabulation of gas pressure profile.\n");
 if (ptsz->pressure_profile == 3)
 tabulate_gas_pressure_profile_gNFW(pba,ptsz);
 else if (ptsz->pressure_profile == 4){
if (ptsz->use_fft_for_profiles_transform){
 tabulate_gas_pressure_profile_B12_fft(pba,ptsz);
 // exit(0);
}
else{
 tabulate_gas_pressure_profile_B12(pba,ptsz);
}
}

}


if (ptsz->has_dcib0dz){
  tabulate_dcib0dz(pba,pnl,ppm,ptsz);
  // printf("%.8e\n",get_dcib0dz_at_z_and_nu(1.,500.,ptsz));
}

if (ptsz->has_dydz){
  tabulate_dydz(pba,pnl,ppm,ptsz);
  // printf("%.8e\n",get_dydz_at_z(1.,ptsz));
}


// tabulate lensing magnificaion integral, only when requested
tabulate_redshift_int_lensmag(ptsz,pba);

tabulate_redshift_int_gallens_sources(ptsz,pba);

// only when requested:
load_unbinned_nl_yy(ptsz);

if (ptsz->has_kSZ_kSZ_gal_1h_fft
 || ptsz->has_kSZ_kSZ_gal_2h_fft
 || ptsz->has_kSZ_kSZ_gal_3h_fft
 || ptsz->has_kSZ_kSZ_gal_3h
){
tabulate_psi_b1g(pba,pnl,ppm,ppt,ptsz);
tabulate_psi_b2g(pba,pnl,ppm,ptsz);
tabulate_psi_b1gt(pba,pnl,ppm,ppt,ptsz);
}


if (ptsz->has_kSZ_kSZ_gallens_1h_fft
 || ptsz->has_kSZ_kSZ_gallens_2h_fft
 || ptsz->has_kSZ_kSZ_gallens_3h_fft
 || ptsz->has_kSZ_kSZ_lens_1h_fft
 || ptsz->has_kSZ_kSZ_lens_2h_fft
 || ptsz->has_kSZ_kSZ_lens_3h_fft
){

tabulate_psi_b1kg(pba,pnl,ppm,ppt,ptsz);
tabulate_psi_b2kg(pba,pnl,ppm,ptsz);
tabulate_psi_b1kgt(pba,pnl,ppm,ppt,ptsz);
}




if (ptsz->has_n5k){
  load_n5k_pk_zk(ptsz);
  load_n5k_cl_K1(ptsz);
  load_n5k_z_of_chi(ptsz);
  // load_n5k_cl_K2(ptsz);
  // double pk_interp = get_n5k_pk_at_z_and_k(3.428571428571428381e+00,9.329304026284677320e+01,ptsz);
  // printf("%.5e\n",pk_interp);
  // double K1_interp = get_n5k_cl_K1_at_chi(6.565659574734836212e+01,ptsz);
  // printf("%.5e\n",K1_interp);
  // double z_interp = get_n5k_z_of_chi(33.,ptsz);
  // printf("%.5e\n",z_interp);
  tabulate_n5k_F1(pba,pnl,ppm,ptsz);
  int index_l;
  printf("\n");
  // printf("ell\n");
  // for (index_l=0; index_l<ptsz->n_l_n5k; index_l++){
  //
  //   printf("%d,\n",ptsz->array_n5k_F1_l[index_l]);
  // }
  printf("\n");
  printf("cls\n");
  char Filepath[_ARGUMENT_LENGTH_MAX_];
  FILE *fp;

  sprintf(Filepath,"%s%s%s",ptsz->root,"n5k_F",".txt");
  fp=fopen(Filepath, "w");
  // char Filepath[_ARGUMENT_LENGTH_MAX_];
  for (index_l=0; index_l<ptsz->n_l_n5k; index_l++){

    fprintf(fp,"%.5e\n",ptsz->array_n5k_F1_F[index_l]);
  }
  fclose(fp);

}
// exit(0);
  // printf("-> start tabulation of gas pressure profile2h. %d\n",ptsz->has_gas_pressure_profile_2h);
if (ptsz->has_gas_pressure_profile_2h){
// printf("-> starting tabulation of pressure profile 2h\n");

tabulate_gas_pressure_profile_2h(pba,pnl,ppm,ppt,ptsz);
// double k_test = 0.36e-1;
// double z_test = 1.51;
// double rho_test = get_gas_pressure_2h_at_k_and_z(k_test,z_test,ptsz);
// printf("k_test = %.3e, z_test = %.3e, rho_test = %.8e\n",
//         k_test,z_test,rho_test);
// printf("-> starting tabulation of pressure profile fft 2h\n");

tabulate_gas_pressure_profile_2h_fft_at_z_and_r(pba,pnl,ppm,ptsz);


// double r_test =  2.42013e-01;
// double z_test = 1.20000e+00;
// double m_test = 3.5e13;
// double rho_test = get_gas_pressure_2h_at_r_and_m_and_z(r_test,m_test,z_test,ptsz,pba);
// printf("r_test = %.5e, z_test = %.5e, rho_test = %.8e\n",
//         r_test,z_test,rho_test);
}
// exit(0);

if (ptsz->has_pk_b_at_z_2h
   +ptsz->has_gas_density_profile_2h){
tabulate_gas_density_profile_2h(pba,pnl,ppm,ppt,ptsz);
// double k_test = 0.36e-1;
// double z_test = 1.51;
// double rho_test = get_rho_2h_at_k_and_z(k_test,z_test,ptsz);
// printf("k_test = %.3e, z_test = %.3e, rho_test = %.8e\n",
//         k_test,z_test,rho_test);
//
// k_test =  2.42013e-01;
// z_test = 3.00000e+00;
// rho_test = get_rho_2h_at_k_and_z(k_test,z_test,ptsz);
// printf("k_test = %.5e, z_test = %.5e, rho_test = %.8e\n",
//         k_test,z_test,rho_test);


tabulate_gas_density_profile_2h_fft_at_z_and_r(pba,pnl,ppm,ptsz);


// double r_test =  2.42013e-01;
// double z_test = 1.20000e+00;
// double m_test = 3.5e13;
// double rho_test = get_rho_2h_at_r_and_m_and_z(r_test,m_test,z_test,ptsz,pba);
// //
// printf("r_test = %.5e, z_test = %.5e, rho_test = %.12e\n",
//         r_test,z_test,rho_test);
// exit(0);

// double k_min = ptsz->k_min_samp_fftw;
// double k_max = ptsz->k_max_samp_fftw; // this is a precision parameter
// // tabulate the integrand in the "l" dimension:
// const int N = ptsz->N_samp_fftw;
//
//
// class_alloc(ptsz->array_profile_rho_2h_at_r_and_z,
//             N*ptsz->n_z_density_profile*sizeof(double),
//             ptsz->error_message);
// class_alloc(ptsz->array_profile_ln_r,
//             N*sizeof(double),
//             ptsz->error_message);
//
// int index_z;
// for (index_z=0; index_z<ptsz->n_z_density_profile; index_z++){
//   double z = exp(ptsz->array_profile_ln_1pz[index_z])-1.;
//   double k[N], Pk1[N];
//   int index_k;
//   for (index_k=0; index_k<N; index_k++)
//   {
//
//     k[index_k] = exp(log(k_min)+index_k/(N-1.)*(log(k_max)-log(k_min)));
//     Pk1[index_k] = get_rho_2h_at_k_and_z(k[index_k],z,ptsz);
//     Pk1[index_k] *= get_pk_lin_at_k_and_z(k[index_k],z,pba,ppm,pnl,ptsz);
//     printf("z = %.3e k = %.5e pk1 = %.5e\n",z,k[index_k],Pk1[index_k]);
//   }
//
//   double rp[N], xi1[N];
//   xi2pk(N,k,Pk1,rp,xi1,ptsz);
//   printf("\n##############\n");
//
//   for (index_k=0; index_k<N; index_k++){
//     int index_k_z = index_k * ptsz->n_z_density_profile + index_z;
//     ptsz->array_profile_rho_2h_at_r_and_z[index_k_z] = xi1[index_k];
//     ptsz->array_profile_ln_r[index_k] = log(rp[index_k]);
//     printf("z = %.3e r = %.5e xi1 = %.5e\n",z,rp[index_k],xi1[index_k]);
//    }
//
//
// }

// int index_k_z = index_k * n_z + index_z;
// exit(0);
}

if (ptsz->has_kSZ_kSZ_gal_1h_fft
 || ptsz->has_kSZ_kSZ_gal_2h_fft
 || ptsz->has_kSZ_kSZ_gal_3h_fft
 || ptsz->has_kSZ_kSZ_gal_3h
 || ptsz->has_kSZ_kSZ_gallens_1h_fft
 || ptsz->has_kSZ_kSZ_gallens_2h_fft
 || ptsz->has_kSZ_kSZ_gallens_3h_fft
 || ptsz->has_kSZ_kSZ_lens_1h_fft
 || ptsz->has_kSZ_kSZ_lens_2h_fft
 || ptsz->has_kSZ_kSZ_lens_3h_fft
){
tabulate_psi_b1t(pba,pnl,ppm,ppt,ptsz);
tabulate_psi_b2t(pba,pnl,ppm,ptsz);
}

if (ptsz->sz_verbose>0)
  printf("-> Starting main parallel block.\n");

   double * Pvecback;
   double * Pvectsz;
   // double * b_l1_l2_l_1d;
   int index_integrand;

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

//printf("number_of_threads = %d\n",number_of_threads);
// number_of_threads= 1;
int id;
omp_lock_t lock;


#pragma omp parallel \
   shared(abort,pba,ptsz,ppm,pnl,lock)\
   private(tstart,tstop,Pvectsz,Pvecback,index_integrand,id)\
   num_threads(number_of_threads)
	 {

#ifdef _OPENMP
	   tstart = omp_get_wtime();
#endif

	   class_alloc_parallel(Pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
       int i;
       for(i = 0; i<ptsz->tsz_size;i++) Pvectsz[i] = 0.;

	   class_alloc_parallel(Pvecback,pba->bg_size*sizeof(double),ptsz->error_message);


     // class_alloc_parallel(b_l1_l2_l_1d,
     //                       sizeof(double *)*ptsz->N_kSZ2_gal_theta_grid*ptsz->N_kSZ2_gal_multipole_grid,
     //                       ptsz->error_message);




//Loop over integrands
//the computation is parallelized with respect to the integrands

#pragma omp for schedule (dynamic)
for (index_integrand=0;index_integrand<ptsz->number_of_integrands;index_integrand++)
	     {
#pragma omp flush(abort)

       Pvectsz[ptsz->index_integrand_id] = index_integrand;

       // class_call_parallel(compute_sz(pba,
       //                                pnl,
       //                                ppm,
       //                                ppt,
       //                                ptsz,
       //                                Pvecback,
       //                                Pvectsz),
       //                               ptsz->error_message,
       //                               ptsz->error_message);
                          compute_sz(pba,
                                      pnl,
                                      ppm,
                                      ppt,
                                      ptsz,
                                      Pvecback,
                                      Pvectsz);



          }
#ifdef _OPENMP
      tstop = omp_get_wtime();
      if (ptsz->sz_verbose > 0)
         printf("In %s: time spent in parallel region (loop over class_sz integrals) = %e s for thread %d\n",
                   __func__,tstop-tstart,omp_get_thread_num());


#endif
   free(Pvecback);
   free(Pvectsz);
   // free(b_l1_l2_l_1d);
	} //end of parallel region

   if (abort == _TRUE_) return _FAILURE_;

   ////////////////end - cl
   // printf("showing res\n");
   // if (ptsz->sz_verbose>1) show_results(pba,pnl,ppm,ptsz);
     }
// printf("printing results\n");

// double z=1.;
// double m200m_pivot = 3e14*0.7;
//
// double r = get_m200m_to_m500c_at_z_and_M(z,m200m_pivot,ptsz)/m200m_pivot;
// printf("r = %.5e\n",r);
// exit(0);


if (ptsz->has_sz_rates){
  ptsz->szunbinned_loglike = - ptsz->hmf_int*4.*_PI_*ptsz->fsky_from_skyfracs;
  int index_rate;
  for (index_rate=0;index_rate<ptsz->szcat_size;index_rate++){
    ptsz->szunbinned_loglike += log(ptsz->szrate[index_rate]);
  // printf("cluster id = %d\trate = %.4e \n",index_rate,ptsz->szrate[index_rate]);
  }
  if (ptsz->sz_verbose >= 1 ) printf("loglike unbinned cc = %.3e\n",ptsz->szunbinned_loglike);
  if (ptsz->sz_verbose >= 1 ) printf("Ntot = %.3e\n",ptsz->hmf_int*4.*_PI_*ptsz->fsky_from_skyfracs);
}


  if ( (ptsz->has_gal_gallens_1h || ptsz->has_gal_gallens_2h) && ptsz->convert_cls_to_gamma){
    // shear
    if (ptsz->sz_verbose > 0) printf("converting cls to gamma\n");
    class_alloc(ptsz->thetas_arcmin,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
    class_alloc(ptsz->gamma_gal_gallens_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
    class_alloc(ptsz->gamma_gal_gallens_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);


    int i;
    double lnl[ptsz->nlSZ],lncl_1h[ptsz->nlSZ],lncl_2h[ptsz->nlSZ];
    for (i=0;i<ptsz->nlSZ;i++){
      double fac = ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_);
      lnl[i] = log(ptsz->ell[i]);
      lncl_1h[i] = log(ptsz->cl_gal_gallens_1h[i]/fac);
      lncl_2h[i] = log(ptsz->cl_gal_gallens_2h[i]/fac);

    }

    const int N = ptsz->N_samp_fftw;
    double l[N],thetas[N], gamma_t_1h[N], gamma_t_2h[N], cl_1h[N], cl_2h[N], l_min, l_max;
    l_min = ptsz->l_min_samp_fftw;
    l_max = ptsz->l_max_samp_fftw;

    for (i=0;i<N;i++){
    l[i] = exp(log(l_min)+i/(N-1.)*(log(l_max)-log(l_min)));

    if ((l[i]<ptsz->ell[0]) || (l[i]>ptsz->ell[ptsz->nlSZ-1])){
      cl_1h[i] = 0.;
      cl_2h[i] = 0.;
    }
    else{
    cl_1h[i] = exp(pwl_value_1d(ptsz->nlSZ,lnl,lncl_1h,log(l[i])));
    cl_2h[i] = exp(pwl_value_1d(ptsz->nlSZ,lnl,lncl_2h,log(l[i])));
  }
    }

    cl2gamma(N,l,cl_1h,thetas,gamma_t_1h,ptsz);
    cl2gamma(N,l,cl_2h,thetas,gamma_t_2h,ptsz);
    // double ltest[N],cltest[N];
    // gamma2cl(N,thetas,gamma_t_1h,ltest,cltest,ptsz);
    // for (i=0;i<N;i++){
    //   printf("%.5e %.5e\n",cltest[i],cl_1h[i]);
    // }
    // gamma2cl(N,l,cl_2h,thetas,gamma_t_2h,ptsz);
    for (i=0;i<ptsz->nlSZ;i++){
      ptsz->thetas_arcmin[i] = 1./ptsz->ell[ptsz->nlSZ-i-1]*(60.*180.)/_PI_;
      ptsz->gamma_gal_gallens_1h[i] = pwl_value_1d(N,thetas,gamma_t_1h,1./ptsz->ell[ptsz->nlSZ-i-1]);
      ptsz->gamma_gal_gallens_2h[i] = pwl_value_1d(N,thetas,gamma_t_2h,1./ptsz->ell[ptsz->nlSZ-i-1]);
if (ptsz->sz_verbose > 0){
      printf("thetas = %.5e gamma_t_1h = %.5e gamma_t_2h = %.5e\n",ptsz->thetas_arcmin[i],ptsz->gamma_gal_gallens_1h[i],ptsz->gamma_gal_gallens_2h[i]);
    }
  }
  }

if (ptsz->need_ksz_template){
   load_cl_ksz_template(ptsz);
 }

 if (ptsz->need_tt_noise){
   load_unbinned_nl_tt(ptsz);
 }

 if (ptsz->need_lensing_noise){
   load_nl_lensing_noise(ptsz);
 }


if (ptsz->has_kSZ_kSZ_gal_covmat
||  ptsz->has_kSZ_kSZ_gallens_covmat
||  ptsz->has_kSZ_kSZ_lens_covmat){

///////////////////
//tabulate cl_ttf//
///////////////////
double Delta_t2 = 0.; // in radians
double theta_fwhm = 1.; // in radians
//collect lensed cl's:
double l_min = ptsz->l_min_samp_fftw;
double l_max = ptsz->l_max_samp_fftw; // this is a precision parameter
const int N = ptsz->N_samp_fftw;
double l[N];
double cl_tt_lensed_fft[N];
double cl_ksz_fft[N];
double cl_noise_fft[N];
double Fl_bl;
double cl_ttf[N];
int i;
double * cl_tot;
class_alloc(cl_tot,
            psp->ct_size*sizeof(double),
            ptsz->error_message);

double cl_tt_lensed[ple->l_lensed_max-2];
double l_lensed[ple->l_lensed_max-2];
// printf("ct_size = %d\n",psp->ct_size);
// printf("l_max lensed %d\n",ple->l_lensed_max);
// exit(0);
for (i=2;i<ple->l_lensed_max;i++){
  lensing_cl_at_l(ple,i,cl_tot);
  l_lensed[i-2] = i;
  cl_tt_lensed[i-2] = cl_tot[0]; //dimensionless total lensed cl TT
//   if ((i<20) || (i>(ple->l_lensed_max-10))){
//   printf("i %d l_lensed[i] %e dl_lensed[i] %.8e cl_lensed[i] %.8e\n",i,l_lensed[i-2],l_lensed[i-2]*(l_lensed[i-2]+1.)/2./_PI_*cl_tt_lensed[i-2],cl_tt_lensed[i-2]);
// }

}
free(cl_tot);

// load_cl_ksz_template(ptsz);
// load_unbinned_nl_tt(ptsz);
// exit(0);
// printf("%.3e %.8e\n",ptsz->unbinned_nl_tt_ell[0],ptsz->unbinned_nl_tt_n_ell[0]);
// printf("%.3e %.8e\n",ptsz->unbinned_nl_tt_ell[1],ptsz->unbinned_nl_tt_n_ell[1]);
// exit(0);
for (i=0;i<N;i++){
l[i] = exp(log(l_min)+i/(N-1.)*(log(l_max)-log(l_min)));
cl_tt_lensed_fft[i] = pwl_value_1d(ple->l_lensed_max-2,l_lensed,cl_tt_lensed,l[i]);
cl_ksz_fft[i] = pwl_value_1d(ptsz->ksz_template_size,ptsz->l_ksz_template,ptsz->cl_ksz_template,l[i]);

if (ptsz->no_tt_noise_in_kSZ2X_cov){
cl_noise_fft[i] = 1e-100;
}
else{
cl_noise_fft[i] = pwl_value_1d(ptsz->unbinned_nl_tt_size,ptsz->unbinned_nl_tt_ell,ptsz->unbinned_nl_tt_n_ell,l[i]);
}

cl_ksz_fft[i] = cl_ksz_fft[i]/l[i]/(l[i]+1.)*2.*_PI_*1e-12/pba->T_cmb/pba->T_cmb;
cl_noise_fft[i] = cl_noise_fft[i]*1e-12/pba->T_cmb/pba->T_cmb;
Fl_bl = get_ksz_filter_at_l(l[i],ptsz);


double result_ttf;
if (ptsz->compute_ksz2ksz2==1){
result_ttf = Fl_bl*Fl_bl*(cl_ksz_fft[i])*sqrt(2./(2.*_PI_)/(2.*_PI_)); // divide by 2pi^2 and multiply by 2 -- see formula

}
else{
result_ttf = Fl_bl*Fl_bl*(cl_tt_lensed_fft[i]+cl_ksz_fft[i]+cl_noise_fft[i])*sqrt(2./(2.*_PI_)/(2.*_PI_)); // divide by 2pi^2 and multiply by 2 -- see formula
}

if (l[i]<l_lensed[0]){
  result_ttf =0.;
  cl_tt_lensed_fft[i] = 0.;}
if (l[i]>l_lensed[ple->l_lensed_max-3]){
  result_ttf =0.;
  cl_tt_lensed_fft[i] = 0.;}

if (l[i]<ptsz->l_ksz_template[0]){
  result_ttf =0.;
  cl_ksz_fft[i] = 0.;}
if (l[i]>ptsz->l_ksz_template[ptsz->ksz_template_size-1]){
  result_ttf =0.;
  cl_ksz_fft[i] = 0.;}

if (l[i]<ptsz->unbinned_nl_tt_ell[0]){
  result_ttf =0.;
  cl_noise_fft[i] = 0.;}
if (l[i]>ptsz->unbinned_nl_tt_ell[ptsz->unbinned_nl_tt_size-1]){
  result_ttf =0.;
  cl_noise_fft[i] = 0.;}


cl_ttf[i] = result_ttf;

// cl_noise_fft[i] = Delta_t2*exp(theta_fwhm*theta_fwhm*l[i]*l[i]/8./log(2.));

// if (i<10)
// printf(" l = %.3e cl tt = %.8e  ksz = %.8e flbl = %.8e ttf = %.8e\n",l[i],cl_tt_lensed_fft[i],cl_ksz_fft[i],Fl_bl,cl_ttf[i]);
// lensing_cl_at_l(ple,l[i])
}

double r[N], xi[N], cl_t2t2f[N];
xi2pk(N,l,cl_ttf,r,xi,ptsz);
for (i=0; i<N; i++){
// convolution:
xi[i] = pow(xi[i],2.);
}
pk2xi(N,r,xi,l,cl_t2t2f,ptsz);

double clp_t2t2f;


for (i=0;i<ptsz->nlSZ;i++){

clp_t2t2f = pwl_value_1d(N,l,cl_t2t2f,ptsz->ell[i]);
ptsz->cl_t2t2f[i] = clp_t2t2f;


if (ptsz->has_kSZ_kSZ_gal_covmat){
  double cl_gg,nl_gg,cl_kSZ2g;

  if ((ptsz->has_gal_gal_hf==1) && ((ptsz->has_gal_gal_2h+ptsz->has_gal_gal_1h)==0)){
  cl_gg = (ptsz->cl_gal_gal_hf[i])
          /(ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_));
  }
  else{
  cl_gg = (ptsz->cl_gal_gal_2h[i] + ptsz->cl_gal_gal_1h[i])
          /(ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_));
  }

  nl_gg = ptsz->cl_gal_gal_A_sn; //galaxy shot-noise

  if ((ptsz->has_kSZ_kSZ_gal_hf==1) && ((ptsz->has_kSZ_kSZ_gal_3h+ptsz->has_kSZ_kSZ_gal_2h+ptsz->has_kSZ_kSZ_gal_1h)==0)){
  cl_kSZ2g = ptsz->cl_kSZ_kSZ_gal_hf[i];
    }
  else{
    cl_kSZ2g = ptsz->cl_kSZ_kSZ_gal_1h_fft[i]
              +ptsz->cl_kSZ_kSZ_gal_2h_fft[i]
              +ptsz->cl_kSZ_kSZ_gal_3h_fft[i];
  }

  ptsz->cov_ll_kSZ_kSZ_gal[i] = (clp_t2t2f*(cl_gg+nl_gg)+cl_kSZ2g*cl_kSZ2g)
                                *1./(2.*ptsz->ell[i]+1.)
                                /ptsz->f_sky;

}

if (ptsz->has_kSZ_kSZ_gallens_covmat){
  double cl_kgkg,nl_kgkg,cl_kSZ2kg;

  cl_kgkg = (ptsz->cl_gallens_gallens_2h[i] + ptsz->cl_gallens_gallens_1h[i])
          /(ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_));


  double arcmin_to_radians = _PI_/(60.*180.);
  nl_kgkg = ptsz->shape_noise_siggamma2/(ptsz->ns_gal_per_arcmin2/arcmin_to_radians/arcmin_to_radians);

  cl_kSZ2kg = ptsz->cl_kSZ_kSZ_gallens_1h_fft[i]
            +ptsz->cl_kSZ_kSZ_gallens_2h_fft[i]
            +ptsz->cl_kSZ_kSZ_gallens_3h_fft[i];
  ptsz->cov_ll_kSZ_kSZ_gallens[i] = (clp_t2t2f*(cl_kgkg+nl_kgkg)+cl_kSZ2kg*cl_kSZ2kg)
                                *1./(2.*ptsz->ell[i]+1.)
                                /ptsz->f_sky;
}

if (ptsz->has_kSZ_kSZ_lens_covmat){

  double cl_kgkg,nl_kcmb_kcmb,cl_kSZ2kg;

  cl_kgkg = (ptsz->cl_lens_lens_2h[i] + ptsz->cl_lens_lens_1h[i])
          /(ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_));

  cl_kSZ2kg = ptsz->cl_kSZ_kSZ_lens_1h_fft[i]
            +ptsz->cl_kSZ_kSZ_lens_2h_fft[i]
            +ptsz->cl_kSZ_kSZ_lens_3h_fft[i];

  nl_kcmb_kcmb = get_lensing_noise_at_ell(ptsz->ell[i],ptsz);
  // printf("i=%d l=%.5e nli=%.5e\n",i,ptsz->ell[i],nl_kcmb_kcmb);

  ptsz->cov_ll_kSZ_kSZ_lens[i] = (clp_t2t2f*(cl_kgkg+nl_kcmb_kcmb)+cl_kSZ2kg*cl_kSZ2kg)
                                *1./(2.*ptsz->ell[i]+1.)
                                /ptsz->f_sky;
}


double ln_ell_up,ln_ell_down,ln_ell_max,ln_ell_min,n_modes;
double ell = ptsz->ell[i];
if (ptsz->dell == 0){
 if (i == 0){
    ln_ell_up = log(ptsz->ell[i+1]);
    ln_ell_max = log(ell) + 0.5*(ln_ell_up-log(ell));
    ln_ell_min = log(ell) - 0.5*(ln_ell_up-log(ell));
    n_modes = exp(ln_ell_max)-exp(ln_ell_min);
 }
 else if (i == ptsz->nlSZ -1){
    ln_ell_down = log(ptsz->ell[i-1]);
    ln_ell_min = log(ell) - 0.5*(log(ell)-ln_ell_down);
    ln_ell_max = log(ell) + 0.5*(log(ell)-ln_ell_down);
    n_modes = exp(ln_ell_max)-exp(ln_ell_min);
 }
 else {
    ln_ell_down = log(ptsz->ell[i-1]);
    ln_ell_up = log(ptsz->ell[i+1]);
    ln_ell_min = log(ell) - 0.5*(log(ell)-ln_ell_down);
    ln_ell_max = log(ell) + 0.5*(ln_ell_up-log(ell));
    n_modes = exp(ln_ell_max)-exp(ln_ell_min);
 }
}
else{
 if (i == 0){
    ln_ell_up = ptsz->ell[i+1];
    ln_ell_max = ell + 0.5*(ln_ell_up-ell);
    ln_ell_min = ell - 0.5*(ln_ell_up-ell);
    n_modes = ln_ell_max-ln_ell_min;
 }
 else if (i == ptsz->nlSZ -1){
    ln_ell_down = ptsz->ell[i-1];
    ln_ell_min = ell - 0.5*(ell-ln_ell_down);
    ln_ell_max = ell + 0.5*(ell-ln_ell_down);
    n_modes = ln_ell_max-ln_ell_min;
 }
 else {
    ln_ell_down = ptsz->ell[i-1];
    ln_ell_up = ptsz->ell[i+1];
    ln_ell_min = ell - 0.5*(ell-ln_ell_down);
    ln_ell_max = ell + 0.5*(ln_ell_up-ell);
    n_modes = ln_ell_max-ln_ell_min;
 }
}
if (ptsz->has_kSZ_kSZ_gal_covmat){
ptsz->cov_ll_kSZ_kSZ_gal[i] *= 1./n_modes;
}
if (ptsz->has_kSZ_kSZ_gallens_covmat){
ptsz->cov_ll_kSZ_kSZ_gallens[i] *= 1./n_modes;
}
if (ptsz->has_kSZ_kSZ_lens_covmat){
ptsz->cov_ll_kSZ_kSZ_lens[i] *= 1./n_modes;
}



}

}


//
// if (ptsz->has_kSZ_kSZ_gal_covmat){
// if (ptsz->sz_verbose > 0){
// printf("starting TTG covmat computation.\n");
// }
//
// ///////////////////
// //tabulate cl_ttf//
// ///////////////////
// double Delta_t2 = 0.; // in radians
// double theta_fwhm = 1.; // in radians
// //collect lensed cl's:
// double l_min = ptsz->l_min_samp_fftw;
// double l_max = ptsz->l_max_samp_fftw; // this is a precision parameter
// const int N = ptsz->N_samp_fftw;
// // double l[N];
// // double cl_tt_lensed_fft[N];
// // double cl_ksz_fft[N];
// // double cl_noise_fft[N];
// // double Fl_bl;
// // double cl_ttf[N];
// int i;
// double * cl_tot;
// class_alloc(cl_tot,
//             psp->ct_size*sizeof(double),
//             ptsz->error_message);
// //
// // double cl_tt_lensed[ple->l_lensed_max-2];
// // double l_lensed[ple->l_lensed_max-2];
// // // printf("ct_size = %d\n",psp->ct_size);
// // // printf("l_max lensed %d\n",ple->l_lensed_max);
// // // exit(0);
// // for (i=2;i<ple->l_lensed_max;i++){
// //   lensing_cl_at_l(ple,i,cl_tot);
// //   l_lensed[i-2] = i;
// //   cl_tt_lensed[i-2] = cl_tot[0]; //dimensionless total lensed cl TT
// // //   if ((i<20) || (i>(ple->l_lensed_max-10))){
// // //   printf("i %d l_lensed[i] %e dl_lensed[i] %.8e cl_lensed[i] %.8e\n",i,l_lensed[i-2],l_lensed[i-2]*(l_lensed[i-2]+1.)/2./_PI_*cl_tt_lensed[i-2],cl_tt_lensed[i-2]);
// // // }
// //
// // }
// //
// // // load_cl_ksz_template(ptsz);
// // // load_unbinned_nl_tt(ptsz);
// // // exit(0);
// // // printf("%.3e %.8e\n",ptsz->unbinned_nl_tt_ell[0],ptsz->unbinned_nl_tt_n_ell[0]);
// // // printf("%.3e %.8e\n",ptsz->unbinned_nl_tt_ell[1],ptsz->unbinned_nl_tt_n_ell[1]);
// // // exit(0);
// // for (i=0;i<N;i++){
// // l[i] = exp(log(l_min)+i/(N-1.)*(log(l_max)-log(l_min)));
// // cl_tt_lensed_fft[i] = pwl_value_1d(ple->l_lensed_max-2,l_lensed,cl_tt_lensed,l[i]);
// // cl_ksz_fft[i] = pwl_value_1d(ptsz->ksz_template_size,ptsz->l_ksz_template,ptsz->cl_ksz_template,l[i]);
// // cl_noise_fft[i] = pwl_value_1d(ptsz->unbinned_nl_tt_size,ptsz->unbinned_nl_tt_ell,ptsz->unbinned_nl_tt_n_ell,l[i]);
// //
// // cl_ksz_fft[i] = cl_ksz_fft[i]/l[i]/(l[i]+1.)*2.*_PI_*1e-12/pba->T_cmb/pba->T_cmb;
// // cl_noise_fft[i] = cl_noise_fft[i]*1e-12/pba->T_cmb/pba->T_cmb;
// // Fl_bl = get_ksz_filter_at_l(l[i],ptsz);
// //
// //
// // if (l[i]<l_lensed[0])
// //   cl_tt_lensed_fft[i] = 0.;
// // if (l[i]>l_lensed[ple->l_lensed_max-3])
// //   cl_tt_lensed_fft[i] = 0.;
// //
// // if (l[i]<ptsz->l_ksz_template[0])
// //   cl_ksz_fft[i] = 0.;
// // if (l[i]>ptsz->l_ksz_template[ptsz->ksz_template_size-1])
// //   cl_ksz_fft[i] = 0.;
// //
// // if (l[i]<ptsz->unbinned_nl_tt_ell[0])
// //   cl_noise_fft[i] = 0.;
// // if (l[i]>ptsz->unbinned_nl_tt_ell[ptsz->unbinned_nl_tt_size-1])
// //   cl_noise_fft[i] = 0.;
// //
// //
// // // cl_noise_fft[i] = Delta_t2*exp(theta_fwhm*theta_fwhm*l[i]*l[i]/8./log(2.));
// //
// // cl_ttf[i] = Fl_bl*Fl_bl*(cl_tt_lensed_fft[i]+cl_ksz_fft[i]+cl_noise_fft[i])*sqrt(2./(2.*_PI_)/(2.*_PI_)); // divide by 2pi^2 and multiply by 2 -- see formula
// // // if (i<10)
// // // printf(" l = %.3e cl tt = %.8e  ksz = %.8e flbl = %.8e ttf = %.8e\n",l[i],cl_tt_lensed_fft[i],cl_ksz_fft[i],Fl_bl,cl_ttf[i]);
// // // lensing_cl_at_l(ple,l[i])
// // }
// //
// // double r[N], xi[N], cl_t2t2f[N];
// // xi2pk(N,l,cl_ttf,r,xi,ptsz);
// // for (i=0; i<N; i++){
// // // convolution:
// // xi[i] = pow(xi[i],2.);
// // }
// // pk2xi(N,r,xi,l,cl_t2t2f,ptsz);
//
// // for (i=0;i<N;i++){
// // printf("l = %.3e t2t2f = %.5e\n",l[i],cl_t2t2f[i]);
// // }
// // double cl_gg,clp_t2t2f,cl_kSZ2g;
// // for (i=0;i<ptsz->nlSZ;i++){
// // cl_gg = (ptsz->cl_gal_gal_2h[i] + ptsz->cl_gal_gal_1h[i])
// //         /(ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_));
// // clp_t2t2f = pwl_value_1d(N,l,cl_t2t2f,ptsz->ell[i]);
//
// cl_kSZ2g = ptsz->cl_kSZ_kSZ_gal_1h_fft[i]
//           +ptsz->cl_kSZ_kSZ_gal_2h_fft[i]
//           +ptsz->cl_kSZ_kSZ_gal_3h_fft[i];
// ptsz->cov_ll_kSZ_kSZ_gal[i] = (clp_t2t2f*cl_gg+cl_kSZ2g*cl_kSZ2g)
//                               *1./(2.*ptsz->ell[i]+1.)
//                               /ptsz->f_sky;
// double ln_ell_up,ln_ell_down,ln_ell_max,ln_ell_min,n_modes;
// double ell = ptsz->ell[i];
// if (ptsz->dell == 0){
//  if (i == 0){
//     ln_ell_up = log(ptsz->ell[i+1]);
//     ln_ell_max = log(ell) + 0.5*(ln_ell_up-log(ell));
//     ln_ell_min = log(ell) - 0.5*(ln_ell_up-log(ell));
//     n_modes = exp(ln_ell_max)-exp(ln_ell_min);
//  }
//  else if (i == ptsz->nlSZ -1){
//     ln_ell_down = log(ptsz->ell[i-1]);
//     ln_ell_min = log(ell) - 0.5*(log(ell)-ln_ell_down);
//     ln_ell_max = log(ell) + 0.5*(log(ell)-ln_ell_down);
//     n_modes = exp(ln_ell_max)-exp(ln_ell_min);
//  }
//  else {
//     ln_ell_down = log(ptsz->ell[i-1]);
//     ln_ell_up = log(ptsz->ell[i+1]);
//     ln_ell_min = log(ell) - 0.5*(log(ell)-ln_ell_down);
//     ln_ell_max = log(ell) + 0.5*(ln_ell_up-log(ell));
//     n_modes = exp(ln_ell_max)-exp(ln_ell_min);
//  }
// }
// else{
//  if (i == 0){
//     ln_ell_up = ptsz->ell[i+1];
//     ln_ell_max = ell + 0.5*(ln_ell_up-ell);
//     ln_ell_min = ell - 0.5*(ln_ell_up-ell);
//     n_modes = ln_ell_max-ln_ell_min;
//  }
//  else if (i == ptsz->nlSZ -1){
//     ln_ell_down = ptsz->ell[i-1];
//     ln_ell_min = ell - 0.5*(ell-ln_ell_down);
//     ln_ell_max = ell + 0.5*(ell-ln_ell_down);
//     n_modes = ln_ell_max-ln_ell_min;
//  }
//  else {
//     ln_ell_down = ptsz->ell[i-1];
//     ln_ell_up = ptsz->ell[i+1];
//     ln_ell_min = ell - 0.5*(ell-ln_ell_down);
//     ln_ell_max = ell + 0.5*(ln_ell_up-ell);
//     n_modes = ln_ell_max-ln_ell_min;
//  }
// }
// ptsz->cov_ll_kSZ_kSZ_gal[i] *= 1./n_modes;


// printf("l = %.3e t2t2f = %.5e gg = %.5e t2t2gg = %.5e t2gt2g = %.5e  cov = %.5e\n",ptsz->ell[i],clp_t2t2f,cl_gg,clp_t2t2f*cl_gg,cl_kSZ2g*cl_kSZ2g,ptsz->cov_ll_kSZ_kSZ_gal[i]);
// }

// free(cl_tot);
// }



if (ptsz->has_kSZ_kSZ_gal_lensing_term
||  ptsz->has_kSZ_kSZ_gallens_lensing_term
||  ptsz->has_kSZ_kSZ_lens_lensing_term){
if (ptsz->sz_verbose > 0){
printf("starting TTG lensing term computation.\n");
}
int i;

//collect unlensed cls at lprime:
// spectra_cl_at_l(psp,(double)l,cl,cl_md,cl_md_ic);
double cl_tt_unlensed[ppt->l_scalar_max-2];
double l_unlensed[ppt->l_scalar_max-2];
double * cl_tot;
class_alloc(cl_tot,
            psp->ct_size*sizeof(double),
            ptsz->error_message);
  double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*psp->ct_size+index_ct] */

  double ** cl_md;    /* array with argument
                         cl_md[index_md][index_ct] */
  int index_md;

    class_alloc(cl_md_ic,
                psp->md_size*sizeof(double *),
                ptsz->error_message);

    class_alloc(cl_md,
                psp->md_size*sizeof(double *),
                ptsz->error_message);
    for (index_md = 0; index_md < psp->md_size; index_md++) {

      if (psp->md_size > 1)

        class_alloc(cl_md[index_md],
                    psp->ct_size*sizeof(double),
                    ptsz->error_message);

      if (psp->ic_size[index_md] > 1)

        class_alloc(cl_md_ic[index_md],
                    psp->ic_ic_size[index_md]*psp->ct_size*sizeof(double),
                    ptsz->error_message);
    }


for (i=2;i<ppt->l_scalar_max;i++){
    l_unlensed[i-2] = i;

    class_call(spectra_cl_at_l(psp,
                               i,
                               cl_tot,
                               cl_md,
                               cl_md_ic),
               psp->error_message,
               ptsz->error_message);
    cl_tt_unlensed[i-2] = cl_tot[0];
//   if ((i<20) || (i>(ppt->l_scalar_max-10))){
//   printf("i %d l_unlensed[i] %e dl_unlensed[i] %.8e cl_unlensed[i] %.8e\n",i,l_unlensed[i-2],l_unlensed[i-2]*(l_unlensed[i-2]+1.)/2./_PI_*cl_tt_unlensed[i-2],cl_tt_unlensed[i-2]);
// }
}

free(cl_tot);
  for (index_md = 0; index_md < psp->md_size; index_md++) {
    if (psp->md_size > 1)
      free(cl_md[index_md]);
    if (psp->ic_size[index_md] > 1)
       free(cl_md_ic[index_md]);
}
free(cl_md);
free(cl_md_ic);

// exit(0);
for (i=0;i<ptsz->nlSZ;i++){
// at each l, we tabulate the integrand:
double * integrand_l_lprime_phi;
class_alloc(integrand_l_lprime_phi,
            sizeof(double *)*ptsz->N_kSZ2_gal_theta_grid*ptsz->N_kSZ2_gal_multipole_grid,
            ptsz->error_message);
int index_lprime = 0;
int index_theta = 0;
int index_lprime_theta = 0;
for (index_lprime=0;index_lprime<ptsz->N_kSZ2_gal_multipole_grid;index_lprime++){
for (index_theta=0;index_theta<ptsz->N_kSZ2_gal_theta_grid;index_theta++){
double lprime = exp(ptsz->ell_kSZ2_gal_multipole_grid[index_lprime]);
double theta = ptsz->theta_kSZ2_gal_theta_grid[index_theta];
double l = ptsz->ell[i];

double flprime = get_ksz_filter_at_l(lprime,ptsz);
double abs_l_plus_lprime = sqrt(l*l+lprime*lprime+2.*l*lprime*cos(theta));
double fl_plus_lprime = get_ksz_filter_at_l(abs_l_plus_lprime,ptsz);

double clprime_tt_unlensed = pwl_value_1d(ppt->l_scalar_max-2,l_unlensed,cl_tt_unlensed,lprime);
if (lprime<l_unlensed[0])
  clprime_tt_unlensed = 0.;
if (lprime>l_unlensed[ppt->l_scalar_max-3])
  clprime_tt_unlensed = 0.;


integrand_l_lprime_phi[index_lprime_theta] = lprime*lprime*flprime*clprime_tt_unlensed*fl_plus_lprime*cos(theta);
integrand_l_lprime_phi[index_lprime_theta] *= lprime; //for integration in log(ell)
// printf("l = %.4e lp = %.4e th = %.3e flprime = %.5e fl_plus_lprime = %.5e clprime = %.5e integrand_l_lprime_phi = %.5e\n",
// l,lprime,theta,flprime,fl_plus_lprime,clprime_tt_unlensed,
// integrand_l_lprime_phi[index_lprime_theta]);
index_lprime_theta += 1;
}
}

// tabulation done
   double epsrel= 1.e-6;//ptsz->redshift_epsrel;//ptsz->patterson_epsrel;
   double epsabs= 1.e-50;//ptsz->redshift_epsabs;//ptsz->patterson_epsabs;
   int show_neval = 0;//ptsz->patterson_show_neval;

   struct Parameters_for_integrand_kSZ2_X_lensing_term V;
   V.pnl= pnl;
   V.ppm= ppm;
   V.ptsz = ptsz;
   V.pba = pba;
   // V.Pvecback = Pvecback;
   // V.Pvectsz = Pvectsz;
   V.ln_ellprime = ptsz->ell_kSZ2_gal_multipole_grid;//ln_ell;
   V.index_ell = i;
   V.integrand_l_lprime_phi = integrand_l_lprime_phi;
   void * params;
   params = &V;

double r_lens = Integrate_using_Patterson_adaptive(0., 2.*_PI_,
                                            epsrel, epsabs,
                                            integrand_kSZ2_X_lensing_term,
                                            params,show_neval);
// printf("r_lens = %.5e\n",r_lens);

if (ptsz->has_kSZ_kSZ_gal_lensing_term){
ptsz->cl_kSZ_kSZ_gal_lensing_term[i] = r_lens;

double cl_gk;// = (ptsz->cl_gal_lens_2h[i] + ptsz->cl_gal_lens_1h[i])
              //  /(ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_));

if ((ptsz->has_gal_lens_1h == _FALSE_) && (ptsz->has_gal_lens_2h == _FALSE_)){
  cl_gk = (ptsz->cl_gal_lens_hf[i])
                  /(ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_));
}
else{
  cl_gk = (ptsz->cl_gal_lens_2h[i] + ptsz->cl_gal_lens_1h[i])
                  /(ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_));

}

double cl_gp = 2./(ptsz->ell[i]*(ptsz->ell[i]+1.))*cl_gk;
ptsz->cl_kSZ_kSZ_gal_lensing_term[i] *= -2./(2.*_PI_)/(2.*_PI_)*ptsz->ell[i]*cl_gp;
}

if (ptsz->has_kSZ_kSZ_gallens_lensing_term){
ptsz->cl_kSZ_kSZ_gallens_lensing_term[i] = r_lens;
double cl_gk = (ptsz->cl_gallens_lens_2h[i] + ptsz->cl_gallens_lens_1h[i])
                /(ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_));
double cl_gp = 2./(ptsz->ell[i]*(ptsz->ell[i]+1.))*cl_gk;
ptsz->cl_kSZ_kSZ_gallens_lensing_term[i] *= -2./(2.*_PI_)/(2.*_PI_)*ptsz->ell[i]*cl_gp;
}

if (ptsz->has_kSZ_kSZ_lens_lensing_term){
ptsz->cl_kSZ_kSZ_lens_lensing_term[i] = r_lens;
// double cl_gk = (ptsz->cl_lens_lens_2h[i] + ptsz->cl_lens_lens_1h[i])
//                 /(ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_));
double cl_kk;


if ((ptsz->has_lens_lens_1h == _FALSE_) && (ptsz->has_lens_lens_2h == _FALSE_)){
  cl_kk = (ptsz->cl_lens_lens_hf[i])
                  /(ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_));
}
else{
  cl_kk = (ptsz->cl_lens_lens_2h[i] + ptsz->cl_lens_lens_1h[i])
                  /(ptsz->ell[i]*(ptsz->ell[i]+1.)/(2*_PI_));

}



double cl_kp = pow(2./(ptsz->ell[i]*(ptsz->ell[i]+1.)),1.)*cl_kk;
// printf("cl_kSZ_kSZ_lens_lensing_term[i] = %.5e cl_kp = %.5e\n",ptsz->cl_kSZ_kSZ_lens_lensing_term[i],cl_kp);
ptsz->cl_kSZ_kSZ_lens_lensing_term[i] *= -2./(2.*_PI_)/(2.*_PI_)*ptsz->ell[i]*cl_kp;
}


free(integrand_l_lprime_phi);
}

}

   if (ptsz->sz_verbose>1) show_results(pba,pnl,ppm,ptsz);
   if (ptsz->write_sz>0 || ptsz->create_ref_trispectrum_for_cobaya){

   write_output_to_files_cl(pnl,pba,ppm,ptsz);
   write_output_to_files_ell_indep_ints(pnl,pba,ptsz);
   write_redshift_dependent_quantities(pba,ptsz);
   if (ptsz->sz_verbose>1)
      printf("done writing things.\n");
 }
   return _SUCCESS_;
}




int szpowerspectrum_free(struct tszspectrum *ptsz)
{
  if (ptsz->sz_verbose>1) printf("-> freeing memory.\n");

    int all_comps = ptsz->has_sz_ps
      + ptsz->has_hmf
      + ptsz->has_pk_at_z_1h
      + ptsz->has_pk_at_z_2h
      + ptsz->has_pk_gg_at_z_1h
      + ptsz->has_pk_gg_at_z_2h
      + ptsz->has_pk_bb_at_z_1h
      + ptsz->has_pk_bb_at_z_2h
      + ptsz->has_pk_b_at_z_2h
      + ptsz->has_gas_density_profile_2h
      + ptsz->has_gas_pressure_profile_2h
      + ptsz->has_pk_em_at_z_1h
      + ptsz->has_pk_em_at_z_2h
      + ptsz->has_pk_HI_at_z_1h
      + ptsz->has_pk_HI_at_z_2h
      + ptsz->has_bk_at_z_1h
      + ptsz->has_bk_at_z_2h
      + ptsz->has_bk_at_z_3h
      + ptsz->has_bk_ttg_at_z_1h
      + ptsz->has_bk_ttg_at_z_2h
      + ptsz->has_bk_ttg_at_z_3h
      + ptsz->has_mean_y
      + ptsz->has_cib_monopole
      + ptsz->has_cib_shotnoise
      + ptsz->has_dcib0dz
      + ptsz->has_dydz
      + ptsz->has_sz_2halo
      + ptsz->has_sz_trispec
      + ptsz->has_sz_m_y_y_1h
      + ptsz->has_sz_m_y_y_2h
      + ptsz->has_sz_te_y_y
      + ptsz->has_sz_cov_N_N
      + ptsz->has_tSZ_tSZ_tSZ_1halo
      + ptsz->has_tSZ_tSZ_tSZ_2h
      + ptsz->has_tSZ_tSZ_tSZ_3h
      + ptsz->has_kSZ_kSZ_1h
      + ptsz->has_kSZ_kSZ_2h
      + ptsz->has_kSZ_kSZ_tSZ_1h
      + ptsz->has_kSZ_kSZ_tSZ_2h
      + ptsz->has_kSZ_kSZ_tSZ_3h
      + ptsz->has_kSZ_kSZ_gal_1h
      + ptsz->has_kSZ_kSZ_gal_1h_fft
      + ptsz->has_kSZ_kSZ_gal_2h_fft
      + ptsz->has_kSZ_kSZ_gal_3h_fft
      + ptsz->has_kSZ_kSZ_gal_2h
      + ptsz->has_kSZ_kSZ_gal_3h
      + ptsz->has_kSZ_kSZ_gal_hf
      + ptsz->has_kSZ_kSZ_lensmag_1halo
      + ptsz->has_tSZ_gal_1h
      + ptsz->has_tSZ_gal_2h
      + ptsz->has_tSZ_lensmag_1h
      + ptsz->has_tSZ_lensmag_2h
      + ptsz->has_tSZ_cib_1h
      + ptsz->has_tSZ_cib_2h
      + ptsz->has_lens_cib_1h
      + ptsz->has_lens_cib_2h
      + ptsz->has_gal_cib_1h
      + ptsz->has_gal_cib_2h
      + ptsz->has_ngal_ngal_1h
      + ptsz->has_ngal_ngal_2h
      + ptsz->has_ngal_ngal_hf
      + ptsz->has_ngal_lens_1h
      + ptsz->has_ngal_lens_2h
      + ptsz->has_ngal_lens_hf
      + ptsz->has_cib_cib_1h
      + ptsz->has_cib_cib_2h
      + ptsz->has_gal_gal_1h
      + ptsz->has_gal_gal_2h
      + ptsz->has_gal_gal_hf
      + ptsz->has_gal_lens_1h
      + ptsz->has_gal_lens_2h
      + ptsz->has_gal_lens_hf
      + ptsz->has_gal_lensmag_1h
      + ptsz->has_gal_lensmag_2h
      + ptsz->has_gal_gallens_1h
      + ptsz->has_gal_gallens_2h
      + ptsz->has_gal_lensmag_hf
      + ptsz->has_lensmag_lensmag_1h
      + ptsz->has_lensmag_lensmag_2h
      + ptsz->has_lensmag_lensmag_hf
      + ptsz->has_lens_lensmag_1h
      + ptsz->has_lens_lensmag_2h
      + ptsz->has_lens_lensmag_hf
      + ptsz->has_lens_lens_1h
      + ptsz->has_lens_lens_2h
      + ptsz->has_lens_lens_hf
      + ptsz->has_tSZ_lens_1h
      + ptsz->has_tSZ_lens_2h
      + ptsz->has_isw_lens
      + ptsz->has_isw_tsz
      + ptsz->has_isw_auto
      + ptsz->has_dndlnM
      + ptsz->has_sz_counts
      + ptsz->has_sz_rates
      + ptsz->tabulate_rhob_xout_at_m_and_z
      + ptsz->need_ng_bias;


  int mass_conversions = ptsz->need_m200m_to_m200c
                        + ptsz->need_m200c_to_m200m
                        + ptsz->need_m200m_to_m500c
                        + ptsz->need_m200c_to_m500c
                        + ptsz->need_m500c_to_m200c;

  int electron_pressure_comps = ptsz->has_sz_ps
      + ptsz->has_mean_y
      + ptsz->has_dydz
      + ptsz->has_sz_2halo
      + ptsz->has_sz_trispec
      + ptsz->has_sz_m_y_y_1h
      + ptsz->has_sz_m_y_y_2h
      + ptsz->has_sz_te_y_y
      + ptsz->has_tSZ_tSZ_tSZ_1halo
      + ptsz->has_tSZ_tSZ_tSZ_2h
      + ptsz->has_tSZ_tSZ_tSZ_3h
      + ptsz->has_kSZ_kSZ_tSZ_1h
      + ptsz->has_kSZ_kSZ_tSZ_2h
      + ptsz->has_kSZ_kSZ_tSZ_3h
      + ptsz->has_tSZ_gal_1h
      + ptsz->has_tSZ_gal_2h
      + ptsz->has_tSZ_lensmag_1h
      + ptsz->has_tSZ_lensmag_2h
      + ptsz->has_tSZ_cib_1h
      + ptsz->has_tSZ_cib_2h
      + ptsz->has_tSZ_lens_1h
      + ptsz->has_tSZ_lens_2h
      + ptsz->has_isw_tsz
      + ptsz->has_gas_pressure_profile_2h;


  if (all_comps + mass_conversions == _FALSE_){
    if (ptsz->sz_verbose>10) printf("-> freeing nothing.\n");
    return  _SUCCESS_;
  }


  if (ptsz->sz_verbose>10) printf("-> freeing l and cl's.\n");

   free(ptsz->ell);
   if (ptsz->sz_verbose>10) printf("-> freeing l and cl's sz\n");

   free(ptsz->cl_sz_1h);
   if (ptsz->sz_verbose>10) printf("-> freeing l and cl's isw_lens\n");

   free(ptsz->frequencies_for_cib);
   free(ptsz->cib_monopole);
   free(ptsz->cib_shotnoise);
   free(ptsz->cl_isw_lens);
   free(ptsz->cl_isw_tsz);
   free(ptsz->cl_isw_auto);
   free(ptsz->cl_gal_gal_1h);
   free(ptsz->cl_gal_gal_2h);
   free(ptsz->cl_gal_gal_hf);
   free(ptsz->cl_gal_lens_1h);
   free(ptsz->cl_gal_lens_2h);
   free(ptsz->cl_gal_lens_hf);
   free(ptsz->cl_gal_lensmag_1h);
   free(ptsz->cl_gal_lensmag_2h);
   free(ptsz->cl_gal_gallens_1h);
   free(ptsz->cl_gal_gallens_2h);
   if (ptsz->sz_verbose>10) printf("-> freeing l and cl's gallens_gallens\n");

   free(ptsz->cl_gallens_gallens_1h);
   free(ptsz->cl_gallens_gallens_2h);
   free(ptsz->cl_gallens_lens_1h);
   free(ptsz->cl_gallens_lens_2h);
     if ( (ptsz->has_gal_gallens_1h || ptsz->has_gal_gallens_2h) && ptsz->convert_cls_to_gamma){
       free(ptsz->thetas_arcmin);
       free(ptsz->gamma_gal_gallens_1h);
       free(ptsz->gamma_gal_gallens_2h);
     }
   free(ptsz->cl_gal_lensmag_hf);
   free(ptsz->cl_tSZ_lensmag_1h);
   free(ptsz->cl_tSZ_lensmag_2h);
   free(ptsz->cl_lensmag_lensmag_1h);
   free(ptsz->cl_lensmag_lensmag_2h);
   free(ptsz->cl_lensmag_lensmag_hf);
   free(ptsz->cl_lens_lensmag_1h);
   free(ptsz->cl_lens_lensmag_2h);
   free(ptsz->cl_lens_lensmag_hf);
   free(ptsz->cl_tSZ_gal_1h);
   free(ptsz->cl_tSZ_gal_2h);
  if (ptsz->sz_verbose>10) printf("-> freeing l and cl's cib.\n");
   int index_nu;
   int index_nu_prime;
   for (index_nu=0;index_nu<ptsz->cib_frequency_list_num;index_nu++){
       for (index_nu_prime=0;index_nu_prime<ptsz->cib_frequency_list_num;index_nu_prime++){
         free(ptsz->cl_cib_cib_1h[index_nu][index_nu_prime]);
         free(ptsz->cl_cib_cib_2h[index_nu][index_nu_prime]);
       }

     free(ptsz->cl_cib_cib_1h[index_nu]);
     free(ptsz->cl_cib_cib_2h[index_nu]);
     free(ptsz->cl_tSZ_cib_1h[index_nu]);
     free(ptsz->cl_tSZ_cib_2h[index_nu]);
     free(ptsz->cl_lens_cib_1h[index_nu]);
     free(ptsz->cl_lens_cib_2h[index_nu]);
     free(ptsz->cl_gal_cib_1h[index_nu]);
     free(ptsz->cl_gal_cib_2h[index_nu]);
   }
    if (ptsz->has_cib_cib_1h
      + ptsz->has_cib_cib_2h
      + ptsz->has_cib_shotnoise
      + ptsz->has_tSZ_cib_1h
      + ptsz->has_tSZ_cib_2h
      + ptsz->has_gal_cib_1h
      + ptsz->has_gal_cib_2h
      + ptsz->has_lens_cib_1h
      + ptsz->has_lens_cib_2h
      != _FALSE_){
    free(ptsz->cib_frequency_list);
      }


      if (ptsz->has_ngal_ngal_1h
        + ptsz->has_ngal_ngal_2h
        + ptsz->has_ngal_ngal_hf
        + ptsz->has_ngal_lens_1h
        + ptsz->has_ngal_lens_2h
        + ptsz->has_ngal_lens_hf
      ){

  if (ptsz->sz_verbose>10) printf("-> freeing l and cl's ngals.\n");
   int index_g;
   int index_g_prime;
   for (index_g=0;index_g<ptsz->galaxy_samples_list_num;index_g++){
       for (index_g_prime=0;index_g_prime<ptsz->galaxy_samples_list_num;index_g_prime++){
         free(ptsz->cl_ngal_ngal_1h[index_g][index_g_prime]);
         free(ptsz->cl_ngal_ngal_2h[index_g][index_g_prime]);
         free(ptsz->cl_ngal_ngal_hf[index_g][index_g_prime]);
       }
     free(ptsz->cl_ngal_ngal_1h[index_g]);
     free(ptsz->cl_ngal_ngal_2h[index_g]);
     free(ptsz->cl_ngal_ngal_hf[index_g]);
     free(ptsz->cl_ngal_lens_1h[index_g]);
     free(ptsz->cl_ngal_lens_2h[index_g]);
     free(ptsz->cl_ngal_lens_hf[index_g]);

   }
 }

   if (ptsz->has_ngal_ngal_1h
     + ptsz->has_ngal_ngal_2h
     + ptsz->has_ngal_ngal_hf
     + ptsz->has_ngal_lens_1h
     + ptsz->has_ngal_lens_2h
     + ptsz->has_ngal_lens_hf
   ){
     int index_g;
     int index_g_prime;
     for (index_g=0;index_g<ptsz->galaxy_samples_list_num;index_g++){
       if (ptsz->has_ngal_ngal_1h
         + ptsz->has_ngal_ngal_2h
         + ptsz->has_ngal_lens_1h
         + ptsz->has_ngal_lens_2h)
         free(ptsz->array_mean_galaxy_number_density_ngal[index_g]);
       free(ptsz->normalized_dndz_ngal_z[index_g]);
       free(ptsz->normalized_dndz_ngal_phig[index_g]);
     }
     if (ptsz->has_ngal_ngal_1h
       + ptsz->has_ngal_ngal_2h
       + ptsz->has_ngal_lens_1h
       + ptsz->has_ngal_lens_2h)
       free(ptsz->array_mean_galaxy_number_density_ngal);
     free(ptsz->normalized_dndz_ngal_z);
     free(ptsz->normalized_dndz_ngal_phig);

     free(ptsz->galaxy_samples_list);

    if (ptsz->has_ngal_ngal_hf
       +ptsz->has_ngal_lens_hf)
      free(ptsz->effective_galaxy_bias_ngal);
    if (ptsz->has_ngal_ngal_1h
      + ptsz->has_ngal_ngal_2h
      + ptsz->has_ngal_lens_1h
      + ptsz->has_ngal_lens_2h){
          free(ptsz->sigma_log10M_HOD_ngal);//[index_g] = 1.;
          free(ptsz->alpha_s_HOD_ngal);//[index_g] = 1.;
          free(ptsz->M1_prime_HOD_ngal);//[index_g] = 1.;
          free(ptsz->M_min_HOD_ngal);//[index_g] = 1e11;
          free(ptsz->M0_HOD_ngal);//[index_g] = 1e11;
          free(ptsz->x_out_truncated_nfw_profile_satellite_galaxies_ngal);//[index_g] = 1.;
          free(ptsz->f_cen_HOD_ngal);//[index_g] = 1.;
        }
}


   free(ptsz->cl_tSZ_cib_1h);
   free(ptsz->cl_tSZ_cib_2h);
   free(ptsz->cl_lens_cib_1h);
   free(ptsz->cl_lens_cib_2h);
   free(ptsz->cl_gal_cib_1h);
   free(ptsz->cl_gal_cib_2h);
   if (ptsz->has_ngal_ngal_1h
     + ptsz->has_ngal_ngal_2h
     + ptsz->has_ngal_ngal_hf
     + ptsz->has_ngal_lens_1h
     + ptsz->has_ngal_lens_2h
     + ptsz->has_ngal_lens_hf
   ){
   free(ptsz->cl_ngal_ngal_1h);
   free(ptsz->cl_ngal_ngal_2h);
   free(ptsz->cl_ngal_ngal_hf);
   free(ptsz->cl_ngal_lens_1h);
   free(ptsz->cl_ngal_lens_2h);
   free(ptsz->cl_ngal_lens_hf);
 }
   free(ptsz->cl_cib_cib_1h);
   free(ptsz->cl_cib_cib_2h);
   free(ptsz->cl_lens_lens_1h);
   free(ptsz->cl_lens_lens_2h);
   free(ptsz->cl_lens_lens_hf);
   free(ptsz->cl_tSZ_lens_1h);
   free(ptsz->cl_tSZ_lens_2h);
   // free(ptsz->szrate);
   free(ptsz->cl_kSZ_kSZ_gal_1h);
   free(ptsz->cl_kSZ_kSZ_gal_1h_fft);
   free(ptsz->cl_kSZ_kSZ_gal_2h_fft);
   free(ptsz->cl_kSZ_kSZ_gal_3h_fft);
   free(ptsz->cl_kSZ_kSZ_gallens_1h_fft);
   free(ptsz->cl_kSZ_kSZ_gallens_2h_fft);
   free(ptsz->cl_kSZ_kSZ_gallens_3h_fft);
   free(ptsz->cl_kSZ_kSZ_lens_1h_fft);
   free(ptsz->cl_kSZ_kSZ_lens_2h_fft);
   free(ptsz->cl_kSZ_kSZ_lens_3h_fft);
   free(ptsz->cov_ll_kSZ_kSZ_gal);
   free(ptsz->cl_t2t2f);
   free(ptsz->cov_ll_kSZ_kSZ_gallens);
   free(ptsz->cov_ll_kSZ_kSZ_lens);
   free(ptsz->cl_kSZ_kSZ_gal_lensing_term);
   free(ptsz->cl_kSZ_kSZ_gallens_lensing_term);
   free(ptsz->cl_kSZ_kSZ_lens_lensing_term);
   free(ptsz->cl_kSZ_kSZ_gal_2h);
   free(ptsz->cl_kSZ_kSZ_gal_3h);
   free(ptsz->cl_kSZ_kSZ_gal_hf);
   free(ptsz->cl_kSZ_kSZ_gallens_hf);
   free(ptsz->cl_kSZ_kSZ_lens_hf);
   free(ptsz->cl_kSZ_kSZ_lensmag_1h);
   free(ptsz->b_tSZ_tSZ_tSZ_1halo);
   free(ptsz->b_tSZ_tSZ_tSZ_2h);
   free(ptsz->b_tSZ_tSZ_tSZ_3h);
   free(ptsz->b_kSZ_kSZ_tSZ_1h);
   free(ptsz->b_kSZ_kSZ_tSZ_2h);
   free(ptsz->b_kSZ_kSZ_tSZ_3h);
   free(ptsz->cl_te_y_y);
   free(ptsz->m_y_y_1h);
   free(ptsz->m_y_y_2h);
   free(ptsz->cov_cl_cl);
   free(ptsz->sig_cl_squared_binned);
   free(ptsz->cl_sz_2h);

   int index_l;
   for (index_l=0;
        index_l<ptsz->nlSZ;
        index_l++){
          free(ptsz->tllprime_sz[index_l]);
          free(ptsz->trispectrum_ref[index_l]);
          free(ptsz->r_cl_clp[index_l]);
          free(ptsz->cov_Y_N[index_l]);
          free(ptsz->cov_Y_N_next_order[index_l]);
          free(ptsz->r_Y_N[index_l]);
        }


   free(ptsz->tllprime_sz);
   free(ptsz->cov_Y_N);
   free(ptsz->cov_Y_N_next_order);
   int index_multipole_1;
   for (index_multipole_1 = 0; index_multipole_1<ptsz->nlSZ; index_multipole_1 ++){
      free(ptsz->cov_Y_Y_ssc[index_multipole_1]);
}
   free(ptsz->cov_Y_Y_ssc);
   free(ptsz->cov_N_N);
   int index_M_bins_1;
   for (index_M_bins_1 = 0; index_M_bins_1<ptsz->nbins_M; index_M_bins_1 ++){
     free(ptsz->cov_N_N_hsv[index_M_bins_1]);
   }
   free(ptsz->cov_N_N_hsv);
   free(ptsz->r_Y_N);
   free(ptsz->r_cl_clp);
   free(ptsz->trispectrum_ref);
   free(ptsz->ln_k_for_tSZ); //BB: added for class_sz

// printf("free 1\n");

  if (ptsz->has_sz_rates || ptsz->has_sz_counts_fft){
    free(ptsz->szrate);
  }

   if ((ptsz->has_completeness_for_ps_SZ == 1)
    || (ptsz->has_sz_counts  == 1)
    || ptsz->has_sz_counts_fft == 1){
   free(ptsz->thetas);
   free(ptsz->skyfracs);
   free(ptsz->szcat_z);
   free(ptsz->szcat_snr);

   int index_patches;
   for (index_patches=0;
        index_patches<ptsz->nskyfracs;
        index_patches++){
          free(ptsz->ylims[index_patches]);
        }
   free(ptsz->ylims);

   free(ptsz->sky_averaged_ylims);
   free(ptsz->erfs_2d_to_1d_th_array);

 }
  if (ptsz->sz_verbose>10) printf("-> freeing kSZ2X.\n");

if( ptsz->has_kSZ_kSZ_gal_1h
 || ptsz->has_kSZ_kSZ_gal_2h
 || ptsz->has_kSZ_kSZ_gal_3h
 || ptsz->has_kSZ_kSZ_gal_hf
 || ptsz->has_kSZ_kSZ_gallens_hf
 || ptsz->has_kSZ_kSZ_lens_hf
 // || ptsz->has_kSZ_kSZ_gal_covmat // not needed for this one...
 || ptsz->has_kSZ_kSZ_gal_lensing_term
 || ptsz->has_kSZ_kSZ_gallens_lensing_term
 || ptsz->has_kSZ_kSZ_lens_lensing_term
 || ptsz->has_kSZ_kSZ_lensmag_1halo
){
  free(ptsz->ell_kSZ2_gal_multipole_grid);
  free(ptsz->theta_kSZ2_gal_theta_grid);
  // free(ptsz->l_unwise_filter);
  // free(ptsz->f_unwise_filter);
}


if(ptsz->has_kSZ_kSZ_gal_1h
|| ptsz->has_kSZ_kSZ_gal_1h_fft
|| ptsz->has_kSZ_kSZ_gal_2h_fft
|| ptsz->has_kSZ_kSZ_gal_3h_fft
|| ptsz->has_kSZ_kSZ_gal_covmat
|| ptsz->has_kSZ_kSZ_gal_lensing_term
|| ptsz->has_kSZ_kSZ_gallens_1h_fft
|| ptsz->has_kSZ_kSZ_gallens_2h_fft
|| ptsz->has_kSZ_kSZ_gallens_3h_fft
|| ptsz->has_kSZ_kSZ_lens_1h_fft
|| ptsz->has_kSZ_kSZ_lens_2h_fft
|| ptsz->has_kSZ_kSZ_lens_3h_fft
|| ptsz->has_kSZ_kSZ_gallens_covmat
|| ptsz->has_kSZ_kSZ_gallens_lensing_term
|| ptsz->has_kSZ_kSZ_gallens_hf
|| ptsz->has_kSZ_kSZ_lens_covmat
|| ptsz->has_kSZ_kSZ_lens_lensing_term
|| ptsz->has_kSZ_kSZ_lens_hf
|| ptsz->has_kSZ_kSZ_lensmag_1halo
|| ptsz->has_kSZ_kSZ_gal_2h
|| ptsz->has_kSZ_kSZ_gal_3h
|| ptsz->has_kSZ_kSZ_gal_hf
){
  // free(ptsz->ell_kSZ2_gal_multipole_grid);
  // free(ptsz->theta_kSZ2_gal_theta_grid);
  free(ptsz->l_unwise_filter);
  free(ptsz->f_unwise_filter);
}
if(ptsz->has_kSZ_kSZ_gal_1h
|| ptsz->has_kSZ_kSZ_gal_2h
|| ptsz->has_kSZ_kSZ_gal_3h
|| ptsz->has_kSZ_kSZ_gal_1h_fft
|| ptsz->has_kSZ_kSZ_gal_2h_fft
|| ptsz->has_kSZ_kSZ_gal_3h_fft
|| ptsz->has_kSZ_kSZ_gallens_1h_fft
|| ptsz->has_kSZ_kSZ_gallens_2h_fft
|| ptsz->has_kSZ_kSZ_gallens_3h_fft
|| ptsz->has_kSZ_kSZ_lens_1h_fft
|| ptsz->has_kSZ_kSZ_lens_2h_fft
|| ptsz->has_kSZ_kSZ_lens_3h_fft
|| ptsz->has_kSZ_kSZ_1h
|| ptsz->has_kSZ_kSZ_2h
|| ptsz->has_kSZ_kSZ_tSZ_1h
|| ptsz->has_kSZ_kSZ_tSZ_2h
|| ptsz->has_kSZ_kSZ_tSZ_3h
|| ptsz->has_pk_bb_at_z_1h
|| ptsz->has_pk_bb_at_z_2h
|| ptsz->has_pk_b_at_z_2h
|| ptsz->has_gas_density_profile_2h
|| ptsz->has_pk_em_at_z_1h
|| ptsz->has_pk_em_at_z_2h
|| ptsz->has_bk_ttg_at_z_1h
|| ptsz->has_bk_ttg_at_z_2h
|| ptsz->has_bk_ttg_at_z_3h
|| ptsz->has_kSZ_kSZ_lensmag_1halo

){

  if (ptsz->sz_verbose>10) printf("-> freeing density profile kmz.\n");
  free(ptsz->array_profile_ln_k);
  free(ptsz->array_profile_ln_m);
  free(ptsz->array_profile_ln_1pz);

  if (ptsz->sz_verbose>10) printf("-> freeing density profile 2h.\n");
  free(ptsz->array_profile_ln_rho_2h_at_k_and_z);
  free(ptsz->array_profile_rho_2h_at_r_and_z);
  free(ptsz->array_profile_ln_r);


if (ptsz->sz_verbose>10) printf("-> freeing density profile.\n");
 int n_ell = ptsz->n_k_density_profile; //hard coded
 int index_l;
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
  free(ptsz->array_profile_ln_rho_at_lnk_lnM_z[index_l]);
}
free(ptsz->array_profile_ln_rho_at_lnk_lnM_z); //here jump
}

if (ptsz->sz_verbose>10) printf("-> freeing more kSZ2X.\n");

if (ptsz->has_kSZ_kSZ_lensmag_1halo
|| ptsz->has_lensmag_lensmag_1h
|| ptsz->has_lensmag_lensmag_2h
|| ptsz->has_lensmag_lensmag_hf
|| ptsz->has_lens_lensmag_1h
|| ptsz->has_lens_lensmag_2h
|| ptsz->has_gal_lensmag_1h
|| ptsz->has_gal_lensmag_2h
|| ptsz->has_gal_lensmag_hf
|| ptsz->has_tSZ_lensmag_1h
|| ptsz->has_tSZ_lensmag_2h
){
  free(ptsz->array_W_lensmag);
  free(ptsz->array_z_W_lensmag);
}

if (ptsz->sz_verbose>10) printf("-> freeing kappa_g n(z).\n");

if (
 ptsz->has_gal_gallens_1h
|| ptsz->has_gal_gallens_2h
|| ptsz->has_gallens_gallens_1h
|| ptsz->has_gallens_gallens_2h
|| ptsz->has_gallens_lens_1h
|| ptsz->has_gallens_lens_2h
){
  free(ptsz->array_W_gallens_sources);
  free(ptsz->array_z_W_gallens_sources);
}

if (ptsz->sz_verbose>10) printf("-> freeing cib.\n");
if (ptsz->has_cib_cib_1h
  ||ptsz->has_cib_cib_2h
  ||ptsz->has_cib_monopole
  ||ptsz->has_cib_shotnoise
  ||ptsz->has_dcib0dz
  ||ptsz->has_tSZ_cib_1h
  ||ptsz->has_tSZ_cib_2h
  ||ptsz->has_lens_cib_1h
  ||ptsz->has_lens_cib_2h
  ||ptsz->has_gal_cib_1h
  ||ptsz->has_gal_cib_2h
  ){

free(ptsz->array_m_L_sat);

  }

if (ptsz->has_cib_cib_1h
  ||ptsz->has_cib_cib_2h
  ||ptsz->has_cib_monopole
  ||ptsz->has_cib_shotnoise
  ||ptsz->has_dcib0dz
  ||ptsz->has_tSZ_cib_1h
  ||ptsz->has_tSZ_cib_2h
  ||ptsz->has_lens_cib_1h
  ||ptsz->has_lens_cib_2h
  ||ptsz->has_gal_cib_1h
  ||ptsz->has_gal_cib_2h
  ){
free(ptsz->array_z_L_sat);

  }

if (ptsz->has_cib_cib_1h
  ||ptsz->has_cib_cib_2h
  ||ptsz->has_cib_monopole
  ||ptsz->has_cib_shotnoise
  ||ptsz->has_dcib0dz
  ||ptsz->has_tSZ_cib_1h
  ||ptsz->has_tSZ_cib_2h
  ||ptsz->has_lens_cib_1h
  ||ptsz->has_lens_cib_2h
  ||ptsz->has_gal_cib_1h
  ||ptsz->has_gal_cib_2h
  ){

// int index_nu;
// for (index_nu=0;index_nu<ptsz->cib_frequency_list_num;index_nu++)
// {
//   free(ptsz->array_L_sat_at_z_and_M_at_nu[index_nu]);
// }
// free(ptsz->array_L_sat_at_z_and_M_at_nu);
// //free(ptsz->array_L_sat_at_z_and_M_at_nu_prime);
//
//   }
//
// if (ptsz->has_cib_monopole || ptsz->has_dcib0dz){
int index_nu;
for (index_nu=0;index_nu<ptsz->n_nu_L_sat;index_nu++)
{
  free(ptsz->array_L_sat_at_M_z_nu[index_nu]);
}
free(ptsz->array_L_sat_at_M_z_nu);
free(ptsz->array_nu_L_sat);
}

if (ptsz->need_ksz_template){
free(ptsz->l_ksz_template);
free(ptsz->cl_ksz_template);
}

if (ptsz->need_ng_bias){
free(ptsz->array_ln_1pz_ng_bias);
free(ptsz->array_ln_k_ng_bias);
free(ptsz->array_ln_ng_bias_at_z_and_k);
}


if (electron_pressure_comps != _FALSE_){
    if (ptsz->sz_verbose>5) printf("-> freeing pressure.\n");
    if (ptsz->pressure_profile == 3){
       free(ptsz->array_profile_ln_l_over_ls);
       free(ptsz->array_profile_ln_PgNFW_at_lnl_over_ls);
       }

    if(ptsz->pressure_profile == 4){
      if (ptsz->sz_verbose>5) printf("-> freeing pressure B12.\n");
      free(ptsz->array_pressure_profile_ln_k);
      free(ptsz->array_pressure_profile_ln_m);
      free(ptsz->array_pressure_profile_ln_1pz);

       int n_k = ptsz->n_k_pressure_profile; //hard coded
       int index_k;
      for (index_k=0;
           index_k<n_k;
           index_k++)
        {
          free(ptsz->array_pressure_profile_ln_p_at_lnk_lnm_z[index_k]);
        }

      free(ptsz->array_pressure_profile_ln_p_at_lnk_lnm_z);

      }

    if (ptsz->has_gas_pressure_profile_2h){
      free(ptsz->array_pressure_profile_ln_r);
      free(ptsz->array_pressure_profile_2h_ln_k);
      free(ptsz->array_pressure_profile_ln_pressure_2h_at_k_and_z);
      free(ptsz->array_pressure_profile_pressure_2h_at_r_and_z);
    }
  }

if (ptsz->sz_verbose>10) printf("-> freeing miscellaneous.\n");

if (ptsz->has_dcib0dz){
   free(ptsz->array_dcib0dz_redshift);
   free(ptsz->array_dcib0dz_nu);
   free(ptsz->array_dcib0dz_at_z_nu);
   }

if (ptsz->has_electron_density == 1 || ptsz->tabulate_rhob_xout_at_m_and_z ==  1){
if (ptsz->sz_verbose>10) printf("-> freeing xout.\n");
if(ptsz->use_xout_in_density_profile_from_enclosed_mass){
  free(ptsz->array_m_to_xout_redshift);
  free(ptsz->array_m_to_xout_mass);
  free(ptsz->array_m_to_xout_at_z_m);
}
}
if (ptsz->sz_verbose>10) printf("-> freeing dydz.\n");
if (ptsz->has_dydz){
   free(ptsz->array_dydz_redshift);
   free(ptsz->array_dydz_at_z);
   }
if (ptsz->sz_verbose>10) printf("-> freeing flag 12.\n");

if (ptsz->need_m200c_to_m200m){

if (ptsz->sz_verbose>10) printf("-> freeing m200c_to_m200m.\n");
   free(ptsz->array_m_m200c_to_m200m);
   free(ptsz->array_ln_1pz_m200c_to_m200m);
   free(ptsz->array_m200c_to_m200m_at_z_and_M);
   }


if (ptsz->need_m200m_to_m200c){
   free(ptsz->array_m_m200m_to_m200c);
   free(ptsz->array_ln_1pz_m200m_to_m200c);
   free(ptsz->array_m200m_to_m200c_at_z_and_M);
   }


if (ptsz->need_m200m_to_m500c){
   free(ptsz->array_m_m200m_to_m500c);
   free(ptsz->array_ln_1pz_m200m_to_m500c);
   free(ptsz->array_m200m_to_m500c_at_z_and_M);
   }


if (ptsz->need_m200c_to_m500c){
  free(ptsz->array_m_m200c_to_m500c);
  free(ptsz->array_ln_1pz_m200c_to_m500c);
  free(ptsz->array_m200c_to_m500c_at_z_and_M);
  }

if (ptsz->need_m500c_to_m200c){
  free(ptsz->array_m_m500c_to_m200c);
  free(ptsz->array_ln_1pz_m500c_to_m200c);
  free(ptsz->array_m500c_to_m200c_at_z_and_M);
  }

if (ptsz->sz_verbose>10) printf("-> freeing flag 11.\n");

 int index_z;
   for (index_z = 0; index_z<ptsz->N_redshift_dndlnM;index_z ++){
     free(ptsz->dndlnM_at_z_and_M[index_z]);
   }
   free(ptsz->dndlnM_at_z_and_M);

   free(ptsz->dndlnM_array_z);
   free(ptsz->dndlnM_array_m);
// printf("free 2\n");

// if (ptsz->need_hmf + ptsz->has_vrms2 >= 1)

if (ptsz->need_sigma == 1
 || ptsz->has_vrms2){
   free(ptsz->array_redshift);

 }

if (ptsz->need_sigma == 1 ){
   //free(ptsz->array_redshift);

if (ptsz->sz_verbose>10) printf("-> freeing sigma(M,r).\n");

   free(ptsz->array_radius);

   free(ptsz->array_sigma_at_z_and_R);
   free(ptsz->array_dsigma2dR_at_z_and_R);
 }

// if (ptsz->use_class_sz_fast_mode){
  free(ptsz->array_pkl_at_z_and_k);
  free(ptsz->array_pknl_at_z_and_k);
  free(ptsz->array_lnk);
// }


if (ptsz->sz_verbose>10) printf("-> freeing flag 1.\n");

if (ptsz->has_kSZ_kSZ_gal_1h_fft
   || ptsz->has_kSZ_kSZ_gal_2h_fft
   || ptsz->has_kSZ_kSZ_gal_3h_fft
   || ptsz->has_kSZ_kSZ_gal_covmat
   || ptsz->has_kSZ_kSZ_gallens_1h_fft
   || ptsz->has_kSZ_kSZ_gallens_2h_fft
   || ptsz->has_kSZ_kSZ_gallens_3h_fft
   || ptsz->has_kSZ_kSZ_gallens_covmat
   || ptsz->has_kSZ_kSZ_lens_1h_fft
   || ptsz->has_kSZ_kSZ_lens_2h_fft
   || ptsz->has_kSZ_kSZ_lens_3h_fft
   || ptsz->has_kSZ_kSZ_lens_covmat
   || ptsz->convert_cls_to_gamma
   || ptsz->has_pk_b_at_z_2h
   || ptsz->has_sz_counts_fft
   || ptsz->has_gas_density_profile_2h
   || ptsz->has_gas_pressure_profile_2h
   || ptsz->use_fft_for_profiles_transform
   || ptsz->has_n5k){
     if (ptsz->sz_verbose>10) printf("-> destroying fft plans freeing flag 1.\n");
  fftw_destroy_plan(ptsz->forward_plan);
  fftw_destroy_plan(ptsz->reverse_plan);
}

if (ptsz->has_dndlnM == 1
 || ptsz->has_sz_counts
 || ptsz->has_kSZ_kSZ_gal_1h
 || ptsz->has_kSZ_kSZ_gal_1h_fft
 || ptsz->has_kSZ_kSZ_gal_2h_fft
 || ptsz->has_kSZ_kSZ_gal_3h_fft){
   if (ptsz->sz_verbose>10) printf("freeing dndlnM.\n");
   free(ptsz->array_m_dndlnM);
   free(ptsz->array_z_dndlnM);
   free(ptsz->array_dndlnM_at_z_and_M);
 }
if (ptsz->need_hmf == 1){
  if (ptsz->hm_consistency == 1){
free(ptsz->array_hmf_counter_terms_nmin);
free(ptsz->array_redshift_hmf_counter_terms);
free(ptsz->array_hmf_counter_terms_b1min);


  // free(ptsz->array_hmf_counter_terms_nmin);
  // free(ptsz->array_redshift_hmf_counter_terms);
   // free(ptsz->array_hmf_counter_terms_b1min);
   free(ptsz->array_hmf_counter_terms_b2min);
 }
 }

if (ptsz->has_vrms2)
   free(ptsz->array_vrms2_at_z);
if (ptsz->has_knl)
  free(ptsz->array_knl_at_z);
if (ptsz->has_nl_index){
  free(ptsz->array_nl_index_at_z_and_k);
  free(ptsz->array_nl_index_at_z_and_k_no_wiggles);
}


if (ptsz->sz_verbose>10) printf("-> freeing flag 2.\n");
if (ptsz->has_tSZ_gal_1h
   || ptsz->has_tSZ_gal_2h

   || ptsz->has_kSZ_kSZ_gal_1h_fft
   || ptsz->has_kSZ_kSZ_gal_2h_fft
   || ptsz->has_kSZ_kSZ_gal_3h_fft
   || ptsz->has_kSZ_kSZ_gal_1h
   || ptsz->has_kSZ_kSZ_gal_2h
   || ptsz->has_kSZ_kSZ_gal_3h
   || ptsz->has_bk_ttg_at_z_1h
   || ptsz->has_bk_ttg_at_z_2h
   || ptsz->has_bk_ttg_at_z_3h
   || ptsz->has_kSZ_kSZ_gal_hf
   || ptsz->has_kSZ_kSZ_lensmag_1halo
   || ptsz->has_gal_gal_1h
   || ptsz->has_gal_gal_2h
   || ptsz->has_pk_gg_at_z_1h
   || ptsz->has_pk_gg_at_z_2h
   || ptsz->has_gal_lens_1h
   || ptsz->has_gal_lens_2h
   || ptsz->has_gal_cib_1h
   || ptsz->has_gal_cib_2h
   || ptsz->has_gal_lensmag_1h
   || ptsz->has_gal_lensmag_2h
   || ptsz->has_gal_gallens_1h
   || ptsz->has_gal_gallens_2h
   // || ptsz->has_tSZ_lensmag_1h
   // || ptsz->has_tSZ_lensmag_2h
   // || ptsz->has_lensmag_lensmag_1h
   // || ptsz->has_lensmag_lensmag_2h
   // || ptsz->has_lens_lensmag_1h
   // || ptsz->has_lens_lensmag_2h
   ){
   free(ptsz->array_mean_galaxy_number_density);

  }

  if (ptsz->has_mean_galaxy_bias){
    free(ptsz->array_mean_galaxy_bias);
  }


if (ptsz->has_tSZ_gal_1h
   || ptsz->has_tSZ_gal_2h
   || ptsz->has_kSZ_kSZ_gal_1h
   || ptsz->has_kSZ_kSZ_gal_1h_fft
   || ptsz->has_kSZ_kSZ_gal_2h_fft
   || ptsz->has_kSZ_kSZ_gal_3h_fft
   || ptsz->has_kSZ_kSZ_gal_2h
   || ptsz->has_kSZ_kSZ_gal_3h
   || ptsz->has_kSZ_kSZ_gal_hf
   || ptsz->has_bk_ttg_at_z_2h
   || ptsz->has_bk_ttg_at_z_3h
   || ptsz->has_bk_ttg_at_z_1h
   || ptsz->has_kSZ_kSZ_lensmag_1halo
   || ptsz->has_gal_gal_1h
   || ptsz->has_gal_gal_2h
   || ptsz->has_gal_cib_1h
   || ptsz->has_gal_cib_2h
   || ptsz->has_gal_gal_hf
   || ptsz->has_gal_lens_1h
   || ptsz->has_gal_lens_2h
   || ptsz->has_gal_lens_hf
   || ptsz->has_gal_lensmag_1h
   || ptsz->has_gal_lensmag_2h
   || ptsz->has_gal_gallens_1h
   || ptsz->has_gal_gallens_2h
   || ptsz->has_gal_lensmag_hf
   || ptsz->has_tSZ_lensmag_1h
   || ptsz->has_tSZ_lensmag_2h
   || ptsz->has_lensmag_lensmag_1h
   || ptsz->has_lensmag_lensmag_2h
   || ptsz->has_lensmag_lensmag_hf
   || ptsz->has_lens_lensmag_1h
   || ptsz->has_lens_lensmag_2h
   || ptsz->has_lens_lensmag_hf
   ){
   // free(ptsz->array_mean_galaxy_number_density);
   free(ptsz->normalized_dndz_z);
   free(ptsz->normalized_dndz_phig);

   // if ( ptsz->has_gal_gallens_1h
   //    || ptsz->has_gal_gallens_2h){
   // free(ptsz->normalized_source_dndz_z);
   // free(ptsz->normalized_source_dndz_phig);
   //    }
  // unwise
  if (ptsz->galaxy_sample ==  1){
    free(ptsz->normalized_fdndz_z);
    free(ptsz->normalized_fdndz_phig);
  //
    free(ptsz->normalized_cosmos_dndz_z);
    free(ptsz->normalized_cosmos_dndz_phig);
   }
  }

if (ptsz->sz_verbose>10) printf("-> freeing flag 3.\n");
if (
    ptsz->has_gal_gallens_1h
   || ptsz->has_gal_gallens_2h
   || ptsz->has_gallens_gallens_1h
   || ptsz->has_gallens_gallens_2h
   || ptsz->has_gallens_lens_1h
   || ptsz->has_gallens_lens_2h
   || ptsz->has_kSZ_kSZ_gallens_1h_fft
   || ptsz->has_kSZ_kSZ_gallens_2h_fft
   || ptsz->has_kSZ_kSZ_gallens_3h_fft
   || ptsz->has_kSZ_kSZ_gallens_hf
   ){

   free(ptsz->normalized_source_dndz_z);
   free(ptsz->normalized_source_dndz_phig);
  }


if (ptsz->sz_verbose>10) printf("-> freeing flag 4.\n");
if (ptsz->include_noise_cov_y_y==1){
   free(ptsz->unbinned_nl_yy_ell);
   free(ptsz->unbinned_nl_yy_n_ell);
}
   if (ptsz->has_sigma2_hsv)
    free(ptsz->array_sigma2_hsv_at_z);

if(ptsz->need_tt_noise){
  free(ptsz->unbinned_nl_tt_ell);
  free(ptsz->unbinned_nl_tt_n_ell);
}

if(ptsz->need_lensing_noise){
if (ptsz->sz_verbose>10) printf("-> freeing lensing noise.\n");
  free(ptsz->nl_lensing_noise);
  free(ptsz->l_lensing_noise);
}

if(ptsz->use_redshift_dependent_M_min){
if (ptsz->sz_verbose>10) printf("-> freeing redshift_dependent_M_min.\n");
  free(ptsz->M_min_of_z_z);
  free(ptsz->M_min_of_z_M_min);
}


   if (ptsz->has_tSZ_gal_1h
      +ptsz->has_tSZ_gal_2h
      +ptsz->has_sz_te_y_y
      +ptsz->has_sz_trispec
      +ptsz->has_sz_m_y_y_1h
      +ptsz->has_sz_m_y_y_2h
      +ptsz->has_sz_cov_Y_N
      +ptsz->has_sz_cov_Y_Y_ssc
      +ptsz->has_sz_cov_Y_N_next_order
      +ptsz->has_tSZ_lensmag_2h
      +ptsz->has_tSZ_lensmag_1h
      +ptsz->has_tSZ_gal_1h
      +ptsz->has_tSZ_gal_2h
      +ptsz->has_tSZ_cib_1h
      +ptsz->has_tSZ_cib_2h
      +ptsz->has_tSZ_lens_1h
      +ptsz->has_tSZ_lens_2h
      +ptsz->has_tSZ_tSZ_tSZ_1halo
      +ptsz->has_tSZ_tSZ_tSZ_2h
      +ptsz->has_tSZ_tSZ_tSZ_3h
      +ptsz->has_kSZ_kSZ_tSZ_1h
      +ptsz->has_kSZ_kSZ_tSZ_2h
      +ptsz->has_kSZ_kSZ_tSZ_3h
      +ptsz->has_sz_ps
      +ptsz->has_mean_y
      +ptsz->has_dydz
      +ptsz->has_sz_2halo
      != 0){
if (ptsz->pressure_profile == 0 || ptsz->pressure_profile == 2 )
{

if (ptsz->sz_verbose>10) printf("-> freeing pressure profile.\n");
   free(ptsz->PP_lnx);
   free(ptsz->PP_lnI);
   free(ptsz->PP_d2lnI);}

 }

if( ptsz->has_pk_at_z_1h
   + ptsz->has_pk_at_z_2h
   + ptsz->has_pk_gg_at_z_1h
   + ptsz->has_pk_gg_at_z_2h
   + ptsz->has_pk_bb_at_z_1h
   + ptsz->has_pk_bb_at_z_2h
   + ptsz->has_pk_b_at_z_2h
   + ptsz->has_pk_em_at_z_1h
   + ptsz->has_pk_em_at_z_2h
   + ptsz->has_pk_HI_at_z_1h
   + ptsz->has_pk_HI_at_z_2h
   + ptsz->has_bk_at_z_1h
   + ptsz->has_bk_at_z_2h
   + ptsz->has_bk_at_z_3h
   + ptsz->has_bk_ttg_at_z_1h
   + ptsz->has_bk_ttg_at_z_2h
   + ptsz->has_bk_ttg_at_z_3h

    >= _TRUE_){
if (ptsz->sz_verbose>10) printf("-> freeing pk's and bk's.\n");
free(ptsz->k_for_pk_hm);
free(ptsz->pk_at_z_1h);
free(ptsz->pk_at_z_2h);
free(ptsz->pk_gg_at_z_1h);
free(ptsz->pk_gg_at_z_2h);
free(ptsz->pk_bb_at_z_1h);
free(ptsz->pk_bb_at_z_2h);
free(ptsz->pk_b_at_z_2h);
free(ptsz->pk_em_at_z_1h);
free(ptsz->pk_em_at_z_2h);
free(ptsz->pk_HI_at_z_1h);
free(ptsz->pk_HI_at_z_2h);
free(ptsz->bk_at_z_1h);
free(ptsz->bk_at_z_2h);
free(ptsz->bk_at_z_3h);
free(ptsz->bk_ttg_at_z_1h);
free(ptsz->bk_ttg_at_z_2h);
free(ptsz->bk_ttg_at_z_3h);
}
//



  // if (ptsz->MF==1 && ptsz->hm_consistency==2){
  if (ptsz->MF==1){
  if (ptsz->T10_alpha_fixed==0){
if (ptsz->sz_verbose>10) printf("-> freeing HMF alpha (T10).\n");
    free(ptsz->T10_ln1pz);
    free(ptsz->T10_lnalpha);
  }
}

if (ptsz->sz_verbose>10) printf("-> freeing c-m relation.\n");
  // }
if (ptsz->concentration_parameter == 4){
  free(ptsz->CM_redshift);
  free(ptsz->CM_logM);
  free(ptsz->CM_logC);
}

  free(ptsz->M_bins);
  free(ptsz->cov_Y_N_mass_bin_edges);

  free(ptsz->ln_x_for_pp);
  free(ptsz->x_for_pp);

if (ptsz->sz_verbose>10) printf("-> freeing cluster counts.\n");
if (ptsz->has_sz_counts == 1){
  free(ptsz->steps_z);
  free(ptsz->steps_m);
  free(ptsz->erfs_2d_to_1d_y_array);

if (ptsz->sz_verbose>10) printf("-> freeing cluster counts freed.\n");
}

if (ptsz->sz_verbose>10) printf("-> freeing more kSZ2X.\n");
if ( ptsz->has_kSZ_kSZ_gal_1h_fft
  || ptsz->has_kSZ_kSZ_gal_2h_fft
  || ptsz->has_kSZ_kSZ_gal_3h_fft
  || ptsz->has_kSZ_kSZ_gal_3h
  || ptsz->has_kSZ_kSZ_gallens_1h_fft
  || ptsz->has_kSZ_kSZ_gallens_2h_fft
  || ptsz->has_kSZ_kSZ_gallens_3h_fft
  || ptsz->has_kSZ_kSZ_lens_1h_fft
  || ptsz->has_kSZ_kSZ_lens_2h_fft
  || ptsz->has_kSZ_kSZ_lens_3h_fft
   ){
  free(ptsz->array_psi_b2t_redshift);
  free(ptsz->array_psi_b2t_multipole);
  free(ptsz->array_psi_b2t_psi);

  free(ptsz->array_psi_b1t_redshift);
  free(ptsz->array_psi_b1t_multipole);
  free(ptsz->array_psi_b1t_psi);
}

if (ptsz->has_n5k){
  free(ptsz->array_n5k_F1_F);
  free(ptsz->array_n5k_F1_k);
  free(ptsz->array_n5k_F1_l);
  free(ptsz->n5k_pk_z);
  free(ptsz->n5k_pk_k);
  free(ptsz->n5k_pk_pk);
  free(ptsz->n5k_cl_K1_K1);
  free(ptsz->n5k_cl_K1_chi);
  free(ptsz->n5k_z_of_chi_z);
  free(ptsz->n5k_z_of_chi_chi);
}

if (ptsz->use_maniyar_cib_model && ptsz->has_cib){
  free(ptsz->cib_Snu_z);
  free(ptsz->cib_Snu_nu);
  free(ptsz->cib_Snu_snu);
}

if (ptsz->sz_verbose>10) printf("-> freeing more kSZ2X.\n");
if ( ptsz->has_kSZ_kSZ_gal_1h_fft
  || ptsz->has_kSZ_kSZ_gal_2h_fft
  || ptsz->has_kSZ_kSZ_gal_3h_fft
  || ptsz->has_kSZ_kSZ_gal_3h
   ){
  free(ptsz->array_psi_b1g_redshift);
  free(ptsz->array_psi_b1g_multipole);
  free(ptsz->array_psi_b1g_psi);

  free(ptsz->array_psi_b2g_redshift);
  free(ptsz->array_psi_b2g_multipole);
  free(ptsz->array_psi_b2g_psi);

  free(ptsz->array_psi_b1gt_redshift);
  free(ptsz->array_psi_b1gt_multipole);
  int iz;
  for (iz=0;iz<ptsz->n_z_psi_b1gt;iz++)
    free(ptsz->array_psi_b1gt_psi[iz]);
  free(ptsz->array_psi_b1gt_psi);
}

  if (ptsz->sz_verbose>10) printf("-> freeing more kSZ2X.\n");
if ( ptsz->has_kSZ_kSZ_gallens_1h_fft
  || ptsz->has_kSZ_kSZ_gallens_2h_fft
  || ptsz->has_kSZ_kSZ_gallens_3h_fft
  || ptsz->has_kSZ_kSZ_lens_1h_fft
  || ptsz->has_kSZ_kSZ_lens_2h_fft
  || ptsz->has_kSZ_kSZ_lens_3h_fft
   ){
  free(ptsz->array_psi_b1kg_redshift);
  free(ptsz->array_psi_b1kg_multipole);
  free(ptsz->array_psi_b1kg_psi);

  free(ptsz->array_psi_b2kg_redshift);
  free(ptsz->array_psi_b2kg_multipole);
  free(ptsz->array_psi_b2kg_psi);

  free(ptsz->array_psi_b1kgt_redshift);
  free(ptsz->array_psi_b1kgt_multipole);
  int iz;
  for (iz=0;iz<ptsz->n_z_psi_b1kgt;iz++)
    free(ptsz->array_psi_b1kgt_psi[iz]);
  free(ptsz->array_psi_b1kgt_psi);
}

// printf("free 3\n");
  if (ptsz->sz_verbose>1) printf("-> memory freed.\n");
  if (ptsz->sz_verbose>1) printf("-> Exiting class_sz.\n");
return _SUCCESS_;
}


int load_cl_ksz_template(struct tszspectrum * ptsz)
{

  if (ptsz->sz_verbose >= 1)
    printf("-> loading the ksz template\n");


  class_alloc(ptsz->l_ksz_template,sizeof(double *)*100,ptsz->error_message);
  class_alloc(ptsz->cl_ksz_template,sizeof(double *)*100,ptsz->error_message);

  char line[_LINE_LENGTH_MAX_];
  //char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  FILE *process,*process2;
  int n_data_guess, n_data = 0;
  int n_data_guess2, n_data2 = 0;
  double *lnx = NULL, *lnI = NULL, *lnJ = NULL,  *tmp = NULL;
  double this_lnx, this_lnI;
  int status;
  int index_x;


  /** 1. Initialization */
  /* Prepare the data (with some initial size) */
  n_data_guess = 100;
  n_data_guess2 = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));
  lnI = (double *)malloc(n_data_guess*sizeof(double));
  lnJ = (double *)malloc(n_data_guess*sizeof(double));


  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];

  class_open(process,ptsz->ksz_template_file, "r",ptsz->error_message);
  if (ptsz->sz_verbose >= 1)
    printf("-> File Name: %s\n",ptsz->ksz_template_file);

  class_open(process2,ptsz->ksz_reio_template_file, "r",ptsz->error_message);
  if (ptsz->sz_verbose >= 1)
    printf("-> File Name: %s\n",ptsz->ksz_reio_template_file);

  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf %lf", &this_lnx, &this_lnI);
    //printf("lnx = %e\n",this_lnx);




    /* Standard technique in C:
     /*if too many data, double the size of the vectors */
    /* (it is faster and safer that reallocating every new line) */
    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 ptsz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 ptsz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnI = tmp;
    };
    /* Store */
    lnx[n_data]   = this_lnx;
    lnI[n_data]   = this_lnI;

    n_data++;
    /* Check ascending order of the k's */
    if(n_data>1) {
      class_test(lnx[n_data-1] <= lnx[n_data-2],
                 ptsz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");


  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process2) != NULL) {
    sscanf(line, "%lf %lf", &this_lnx, &this_lnI);
    //printf("lnx = %e\n",this_lnx);




    /* Standard technique in C:
     /*if too many data, double the size of the vectors */
    /* (it is faster and safer that reallocating every new line) */
    if((n_data2+1) > n_data_guess2) {
      n_data_guess2 *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess2*sizeof(double));
      class_test(tmp == NULL,
                 ptsz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnJ, n_data_guess2*sizeof(double));
      class_test(tmp == NULL,
                 ptsz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnJ = tmp;
    };
    /* Store */
    lnx[n_data2]   = this_lnx;
    lnJ[n_data2]   = this_lnI;

    n_data2++;
    /* Check ascending order of the k's */
    if(n_data2>1) {
      class_test(lnx[n_data2-1] <= lnx[n_data2-2],
                 ptsz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  // status = pclose(process);
  status = fclose(process2);
  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");




  /** 3. Store the read results into CLASS structures */
  ptsz->ksz_template_size = n_data;
  /** Make room */

  class_realloc(ptsz->l_ksz_template,
                ptsz->l_ksz_template,
                ptsz->ksz_template_size*sizeof(double),
                ptsz->error_message);
  class_realloc(ptsz->cl_ksz_template,
                ptsz->cl_ksz_template,
                ptsz->ksz_template_size*sizeof(double),
                ptsz->error_message);



  /** Store them */
  for (index_x=0; index_x<ptsz->ksz_template_size; index_x++) {
    ptsz->l_ksz_template[index_x] = lnx[index_x];
    ptsz->cl_ksz_template[index_x] = lnI[index_x]+lnJ[index_x];
// printf("%.3e  %.3e  %.3e  %.3e\n",
// ptsz->l_ksz_template[index_x],
// ptsz->cl_ksz_template[index_x],
// lnI[index_x],
// lnJ[index_x]
// );
  };


  /** Release the memory used locally */
  free(lnx);
  free(lnI);
  free(lnJ);

  return _SUCCESS_;
}





int compute_sz(struct background * pba,
                struct nonlinear * pnl,
                struct primordial * ppm,
                struct perturbs * ppt,
                struct tszspectrum * ptsz,
                double * Pvecback,
                double * Pvectsz){



   int index_integrand = (int) Pvectsz[ptsz->index_integrand_id];


   // printf("index_integrand = %d %d\n",index_integrand,ppt->has_density_transfers);
   //Pvectsz[ptsz->index_multipole] = 0.;
   if (index_integrand>=ptsz->index_integrand_id_dndlnM_first && index_integrand <= ptsz->index_integrand_id_dndlnM_last && ptsz->has_dndlnM){
      Pvectsz[ptsz->index_md] = ptsz->index_md_dndlnM;

      int index_redshift_mass = (int) (index_integrand - ptsz->index_integrand_id_dndlnM_first);
      int index_redshift = index_redshift_mass / ptsz->N_mass_dndlnM;
      int index_mass = index_redshift_mass % ptsz->N_mass_dndlnM;
      Pvectsz[ptsz->index_redshift_for_dndlnM] = (double)  (index_redshift);
      Pvectsz[ptsz->index_mass_for_dndlnM] = (double) (index_mass);


      if (ptsz->sz_verbose > 0 && index_integrand==ptsz->index_integrand_id_dndlnM_first) printf("computing dndlnM @ redshift_id = %.0f and mass_id = %.0f\n",Pvectsz[ptsz->index_redshift_for_dndlnM], Pvectsz[ptsz->index_mass_for_dndlnM]);
      if (ptsz->sz_verbose > 0 && index_integrand==ptsz->index_integrand_id_dndlnM_last) printf("computing dndlnM @ redshift_id = %.0f and mass_id = %.0f\n",Pvectsz[ptsz->index_redshift_for_dndlnM], Pvectsz[ptsz->index_mass_for_dndlnM]);
    }
   else if(index_integrand == ptsz->index_integrand_id_hmf && ptsz->has_hmf){
      Pvectsz[ptsz->index_md] = ptsz->index_md_hmf;
      if (ptsz->sz_verbose > 0) printf("computing hmf int\n");
   }
   else if (index_integrand>=ptsz->index_integrand_id_szrates_first && index_integrand <= ptsz->index_integrand_id_szrates_last && ptsz->has_sz_rates){
      Pvectsz[ptsz->index_md] = ptsz->index_md_szrates;

      Pvectsz[ptsz->index_szrate] = (double) (index_integrand - ptsz->index_integrand_id_szrates_first);

      if (ptsz->sz_verbose > 1) printf("computing szrates @ cluster id = %.0f\n",Pvectsz[ptsz->index_szrate]);
    }

   else if (index_integrand == ptsz->index_integrand_id_mean_y && ptsz->has_mean_y) {
      Pvectsz[ptsz->index_md] = ptsz->index_md_mean_y;
      Pvectsz[ptsz->index_has_electron_pressure] = 1;
      if (ptsz->sz_verbose > 0) printf("computing mean y\n");
   }
   else if (index_integrand>=ptsz->index_integrand_id_sz_ps_first && index_integrand <= ptsz->index_integrand_id_sz_ps_last && ptsz->has_sz_ps){
      Pvectsz[ptsz->index_md] = ptsz->index_md_sz_ps;
      Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_sz_ps_first);
      Pvectsz[ptsz->index_has_electron_pressure] = 1;
      if (ptsz->sz_verbose == 1 && index_integrand ==ptsz->index_integrand_id_sz_ps_first ) printf("computing cl^yy @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      if (ptsz->sz_verbose == 1 && index_integrand ==ptsz->index_integrand_id_sz_ps_last ) printf("computing cl^yy @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      if (ptsz->sz_verbose >1) printf("computing cl^yy @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
    }
   else if (index_integrand>=ptsz->index_integrand_id_cib_monopole_first && index_integrand <= ptsz->index_integrand_id_cib_monopole_last && ptsz->has_cib_monopole){
      Pvectsz[ptsz->index_md] = ptsz->index_md_cib_monopole;
      Pvectsz[ptsz->index_frequency_for_cib_profile] = (double) (index_integrand - ptsz->index_integrand_id_cib_monopole_first);
      // Pvectsz[ptsz->index_has_electron_pressure] = 1;
      Pvectsz[ptsz->index_has_cib] = 1;
      if (ptsz->sz_verbose == 1 && index_integrand ==ptsz->index_integrand_id_cib_monopole_first ) printf("computing CIB Inu @ nu_id = %.0f\n",Pvectsz[ptsz->index_frequency_for_cib_profile]);
      if (ptsz->sz_verbose == 1 && index_integrand ==ptsz->index_integrand_id_cib_monopole_last ) printf("computing CIB Inu @ nu_id = %.0f\n",Pvectsz[ptsz->index_frequency_for_cib_profile]);
      if (ptsz->sz_verbose >1) printf("computing cib_monopole @ nu_id = %.0f\n",Pvectsz[ptsz->index_frequency_for_cib_profile]);
    }
   else if (index_integrand>=ptsz->index_integrand_id_cib_shotnoise_first && index_integrand <= ptsz->index_integrand_id_cib_shotnoise_last && ptsz->has_cib_shotnoise){
      Pvectsz[ptsz->index_md] = ptsz->index_md_cib_shotnoise;
      Pvectsz[ptsz->index_frequency_for_cib_profile] = (double) (index_integrand - ptsz->index_integrand_id_cib_shotnoise_first);
      // Pvectsz[ptsz->index_has_electron_pressure] = 1;
      Pvectsz[ptsz->index_has_cib] = 1;
      if (ptsz->sz_verbose == 1 && index_integrand ==ptsz->index_integrand_id_cib_shotnoise_first ) printf("computing CIB Inu @ nu_id = %.0f\n",Pvectsz[ptsz->index_frequency_for_cib_profile]);
      if (ptsz->sz_verbose == 1 && index_integrand ==ptsz->index_integrand_id_cib_shotnoise_last ) printf("computing CIB Inu @ nu_id = %.0f\n",Pvectsz[ptsz->index_frequency_for_cib_profile]);
      if (ptsz->sz_verbose >1) printf("computing cib_shotnoise @ nu_id = %.0f\n",Pvectsz[ptsz->index_frequency_for_cib_profile]);
    }

    else if (index_integrand>=ptsz->index_integrand_id_trispectrum_first && index_integrand <= ptsz->index_integrand_id_trispectrum_last && ptsz->has_sz_trispec){
       Pvectsz[ptsz->index_md] = ptsz->index_md_trispectrum;
       int index_ell_ell_prime = (int) (index_integrand - ptsz->index_integrand_id_trispectrum_first);
       int n = (-1.+sqrt(1. + 4.*2.*index_ell_ell_prime))/2.;
       int index_ell = floor(n);
       int index_ell_prime = index_ell_ell_prime -index_ell*(index_ell+1)/2;
       Pvectsz[ptsz->index_multipole] = (double) (index_ell);
       Pvectsz[ptsz->index_multipole_prime] = (double) (index_ell_prime);
       Pvectsz[ptsz->index_has_electron_pressure] = 1;
       if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_trispectrum_first) printf("computing trispectrum @ ell_id = %.0f, ell_id_prime = %.0f\n",Pvectsz[ptsz->index_multipole],Pvectsz[ptsz->index_multipole_prime]);
       if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_trispectrum_last) printf("computing trispectrum @ ell_id = %.0f, ell_id_prime = %.0f\n",Pvectsz[ptsz->index_multipole],Pvectsz[ptsz->index_multipole_prime]);
       if (ptsz->sz_verbose > 1) printf("computing trispectrum @ ell_id = %.0f, ell_id_prime = %.0f\n",Pvectsz[ptsz->index_multipole],Pvectsz[ptsz->index_multipole_prime]);
     }
   else if (index_integrand>=ptsz->index_integrand_id_sz_ps_2halo_first && index_integrand <= ptsz->index_integrand_id_sz_ps_2halo_last && ptsz->has_sz_2halo){
      Pvectsz[ptsz->index_md] = ptsz->index_md_2halo;
      Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_sz_ps_2halo_first);
      Pvectsz[ptsz->index_has_electron_pressure] = 1;
      if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_sz_ps_2halo_first) printf("computing cl^yy 2-halo term @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_sz_ps_2halo_last) printf("computing cl^yy 2-halo term @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      if (ptsz->sz_verbose > 1) printf("computing cl^yy 2-halo term @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
    }


   else if (index_integrand>=ptsz->index_integrand_id_sz_ps_te_y_y_first && index_integrand <= ptsz->index_integrand_id_sz_ps_te_y_y_last && ptsz->has_sz_te_y_y){
      Pvectsz[ptsz->index_md] = ptsz->index_md_te_y_y;
      Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_sz_ps_te_y_y_first);
      Pvectsz[ptsz->index_has_electron_pressure] = 1;
      if (ptsz->sz_verbose > 0) printf("computing cl^Teyy @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
    }

    else if (index_integrand>=ptsz->index_integrand_id_sz_ps_m_y_y_1h_first && index_integrand <= ptsz->index_integrand_id_sz_ps_m_y_y_1h_last && ptsz->has_sz_m_y_y_1h){
       Pvectsz[ptsz->index_md] = ptsz->index_md_m_y_y_1h;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_sz_ps_m_y_y_1h_first);
       Pvectsz[ptsz->index_has_electron_pressure] = 1;
       if (ptsz->sz_verbose > 0) printf("computing m_y_y_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
     else if (index_integrand>=ptsz->index_integrand_id_sz_ps_m_y_y_2h_first && index_integrand <= ptsz->index_integrand_id_sz_ps_m_y_y_2h_last && ptsz->has_sz_m_y_y_2h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_m_y_y_2h;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_sz_ps_m_y_y_2h_first);
        Pvectsz[ptsz->index_has_electron_pressure] = 1;
        if (ptsz->sz_verbose > 0) printf("computing m_y_y_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
    else if (index_integrand>=ptsz->index_integrand_id_cov_Y_N_first && index_integrand <= ptsz->index_integrand_id_cov_Y_N_last && ptsz->has_sz_cov_Y_N){
       Pvectsz[ptsz->index_md] = ptsz->index_md_cov_Y_N;
       int index_ell_mass = (int) (index_integrand - ptsz->index_integrand_id_cov_Y_N_first);
       int index_ell = index_ell_mass / ptsz->nbins_M;
       int index_mass = index_ell_mass % ptsz->nbins_M;
       Pvectsz[ptsz->index_multipole] = (double) (index_ell);
       Pvectsz[ptsz->index_mass_bin_1] = (double) (index_mass);
       Pvectsz[ptsz->index_has_electron_pressure] = 1;
       if (ptsz->sz_verbose > 0 && index_integrand==ptsz->index_integrand_id_cov_Y_N_first) printf("computing cov(Y,N) @ ell_id = %.0f, mass_bin_id = %.0f\n",Pvectsz[ptsz->index_multipole],Pvectsz[ptsz->index_mass_bin_1]);
       if (ptsz->sz_verbose > 0 && index_integrand==ptsz->index_integrand_id_cov_Y_N_last) printf("computing cov(Y,N) @ ell_id = %.0f, mass_bin_id = %.0f\n",Pvectsz[ptsz->index_multipole],Pvectsz[ptsz->index_mass_bin_1]);
     }
    else if (index_integrand>=ptsz->index_integrand_id_cov_N_N_first && index_integrand <= ptsz->index_integrand_id_cov_N_N_last && ptsz->has_sz_cov_N_N){

        Pvectsz[ptsz->index_md] = ptsz->index_md_cov_N_N;
        int index_mass_bin_1 = (int) (index_integrand - ptsz->index_integrand_id_cov_N_N_first);
        Pvectsz[ptsz->index_mass_bin_1] = (double) (index_mass_bin_1);

        if (ptsz->sz_verbose > 0) printf("computing cov(N,N) @ mass_bin_id = %.0f\n",Pvectsz[ptsz->index_mass_bin_1]);
      }
    else if (index_integrand>=ptsz->index_integrand_id_cov_N_N_hsv_first && index_integrand <= ptsz->index_integrand_id_cov_N_N_hsv_last && ptsz->has_sz_cov_N_N_hsv){

       Pvectsz[ptsz->index_md] = ptsz->index_md_cov_N_N_hsv;
       int index_mass_bin_1_mass_bin_2 = (int) (index_integrand - ptsz->index_integrand_id_cov_N_N_hsv_first);
       int n = (-1.+sqrt(1. + 4.*2.*index_mass_bin_1_mass_bin_2))/2.;
       int index_mass_bin_1 = floor(n);
       int index_mass_bin_2 = index_mass_bin_1_mass_bin_2 -index_mass_bin_1*(index_mass_bin_1+1)/2;
       Pvectsz[ptsz->index_mass_bin_1] = (double) (index_mass_bin_1);
       Pvectsz[ptsz->index_mass_bin_2] = (double) (index_mass_bin_2);

       if (ptsz->sz_verbose > 0) printf("computing cov(N,N) [hsv] @ mass_bin_1_id = %.0f, mass_bin_2_id = %.0f\n",Pvectsz[ptsz->index_mass_bin_1],Pvectsz[ptsz->index_mass_bin_2]);
     }
     else if (index_integrand>=ptsz->index_integrand_id_cov_Y_N_next_order_first && index_integrand <= ptsz->index_integrand_id_cov_Y_N_next_order_last && ptsz->has_sz_cov_Y_N_next_order){
        Pvectsz[ptsz->index_md] = ptsz->index_md_cov_Y_N_next_order;
        int index_ell_mass = (int) (index_integrand - ptsz->index_integrand_id_cov_Y_N_next_order_first);
        int index_ell = index_ell_mass / ptsz->nbins_M;
        int index_mass = index_ell_mass % ptsz->nbins_M;
        Pvectsz[ptsz->index_multipole] = (double) (index_ell);
        Pvectsz[ptsz->index_mass_bin_1] = (double) (index_mass);
        Pvectsz[ptsz->index_has_electron_pressure] = 1;
        if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_cov_Y_N_next_order_first) printf("computing cov(Y,N) [ssc] @ ell_id = %.0f, mass_bin_id = %.0f\n",Pvectsz[ptsz->index_multipole],Pvectsz[ptsz->index_mass_bin_1]);
        if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_cov_Y_N_next_order_last) printf("computing cov(Y,N) [ssc] @ ell_id = %.0f, mass_bin_id = %.0f\n",Pvectsz[ptsz->index_multipole],Pvectsz[ptsz->index_mass_bin_1]);
        if (ptsz->sz_verbose > 1) printf("computing cov(Y,N) [hsv] @ ell_id = %.0f, mass_bin_id = %.0f\n",Pvectsz[ptsz->index_multipole],Pvectsz[ptsz->index_mass_bin_1]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_cov_Y_Y_ssc_first && index_integrand <= ptsz->index_integrand_id_cov_Y_Y_ssc_last && ptsz->has_sz_cov_Y_Y_ssc){
         Pvectsz[ptsz->index_md] = ptsz->index_md_cov_Y_Y_ssc;
         int index_multipole_1_multipole_2 = (int) (index_integrand - ptsz->index_integrand_id_cov_Y_Y_ssc_first);
         int n = (-1.+sqrt(1. + 4.*2.*index_multipole_1_multipole_2))/2.;
         int index_multipole_1 = floor(n);
         int index_multipole_2 = index_multipole_1_multipole_2 -index_multipole_1*(index_multipole_1+1)/2;
         Pvectsz[ptsz->index_multipole_1] = (double) (index_multipole_1);
         Pvectsz[ptsz->index_multipole_2] = (double) (index_multipole_2);
         Pvectsz[ptsz->index_has_electron_pressure] = 1;
         if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_cov_Y_Y_ssc_first) printf("computing cov(Y,Y) [ssc] @ ell_id = %.0f, ell_id = %.0f\n",Pvectsz[ptsz->index_multipole_1],Pvectsz[ptsz->index_multipole_2]);
         if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_cov_Y_Y_ssc_last) printf("computing cov(Y,Y) [ssc] @ ell_id = %.0f, ell_id = %.0f\n",Pvectsz[ptsz->index_multipole_1],Pvectsz[ptsz->index_multipole_2]);
         if (ptsz->sz_verbose > 1) printf("computing cov(Y,Y) [ssc] @ ell_id = %.0f, ell_id = %.0f\n",Pvectsz[ptsz->index_multipole_1],Pvectsz[ptsz->index_multipole_2]);
       }

    else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_gal_1h_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_gal_1h_last && ptsz->has_kSZ_kSZ_gal_1h){
       Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gal_1h;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_gal_1h_first);
       Pvectsz[ptsz->index_has_electron_density] = 1;
       Pvectsz[ptsz->index_has_galaxy] = 1;
       if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-gal (1h) @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
    else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_gal_1h_fft_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_gal_1h_fft_last && ptsz->has_kSZ_kSZ_gal_1h_fft){
       Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gal_1h_fft;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_gal_1h_fft_first);
       Pvectsz[ptsz->index_has_electron_density] = 1;
       Pvectsz[ptsz->index_has_galaxy] = 1;
       if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-gal (1h) FFT @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
     else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_gal_2h_fft_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_gal_2h_fft_last && ptsz->has_kSZ_kSZ_gal_2h_fft){
        Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gal_2h_fft;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_gal_2h_fft_first);
        Pvectsz[ptsz->index_has_electron_density] = 1;
        Pvectsz[ptsz->index_has_galaxy] = 1;
        if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-gal (2h) FFT @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
     else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_gal_3h_fft_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_gal_3h_fft_last && ptsz->has_kSZ_kSZ_gal_3h_fft){
        Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gal_3h_fft;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_gal_3h_fft_first);
        Pvectsz[ptsz->index_has_electron_density] = 1;
        Pvectsz[ptsz->index_has_galaxy] = 1;
        if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-gal (3h) FFT @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
     else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_gal_2h_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_gal_2h_last && ptsz->has_kSZ_kSZ_gal_2h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gal_2h;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_gal_2h_first);
        Pvectsz[ptsz->index_has_electron_density] = 1;
        Pvectsz[ptsz->index_has_galaxy] = 1;
        if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-gal (2h) @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_gal_3h_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_gal_3h_last && ptsz->has_kSZ_kSZ_gal_3h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gal_3h;
         Pvectsz[ptsz->index_has_electron_density] = 1;
         Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_gal_3h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-gal (3h) @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
      else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_gal_hf_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_gal_hf_last && ptsz->has_kSZ_kSZ_gal_hf){
         Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gal_hf;
         // Pvectsz[ptsz->index_has_electron_density] = 1;
         // Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_gal_hf_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-gal (hf) @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_last && ptsz->has_kSZ_kSZ_lensmag_1halo){
        Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_lensmag_1halo;
        Pvectsz[ptsz->index_has_electron_density] = 1;
        Pvectsz[ptsz->index_has_lensing] = 1;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-lensmag @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
    else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_gallens_1h_fft_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_gallens_1h_fft_last && ptsz->has_kSZ_kSZ_gallens_1h_fft){
       Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gallens_1h_fft;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_gallens_1h_fft_first);
       Pvectsz[ptsz->index_has_electron_density] = 1;
       Pvectsz[ptsz->index_has_lensing] = 1;
       if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-gallens (1h) FFT @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
     else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_gallens_2h_fft_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_gallens_2h_fft_last && ptsz->has_kSZ_kSZ_gallens_2h_fft){
        Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gallens_2h_fft;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_gallens_2h_fft_first);
        Pvectsz[ptsz->index_has_electron_density] = 1;
        Pvectsz[ptsz->index_has_lensing] = 1;
        if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-gallens (2h) FFT @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
     else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_gallens_3h_fft_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_gallens_3h_fft_last && ptsz->has_kSZ_kSZ_gallens_3h_fft){
        Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gallens_3h_fft;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_gallens_3h_fft_first);
        Pvectsz[ptsz->index_has_electron_density] = 1;
        Pvectsz[ptsz->index_has_lensing] = 1;
        if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-gallens (3h) FFT @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_gallens_hf_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_gallens_hf_last && ptsz->has_kSZ_kSZ_gallens_hf){
         Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gallens_hf;
         // Pvectsz[ptsz->index_has_electron_density] = 1;
         // Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_gallens_hf_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-gallens (hf) @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
    else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_lens_1h_fft_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_lens_1h_fft_last && ptsz->has_kSZ_kSZ_lens_1h_fft){
       Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_lens_1h_fft;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_lens_1h_fft_first);
       Pvectsz[ptsz->index_has_electron_density] = 1;
       Pvectsz[ptsz->index_has_lensing] = 1;
       if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-lens (1h) FFT @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
     else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_lens_2h_fft_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_lens_2h_fft_last && ptsz->has_kSZ_kSZ_lens_2h_fft){
        Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_lens_2h_fft;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_lens_2h_fft_first);
        Pvectsz[ptsz->index_has_electron_density] = 1;
        Pvectsz[ptsz->index_has_lensing] = 1;
        if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-lens (2h) FFT @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
     else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_lens_3h_fft_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_lens_3h_fft_last && ptsz->has_kSZ_kSZ_lens_3h_fft){
        Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_lens_3h_fft;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_lens_3h_fft_first);
        Pvectsz[ptsz->index_has_electron_density] = 1;
        Pvectsz[ptsz->index_has_lensing] = 1;
        if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-lens (3h) FFT @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_lens_hf_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_lens_hf_last && ptsz->has_kSZ_kSZ_lens_hf){
         Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_lens_hf;
         // Pvectsz[ptsz->index_has_electron_density] = 1;
         // Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_lens_hf_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-lens (hf) @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
    else if (index_integrand>=ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_first && index_integrand <= ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_last && ptsz->has_tSZ_tSZ_tSZ_1halo){
       Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_tSZ_tSZ_1halo;
       Pvectsz[ptsz->index_has_electron_pressure] = 1;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_first);
       if (ptsz->sz_verbose > 0) printf("computing b^y-y-y 1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
    else if (index_integrand>=ptsz->index_integrand_id_tSZ_tSZ_tSZ_2h_first && index_integrand <= ptsz->index_integrand_id_tSZ_tSZ_tSZ_2h_last && ptsz->has_tSZ_tSZ_tSZ_2h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_tSZ_tSZ_2h;
        Pvectsz[ptsz->index_has_electron_pressure] = 1;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_tSZ_tSZ_2h_first);
        if (ptsz->sz_verbose > 0) printf("computing b^y-y-y 2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
  else if (index_integrand>=ptsz->index_integrand_id_tSZ_tSZ_tSZ_3h_first && index_integrand <= ptsz->index_integrand_id_tSZ_tSZ_tSZ_3h_last && ptsz->has_tSZ_tSZ_tSZ_3h){
       Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_tSZ_tSZ_3h;
       Pvectsz[ptsz->index_has_electron_pressure] = 1;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_tSZ_tSZ_3h_first);
       if (ptsz->sz_verbose > 0) printf("computing b^y-y-y  3h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
    else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_tSZ_1h_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_tSZ_1h_last && ptsz->has_kSZ_kSZ_tSZ_1h){
       Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_tSZ_1h;
       Pvectsz[ptsz->index_has_electron_pressure] = 1;
       Pvectsz[ptsz->index_has_electron_density] = 1;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_tSZ_1h_first);
       if (ptsz->sz_verbose > 0) printf("computing b^t-t-y 1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
    else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_tSZ_2h_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_tSZ_2h_last && ptsz->has_kSZ_kSZ_tSZ_2h){
       Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_tSZ_2h;
       Pvectsz[ptsz->index_has_electron_pressure] = 1;
       Pvectsz[ptsz->index_has_electron_density] = 1;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_tSZ_2h_first);
       if (ptsz->sz_verbose > 0) printf("computing b^t-t-y 2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
    else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_1h_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_1h_last && ptsz->has_kSZ_kSZ_1h){
       Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_1h;
       Pvectsz[ptsz->index_has_electron_density] = 1;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_1h_first);
       if (ptsz->sz_verbose > 0) printf("computing cl ksz ksz 1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
    else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_2h_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_2h_last && ptsz->has_kSZ_kSZ_2h){
       Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_2h;
       Pvectsz[ptsz->index_has_electron_density] = 1;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_2h_first);
       if (ptsz->sz_verbose > 0) printf("computing cl ksz ksz 2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
     else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_tSZ_3h_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_tSZ_3h_last && ptsz->has_kSZ_kSZ_tSZ_3h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_tSZ_3h;
        Pvectsz[ptsz->index_has_electron_pressure] = 1;
        Pvectsz[ptsz->index_has_electron_density] = 1;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_tSZ_3h_first);
        if (ptsz->sz_verbose > 0) printf("computing b^t-t-y 3h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
     else if (index_integrand>=ptsz->index_integrand_id_pk_at_z_1h_first && index_integrand <= ptsz->index_integrand_id_pk_at_z_1h_last && ptsz->has_pk_at_z_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_pk_at_z_1h;
        Pvectsz[ptsz->index_has_matter_density] = 1;
        Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_at_z_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing pk^1h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_pk_at_z_2h_first && index_integrand <= ptsz->index_integrand_id_pk_at_z_2h_last && ptsz->has_pk_at_z_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_pk_at_z_2h;
         Pvectsz[ptsz->index_has_matter_density] = 1;
         Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_at_z_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing pk^2h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_pk_gg_at_z_1h_first && index_integrand <= ptsz->index_integrand_id_pk_gg_at_z_1h_last && ptsz->has_pk_gg_at_z_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_pk_gg_at_z_1h;
        Pvectsz[ptsz->index_has_galaxy] = 1;
        Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_gg_at_z_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing pk_gg^1h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_pk_gg_at_z_2h_first && index_integrand <= ptsz->index_integrand_id_pk_gg_at_z_2h_last && ptsz->has_pk_gg_at_z_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_pk_gg_at_z_2h;
         Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_gg_at_z_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing pk_gg^2h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_pk_bb_at_z_1h_first && index_integrand <= ptsz->index_integrand_id_pk_bb_at_z_1h_last && ptsz->has_pk_bb_at_z_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_pk_bb_at_z_1h;
        Pvectsz[ptsz->index_has_electron_density] = 1;
        Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_bb_at_z_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing pk_bb^1h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_pk_bb_at_z_2h_first && index_integrand <= ptsz->index_integrand_id_pk_bb_at_z_2h_last && ptsz->has_pk_bb_at_z_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_pk_bb_at_z_2h;
         Pvectsz[ptsz->index_has_electron_density] = 1;
         Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_bb_at_z_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing pk_bb^2h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
       }
      else if (index_integrand>=ptsz->index_integrand_id_pk_b_at_z_2h_first && index_integrand <= ptsz->index_integrand_id_pk_b_at_z_2h_last && ptsz->has_pk_b_at_z_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_pk_b_at_z_2h;
         Pvectsz[ptsz->index_has_electron_density] = 1;
         Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_b_at_z_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing pk_b^2h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_pk_em_at_z_1h_first && index_integrand <= ptsz->index_integrand_id_pk_em_at_z_1h_last && ptsz->has_pk_em_at_z_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_pk_em_at_z_1h;
        Pvectsz[ptsz->index_has_electron_density] = 1;
        Pvectsz[ptsz->index_has_matter_density] = 1;
        Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_em_at_z_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing pk_bb^1h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_pk_em_at_z_2h_first && index_integrand <= ptsz->index_integrand_id_pk_em_at_z_2h_last && ptsz->has_pk_em_at_z_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_pk_em_at_z_2h;
         Pvectsz[ptsz->index_has_electron_density] = 1;
         Pvectsz[ptsz->index_has_matter_density] = 1;
         Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_em_at_z_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing pk_em^2h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_pk_HI_at_z_1h_first && index_integrand <= ptsz->index_integrand_id_pk_HI_at_z_1h_last && ptsz->has_pk_HI_at_z_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_pk_HI_at_z_1h;
        Pvectsz[ptsz->index_has_HI_density] = 1;
        Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_HI_at_z_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing pk_HI^1h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_pk_HI_at_z_2h_first && index_integrand <= ptsz->index_integrand_id_pk_HI_at_z_2h_last && ptsz->has_pk_HI_at_z_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_pk_HI_at_z_2h;
         Pvectsz[ptsz->index_has_HI_density] = 1;
         Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_HI_at_z_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing pk_HI^2h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_bk_at_z_1h_first && index_integrand <= ptsz->index_integrand_id_bk_at_z_1h_last && ptsz->has_bk_at_z_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_bk_at_z_1h;
        Pvectsz[ptsz->index_has_matter_density] = 1;
        Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_bk_at_z_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing bk^1h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_bk_at_z_2h_first && index_integrand <= ptsz->index_integrand_id_bk_at_z_2h_last && ptsz->has_bk_at_z_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_bk_at_z_2h;
        Pvectsz[ptsz->index_has_matter_density] = 1;
         Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_bk_at_z_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing bk^2h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
       }
      else if (index_integrand>=ptsz->index_integrand_id_bk_at_z_3h_first && index_integrand <= ptsz->index_integrand_id_bk_at_z_3h_last && ptsz->has_bk_at_z_3h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_bk_at_z_3h;
         Pvectsz[ptsz->index_has_matter_density] = 1;
         Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_bk_at_z_3h_first);
         if (ptsz->sz_verbose > 0) printf("computing bk^3h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_bk_ttg_at_z_1h_first && index_integrand <= ptsz->index_integrand_id_bk_ttg_at_z_1h_last && ptsz->has_bk_ttg_at_z_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_bk_ttg_at_z_1h;
        Pvectsz[ptsz->index_has_electron_density] = 1;
        Pvectsz[ptsz->index_has_galaxy] = 1;
        Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_bk_ttg_at_z_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing bk^1h (ttg) @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_bk_ttg_at_z_2h_first && index_integrand <= ptsz->index_integrand_id_bk_ttg_at_z_2h_last && ptsz->has_bk_ttg_at_z_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_bk_ttg_at_z_2h;
         Pvectsz[ptsz->index_has_electron_density] = 1;
         Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_bk_ttg_at_z_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing bk^2h (ttg) @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
       }
      else if (index_integrand>=ptsz->index_integrand_id_bk_ttg_at_z_3h_first && index_integrand <= ptsz->index_integrand_id_bk_ttg_at_z_3h_last && ptsz->has_bk_ttg_at_z_3h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_bk_ttg_at_z_3h;
         Pvectsz[ptsz->index_has_electron_density] = 1;
         Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_bk_ttg_at_z_3h_first);
         if (ptsz->sz_verbose > 0) printf("computing bk^3h (ttg) @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_gal_gal_1h_first && index_integrand <= ptsz->index_integrand_id_gal_gal_1h_last && ptsz->has_gal_gal_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_gal_gal_1h;
        Pvectsz[ptsz->index_has_galaxy] = 1;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_gal_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^gal-gal_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
     else if (index_integrand>=ptsz->index_integrand_id_gal_gal_2h_first && index_integrand <= ptsz->index_integrand_id_gal_gal_2h_last && ptsz->has_gal_gal_2h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_gal_gal_2h;
        Pvectsz[ptsz->index_has_galaxy] = 1;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_gal_2h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^gal-gal_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_gal_gal_hf_first && index_integrand <= ptsz->index_integrand_id_gal_gal_hf_last && ptsz->has_gal_gal_hf){
         Pvectsz[ptsz->index_md] = ptsz->index_md_gal_gal_hf;
         // Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_gal_hf_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^gal-gal_hf @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
      else if (index_integrand>=ptsz->index_integrand_id_gal_lens_1h_first && index_integrand <= ptsz->index_integrand_id_gal_lens_1h_last && ptsz->has_gal_lens_1h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_gal_lens_1h;
         Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_has_lensing] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_lens_1h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^gal-lens_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
       else if (index_integrand>=ptsz->index_integrand_id_gal_lens_2h_first && index_integrand <= ptsz->index_integrand_id_gal_lens_2h_last && ptsz->has_gal_lens_2h){
          Pvectsz[ptsz->index_md] = ptsz->index_md_gal_lens_2h;
          Pvectsz[ptsz->index_has_galaxy] = 1;
          Pvectsz[ptsz->index_has_lensing] = 1;
          Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_lens_2h_first);
          if (ptsz->sz_verbose > 0) printf("computing cl^gal-lens_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
        }
       else if (index_integrand>=ptsz->index_integrand_id_gal_lens_hf_first && index_integrand <= ptsz->index_integrand_id_gal_lens_hf_last && ptsz->has_gal_lens_hf){
          Pvectsz[ptsz->index_md] = ptsz->index_md_gal_lens_hf;
          // Pvectsz[ptsz->index_has_galaxy] = 1;
          // Pvectsz[ptsz->index_has_lensing] = 1;
          Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_lens_hf_first);
          if (ptsz->sz_verbose > 0) printf("computing cl^gal-lens_hf @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
        }
      else if (index_integrand>=ptsz->index_integrand_id_gal_lensmag_1h_first && index_integrand <= ptsz->index_integrand_id_gal_lensmag_1h_last && ptsz->has_gal_lensmag_1h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_gal_lensmag_1h;
         Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_has_lensing] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_lensmag_1h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^gal-lensmag_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
       else if (index_integrand>=ptsz->index_integrand_id_gal_lensmag_2h_first && index_integrand <= ptsz->index_integrand_id_gal_lensmag_2h_last && ptsz->has_gal_lensmag_2h){
          Pvectsz[ptsz->index_md] = ptsz->index_md_gal_lensmag_2h;
          Pvectsz[ptsz->index_has_galaxy] = 1;
          Pvectsz[ptsz->index_has_lensing] = 1;
          Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_lensmag_2h_first);
          if (ptsz->sz_verbose > 0) printf("computing cl^gal-lensmag_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
        }
      else if (index_integrand>=ptsz->index_integrand_id_gal_gallens_1h_first && index_integrand <= ptsz->index_integrand_id_gal_gallens_1h_last && ptsz->has_gal_gallens_1h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_gal_gallens_1h;
         Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_has_lensing] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_gallens_1h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^gal-gallens_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
       else if (index_integrand>=ptsz->index_integrand_id_gal_gallens_2h_first && index_integrand <= ptsz->index_integrand_id_gal_gallens_2h_last && ptsz->has_gal_gallens_2h){
          Pvectsz[ptsz->index_md] = ptsz->index_md_gal_gallens_2h;
          Pvectsz[ptsz->index_has_galaxy] = 1;
          Pvectsz[ptsz->index_has_lensing] = 1;
          Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_gallens_2h_first);
          if (ptsz->sz_verbose > 0) printf("computing cl^gal-gallens_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
        }
      else if (index_integrand>=ptsz->index_integrand_id_gallens_gallens_1h_first && index_integrand <= ptsz->index_integrand_id_gallens_gallens_1h_last && ptsz->has_gallens_gallens_1h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_gallens_gallens_1h;
         // Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_has_lensing] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gallens_gallens_1h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^gallens-gallens_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
       else if (index_integrand>=ptsz->index_integrand_id_gallens_gallens_2h_first && index_integrand <= ptsz->index_integrand_id_gallens_gallens_2h_last && ptsz->has_gallens_gallens_2h){
          Pvectsz[ptsz->index_md] = ptsz->index_md_gallens_gallens_2h;
          // Pvectsz[ptsz->index_has_galaxy] = 1;
          Pvectsz[ptsz->index_has_lensing] = 1;
          Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gallens_gallens_2h_first);
          if (ptsz->sz_verbose > 0) printf("computing cl^gallens-gallens_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
        }
      else if (index_integrand>=ptsz->index_integrand_id_gallens_lens_1h_first && index_integrand <= ptsz->index_integrand_id_gallens_lens_1h_last && ptsz->has_gallens_lens_1h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_gallens_lens_1h;
         // Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_has_lensing] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gallens_lens_1h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^gallens-lens_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
       else if (index_integrand>=ptsz->index_integrand_id_gallens_lens_2h_first && index_integrand <= ptsz->index_integrand_id_gallens_lens_2h_last && ptsz->has_gallens_lens_2h){
          Pvectsz[ptsz->index_md] = ptsz->index_md_gallens_lens_2h;
          // Pvectsz[ptsz->index_has_galaxy] = 1;
          Pvectsz[ptsz->index_has_lensing] = 1;
          Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gallens_lens_2h_first);
          if (ptsz->sz_verbose > 0) printf("computing cl^gallens-lens_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
        }
       else if (index_integrand>=ptsz->index_integrand_id_gal_lensmag_hf_first && index_integrand <= ptsz->index_integrand_id_gal_lensmag_hf_last && ptsz->has_gal_lensmag_hf){
          Pvectsz[ptsz->index_md] = ptsz->index_md_gal_lensmag_hf;
          // Pvectsz[ptsz->index_has_galaxy] = 1;
          // Pvectsz[ptsz->index_has_lensing] = 1;
          Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_lensmag_hf_first);
          if (ptsz->sz_verbose > 0) printf("computing cl^gal-lensmag_hf @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
        }
      else if (index_integrand>=ptsz->index_integrand_id_lensmag_lensmag_1h_first && index_integrand <= ptsz->index_integrand_id_lensmag_lensmag_1h_last && ptsz->has_lensmag_lensmag_1h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_lensmag_lensmag_1h;
         Pvectsz[ptsz->index_has_lensing] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_lensmag_lensmag_1h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^lensmag-lensmag_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
      else if (index_integrand>=ptsz->index_integrand_id_lensmag_lensmag_hf_first && index_integrand <= ptsz->index_integrand_id_lensmag_lensmag_hf_last && ptsz->has_lensmag_lensmag_hf){
         Pvectsz[ptsz->index_md] = ptsz->index_md_lensmag_lensmag_hf;
         // Pvectsz[ptsz->index_has_galaxy] = 1;
         // Pvectsz[ptsz->index_has_lensing] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_lensmag_lensmag_hf_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^lensmag-lensmag_hf @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
       else if (index_integrand>=ptsz->index_integrand_id_lensmag_lensmag_2h_first && index_integrand <= ptsz->index_integrand_id_lensmag_lensmag_2h_last && ptsz->has_lensmag_lensmag_2h){
          Pvectsz[ptsz->index_md] = ptsz->index_md_lensmag_lensmag_2h;
          Pvectsz[ptsz->index_has_lensing] = 1;
          Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_lensmag_lensmag_2h_first);
          if (ptsz->sz_verbose > 0) printf("computing cl^lensmag-lensmag_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
        }
        else if (index_integrand>=ptsz->index_integrand_id_lens_lensmag_1h_first && index_integrand <= ptsz->index_integrand_id_lens_lensmag_1h_last && ptsz->has_lens_lensmag_1h){
           Pvectsz[ptsz->index_md] = ptsz->index_md_lens_lensmag_1h;
           Pvectsz[ptsz->index_has_lensing] = 1;
           Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_lens_lensmag_1h_first);
           if (ptsz->sz_verbose > 0) printf("computing cl^lens-lensmag_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
         }
       else if (index_integrand>=ptsz->index_integrand_id_lens_lensmag_2h_first && index_integrand <= ptsz->index_integrand_id_lens_lensmag_2h_last && ptsz->has_lens_lensmag_2h){
          Pvectsz[ptsz->index_md] = ptsz->index_md_lens_lensmag_2h;
          Pvectsz[ptsz->index_has_lensing] = 1;
          Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_lens_lensmag_2h_first);
          if (ptsz->sz_verbose > 0) printf("computing cl^lens-lensmag_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
        }
       else if (index_integrand>=ptsz->index_integrand_id_lens_lensmag_hf_first && index_integrand <= ptsz->index_integrand_id_lens_lensmag_hf_last && ptsz->has_lens_lensmag_hf){
          Pvectsz[ptsz->index_md] = ptsz->index_md_lens_lensmag_hf;
          // Pvectsz[ptsz->index_has_lensing] = 1;
          Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_lens_lensmag_hf_first);
          if (ptsz->sz_verbose > 0) printf("computing cl^lens-lensmag_hf @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
        }
      else if (index_integrand>=ptsz->index_integrand_id_tSZ_cib_1h_first && index_integrand <= ptsz->index_integrand_id_tSZ_cib_1h_last && ptsz->has_tSZ_cib_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_cib_1h;
        int index_multipole_cib1 = (double) (index_integrand - ptsz->index_integrand_id_tSZ_cib_1h_first);
        int index_cib1 = index_multipole_cib1 / ptsz->nlSZ;
        int index_multipole = index_multipole_cib1 % ptsz->nlSZ;
        Pvectsz[ptsz->index_multipole] = (double) index_multipole ;
        Pvectsz[ptsz->index_frequency_for_cib_profile] = (double) index_cib1;
        Pvectsz[ptsz->index_has_cib] = 1;
        Pvectsz[ptsz->index_has_electron_pressure] = 1;
        if (ptsz->sz_verbose > 0) printf("computing cl^y-cib_1h @ frequency_id = %.0f, ell_id = %.0f\n",
                                         Pvectsz[ptsz->index_frequency_for_cib_profile],
                                         Pvectsz[ptsz->index_multipole]);
        }
       else if (index_integrand>=ptsz->index_integrand_id_tSZ_cib_2h_first && index_integrand <= ptsz->index_integrand_id_tSZ_cib_2h_last && ptsz->has_tSZ_cib_2h){
          Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_cib_2h;
          int index_multipole_cib1 = (double) (index_integrand - ptsz->index_integrand_id_tSZ_cib_2h_first);
          int index_cib1 = index_multipole_cib1 / ptsz->nlSZ;
          int index_multipole = index_multipole_cib1 % ptsz->nlSZ;
          Pvectsz[ptsz->index_multipole] = (double) index_multipole ;
          Pvectsz[ptsz->index_frequency_for_cib_profile] = (double) index_cib1;
          Pvectsz[ptsz->index_has_cib] = 1;
          Pvectsz[ptsz->index_has_electron_pressure] = 1;
          if (ptsz->sz_verbose > 0) printf("computing cl^y-cib_2h @ frequency_id = %.0f, ell_id = %.0f\n",
                                           Pvectsz[ptsz->index_frequency_for_cib_profile],
                                           Pvectsz[ptsz->index_multipole]);
        }
      else if (index_integrand>=ptsz->index_integrand_id_lens_cib_1h_first && index_integrand <= ptsz->index_integrand_id_lens_cib_1h_last && ptsz->has_lens_cib_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_lens_cib_1h;
        int index_multipole_cib1 = (double) (index_integrand - ptsz->index_integrand_id_lens_cib_1h_first);
        int index_cib1 = index_multipole_cib1 / ptsz->nlSZ;
        int index_multipole = index_multipole_cib1 % ptsz->nlSZ;
        Pvectsz[ptsz->index_multipole] = (double) index_multipole ;
        Pvectsz[ptsz->index_frequency_for_cib_profile] = (double) index_cib1;
        Pvectsz[ptsz->index_has_cib] = 1;
        Pvectsz[ptsz->index_has_lensing] = 1;
        if (ptsz->sz_verbose > 0) printf("computing cl^lens-cib_1h @ frequency_id = %.0f, ell_id = %.0f\n",
                                         Pvectsz[ptsz->index_frequency_for_cib_profile],
                                         Pvectsz[ptsz->index_multipole]);
        }
       else if (index_integrand>=ptsz->index_integrand_id_lens_cib_2h_first && index_integrand <= ptsz->index_integrand_id_lens_cib_2h_last && ptsz->has_lens_cib_2h){
          Pvectsz[ptsz->index_md] = ptsz->index_md_lens_cib_2h;
          int index_multipole_cib1 = (double) (index_integrand - ptsz->index_integrand_id_lens_cib_2h_first);
          int index_cib1 = index_multipole_cib1 / ptsz->nlSZ;
          int index_multipole = index_multipole_cib1 % ptsz->nlSZ;
          Pvectsz[ptsz->index_multipole] = (double) index_multipole ;
          Pvectsz[ptsz->index_frequency_for_cib_profile] = (double) index_cib1;
          Pvectsz[ptsz->index_has_cib] = 1;
          Pvectsz[ptsz->index_has_lensing] = 1;
          if (ptsz->sz_verbose > 0) printf("computing cl^lens-cib_2h @ frequency_id = %.0f, ell_id = %.0f\n",
                                           Pvectsz[ptsz->index_frequency_for_cib_profile],
                                           Pvectsz[ptsz->index_multipole]);
        }
      else if (index_integrand>=ptsz->index_integrand_id_gal_cib_1h_first && index_integrand <= ptsz->index_integrand_id_gal_cib_1h_last && ptsz->has_gal_cib_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_gal_cib_1h;
        int index_multipole_cib1 = (double) (index_integrand - ptsz->index_integrand_id_gal_cib_1h_first);
        int index_cib1 = index_multipole_cib1 / ptsz->nlSZ;
        int index_multipole = index_multipole_cib1 % ptsz->nlSZ;
        Pvectsz[ptsz->index_multipole] = (double) index_multipole ;
        Pvectsz[ptsz->index_frequency_for_cib_profile] = (double) index_cib1;
        Pvectsz[ptsz->index_has_cib] = 1;
        Pvectsz[ptsz->index_has_galaxy] = 1;
        if (ptsz->sz_verbose > 0) printf("computing cl^gal-cib_1h @ frequency_id = %.0f, ell_id = %.0f\n",
                                         Pvectsz[ptsz->index_frequency_for_cib_profile],
                                         Pvectsz[ptsz->index_multipole]);
        }
       else if (index_integrand>=ptsz->index_integrand_id_gal_cib_2h_first && index_integrand <= ptsz->index_integrand_id_gal_cib_2h_last && ptsz->has_gal_cib_2h){
          Pvectsz[ptsz->index_md] = ptsz->index_md_gal_cib_2h;
          int index_multipole_cib1 = (double) (index_integrand - ptsz->index_integrand_id_gal_cib_2h_first);
          int index_cib1 = index_multipole_cib1 / ptsz->nlSZ;
          int index_multipole = index_multipole_cib1 % ptsz->nlSZ;
          Pvectsz[ptsz->index_multipole] = (double) index_multipole ;
          Pvectsz[ptsz->index_frequency_for_cib_profile] = (double) index_cib1;
          Pvectsz[ptsz->index_has_cib] = 1;
          Pvectsz[ptsz->index_has_galaxy] = 1;
          if (ptsz->sz_verbose > 0) printf("computing cl^gal-cib_2h @ frequency_id = %.0f, ell_id = %.0f\n",
                                           Pvectsz[ptsz->index_frequency_for_cib_profile],
                                           Pvectsz[ptsz->index_multipole]);
        }
      else if (index_integrand>=ptsz->index_integrand_id_cib_cib_1h_first && index_integrand <= ptsz->index_integrand_id_cib_cib_1h_last && ptsz->has_cib_cib_1h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_cib_cib_1h;
         int index_multipole_cib1_cib2 = (int) (index_integrand - ptsz->index_integrand_id_cib_cib_1h_first);
         int index_multipole = index_multipole_cib1_cib2 % ptsz->nlSZ;
         int index_cib1_cib2 = index_multipole_cib1_cib2 / ptsz->nlSZ;
         int n = (-1.+sqrt(1. + 4.*2.*index_cib1_cib2))/2.;
         int index_cib1 = floor(n);
         int index_cib2 = index_cib1_cib2 -index_cib1*(index_cib1+1)/2;
         Pvectsz[ptsz->index_frequency_for_cib_profile] = (double) index_cib1;
         Pvectsz[ptsz->index_frequency_prime_for_cib_profile] = (double) index_cib2;
         //int index_multipole = (int) (index_integrand - ptsz->index_integrand_id_cib_cib_1h_first);
         Pvectsz[ptsz->index_multipole] = (double) index_multipole;
         Pvectsz[ptsz->index_has_cib] = 1;
         if (ptsz->sz_verbose > 0) printf("computing cl^cib-cib_1h @ frequency_id = %.0f, frequency_prime_id = %.0f, ell_id = %.0f\n",
                                          Pvectsz[ptsz->index_frequency_for_cib_profile],
                                          Pvectsz[ptsz->index_frequency_prime_for_cib_profile],
                                          Pvectsz[ptsz->index_multipole]);
       }
      else if (index_integrand>=ptsz->index_integrand_id_ngal_ngal_1h_first && index_integrand <= ptsz->index_integrand_id_ngal_ngal_1h_last && ptsz->has_ngal_ngal_1h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_ngal_ngal_1h;
         int index_multipole_ngal1_ngal2 = (int) (index_integrand - ptsz->index_integrand_id_ngal_ngal_1h_first);
         int index_multipole = index_multipole_ngal1_ngal2 % ptsz->nlSZ;
         int index_ngal1_ngal2 = index_multipole_ngal1_ngal2 / ptsz->nlSZ;
         int n = (-1.+sqrt(1. + 4.*2.*index_ngal1_ngal2))/2.;
         int index_ngal1 = floor(n);
         int index_ngal2 = index_ngal1_ngal2 -index_ngal1*(index_ngal1+1)/2;
         Pvectsz[ptsz->index_ngal_for_galaxy_profile] = (double) index_ngal1;
         Pvectsz[ptsz->index_ngal_prime_for_galaxy_profile] = (double) index_ngal2;
         //int index_multipole = (int) (index_integrand - ptsz->index_integrand_id_cib_cib_1h_first);
         Pvectsz[ptsz->index_multipole] = (double) index_multipole;
         Pvectsz[ptsz->index_has_galaxy] = 1;
         if (ptsz->sz_verbose > 0) printf("computing cl^ngal-ngal_1h @ ngal_id = %.0f, ngal_prime_id = %.0f, ell_id = %.0f\n",
                                          Pvectsz[ptsz->index_ngal_for_galaxy_profile],
                                          Pvectsz[ptsz->index_ngal_prime_for_galaxy_profile],
                                          Pvectsz[ptsz->index_multipole]);
       }

      else if (index_integrand>=ptsz->index_integrand_id_ngal_ngal_2h_first && index_integrand <= ptsz->index_integrand_id_ngal_ngal_2h_last && ptsz->has_ngal_ngal_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_ngal_ngal_2h;
         int index_multipole_ngal1_ngal2 = (int) (index_integrand - ptsz->index_integrand_id_ngal_ngal_2h_first);
         int index_multipole = index_multipole_ngal1_ngal2 % ptsz->nlSZ;
         int index_ngal1_ngal2 = index_multipole_ngal1_ngal2 / ptsz->nlSZ;
         int n = (-1.+sqrt(1. + 4.*2.*index_ngal1_ngal2))/2.;
         int index_ngal1 = floor(n);
         int index_ngal2 = index_ngal1_ngal2 -index_ngal1*(index_ngal1+1)/2;
         Pvectsz[ptsz->index_ngal_for_galaxy_profile] = (double) index_ngal1;
         Pvectsz[ptsz->index_ngal_prime_for_galaxy_profile] = (double) index_ngal2;
         //int index_multipole = (int) (index_integrand - ptsz->index_integrand_id_cib_cib_1h_first);
         Pvectsz[ptsz->index_multipole] = (double) index_multipole;
         Pvectsz[ptsz->index_has_galaxy] = 1;
         if (ptsz->sz_verbose > 0) printf("computing cl^ngal-ngal_2h @ ngal_id = %.0f, ngal_prime_id = %.0f, ell_id = %.0f\n",
                                          Pvectsz[ptsz->index_ngal_for_galaxy_profile],
                                          Pvectsz[ptsz->index_ngal_prime_for_galaxy_profile],
                                          Pvectsz[ptsz->index_multipole]);
       }

      else if (index_integrand>=ptsz->index_integrand_id_ngal_ngal_hf_first && index_integrand <= ptsz->index_integrand_id_ngal_ngal_hf_last && ptsz->has_ngal_ngal_hf){
         Pvectsz[ptsz->index_md] = ptsz->index_md_ngal_ngal_hf;
         int index_multipole_ngal1_ngal2 = (int) (index_integrand - ptsz->index_integrand_id_ngal_ngal_hf_first);
         int index_multipole = index_multipole_ngal1_ngal2 % ptsz->nlSZ;
         int index_ngal1_ngal2 = index_multipole_ngal1_ngal2 / ptsz->nlSZ;
         int n = (-1.+sqrt(1. + 4.*2.*index_ngal1_ngal2))/2.;
         int index_ngal1 = floor(n);
         int index_ngal2 = index_ngal1_ngal2 -index_ngal1*(index_ngal1+1)/2;
         Pvectsz[ptsz->index_ngal_for_galaxy_profile] = (double) index_ngal1;
         Pvectsz[ptsz->index_ngal_prime_for_galaxy_profile] = (double) index_ngal2;
         //int index_multipole = (int) (index_integrand - ptsz->index_integrand_id_cib_cib_1h_first);
         Pvectsz[ptsz->index_multipole] = (double) index_multipole;
         // Pvectsz[ptsz->index_has_galaxy] = 1;
         if (ptsz->sz_verbose > 0) printf("computing cl^ngal-ngal_hf @ ngal_id = %.0f, ngal_prime_id = %.0f, ell_id = %.0f\n",
                                          Pvectsz[ptsz->index_ngal_for_galaxy_profile],
                                          Pvectsz[ptsz->index_ngal_prime_for_galaxy_profile],
                                          Pvectsz[ptsz->index_multipole]);
       }
      else if (index_integrand>=ptsz->index_integrand_id_ngal_lens_1h_first && index_integrand <= ptsz->index_integrand_id_ngal_lens_1h_last && ptsz->has_ngal_lens_1h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_ngal_lens_1h;
         int index_multipole_ngal1 = (int) (index_integrand - ptsz->index_integrand_id_ngal_lens_1h_first);
         int index_ngal1 = index_multipole_ngal1 / ptsz->nlSZ;
         int index_multipole = index_multipole_ngal1 % ptsz->nlSZ;
         Pvectsz[ptsz->index_ngal_for_galaxy_profile] = (double) index_ngal1;
         Pvectsz[ptsz->index_multipole] = (double) index_multipole;
         Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_has_lensing] = 1;
         if (ptsz->sz_verbose > 0) printf("computing cl^ngal-lens_1h @ ngal_id = %.0f, ell_id = %.0f\n",
                                          Pvectsz[ptsz->index_ngal_for_galaxy_profile],
                                          Pvectsz[ptsz->index_multipole]);

       }

      else if (index_integrand>=ptsz->index_integrand_id_ngal_lens_2h_first && index_integrand <= ptsz->index_integrand_id_ngal_lens_2h_last && ptsz->has_ngal_lens_2h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_ngal_lens_2h;
        int index_multipole_ngal1 = (int) (index_integrand - ptsz->index_integrand_id_ngal_lens_2h_first);
        int index_ngal1 = index_multipole_ngal1 / ptsz->nlSZ;
        int index_multipole = index_multipole_ngal1 % ptsz->nlSZ;
        Pvectsz[ptsz->index_ngal_for_galaxy_profile] = (double) index_ngal1;
        Pvectsz[ptsz->index_multipole] = (double) index_multipole;
        Pvectsz[ptsz->index_has_galaxy] = 1;
        Pvectsz[ptsz->index_has_lensing] = 1;
        if (ptsz->sz_verbose > 0) printf("computing cl^ngal-lens_2h @ ngal_id = %.0f, ell_id = %.0f\n",
                                         Pvectsz[ptsz->index_ngal_for_galaxy_profile],
                                         Pvectsz[ptsz->index_multipole]);
       }

      else if (index_integrand>=ptsz->index_integrand_id_ngal_lens_hf_first && index_integrand <= ptsz->index_integrand_id_ngal_lens_hf_last && ptsz->has_ngal_lens_hf){
        Pvectsz[ptsz->index_md] = ptsz->index_md_ngal_lens_hf;
        int index_multipole_ngal1 = (int) (index_integrand - ptsz->index_integrand_id_ngal_lens_hf_first);
        int index_ngal1 = index_multipole_ngal1 / ptsz->nlSZ;
        int index_multipole = index_multipole_ngal1 % ptsz->nlSZ;
        Pvectsz[ptsz->index_ngal_for_galaxy_profile] = (double) index_ngal1;
        Pvectsz[ptsz->index_multipole] = (double) index_multipole;
        // Pvectsz[ptsz->index_has_galaxy] = 1;
        if (ptsz->sz_verbose > 0) printf("computing cl^ngal-lens_hf @ ngal_id = %.0f, ell_id = %.0f\n",
                                          Pvectsz[ptsz->index_ngal_for_galaxy_profile],
                                          Pvectsz[ptsz->index_multipole]);

       }

      else if (index_integrand>=ptsz->index_integrand_id_cib_cib_2h_first && index_integrand <= ptsz->index_integrand_id_cib_cib_2h_last && ptsz->has_cib_cib_2h){
          //Pvectsz[ptsz->index_md] = ptsz->index_md_cib_cib_2h;
          //Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_cib_cib_2h_first);
         Pvectsz[ptsz->index_md] = ptsz->index_md_cib_cib_2h;
         int index_multipole_cib1_cib2 = (int) (index_integrand - ptsz->index_integrand_id_cib_cib_2h_first);
         int index_multipole = index_multipole_cib1_cib2 % ptsz->nlSZ;
         int index_cib1_cib2 = index_multipole_cib1_cib2 / ptsz->nlSZ;
         int n = (-1.+sqrt(1. + 4.*2.*index_cib1_cib2))/2.;
         int index_cib1 = floor(n);
         int index_cib2 = index_cib1_cib2 -index_cib1*(index_cib1+1)/2;
         Pvectsz[ptsz->index_frequency_for_cib_profile] = (double) index_cib1;
         Pvectsz[ptsz->index_frequency_prime_for_cib_profile] = (double) index_cib2;
         //int index_multipole = (int) (index_integrand - ptsz->index_integrand_id_cib_cib_1h_first);
         Pvectsz[ptsz->index_multipole] = (double) index_multipole;
         Pvectsz[ptsz->index_has_cib] = 1;
         //  if (ptsz->sz_verbose > 0) printf("computing cl^cib-cib_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
         if (ptsz->sz_verbose > 0) printf("computing cl^cib-cib_2h @ frequency_id = %.0f, frequency_prime_id = %.0f, ell_id = %.0f\n",
                                Pvectsz[ptsz->index_frequency_for_cib_profile],
                                Pvectsz[ptsz->index_frequency_prime_for_cib_profile],
                                Pvectsz[ptsz->index_multipole]);
        }
     else if (index_integrand>=ptsz->index_integrand_id_tSZ_gal_1h_first && index_integrand <= ptsz->index_integrand_id_tSZ_gal_1h_last && ptsz->has_tSZ_gal_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_gal_1h;
        Pvectsz[ptsz->index_has_electron_pressure] = 1;
        Pvectsz[ptsz->index_has_galaxy] = 1;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_gal_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^y-gal_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_tSZ_gal_2h_first && index_integrand <= ptsz->index_integrand_id_tSZ_gal_2h_last && ptsz->has_tSZ_gal_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_gal_2h;
         Pvectsz[ptsz->index_has_electron_pressure] = 1;
         Pvectsz[ptsz->index_has_galaxy] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_gal_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^y-gal_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }

      else if (index_integrand>=ptsz->index_integrand_id_tSZ_lensmag_1h_first && index_integrand <= ptsz->index_integrand_id_tSZ_lensmag_1h_last && ptsz->has_tSZ_lensmag_1h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_lensmag_1h;
         Pvectsz[ptsz->index_has_electron_pressure] = 1;
         Pvectsz[ptsz->index_has_lensing] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_lensmag_1h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^tSZ-lensmag_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
       else if (index_integrand>=ptsz->index_integrand_id_tSZ_lensmag_2h_first && index_integrand <= ptsz->index_integrand_id_tSZ_lensmag_2h_last && ptsz->has_tSZ_lensmag_2h){
          Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_lensmag_2h;
          Pvectsz[ptsz->index_has_electron_pressure] = 1;
          Pvectsz[ptsz->index_has_lensing] = 1;
          Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_lensmag_2h_first);
          if (ptsz->sz_verbose > 0) printf("computing cl^tSZ-lensmag_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
        }
     else if (index_integrand>=ptsz->index_integrand_id_lens_lens_1h_first && index_integrand <= ptsz->index_integrand_id_lens_lens_1h_last && ptsz->has_lens_lens_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_lens_lens_1h;
        Pvectsz[ptsz->index_has_lensing] = 1;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_lens_lens_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^lens-lens_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_lens_lens_2h_first && index_integrand <= ptsz->index_integrand_id_lens_lens_2h_last && ptsz->has_lens_lens_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_lens_lens_2h;
         Pvectsz[ptsz->index_has_lensing] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_lens_lens_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^lens-lens_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
       else if (index_integrand>=ptsz->index_integrand_id_lens_lens_hf_first && index_integrand <= ptsz->index_integrand_id_lens_lens_hf_last && ptsz->has_lens_lens_hf){
          Pvectsz[ptsz->index_md] = ptsz->index_md_lens_lens_hf;
          Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_lens_lens_hf_first);
          if (ptsz->sz_verbose > 0) printf("computing cl^lens-lens_hf @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
        }
     else if (index_integrand>=ptsz->index_integrand_id_tSZ_lens_1h_first && index_integrand <= ptsz->index_integrand_id_tSZ_lens_1h_last && ptsz->has_tSZ_lens_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_lens_1h;
        Pvectsz[ptsz->index_has_lensing] = 1;
        Pvectsz[ptsz->index_has_electron_pressure] = 1;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_lens_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^y-phi_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
     else if (index_integrand>=ptsz->index_integrand_id_tSZ_lens_2h_first && index_integrand <= ptsz->index_integrand_id_tSZ_lens_2h_last && ptsz->has_tSZ_lens_2h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_lens_2h;
        Pvectsz[ptsz->index_has_lensing] = 1;
        Pvectsz[ptsz->index_has_electron_pressure] = 1;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_lens_2h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^y-phi_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_isw_lens_first && index_integrand <= ptsz->index_integrand_id_isw_lens_last && ptsz->has_isw_lens){
         Pvectsz[ptsz->index_md] = ptsz->index_md_isw_lens;
         Pvectsz[ptsz->index_has_lensing] = 1;
         Pvectsz[ptsz->index_has_isw] = 1;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_isw_lens_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^isw-phi @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_isw_tsz_first && index_integrand <= ptsz->index_integrand_id_isw_tsz_last && ptsz->has_isw_tsz){
        Pvectsz[ptsz->index_md] = ptsz->index_md_isw_tsz;
        Pvectsz[ptsz->index_has_isw] = 1;
        Pvectsz[ptsz->index_has_electron_pressure] = 1;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_isw_tsz_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^isw-y @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
     else if (index_integrand>=ptsz->index_integrand_id_isw_auto_first && index_integrand <= ptsz->index_integrand_id_isw_auto_last && ptsz->has_isw_auto){
        Pvectsz[ptsz->index_md] = ptsz->index_md_isw_auto;
        Pvectsz[ptsz->index_has_isw] = 1;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_isw_auto_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^isw-isw @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }


     else
     {
       // printf("id not found. index_integrand = %d \n",index_integrand);
     return _SUCCESS_;
}



 //return 0;

 int index_md = (int) Pvectsz[ptsz->index_md];


 if (_dndlnM_){
   int index_z = (int) Pvectsz[ptsz->index_redshift_for_dndlnM];
   int index_m = (int) Pvectsz[ptsz->index_mass_for_dndlnM];

   double z_asked = ptsz->dndlnM_array_z[index_z];
   double m_asked = ptsz->dndlnM_array_m[index_m];

  ptsz->dndlnM_at_z_and_M[index_z][index_m] = get_dndlnM_at_z_and_M(z_asked,m_asked,ptsz);


 }
 else {

   if (_kSZ_kSZ_gal_1h_
    || _kSZ_kSZ_lensmag_1halo_
    || _kSZ_kSZ_gal_2h_
    || _kSZ_kSZ_gal_3h_
    || _kSZ_kSZ_gal_hf_
    || _kSZ_kSZ_gallens_hf_
    || _kSZ_kSZ_lens_hf_){
     // loop over l1,l2 for each ell
     int index_ell_1 = 0;
     int index_theta_1 = 0;
     int index_ell_2 = 0;
     int index_ell_3 = (int) Pvectsz[ptsz->index_multipole]; // ell_3 of the bispectrum is always ell of the power spectrum


    // put bispectrum in 1d format for 2d interpolation
    int index_l1_l2 = 0;
    double * b_l1_l2_l_1d;
    // double * ln_ell;
    class_alloc(b_l1_l2_l_1d,
                sizeof(double *)*ptsz->N_kSZ2_gal_theta_grid*ptsz->N_kSZ2_gal_multipole_grid,
                ptsz->error_message);
    // class_alloc(ln_ell,
    //             sizeof(double *)*ptsz->N_kSZ2_gal_multipole_grid,
    //             ptsz->error_message);


 for (index_ell_2=0;index_ell_2<ptsz->N_kSZ2_gal_multipole_grid;index_ell_2++){
//  // // ln_ell[index_ell_2] = log(ptsz->ell_kSZ2_gal_multipole_grid[index_ell_2]);
      for (index_theta_1=0;index_theta_1<ptsz->N_kSZ2_gal_theta_grid;index_theta_1++){
//
        if (ptsz->sz_verbose > 100){
         if (_kSZ_kSZ_gal_1h_)
          printf("computing b_kSZ_kSZ_g_1h @ l3_id = %d and (th1,l2) = (%d,%d)\n", index_ell_3, index_theta_1, index_ell_2);
         if (_kSZ_kSZ_gal_2h_)
          printf("computing b_kSZ_kSZ_g_2h @ l3_id = %d and (th1,l2) = (%d,%d)\n", index_ell_3, index_theta_1, index_ell_2);
         if (_kSZ_kSZ_gal_3h_)
          printf("computing b_kSZ_kSZ_g_3h @ l3_id = %d and (th1,l2) = (%d,%d)\n", index_ell_3, index_theta_1, index_ell_2);
         if (_kSZ_kSZ_gal_hf_)
          printf("computing b_kSZ_kSZ_g_hf @ l3_id = %d and (th1,l2) = (%d,%d)\n", index_ell_3, index_theta_1, index_ell_2);
         if (_kSZ_kSZ_gallens_hf_)
          printf("computing b_kSZ_kSZ_kg_hf @ l3_id = %d and (th1,l2) = (%d,%d)\n", index_ell_3, index_theta_1, index_ell_2);
         if (_kSZ_kSZ_lens_hf_)
           printf("computing b_kSZ_kSZ_kcmb_hf @ l3_id = %d and (th1,l2) = (%d,%d)\n", index_ell_3, index_theta_1, index_ell_2);
         if (_kSZ_kSZ_lensmag_1halo_)
          printf("computing b_kSZ_kSZ_mu_1h @ l3_id = %d and (th1,l2) = (%d,%d)\n", index_ell_3, index_theta_1, index_ell_2);
          }
          Pvectsz[ptsz->index_multipole_1] = index_theta_1;
          Pvectsz[ptsz->index_multipole_2] = index_ell_2;
          Pvectsz[ptsz->index_multipole_3] = index_ell_3;

          class_call(integrate_over_redshift(pba,
                                             pnl,
                                             ppm,
                                             ppt,
                                             ptsz,
                                             Pvecback,
                                             Pvectsz),
                          ptsz->error_message,
                          ptsz->error_message);
//
       b_l1_l2_l_1d[index_l1_l2] = Pvectsz[ptsz->index_integral];
//
//

double db = b_l1_l2_l_1d[index_l1_l2];
// printf("db = %.5e\n",db);
if (isnan(db) || isinf(db)){
  // db = 0.;
if (isnan(db))
printf("found nan in grid of b_l1_l2_l_1d\n");
if (isinf(db))
printf("found inf in grid of b_l1_l2_l_1d\n");
printf("id_theta = %d \t id_l2 = %d \n",index_theta_1,index_ell_2);

printf("\n\n");
exit(0);
}

       index_l1_l2 += 1;
//

// printf("ok %d %d\n",index_ell_2,index_theta_1);
       }


     }
   // now we integrate the bispectrum to compute power spectrum
   double cl_kSZ2_gal = 0.;
   double r = 0.; //result of the integral

   struct Parameters_for_integrand_kSZ2_X V;
   V.pnl= pnl;
   V.ppm= ppm;
   V.ptsz = ptsz;
   V.pba = pba;
   V.Pvecback = Pvecback;
   V.Pvectsz = Pvectsz;
   V.ln_ell = ptsz->ell_kSZ2_gal_multipole_grid;//ln_ell;
   V.index_ell_3 = index_ell_3;
   V.b_l1_l2_l_1d = b_l1_l2_l_1d;
   void * params;


   double epsrel= 1.e-6;//ptsz->redshift_epsrel;//ptsz->patterson_epsrel;
   double epsabs= 1.e-50;//ptsz->redshift_epsabs;//ptsz->patterson_epsabs;
   int show_neval = 0;//ptsz->patterson_show_neval;

    params = &V;

    // integral is symmetric (triangular configurations):
    // int(0,2*PI) = 2*int(0,PI)
    r = 2.*Integrate_using_Patterson_adaptive(0., _PI_,
                                            epsrel, epsabs,
                                            integrand_kSZ2_X,
                                            params,show_neval);



// //ROMBERG
// gsl_function F;
// double result_gsl, error;
// F.function = &integrand_kSZ2_X;
// F.params = params;
// int n_subintervals_gsl = 100;
// gsl_integration_romberg_workspace * w = gsl_integration_romberg_alloc (n_subintervals_gsl);
//
// size_t neval;
// double xin = 0.;
// double xout = _PI_;
// gsl_integration_romberg(&F,xin,xout,epsabs,epsrel,&result_gsl,&neval,w);
// gsl_integration_romberg_free(w);
// // *result = result_gsl;
// r = 2.*result_gsl;
// printf("r 0-pi = %.8e \n",r);

    cl_kSZ2_gal = r;

   free(b_l1_l2_l_1d);


   int index_l = index_ell_3;
   // double cl_kSZ2_gal = 0.;

  if (_kSZ_kSZ_gal_1h_)
  ptsz->cl_kSZ_kSZ_gal_1h[index_l] = cl_kSZ2_gal;
  else if(_kSZ_kSZ_lensmag_1halo_)
  ptsz->cl_kSZ_kSZ_lensmag_1h[index_l] = cl_kSZ2_gal;
  else if(_kSZ_kSZ_gal_2h_)
  ptsz->cl_kSZ_kSZ_gal_2h[index_l] = cl_kSZ2_gal;
  else if(_kSZ_kSZ_gal_3h_)
  ptsz->cl_kSZ_kSZ_gal_3h[index_l] = cl_kSZ2_gal;
  else if(_kSZ_kSZ_gal_hf_)
  ptsz->cl_kSZ_kSZ_gal_hf[index_l] = cl_kSZ2_gal;
  else if(_kSZ_kSZ_gallens_hf_)
  ptsz->cl_kSZ_kSZ_gallens_hf[index_l] = cl_kSZ2_gal;
  else if(_kSZ_kSZ_lens_hf_)
  ptsz->cl_kSZ_kSZ_lens_hf[index_l] = cl_kSZ2_gal;

   }

   else {
  // printf("integrating over redshift\n");
   class_call(integrate_over_redshift(pba,
                                      pnl,
                                      ppm,
                                      ppt,
                                      ptsz,
                                      Pvecback,
                                      Pvectsz),
                   ptsz->error_message,
                   ptsz->error_message);
 // Once the integration is done, the results of the integral
 // is stored in:
 // Pvectsz[ptsz->index_integral]
 // We collect the results hereafter...
          }

 }


   if (_hmf_){

      ptsz->hmf_int = Pvectsz[ptsz->index_integral];

   }

   if (_mean_y_){
      ptsz->y_monopole = Pvectsz[ptsz->index_integral]/pow(ptsz->Tcmb_gNU,1)/1.e6; //1e6 to convert Tcmb in micro Kelvins

   }
   if (_cib_monopole_){
    int index_cib1 = (int) Pvectsz[ptsz->index_frequency_for_cib_profile];
    ptsz->cib_monopole[index_cib1] = Pvectsz[ptsz->index_integral];

    }
   if (_cib_shotnoise_){
    int index_cib1 = (int) Pvectsz[ptsz->index_frequency_for_cib_profile];
    // int index_cib1 = (int) Pvectsz[ptsz->index_frequency_for_cib_profile];
    ptsz->cib_shotnoise[index_cib1] = Pvectsz[ptsz->index_integral];

    }
   if (_kSZ_kSZ_gal_1h_fft_){
   // if (1==0){
     int index_l = (int) Pvectsz[ptsz->index_multipole];
     ptsz->cl_kSZ_kSZ_gal_1h_fft[index_l] = Pvectsz[ptsz->index_integral]/(2.*_PI_)/(2.*_PI_);
   }

   if (_kSZ_kSZ_gal_2h_fft_){
   // if (1==0){
     int index_l = (int) Pvectsz[ptsz->index_multipole];
     ptsz->cl_kSZ_kSZ_gal_2h_fft[index_l] = Pvectsz[ptsz->index_integral]/(2.*_PI_)/(2.*_PI_);
   }

   if (_kSZ_kSZ_gal_3h_fft_){
   // if (1==0){
     int index_l = (int) Pvectsz[ptsz->index_multipole];
     ptsz->cl_kSZ_kSZ_gal_3h_fft[index_l] = Pvectsz[ptsz->index_integral]/(2.*_PI_)/(2.*_PI_);
   }

   if (_kSZ_kSZ_gallens_1h_fft_){
   // if (1==0){
     int index_l = (int) Pvectsz[ptsz->index_multipole];
     ptsz->cl_kSZ_kSZ_gallens_1h_fft[index_l] = Pvectsz[ptsz->index_integral]/(2.*_PI_)/(2.*_PI_);
   }

   if (_kSZ_kSZ_gallens_2h_fft_){
   // if (1==0){
     int index_l = (int) Pvectsz[ptsz->index_multipole];
     ptsz->cl_kSZ_kSZ_gallens_2h_fft[index_l] = Pvectsz[ptsz->index_integral]/(2.*_PI_)/(2.*_PI_);
   }

   if (_kSZ_kSZ_gallens_3h_fft_){
   // if (1==0){
     int index_l = (int) Pvectsz[ptsz->index_multipole];
     ptsz->cl_kSZ_kSZ_gallens_3h_fft[index_l] = Pvectsz[ptsz->index_integral]/(2.*_PI_)/(2.*_PI_);
   }

   if (_kSZ_kSZ_lens_1h_fft_){
   // if (1==0){
     int index_l = (int) Pvectsz[ptsz->index_multipole];
     ptsz->cl_kSZ_kSZ_lens_1h_fft[index_l] = Pvectsz[ptsz->index_integral]/(2.*_PI_)/(2.*_PI_);
   }

   if (_kSZ_kSZ_lens_2h_fft_){
   // if (1==0){
     int index_l = (int) Pvectsz[ptsz->index_multipole];
     ptsz->cl_kSZ_kSZ_lens_2h_fft[index_l] = Pvectsz[ptsz->index_integral]/(2.*_PI_)/(2.*_PI_);
   }

   if (_kSZ_kSZ_lens_3h_fft_){
   // if (1==0){
     int index_l = (int) Pvectsz[ptsz->index_multipole];
     ptsz->cl_kSZ_kSZ_lens_3h_fft[index_l] = Pvectsz[ptsz->index_integral]/(2.*_PI_)/(2.*_PI_);
   }

   if (_tSZ_power_spectrum_){
       int index_l = (int) Pvectsz[ptsz->index_multipole];
       ptsz->cl_sz_1h[index_l] = Pvectsz[ptsz->index_integral]
                                  *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                  /(2*_PI_*pow(ptsz->Tcmb_gNU,ptsz->exponent_unit));

   }

   if (_trispectrum_){
     int index_l = (int) Pvectsz[ptsz->index_multipole];
     int index_l_prime = (int) Pvectsz[ptsz->index_multipole_prime];

     ptsz->tllprime_sz[index_l][index_l_prime] = Pvectsz[ptsz->index_integral]
                                                 /pow(ptsz->Tcmb_gNU,2.*ptsz->exponent_unit);
  }
   //
   if (_2halo_){
    int index_l = (int) Pvectsz[ptsz->index_multipole];
    ptsz->cl_sz_2h[index_l] = Pvectsz[ptsz->index_integral]
                              *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                              /(2*_PI_*pow(ptsz->Tcmb_gNU,ptsz->exponent_unit));
    }

    if (_te_y_y_){
    int index_l = (int) Pvectsz[ptsz->index_multipole];
    ptsz->cl_te_y_y[index_l] = Pvectsz[ptsz->index_integral]
                                *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                /(2*_PI_*pow(ptsz->Tcmb_gNU,ptsz->exponent_unit));

    }
    if (_m_y_y_1h_){
    int index_l = (int) Pvectsz[ptsz->index_multipole];
    ptsz->m_y_y_1h[index_l] = Pvectsz[ptsz->index_integral]
                                *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                /(2*_PI_*pow(ptsz->Tcmb_gNU,ptsz->exponent_unit));

    }
    if (_m_y_y_2h_){
    int index_l = (int) Pvectsz[ptsz->index_multipole];
    ptsz->m_y_y_2h[index_l] = Pvectsz[ptsz->index_integral]
                                *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                /(2*_PI_*pow(ptsz->Tcmb_gNU,ptsz->exponent_unit));

    }
   if (_cov_Y_N_){
     int index_l = (int) Pvectsz[ptsz->index_multipole];
     int index_m = (int) Pvectsz[ptsz->index_mass_bin_1];


    ptsz->cov_Y_N[index_l][index_m] = Pvectsz[ptsz->index_integral]
                                      /pow(ptsz->Tcmb_gNU,ptsz->exponent_unit);


  }


   if (_cov_N_N_){
     int index_m_1 = (int) Pvectsz[ptsz->index_mass_bin_1];
    ptsz->cov_N_N[index_m_1] = Pvectsz[ptsz->index_integral];

  }

  if (_cov_N_N_hsv_){
    int index_m_1 = (int) Pvectsz[ptsz->index_mass_bin_1];
    int index_m_2 = (int) Pvectsz[ptsz->index_mass_bin_2];

   ptsz->cov_N_N_hsv[index_m_1][index_m_2] = Pvectsz[ptsz->index_integral];
   ptsz->cov_N_N_hsv[index_m_2][index_m_1] = ptsz->cov_N_N_hsv[index_m_1][index_m_2];

 }
  if (_cov_Y_N_next_order_){
    int index_l = (int) Pvectsz[ptsz->index_multipole];
    int index_m = (int) Pvectsz[ptsz->index_mass_bin_1];


   ptsz->cov_Y_N_next_order[index_l][index_m] = Pvectsz[ptsz->index_integral]
                                                /pow(ptsz->Tcmb_gNU,ptsz->exponent_unit);


 }
 if (_cov_Y_Y_ssc_){
   int index_multipole_1 = (int) Pvectsz[ptsz->index_multipole_1];
   int index_multipole_2 = (int) Pvectsz[ptsz->index_multipole_2];


  ptsz->cov_Y_Y_ssc[index_multipole_1][index_multipole_2] = Pvectsz[ptsz->index_integral]
                                                            /pow(ptsz->Tcmb_gNU,2.*ptsz->exponent_unit);
  ptsz->cov_Y_Y_ssc[index_multipole_2][index_multipole_1] = ptsz->cov_Y_Y_ssc[index_multipole_1][index_multipole_2];


}


  if (_kSZ_kSZ_1h_){
    int index_l = (int) Pvectsz[ptsz->index_multipole];
    ptsz->cl_kSZ_kSZ_1h[index_l] = Pvectsz[ptsz->index_integral] // dimensionless
                                  *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                  /(2*_PI_);
 }

  if (_kSZ_kSZ_2h_){
    int index_l = (int) Pvectsz[ptsz->index_multipole];
    ptsz->cl_kSZ_kSZ_2h[index_l] = Pvectsz[ptsz->index_integral] // dimensionless
                                  *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                  /(2*_PI_);
 }


 if (_tSZ_tSZ_tSZ_1halo_){
   int index_l = (int) Pvectsz[ptsz->index_multipole];
  ptsz->b_tSZ_tSZ_tSZ_1halo[index_l] = Pvectsz[ptsz->index_integral]/pow(ptsz->Tcmb_gNU,3)/1.e18; // dimensionless
}

 if (_tSZ_tSZ_tSZ_2h_){
   int index_l = (int) Pvectsz[ptsz->index_multipole];
  ptsz->b_tSZ_tSZ_tSZ_2h[index_l] = Pvectsz[ptsz->index_integral]/pow(ptsz->Tcmb_gNU,3)/1.e18; // dimensionless
}

 if (_tSZ_tSZ_tSZ_3h_){
   int index_l = (int) Pvectsz[ptsz->index_multipole];
  ptsz->b_tSZ_tSZ_tSZ_3h[index_l] = Pvectsz[ptsz->index_integral]/pow(ptsz->Tcmb_gNU,3)/1.e18; // dimensionless
}

  if (_kSZ_kSZ_tSZ_1h_){
    int index_l = (int) Pvectsz[ptsz->index_multipole];
    ptsz->b_kSZ_kSZ_tSZ_1h[index_l] = Pvectsz[ptsz->index_integral]/pow(ptsz->Tcmb_gNU,0.)/1.e0; // dimensionless
 }

  if (_kSZ_kSZ_tSZ_2h_){
    int index_l = (int) Pvectsz[ptsz->index_multipole];
    ptsz->b_kSZ_kSZ_tSZ_2h[index_l] = Pvectsz[ptsz->index_integral]/pow(ptsz->Tcmb_gNU,0.)/1.e0; // dimensionless
 }

  if (_kSZ_kSZ_tSZ_3h_){
    int index_l = (int) Pvectsz[ptsz->index_multipole];
    ptsz->b_kSZ_kSZ_tSZ_3h[index_l] = Pvectsz[ptsz->index_integral]/pow(ptsz->Tcmb_gNU,0.)/1.e0; // dimensionless
 }

 if (_pk_at_z_1h_){
  int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
  ptsz->pk_at_z_1h[index_k] = Pvectsz[ptsz->index_integral];

}

if (_pk_at_z_2h_){
 int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
 ptsz->pk_at_z_2h[index_k] = Pvectsz[ptsz->index_integral];

}
 if (_pk_gg_at_z_1h_){
  int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
  ptsz->pk_gg_at_z_1h[index_k] = Pvectsz[ptsz->index_integral];

}

if (_pk_gg_at_z_2h_){
 int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
 ptsz->pk_gg_at_z_2h[index_k] = Pvectsz[ptsz->index_integral];

}
 if (_pk_bb_at_z_1h_){
  int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
  ptsz->pk_bb_at_z_1h[index_k] = Pvectsz[ptsz->index_integral];

}

if (_pk_bb_at_z_2h_){
 int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
 ptsz->pk_bb_at_z_2h[index_k] = Pvectsz[ptsz->index_integral];

}
if (_pk_b_at_z_2h_){
 int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
 ptsz->pk_b_at_z_2h[index_k] = Pvectsz[ptsz->index_integral];

}
 if (_pk_em_at_z_1h_){
  int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
  ptsz->pk_em_at_z_1h[index_k] = Pvectsz[ptsz->index_integral];

}

if (_pk_em_at_z_2h_){
 int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
 ptsz->pk_em_at_z_2h[index_k] = Pvectsz[ptsz->index_integral];

}

 if (_pk_HI_at_z_1h_){
  int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
  ptsz->pk_HI_at_z_1h[index_k] = Pvectsz[ptsz->index_integral];

}

if (_pk_HI_at_z_2h_){
 int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
 ptsz->pk_HI_at_z_2h[index_k] = Pvectsz[ptsz->index_integral];

}

 if (_bk_at_z_1h_){
  int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
  ptsz->bk_at_z_1h[index_k] = Pvectsz[ptsz->index_integral];
}

if (_bk_at_z_2h_){
 int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
 ptsz->bk_at_z_2h[index_k] = Pvectsz[ptsz->index_integral];
}

if (_bk_at_z_3h_){
 int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
 ptsz->bk_at_z_3h[index_k] = Pvectsz[ptsz->index_integral];
}

 if (_bk_ttg_at_z_1h_){
  int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
  ptsz->bk_ttg_at_z_1h[index_k] = Pvectsz[ptsz->index_integral];
}

if (_bk_ttg_at_z_2h_){
 int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
 ptsz->bk_ttg_at_z_2h[index_k] = Pvectsz[ptsz->index_integral];
}

if (_bk_ttg_at_z_3h_){
 int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
 ptsz->bk_ttg_at_z_3h[index_k] = Pvectsz[ptsz->index_integral];
}



 // Collect gxg 1-halo at each multipole:
 // result in y-units (dimensionless)
 // [l(l+1)/2pi]*cl
 if (_gal_gal_1h_){
  int index_l = (int) Pvectsz[ptsz->index_multipole];
  ptsz->cl_gal_gal_1h[index_l] = Pvectsz[ptsz->index_integral]
                                  *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                  /(2*_PI_);

}


 // Collect gxg 2-halo at each multipole:
 // result in y-units (dimensionless)
 // [l(l+1)/2pi]*cl
 if (_gal_gal_2h_){
  int index_l = (int) Pvectsz[ptsz->index_multipole];
  ptsz->cl_gal_gal_2h[index_l] = Pvectsz[ptsz->index_integral]
                                  *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                  /(2*_PI_);

}


 // Collect gxg effective approach (halofit/hmcode) at each multipole:
 // result in y-units (dimensionless)
 // [l(l+1)/2pi]*cl
 if (_gal_gal_hf_){

  int index_l = (int) Pvectsz[ptsz->index_multipole];
  ptsz->cl_gal_gal_hf[index_l] = Pvectsz[ptsz->index_integral]
                                  *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                  /(2*_PI_);

}

// Collect gxlens 1-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_gal_lens_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_gal_lens_1h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}

// Collect gxlens 2-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_gal_lens_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_gal_lens_2h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}

// Collect gxlens (effective approach) at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_gal_lens_hf_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_gal_lens_hf[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}
// Collect gxlensmag 1-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_gal_lensmag_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_gal_lensmag_1h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}
// Collect gxlensmag 1-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_gal_lensmag_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_gal_lensmag_2h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}
// Collect gxlensmag 1-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_gal_gallens_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_gal_gallens_1h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);
  // printf("%.3e\n",ptsz->cl_gal_gallens_1h[index_l]);
  if (isnan(ptsz->cl_gal_gallens_1h[index_l])||isinf(ptsz->cl_gal_gallens_1h[index_l])){
  printf("nan or inf in cls 1h\n");
  exit(0);
  }
}
// Collect gxlensmag 1-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_gal_gallens_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_gal_gallens_2h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}
// Collect gxlensmag 1-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_gallens_gallens_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_gallens_gallens_1h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);
  // printf("%.3e\n",ptsz->cl_gal_gallens_1h[index_l]);
  if (isnan(ptsz->cl_gallens_gallens_1h[index_l])||isinf(ptsz->cl_gallens_gallens_1h[index_l])){
  printf("nan or inf in cls 1h\n");
  exit(0);
  }
}
// Collect gxlensmag 1-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_gallens_gallens_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_gallens_gallens_2h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}
// Collect gxlensmag 1-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_gallens_lens_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_gallens_lens_1h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);
  // printf("%.3e\n",ptsz->cl_gal_gallens_1h[index_l]);
  if (isnan(ptsz->cl_gallens_lens_1h[index_l])||isinf(ptsz->cl_gallens_lens_1h[index_l])){
  printf("nan or inf in cls 1h\n");
  exit(0);
  }
}
// Collect gxlensmag 1-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_gallens_lens_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_gallens_lens_2h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}

// Collect gxlensmag effective approach at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_gal_lensmag_hf_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_gal_lensmag_hf[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}

// Collect gxlensmag 1-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_lensmag_lensmag_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_lensmag_lensmag_1h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}

// Collect gxlensmag 2-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_lensmag_lensmag_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_lensmag_lensmag_2h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}


// Collect gxlensmag 2-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_lensmag_lensmag_hf_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_lensmag_lensmag_hf[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}

// Collect gxlensmag 1-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_lens_lensmag_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_lens_lensmag_1h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}

// Collect gxlensmag 2-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_lens_lensmag_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_lens_lensmag_2h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);
 //printf("cl_lens_lensmag = %.3e\n",ptsz->cl_lens_lensmag_2h[index_l]);

}

// Collect gxlensmag 2-halo at each multipole:
// result in y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_lens_lensmag_hf_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_lens_lensmag_hf[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);
 //printf("cl_lens_lensmag = %.3e\n",ptsz->cl_lens_lensmag_2h[index_l]);

}


 // Collect Yxgal 1-halo at each multipole:
 // result in 10^-6 y-units (dimensionless)
 // [l(l+1)/2pi]*cl
 if (_tSZ_gal_1h_){
  int index_l = (int) Pvectsz[ptsz->index_multipole];
  ptsz->cl_tSZ_gal_1h[index_l] = Pvectsz[ptsz->index_integral]
                                  *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                  /(2*_PI_)
                                  /pow(ptsz->Tcmb_gNU,1);
  //printf("ell = %.3e and cl_yg = %.3e\n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_1h[index_l]);

}


if (_tSZ_gal_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_tSZ_gal_2h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_)
                                 /pow(ptsz->Tcmb_gNU,1);
 //printf("ell = %.3e and cl_yg = %.3e\n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_1h[index_l]);

}

 // Collect Yxmu 1-halo at each multipole:
 // result in 10^-6 y-units (dimensionless)
 // [l(l+1)/2pi]*cl
 if (_tSZ_lensmag_1h_){
  int index_l = (int) Pvectsz[ptsz->index_multipole];
  ptsz->cl_tSZ_lensmag_1h[index_l] = Pvectsz[ptsz->index_integral]
                                  *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                  /(2*_PI_)
                                  /pow(ptsz->Tcmb_gNU,1);
  //printf("ell = %.3e and cl_yg = %.3e\n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_1h[index_l]);

}


if (_tSZ_lensmag_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_tSZ_lensmag_2h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_)
                                 /pow(ptsz->Tcmb_gNU,1);
 //printf("ell = %.3e and cl_yg = %.3e\n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_1h[index_l]);

}



// Collect Yxcib 1-halo at each multipole:
// result in 10^-6 y-units (dimensionless)
// [l(l+1)/2pi]*cl
if (_tSZ_cib_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 int index_cib1 = (int) Pvectsz[ptsz->index_frequency_for_cib_profile];
 ptsz->cl_tSZ_cib_1h[index_cib1][index_l] = Pvectsz[ptsz->index_integral]
                                             *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                             /(2*_PI_)
                                             /pow(ptsz->Tcmb_gNU,1);
 //printf("ell = %.3e and cl_yg = %.3e\n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_1h[index_l]);

}


if (_tSZ_cib_2h_){
int index_l = (int) Pvectsz[ptsz->index_multipole];
int index_cib1 = (int) Pvectsz[ptsz->index_frequency_for_cib_profile];
ptsz->cl_tSZ_cib_2h[index_cib1][index_l] = Pvectsz[ptsz->index_integral]
                                            *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                            /(2*_PI_)
                                            /pow(ptsz->Tcmb_gNU,1);
//printf("ell = %.3e and cl_yg = %.3e\n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_1h[index_l]);

}

if (_lens_cib_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 int index_cib1 = (int) Pvectsz[ptsz->index_frequency_for_cib_profile];
 ptsz->cl_lens_cib_1h[index_cib1][index_l] = Pvectsz[ptsz->index_integral]
                                             *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                             /(2*_PI_);
 //printf("ell = %.3e and cl_yg = %.3e\n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_1h[index_l]);

}


if (_lens_cib_2h_){
int index_l = (int) Pvectsz[ptsz->index_multipole];
int index_cib1 = (int) Pvectsz[ptsz->index_frequency_for_cib_profile];
ptsz->cl_lens_cib_2h[index_cib1][index_l] = Pvectsz[ptsz->index_integral]
                                            *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                            /(2*_PI_);
//printf("ell = %.3e and cl_yg = %.3e\n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_1h[index_l]);

}

if (_gal_cib_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 int index_cib1 = (int) Pvectsz[ptsz->index_frequency_for_cib_profile];
 ptsz->cl_gal_cib_1h[index_cib1][index_l] = Pvectsz[ptsz->index_integral]
                                             *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                             /(2*_PI_);
 //printf("ell = %.3e and cl_yg = %.3e\n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_1h[index_l]);

}


if (_gal_cib_2h_){
int index_l = (int) Pvectsz[ptsz->index_multipole];
int index_cib1 = (int) Pvectsz[ptsz->index_frequency_for_cib_profile];
ptsz->cl_gal_cib_2h[index_cib1][index_l] = Pvectsz[ptsz->index_integral]
                                           *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                           /(2*_PI_);
//printf("ell = %.3e and cl_yg = %.3e\n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_1h[index_l]);

}


// Collect ngalxngal 1-halo at each multipole:
// [l(l+1)/2pi]*cl
if (_ngal_ngal_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 int index_ngal1 = (int) Pvectsz[ptsz->index_ngal_for_galaxy_profile];
 int index_ngal2 = (int) Pvectsz[ptsz->index_ngal_prime_for_galaxy_profile];

 ptsz->cl_ngal_ngal_1h[index_ngal1][index_ngal2][index_l] = Pvectsz[ptsz->index_integral]
                                                            *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                                            /(2*_PI_);

 ptsz->cl_ngal_ngal_1h[index_ngal2][index_ngal1][index_l] = ptsz->cl_ngal_ngal_1h[index_ngal1][index_ngal2][index_l];

}

// Collect ngalxngal 1-halo at each multipole:
// [l(l+1)/2pi]*cl
if (_ngal_ngal_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 int index_ngal1 = (int) Pvectsz[ptsz->index_ngal_for_galaxy_profile];
 int index_ngal2 = (int) Pvectsz[ptsz->index_ngal_prime_for_galaxy_profile];

 ptsz->cl_ngal_ngal_2h[index_ngal1][index_ngal2][index_l] = Pvectsz[ptsz->index_integral]
                                                            *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                                            /(2*_PI_);

 ptsz->cl_ngal_ngal_2h[index_ngal2][index_ngal1][index_l] = ptsz->cl_ngal_ngal_2h[index_ngal1][index_ngal2][index_l];


}

// Collect ngalxngal hf at each multipole:
// [l(l+1)/2pi]*cl
if (_ngal_ngal_hf_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 int index_ngal1 = (int) Pvectsz[ptsz->index_ngal_for_galaxy_profile];
 int index_ngal2 = (int) Pvectsz[ptsz->index_ngal_prime_for_galaxy_profile];

 ptsz->cl_ngal_ngal_hf[index_ngal1][index_ngal2][index_l] = Pvectsz[ptsz->index_integral]
                                                            *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                                            /(2*_PI_);

 ptsz->cl_ngal_ngal_hf[index_ngal2][index_ngal1][index_l] = ptsz->cl_ngal_ngal_hf[index_ngal1][index_ngal2][index_l];

}



// Collect ngalxngal 1-halo at each multipole:
// [l(l+1)/2pi]*cl
if (_ngal_lens_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 int index_ngal1 = (int) Pvectsz[ptsz->index_ngal_for_galaxy_profile];


 ptsz->cl_ngal_lens_1h[index_ngal1][index_l] = Pvectsz[ptsz->index_integral]
                                                            *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                                            /(2*_PI_);

}

// Collect ngalxngal 1-halo at each multipole:
// [l(l+1)/2pi]*cl
if (_ngal_lens_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 int index_ngal1 = (int) Pvectsz[ptsz->index_ngal_for_galaxy_profile];


 ptsz->cl_ngal_lens_2h[index_ngal1][index_l] = Pvectsz[ptsz->index_integral]
                                                            *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                                            /(2*_PI_);


}

// Collect ngalxngal hf at each multipole:
// [l(l+1)/2pi]*cl
if (_ngal_lens_hf_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 int index_ngal1 = (int) Pvectsz[ptsz->index_ngal_for_galaxy_profile];


 ptsz->cl_ngal_lens_hf[index_ngal1][index_l] = Pvectsz[ptsz->index_integral]
                                               *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                               /(2*_PI_);

}



// Collect cibxcib 1-halo at each multipole:
// [l(l+1)/2pi]*cl
if (_cib_cib_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 int index_cib1 = (int) Pvectsz[ptsz->index_frequency_for_cib_profile];
 int index_cib2 = (int) Pvectsz[ptsz->index_frequency_prime_for_cib_profile];

 ptsz->cl_cib_cib_1h[index_cib1][index_cib2][index_l] = Pvectsz[ptsz->index_integral]
                                                        *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                                        /(2*_PI_);

ptsz->cl_cib_cib_1h[index_cib2][index_cib1][index_l] = ptsz->cl_cib_cib_1h[index_cib1][index_cib2][index_l];

 //printf("ell = %.3e and cl_yg = %.3e\n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_1h[index_l]);

}


if (_cib_cib_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 int index_cib1 = (int) Pvectsz[ptsz->index_frequency_for_cib_profile];
 int index_cib2 = (int) Pvectsz[ptsz->index_frequency_prime_for_cib_profile];

ptsz->cl_cib_cib_2h[index_cib1][index_cib2][index_l] = Pvectsz[ptsz->index_integral]
                                                       *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                                       /(2*_PI_);

ptsz->cl_cib_cib_2h[index_cib2][index_cib1][index_l] = ptsz->cl_cib_cib_2h[index_cib1][index_cib2][index_l];
//printf("ell = %.3e and cl_cib_cib = %.3e\n",ptsz->ell[index_l],ptsz->cl_cib_cib_2h[index_l]);

}

// Collect kappaxkappa 1-halo at each multipole:
// [l(l+1)/2pi]*cl
if (_lens_lens_1h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_lens_lens_1h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);
}

// Collect kappaxkappa 2-halo at each multipole:
// [l(l+1)/2pi]*cl
if (_lens_lens_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_lens_lens_2h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);
}


if (_lens_lens_hf_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_lens_lens_hf[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_);

}
 // Collect YxPhi 1-halo at each multipole:
 // result in 10^-6 y-units (dimensionless)
 // [l^2(l+1)/2pi]*cl
 if (_tSZ_lens_1h_){
  int index_l = (int) Pvectsz[ptsz->index_multipole];
  ptsz->cl_tSZ_lens_1h[index_l] = Pvectsz[ptsz->index_integral]
                                  *ptsz->ell[index_l]*ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                  /(2*_PI_)
                                  /pow(ptsz->Tcmb_gNU,1);
}

// Collect YxPhi 2-halo at each multipole:
// result in 10^-6 y-units (dimensionless)
// [l^2(l+1)/2pi]*cl
if (_tSZ_lens_2h_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_tSZ_lens_2h[index_l] = Pvectsz[ptsz->index_integral]
                                 *ptsz->ell[index_l]*ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                 /(2*_PI_)
                                 /pow(ptsz->Tcmb_gNU,1);
}

// Collect ISWxPhi at each multipole:
// [l(l+1)/2pi]*cl
// result in y-units (dimensionless)
if (_isw_lens_){
int index_l = (int) Pvectsz[ptsz->index_multipole];
ptsz->cl_isw_lens[index_l] = Pvectsz[ptsz->index_integral]
                             *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                             /(2*_PI_);
}

// Collect YxISW at each multipole:
// [l(l+1)/2pi]*cl
// result in y-units (dimensionless)
if (_isw_tsz_){
 int index_l = (int) Pvectsz[ptsz->index_multipole];
 ptsz->cl_isw_tsz[index_l] = Pvectsz[ptsz->index_integral]
                             *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                             /(2*_PI_*pow(ptsz->Tcmb_gNU*1e6,1.));
}

if (_isw_auto_){
  int index_l = (int) Pvectsz[ptsz->index_multipole];
  //int index_m = (int) Pvectsz[ptsz->index_mass_bin_1];
 ptsz->cl_isw_auto[index_l] =  Pvectsz[ptsz->index_integral]
                               *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                               /(2*_PI_); //dimensionnless
}


if (_szrates_){
  int index_rate = (int) Pvectsz[ptsz->index_szrate];
  ptsz->szrate[index_rate] = Pvectsz[ptsz->index_integral]*4.*_PI_*ptsz->fsky_from_skyfracs;;
}

return _SUCCESS_;
}





double integrand_at_m_and_z(double logM,
                             double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct primordial * ppm,
                             struct nonlinear * pnl,
                             struct perturbs * ppt,
                             struct tszspectrum * ptsz)
    {

   int index_md = (int) pvectsz[ptsz->index_md];

   int index_l = (int) pvectsz[ptsz->index_multipole];





   double z = pvectsz[ptsz->index_z];




   double chi = sqrt(pvectsz[ptsz->index_chi2]);

   //collect the necessary overdensity masses
   //for the different fields
   do_mass_conversions(logM,z,pvecback,pvectsz,pba,ptsz);


  // printf("mhmf = %.3e m500c = %.3e m = %.3e z=%.3e\n",pvectsz[ptsz->index_mass_for_hmf],pvectsz[ptsz->index_m500c],exp(logM),z);
   //Return the HMF - dn/dlogM in units of h^3 Mpc^-3
   //result stored in pvectsz[ptsz->index_hmf]

   evaluate_HMF_at_logM_and_z(log(pvectsz[ptsz->index_mass_for_hmf]),z,pvecback,pvectsz,pba,pnl,ptsz);

   pvectsz[ptsz->index_dlnMdeltadlnM]= 1.;//evaluate_dlnMdeltadlnM(log(pvectsz[ptsz->index_mVIR]),
   //                                                            pvecback,
   //                                                            pvectsz,
   //                                                            pba,
   //                                                            pnl,
   //                                                            ptsz);
   //                                                        //  }

   double r_delta_gal, c_delta_gal, m_delta_gal;
   double r_delta_electron_pressure, c_delta_electron_pressure, m_delta_electron_pressure; //(tsz)
   double r_delta_matter, c_delta_matter, m_delta_matter;
   double r_delta_lensing, c_delta_lensing, m_delta_lensing;
   double r_delta_cib, c_delta_cib, m_delta_cib;
   double r_delta_electron_density, c_delta_electron_density, m_delta_electron_density; //(ksz)
   double r_delta_HI_density, c_delta_HI_density, m_delta_HI_density; //(HI)

   m_delta_gal = pvectsz[ptsz->index_mass_for_galaxies];
   r_delta_gal = pvectsz[ptsz->index_radius_for_galaxies];
   c_delta_gal = pvectsz[ptsz->index_concentration_for_galaxies];

   m_delta_electron_pressure = pvectsz[ptsz->index_mass_for_electron_pressure];
   r_delta_electron_pressure = pvectsz[ptsz->index_radius_for_electron_pressure];
   c_delta_electron_pressure = pvectsz[ptsz->index_concentration_for_electron_pressure];

   m_delta_electron_density = pvectsz[ptsz->index_mass_for_electron_density];
   r_delta_electron_density = pvectsz[ptsz->index_radius_for_electron_density];
   c_delta_electron_density = pvectsz[ptsz->index_concentration_for_electron_density];

   m_delta_HI_density = pvectsz[ptsz->index_mass_for_HI_density];
   r_delta_HI_density = pvectsz[ptsz->index_radius_for_HI_density];
   c_delta_HI_density = pvectsz[ptsz->index_concentration_for_HI_density];

   m_delta_matter = pvectsz[ptsz->index_mass_for_matter_density];
   r_delta_matter = pvectsz[ptsz->index_radius_for_matter_density];
   c_delta_matter = pvectsz[ptsz->index_concentration_for_matter_density];

   m_delta_lensing = pvectsz[ptsz->index_mass_for_matter_density];
   r_delta_lensing = pvectsz[ptsz->index_radius_for_matter_density];
   c_delta_lensing = pvectsz[ptsz->index_concentration_for_matter_density];

   m_delta_cib = pvectsz[ptsz->index_mass_for_cib];
   r_delta_cib = pvectsz[ptsz->index_radius_for_cib];
   c_delta_cib = pvectsz[ptsz->index_concentration_for_cib];

   // evaluate_completeness(pvecback,pvectsz,pba,ptsz);
   pvectsz[ptsz->index_completeness] = 1.;



   // double z = pvectsz[ptsz->index_z];
   double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi
   double kl;

   double kl1;
   double kl2;
   double kl3;
   double l1 = ptsz->ell[index_l];
   double l2 = ptsz->bispectrum_lambda_k2*ptsz->ell[index_l];
   double l3 = ptsz->bispectrum_lambda_k3*ptsz->ell[index_l];

   if (_pk_at_z_1h_
    || _pk_gg_at_z_1h_
    || _pk_at_z_2h_
    || _pk_gg_at_z_2h_
    || _pk_bb_at_z_1h_
    || _pk_bb_at_z_2h_
    || _pk_b_at_z_2h_
    || _pk_em_at_z_1h_
    || _pk_em_at_z_2h_
    || _pk_HI_at_z_1h_
    || _pk_HI_at_z_2h_
    || _bk_at_z_1h_
    || _bk_at_z_2h_
    || _bk_at_z_3h_
    || _bk_ttg_at_z_1h_
    || _bk_ttg_at_z_2h_
    || _bk_ttg_at_z_3h_
  ){
     int index_k = (int) pvectsz[ptsz->index_k_for_pk_hm];
     kl = ptsz->k_for_pk_hm[index_k];
      }
   else{
     kl = (ptsz->ell[index_l]+0.5)/d_A;
     kl1 = (l1+0.5)/d_A;
     kl2 = (l2+0.5)/d_A;
     kl3 = (l3+0.5)/d_A;
    }


double damping_1h_term;
   if (ptsz->damping_1h_term ==  1){
   damping_1h_term = (1.-exp(-pow(kl*pba->h/ptsz->kstar_damping_1h_term_Mpc,2.)));
 }
   else{
   damping_1h_term = 1.;
   }


   if (_szrates_){
   int index_szrate = (int) pvectsz[ptsz->index_szrate];
     if (ptsz->sz_verbose>10)
   printf("computing integrand at m = %.3e and z = %.3e for rate id = %d.\n",exp(logM),z,index_szrate);
   double thetap = get_theta_at_m_and_z(pvectsz[ptsz->index_m500c],z,ptsz,pba);
   double yp = get_y_at_m_and_z(pvectsz[ptsz->index_m500c],z,ptsz,pba);
   double noisep;
   if ((thetap < ptsz->thetas[0]) || (thetap > ptsz->thetas[ptsz->nthetas-1])){
     noisep = 1e100;
   }
   else {
   noisep = pwl_value_1d(ptsz->nthetas,
                                ptsz->thetas,
                                ptsz->sky_averaged_ylims,
                                thetap);
   }
   double snrp = yp/noisep;

   if (ptsz->sz_verbose>10) printf("predicted snr = %.3e cat snr = %.3e.\n",snrp,ptsz->szcat_snr[index_szrate]);
   double arg = (snrp - ptsz->szcat_snr[index_szrate])/sqrt(2.);
   double prate = 1./sqrt(2.*_PI_)*exp(-arg*arg);
   pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]*prate;

   arg = (snrp - ptsz->sn_cutoff)/sqrt(2.);
   prate = (erf(arg) + 1.)/2.;
   pvectsz[ptsz->index_integrand] *= prate;


 }


   if (_hmf_){


// printf("z=%.3e m=%.3e\n",z,pvectsz[ptsz->index_m500c]);

   double thetap = get_theta_at_m_and_z(pvectsz[ptsz->index_m500c],z,ptsz,pba);
   double ypp = get_y_at_m_and_z(pvectsz[ptsz->index_m500c],z,ptsz,pba);
   double noisep;
   double prate;
   if ((thetap < ptsz->thetas[0]) || (thetap > ptsz->thetas[ptsz->nthetas-1])){
     prate = 0.;
   }
   else {
   noisep = pwl_value_1d(ptsz->nthetas,
                                ptsz->thetas,
                                ptsz->sky_averaged_ylims,
                                thetap);

  double snrp = ypp/noisep;

  double arg = (snrp - ptsz->sn_cutoff)/sqrt(2.);
    prate = (erf(arg) + 1.)/2.;


      // printf("z=%.3e m=%.3e noisep = %.3e snrp = %.3e prate = %.3e\n",z,pvectsz[ptsz->index_m500c],noisep,snrp,prate);
   }


   pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                     *prate;

  return  pvectsz[ptsz->index_integrand];
   }




   else if (_mean_y_){
     // pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];
     evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
     double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];
     pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                       pvectsz[ptsz->index_hmf]
                                       //*pvectsz[ptsz->index_completeness]
                                       *pow(pressure_profile_at_ell,1.);

   }
   else if (_cib_monopole_){

     evaluate_cib_profile(m_delta_cib,r_delta_cib,c_delta_cib,pvecback,pvectsz,pba,ptsz);
     double cib_profile = pvectsz[ptsz->index_cib_profile];
     pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                       *pow(cib_profile,1.);

   }

   else if (_cib_shotnoise_){

     evaluate_cib_profile(m_delta_cib,r_delta_cib,c_delta_cib,pvecback,pvectsz,pba,ptsz);
     double cib_profile = pvectsz[ptsz->index_cib_profile];
     pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                       *pow(cib_profile,1.);

   }


   else if (_tSZ_power_spectrum_){

      // pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];

      evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
      double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_completeness]
                                           *pow(pressure_profile_at_ell,2.)
                                           *damping_1h_term;

      if (isinf(pvectsz[ptsz->index_integrand])){
      printf("%.8e %.8e\n",pvectsz[ptsz->index_integrand],
      pressure_profile_at_ell);
      exit(0);
    }

    }


   else if (_trispectrum_){
     // pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];
     evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
     double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];


     int index_l_prime = (int) pvectsz[ptsz->index_multipole_prime];
     // pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l_prime];
     evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
     double pressure_profile_at_ell_prime = pvectsz[ptsz->index_pressure_profile];

     pvectsz[ptsz->index_integrand] = pvectsz[ptsz->index_hmf]
                                     *pvectsz[ptsz->index_completeness]
                                     *pow(pressure_profile_at_ell,2.)
                                     *pow(pressure_profile_at_ell_prime,2.);
 }

 else if (_2halo_){
     // pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];
     evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
     double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];
     evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
     //this integrand is squared afterward
     pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                        *pvectsz[ptsz->index_halo_bias]
                                        *pvectsz[ptsz->index_completeness]
                                        *pow(pressure_profile_at_ell,1.);
 }

 else if  (_te_y_y_){
        // pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];
        evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
        double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];
        //int index_l = (int) pvectsz[ptsz->index_multipole];
        evaluate_temperature_mass_relation(pvecback,pvectsz,pba,ptsz);

        pvectsz[ptsz->index_integrand] = pvectsz[ptsz->index_te_of_m]
                                         //*pvectsz[ptsz->index_chi2]
                                         *pvectsz[ptsz->index_hmf]
                                         // *pvectsz[ptsz->index_completeness]
                                         *pow(pvectsz[ptsz->index_pressure_profile],2.);


        }
 else if  (_m_y_y_1h_){

          //int index_l = (int) pvectsz[ptsz->index_multipole];
           //evaluate_temperature_mass_relation(pvecback,pvectsz,pba,ptsz);
           // pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];
           evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
           double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];
           pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                             pvectsz[ptsz->index_mass_for_hmf]
                                             *pvectsz[ptsz->index_hmf]
                                             //*pvectsz[ptsz->index_dlnMdeltadlnM]
                                             *pvectsz[ptsz->index_completeness]
                                             *pow(pressure_profile_at_ell,2.)
                                             *damping_1h_term;
        }
 else if  (_m_y_y_2h_){
         // pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];
         evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
         double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];


           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

           //this integrand is squared afterward
           pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                              *pvectsz[ptsz->index_mass_for_hmf]
                                              //*pvectsz[ptsz->index_dlnMdeltadlnM]
                                              *pvectsz[ptsz->index_halo_bias]
                                              *pvectsz[ptsz->index_completeness]
                                              *pow(pressure_profile_at_ell,1.);
                                              //*pow(pvectsz[ptsz->index_pk_for_halo_bias],0.5);

        }
 else if  (_cov_Y_N_){
           // pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];
           evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
           double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];

           pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                             pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness]
                                             *pow(pvectsz[ptsz->index_pressure_profile],2.);


        }

 else if  (_cov_N_N_){

           pvectsz[ptsz->index_integrand] =  ptsz->Omega_survey
                                             //*pvectsz[ptsz->index_chi2]
                                             *pvectsz[ptsz->index_hmf];
                                             // *pvectsz[ptsz->index_completeness];
                      }

 else if  (_cov_N_N_hsv_){


           if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

           // evaluate_sigma2_hsv(pvecback,pvectsz,pba,pnl,ptsz);
           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
           pvectsz[ptsz->index_halo_bias] = 1.;

           pvectsz[ptsz->index_integrand] =  ptsz->Omega_survey
                                             //*pvectsz[ptsz->index_chi2]
                                             *pvectsz[ptsz->index_hmf]
                                             // *pvectsz[ptsz->index_completeness]
                                             *pvectsz[ptsz->index_halo_bias];
                                           }

           if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
           pvectsz[ptsz->index_halo_bias] = 1.;

           pvectsz[ptsz->index_integrand] =  ptsz->Omega_survey
                                             //*pvectsz[ptsz->index_chi2]
                                             *pvectsz[ptsz->index_hmf]
                                             // *pvectsz[ptsz->index_completeness]
                                             *pvectsz[ptsz->index_halo_bias];

                                           }

        }


 else if  (_cov_Y_N_next_order_){

           if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

           evaluate_sigma2_hsv(pvecback,pvectsz,pba,pnl,ptsz);
           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

           pvectsz[ptsz->index_integrand] =  ptsz->Omega_survey
                                             //*pvectsz[ptsz->index_chi2]
                                             *pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness]
                                             *pvectsz[ptsz->index_halo_bias]
                                             *pvectsz[ptsz->index_sigma2_hsv];
                                           }

           if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {
         // pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];
         evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
         double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];

           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

           pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                             pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness]
                                             *pvectsz[ptsz->index_halo_bias]
                                             *pow(pressure_profile_at_ell,2.);

                                           }
        }
 else if  (_cov_Y_Y_ssc_){
   // pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];
   evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];

           if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

           evaluate_sigma2_hsv(pvecback,pvectsz,pba,pnl,ptsz);
           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

           pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                             pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness]
                                             *pvectsz[ptsz->index_halo_bias]
                                             *pow(pressure_profile_at_ell,2.)
                                             *pvectsz[ptsz->index_sigma2_hsv];
                                           }

           if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

           pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                             pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness]
                                             *pvectsz[ptsz->index_halo_bias]
                                             *pow(pressure_profile_at_ell,2.);

                                           }
        }

 else if (_kSZ_kSZ_1h_){
   // int index_l = (int) pvectsz[ptsz->index_multipole];
   // double l1 = ptsz->ell[index_l];

   evaluate_tau_profile(kl,pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];

   pvectsz[ptsz->index_integrand] = pvectsz[ptsz->index_hmf]
                                    *tau_profile_at_ell_1
                                    *tau_profile_at_ell_1;
  }
 else if (_kSZ_kSZ_2h_){
   // int index_l = (int) pvectsz[ptsz->index_multipole];
   // double l1 = ptsz->ell[index_l];
   evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
   evaluate_tau_profile(kl,pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];

   pvectsz[ptsz->index_integrand] = pvectsz[ptsz->index_hmf]
                                    *pvectsz[ptsz->index_halo_bias]
                                    *tau_profile_at_ell_1;
  }


 else if (_tSZ_tSZ_tSZ_1halo_){
   // int index_l = (int) pvectsz[ptsz->index_multipole];
   // double l1 = ptsz->ell[index_l];
   // double l2 = ptsz->bispectrum_lambda_k2*ptsz->ell[index_l];
   // double l3 = ptsz->bispectrum_lambda_k3*ptsz->ell[index_l];

   // pvectsz[ptsz->index_multipole_for_pressure_profile] = l1;
   // double kl1 =
   evaluate_pressure_profile(kl1,pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_1 = pvectsz[ptsz->index_pressure_profile];

   // pvectsz[ptsz->index_multipole_for_pressure_profile] = l2;
   evaluate_pressure_profile(kl2,pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_2 = pvectsz[ptsz->index_pressure_profile];

   // pvectsz[ptsz->index_multipole_for_pressure_profile] = l3;
   evaluate_pressure_profile(kl3,pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_3 = pvectsz[ptsz->index_pressure_profile];

       pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                         pvectsz[ptsz->index_hmf]
                                         // *pvectsz[ptsz->index_dlnMdeltadlnM]
                                         // //*pvectsz[ptsz->index_completeness]
                                         *pressure_profile_at_ell_1
                                         *pressure_profile_at_ell_2
                                         *pressure_profile_at_ell_3;
     // printf("l1 = %.3e\n",pressure_profile_at_ell_1);
     // printf("l2 = %.3e\n",pressure_profile_at_ell_2);
     // printf("l3 = %.3e\n",pressure_profile_at_ell_3);


  }


 else if (_tSZ_tSZ_tSZ_2h_){
   // int index_l = (int) pvectsz[ptsz->index_multipole];
   // double l1 = ptsz->ell[index_l];
   // double l2 = ptsz->bispectrum_lambda_k2*ptsz->ell[index_l];
   // double l3 = ptsz->bispectrum_lambda_k3*ptsz->ell[index_l];

   // pvectsz[ptsz->index_multipole_for_pressure_profile] = l1;
   evaluate_pressure_profile(kl1,pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_1 = pvectsz[ptsz->index_pressure_profile];
   // pvectsz[ptsz->index_multipole_for_pressure_profile] = l2;
   evaluate_pressure_profile(kl2,pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_2 = pvectsz[ptsz->index_pressure_profile];
   // pvectsz[ptsz->index_multipole_for_pressure_profile] = l3;
   evaluate_pressure_profile(kl3,pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_3 = pvectsz[ptsz->index_pressure_profile];

   evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {
       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_halo_bias]
                                         *pressure_profile_at_ell_1
                                         *pressure_profile_at_ell_2;
                                       }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {
       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_halo_bias]
                                         *pressure_profile_at_ell_3;
                                       }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  3) {
       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_halo_bias]
                                         *pressure_profile_at_ell_1
                                         *pressure_profile_at_ell_3;

                                       }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  4) {
       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_halo_bias]
                                         *pressure_profile_at_ell_2;
                                       }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  5) {
       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_halo_bias]
                                         *pressure_profile_at_ell_2
                                         *pressure_profile_at_ell_3;

                                       }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  6) {
       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_halo_bias]
                                         *pressure_profile_at_ell_1;
                                       }

  }

   else if (_tSZ_tSZ_tSZ_3h_){
     int index_l = (int) pvectsz[ptsz->index_multipole];
     // double l1 = ptsz->ell[index_l];
     // double l2 = ptsz->bispectrum_lambda_k2*ptsz->ell[index_l];
     // double l3 = ptsz->bispectrum_lambda_k3*ptsz->ell[index_l];

     // pvectsz[ptsz->index_multipole_for_pressure_profile] = l1;
     evaluate_pressure_profile(kl1,pvecback,pvectsz,pba,ptsz);
     double pressure_profile_at_ell_1 = pvectsz[ptsz->index_pressure_profile];
     // pvectsz[ptsz->index_multipole_for_pressure_profile] = l2;
     evaluate_pressure_profile(kl2,pvecback,pvectsz,pba,ptsz);
     double pressure_profile_at_ell_2 = pvectsz[ptsz->index_pressure_profile];
     // pvectsz[ptsz->index_multipole_for_pressure_profile] = l3;
     evaluate_pressure_profile(kl3,pvecback,pvectsz,pba,ptsz);
     double pressure_profile_at_ell_3 = pvectsz[ptsz->index_pressure_profile];


     evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
     evaluate_halo_bias_b2(pvecback,pvectsz,pba,ppm,pnl,ptsz);

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias]
                                           *pressure_profile_at_ell_1;
                                         }

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias]
                                           *pressure_profile_at_ell_2;
                                         }

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  3) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias]
                                           *pressure_profile_at_ell_3;

                                         }

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  4) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias_b2]
                                           *pressure_profile_at_ell_1;
                                         }


  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  5) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias_b2]
                                           *pressure_profile_at_ell_2;

                                         }

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  6) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias_b2]
                                           *pressure_profile_at_ell_2;

                                         }

    }



 else if (_kSZ_kSZ_tSZ_1h_){
   // int index_l = (int) pvectsz[ptsz->index_multipole];
   // double l1 = ptsz->ell[index_l];
   // double l2 = ptsz->bispectrum_lambda_k2*ptsz->ell[index_l];
   // double l3 = ptsz->bispectrum_lambda_k3*ptsz->ell[index_l];

   // pvectsz[ptsz->index_multipole_for_pressure_profile] = l1;//the actual multipole
   evaluate_pressure_profile(kl1,pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_1 = pvectsz[ptsz->index_pressure_profile];

   evaluate_tau_profile(kl2,pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];

   evaluate_tau_profile(kl3,pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_3 = pvectsz[ptsz->index_tau_profile];

       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pressure_profile_at_ell_1
                                         *tau_profile_at_ell_2
                                         *tau_profile_at_ell_3;
  }

 else if (_kSZ_kSZ_tSZ_2h_){
   // int index_l = (int) pvectsz[ptsz->index_multipole];
   // double l1 = ptsz->ell[index_l];
   // double l2 = ptsz->bispectrum_lambda_k2*ptsz->ell[index_l];
   // double l3 = ptsz->bispectrum_lambda_k3*ptsz->ell[index_l];

   // pvectsz[ptsz->index_multipole_for_pressure_profile] = l1;
   evaluate_pressure_profile(kl1,pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_1 = pvectsz[ptsz->index_pressure_profile];
   // pvectsz[ptsz->index_multipole_for_pressure_profile] = l2;
   evaluate_pressure_profile(kl2,pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_2 = pvectsz[ptsz->index_pressure_profile];
   // pvectsz[ptsz->index_multipole_for_pressure_profile] = l3;
   evaluate_pressure_profile(kl3,pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_3 = pvectsz[ptsz->index_pressure_profile];

   evaluate_tau_profile(kl1,pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];
   evaluate_tau_profile(kl2,pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];
   evaluate_tau_profile(kl3,pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_3 = pvectsz[ptsz->index_tau_profile];


   evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {
       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_halo_bias]
                                         *pressure_profile_at_ell_1
                                         *tau_profile_at_ell_2;
                                       }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {
       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_halo_bias]
                                         *tau_profile_at_ell_3;
                                       }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  3) {
       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_halo_bias]
                                         *tau_profile_at_ell_1
                                         *pressure_profile_at_ell_3;

                                       }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  4) {
       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_halo_bias]
                                         *tau_profile_at_ell_2;
                                       }


if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  5) {
       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_halo_bias]
                                         *tau_profile_at_ell_2
                                         *tau_profile_at_ell_1;

                                       }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  6) {
       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_halo_bias]
                                         *pressure_profile_at_ell_3;
                                       }

  }


   else if (_kSZ_kSZ_tSZ_3h_){
     // int index_l = (int) pvectsz[ptsz->index_multipole];
     // double l1 = ptsz->ell[index_l];
     // double l2 = ptsz->bispectrum_lambda_k2*ptsz->ell[index_l];
     // double l3 = ptsz->bispectrum_lambda_k3*ptsz->ell[index_l];

     // pvectsz[ptsz->index_multipole_for_pressure_profile] = l1;
     evaluate_pressure_profile(kl1,pvecback,pvectsz,pba,ptsz);
     double pressure_profile_at_ell_1 = pvectsz[ptsz->index_pressure_profile];
     // pvectsz[ptsz->index_multipole_for_pressure_profile] = l2;
     evaluate_pressure_profile(kl2,pvecback,pvectsz,pba,ptsz);
     double pressure_profile_at_ell_2 = pvectsz[ptsz->index_pressure_profile];
     // pvectsz[ptsz->index_multipole_for_pressure_profile] = l3;
     evaluate_pressure_profile(kl3,pvecback,pvectsz,pba,ptsz);
     double pressure_profile_at_ell_3 = pvectsz[ptsz->index_pressure_profile];

     evaluate_tau_profile(kl1,pvecback,pvectsz,pba,ptsz);
     double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];
     evaluate_tau_profile(kl2,pvecback,pvectsz,pba,ptsz);
     double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];
     evaluate_tau_profile(kl3,pvecback,pvectsz,pba,ptsz);
     double tau_profile_at_ell_3 = pvectsz[ptsz->index_tau_profile];


     evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
     evaluate_halo_bias_b2(pvecback,pvectsz,pba,ppm,pnl,ptsz);

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias]
                                           *pressure_profile_at_ell_1;
                                         }

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias]
                                           *pressure_profile_at_ell_2;
                                         }

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  3) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias]
                                           *pressure_profile_at_ell_3;

                                         }

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  4) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias]
                                           *tau_profile_at_ell_1;
                                         }


  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  5) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias]
                                           *tau_profile_at_ell_2;

                                         }

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  6) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias]
                                           *tau_profile_at_ell_3;
                                         }


  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  7) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias_b2]
                                           *pressure_profile_at_ell_1;
                                         }

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  8) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias_b2]
                                           *pressure_profile_at_ell_2;
                                         }

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  9) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias_b2]
                                           *pressure_profile_at_ell_3;

                                         }

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  10) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias_b2]
                                           *tau_profile_at_ell_1;
                                         }


  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  11) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias_b2]
                                           *tau_profile_at_ell_2;

                                         }

  if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  12) {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias_b2]
                                           *tau_profile_at_ell_3;
                                         }

    }



  else if (_kSZ_kSZ_gal_1h_fft_){

    double l_min = ptsz->l_min_samp_fftw;
    double l_max = ptsz->l_max_samp_fftw; // this is a precision parameter
        // tabulate the integrand in the "l" dimension:
        const int N = ptsz->N_samp_fftw;
        double k[N], Pk[N],Pko[N];
        double lnk[N],lnpk[N];
        int ik;
        double fl;
        double taul;
        double l;
        double m = exp(logM);

        for (ik=0; ik<N; ik++){
        k[ik] = exp(log(l_min)+ik/(N-1.)*(log(l_max)-log(l_min)));
        lnk[ik] = log(l_min)+ik/(N-1.)*(log(l_max)-log(l_min));
        l = k[ik];
        fl = get_ksz_filter_at_l(l,ptsz);

        // pvectsz[ptsz->index_multipole_for_tau_profile] = l;
        evaluate_tau_profile((l+0.5)/chi,pvecback,pvectsz,pba,ptsz);
        taul = pvectsz[ptsz->index_tau_profile];
        Pk[ik] = fl*taul;
        Pko[ik] = fl*taul;
        }

        double r[N], xi[N];

        xi2pk(N,k,Pk,r,xi,ptsz);
        for (ik=0; ik<N; ik++){
        // convolution:
        xi[ik] = pow(xi[ik],2.);
        }

        pk2xi(N,r,xi,k,Pk,ptsz);

        // evaluate at ell:
        int index_l = (int) pvectsz[ptsz->index_multipole];
        double ell = ptsz->ell[index_l];
        double f_tau_f_tau = pwl_value_1d(N,lnk,Pk,log(ell));

        // galaxy term:
        double g; // this has 1/ng_bar
        // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ell;

        double kp = (ell+0.5)/chi;
        evaluate_galaxy_profile_2h(kp,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
        g =   pvectsz[ptsz->index_galaxy_profile];


        double dn = pvectsz[ptsz->index_hmf];
        double integrand_projected_fields = dn*f_tau_f_tau*g;


        if (isnan(integrand_projected_fields) || isinf(integrand_projected_fields)){
        integrand_projected_fields = 0.;
        }

        pvectsz[ptsz->index_integrand] = integrand_projected_fields;
  }




  else if (_kSZ_kSZ_gallens_1h_fft_ || _kSZ_kSZ_lens_1h_fft_){

    double l_min = ptsz->l_min_samp_fftw;
    double l_max = ptsz->l_max_samp_fftw; // this is a precision parameter
        // tabulate the integrand in the "l" dimension:
        const int N = ptsz->N_samp_fftw;
        double k[N], Pk[N],Pko[N];
        double lnk[N],lnpk[N];
        int ik;
        double fl;
        double taul;
        double l;
        double m = exp(logM);

        for (ik=0; ik<N; ik++){
        k[ik] = exp(log(l_min)+ik/(N-1.)*(log(l_max)-log(l_min)));
        lnk[ik] = log(l_min)+ik/(N-1.)*(log(l_max)-log(l_min));
        l = k[ik];
        fl = get_ksz_filter_at_l(l,ptsz);

        // pvectsz[ptsz->index_multipole_for_tau_profile] = l;
        evaluate_tau_profile((l+0.5)/chi,pvecback,pvectsz,pba,ptsz);
        taul = pvectsz[ptsz->index_tau_profile];
        Pk[ik] = fl*taul;
        Pko[ik] = fl*taul;
        }

        double r[N], xi[N];

        xi2pk(N,k,Pk,r,xi,ptsz);
        for (ik=0; ik<N; ik++){
        // convolution:
        xi[ik] = pow(xi[ik],2.);
        }

        pk2xi(N,r,xi,k,Pk,ptsz);

        // evaluate at ell:
        int index_l = (int) pvectsz[ptsz->index_multipole];
        double ell = ptsz->ell[index_l];
        double f_tau_f_tau = pwl_value_1d(N,lnk,Pk,log(ell));

        // galaxy term:
        double g; // this has 1/ng_bar
        // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ell;

        double kp = (ell+0.5)/chi;
        evaluate_lensing_profile(kp,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
        g =   pvectsz[ptsz->index_lensing_profile];


        double dn = pvectsz[ptsz->index_hmf];
        double integrand_projected_fields = dn*f_tau_f_tau*g;


        if (isnan(integrand_projected_fields) || isinf(integrand_projected_fields)){
        integrand_projected_fields = 0.;
        }

        pvectsz[ptsz->index_integrand] = integrand_projected_fields;
  }

  else if (_kSZ_kSZ_gal_1h_){

   int index_theta_1 = (int) pvectsz[ptsz->index_multipole_1];
   double theta_1 = ptsz->theta_kSZ2_gal_theta_grid[index_theta_1];
   int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
   int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
   double l2 = exp(ptsz->ell_kSZ2_gal_multipole_grid[index_l_2]);
   double l3 = ptsz->ell[index_l_3];
   double ell = l3;
   double ell_prime = l2;
   double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos(theta_1));

   // pvectsz[ptsz->index_multipole_for_tau_profile] = ptsz->ell_kSZ2_gal_multipole_grid[index_l_1];
   // pvectsz[ptsz->index_multipole_for_tau_profile] = l1;//ptsz->ell_kSZ2_gal_multipole_grid[index_l_1];
   double k1p,k2p,k3p;
   k1p = (l1+0.5)/chi;
   k2p = (l2+0.5)/chi;
   k3p = (l3+0.5)/chi;

   evaluate_tau_profile((l1+0.5)/chi,pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];

   // int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
   // pvectsz[ptsz->index_multipole_for_tau_profile] = l2;//ptsz->ell_kSZ2_gal_multipole_grid[index_l_2];
   evaluate_tau_profile((l2+0.5)/chi,pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];

   // int index_l_3 = (int) pvectsz[ptsz->index_multipole_3]; // this is the ell of the power spectrum
   pvectsz[ptsz->index_multipole_for_galaxy_profile] = l3;//ptsz->ell[index_l_3];
   // double chi = sqrt(pvectsz[ptsz->index_chi2]);
   double kp = (l3+0.5)/chi;
   evaluate_galaxy_profile_2h(kp,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
   double galaxy_profile_at_ell_3 = pvectsz[ptsz->index_galaxy_profile];



   pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                     *tau_profile_at_ell_1
                                     *tau_profile_at_ell_2
                                     *galaxy_profile_at_ell_3;


////configurations::
if (ptsz->bispec_conf_id!=0){
double lambda1p  = fabs(k1p/k3p);
// k1 = lambda1*k3
double lambda2p  = fabs(k2p/k3p);
double flag_conf = 0.;
if (ptsz->bispec_conf_id==1){ // equi
// just equilateral:
if ((fabs(lambda1p-1.)<0.3) && (fabs(lambda2p-1.)<0.3)){
flag_conf = 1.;
}
}
else if (ptsz->bispec_conf_id==2){ //squeezed
// // just flattened:
// if ((fabs(lambda1p-(1.-lambda2p))<0.1)){
// flag_conf = 1.;
// }
// just squeezed:
if ((fabs(lambda1p)>3.) && (fabs(lambda2p)>3.) && fabs((lambda1p/lambda2p)-1.)<0.1){
flag_conf = 1.;
}
}
else if (ptsz->bispec_conf_id==3){ //anti-squeezed
flag_conf = 1.;
// // just flattened:
// if ((fabs(lambda1p-(1.-lambda2p))<0.1)){
// flag_conf = 1.;
// }
// just squeezed:
if ((fabs(lambda1p)>3.) && (fabs(lambda2p)>3.) && fabs((lambda1p/lambda2p)-1.)<0.1){
flag_conf = 0.;
}
}

pvectsz[ptsz->index_integrand] *= flag_conf;
////end configurations.
}


  }

  else if (_kSZ_kSZ_gal_2h_fft_){
  // else if (1==0){
// double l_min = 1e-2;
// double l_max = 2e5; // this is a precision parameter
double l_min = ptsz->l_min_samp_fftw;
double l_max = ptsz->l_max_samp_fftw; // this is a precision parameter
// tabulate the integrand in the "l" dimension:
const int N = ptsz->N_samp_fftw;
double k[N], Pk[N],Pko[N];
double lnk[N],lnpk[N];
int ik;
double fl;
double taul;
double l;
// double m = exp(logM);

for (ik=0; ik<N; ik++){
k[ik] = exp(log(l_min)+ik/(N-1.)*(log(l_max)-log(l_min)));
lnk[ik] = log(l_min)+ik/(N-1.)*(log(l_max)-log(l_min));
l = k[ik];
fl = get_ksz_filter_at_l(l,ptsz);

// pvectsz[ptsz->index_multipole_for_tau_profile] = l;//ptsz->ell_kSZ2_gal_multipole_grid[index_l_2];
evaluate_tau_profile((l+0.5)/chi,pvecback,pvectsz,pba,ptsz);
taul = pvectsz[ptsz->index_tau_profile];//get_tau_profile_at_z_m_l(z,m,l,ptsz,pba);
Pk[ik] = fl*taul;
Pko[ik] = fl*taul;
// if(l>3e3)
//   printf("k = %.5e pk = %.5e\n",l,Pk[ik]);
}

double r[N], xi[N];

// printf("\n\n####################\n");
// printf("To position space:\n");
// printf("####################\n\n");
// pk2xi(N,k,Pk,r,xi,ptsz);
xi2pk(N,k,Pk,r,xi,ptsz);
for (ik=0; ik<N; ik++){
// printf("r = %.5e xi = %.5e\n",r[ik],xi[ik]);
// convolution:
xi[ik] = pow(xi[ik],2.);
}

// printf("\n\n####################\n");
// printf("Back to k space:\n");
// printf("####################\n\n");
pk2xi(N,r,xi,k,Pk,ptsz);
// xi2pk(N,r,xi,k,Pk,ptsz);
// xi2pk(N,r,xi,k,Pk,ptsz);

// for (ik=0; ik<N; ik++){
// // printf("k = %.5e pk rec = %.5e pk = %.5e\n",k[ik],Pk[ik],Pk[ik]/Pko[ik]);
// lnpk[ik] = log(fabs(Pk[ik]));
// }
   // evaluate at ell:
   int index_l = (int) pvectsz[ptsz->index_multipole];
   double ell = ptsz->ell[index_l];
   // printf("ok1\n");

   // double f_tau_f_tau = exp(pwl_value_1d(N,lnk,lnpk,log(ell)));
   double f_tau_f_tau = pwl_value_1d(N,lnk,Pk,log(ell));
   // printf("ok2\n");
   // galaxy term:
   double g;// = get_galaxy_profile_at_z_m_l_1h(z,m,ell,ptsz,pba); // this has 1/ng_bar
   // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ell;//ptsz->ell[index_l_3];
   // double chi = sqrt(pvectsz[ptsz->index_chi2]);
   double kp = (ell+0.5)/chi;
   evaluate_galaxy_profile_2h(kp,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
   g =   pvectsz[ptsz->index_galaxy_profile];
      // printf("ok3\n");
   // double ell = l3;
   double dn = pvectsz[ptsz->index_hmf];//get_dndlnM_at_z_and_M(z,m,ptsz);
      // printf("ok4\n");
   // for (ik=0; ik<N; ik++){

   // lnpk[ik] = log(Pk[ik]);
   // }
   evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
   double b1 = pvectsz[ptsz->index_halo_bias];
   // printf("k = %.5e f_tau_f_tau = %.5e g = %.5e dn = %.5e b1 = %.5e\n",ell,f_tau_f_tau,g,dn,b1);
   double integrand_projected_fields = dn*f_tau_f_tau*b1;





// --> integrand wrt ln(1+z), ln(m)

// printf("l = %.5e m = %.5e z = %.5e integrand = %.5e\n",ell,m,z,integrand_projected_fields);

if (isnan(integrand_projected_fields) || isinf(integrand_projected_fields)){
// printf("nans\n");
integrand_projected_fields = 0.;
}

pvectsz[ptsz->index_integrand] = integrand_projected_fields;
  }

  else if (_kSZ_kSZ_gallens_2h_fft_ || _kSZ_kSZ_lens_2h_fft_){
  // else if (1==0){
// double l_min = 1e-2;
// double l_max = 2e5; // this is a precision parameter
double l_min = ptsz->l_min_samp_fftw;
double l_max = ptsz->l_max_samp_fftw; // this is a precision parameter
// tabulate the integrand in the "l" dimension:
const int N = ptsz->N_samp_fftw;
double k[N], Pk[N],Pko[N];
double lnk[N],lnpk[N];
int ik;
double fl;
double taul;
double l;
// double m = exp(logM);

for (ik=0; ik<N; ik++){
k[ik] = exp(log(l_min)+ik/(N-1.)*(log(l_max)-log(l_min)));
lnk[ik] = log(l_min)+ik/(N-1.)*(log(l_max)-log(l_min));
l = k[ik];
fl = get_ksz_filter_at_l(l,ptsz);

// pvectsz[ptsz->index_multipole_for_tau_profile] = l;//ptsz->ell_kSZ2_gal_multipole_grid[index_l_2];
evaluate_tau_profile((l+0.5)/chi,pvecback,pvectsz,pba,ptsz);
taul = pvectsz[ptsz->index_tau_profile];//get_tau_profile_at_z_m_l(z,m,l,ptsz,pba);
Pk[ik] = fl*taul;
Pko[ik] = fl*taul;
// if(l>3e3)
//   printf("k = %.5e pk = %.5e\n",l,Pk[ik]);
}

double r[N], xi[N];

// printf("\n\n####################\n");
// printf("To position space:\n");
// printf("####################\n\n");
// pk2xi(N,k,Pk,r,xi,ptsz);
xi2pk(N,k,Pk,r,xi,ptsz);
for (ik=0; ik<N; ik++){
// printf("r = %.5e xi = %.5e\n",r[ik],xi[ik]);
// convolution:
xi[ik] = pow(xi[ik],2.);
}

// printf("\n\n####################\n");
// printf("Back to k space:\n");
// printf("####################\n\n");
pk2xi(N,r,xi,k,Pk,ptsz);
// xi2pk(N,r,xi,k,Pk,ptsz);
// xi2pk(N,r,xi,k,Pk,ptsz);

// for (ik=0; ik<N; ik++){
// // printf("k = %.5e pk rec = %.5e pk = %.5e\n",k[ik],Pk[ik],Pk[ik]/Pko[ik]);
// lnpk[ik] = log(fabs(Pk[ik]));
// }
   // evaluate at ell:
   int index_l = (int) pvectsz[ptsz->index_multipole];
   double ell = ptsz->ell[index_l];
   // printf("ok1\n");

   // double f_tau_f_tau = exp(pwl_value_1d(N,lnk,lnpk,log(ell)));
   double f_tau_f_tau = pwl_value_1d(N,lnk,Pk,log(ell));
   // printf("ok2\n");
   // galaxy term:
   double g;// = get_galaxy_profile_at_z_m_l_1h(z,m,ell,ptsz,pba); // this has 1/ng_bar
   // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ell;//ptsz->ell[index_l_3];
   // double chi = sqrt(pvectsz[ptsz->index_chi2]);
   double kp = (ell+0.5)/chi;
   evaluate_lensing_profile(kp,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
   g =   pvectsz[ptsz->index_lensing_profile];
      // printf("ok3\n");
   // double ell = l3;
   double dn = pvectsz[ptsz->index_hmf];//get_dndlnM_at_z_and_M(z,m,ptsz);
      // printf("ok4\n");
   // for (ik=0; ik<N; ik++){

   // lnpk[ik] = log(Pk[ik]);
   // }
   evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
   double b1 = pvectsz[ptsz->index_halo_bias];
   // printf("k = %.5e f_tau_f_tau = %.5e g = %.5e dn = %.5e b1 = %.5e\n",ell,f_tau_f_tau,g,dn,b1);
   double integrand_projected_fields = dn*f_tau_f_tau*b1;





// --> integrand wrt ln(1+z), ln(m)

// printf("l = %.5e m = %.5e z = %.5e integrand = %.5e\n",ell,m,z,integrand_projected_fields);

if (isnan(integrand_projected_fields) || isinf(integrand_projected_fields)){
// printf("nans\n");
integrand_projected_fields = 0.;
}

pvectsz[ptsz->index_integrand] = integrand_projected_fields;
  }

  else if (_kSZ_kSZ_gal_2h_){
   int index_theta_1 = (int) pvectsz[ptsz->index_multipole_1];
   double theta_1 = ptsz->theta_kSZ2_gal_theta_grid[index_theta_1];
   int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
   int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
   double l2 = exp(ptsz->ell_kSZ2_gal_multipole_grid[index_l_2]);
   double l3 = ptsz->ell[index_l_3];
   double ell = l3;
   double ell_prime = l2;
   double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos(theta_1));

   double k1p,k2p,k3p;
   k1p = (l1+0.5)/chi;
   k2p = (l2+0.5)/chi;
   k3p = (l3+0.5)/chi;


    // g3 - t1t2
    if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    // int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
    // pvectsz[ptsz->index_multipole_for_tau_profile] = l2;//ptsz->ell_kSZ2_gal_multipole_grid[index_l_2];
    evaluate_tau_profile((l2+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];

    // int index_l_1 = (int) pvectsz[ptsz->index_multipole_1];
    pvectsz[ptsz->index_multipole_for_tau_profile] = l1;//ptsz->ell_kSZ2_gal_multipole_grid[index_l_1];
    evaluate_tau_profile((l1+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *tau_profile_at_ell_2
                                      *tau_profile_at_ell_1;
    }

    if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    // int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
    pvectsz[ptsz->index_multipole_for_galaxy_profile] = l3;//ptsz->ell[index_l_3];
    // double chi = sqrt(pvectsz[ptsz->index_chi2]);
    double kp = (l3+0.5)/chi;
    evaluate_galaxy_profile_2h(kp,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_ell_3 = pvectsz[ptsz->index_galaxy_profile];


    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *galaxy_profile_at_ell_3;
    }

    // t2 - g3t1
    if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  3) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    // int index_l_2 = (int) pvectsz[ptsz->index_multipole_3];
    pvectsz[ptsz->index_multipole_for_tau_profile] = l2;//ptsz->ell_kSZ2_gal_multipole_grid[index_l_2];
    evaluate_tau_profile((l2+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *tau_profile_at_ell_2;
    }

    if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  4) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    // int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
    pvectsz[ptsz->index_multipole_for_galaxy_profile] = l3;//ptsz->ell[index_l_3];
    // double chi = sqrt(pvectsz[ptsz->index_chi2]);
    double kp = (l3+0.5)/chi;
    evaluate_galaxy_profile_2h(kp,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_ell_3 = pvectsz[ptsz->index_galaxy_profile];

    // int index_l_1 = (int) pvectsz[ptsz->index_multipole_1];
    pvectsz[ptsz->index_multipole_for_tau_profile] = l1;// ptsz->ell_kSZ2_gal_multipole_grid[index_l_1];
    evaluate_tau_profile((l1+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];


    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *galaxy_profile_at_ell_3
                                      *tau_profile_at_ell_1;
    }

    // t1 - t2g3
    if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  5) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    // pvectsz[ptsz->index_multipole_for_tau_profile] = l1;
    evaluate_tau_profile((l1+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *tau_profile_at_ell_1;
    }

    if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  6) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    // pvectsz[ptsz->index_multipole_for_tau_profile] = l2;
    evaluate_tau_profile((l2+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];


    // pvectsz[ptsz->index_multipole_for_galaxy_profile] = l3;
    // double chi = sqrt(pvectsz[ptsz->index_chi2]);
    double kp = (l3+0.5)/chi;
    evaluate_galaxy_profile_2h(kp,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_ell_3 = pvectsz[ptsz->index_galaxy_profile];


    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *galaxy_profile_at_ell_3
                                      *tau_profile_at_ell_2;
    }


////configurations::
if (ptsz->bispec_conf_id!=0){
double lambda1p  = fabs(k1p/k3p);
// k1 = lambda1*k3
double lambda2p  = fabs(k2p/k3p);
double flag_conf = 0.;
if (ptsz->bispec_conf_id==1){ // equi
// just equilateral:
if ((fabs(lambda1p-1.)<0.3) && (fabs(lambda2p-1.)<0.3)){
flag_conf = 1.;
}
}
else if (ptsz->bispec_conf_id==2){ //squeezed
// // just flattened:
// if ((fabs(lambda1p-(1.-lambda2p))<0.1)){
// flag_conf = 1.;
// }
// just squeezed:
if ((fabs(lambda1p)>3.) && (fabs(lambda2p)>3.) && fabs((lambda1p/lambda2p)-1.)<0.1){
flag_conf = 1.;
}
}
else if (ptsz->bispec_conf_id==3){ //anti-squeezed
flag_conf = 1.;
// // just flattened:
// if ((fabs(lambda1p-(1.-lambda2p))<0.1)){
// flag_conf = 1.;
// }
// just squeezed:
if ((fabs(lambda1p)>3.) && (fabs(lambda2p)>3.) && fabs((lambda1p/lambda2p)-1.)<0.1){
flag_conf = 0.;
}
}

pvectsz[ptsz->index_integrand] *= flag_conf;
////end configurations.
}




    }



  else if (_kSZ_kSZ_gal_3h_){
   int index_theta_1 = (int) pvectsz[ptsz->index_multipole_1];
   double theta_1 = ptsz->theta_kSZ2_gal_theta_grid[index_theta_1];
   // double cos_theta_1 = ptsz->theta_kSZ2_gal_theta_grid[index_theta_1];
   int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
   int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
   double l2 = exp(ptsz->ell_kSZ2_gal_multipole_grid[index_l_2]);
   double l3 = ptsz->ell[index_l_3];
   double ell = l3;
   double ell_prime = l2;
   double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos(theta_1));
   // double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos_theta_1);
   double k1p,k2p,k3p;
   k1p = (l1+0.5)/chi;
   k2p = (l2+0.5)/chi;
   k3p = (l3+0.5)/chi;
    // b1t1
    if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    // int index_l_1 = (int) pvectsz[ptsz->index_multipole_1];
    // pvectsz[ptsz->index_multipole_for_tau_profile] = l1;
    evaluate_tau_profile((l1+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *tau_profile_at_ell_1;

    // printf("b1t1 ell = %.3e tau = %.3e\n",pvectsz[ptsz->index_multipole_for_tau_profile],tau_profile_at_ell_1);
    }

    // b1t2
    if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    // int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
    // pvectsz[ptsz->index_multipole_for_tau_profile] = l2;
    evaluate_tau_profile((l2+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];


    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *tau_profile_at_ell_2;
    // printf("b1t2 ell = %.3e tau = %.3e\n",pvectsz[ptsz->index_multipole_for_tau_profile],tau_profile_at_ell_2);

    }

    // b1g3
    if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  3) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    // int index_l_3 = (int) pvectsz[ptsz->index_multipole_3]; // ell_3 is the ell of the power spectrum
    pvectsz[ptsz->index_multipole_for_galaxy_profile] = l3;
    // double chi = sqrt(pvectsz[ptsz->index_chi2]);
    double kp = (l3+0.5)/chi;
    evaluate_galaxy_profile_2h(kp,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_ell_3 = pvectsz[ptsz->index_galaxy_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *galaxy_profile_at_ell_3;

    // printf("b1g3 ell = %.3e tau = %.3e\n",pvectsz[ptsz->index_multipole_for_galaxy_profile],galaxy_profile_at_ell_3);

    }

    // b2g3
    if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  4) {
    evaluate_halo_bias_b2(pvecback,pvectsz,pba,ppm,pnl,ptsz);

    // int index_l_3 = (int) pvectsz[ptsz->index_multipole_3]; // ell_3 is the ell of the power spectrum
    pvectsz[ptsz->index_multipole_for_galaxy_profile] = l3;
    // double chi = sqrt(pvectsz[ptsz->index_chi2]);
    double kp = (l3+0.5)/chi;
    evaluate_galaxy_profile_2h(kp,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_ell_3 = pvectsz[ptsz->index_galaxy_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias_b2]
                                      *galaxy_profile_at_ell_3;

    //printf("b2g3 ell = %.3e tau = %.3e\n",pvectsz[ptsz->index_multipole_for_galaxy_profile],galaxy_profile_at_ell_3);

      }


////configurations::
if (ptsz->bispec_conf_id!=0){
double lambda1p  = fabs(k1p/k3p);
// k1 = lambda1*k3
double lambda2p  = fabs(k2p/k3p);
double flag_conf = 0.;
if (ptsz->bispec_conf_id==1){ // equi
// just equilateral:
if ((fabs(lambda1p-1.)<0.3) && (fabs(lambda2p-1.)<0.3)){
flag_conf = 1.;
}
}
else if (ptsz->bispec_conf_id==2){ //squeezed
// // just flattened:
// if ((fabs(lambda1p-(1.-lambda2p))<0.1)){
// flag_conf = 1.;
// }
// just squeezed:
if ((fabs(lambda1p)>3.) && (fabs(lambda2p)>3.) && fabs((lambda1p/lambda2p)-1.)<0.1){
flag_conf = 1.;
}
}
else if (ptsz->bispec_conf_id==3){ //anti-squeezed
flag_conf = 1.;
// // just flattened:
// if ((fabs(lambda1p-(1.-lambda2p))<0.1)){
// flag_conf = 1.;
// }
// just squeezed:
if ((fabs(lambda1p)>3.) && (fabs(lambda2p)>3.) && fabs((lambda1p/lambda2p)-1.)<0.1){
flag_conf = 0.;
}
}

pvectsz[ptsz->index_integrand] *= flag_conf;
////end configurations.
}


    }




 else if (_kSZ_kSZ_lensmag_1halo_){
   int index_theta_1 = (int) pvectsz[ptsz->index_multipole_1];
   double theta_1 = ptsz->theta_kSZ2_gal_theta_grid[index_theta_1];
   int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
   int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
   double l2 = exp(ptsz->ell_kSZ2_gal_multipole_grid[index_l_2]);
   double l3 = ptsz->ell[index_l_3];
   double ell = l3;
   double ell_prime = l2;
   double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos(theta_1));

   // int index_l_1 = (int) pvectsz[ptsz->index_multipole_1];
   // pvectsz[ptsz->index_multipole_for_tau_profile] = l1;//ptsz->ell_kSZ2_gal_multipole_grid[index_l_1];//the actual multipole
   evaluate_tau_profile((l1+0.5)/chi,pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];

   // int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
   // pvectsz[ptsz->index_multipole_for_tau_profile] = l2;//ptsz->ell_kSZ2_gal_multipole_grid[index_l_2];
   evaluate_tau_profile((l2+0.5)/chi,pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];

   // int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
   // pvectsz[ptsz->index_multipole_for_lensing_profile] = l3;//ptsz->ell[index_l_3];
   evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
   double lensing_profile_at_ell_3 = pvectsz[ptsz->index_lensing_profile];

   // //velocity dispersion (kSZ)
   // evaluate_vrms2(pvecback,pvectsz,pba,pnl,ptsz);

       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *tau_profile_at_ell_1
                                         *tau_profile_at_ell_2
                                         *lensing_profile_at_ell_3;
   // printf("intg = %.3e\n",pvectsz[ptsz->index_integrand]);
   // printf("hmf = %.3e\n",pvectsz[ptsz->index_hmf]);
   // printf("tau1 = %.3e\n",tau_profile_at_ell_1);
   // printf("tau2 = %.3e\n",tau_profile_at_ell_2);
   // printf("lens3 = %.3e\n",lensing_profile_at_ell_3);


  }


 else if (_tSZ_lensmag_1h_){

   // int index_l = (int) pvectsz[ptsz->index_multipole];
   // pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
   evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_1 = pvectsz[ptsz->index_pressure_profile];

   // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
   evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
   double lensing_profile_at_ell_2 = pvectsz[ptsz->index_lensing_profile];

       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *pressure_profile_at_ell_1
                                         *lensing_profile_at_ell_2
                                         *damping_1h_term;
  }

   else if  (_tSZ_lensmag_2h_){



             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
             evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_pressure_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {


             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

          }



  else if (_pk_gg_at_z_1h_){

    evaluate_galaxy_profile_1h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_k_1 = pvectsz[ptsz->index_galaxy_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *galaxy_profile_at_k_1
                                      *galaxy_profile_at_k_1
                                      *damping_1h_term;
   }

  else if (_pk_gg_at_z_2h_){

    evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_k_1 = pvectsz[ptsz->index_galaxy_profile];
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *galaxy_profile_at_k_1;
   }

  else if (_pk_bb_at_z_1h_){

    double gas_profile_at_k_1 = get_gas_density_profile_at_k_M_z(kl,m_delta_electron_density,z,ptsz);

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *gas_profile_at_k_1
                                      *gas_profile_at_k_1
                                      *pow((pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0,-2)
                                      *damping_1h_term;
   }

  else if (_pk_bb_at_z_2h_){

    double gas_profile_at_k_1 = get_gas_density_profile_at_k_M_z(kl,m_delta_electron_density,z,ptsz);
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *gas_profile_at_k_1
                                      *pow((pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0,-1);
   }

  else if (_pk_b_at_z_2h_){

    double gas_profile_at_k_1 = get_gas_density_profile_at_k_M_z(kl,m_delta_electron_density,z,ptsz);
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *gas_profile_at_k_1
                                      *pow((pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0,0);
   }


  else if (_pk_em_at_z_1h_){

    double gas_profile_at_k_1 = get_gas_density_profile_at_k_M_z(kl,m_delta_electron_density,z,ptsz);

    evaluate_matter_density_profile(kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_1 = pvectsz[ptsz->index_density_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *density_profile_at_k_1
                                      *gas_profile_at_k_1
                                      *ptsz->f_free
                                      *pow((pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0,-1)
                                      *damping_1h_term;
   }

  else if (_pk_em_at_z_2h_){
 if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {
    double gas_profile_at_k_1 = get_gas_density_profile_at_k_M_z(kl,m_delta_electron_density,z,ptsz);
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *gas_profile_at_k_1
                                      *ptsz->f_free
                                      *pow((pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0,-1);
                                    }
if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {
    evaluate_matter_density_profile(kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_1 = pvectsz[ptsz->index_density_profile];
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *density_profile_at_k_1;

}

   }

  else if (_pk_HI_at_z_1h_){

    double HI_profile_at_k_1 = get_HI_density_profile_at_k_M_z(kl,m_delta_HI_density,z,ptsz);
    // printf("%.5e %.5e %.5e\n",pvectsz[ptsz->index_hmf],HI_profile_at_k_1,damping_1h_term);
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *HI_profile_at_k_1
                                      *HI_profile_at_k_1
                                      *pow((pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0,-2)
                                      *damping_1h_term;
   }

  else if (_pk_HI_at_z_2h_){

    double HI_profile_at_k_1 = get_HI_density_profile_at_k_M_z(kl,m_delta_HI_density,z,ptsz);
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *HI_profile_at_k_1
                                      *pow((pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0,-1);
   }



  else if (_pk_at_z_1h_){

    evaluate_matter_density_profile(kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_1 = pvectsz[ptsz->index_density_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      // *pvectsz[ptsz->index_mass_for_hmf]
                                      // *pvectsz[ptsz->index_mass_for_hmf]
                                      *density_profile_at_k_1
                                      *density_profile_at_k_1
                                      *damping_1h_term;
   }

  else if (_pk_at_z_2h_){

    evaluate_matter_density_profile(kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_1 = pvectsz[ptsz->index_density_profile];
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      // *pvectsz[ptsz->index_mass_for_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *density_profile_at_k_1;
   }
  else if (_bk_at_z_1h_){

    evaluate_matter_density_profile(kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_1 = pvectsz[ptsz->index_density_profile];
    evaluate_matter_density_profile(ptsz->bispectrum_lambda_k2*kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_2 = pvectsz[ptsz->index_density_profile];
    evaluate_matter_density_profile(ptsz->bispectrum_lambda_k3*kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_3 = pvectsz[ptsz->index_density_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *density_profile_at_k_1
                                      *density_profile_at_k_2
                                      *density_profile_at_k_3
                                      *damping_1h_term;
   }

  else if (_bk_at_z_2h_){

    evaluate_matter_density_profile(kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_1 = pvectsz[ptsz->index_density_profile];
    evaluate_matter_density_profile(ptsz->bispectrum_lambda_k2*kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_2 = pvectsz[ptsz->index_density_profile];
    evaluate_matter_density_profile(ptsz->bispectrum_lambda_k3*kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_3 = pvectsz[ptsz->index_density_profile];

    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *density_profile_at_k_1;
                                    }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *density_profile_at_k_2
                                      *density_profile_at_k_3;
                                    }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  3) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *density_profile_at_k_3;
                                    }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  4) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *density_profile_at_k_1
                                      *density_profile_at_k_2;
                                    }
if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  5) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *density_profile_at_k_2;
                                    }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  6) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *density_profile_at_k_3
                                      *density_profile_at_k_1;
                                    }


   }

  else if (_bk_at_z_3h_){

    evaluate_matter_density_profile(kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_1 = pvectsz[ptsz->index_density_profile];
    evaluate_matter_density_profile(ptsz->bispectrum_lambda_k2*kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_2 = pvectsz[ptsz->index_density_profile];
    evaluate_matter_density_profile(ptsz->bispectrum_lambda_k3*kl,r_delta_matter,c_delta_matter,pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_3 = pvectsz[ptsz->index_density_profile];

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *density_profile_at_k_1;
}

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {
    evaluate_halo_bias_b2(pvecback,pvectsz,pba,ppm,pnl,ptsz);
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias_b2]
                                      *density_profile_at_k_1;
}

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  6) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *density_profile_at_k_2;
}

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  7) {
    evaluate_halo_bias_b2(pvecback,pvectsz,pba,ppm,pnl,ptsz);
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias_b2]
                                      *density_profile_at_k_2;
}

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  8) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *density_profile_at_k_3;
}

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  9) {
    evaluate_halo_bias_b2(pvecback,pvectsz,pba,ppm,pnl,ptsz);
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias_b2]
                                      *density_profile_at_k_3;
}


if (ptsz->check_consistency_conditions == 1){
// mass consistency:
if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  3) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *density_profile_at_k_1;
}
// b1 consistency:
if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  4) {
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *density_profile_at_k_1;
}

// b2 consistency:
if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  5) {
    evaluate_halo_bias_b2(pvecback,pvectsz,pba,ppm,pnl,ptsz);
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias_b2]
                                      *density_profile_at_k_1;
}
}
   }

  else if (_bk_ttg_at_z_1h_){
    // double chi = sqrt(pvectsz[ptsz->index_chi2]);
    double l1,l2,l3;
    l1 = kl*chi-0.5;
    l2 = ptsz->bispectrum_lambda_k2*kl*chi-0.5;
    l3 = ptsz->bispectrum_lambda_k3*kl*chi-0.5;

    // if (l1<0.) l1 =1e-100;
    // if (l2<0.) l2 =1e-100;
    // if (l3<0.) l3 =1e-100;
    //


    // pvectsz[ptsz->index_multipole_for_tau_profile] = l1;
    evaluate_tau_profile((l1+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];


    // pvectsz[ptsz->index_multipole_for_tau_profile] = l2;
    evaluate_tau_profile((l2+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];


    // pvectsz[ptsz->index_multipole_for_galaxy_profile] = l3;

    double kp = (l3+0.5)/chi;
    evaluate_galaxy_profile_2h(kp,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_ell_3 = pvectsz[ptsz->index_galaxy_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *tau_profile_at_ell_1
                                      *tau_profile_at_ell_2
                                      *galaxy_profile_at_ell_3
                                      *damping_1h_term;

    if (isnan(pvectsz[ptsz->index_integrand]) || isinf(pvectsz[ptsz->index_integrand])){
    printf("z %.8e rg %.8e cg %.8e t %.8e t %.8e g %.8e\n",pvectsz[ptsz->index_z],r_delta_gal,c_delta_gal,tau_profile_at_ell_1,tau_profile_at_ell_2,galaxy_profile_at_ell_3);
    exit(0);
  }
   }

  else if (_bk_ttg_at_z_2h_){
    // double chi = sqrt(pvectsz[ptsz->index_chi2]);
    double l1,l2,l3;
    l1 = kl*chi-0.5;
    l2 = ptsz->bispectrum_lambda_k2*kl*chi-0.5;
    l3 = ptsz->bispectrum_lambda_k3*kl*chi-0.5;
    // if (l1<0.) l1 =1e-100;
    // if (l2<0.) l2 =1e-100;
    // if (l3<0.) l3 =1e-100;
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

    // pvectsz[ptsz->index_multipole_for_tau_profile] = l1;
    evaluate_tau_profile((l1+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];

    // pvectsz[ptsz->index_multipole_for_tau_profile] = l2;
    evaluate_tau_profile((l2+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];


    pvectsz[ptsz->index_multipole_for_galaxy_profile] = l3;
    double kp = (l3+0.5)/chi;
    evaluate_galaxy_profile_2h(kp,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_ell_3 = pvectsz[ptsz->index_galaxy_profile];


if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *tau_profile_at_ell_1;
                                    }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *tau_profile_at_ell_2;
                                    }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  3) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *tau_profile_at_ell_1
                                      *galaxy_profile_at_ell_3;
                                    }
if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  4) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *tau_profile_at_ell_2
                                      *galaxy_profile_at_ell_3;
                                    }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  5) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *galaxy_profile_at_ell_3;
                                    }

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  6) {
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *tau_profile_at_ell_1
                                      *tau_profile_at_ell_2;
                                    }

   }

  else if (_bk_ttg_at_z_3h_){
    // double chi = sqrt(pvectsz[ptsz->index_chi2]);
    double l1,l2,l3;
    l1 = kl*chi-0.5;
    l2 = ptsz->bispectrum_lambda_k2*kl*chi-0.5;
    l3 = ptsz->bispectrum_lambda_k3*kl*chi-0.5;
    // if (l1<0.) l1 =1e-100;
    // if (l2<0.) l2 =1e-100;
    // if (l3<0.) l3 =1e-100;
    // pvectsz[ptsz->index_multipole_for_tau_profile] = l1;
    evaluate_tau_profile((l1+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];

    // pvectsz[ptsz->index_multipole_for_tau_profile] = l2;
    evaluate_tau_profile((l2+0.5)/chi,pvecback,pvectsz,pba,ptsz);
    double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];


    pvectsz[ptsz->index_multipole_for_galaxy_profile] = l3;
    double kp = (l3+0.5)/chi;
    evaluate_galaxy_profile_2h(kp,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_ell_3 = pvectsz[ptsz->index_galaxy_profile];

    evaluate_halo_bias_b2(pvecback,pvectsz,pba,ppm,pnl,ptsz);
    evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
    // pvectsz[ptsz->index_halo_bias] = 1.;

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *tau_profile_at_ell_1;
}

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *tau_profile_at_ell_2;
}

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  3) {

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias]
                                      *galaxy_profile_at_ell_3;
}

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  4) {

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias_b2]
                                      *tau_profile_at_ell_1;
}

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  5) {

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias_b2]
                                      *tau_profile_at_ell_2;
}

if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  6) {

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_halo_bias_b2]
                                      *galaxy_profile_at_ell_3;
}
   }



  else if (_gal_gal_1h_){

    // int index_l = (int) pvectsz[ptsz->index_multipole];
    // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
    evaluate_galaxy_profile_1h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_ell_1 = pvectsz[ptsz->index_galaxy_profile];

        pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                          *galaxy_profile_at_ell_1
                                          *galaxy_profile_at_ell_1
                                          *damping_1h_term;

  // printf("hmf = %.5e ug = %.5e\n",pvectsz[ptsz->index_hmf],galaxy_profile_at_ell_1);
   }

   else if (_gal_gal_2h_){

     // int index_l = (int) pvectsz[ptsz->index_multipole];
     // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
     evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
     double galaxy_profile_at_ell_1 = pvectsz[ptsz->index_galaxy_profile];

     evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias] // BB: commented for debug
                                           *galaxy_profile_at_ell_1;
    // printf("galaxy_profile = %.3e\n",galaxy_profile_at_ell_1); // BB debug
    }

   else if (_gal_lens_1h_){
    // kl = 1e-1;
    // m_delta_gal = 5e13;
    // r_delta_gal = 0.3;
    // c_delta_gal = 8.;
    // pvectsz[ptsz->index_z] = 0.3;
    // pvectsz[ptsz->index_chi2] = 13.4;
    // pvecback[pba->index_bg_ang_distance] = 789.9;

             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
             // pvectsz[ptsz->index_multipole_for_galaxy_profile] =  ptsz->ell[index_l];
             evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    // printf("lens = %.8e gal = %.8e\n",pvectsz[ptsz->index_lensing_profile],pvectsz[ptsz->index_galaxy_profile]);
    //
    //          exit(0);
             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_galaxy_profile]
                                               *damping_1h_term;
          }
   else if (_gal_lens_2h_){


            if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

            evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
            int index_l = (int) pvectsz[ptsz->index_multipole];
            pvectsz[ptsz->index_multipole_for_galaxy_profile] =  ptsz->ell[index_l];
            evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
            pvectsz[ptsz->index_integrand] = pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_galaxy_profile]
                                             *pvectsz[ptsz->index_halo_bias];


            }

            if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

            evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

            // int index_l = (int) pvectsz[ptsz->index_multipole];


            // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
            evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);

            pvectsz[ptsz->index_integrand] = pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_lensing_profile]
                                             *pvectsz[ptsz->index_halo_bias];
            }
        }


  else if (_ngal_lens_1h_){
    // kl = 1e-1;
    // m_delta_gal = 5e13;
    // r_delta_gal = 0.3;
    // c_delta_gal = 8.;
    // pvectsz[ptsz->index_z] = 0.3;
    // pvectsz[ptsz->index_chi2] = 13.4;
    // pvecback[pba->index_bg_ang_distance] = 789.9;
    evaluate_galaxy_profile_ngal(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
    // printf("lens = %.8e gal = %.8e\n",pvectsz[ptsz->index_lensing_profile],pvectsz[ptsz->index_galaxy_profile]);
    // exit(0);
    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_lensing_profile]
                                      *pvectsz[ptsz->index_galaxy_profile]
                                      *damping_1h_term;
   }

  else if (_ngal_lens_2h_){

            // pvectsz[ptsz->index_part_id_cov_hsv] = 1;
            if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

            evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
            evaluate_galaxy_profile_ngal(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
            pvectsz[ptsz->index_integrand] = pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_galaxy_profile]
                                             *pvectsz[ptsz->index_halo_bias];


            }

            // pvectsz[ptsz->index_part_id_cov_hsv] =  2

            if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

            evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
            evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);

            pvectsz[ptsz->index_integrand] = pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_lensing_profile]
                                             *pvectsz[ptsz->index_halo_bias];
            }

          }
  else if (_ngal_ngal_2h_){

            // pvectsz[ptsz->index_part_id_cov_hsv] = 1;
            // if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

            evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
            evaluate_galaxy_profile_ngal(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
            pvectsz[ptsz->index_integrand] = pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_galaxy_profile]
                                             *pvectsz[ptsz->index_halo_bias];


            // // }
            //
            // pvectsz[ptsz->index_part_id_cov_hsv] =  2
            //
            // // if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {
            //
            // evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
            // evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
            //
            // pvectsz[ptsz->index_integrand] = pvectsz[ptsz->index_hmf]
            //                                  *pvectsz[ptsz->index_lensing_profile]
            //                                  *pvectsz[ptsz->index_halo_bias];
            // // }

          }





           else if (_gal_lensmag_1h_){

             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
             evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
             double galaxy_profile_at_ell_1 = pvectsz[ptsz->index_galaxy_profile];

             // pvectsz[ptsz->index_multipole_for_lensing_profile] =  ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
             double lensing_profile_at_ell_2 = pvectsz[ptsz->index_lensing_profile];

                 pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                                   *galaxy_profile_at_ell_1
                                                   *lensing_profile_at_ell_2
                                                   *damping_1h_term;

                  }
           else if (_gal_lensmag_2h_){


// printf("tsz lensmag 2h %.3f\n",1.);
             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
             evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_galaxy_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

                }



           else if (_gal_gallens_1h_){

             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
             evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
             double galaxy_profile_at_ell_1 = pvectsz[ptsz->index_galaxy_profile];

             // pvectsz[ptsz->index_multipole_for_lensing_profile] =  ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
             double lensing_profile_at_ell_2 = pvectsz[ptsz->index_lensing_profile];

                 pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                                   *galaxy_profile_at_ell_1
                                                   *lensing_profile_at_ell_2
                                                   *damping_1h_term;
            if (isnan(pvectsz[ptsz->index_integrand])||isinf(pvectsz[ptsz->index_integrand])){
            printf("nan or inf in integrand gallens 1h\n");
            printf("%.5e %.5e %.5e\n",pvectsz[ptsz->index_hmf],galaxy_profile_at_ell_1,lensing_profile_at_ell_2);
            exit(0);
            }

                  }
           else if (_gal_gallens_2h_){


// printf("tsz lensmag 2h %.3f\n",1.);
             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
             evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_galaxy_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

                }


           else if (_gallens_lens_1h_){

             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
             // evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
             // double galaxy_profile_at_ell_1 = pvectsz[ptsz->index_galaxy_profile];

             // pvectsz[ptsz->index_multipole_for_lensing_profile] =  ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
             double lensing_profile_at_ell_2 = pvectsz[ptsz->index_lensing_profile];

                 pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                                   *lensing_profile_at_ell_2
                                                   *lensing_profile_at_ell_2
                                                   *damping_1h_term;
            if (isnan(pvectsz[ptsz->index_integrand])||isinf(pvectsz[ptsz->index_integrand])){
            printf("nan or inf in integrand gallens lens 1h\n");
            printf("%.5e %.5e %.5e\n",pvectsz[ptsz->index_hmf],lensing_profile_at_ell_2,lensing_profile_at_ell_2);
            exit(0);
            }

                  }
           else if (_gallens_lens_2h_){


             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_halo_bias];


                }


           else if (_gallens_gallens_1h_){

             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
             // evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
             // double galaxy_profile_at_ell_1 = pvectsz[ptsz->index_galaxy_profile];

             // pvectsz[ptsz->index_multipole_for_lensing_profile] =  ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
             double lensing_profile_at_ell_2 = pvectsz[ptsz->index_lensing_profile];

                 pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                                   *lensing_profile_at_ell_2
                                                   *lensing_profile_at_ell_2
                                                   *damping_1h_term;
            if (isnan(pvectsz[ptsz->index_integrand])||isinf(pvectsz[ptsz->index_integrand])){
            printf("nan or inf in integrand gallens gallens 1h\n");
            printf("%.5e %.5e %.5e\n",pvectsz[ptsz->index_hmf],lensing_profile_at_ell_2,lensing_profile_at_ell_2);
            exit(0);
            }

                  }
           else if (_gallens_gallens_2h_){


             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_halo_bias];


                }


           else if (_lensmag_lensmag_1h_){
             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
             double lensing_profile_at_ell_1 = pvectsz[ptsz->index_lensing_profile];
             double lensing_profile_at_ell_2 = lensing_profile_at_ell_1;

                 pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                                   *lensing_profile_at_ell_1
                                                   *lensing_profile_at_ell_2
                                                   *damping_1h_term;
                  }
           else if (_lensmag_lensmag_2h_){

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_halo_bias];

                }


           else if (_lens_lensmag_1h_){

             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
             double lensing_profile_at_ell_1 = pvectsz[ptsz->index_lensing_profile];
             double lensing_profile_at_ell_2 = lensing_profile_at_ell_1;

                 pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                                   *lensing_profile_at_ell_1
                                                   *lensing_profile_at_ell_2
                                                   *damping_1h_term;


                  }
           else if (_lens_lensmag_2h_){


             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_halo_bias];

                }


  else if (_ngal_ngal_1h_){

    // int index_l = (int) pvectsz[ptsz->index_multipole];
    // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];


    // evaluate_galaxy_profile_1h(m_delta_cib,r_delta_cib,c_delta_cib,pvecback,pvectsz,pba,ptsz);
    evaluate_galaxy_profile_ngal(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);

    double ngal_profile_at_ell_1 = pvectsz[ptsz->index_galaxy_profile];

        pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                          *ngal_profile_at_ell_1
                                          *ngal_profile_at_ell_1
                                          *damping_1h_term;
   }

  else if (_cib_cib_1h_){

    int index_l = (int) pvectsz[ptsz->index_multipole];
    pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];


    evaluate_cib_profile(m_delta_cib,r_delta_cib,c_delta_cib,pvecback,pvectsz,pba,ptsz);
    double cib_profile_at_ell_1 = pvectsz[ptsz->index_cib_profile];

        pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                          *cib_profile_at_ell_1
                                          *cib_profile_at_ell_1
                                          *damping_1h_term;
   }

   else if (_cib_cib_2h_){
             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];
              evaluate_cib_profile(m_delta_cib,r_delta_cib,c_delta_cib,pvecback,pvectsz,pba,ptsz);


             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_cib_profile]
                                               *pvectsz[ptsz->index_halo_bias];

    }

  else if (_tSZ_cib_1h_){

    int index_l = (int) pvectsz[ptsz->index_multipole];
    pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];
    evaluate_cib_profile(m_delta_cib,r_delta_cib,c_delta_cib,pvecback,pvectsz,pba,ptsz);
    double cib_profile_at_ell_1 = pvectsz[ptsz->index_cib_profile];
    // pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
    evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
    double pressure_profile_at_ell_2 = pvectsz[ptsz->index_pressure_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *cib_profile_at_ell_1
                                      *pressure_profile_at_ell_2
                                      *damping_1h_term;
   }
   else if  (_tSZ_cib_2h_){


             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
             evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_pressure_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];
             evaluate_cib_profile(m_delta_cib,r_delta_cib,c_delta_cib,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_cib_profile]
                                               *pvectsz[ptsz->index_halo_bias];

                                             }

          }

  else if (_lens_cib_1h_){

    int index_l = (int) pvectsz[ptsz->index_multipole];
    pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];

     evaluate_cib_profile(m_delta_cib,r_delta_cib,c_delta_cib,pvecback,pvectsz,pba,ptsz);

    double cib_profile_at_ell_1 = pvectsz[ptsz->index_cib_profile];
    // pvectsz[ptsz->index_multipole_for_lensing_profile] =  ptsz->ell[index_l];
    evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
    double lensing_profile_at_ell_2 = pvectsz[ptsz->index_lensing_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *cib_profile_at_ell_1
                                      *lensing_profile_at_ell_2
                                      *damping_1h_term;
   }
   else if  (_lens_cib_2h_){


             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_lensing_profile] =  ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];
             evaluate_cib_profile(m_delta_cib,r_delta_cib,c_delta_cib,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_cib_profile]
                                               *pvectsz[ptsz->index_halo_bias];

                                             }

          }






  else if (_gal_cib_1h_){

    int index_l = (int) pvectsz[ptsz->index_multipole];
    pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];

    evaluate_cib_profile(m_delta_cib,r_delta_cib,c_delta_cib,pvecback,pvectsz,pba,ptsz);

    double cib_profile_at_ell_1 = pvectsz[ptsz->index_cib_profile];
    pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
    evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *cib_profile_at_ell_1
                                      *pvectsz[ptsz->index_galaxy_profile]
                                      *damping_1h_term;
   }
   else if  (_gal_cib_2h_){


             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

             int index_l = (int) pvectsz[ptsz->index_multipole];
            pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
            evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_galaxy_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];

            evaluate_cib_profile(m_delta_cib,r_delta_cib,c_delta_cib,pvecback,pvectsz,pba,ptsz);


             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_cib_profile]
                                               *pvectsz[ptsz->index_halo_bias];

                                             }

          }


  else if (_tSZ_gal_1h_){

    // int index_l = (int) pvectsz[ptsz->index_multipole];
    // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
    evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_ell_1 = pvectsz[ptsz->index_galaxy_profile];
    // pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
    evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);

    double pressure_profile_at_ell_2 = pvectsz[ptsz->index_pressure_profile];

        pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                          *galaxy_profile_at_ell_1
                                          *pressure_profile_at_ell_2
                                          *damping_1h_term;


   }

   else if  (_tSZ_gal_2h_){


             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
             evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_dlnMdeltadlnM]
                                               *pvectsz[ptsz->index_completeness]
                                               *pvectsz[ptsz->index_pressure_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
             evaluate_galaxy_profile_2h(kl,m_delta_gal,r_delta_gal,c_delta_gal,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_dlnMdeltadlnM]
                                               *pvectsz[ptsz->index_completeness]
                                               *pvectsz[ptsz->index_galaxy_profile]
                                               *pvectsz[ptsz->index_halo_bias];

                                             }

          }
  else if (_lens_lens_1h_){

    // int index_l = (int) pvectsz[ptsz->index_multipole];
    // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
    evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
    double lensing_profile_at_ell_1 = pvectsz[ptsz->index_lensing_profile];

        pvectsz[ptsz->index_integrand] =
                                          pvectsz[ptsz->index_hmf]
                                          *pvectsz[ptsz->index_dlnMdeltadlnM]
                                          *pvectsz[ptsz->index_completeness]
                                          *lensing_profile_at_ell_1
                                          *lensing_profile_at_ell_1
                                          *damping_1h_term;
   }

   else if (_lens_lens_2h_){

     // int index_l = (int) pvectsz[ptsz->index_multipole];
     // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
     evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
     double lensing_profile_at_ell_1 = pvectsz[ptsz->index_lensing_profile];

     evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);

         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias]
                                           *lensing_profile_at_ell_1;

    }


  else if (_tSZ_lens_1h_){

    // int index_l = (int) pvectsz[ptsz->index_multipole];
    // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
    evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);
    double lensing_profile_at_ell_1 = pvectsz[ptsz->index_lensing_profile];
    // pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
    evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
    double pressure_profile_at_ell_2 = pvectsz[ptsz->index_pressure_profile];

        pvectsz[ptsz->index_integrand] =
                                          pvectsz[ptsz->index_hmf]
                                          *pvectsz[ptsz->index_dlnMdeltadlnM]
                                          *pvectsz[ptsz->index_completeness]
                                          *lensing_profile_at_ell_1
                                          *pressure_profile_at_ell_2
                                          *damping_1h_term;
   }

   else if  (_tSZ_lens_2h_){


             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
             evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_dlnMdeltadlnM]
                                               *pvectsz[ptsz->index_completeness]
                                               *pvectsz[ptsz->index_pressure_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
             // int index_l = (int) pvectsz[ptsz->index_multipole];
             // pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
             evaluate_lensing_profile(kl,m_delta_lensing,r_delta_lensing,c_delta_lensing,pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_dlnMdeltadlnM]
                                               *pvectsz[ptsz->index_completeness]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_halo_bias];

                                             }

          }



   else if (_isw_tsz_){
     // pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];
     evaluate_pressure_profile(kl,pvecback,pvectsz,pba,ptsz);
     double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];
     evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ppt,ptsz);
     pvectsz[ptsz->index_integrand] =
                                       pvectsz[ptsz->index_hmf]
                                       *pvecback[pba->index_bg_D]
                                       *pvectsz[ptsz->index_dlnMdeltadlnM]
                                       *pvectsz[ptsz->index_completeness]
                                       *pvectsz[ptsz->index_halo_bias]
                                       *pressure_profile_at_ell;


   }
   return pvectsz[ptsz->index_integrand];
}


double get_te_of_m500c_at_z_arnaud(double m, double z, struct background * pba,struct tszspectrum * ptsz){
  double tau;
  int first_index_back = 0;
  double * pvecback;
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              ptsz->error_message);

  class_call(background_tau_of_z(pba,z,&tau),
             ptsz->error_message,
             ptsz->error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &first_index_back,
                               pvecback),
             ptsz->error_message,
             ptsz->error_message);

double Eh = pvecback[pba->index_bg_H]/pba->H0;

free(pvecback);

double mass = m;
return  5.*pow(Eh*mass/3.e14,2./3.); //kB*Te in keV

}

double get_te_of_m500c_at_z_lee(double m, double z, struct background * pba,struct tszspectrum * ptsz){
  double tau;
  int first_index_back = 0;
  double * pvecback;
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              ptsz->error_message);

  class_call(background_tau_of_z(pba,z,&tau),
             ptsz->error_message,
             ptsz->error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &first_index_back,
                               pvecback),
             ptsz->error_message,
             ptsz->error_message);

double Eh = pvecback[pba->index_bg_H]/pba->H0;

free(pvecback);

      double A[3] = {4.763,4.353,3.997};
      double B[3] = {0.581,0.571,0.593};
      double C[3] = {0.013,0.008,0.009};
      double z_interp[3] = {0.0,0.5,1.0};

      double Ap,Bp,Cp,zp;

      zp = z;
      if (zp <= 1.){
      Ap = pwl_value_1d(3,z_interp,A,zp);
      Bp = pwl_value_1d(3,z_interp,B,zp);
      Cp = pwl_value_1d(3,z_interp,C,zp);}

      else {
        Ap = A[2];
        Bp = B[2];
        Cp = C[2];
      }



   double Mfid = 3.e14 ; //Msun/h

   return Ap*pow(Eh,2./3.)*pow(m/Mfid,Bp+Cp*log(m/Mfid)); //kB*Te in keV

}


int evaluate_temperature_mass_relation( double * pvecback,
                                        double * pvectsz,
                                        struct background * pba,
                                        struct tszspectrum * ptsz)
  {


double mass = pvectsz[ptsz->index_m500]; //biased mass = M/B (X-ray mass)

if (ptsz->temperature_mass_relation == 0){
   // pvectsz[ptsz->index_te_of_m] = 5.*pow(Eh*mass //biased mass
   //                                       /3.e14,2./3.); //kB*Te in keV
pvectsz[ptsz->index_te_of_m] = get_te_of_m500c_at_z_arnaud(mass,pvectsz[ptsz->index_z],pba,ptsz);
      }
   //lee et al [1912.07924v1]
else if (ptsz->temperature_mass_relation == 1){


pvectsz[ptsz->index_te_of_m] = get_te_of_m500c_at_z_lee(mass*ptsz->HSEbias,pvectsz[ptsz->index_z],pba,ptsz);

   }


   return _SUCCESS_;
}


int evaluate_tau_profile(
                        double k,
                        double * pvecback,
                        double * pvectsz,
                        struct background * pba,
                        struct tszspectrum * ptsz)
{

  double m_asked;
   m_asked = pvectsz[ptsz->index_mass_for_electron_density]; // in Msun/h
   // m_asked = pvectsz[ptsz->index_m200m]; // in Msun/h


   double result;
   double tau_normalisation;
   double rho0;


   // double l_asked = pvectsz[ptsz->index_multipole_for_tau_profile];
   double chi = sqrt(pvectsz[ptsz->index_chi2]);
   // double k = (l_asked+0.5)/chi;
   double z_asked = pvectsz[ptsz->index_z];
   // printf("l_asked %.8e")
   result = get_gas_density_profile_at_k_M_z(k,m_asked,z_asked,ptsz);
   // printf("%.5e %.5e %.5e %.5e %.5e\n",
   // z_asked,
   // m_asked,
   // get_gas_density_profile_at_k_M_z(1e-3,m_asked,z_asked,ptsz)/m_asked,
   // get_gas_density_profile_at_k_M_z(1e-2,m_asked,z_asked,ptsz)/m_asked,
   // ptsz->f_b_gas);



   // double m_delta = m_asked;
   // double xout = 1;
   // double r_delta = pow(3.*m_delta/(4.*_PI_*200.*pvecback[pba->index_bg_Omega_m]*pvectsz[ptsz->index_Rho_crit]),1./3.);
   // double c_delta = get_c200m_at_m_and_z(m_delta,z_asked,pba,ptsz);
   // double result_ex =  evaluate_truncated_nfw_profile(k,r_delta,c_delta,xout,pvectsz,pba,ptsz);

   // double f_b = pba->Omega0_b/ptsz->Omega_m_0;
   // result_ex *= f_b*m_delta*pow((pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0,-1);

   // printf("result = %.8e result_ex = %.8e\n",result,result_ex);

   // pvectsz[ptsz->index_tau_profile] = 1./ptsz->mu_e*ptsz->f_free*result;
   // rho0 = 1.;
  double sigmaT_over_mp = 8.305907197761162e-17 * pow(pba->h,2)/pba->h; // !this is sigmaT / m_prot in (Mpc/h)**2/(Msun/h)
  double z = pvectsz[ptsz->index_z];
  double a = 1. / (1. + z);
  double tau_fac = a*sigmaT_over_mp/ptsz->mu_e*ptsz->f_free;

  // pvectsz[ptsz->index_tau_profile] = tau_fac*result*pow(1.+z,3.)/pow(chi,2.);
  pvectsz[ptsz->index_tau_profile] = tau_fac*result*pow(1.+z,3.)/pow(chi,2.);

                                    //  //*pow(pvecback[pba->index_bg_ang_distance]*pba->h,-2.); //(rs*ls)^2 in [Mpc/h]^2


   return _SUCCESS_;


  }

double get_ksz_filter_at_l(double ell,
                           struct tszspectrum * ptsz){
double fl = 1.;
if  (ell <= ptsz->l_unwise_filter[0] || ell >= ptsz->l_unwise_filter[ptsz->unwise_filter_size-1])
  fl = 0.;
else
  fl = pwl_value_1d(ptsz->unwise_filter_size,
                    ptsz->l_unwise_filter,
                    ptsz->f_unwise_filter,
                    ell);
return fl;
                           }

double get_M_min_of_z(double ell,
                      struct tszspectrum * ptsz){
double fl;
if  (ell <= ptsz->M_min_of_z_z[0])
  fl = ptsz->M_min_of_z_M_min[0];
else if  (ell >= ptsz->M_min_of_z_z[ptsz->M_min_of_z_size-1])
  fl = ptsz->M_min_of_z_M_min[ptsz->M_min_of_z_size-1];
else
  fl = pwl_value_1d(ptsz->M_min_of_z_size,
                    ptsz->M_min_of_z_z,
                    ptsz->M_min_of_z_M_min,
                    ell);
return fl;
                           }


double get_tau_profile_at_z_m_l(double z,
                                double m,
                                double l,
                                struct tszspectrum * ptsz,
                                struct background * pba){

  double tau;
  int first_index_back = 0;
  double * pvecback;
  class_alloc(pvecback,
        pba->bg_size*sizeof(double),
        ptsz->error_message);

  class_call(background_tau_of_z(pba,z,&tau),
       ptsz->error_message,
       ptsz->error_message);

  class_call(background_at_tau(pba,
                         tau,
                         pba->long_info,
                         pba->inter_normal,
                         &first_index_back,
                         pvecback),
       ptsz->error_message,
       ptsz->error_message);


double chi = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi
// double k1 = (l+0.5)/chi;
// double l = chi*k-0.5;

double tau_e  = get_gas_density_profile_at_k_M_z(l,m,z,ptsz);
// printf("%.5e %.5e %.5e %.5e\n",z,m,get_gas_density_profile_at_k_M_z(1e-3,m,z,ptsz)/m,ptsz->f_b_gas);

tau_e *= 1./ptsz->mu_e*ptsz->f_free;
double sigmaT_over_mp = 8.305907197761162e-17 * pow(pba->h,2)/pba->h; // !this is sigmaT / m_prot in (Mpc/h)**2/(Msun/h)
// double z = pvectsz[ptsz->index_z];
double a = 1. / (1. + z);
double rho0 = 1.;
tau_e *= sigmaT_over_mp*a*rho0*pow(pvecback[pba->index_bg_ang_distance]*pba->h,-2.);

free(pvecback);
return tau_e;
}


int evaluate_matter_density_profile(double k,
                                    double r_delta,
                                    double c_delta,
                                    double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct tszspectrum * ptsz)
{

    // int index_k = (int) pvectsz[ptsz->index_k_for_pk_hm];
    // double k = ptsz->k_for_pk_hm[index_k];

    double result;
    double characteristic_radius;
    //double characteristic_multipole;
    double density_normalisation;

    double xout = ptsz->x_out_truncated_nfw_profile;
    double result_trunc =  evaluate_truncated_nfw_profile(pvectsz[ptsz->index_z],k,r_delta,c_delta,xout);

    pvectsz[ptsz->index_density_profile] = result_trunc;



    density_normalisation = 1.;

    double rho0 = 1.;


    pvectsz[ptsz->index_density_profile] = density_normalisation
                                         *rho0
                                         *pvectsz[ptsz->index_mass_for_matter_density]
                                         *pow((pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0,-1)
                                         *pvectsz[ptsz->index_density_profile]; //rs in Mpc/h

   return _SUCCESS_;
}



int evaluate_lensing_profile(double kl,
                             double m_delta,
                             double r_delta,
                             double c_delta,
                             double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct tszspectrum * ptsz)
{

  double mass_nfw = m_delta;

   double rs = r_delta/c_delta;
   // pvectsz[ptsz->index_rs] = rs;
   // pvectsz[ptsz->index_ls] = sqrt(pvectsz[ptsz->index_chi2])/(1.+pvectsz[ptsz->index_z])/pvectsz[ptsz->index_rs];



   double result;
   double characteristic_radius;
   double characteristic_multipole;
   double lensing_normalisation;


   // characteristic_radius = pvectsz[ptsz->index_rs]; // in Mpc/h
   // characteristic_multipole = pvectsz[ptsz->index_ls];

   // pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile] = characteristic_multipole;
   // pvectsz[ptsz->index_multipole_for_nfw_profile] = pvectsz[ptsz->index_multipole_for_lensing_profile];


   // pvectsz[ptsz->index_multipole_for_truncated_nfw_profile] = pvectsz[ptsz->index_multipole_for_lensing_profile];  // multipole going into truncated nfw, unfortunate name
   double xout = ptsz->x_out_truncated_nfw_profile;
   // double l = pvectsz[ptsz->index_multipole_for_truncated_nfw_profile];
   double chi = sqrt(pvectsz[ptsz->index_chi2]);
   double z = pvectsz[ptsz->index_z];
   // double k = (l+0.5)/chi;
   // double l = chi*kl-0.5;
   // if (l<0.)
   //  result = 1e-100;
   // else
    result =  evaluate_truncated_nfw_profile(pvectsz[ptsz->index_z],kl,r_delta,c_delta,xout);



   pvectsz[ptsz->index_lensing_profile] = result;

  // lensing normalisation is dimensionless
  // we write down things in terms of phi here
  // (\kappa = l(l+1)\phi_l/2, see Hill & Spergel 1312.4525)

  int index_md = (int) pvectsz[ptsz->index_md];
  // if ( _gal_lens_2h_
  //   || _gal_lens_1h_
  //   || _gal_lensmag_1h_
  //   || _gal_lensmag_2h_
  //   || _gal_gallens_1h_
  //   || _gal_gallens_2h_
  //   || _kSZ_kSZ_gallens_1h_fft_
  //   || _kSZ_kSZ_gallens_2h_fft_
  //   || _kSZ_kSZ_gallens_3h_fft_
  //   || _gallens_gallens_1h_
  //   || _gallens_gallens_2h_
  //   || _gallens_lens_1h_
  //   || _gallens_lens_2h_
  //   || _lensmag_lensmag_1h_
  //   || _lensmag_lensmag_2h_
  //   || _lens_lensmag_1h_
  //   || _lens_lensmag_2h_
  //   //|| _cib_lens_1h_
  //   //|| _cib_lens_1h_
  //   || _tSZ_lens_1h_
  //   || _tSZ_lens_2h_
  //   || _tSZ_lensmag_1h_
  //   || _tSZ_lensmag_2h_
  //   || _kSZ_kSZ_lensmag_1halo_
  //   || _lens_lens_1h_
  //   || _lens_lens_2h_)
  lensing_normalisation = 1.; // kappa
  // else // phi
  // lensing_normalisation = 2./(l*(l+1.));


  // double rho0;
  // rho0 = mass_nfw/(4.*_PI_*pow(rs,3.));


   pvectsz[ptsz->index_lensing_profile] =  lensing_normalisation // dim less
                                           *mass_nfw // M
                                           *pvectsz[ptsz->index_lensing_profile]
                                           /(pow(3.*pow(pba->H0/pba->h,2)/2./ptsz->Rho_crit_0,-1)*pow((1.+z),1.)/chi)
                                           ///pvectsz[ptsz->index_lensing_Sigma_crit] // Sigma_crit is in M/L^2
                                           *pow(pvecback[pba->index_bg_ang_distance]*pba->h,-2.); // 1/L^2



   return _SUCCESS_;
}

int evaluate_pressure_profile(double kl,
                              double * pvecback,
                              double * pvectsz,
                              struct background * pba,
                              struct tszspectrum * ptsz)
{
  double z = pvectsz[ptsz->index_z];
  double chi =  sqrt(pvectsz[ptsz->index_chi2]);
   //Fourier tranform
   //of the  pressure profile

   //If customized profile (!=P13 or A10)
   //perform the Fourier transform
   //at each ell, M, z...
   //time consuming.
   //(but still manageable for
   //MCMC runs)

   //Else, read the FT from tabulated
   //values and interpolate
   //at the desired ell/ell500.
   int index_md = (int) pvectsz[ptsz->index_md];
   double lnx_asked;

   double result = 0.;

   //Battaglia et al 2012
   if (ptsz->pressure_profile == 4 ){

         double m_asked = pvectsz[ptsz->index_m200c];
         double k_asked = kl*(1.+z)*pvectsz[ptsz->index_r200c]; // this is (l+1/2)/l200c
         // double m_asked = pvectsz[ptsz->index_m200m]; // in Msun/h
         double z_asked = pvectsz[ptsz->index_z];


         if( (log(k_asked)<ptsz->array_pressure_profile_ln_k[0]) || _mean_y_ || _dydz_){ // get large scale limit
           result = get_gas_pressure_profile_at_k_m_z(exp(ptsz->array_pressure_profile_ln_k[0]),m_asked,z_asked,ptsz);
         }
         // else if(log(k_asked)<ptsz->array_pressure_profile_ln_k[0]){
         //   // also large scale limit in this case
         //   result = get_gas_pressure_profile_at_k_m_z(exp(ptsz->array_pressure_profile_ln_k[0]),m_asked,z_asked,ptsz);
         // }
         else{
           double result_tabulated = get_gas_pressure_profile_at_k_m_z(k_asked,m_asked,z_asked,ptsz);
           result = result_tabulated;
         }
         // printf("%.7e %.7e %.7e\n",result,result_tabulated,m_asked);

   }
   //custom gNFW pressure profile at 500c
   else if(ptsz->pressure_profile == 3){

     // printf("ptsz->delta_def_electron_pressure = %d\n",ptsz->delta_def_electron_pressure);
          // double l_delta = 0.;
       if ( ptsz->delta_def_electron_pressure == 2){
         if (ptsz->mass_dependent_bias == 1)
            ptsz->HSEbias = 1./(0.8/(1.+ ptsz->Ap*pow(pvectsz[ptsz->index_m500]/3.e14,ptsz->alpha_b)));
        //m500 X-ray for the pressure profiles A10 and P13
         pvectsz[ptsz->index_m500] = pvectsz[ptsz->index_m500c]/ptsz->HSEbias;
         // printf("m500 = %.3e\n",pvectsz[ptsz->index_m500]);
         //r500 X-ray for the pressure profiles A10 and P13
         pvectsz[ptsz->index_r500] = pow(3.*pvectsz[ptsz->index_m500]/(4.*_PI_*500.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc
         // printf("r500 = %.3e\n",pvectsz[ptsz->index_r500]);
         pvectsz[ptsz->index_l500] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_r500];
         // l_delta = pvectsz[ptsz->index_l500];
       }
       else if (ptsz->delta_def_electron_pressure == 1){ // 200c
         pvectsz[ptsz->index_m500] = pvectsz[ptsz->index_m200c];
         // printf("m500 = %.3e\n",pvectsz[ptsz->index_m500]);
         //r500 X-ray for the pressure profiles A10 and P13
         pvectsz[ptsz->index_r500] = pow(3.*pvectsz[ptsz->index_m500]/(4.*_PI_*200.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc
         // printf("r500 = %.3e\n",pvectsz[ptsz->index_r500]);
         pvectsz[ptsz->index_l500] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_r500];
       }
       else{
         printf("This delta definition for electron pressure is not implemented yet.\n");
         exit(0);
       }

      // lnx_asked = log(kl*chi/pvectsz[ptsz->index_l500]);
      lnx_asked = log(kl*(1.+z)*pvectsz[ptsz->index_r500]);

      if(lnx_asked<ptsz->array_profile_ln_l_over_ls[0] || _mean_y_ || _dydz_) // large scale limit
         result = ptsz->array_profile_ln_PgNFW_at_lnl_over_ls[0];
      else if (lnx_asked>ptsz->array_profile_ln_l_over_ls[ptsz->array_profile_ln_PgNFW_at_lnl_over_ls_size-1])
         result = -100.;
      else
        result = pwl_value_1d(ptsz->array_profile_ln_PgNFW_at_lnl_over_ls_size,
                              ptsz->array_profile_ln_l_over_ls,
                              ptsz->array_profile_ln_PgNFW_at_lnl_over_ls,
                              lnx_asked);

      result = exp(result);
      // printf("lnx_asked = %.3e pp=%.3e\n",lnx_asked ,result);

   }
  // tabulated A10 and P13 profiles
  else {

   if (ptsz->mass_dependent_bias == 1)
      ptsz->HSEbias = 1./(0.8/(1.+ ptsz->Ap*pow(pvectsz[ptsz->index_m500]/3.e14,ptsz->alpha_b)));


  //m500 X-ray for the pressure profiles A10 and P13
   pvectsz[ptsz->index_m500] = pvectsz[ptsz->index_m500c]/ptsz->HSEbias;
   // printf("m500 = %.3e\n",pvectsz[ptsz->index_m500]);

   //r500 X-ray for the pressure profiles A10 and P13
   pvectsz[ptsz->index_r500] = pow(3.*pvectsz[ptsz->index_m500]/(4.*_PI_*500.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc
   // printf("r500 = %.3e\n",pvectsz[ptsz->index_r500]);

   pvectsz[ptsz->index_l500] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_r500];


      // lnx_asked = log(kl*chi/pvectsz[ptsz->index_l500]);
      lnx_asked = log(kl*(1.+z)*pvectsz[ptsz->index_r500]);
      if(lnx_asked<ptsz->PP_lnx[0] || _mean_y_ || _dydz_)
         result = ptsz->PP_lnI[0];
      else if (lnx_asked>ptsz->PP_lnx[ptsz->PP_lnx_size-1])
         result = -100.;
      else splint(ptsz->PP_lnx,
                  ptsz->PP_lnI,
                  ptsz->PP_d2lnI,
                  ptsz->PP_lnx_size,
                  lnx_asked,
                  &result);

      result = exp(result);
   }



   pvectsz[ptsz->index_pressure_profile] = result;

    // in units of Mpc^-1*micro Kelvins
    // old version (szfast)
    // double sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb = 283./0.5176; //1./0.5176=1.932=(5Xh+3)/2(Xh+1) with Xh = 0.76 and Pth=1.932Pe
    // (Xh is the primodial hydrogen mass fraction)
    // more accurate version (see explanation below):
    // in units of Mpc^-1*micro Kelvins
    double sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb = 283.2980000259841/0.5176*pba->T_cmb/2.725;

    // Explanation of the above factors:
    // sigma_thomson_in_m2 = 6.6524587321e-29
    // me_in_eV = 0.510998946e6
    // Mpc_over_m =  3.085677581282e22
    // factor_in_inverse_Mpc =  sigma_thomson_in_m2/me_in_eV*50.*1e6*Mpc_over_m*0.5176
    // print(factor_in_inverse_Mpc*2.725e6)
    // 283.2980000259841
    // we divide by 0.5176 because it was included in the 283. It is necessary when working with
    // thermal pressure, but not when we work with electron pressure:  Pth=1.932*Pe, Pe = 0.5176*Pth


   double characteristic_radius;
   double characteristic_multipole;
   double pressure_normalisation;


   //Battaglia et al 2012 pressure profile
   if (ptsz->pressure_profile == 4) {

      double P_200; //in units of eV/cm^3, corresponds to Pth


      double R_200crit = pvectsz[ptsz->index_r200c]; //in units of h^-1 Mpc
      double f_b = ptsz->f_b_gas;//pba->Omega0_b/ptsz->Omega_m_0;

      double Eh = pvecback[pba->index_bg_H]/pba->H0;

      //double rho_crit_at_z = pvectsz[ptsz->index_Rho_crit]; //in units of h^2 M_sun/Mpc^3
      //double _G_in_eV_Mpc_over_Msun2 = _G_/(_eV_ *_Mpc_over_m_ /_M_sun_/_M_sun_);
      //double P_200_boris = _G_in_eV_Mpc_over_Msun2*pvectsz[ptsz->index_m200c]
      //            /pba->h*200.*rho_crit_at_z*pba->h*pba->h
      //            *f_b/2./(R_200crit/pba->h)/pow(_Mpc_over_m_,3.)/1e6;

      P_200 = pvectsz[ptsz->index_m200c]/(R_200crit)*f_b
              *2.61051e-18*pow(100.*pba->h*Eh,2.);


     // NB JCH implementation:
     // transform2d=4.0d0*pi*(r500/h0)/(l500**2d0)* &
     //      2.61051d-18*(obh2/om0/h0**2d0)*(100.0d0*h0*Ez(z))**2d0* &
     //      m500/r500*2.5d0*rombint3(padia,xin,xoutpress,tol,(ell+0.5d0)/l500,m500,z)
     // ! the 10.94d0 factor below is sigmaT/m_e/c^2*Mpcincm*Tcmb*10^6
     // ! in units that agree with the eV/cm^3 in transform2d
     // Tsz=-2d0*10.94d0*transform2d ! uK, Rayleigh-Jeans limit

     // link with class_sz: (283./0.5176/50.)=10.935085007727976

     // (2.61....m500/r500) is the pressure normalization for the Battaglia pressure profile in eV/cm^3




      pressure_normalisation = P_200;///1.932; //1./0.5176=1.932=(5Xh+3)/2(Xh+1) with Xh = 0.76
      //divide by 1.932 to obtain Pe
      characteristic_radius = pvectsz[ptsz->index_r200c]/pba->h; //in Mpc
      characteristic_multipole = pvectsz[ptsz->index_l200c];



   }

   else {

     // formula D1 of WMAP 7 year paper (normalisation of pressure profile)
     // see also formula 5 of Planck 2013 profile paper, P500:
      double C_pressure = 1.65*pow(pba->h/0.7,2)
                          *pow(pvecback[pba->index_bg_H]/pba->H0,8./3.)
                          *pow(pvectsz[ptsz->index_m500]/(3.e14*0.7),2./3.+ptsz->alpha_p)
                          *pow(pvectsz[ptsz->index_m500]/3.e14, ptsz->delta_alpha);



      //A10:
      if (ptsz->pressure_profile == 2)
      pressure_normalisation = C_pressure
                               *ptsz->P0GNFW
                               *pow(0.7/pba->h, 1.5); // as found by dimensional analysis (X-ray data, see email with E. Komatsu and R. Makya)
      //P13:
      else if (ptsz->pressure_profile == 0)
      pressure_normalisation = C_pressure
                               *ptsz->P0GNFW
                               *pow(0.7/pba->h, 1.); // as found by dimensional analysis (sz data, see email with E. Komatsu and R. Makya)

      //Custom. GNFW
      else if (ptsz->pressure_profile == 3)
      pressure_normalisation = C_pressure
                               *ptsz->P0GNFW
                               *pow(0.7/pba->h, 1.5); // assuming X-ray data based pressure profile
     // //Custom. GNFW
     // else if (ptsz->pressure_profile == 3)
     // pressure_normalisation = C_pressure
     //                          *ptsz->P0GNFW
     //                          *pow(0.7/pba->h, 1.); // assuming SZ data based pressure profile


      characteristic_radius = pvectsz[ptsz->index_r500]/pba->h; // in Mpc
      characteristic_multipole = pvectsz[ptsz->index_l500];


   }


   //(see comments after to link this way of writing with the szfast implementation)
   pvectsz[ptsz->index_pressure_profile] = sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb // here Tcmb is in muK
                                           /50. // to cancel the factor 50 above 50eV/cm^3
                                           /pba->T_cmb
                                           *pressure_normalisation
                                           *pvectsz[ptsz->index_pressure_profile]
                                           *(4*_PI_)
                                           *pow(characteristic_multipole,-2)
                                           *characteristic_radius //rs in Mpc
                                           *ptsz->Tcmb_gNU;
                                           // now in units Tcmb*gNU at the requested frequency
                                           // Result in muK due to the units in sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb
                                           // which is in units of Mpc^-1*micro Kelvins
                                           // (but Tcmb is in K in Tcmb_gNU)

   // Ancient way of writing (like in szfast):
   // pvectsz[ptsz->index_pressure_profile] = (ptsz->Tcmb_gNU_at_150GHz/pba->T_cmb)      //gnu at 150 GHz (was -0.953652 in sz_fast), we replace it with the exact formula
   //                                         *sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb
   //                                         *pressure_normalisation
   //                                         *pvectsz[ptsz->index_pressure_profile]
   //                                         *(4*_PI_)
   //                                         *pow(characteristic_multipole,-2)
   //                                         *characteristic_radius //rs in Mpc
   //                                         /50. //pressure normalised to 50eV/cm^3
   //                                         *ptsz->Tcmb_gNU/ptsz->Tcmb_gNU_at_150GHz;
   //                                         // now in units Tcmb*gNU at the requested frequency
   //                                         // Result in muK due to the units in sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb
   //                                         // which is in units of Mpc^-1*micro Kelvins
   //                                         // (but Tcmb is in K in Tcmb_gNU)

   return _SUCCESS_;
}







// normaized pressure profile battaglia et al 2012
double get_pressure_P_over_P_delta_at_x_M_z_b12_200c(double x_asked,
                                                     double m_asked,
                                                     double z_asked,
                                                     double c_asked,
                                                     double A_P0,
                                                     double A_xc,
                                                     double A_beta,
                                                     double alpha_m_P0,
                                                     double alpha_m_xc,
                                                     double alpha_m_beta,
                                                     double alpha_z_P0,
                                                     double alpha_z_xc,
                                                     double alpha_z_beta,
                                                     // break model
                                  							     double mcut,
                                  							     double alphap_m_P0,
                                  							     double alphap_m_xc,
                                  							     double alphap_m_beta,
                                  							     double alpha_c_P0,
                                  							     double alpha_c_xc,
                                  							     double alpha_c_beta,
                                                     // end break model
                                                     double alpha,
                                                     double gamma,
                                                     struct background * pba,
                                                     struct tszspectrum * ptsz){

  /// NORMALIZED PRESSURE PROFILE PART

  double xc;
  double beta;
  double P0;

  double m200_over_msol = m_asked/pba->h; // convert to Msun
  double x = x_asked;
  double z = z_asked;

  if (m200_over_msol > mcut) {

  // P0 = A_P0*pow(m200_over_msol/1e14,alpha_m_P0)*pow(1.+z,alpha_z_P0);
  // xc = A_xc*pow(m200_over_msol/1e14,alpha_m_xc)*pow(1.+z,alpha_z_xc);
  // beta = A_beta*pow(m200_over_msol/1e14,alpha_m_beta)*pow(1.+z,alpha_z_beta);

    P0 = A_P0*pow(m200_over_msol/mcut,alpha_m_P0)*pow(1.+z,alpha_z_P0)*pow(1.+c_asked,alpha_c_P0);
    xc = A_xc*pow(m200_over_msol/mcut,alpha_m_xc)*pow(1.+z,alpha_z_xc)*pow(1.+c_asked,alpha_c_xc);
    beta = A_beta*pow(m200_over_msol/mcut,alpha_m_beta)*pow(1.+z,alpha_z_beta)*pow(1.+c_asked,alpha_c_beta);

}
else {
  P0 = A_P0*pow(m200_over_msol/mcut,alphap_m_P0)*pow(1.+z,alpha_z_P0)*pow(1.+c_asked,alpha_c_P0);
  xc = A_xc*pow(m200_over_msol/mcut,alphap_m_xc)*pow(1.+z,alpha_z_xc)*pow(1.+c_asked,alpha_c_xc);
  beta = A_beta*pow(m200_over_msol/mcut,alphap_m_beta)*pow(1.+z,alpha_z_beta)*pow(1.+c_asked,alpha_c_beta);
}
  // double gamma = -0.3;
  // double alpha = 1.0;

  double plc_x = P0*pow(x/xc,gamma)*pow(1.+ pow(x/xc,alpha),-beta);


  //putting everything together
  double result =  plc_x;

  return result;

}


// this is r_200c*P_200c
double get_1e6xdy_from_battaglia_pressure_at_x_z_and_m200c(double x,
                                                           double z,
                                                           double m,
                                                           struct background * pba,
                                                           struct tszspectrum * ptsz){

  double tau;
  int first_index_back = 0;
  double * pvecback;
  class_alloc(pvecback,
        pba->bg_size*sizeof(double),
        ptsz->error_message);

  class_call(background_tau_of_z(pba,z,&tau),
       ptsz->error_message,
       ptsz->error_message);

  class_call(background_at_tau(pba,
                         tau,
                         pba->long_info,
                         pba->inter_normal,
                         &first_index_back,
                         pvecback),
       ptsz->error_message,
       ptsz->error_message);

  double sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb = 283.2980000259841/0.5176*pba->T_cmb/2.725;


  double rho_crit = (3./(8.*_PI_*_G_*_M_sun_))
                  *pow(_Mpc_over_m_,1)
                  *pow(_c_,2)
                  *pvecback[pba->index_bg_rho_crit]
                  /pow(pba->h,2);

  double r200c =  pow(3.*m/(4.*_PI_*200.*rho_crit),1./3.);


  double f_b = ptsz->f_b_gas;//pba->Omega0_b/ptsz->Omega_m_0;
  double Eh = pvecback[pba->index_bg_H]/pba->H0;
  double d_A = pvecback[pba->index_bg_ang_distance]*pba->h; // in Mpc/h


  double P_200 = m/r200c*f_b*2.61051e-18*pow(100.*pba->h*Eh,2.);


  /// NORMALIZED PRESSURE PROFILE PART

  // double xc;
  // double beta;
  // double P0;
  //
  double m200_over_msol = m/pba->h; // convert to Msun
  //
  //
  // P0 = ptsz->P0_B12*pow(m200_over_msol/1e14,ptsz->alpha_m_P0_B12)*pow(1+z,ptsz->alpha_z_P0_B12);
  // xc = ptsz->xc_B12*pow(m200_over_msol/1e14,ptsz->alpha_m_xc_B12)*pow(1+z,ptsz->alpha_z_xc_B12);
  // beta = ptsz->beta_B12*pow(m200_over_msol/1e14,ptsz->alpha_m_beta_B12)*pow(1+z,ptsz->alpha_z_beta_B12);
  //
  // double gamma = ptsz->gamma_B12;
  // double alpha = ptsz->alpha_B12;
  //
  // double plc_x = P0*pow(x/xc,gamma)*pow(1.+ pow(x/xc,alpha),-beta);

double c_asked = ptsz->c_B12;
double Px = get_pressure_P_over_P_delta_at_x_M_z_b12_200c(x,m200_over_msol,z,
                                              c_asked,ptsz->P0_B12,
                                              ptsz->xc_B12,ptsz->beta_B12,
                                              ptsz->alpha_m_P0_B12,ptsz->alpha_m_xc_B12,
                                              ptsz->alpha_m_beta_B12,ptsz->alpha_z_P0_B12,
                                              ptsz->alpha_z_xc_B12,ptsz->alpha_z_beta_B12,
                                              // break model
                                  						ptsz->mcut_B12,ptsz->alphap_m_P0_B12,
                                  						ptsz->alphap_m_xc_B12,ptsz->alphap_m_beta_B12,
                                  						ptsz->alpha_c_P0_B12,
                                  						ptsz->alpha_c_xc_B12,
                                  						ptsz->alpha_c_beta_B12,
                                                     // end break model
                                              ptsz->alpha_B12,
                                              ptsz->gamma_B12,
                                              pba,ptsz);
double plc_x = Px;
  //putting everything together
  double result = sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb // here Tcmb is in muK
                   /50. // to cancel the factor 50 above 50eV/cm^3
                   /pba->T_cmb
                   *P_200
                   *plc_x
                   *r200c/pba->h; //in Mpc


  free(pvecback);
  return result;

}

// this is r_500c*P_500c
double get_1e6xdy_from_gnfw_pressure_at_x_z_and_m500c(double x,
                                                      double z,
                                                      double m,
                                                      double delta,
                                                      struct background * pba,
                                                      struct tszspectrum * ptsz){

  double tau;
  int first_index_back = 0;
  double * pvecback;
  class_alloc(pvecback,
        pba->bg_size*sizeof(double),
        ptsz->error_message);

  class_call(background_tau_of_z(pba,z,&tau),
       ptsz->error_message,
       ptsz->error_message);

  class_call(background_at_tau(pba,
                         tau,
                         pba->long_info,
                         pba->inter_normal,
                         &first_index_back,
                         pvecback),
       ptsz->error_message,
       ptsz->error_message);

  double sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb = 283.2980000259841/0.5176*pba->T_cmb/2.725;


  double rho_crit = (3./(8.*_PI_*_G_*_M_sun_))
                  *pow(_Mpc_over_m_,1)
                  *pow(_c_,2)
                  *pvecback[pba->index_bg_rho_crit]
                  /pow(pba->h,2);

  double r500c =  pow(3.*m/(4.*_PI_*delta*rho_crit),1./3.);

  double Eh = pvecback[pba->index_bg_H]/pba->H0;

  // formula D1 of WMAP 7 year paper (normalisation of pressure profile)
  // see also formula 5 of Planck 2013 profile paper, P500:
   // double C_pressure = 1.65*pow(pba->h/0.7,2)
   //                     *pow(Eh,8./3.)
   //                     *pow(m/(3.e14*0.7),2./3.+ptsz->alpha_p);
                       //*pow(m/3.e14, ptsz->delta_alpha);

  // hasselfield et al 2013:
   double C_pressure = 1.65*pow(pba->h/0.7,2)
                       *pow(Eh,8./3.)
                       *pow(m/(3.e14*0.7),2./3.+ptsz->alpha_p);
                       // *pow(m/(3.e14*0.7),0.22/(1.+8.*pow(x,3.))); hasselfield et al has this part!


   //A10:
  double pressure_normalisation = C_pressure
                                  *ptsz->P0GNFW
                                  *pow(0.7/pba->h, 1.5); // as found by dimensional analysis (X-ray data, see email with E. Komatsu and R. Makya)




  //putting everything together



  double plc_x =  (1./(pow(ptsz->c500*x,ptsz->gammaGNFW)
                  *pow(1.+ pow(ptsz->c500*x,ptsz->alphaGNFW),
                  (ptsz->betaGNFW-ptsz->gammaGNFW)/ptsz->alphaGNFW)));


  // double result = plc_x*ptsz->P0GNFW*pow(0.7/pba->h, 1.5);


  double result = sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb // here Tcmb is in muK
                  /50. // to cancel the factor 50 above 50eV/cm^3
                  /pba->T_cmb
                  *pressure_normalisation
                  *plc_x
                  *r500c/pba->h; //in Mpc
  free(pvecback);
  return result;

}


// normaized pressure profile gnfw Arnaud et al 2010
double get_pressure_P_over_P_delta_at_x_gnfw_500c(double x_asked,
                                                  double P0GNFW,
                                                  double alphaGNFW,
                                                  double betaGNFW,
                                                  double gammaGNFW,
                                                  double c500,
                                                  struct background * pba,
                                                  struct tszspectrum * ptsz){

  //  //A10:
  // double pressure_normalisation = C_pressure
  //                                 *ptsz->P0GNFW
  //                                 *pow(0.7/pba->h, 1.5); // as found by dimensional analysis (X-ray data, see email with E. Komatsu and R. Makya)
  //



  //putting everything together
  double x = x_asked;


  double plc_x =  (1./(pow(c500*x,gammaGNFW)
                  *pow(1.+ pow(c500*x,alphaGNFW),
                  (betaGNFW-gammaGNFW)/alphaGNFW)));


  // double result = plc_x*ptsz->P0GNFW*pow(0.7/pba->h, 1.5);


  double result = P0GNFW*plc_x; //in Mpc
  return result;

}



int evaluate_completeness(double * pvecback,
                          double * pvectsz,
                          struct background * pba,
                          struct tszspectrum * ptsz){

    double comp_at_M_and_z = 0.;

printf("check units of dA!!.\n");
exit(0);

    if (ptsz->has_completeness_for_ps_SZ == 1){

    comp_at_M_and_z = 0.;

    double mp_bias = pvectsz[ptsz->index_m500c]/ptsz->HSEbias; //biased mass = M/B
    //double redshift = pvectsz[ptsz->index_z];

    //printf("mass m500 = %e\n",mp_bias);
    //printf("mass m200 = %e\n",pvectsz[ptsz->index_m200]); //true mass
    //printf("bias = %e\n",ptsz->HSEbias);
    //printf("redshift = %e\n",redshift);
    double Eh = pvecback[pba->index_bg_H]/pba->H0;
    double d_A = pvecback[pba->index_bg_ang_distance]; //units Mpc

    //! szcounts.f90: angular diameter distance in units of h^-1 Mpc


    double H0 = pba->h*100.;


    double thetastar = 6.997;
    double alpha_theta = 1./3.;
    double alpha = 1.78;//pcsz->alpha; //1.78;
    double beta = 0.66;

    double thetastar2 = thetastar * pow(H0/70.,-2./3.);
    double theta500_for_mp_at_zp =  thetastar2 * pow(mp_bias/3.e14* (100./H0),alpha_theta);
    theta500_for_mp_at_zp *=    pow(Eh,-2./3) *pow(100.*d_A/(500.0*H0),-1.);
    double thp = theta500_for_mp_at_zp;
    double ystar = pow(10.,-0.19)/pow(2., alpha)*0.00472724;
    double ystar2 = ystar;
    ystar2 *=  pow(H0/70.,-2.+alpha);
    double y500_for_mp_at_zp =  ystar2 * pow(mp_bias/3.e14* (100./H0),alpha);
    y500_for_mp_at_zp *=   pow(Eh,beta) *pow(100.*d_A/(500.0*H0),-2.);
    double yp = y500_for_mp_at_zp;
    double sn_cutoff = ptsz->sn_cutoff;
    //printf("with cutoff = %.3e\n",sn_cutoff); // BB debug
    double y;
    int l1,l2;
    double th1,th2;


  if (thp > ptsz->theta_bin_max){
          l1 = ptsz->nthetas - 1;
          l2 = ptsz->nthetas - 2;
          th1 = ptsz->thetas[l1];
          th2 = ptsz->thetas[l2];
       }

    else if (thp < ptsz->theta_bin_min){
       l1 = 0;
       l2 = 1;
       th1 = ptsz->thetas[l1];
       th2 = ptsz->thetas[l2];
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

    // //Sum over sky pathces:
    // int index_patches;
    // double sum_skyfracs = 0.;
    //
    // for (index_patches =0;
    // index_patches<ptsz->nskyfracs;
    // index_patches++){
    //
    // double y1 = ptsz->ylims[index_patches][l1];
    // double y2 = ptsz->ylims[index_patches][l2];
    // y = y1 + (y2-y1)/(th2-th1)*(thp-th1);
    //
    // double c2 = erf_compl_ps(yp,y,sn_cutoff);
    //
    // comp_at_M_and_z += c2*ptsz->skyfracs[index_patches];
    // sum_skyfracs += ptsz->skyfracs[index_patches];
    // }
    // //Now divide by sky fraction
    // comp_at_M_and_z = comp_at_M_and_z/sum_skyfracs;


    // Using skyaveraged noisemap:
    double y1 = ptsz->sky_averaged_ylims[l1];
    double y2 = ptsz->sky_averaged_ylims[l2];
    y = y1 + (y2-y1)/(th2-th1)*(thp-th1);
    double c2 = erf_compl_ps(yp,y,sn_cutoff);
    comp_at_M_and_z =  c2;

// printf("completensess = %.3e\n",c2); // BB debug
//
// pvectsz[ptsz->index_completeness] = c2; // BB debug
//    return _SUCCESS_; // BB debug

    }


   if (ptsz->which_ps_sz == 0)
      pvectsz[ptsz->index_completeness] = 1.;
   else if (ptsz->which_ps_sz == 1) {//ps resolved
     //printf("completensess = %.3e\n",pvectsz[ptsz->index_completeness]); // BB debug

      pvectsz[ptsz->index_completeness] = comp_at_M_and_z;
       }
   else if (ptsz->which_ps_sz == 2){ //ps unresolved

      pvectsz[ptsz->index_completeness] = (1.-comp_at_M_and_z);
      //printf("completensess = %.3e\n",pvectsz[ptsz->index_completeness]); // BB debug
        }


   return _SUCCESS_;
}




int evaluate_effective_galaxy_bias_ngal(
                                   int index_g,
                                   double * pvecback,
                                   double * pvectsz,
                                   struct background * pba,
                                   struct primordial * ppm,
                                   struct nonlinear * pnl,
                                   struct tszspectrum * ptsz)
{


  double b;
  double z = pvectsz[ptsz->index_z];
  int index_md = (int) pvectsz[ptsz->index_md];

    b = ptsz->effective_galaxy_bias_ngal[index_g];
    pvectsz[ptsz->index_halo_bias] = b;

if (ptsz->has_ng_in_bh){
  double bng = 0.;
  int index_l = (int) pvectsz[ptsz->index_multipole];
  double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi
  double kl = (ptsz->ell[index_l]+0.5)/d_A;
  double bh=b;
  bng = get_scale_dependent_bias_at_z_and_k(z,kl,bh,ptsz);
  pvectsz[ptsz->index_halo_bias] += bng;
}



   return _SUCCESS_;
}




int evaluate_effective_galaxy_bias(double * pvecback,
                                   double * pvectsz,
                                   struct background * pba,
                                   struct primordial * ppm,
                                   struct nonlinear * pnl,
                                   struct tszspectrum * ptsz)
{



  // For unWISE galaxies, we use an effective galaxy bias, Eq. 6.2 of KFW20
  // it is simply redshift dependent and depends on the color of the sample:
  // See Eq. 5.6 of KFW20 for the details
  double b;
  double z = pvectsz[ptsz->index_z];
  int index_md = (int) pvectsz[ptsz->index_md];

  // // blue sample
  // if (ptsz->unwise_galaxy_sample_id==3){
  //
  //   if (_gal_gal_hf_ || _gal_lensmag_hf_)
  //   b = 1.74;
  //   else if (_gal_lens_hf_)
  //   //b = 0.8+1.2*z;
  //   b = 1.56;
  //
  // }
  //
  // // green (two greens, depending on the min halo mass cut)
  // else if (ptsz->unwise_galaxy_sample_id==2 || ptsz->unwise_galaxy_sample_id==1){
  //
  //   if (_gal_gal_hf_ || _gal_lensmag_hf_)
  //   b = 2.44;
  //   else if (_gal_lens_hf_)
  //   //b = gsl_max(1.6*z*z, 1.);
  //   b = 2.23;
  // }
  //
  // // red
  // else if (ptsz->unwise_galaxy_sample_id==0){
  //
  //   if (_gal_gal_hf_ || _gal_lensmag_hf_)
  //   b = 3.47;
  //   else if (_gal_lens_hf_)
  //   //b = gsl_max(2.*pow(z,1.5),1.);
  //   b = 3.29;
  // }
  //
  // else{
    b = ptsz->effective_galaxy_bias;
  // }



   pvectsz[ptsz->index_halo_bias] = b;

if (ptsz->has_ng_in_bh){
  double bng = 0.;
  int index_l = (int) pvectsz[ptsz->index_multipole];
  double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi
  double kl = (ptsz->ell[index_l]+0.5)/d_A;
  double bh=b;
  bng = get_scale_dependent_bias_at_z_and_k(z,kl,bh,ptsz);
  pvectsz[ptsz->index_halo_bias] += bng;
}



   return _SUCCESS_;
}


double get_scale_dependent_bias_at_z_and_k(double z_asked,
                                           double k_asked,
                                           double bh,
                                           struct tszspectrum *ptsz){
  double z = log(1.+z_asked);
  double k = log(k_asked);

  double result;

 if (z<ptsz->array_ln_1pz_ng_bias[0]){
// z = ptsz->array_ln_1pz_ng_bias[0];
result = 0.;
 }
else if (z>ptsz->array_ln_1pz_ng_bias[ptsz->nz_ng_bias-1]){
    // z = ptsz->array_ln_1pz_ng_bias[ptsz->nz_ng_bias-1];
result = 0.;
}
else if (k<ptsz->array_ln_k_ng_bias[0]){
    // k = ptsz->array_ln_k_ng_bias[0];
result = 0.;
}
else if (k>ptsz->array_ln_k_ng_bias[ptsz->nk_ng_bias-1]){
    // k =  ptsz->array_ln_k_ng_bias[ptsz->nk_ng_bias-1];
result = 0.;
  }
else{

 result = exp(pwl_interp_2d(ptsz->nz_ng_bias,
                          ptsz->nk_ng_bias,
                          ptsz->array_ln_1pz_ng_bias,
                          ptsz->array_ln_k_ng_bias,
                          ptsz->array_ln_ng_bias_at_z_and_k,
                          1,
                          &z,
                          &k))*(bh-1.);
}
return result;
}

//
// double get_ng_bias_contribution_at_z_and_k(double z,
//                                            double k,
//                                            double bh,
//                                            struct background * pba,
//                                            struct perturbs * ppt,
//                                            struct tszspectrum * ptsz){
//
// double fNL = ptsz->fNL;
// // double bh = get_first_order_bias_at_z_and_nu(z,nu,ptsz);
// double beta_f = 2.*ptsz->delta_cSZ*(bh-1.);
// double alpha_k = 1.;
//
//
// // start collecting transfer functions
// char titles[_MAXTITLESTRINGLENGTH_]={0};
// double * data;
// int size_data;
// int number_of_titles = ptsz->number_of_titles;
// int index_md;
// int index_title;
// int index_d_tot = ptsz->index_d_tot;
// int index_phi = ptsz->index_phi;
//
// // class_call(perturb_output_titles(pba,ppt,class_format,titles),
// //            ptsz->error_message,
// //            ptsz->error_message);
// //
// // // printf("ok titles diones\n");
// //
// // number_of_titles = get_number_of_titles(titles);
// // // printf("number_of_titles  = %d %s\n",
// // // number_of_titles,
// // // titles);
// // char *pch;
// // pch = strtok(titles,"\t");
// // int id = 0;
// //
// // while( pch != NULL ) {
// //       // printf( "%s\n",pch );
// //       // strcpy(string1,pch);
// //       // printf( "%s\n",string1 );
// //       if (strstr(pch,"d_tot") != 0)
// //         index_d_tot = id;
// //       if (strstr(pch,"phi") != 0)
// //         index_phi = id;
// //       pch = strtok(NULL,"\t");
// //       id+=1;
// //    }
//
// // printf("index_d_tot = %d index_phi = %d\n",
// //       index_d_tot,index_phi);
// // exit(0);
// //
// // for (index_title=0; index_title<number_of_titles; index_title++){
// // printf("title %s\n",pch);
// // }
//
// // printf("finnished printing titles\n");
//
// index_md=ppt->index_md_scalars;
// size_data = number_of_titles*ppt->k_size[index_md];
//
// if (ppt->ic_size[index_md] != 1){
//   printf("Please run with only one type of initial conditions to avoid confusion in class_sz.\n");
//   exit(0);
// }
//
// class_alloc(data, sizeof(double)*ppt->ic_size[index_md]*size_data, ptsz->error_message);
// // z = 1.;
// perturb_output_data(pba,
//                     ppt,
//                     class_format,
//                     z,
//                     number_of_titles,
//                     data);
//
// int index_k;
// double alpha_kp;
// double kp;
// double * array_ln_alpha_k;
// double * array_ln_k;
// class_alloc(array_ln_alpha_k,sizeof(double *)*ppt->k_size[index_md],ptsz->error_message);
// class_alloc(array_ln_k,sizeof(double *)*ppt->k_size[index_md],ptsz->error_message);
// // printf("starting loop k_size = %d\n",ppt->k_size[index_md]);
// for (index_k=0; index_k<ppt->k_size[index_md]; index_k++){
// // index_title = index_d_tot;
// // printf("d_tot k = %.5e h/Mpc data = %.5e\n",
// // ppt->k[index_md][index_tau]/pba->h,
// // data[index_tau*number_of_titles+index_title]);
// // index_title = index_phi;
// // printf("phi k = %.5e h/Mpc data = %.5e\n",
// // ppt->k[index_md][index_tau]/pba->h,
// // data[index_tau*number_of_titles+index_title]);
// kp = ppt->k[index_md][index_k]/pba->h;
// alpha_kp = data[index_k*number_of_titles+index_d_tot]/data[index_k*number_of_titles+index_phi];
//
// // printf("%d %d %d %d %d %d\n",
// // index_k*number_of_titles+index_d_tot,
// // index_k*number_of_titles+index_phi,
// // size_data,
// // index_k,
// // index_d_tot,
// // index_phi
// // );
// double lp = data[index_k*number_of_titles+index_d_tot]/data[index_k*number_of_titles+index_phi];
//
// if (isnan(alpha_kp)||isinf(alpha_kp)){
//   printf("alpha_kp = %.5e den = %.5e num = %.5e k = %.5e z = %.5e\n",
//          alpha_kp,
//          data[index_k*number_of_titles+index_phi],
//          data[index_k*number_of_titles+index_d_tot],
//          kp,
//          z
//        );
//   exit(0);
// }
// else{
// if (lp>0){
//   printf("alpha>0\n");
//   exit(0);
// }
// array_ln_alpha_k[index_k] = log(-alpha_kp);
// }
//
// // array_ln_alpha_k[index_k] = log(10);
//
// array_ln_k[index_k] = log(kp);
//
// // printf("alpha k = %.5e h/Mpc data = %.5e\n",
// // array_ln_k[index_k],
// // array_ln_alpha_k[index_k]);
// // printf("\n");
//
//
// }
//
//
//
// // exit(0);
//
// free(data);
//
//
//
//
//
// double b_ng;
//
// if ((log(k)<array_ln_k[0]) || (log(k)>array_ln_k[ppt->k_size[index_md]-1])){
//   if (ptsz->sz_verbose>3)
//     printf("k outside interpolation assigning bng=0\n");
//   b_ng = 0.;
// }
// else{
// alpha_k =  exp(pwl_value_1d(ppt->k_size[index_md],
//                             array_ln_k,
//                             array_ln_alpha_k,
//                             log(k)));
// b_ng = fNL*beta_f/alpha_k;
//     }
//
// free(array_ln_alpha_k);
// free(array_ln_k);
// // printf("z k bng = %.5e %.5e %.5e\n",z,k,b_ng);
//
// double bng_int = get_scale_dependent_bias_at_z_and_k(z,k,bh,ptsz);
// printf("z = %.5e fNL = %.5e beta_f = %.5e k = %.5e alpha_k = %.5e b_ng = %.5e b_ng (interp.) = %.5e\n",z,fNL,beta_f,k,alpha_k,b_ng,bng_int);
// // exit(0);
// return b_ng;
//
//                                            }
//



double get_first_order_bias_at_z_and_nu(double z,
                                        double nu,
                                        struct tszspectrum * ptsz){

   //double nu = exp(pvectsz[ptsz->index_lognu]);
   double nuTink = sqrt(nu); //in T10 paper: nu=delta/sigma while here nu=(delta/sigma)^2

   double Delta;

  if (z>3.) z = 3.;

  // double om0 = ptsz->Omega_m_0;
  // double ol0 = 1.-ptsz->Omega_m_0;
  // double Omega_m_at_z = om0*pow(1.+z,3.)/(om0*pow(1.+z,3.)+ ol0);
  double om0 = ptsz->Omega_m_0;
  double or0 = ptsz->Omega_r_0;
  double ol0 = 1.-om0-or0;
  double Omega_m_at_z = om0*pow(1.+z,3.)/(om0*pow(1.+z,3.)+ ol0 + or0*pow(1.+z,4.));

// Omega_m_0

   //Tinker et al 2008 @ M1600-mean
   if (ptsz->MF==6)
   Delta = 1600.;

   //Jenkins et al 2001
   else if (ptsz->MF==3)
   Delta = 180.;

   //Tinker et al 2008 @ m500
   //Boquet et al 2015
   else if (ptsz->MF==5 || ptsz->MF==7)
   Delta = 500./Omega_m_at_z;
   else if (ptsz->MF==8) // 200c
   Delta = 200./Omega_m_at_z;
   else
   Delta = 200.;


   // Table 2 of Tinker et al 2010:
   double y = log10(Delta);
   double ATT = 1.+0.24*y*exp(-pow(4./y,4.));
   double aTTT = 0.44*y-0.88;
   double BTT = 0.183;
   double bTTT = 1.5;
   double CTT = 0.019+0.107*y+0.19*exp(-pow(4./y,4.));
   double cTTT = 2.4;
   double bh;
   bh = 1.-ATT*(pow(nuTink,aTTT)/(pow(nuTink,aTTT)+pow(ptsz->delta_cSZ,aTTT)))
                                                      +BTT*pow(nuTink,bTTT)+CTT*pow(nuTink,cTTT);

   return bh;
  }



int evaluate_halo_bias(double * pvecback,
                       double * pvectsz,
                       struct background * pba,
                       struct primordial * ppm,
                       struct nonlinear * pnl,
                       struct perturbs * ppt,
                       struct tszspectrum * ptsz)
{

   double nu = exp(pvectsz[ptsz->index_lognu]);
   double nuTink = sqrt(nu); //in T10 paper: nu=delta/sigma while here nu=(delta/sigma)^2
   double z = pvectsz[ptsz->index_z];
   double bh = get_first_order_bias_at_z_and_nu(z,nu,ptsz);
   pvectsz[ptsz->index_halo_bias] = bh;//get_first_order_bias_at_z_and_nu(z,nu,ptsz);



   // if fNL bias:
   if (ptsz->has_ng_in_bh){
    double bng = 0.;
    int index_l = (int) pvectsz[ptsz->index_multipole];
    int index_md = (int) pvectsz[ptsz->index_md];
    double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi
    double kl;
    if (_pk_at_z_1h_
     || _pk_gg_at_z_1h_
     || _pk_at_z_2h_
     || _pk_gg_at_z_2h_
     || _pk_bb_at_z_1h_
     || _pk_bb_at_z_2h_
     || _pk_b_at_z_2h_
     || _pk_em_at_z_1h_
     || _pk_em_at_z_2h_
     || _pk_HI_at_z_1h_
     || _pk_HI_at_z_2h_
     || _bk_at_z_1h_
     || _bk_at_z_2h_
     || _bk_at_z_3h_
     || _bk_ttg_at_z_1h_
     || _bk_ttg_at_z_2h_
     || _bk_ttg_at_z_3h_
   ){
    int index_k = (int) pvectsz[ptsz->index_k_for_pk_hm];
      kl = ptsz->k_for_pk_hm[index_k];
       }
    else{
      kl = (ptsz->ell[index_l]+0.5)/d_A;
     }

    // bng = get_ng_bias_contribution_at_z_and_k(z,kl,bh,pba,ppt,ptsz);
    bng = get_scale_dependent_bias_at_z_and_k(z,kl,bh,ptsz);
    pvectsz[ptsz->index_halo_bias] += bng;
  }

   return _SUCCESS_;
}
double get_sigma8_at_z(double z,
                      struct tszspectrum * ptsz,
                      struct background * pba){
double z_asked = log(1.+z);
double R_asked = log(8./pba->h); //log(R) in Mpc
double sigma8_at_z =  exp(pwl_interp_2d(ptsz->n_arraySZ,
                           ptsz->ndimSZ,
                           ptsz->array_redshift,
                           ptsz->array_radius,
                           ptsz->array_sigma_at_z_and_R,
                           1,
                           &z_asked,
                           &R_asked));
return sigma8_at_z;
                         }

double get_sigma_at_z_and_m(double z,
                            double m,
                            struct tszspectrum * ptsz,
                            struct background * pba){

  double rh;
  if (ptsz->HMF_prescription_NCDM == 0) //Matter
    rh = pow(3.*m/(4*_PI_*(pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0),1./3.);

  else if (ptsz->HMF_prescription_NCDM == 1) //CDM
    rh = pow(3.*m/(4*_PI_*(pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0),1./3.);

  else if (ptsz->HMF_prescription_NCDM == 2) //No-pres
    rh = pow(3.*m/(4*_PI_*ptsz->Omega_m_0*ptsz->Rho_crit_0),1./3.);

   double z_asked = log(1.+z);
   // double R_asked = log(exp(log(rh))/pba->h);
   double R_asked = log(rh/pba->h); // this is in Mpc


 if (z_asked<ptsz->array_redshift[0]){
    // z_asked = ptsz->array_redshift[0];
printf("get_sigm: z_asked<ptsz->array_redshift[0].. check bounds.\n");
exit(0);
  }
 if (z_asked>ptsz->array_redshift[ptsz->n_arraySZ-1]){
    // z_asked = ptsz->array_redshift[ptsz->n_arraySZ-1];
printf("get_sigm: z_asked>ptsz->array_redshift[ptsz->n_arraySZ-1].. check bounds.\n");
printf("z_asked = %.15e ptsz->array_redshift[ptsz->n_arraySZ-1] = %.15e\n",exp(z_asked)-1.,exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.);
exit(0);
  }

 if (R_asked<ptsz->array_radius[0]){
    // R_asked = ptsz->array_radius[0];
printf("get_sigm: R_asked<ptsz->array_radius[0].. check bounds.\n");
exit(0);
  }
      // printf("dealing with mass conversion in hmf3\n");
 if (R_asked>ptsz->array_radius[ptsz->ndimSZ-1]){
   // R_asked =  ptsz->array_radius[ptsz->ndimSZ-1];
printf("get_sigm: R_asked>ptsz->array_radius[ptsz->ndimSZ-1].. check bounds.\n");
exit(0);
 }



   double sigma = exp(pwl_interp_2d(ptsz->n_arraySZ,
                                    ptsz->ndimSZ,
                                    ptsz->array_redshift,
                                    ptsz->array_radius,
                                    ptsz->array_sigma_at_z_and_R,
                                    1,
                                    &z_asked,
                                    &R_asked));
  if (isnan(sigma) || isinf(sigma)){
    printf("failed interpolation of sigma.\n");
    printf("z=%.8e zmin=%.8e m=%.8e\n",z,ptsz->array_redshift[0],m);
    exit(0);
  }

   return sigma;
                            }


double get_dlnsigma_dlnR_at_z_and_m(double z,
                                    double m,
                                    struct tszspectrum * ptsz,
                                    struct background * pba){

  double rh;
  if (ptsz->HMF_prescription_NCDM == 0) //Matter
    rh = pow(3.*m/(4*_PI_*(pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0),1./3.);

  else if (ptsz->HMF_prescription_NCDM == 1) //CDM
    rh = pow(3.*m/(4*_PI_*(pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0),1./3.);

  else if (ptsz->HMF_prescription_NCDM == 2) //No-pres
    rh = pow(3.*m/(4*_PI_*ptsz->Omega_m_0*ptsz->Rho_crit_0),1./3.);

   double z_asked = log(1.+z);
   // double R_asked = log(exp(log(rh))/pba->h);
   double R_asked = log(rh/pba->h); // in Mpc


 // if (z_asked<ptsz->array_redshift[0])
 //    z_asked = ptsz->array_redshift[0];
 // if (z_asked>ptsz->array_redshift[ptsz->n_arraySZ-1])
 //    z_asked = ptsz->array_redshift[ptsz->n_arraySZ-1];
 //
 // if (R_asked<ptsz->array_radius[0])
 //    R_asked = ptsz->array_radius[0];
 //      // printf("dealing with mass conversion in hmf3\n");
 // if (R_asked>ptsz->array_radius[ptsz->ndimSZ-1])
 //    R_asked =  ptsz->array_radius[ptsz->ndimSZ-1];

 if (z_asked<ptsz->array_redshift[0]){
    // z_asked = ptsz->array_redshift[0];
printf("get_dlnsigm: z_asked<ptsz->array_redshift[0].. check bounds.\n");
exit(0);
  }
 if (z_asked>ptsz->array_redshift[ptsz->n_arraySZ-1]){
    // z_asked = ptsz->array_redshift[ptsz->n_arraySZ-1];
printf("get_dlnsigm: z_asked>ptsz->array_redshift[ptsz->n_arraySZ-1].. check bounds.\n");
exit(0);
  }

 if (R_asked<ptsz->array_radius[0]){
    // R_asked = ptsz->array_radius[0];
printf("get_dlnsigm: R_asked<ptsz->array_radius[0].. check bounds.\n");
exit(0);
  }
      // printf("dealing with mass conversion in hmf3\n");
 if (R_asked>ptsz->array_radius[ptsz->ndimSZ-1]){
   // R_asked =  ptsz->array_radius[ptsz->ndimSZ-1];
printf("get_dlnsigm: R_asked>ptsz->array_radius[ptsz->ndimSZ-1].. check bounds.\n");
exit(0);
 }



   double sigma =  pwl_interp_2d(
                   ptsz->n_arraySZ,
                   ptsz->ndimSZ,
                   ptsz->array_redshift,
                   ptsz->array_radius,
                   ptsz->array_dsigma2dR_at_z_and_R,
                   1,
                   &z_asked,
                   &R_asked
                   );
  sigma = sigma/2.;


  if (isnan(sigma) || isinf(sigma)){
    printf("failed interpolation of sigma.\n");
    printf("z=%.8e zmin=%.8e m=%.8e\n",z,ptsz->array_redshift[0],m);
    exit(0);
  }

   return sigma;
                            }




double get_nu_at_z_and_m(double z,
                         double m,
                         struct tszspectrum * ptsz,
                         struct background * pba){

  double * pvecback;
  double tau;
  int first_index_back = 0;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);


  class_call(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             pba->error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &first_index_back,
                               pvecback),
             pba->error_message,
             pba->error_message);


double sigma = get_sigma_at_z_and_m(z,m,ptsz,pba);
double lnsigma2 = log(sigma*sigma);
//Nakamura-Suto or not
double delta_c = ptsz->delta_cSZ;//*(1.+0.012299*log10(pvecback[pba->index_bg_Omega_m]));
// double lnnu = 2.*log(ptsz->delta_cSZ) - lnsigma2;
double lnnu = 2.*log(delta_c) - lnsigma2;
double nu = exp(lnnu);
free(pvecback);
return nu;
                         }



double get_second_order_bias_at_z_and_nu(double z,
                                         double nu,
                                         struct tszspectrum * ptsz,
                                         struct background * pba){
if (z>3.) z=3.;

//  double beta = ptsz->beta0SZ*pow(1.+z,0.2);
//  double gamma = ptsz->gamma0SZ*pow(1.+z,-0.01);
//  double eta = ptsz->eta0SZ*pow(1.+z,0.27);
//  double phi = ptsz->phi0SZ*pow(1.+z,-0.08);
//
//   double * pvecback;
//   double tau;
//   int first_index_back = 0;
//   class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
//
//
//   class_call(background_tau_of_z(pba,z,&tau),
//              pba->error_message,
//              pba->error_message);
//
//   class_call(background_at_tau(pba,
//                                tau,
//                                pba->long_info,
//                                pba->inter_normal,
//                                &first_index_back,
//                                pvecback),
//              pba->error_message,
//              pba->error_message);
//
//
// // Nakamura Suto
//  // double delta_c = ptsz->delta_cSZ*(1.+0.012299*log10(pvecback[pba->index_bg_Omega_m]));
 double delta_c = ptsz->delta_cSZ;
//  // printf("z = %.3e delc = %.3e omega_m = %.3e\n",z, delta_c,pvecback[pba->index_bg_Omega_m]);
//  free(pvecback);
//  double t1_num = 2.*(42.*gamma*nu*phi
//                      +8.*delta_c*phi
//                      -84.*eta*phi
//                      +42.*phi*phi
//                      -21.*phi);
//  double t1_den = 21.*delta_c*delta_c*(1.+pow(beta*sqrt(nu),2.*phi));
//  double t2_num = 21.*gamma*gamma*pow(nu,2.)+8.*gamma*delta_c*nu-84*gamma*eta*nu-63.*gamma*nu;
//  double t2_den = 21.*pow(delta_c,2.);
//  double t3_num = -16.*delta_c*eta-8.*delta_c+84.*eta*eta+42*eta;
//  double t3_den = 21.*delta_c*delta_c;
//
//  return t1_num/t1_den + t2_num/t2_den + t3_num/t3_den;

  // K. Hoffmann, J. Bel and E. Gaztanaga 2015
  // 2(1+a2)(ε1 +E1)+ε2 +E2
  // if (ptsz->no_b2){
  //   return 0.;
  // }
  // else{
  double a2 = -17./21.;
  //νf(ν) = A[1+(bν)^a]ν^d*e^(−cν/2),
  double a = -ptsz->phi0SZ*pow(1.+z,-0.08);
  double b = pow(ptsz->beta0SZ*pow(1.+z,0.2),2.);
  double c = ptsz->gamma0SZ*pow(1.+z,-0.01);
  double d = ptsz->eta0SZ*pow(1.+z,0.27)+0.5;
  double epsilon_1 = (c*nu-2.*d)/delta_c;
  double epsilon_2 = (c*nu*(c*nu-4.*d-1.)+2.*d*(2.*d-1.))/(delta_c*delta_c);
  double E1 = -2.*a/(delta_c*(pow(b*nu,-a)+1.));
  double E2 = E1*(-2.*a+2.*c*nu-4.*d+1.)/delta_c;
  return 2.*(1.+a2)*(epsilon_1+E1)+epsilon_2+E2;
// }
  // return 0.;
                                     }



int evaluate_halo_bias_b2(double * pvecback,
                          double * pvectsz,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct tszspectrum * ptsz)
{

   double nu = exp(pvectsz[ptsz->index_lognu]);
   // double nuTink = sqrt(nu); //in T10 paper: nu=delta/sigma while here nu=(delta/sigma)^2

   double z = pvectsz[ptsz->index_z];
   pvectsz[ptsz->index_halo_bias_b2] = get_second_order_bias_at_z_and_nu(z,nu,ptsz,pba);


   return _SUCCESS_;
}

double get_pk_nonlin_at_k_and_z(double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct tszspectrum * ptsz){

if (ptsz->use_class_sz_fast_mode){
double pk = get_pk_lin_at_k_and_z_fast(k,z,pba,ppm,pnl,ptsz);
return pk;
}
else{
if ((k*pba->h < exp(pnl->ln_k[0])) || (k*pba->h > exp(pnl->ln_k[pnl->k_size-1]))){
  return 0.;
}

double pk;
double * pk_ic = NULL;
          class_call(nonlinear_pk_at_k_and_z(
                                            pba,
                                            ppm,
                                            pnl,
                                            pk_nonlinear,
                                            k*pba->h,
                                            z,
                                            pnl->index_pk_m,
                                            &pk, // number *out_pk_l
                                            pk_ic // array out_pk_ic_l[index_ic_ic]
                                          ),
                                          pnl->error_message,
                                          pnl->error_message);
return pk*pow(pba->h,3.);
}
                          }


double get_pk_lin_at_k_and_z_fast(double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct tszspectrum * ptsz){
  if ((k*pba->h < exp(ptsz->array_lnk[0])) || (k*pba->h > exp(ptsz->array_lnk[ptsz->ndimSZ-1]))){
    return 0.;
  }

   double zp = log(1.+z);
   double kp = log(k*pba->h);
   double pk = pwl_interp_2d(ptsz->n_arraySZ,
                      ptsz->ndimSZ,
                      ptsz->array_redshift,
                      ptsz->array_lnk,
                      ptsz->array_pkl_at_z_and_k,
                      1,
                      &zp,
                      &kp);

return pk*pow(pba->h,3.);
}

double get_pk_nonlin_at_k_and_z_fast(double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct tszspectrum * ptsz){
  if ((k*pba->h < exp(ptsz->array_lnk[0])) || (k*pba->h > exp(ptsz->array_lnk[ptsz->ndimSZ-1]))){
    return 0.;
  }

   double zp = log(1.+z);
   double kp = log(k*pba->h);
   double pk = pwl_interp_2d(ptsz->n_arraySZ,
                      ptsz->ndimSZ,
                      ptsz->array_redshift,
                      ptsz->array_lnk,
                      ptsz->array_pknl_at_z_and_k,
                      1,
                      &zp,
                      &kp);

return pk*pow(pba->h,3.);
}


double get_pk_lin_at_k_and_z(double k, double z,
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear * pnl,
                          struct tszspectrum * ptsz){



if (ptsz->use_class_sz_fast_mode){
double pk = get_pk_lin_at_k_and_z_fast(k,z,pba,ppm,pnl,ptsz);
return pk;
}
else{
double pk;
if ((k*pba->h < exp(pnl->ln_k[0])) || (k*pba->h > exp(pnl->ln_k[pnl->k_size-1]))){
  return 0.;
}

double * pk_ic = NULL;
// printf("before\n");
        class_call(nonlinear_pk_at_k_and_z(
                                            pba,
                                            ppm,
                                            pnl,
                                            pk_linear,
                                            k*pba->h,
                                            z,
                                            pnl->index_pk_m,
                                            &pk, // number *out_pk_l
                                            pk_ic // array out_pk_ic_l[index_ic_ic]
                                          ),
                                          ptsz->error_message,
                                          ptsz->error_message);
// printf("k=%.3e z=%.3e pklin=%.3e\n",k,z,pk*pow(pba->h,3.));
// printf("after\n");

// evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);

return pk*pow(pba->h,3.);
}
                          }


// //
// //Compute P(k) in units of h^-3 Mpc^3
// int evaluate_pk_at_ell_plus_one_half_over_chi(double * pvecback,
//                                               double * pvectsz,
//                                               struct background * pba,
//                                               struct primordial * ppm,
//                                               struct nonlinear * pnl,
//                                               struct tszspectrum * ptsz)
// {
//
//   double k;
//   int index_md = (int) pvectsz[ptsz->index_md];
//   double z = pvectsz[ptsz->index_z];
//
//
// if (_pk_at_z_2h_
//   ||_pk_gg_at_z_2h_
//   || _bk_at_z_2h_
//   || _bk_at_z_3h_
//   || _bk_ttg_at_z_2h_
//   || _bk_ttg_at_z_3h_){
//   int index_k = (int) pvectsz[ptsz->index_k_for_pk_hm];
//   k = ptsz->k_for_pk_hm[index_k];
//   // printf("exiting\n");
//   // exit(0);
//
// }
// else {
//    //int index_l = (int)  pvectsz[ptsz->index_multipole];
//    //identical to sqrt(pvectsz[index_chi2])
//    double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi
//
//    pvectsz[ptsz->index_k_value_for_halo_bias] = (pvectsz[ptsz->index_multipole_for_pk]+0.5)/d_A; //units h/Mpc
//
//
//    k = pvectsz[ptsz->index_k_value_for_halo_bias]; //in h/Mpc
//  }
//
//    double pk;
//    double * pk_ic = NULL;
//
//    //printf("evaluating pk at k=%.3e h/Mpc\n",k);
//
//
//
//    //int index_md = (int) pvectsz[ptsz->index_md];
//    if (
//        ((ptsz->has_gal_gal_hf == _TRUE_) && (index_md == ptsz->index_md_gal_gal_hf)) // WIxSC
//     || ((ptsz->has_lens_lens_2h == _TRUE_) && (index_md == ptsz->index_md_lens_lens_2h) && ptsz->use_hod == 0)
//     || ((ptsz->has_gal_lens_hf == _TRUE_) && (index_md == ptsz->index_md_gal_lens_hf))
//     || ((ptsz->has_gal_lensmag_hf == _TRUE_) && (index_md == ptsz->index_md_gal_lensmag_hf))
//     || ((ptsz->has_lens_lensmag_hf == _TRUE_) && (index_md == ptsz->index_md_lens_lensmag_hf))
//     || ((ptsz->has_lensmag_lensmag_hf == _TRUE_) && (index_md == ptsz->index_md_lensmag_lensmag_hf))
//     || ((ptsz->has_kSZ_kSZ_gal_hf == _TRUE_) && (index_md == ptsz->index_md_kSZ_kSZ_gal_hf))) // unWISE
//
//    {
//        //Input: wavenumber in 1/Mpc
//        //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// // printf("evaluating pk at z=%.3e and k=%.3e\n",z,k);
//
//           // halofit approach: uses non-linear pk
//           class_call(nonlinear_pk_at_k_and_z(
//                                             pba,
//                                             ppm,
//                                             pnl,
//                                             pk_nonlinear,
//                                             //pk_linear,
//                                             k*pba->h,
//                                             z,
//                                             pnl->index_pk_m,
//                                             &pk, // number *out_pk_l
//                                             pk_ic // array out_pk_ic_l[index_ic_ic]
//                                           ),
//                                           pnl->error_message,
//                                           pnl->error_message);
//
//   // printf("evaluating pk at z=%.3e and k=%.3e, getting pk =%.3e\n",z,k,pk);
//
//    }
//    else{
//
//        //Input: wavenumber in 1/Mpc
//        //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
//        // printf("z=%.3e\n",z);
//           class_call(nonlinear_pk_at_k_and_z(
//                                             pba,
//                                             ppm,
//                                             pnl,
//                                             pk_linear,
//                                             k*pba->h,
//                                             z,
//                                             pnl->index_pk_m,
//                                             &pk, // number *out_pk_l
//                                             pk_ic // array out_pk_ic_l[index_ic_ic]
//                                           ),
//                                           pnl->error_message,
//                                           pnl->error_message);
// }
//
//
//
//    //now compute P(k) in units of h^-3 Mpc^3
//    pvectsz[ptsz->index_pk_for_halo_bias] = pk*pow(pba->h,3.); //in units Mpc^3/h^3
//    return _SUCCESS_;
// }


//
//
// double evaluate_pk_halofit_over_pk_linear_at_ell_plus_one_half_over_chi(double * pvecback,
//                                                                      double * pvectsz,
//                                                                      struct background * pba,
//                                                                      struct primordial * ppm,
//                                                                      struct nonlinear * pnl,
//                                                                      struct tszspectrum * ptsz)
// {
//
//
//    int index_l = (int)  pvectsz[ptsz->index_multipole];
//    double z = pvectsz[ptsz->index_z];
//    //identical to sqrt(pvectsz[index_chi2])
//    double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi
//
//    pvectsz[ptsz->index_k_value_for_halo_bias] = (ptsz->ell[index_l]+0.5)/d_A; //units h/Mpc
//
//
//    double k = pvectsz[ptsz->index_k_value_for_halo_bias]; //in h/Mpc
//
//    double pk;
//    double * pk_ic = NULL;
//
//    //printf("evaluating pk at k=%.3e h/Mpc\n",k);
//
//
//
//    int index_md = (int) pvectsz[ptsz->index_md];
//
//
//           class_call(nonlinear_pk_at_k_and_z(
//                                             pba,
//                                             ppm,
//                                             pnl,
//                                             //pk_nonlinear,
//                                             pk_nonlinear,
//                                             k*pba->h,
//                                             z,
//                                             pnl->index_pk_cb,
//                                             &pk, // number *out_pk_l
//                                             pk_ic // array out_pk_ic_l[index_ic_ic]
//                                           ),
//                                           pnl->error_message,
//                                           pnl->error_message);
//
//    double pk_halofit = pk;
//
//
//
//        //Input: wavenumber in 1/Mpc
//        //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
//           class_call(nonlinear_pk_at_k_and_z(
//                                             pba,
//                                             ppm,
//                                             pnl,
//                                             pk_linear,
//                                             k*pba->h,
//                                             z,
//                                             pnl->index_pk_cb,
//                                             &pk, // number *out_pk_l
//                                             pk_ic // array out_pk_ic_l[index_ic_ic]
//                                           ),
//                                           pnl->error_message,
//                                           pnl->error_message);
//
//     double pk_lin  = pk;
//
//
//
//
//    double  result = pk_halofit/pk_lin;
//    return result;
// }
//
//
//


int evaluate_pk_at_ell_plus_one_half_over_chi_today(double * pvecback,
                                                   double * pvectsz,
                                                   struct background * pba,
                                                   struct primordial * ppm,
                                                   struct nonlinear * pnl,
                                                   struct tszspectrum * ptsz)
{

   int index_l = (int)  pvectsz[ptsz->index_multipole];
   double z = pvectsz[ptsz->index_z];
   //identical to sqrt(pvectsz[index_chi2])
   double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi

   pvectsz[ptsz->index_k_value_for_halo_bias] = (ptsz->ell[index_l]+0.5)/d_A; //units h/Mpc


   double k = pvectsz[ptsz->index_k_value_for_halo_bias]; //in h/Mpc

   double pk;
   double * pk_ic = NULL;


  //Input: wavenumber in 1/Mpc
  //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
   class_call(nonlinear_pk_at_k_and_z(
                                     pba,
                                     ppm,
                                     pnl,
                                     pk_linear,
                                     k*pba->h,
                                     0.,
                                     pnl->index_pk_m,
                                     &pk, // number *out_pk_l
                                     pk_ic // array out_pk_ic_l[index_ic_ic]
                                   ),
                                   pnl->error_message,
                                   pnl->error_message);


   //now compute P(k) in units of h^-3 Mpc^3
   pvectsz[ptsz->index_pk_for_halo_bias] = pk*pow(pba->h,3.); //in units Mpc^3/h^3
   return _SUCCESS_;
}


int do_mass_conversions( double logM,
                         double z,
                         double * pvecback,
                         double * pvectsz,
                         struct background * pba,
                         struct tszspectrum * ptsz)
{
//
// if (pvectsz[ptsz->index_has_electron_pressure] == 1){
//         // if battaglia pressure profile need 200c
//         // and virial quantities for integration bound
//       if (ptsz->pressure_profile == 4){
//         pvectsz[ptsz->index_has_200c] = 1;
//         pvectsz[ptsz->index_has_vir] = 1;
//       }
//       else if (ptsz->pressure_profile == 0){ // Planck 2013
//         pvectsz[ptsz->index_has_500c] = 1;
//       }
//       else if (ptsz->pressure_profile == 2){ // Arnaud et al 2010
//         pvectsz[ptsz->index_has_500c] = 1;
//       }
//       else if (ptsz->pressure_profile == 3){ // custom
//         pvectsz[ptsz->index_has_500c] = 1;
//       }
// }
//
// if (pvectsz[ptsz->index_has_electron_density] == 1){
//   // if (ptsz->tau_profile == 1){ // battaglia tau profile, need m200c
//   pvectsz[ptsz->index_has_200c] = 1;
//   // }
// }
//
// if (pvectsz[ptsz->index_has_cib] == 1){
//   pvectsz[ptsz->index_has_200m] = 1;
// }


if (pvectsz[ptsz->index_has_galaxy] == 1){
  if (ptsz->delta_def_galaxies == 0)
    pvectsz[ptsz->index_has_200m] = 1;
  else if (ptsz->delta_def_galaxies == 1)
    pvectsz[ptsz->index_has_200c] = 1;
  else if (ptsz->delta_def_galaxies == 2)
    pvectsz[ptsz->index_has_500c] = 1;
}

if (pvectsz[ptsz->index_has_cib] == 1){
  if (ptsz->delta_def_cib == 0)
    pvectsz[ptsz->index_has_200m] = 1;
  else if (ptsz->delta_def_cib == 1)
    pvectsz[ptsz->index_has_200c] = 1;
  else if (ptsz->delta_def_cib == 2)
    pvectsz[ptsz->index_has_500c] = 1;
}


if (pvectsz[ptsz->index_has_matter_density] == 1 || pvectsz[ptsz->index_has_lensing] == 1){
  if (ptsz->delta_def_matter_density == 0)
    pvectsz[ptsz->index_has_200m] = 1;
  else if (ptsz->delta_def_matter_density == 1)
    pvectsz[ptsz->index_has_200c] = 1;
  else if (ptsz->delta_def_matter_density == 2)
    pvectsz[ptsz->index_has_500c] = 1;
}

if (pvectsz[ptsz->index_has_electron_density] == 1){
  if (ptsz->delta_def_electron_density == 0)
    pvectsz[ptsz->index_has_200m] = 1;
  else if (ptsz->delta_def_electron_density == 1)
    pvectsz[ptsz->index_has_200c] = 1;
  else if (ptsz->delta_def_electron_density == 2)
    pvectsz[ptsz->index_has_500c] = 1;
}

if (pvectsz[ptsz->index_has_HI_density] == 1){
  if (ptsz->delta_def_HI_density == 0)
    pvectsz[ptsz->index_has_200m] = 1;
  else if (ptsz->delta_def_HI_density == 1)
    pvectsz[ptsz->index_has_200c] = 1;
  else if (ptsz->delta_def_HI_density == 2)
    pvectsz[ptsz->index_has_500c] = 1;
}

if (pvectsz[ptsz->index_has_electron_pressure] == 1){
  if (ptsz->delta_def_electron_pressure == 0)
    pvectsz[ptsz->index_has_200m] = 1;
  else if (ptsz->delta_def_electron_pressure == 1)
    pvectsz[ptsz->index_has_200c] = 1;
  else if (ptsz->delta_def_electron_pressure == 2)
    pvectsz[ptsz->index_has_500c] = 1;

if (ptsz->pressure_profile == 4){
      pvectsz[ptsz->index_has_vir] = 1; // cut wrt virial radius
      }
}

// if (ptsz->MF == 1) {
//   pvectsz[ptsz->index_has_200m] = 1;
// }
// //Tinker et al 2008 @ m500
// //Boquet et al 2015
// else if (ptsz->MF==5 || ptsz->MF==7){
//   pvectsz[ptsz->index_has_500c] = 1;
// }


   // deal with mass definitions
   if (ptsz->integrate_wrt_mvir == 1){

    pvectsz[ptsz->index_mVIR] = exp(logM);
    pvectsz[ptsz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[ptsz->index_mVIR],pvectsz[ptsz->index_Delta_c],pvectsz[ptsz->index_Rho_crit],ptsz);
    pvectsz[ptsz->index_cVIR] = evaluate_cvir_of_mvir(pvectsz[ptsz->index_mVIR],z,ptsz,pba);

    if (pvectsz[ptsz->index_has_500c] == 1){
  class_call(mVIR_to_mDEL(pvectsz[ptsz->index_mVIR],
                       pvectsz[ptsz->index_rVIR],
                       pvectsz[ptsz->index_cVIR],
                       500.*(pvectsz[ptsz->index_Rho_crit]),
                       &pvectsz[ptsz->index_m500c],
                       ptsz),
                  ptsz->error_message,
                  ptsz->error_message);
  pvectsz[ptsz->index_r500c] = pow(3.*pvectsz[ptsz->index_m500c]/(4.*_PI_*500.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc
  pvectsz[ptsz->index_c500c] = get_c500c_at_m_and_z(pvectsz[ptsz->index_m500c],z,pba,ptsz);


    }

    if (pvectsz[ptsz->index_has_200c] == 1){

   //compute m200c
  class_call(mVIR_to_mDEL(pvectsz[ptsz->index_mVIR],
                       pvectsz[ptsz->index_rVIR],
                       pvectsz[ptsz->index_cVIR],
                       200.*(pvectsz[ptsz->index_Rho_crit]),
                       &pvectsz[ptsz->index_m200c],
                       ptsz),
                  ptsz->error_message,
                  ptsz->error_message);

    //r200c for the pressure profile B12
    pvectsz[ptsz->index_r200c] = pow(3.*pvectsz[ptsz->index_m200c]/(4.*_PI_*200.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc
    pvectsz[ptsz->index_l200c] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_r200c];
    pvectsz[ptsz->index_c200c] = get_c200c_at_m_and_z(pvectsz[ptsz->index_m200c],z,pba,ptsz);


    }

    if (pvectsz[ptsz->index_has_200m] == 1){
  class_call(mVIR_to_mDEL(pvectsz[ptsz->index_mVIR],
                      pvectsz[ptsz->index_rVIR],
                      pvectsz[ptsz->index_cVIR],
                      200.*pvecback[pba->index_bg_Omega_m]
                      *pvectsz[ptsz->index_Rho_crit],
                      &pvectsz[ptsz->index_m200m],
                      ptsz),
                 ptsz->error_message,
                 ptsz->error_message);
  pvectsz[ptsz->index_r200m] = pow(3.*pvectsz[ptsz->index_m200m]
                                  /(4.*_PI_*200.*pvecback[pba->index_bg_Omega_m]
                                  *pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc


  pvectsz[ptsz->index_c200m] =  get_c200m_at_m_and_z(pvectsz[ptsz->index_m200m],z,pba,ptsz);

    }

  }
   else if (ptsz->integrate_wrt_m500c == 1){

    // printf("entering m500c\n");
    pvectsz[ptsz->index_m500c] = exp(logM);
    // printf("m500c = %.3e\n",pvectsz[ptsz->index_m500c]);
    pvectsz[ptsz->index_r500c] = pow(3.*pvectsz[ptsz->index_m500c]/(4.*_PI_*500.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc
    // printf("r500c = %.3e\n",pvectsz[ptsz->index_r500c]);

    pvectsz[ptsz->index_c500c] = get_c500c_at_m_and_z(pvectsz[ptsz->index_m500c],z,pba,ptsz);


    if (pvectsz[ptsz->index_has_vir] == 1){

    // //convert from m500c to mVIR:
    class_call(mDEL_to_mVIR(pvectsz[ptsz->index_m500c],
                             500.*(pvectsz[ptsz->index_Rho_crit]),
                             pvectsz[ptsz->index_Delta_c],
                             pvectsz[ptsz->index_Rho_crit],
                             z,
                             &pvectsz[ptsz->index_mVIR],
                             ptsz,
                             pba),
                    ptsz->error_message,
                    ptsz->error_message);
    // printf("mvir = %.3e\n",pvectsz[ptsz->index_mVIR]);

    pvectsz[ptsz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[ptsz->index_mVIR],pvectsz[ptsz->index_Delta_c],pvectsz[ptsz->index_Rho_crit],ptsz);
    pvectsz[ ptsz->index_cVIR] = evaluate_cvir_of_mvir(pvectsz[ptsz->index_mVIR],z,ptsz,pba);

    }

    if (pvectsz[ptsz->index_has_200c] == 1){

    pvectsz[ptsz->index_m200c] = get_m500c_to_m200c_at_z_and_M(z,pvectsz[ptsz->index_m500c],ptsz);
        //r200c for the pressure profile B12
    pvectsz[ptsz->index_r200c] = pow(3.*pvectsz[ptsz->index_m200c]/(4.*_PI_*200.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc
    pvectsz[ptsz->index_l200c] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_r200c];
    pvectsz[ptsz->index_c200c] = get_c200c_at_m_and_z(pvectsz[ptsz->index_m200c],z,pba,ptsz);

    }

    if (pvectsz[ptsz->index_has_200m] == 1){
      //exit(0);
    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[ptsz->index_Rho_crit];
    double delrho = 500.*rhoc; // 200m
    double delrho_prime = 200.*omega*rhoc; //500c
    double mdel = pvectsz[ptsz->index_m500c];

    double mdel_prime;
    class_call(mDEL_to_mDELprime(mdel,
                           delrho,
                           delrho_prime,
                           delc,
                           rhoc,
                           z,
                           &mdel_prime,
                           ptsz,
                           pba),
                    ptsz->error_message,
                    ptsz->error_message);
    pvectsz[ptsz->index_m200m] = mdel_prime;
    pvectsz[ptsz->index_r200m]  = pow(3.*pvectsz[ptsz->index_m200m]
                                      /(4.*_PI_*200.*pvecback[pba->index_bg_Omega_m]
                                      *pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc
    pvectsz[ptsz->index_c200m]  = get_c200m_at_m_and_z(pvectsz[ptsz->index_m200m],z,pba,ptsz);


    }


  }
   else if (ptsz->integrate_wrt_m200m == 1){

    pvectsz[ptsz->index_m200m] = exp(logM);
    pvectsz[ptsz->index_r200m]  = pow(3.*pvectsz[ptsz->index_m200m]
                                      /(4.*_PI_*200.*pvecback[pba->index_bg_Omega_m]
                                      *pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc

    pvectsz[ptsz->index_c200m] = get_c200m_at_m_and_z(pvectsz[ptsz->index_m200m],z,pba,ptsz);

    if (pvectsz[ptsz->index_has_vir] == 1){
    class_call(mDEL_to_mVIR(pvectsz[ptsz->index_m200m],
                             200.*pvecback[pba->index_bg_Omega_m]*(pvectsz[ptsz->index_Rho_crit]),
                             pvectsz[ptsz->index_Delta_c],
                             pvectsz[ptsz->index_Rho_crit],
                             z,
                             &pvectsz[ptsz->index_mVIR],
                             ptsz,
                             pba),
                    ptsz->error_message,
                    ptsz->error_message);

    pvectsz[ptsz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[ptsz->index_mVIR],
                                                      pvectsz[ptsz->index_Delta_c],
                                                      pvectsz[ptsz->index_Rho_crit],
                                                      ptsz);
    pvectsz[ ptsz->index_cVIR] = evaluate_cvir_of_mvir(pvectsz[ptsz->index_mVIR],z,ptsz,pba);

    }

    if (pvectsz[ptsz->index_has_200c] == 1){
    pvectsz[ptsz->index_m200c] = get_m200m_to_m200c_at_z_and_M(z,pvectsz[ptsz->index_m200m],ptsz);
        //r200c for the pressure profile B12
    pvectsz[ptsz->index_r200c] = pow(3.*pvectsz[ptsz->index_m200c]/(4.*_PI_*200.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc
    pvectsz[ptsz->index_l200c] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_r200c];
    pvectsz[ptsz->index_c200c] = get_c200c_at_m_and_z(pvectsz[ptsz->index_m200c],z,pba,ptsz);

    }
    if (pvectsz[ptsz->index_has_500c] == 1){

  pvectsz[ptsz->index_m500c] = get_m200m_to_m500c_at_z_and_M(z,pvectsz[ptsz->index_m200m],ptsz);
  pvectsz[ptsz->index_r500c] = pow(3.*pvectsz[ptsz->index_m500c]/(4.*_PI_*500.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc
  pvectsz[ptsz->index_c500c] = get_c500c_at_m_and_z(pvectsz[ptsz->index_m500c],z,pba,ptsz);


    }

  }
 else if (ptsz->integrate_wrt_m200c == 1){

pvectsz[ptsz->index_m200c] = exp(logM);
pvectsz[ptsz->index_r200c]  = pow(3.*pvectsz[ptsz->index_m200c]
                                  /(4.*_PI_*200.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc

pvectsz[ptsz->index_c200c] = get_c200c_at_m_and_z(pvectsz[ptsz->index_m200c],z,pba,ptsz);
pvectsz[ptsz->index_l200c] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_r200c];

    if (pvectsz[ptsz->index_has_200m] == 1){

    pvectsz[ptsz->index_m200m] = get_m200c_to_m200m_at_z_and_M(z,pvectsz[ptsz->index_m200c],ptsz);
    pvectsz[ptsz->index_r200m] = pow(3.*pvectsz[ptsz->index_m200m]/(4.*_PI_*200.*pvecback[pba->index_bg_Omega_m]*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc
    // pvectsz[ptsz->index_l200m] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_r200m];
    pvectsz[ptsz->index_c200m] = get_c200m_at_m_and_z(pvectsz[ptsz->index_m200m],z,pba,ptsz);

    }
    if (pvectsz[ptsz->index_has_500c] == 1){

pvectsz[ptsz->index_m500c] = get_m200c_to_m500c_at_z_and_M(z,pvectsz[ptsz->index_m200c],ptsz);
pvectsz[ptsz->index_r500c]  = pow(3.*pvectsz[ptsz->index_m500c]
                                  /(4.*_PI_*500.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc

pvectsz[ptsz->index_c500c] = get_c500c_at_m_and_z(pvectsz[ptsz->index_m500c],z,pba,ptsz);


    }


}

  double m_for_hmf;

    //Tinker et al 2010
    if (ptsz->MF==1)
      m_for_hmf = pvectsz[ptsz->index_m200m];
      //m_for_hmf = pvectsz[ptsz->index_mVIR]; // BB debug

    //Boquet et al 2015
    if (ptsz->MF==2)
      m_for_hmf = pvectsz[ptsz->index_m200m];

    //Jenkins et al 2001
    if (ptsz->MF==3)
      m_for_hmf = pvectsz[ptsz->index_m180m];

    //Tinker et al 2008
    if (ptsz->MF==4)
      m_for_hmf = pvectsz[ptsz->index_m200m];

     //Tinker et al 2008 @ m500
    if (ptsz->MF==5)
       m_for_hmf = pvectsz[ptsz->index_m500c];

    //Tinker et al 2008 @ M1600-mean
    if (ptsz->MF==6)
      m_for_hmf = pvectsz[ptsz->index_m1600m];

    //Boquet et al 2015
    if (ptsz->MF==7)
      m_for_hmf = pvectsz[ptsz->index_m500c];

    //Tinker et al 2008 @ m200c
    if (ptsz->MF==8)
      m_for_hmf = pvectsz[ptsz->index_m200c];

    pvectsz[ptsz->index_mass_for_hmf] = m_for_hmf;



  double r_delta, c_delta, m_delta;
  if (ptsz->integrate_wrt_mvir == 1){
     m_delta = pvectsz[ptsz->index_mVIR];
     r_delta = pvectsz[ptsz->index_rVIR];
     c_delta = pvectsz[ptsz->index_cVIR];
   }
  else if (ptsz->integrate_wrt_m500c == 1){
    m_delta = pvectsz[ptsz->index_m500c];
    r_delta = pvectsz[ptsz->index_r500c];
    c_delta = pvectsz[ptsz->index_c500c];
  }
  else if (ptsz->integrate_wrt_m200m == 1){
    m_delta = pvectsz[ptsz->index_m200m];
    r_delta = pvectsz[ptsz->index_r200m];
    c_delta = pvectsz[ptsz->index_c200m];
  }
  else if (ptsz->integrate_wrt_m200c == 1){
    m_delta = pvectsz[ptsz->index_m200c];
    r_delta = pvectsz[ptsz->index_r200c];
    c_delta = pvectsz[ptsz->index_c200c];
  }


   double r_delta_gal, c_delta_gal, m_delta_gal;
   double r_delta_electron_pressure, c_delta_electron_pressure, m_delta_electron_pressure; //(tsz)
   double r_delta_matter, c_delta_matter, m_delta_matter;
   double r_delta_lensing, c_delta_lensing, m_delta_lensing;
   double r_delta_cib, c_delta_cib, m_delta_cib;
   double r_delta_electron_density, c_delta_electron_density, m_delta_electron_density; //(ksz)
   double r_delta_HI_density, c_delta_HI_density, m_delta_HI_density;


if (pvectsz[ptsz->index_has_electron_pressure] == 1){

  if (ptsz->delta_def_electron_pressure == 0){
    m_delta_electron_pressure = pvectsz[ptsz->index_m200m];
    r_delta_electron_pressure = pvectsz[ptsz->index_r200m];
    c_delta_electron_pressure = pvectsz[ptsz->index_c200m];
  }
  else if (ptsz->delta_def_electron_pressure == 1){
    m_delta_electron_pressure = pvectsz[ptsz->index_m200c];
    r_delta_electron_pressure = pvectsz[ptsz->index_r200c];
    c_delta_electron_pressure = pvectsz[ptsz->index_c200c];
  }
  else if (ptsz->delta_def_electron_pressure == 2){
    m_delta_electron_pressure = pvectsz[ptsz->index_m500c];
    r_delta_electron_pressure = pvectsz[ptsz->index_r500c];
    c_delta_electron_pressure = pvectsz[ptsz->index_c500c];
  }

  pvectsz[ptsz->index_mass_for_electron_pressure] = m_delta_electron_pressure;
  pvectsz[ptsz->index_radius_for_electron_pressure] = r_delta_electron_pressure;
  pvectsz[ptsz->index_concentration_for_electron_pressure] = c_delta_electron_pressure;

 }// end electron pressure

if (pvectsz[ptsz->index_has_electron_density] == 1){
  if (ptsz->delta_def_electron_density == 0){
    m_delta_electron_density = pvectsz[ptsz->index_m200m];
    r_delta_electron_density = pvectsz[ptsz->index_r200m];
    c_delta_electron_density = pvectsz[ptsz->index_c200m];
  }
  else if (ptsz->delta_def_electron_density == 1){
    m_delta_electron_density = pvectsz[ptsz->index_m200c];
    r_delta_electron_density = pvectsz[ptsz->index_r200c];
    c_delta_electron_density = pvectsz[ptsz->index_c200c];
  }
  else if (ptsz->delta_def_electron_density == 2){
    m_delta_electron_density = pvectsz[ptsz->index_m500c];
    r_delta_electron_density = pvectsz[ptsz->index_r500c];
    c_delta_electron_density = pvectsz[ptsz->index_c500c];
  }

  pvectsz[ptsz->index_mass_for_electron_density] = m_delta_electron_density;
  pvectsz[ptsz->index_radius_for_electron_density] = r_delta_electron_density;
  pvectsz[ptsz->index_concentration_for_electron_density] = c_delta_electron_density;

 }// end electron density

if (pvectsz[ptsz->index_has_HI_density] == 1){
  if (ptsz->delta_def_HI_density == 0){
    m_delta_HI_density = pvectsz[ptsz->index_m200m];
    r_delta_HI_density = pvectsz[ptsz->index_r200m];
    c_delta_HI_density = pvectsz[ptsz->index_c200m];
  }
  else if (ptsz->delta_def_HI_density == 1){
    m_delta_HI_density = pvectsz[ptsz->index_m200c];
    r_delta_HI_density = pvectsz[ptsz->index_r200c];
    c_delta_HI_density = pvectsz[ptsz->index_c200c];
  }
  else if (ptsz->delta_def_HI_density == 2){
    m_delta_HI_density = pvectsz[ptsz->index_m500c];
    r_delta_HI_density = pvectsz[ptsz->index_r500c];
    c_delta_HI_density = pvectsz[ptsz->index_c500c];
  }

  pvectsz[ptsz->index_mass_for_HI_density] = m_delta_HI_density;
  pvectsz[ptsz->index_radius_for_HI_density] = r_delta_HI_density;
  pvectsz[ptsz->index_concentration_for_HI_density] = c_delta_HI_density;

}// end HI density



// matter density
if (pvectsz[ptsz->index_has_matter_density] == 1 || pvectsz[ptsz->index_has_lensing] == 1){
   //nfw profile
  if (ptsz->delta_def_matter_density == 0){
    m_delta_matter = pvectsz[ptsz->index_m200m];
    r_delta_matter = pvectsz[ptsz->index_r200m];
    c_delta_matter = pvectsz[ptsz->index_c200m];
  }
  else if (ptsz->delta_def_matter_density == 1){
    m_delta_matter = pvectsz[ptsz->index_m200c];
    r_delta_matter = pvectsz[ptsz->index_r200c];
    c_delta_matter = pvectsz[ptsz->index_c200c];
  }
  else if (ptsz->delta_def_matter_density == 2){
    m_delta_matter = pvectsz[ptsz->index_m500c];
    r_delta_matter = pvectsz[ptsz->index_r500c];
    c_delta_matter = pvectsz[ptsz->index_c500c];
  }

  pvectsz[ptsz->index_mass_for_matter_density] = m_delta_matter;
  pvectsz[ptsz->index_radius_for_matter_density] = r_delta_matter;
  pvectsz[ptsz->index_concentration_for_matter_density] = c_delta_matter;

  m_delta_lensing = m_delta_matter;
  r_delta_lensing = r_delta_matter;
  c_delta_lensing = c_delta_matter;


 }// end matter density


// cib
if (pvectsz[ptsz->index_has_cib] == 1){
   // cib profile
  if (ptsz->delta_def_cib == 0){
    m_delta_cib = pvectsz[ptsz->index_m200m];
    r_delta_cib = pvectsz[ptsz->index_r200m];
    c_delta_cib = pvectsz[ptsz->index_c200m];
  }
  else if (ptsz->delta_def_cib == 1){
    m_delta_cib = pvectsz[ptsz->index_m200c];
    r_delta_cib = pvectsz[ptsz->index_r200c];
    c_delta_cib = pvectsz[ptsz->index_c200c];
  }
  else if (ptsz->delta_def_cib == 2){
    m_delta_cib = pvectsz[ptsz->index_m500c];
    r_delta_cib = pvectsz[ptsz->index_r500c];
    c_delta_cib = pvectsz[ptsz->index_c500c];
  }

  pvectsz[ptsz->index_mass_for_cib] = m_delta_cib;
  pvectsz[ptsz->index_radius_for_cib] = r_delta_cib;
  pvectsz[ptsz->index_concentration_for_cib] = c_delta_cib;

 }// end cib


// galaxies
if (pvectsz[ptsz->index_has_galaxy] == 1){
   // galaxies
  if (ptsz->delta_def_galaxies == 0){
    m_delta_gal = pvectsz[ptsz->index_m200m];
    r_delta_gal = pvectsz[ptsz->index_r200m];
    c_delta_gal = pvectsz[ptsz->index_c200m];
  }
  else if (ptsz->delta_def_galaxies == 1){
    m_delta_gal = pvectsz[ptsz->index_m200c];
    r_delta_gal = pvectsz[ptsz->index_r200c];
    c_delta_gal = pvectsz[ptsz->index_c200c];
    // printf("setting galaxy's delta\n");
    // printf("c_delta_gal = %.5e\n",c_delta_gal);
  }
  else if (ptsz->delta_def_galaxies == 2){
    m_delta_gal = pvectsz[ptsz->index_m500c];
    r_delta_gal = pvectsz[ptsz->index_r500c];
    c_delta_gal = pvectsz[ptsz->index_c500c];
  }

  pvectsz[ptsz->index_mass_for_galaxies] = m_delta_gal;
  pvectsz[ptsz->index_radius_for_galaxies] = r_delta_gal;
  pvectsz[ptsz->index_concentration_for_galaxies] = ptsz->csat_over_cdm*c_delta_gal;
 }// end galaxies






}




int evaluate_HMF_at_logM_and_z(
                 double logM ,
                 double z,
                 double * pvecback,
                 double * pvectsz,
                 struct background * pba,
                 struct nonlinear * pnl,
                 struct tszspectrum * ptsz)
{

   double m_for_hmf = exp(logM);



  if (ptsz->HMF_prescription_NCDM == 0) //Matter
    pvectsz[ptsz->index_Rh] = pow(3.*m_for_hmf/(4*_PI_*(pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0),1./3.);

  else if (ptsz->HMF_prescription_NCDM == 1) //CDM
    pvectsz[ptsz->index_Rh] = pow(3.*m_for_hmf/(4*_PI_*(pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0),1./3.);

  else if (ptsz->HMF_prescription_NCDM == 2) //No-pres
    pvectsz[ptsz->index_Rh] = pow(3.*m_for_hmf/(4*_PI_*ptsz->Omega_m_0*ptsz->Rho_crit_0),1./3.);


   double z_asked = log(1.+z);

   // double R_asked = log(exp(log(pvectsz[ptsz->index_Rh]))/pba->h);
// printf("dealing with mass conversion in hmf0 %.3e\n",ptsz->array_redshift[0]);
   // if (z<exp(ptsz->array_redshift[0])-1.)
   //    z_asked = ptsz->array_redshift[0];
   //      // printf("dealing with mass conversion in hmf\n");
   // if (z>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
   //    z_asked =  ptsz->array_redshift[ptsz->n_arraySZ-1];
   //      // printf("dealing with mass conversion in hmf2\n");
   // if (log(exp(log(pvectsz[ptsz->index_Rh]))/pba->h)<ptsz->array_radius[0])
   //    R_asked = ptsz->array_radius[0];
   //      // printf("dealing with mass conversion in hmf3\n");
   // if (log(exp(log(pvectsz[ptsz->index_Rh]))/pba->h)>ptsz->array_radius[ptsz->ndimSZ-1])
   //    R_asked =  ptsz->array_radius[ptsz->ndimSZ-1];
   //

  pvectsz[ptsz->index_logSigma2] = 2.*log(get_sigma_at_z_and_m(exp(z_asked)-1.,m_for_hmf,ptsz,pba));

  pvectsz[ptsz->index_lognu] = log(get_nu_at_z_and_m(exp(z_asked)-1.,m_for_hmf,ptsz,pba));

  pvectsz[ptsz->index_dlogSigma2dlogRh] = 2.*get_dlnsigma_dlnR_at_z_and_m(z,m_for_hmf,ptsz,pba);
  // pwl_interp_2d(
  //               ptsz->n_arraySZ,
  //               ptsz->ndimSZ,
  //               ptsz->array_redshift,
  //               ptsz->array_radius,
  //               ptsz->array_dsigma2dR_at_z_and_R,
  //               1,
  //               &z_asked,
  //               &R_asked
  //               );



  pvectsz[ptsz->index_dlogSigma2dlogRh] *= exp(log(pvectsz[ptsz->index_Rh]))/pba->h
                                           /exp(pvectsz[ptsz->index_logSigma2]);


  pvectsz[ptsz->index_dlognudlogRh] = -pvectsz[ptsz->index_dlogSigma2dlogRh];
   //HMF evaluation:
   //Tinker et al 2010
   if (ptsz->MF==1) {
      class_call(
                      MF_T10(
                                 &pvectsz[ptsz->index_mf],
                                 &pvectsz[ptsz->index_lognu],
                                 z,
                                 ptsz
                                 ),
                      ptsz->error_message,
                      ptsz->error_message
                      );
   }

   //HMF evaluation:
   //Boquet et al 2015
   else if (ptsz->MF==2) {
      class_call(
                      MF_B15(
                                 &pvectsz[ptsz->index_mf],
                                 &pvectsz[ptsz->index_lognu],
                                 z,
                                 ptsz
                                 ),
                      ptsz->error_message,
                      ptsz->error_message
                      );
   }

   //HMF evaluation:
   //Jenkins et al 2001
   else if (ptsz->MF==3){
      class_call(
                      MF_J01(
                                 &pvectsz[ptsz->index_mf],
                                 &pvectsz[ptsz->index_lognu],
                                 ptsz
                                 ),
                      ptsz->error_message,
                      ptsz->error_message
                      );
   }
   //HMF evaluation:
   //Tinker et al 2008
   else if (ptsz->MF==4) {
       double zz = z;
       // if(z>3.) zz=3.;
       // double om0 = ptsz->Omega_m_0;
       // double ol0 = 1.-ptsz->Omega_m_0;
       // double Omega_m_z = om0*pow(1.+zz,3.)/(om0*pow(1.+zz,3.)+ ol0);
       double om0 = ptsz->Omega_m_0;
       double or0 = ptsz->Omega_r_0;
       double ol0 = 1.-om0-or0;
       double Omega_m_z = om0*pow(1.+z,3.)/(om0*pow(1.+z,3.)+ ol0 + or0*pow(1.+z,4.));
      class_call(
                      // MF_T08(
                      //            &pvectsz[ptsz->index_mf],
                      //            &pvectsz[ptsz->index_lognu],
                      //            z,
                      //            ptsz
                      //            ),
                      // ptsz->error_message,
                      // ptsz->error_message
                      // );

                      MF_T08_m500(
                                        &pvectsz[ptsz->index_mf],
                                        &pvectsz[ptsz->index_lognu],
                                        z,
                                        200.*Omega_m_z,
                                        ptsz
                                        ),
                      ptsz->error_message,
                      ptsz->error_message
                      );



   }

   //HMF evaluation:
   //Tinker et al 2008 @ m500
   else if (ptsz->MF==5) {
      class_call(
                      MF_T08_m500(
                                        &pvectsz[ptsz->index_mf],
                                        &pvectsz[ptsz->index_lognu],
                                        z,
                                        500.,
                                        ptsz
                                        ),
                      ptsz->error_message,
                      ptsz->error_message
                      );
    // exit(0);
   }


   //HMF evaluation:
   //Tinker et al 2008 @ M1600-mean
   else if (ptsz->MF==6) {
      class_call(
                      MF_T08_M1600m(
                                           &pvectsz[ptsz->index_mf],
                                           &pvectsz[ptsz->index_lognu],
                                           z,
                                           ptsz
                                           ),
                      ptsz->error_message,
                      ptsz->error_message
                      );
   }

   //HMF evaluation:
   //Boquet et al 2015
   else if (ptsz->MF==7) {
      class_call(MF_B15_M500c(&pvectsz[ptsz->index_mf],
                                          &pvectsz[ptsz->index_lognu],
                                          z,ptsz),
                      ptsz->error_message,
                      ptsz->error_message);
   }


   //HMF evaluation:
   //Tinker et al 2008 @ m200c
   else if (ptsz->MF==8) {
      class_call(
                      MF_T08_m500(
                                        &pvectsz[ptsz->index_mf],
                                        &pvectsz[ptsz->index_lognu],
                                        z,
                                        200.,
                                        ptsz
                                        ),
                      ptsz->error_message,
                      ptsz->error_message
                      );
  // printf("calling T08 m200c\n");
   }


pvectsz[ptsz->index_dndlogRh] = 3./(4.*_PI_*pow(pvectsz[ptsz->index_Rh],3))
                                               *pvectsz[ptsz->index_dlognudlogRh]
                                               *pvectsz[ptsz->index_mf];

//Return the HMF - dn/dlogM in units of h^3 Mpc^-3
pvectsz[ptsz->index_hmf] = pvectsz[ptsz->index_dndlogRh]/3.;
// pvectsz[ptsz->index_hmf] = pvectsz[ptsz->index_mf];

return _SUCCESS_;
}



int evaluate_vrms2(double * pvecback,
                   double * pvectsz,
                   struct background * pba,
                   struct nonlinear * pnl,
                   struct tszspectrum * ptsz)
  {

   double z = pvectsz[ptsz->index_z];
   double z_asked = log(1.+z);

   if (z<exp(ptsz->array_redshift[0])-1.)
      z_asked = ptsz->array_redshift[0];
   if (z>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
      z_asked =  ptsz->array_redshift[ptsz->n_arraySZ-1];


   pvectsz[ptsz->index_vrms2] =  exp(pwl_value_1d(ptsz->n_arraySZ,
                                                  ptsz->array_redshift,
                                                  ptsz->array_vrms2_at_z,
                                                  z_asked));

return _SUCCESS_;
}




double get_ttg_bispectrum_at_z_tree_level_PT(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm){

double * pk_ic = NULL;
double result = 0.;
double pk1 = 0.;
double pk2 = 0.;
double pk3 = 0.;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_linear,
                                    k1_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk1, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  );
//now compute P(k) in units of h^-3 Mpc^3
pk1 *= pow(pba->h,3.);
// double pk2;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_linear,
                                    k2_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk2, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  );
//now compute P(k) in units of h^-3 Mpc^3
pk2 *= pow(pba->h,3.);
// double pk3;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_linear,
                                    k3_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk3, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  );
//now compute P(k) in units of h^-3 Mpc^3
pk3 *= pow(pba->h,3.);


  double f2_eff_12 =  bispectrum_f2_kernel(k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc);
  double f2_eff_23 =  bispectrum_f2_kernel(k2_in_h_over_Mpc,k3_in_h_over_Mpc,k1_in_h_over_Mpc);
  double f2_eff_31 =  bispectrum_f2_kernel(k3_in_h_over_Mpc,k1_in_h_over_Mpc,k2_in_h_over_Mpc);

  // printf("f1 = %.3e \t f2 = %.3e \t f3 = %.3e\n",f2_eff_12,f2_eff_23,f2_eff_31);

  double b123 = 2.*pk1*pk2*f2_eff_12 + 2.*pk2*pk3*f2_eff_23 + 2.*pk3*pk1*f2_eff_31;




  //Evaluation of background quantities @ z:
  double tau;
  int first_index_back = 0;
  double * pvecback;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  class_call(background_tau_of_z(pba,z,&tau),
             ptsz->error_message,
             ptsz->error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &first_index_back,
                               pvecback),
             ptsz->error_message,
             ptsz->error_message);


  double H_over_c_in_h_over_Mpc = pvecback[pba->index_bg_H]/pba->h;
  double galaxy_normalisation = 1.;//H_over_c_in_h_over_Mpc; // here is normally also the bias but we take it out and multiply afterward
  //galaxy_normalisation *= pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,-2.);///(1.+z);
  double sigmaT_over_mp = 8.305907197761162e-17 * pow(pba->h,2)/pba->h; // !this is sigmaT / m_prot in (Mpc/h)**2/(Msun/h)
  double a = 1. / (1. + z);
  double tau_fac = a*sigmaT_over_mp*ptsz->f_b_gas/ptsz->mu_e*ptsz->f_free;
  double rho_crit  = (3./(8.*_PI_*_G_*_M_sun_))
                                  *pow(_Mpc_over_m_,1)
                                  *pow(_c_,2)
                                  *pvecback[pba->index_bg_rho_crit]
                                  /pow(pba->h,2);

  double tau_normalisation =tau_fac
                            *pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,-2.)
                            // *pow((pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0,1);
                            *pvecback[pba->index_bg_Omega_m]*rho_crit;

  free(pvecback);
  // tau_normalisation = 1.;
  double bg = 1.;//get_mean_galaxy_bias_at_z(z,ptsz);
  double vrms2 = get_vrms2_at_z(z,ptsz);
  result = bg*tau_normalisation*tau_normalisation*b123*vrms2;
  return result;

}




int evaluate_ttg_bispectrum_at_z_effective_approach( double * r,
                                                     double k1,
                                                     double k2,
                                                     double k3,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm){

*r =  get_ttg_bispectrum_at_z_effective_approach(k1,k2,k3,z,ptsz,pba,pnl,ppm);
return _SUCCESS_;
}

int evaluate_ttg_bispectrum_at_z_tree_level_PT( double * r,
                                                     double k1,
                                                     double k2,
                                                     double k3,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm){

*r =  get_ttg_bispectrum_at_z_tree_level_PT(k1,k2,k3,z,ptsz,pba,pnl,ppm);
return _SUCCESS_;
}

double get_ttg_bispectrum_at_z_effective_approach(double k1_in_h_over_Mpc,
                                                     double k2_in_h_over_Mpc,
                                                     double k3_in_h_over_Mpc,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm){

double * pk_ic = NULL;
double result = 0.;
double pk1 = 0.;
double pk2 = 0.;
double pk3 = 0.;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_nonlinear,
                                    k1_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk1, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  );
//now compute P(k) in units of h^-3 Mpc^3
pk1 *= pow(pba->h,3.);
// if (isnan(pk1) || isinf(pk1)) pk1 = 0.;

//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_nonlinear,
                                    k2_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk2, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  );
//now compute P(k) in units of h^-3 Mpc^3
pk2 *= pow(pba->h,3.);
// if (isnan(pk2) || isinf(pk2)) pk2 = 0.;

//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_nonlinear,
                                    k3_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk3, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  );

//now compute P(k) in units of h^-3 Mpc^3
pk3 *= pow(pba->h,3.);
// if (isnan(pk3) || isinf(pk3)) pk3 = 0.;


  double z_asked = log(1.+z);
  double R_asked = log(8./pba->h); //log(R) in Mpc
  double sigma8_at_z =  exp(pwl_interp_2d(ptsz->n_arraySZ,
                             ptsz->ndimSZ,
                             ptsz->array_redshift,
                             ptsz->array_radius,
                             ptsz->array_sigma_at_z_and_R,
                             1,
                             &z_asked,
                             &R_asked));

  double knl = get_knl_at_z(z,ptsz);

  double n1 = get_nl_index_at_z_and_k_no_wiggles(z,k1_in_h_over_Mpc,ptsz,pnl);
  double n2 = get_nl_index_at_z_and_k_no_wiggles(z,k2_in_h_over_Mpc,ptsz,pnl);
  double n3 = get_nl_index_at_z_and_k_no_wiggles(z,k3_in_h_over_Mpc,ptsz,pnl);

  // printf("n1 = %.3e \t n2 = %.3e \t n3 = %.3e\n",n1,n2,n3);

  double f2_eff_12 =  bispectrum_f2_kernel_eff(k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,n1,n2,sigma8_at_z,knl);
  double f2_eff_23 =  bispectrum_f2_kernel_eff(k2_in_h_over_Mpc,k3_in_h_over_Mpc,k1_in_h_over_Mpc,n2,n3,sigma8_at_z,knl);
  double f2_eff_31 =  bispectrum_f2_kernel_eff(k3_in_h_over_Mpc,k1_in_h_over_Mpc,k2_in_h_over_Mpc,n3,n1,sigma8_at_z,knl);

  // printf("f1 = %.3e \t f2 = %.3e \t f3 = %.3e\n",f2_eff_12,f2_eff_23,f2_eff_31);

  double b123 = 2.*pk1*pk2*f2_eff_12 + 2.*pk2*pk3*f2_eff_23 + 2.*pk3*pk1*f2_eff_31;
  // if (isnan(b123) || isinf(b123) || b123 == 0.){
  //   // printf("f1 = %.3e \t f2 = %.3e \t f3 = %.3e\n",f2_eff_12,f2_eff_23,f2_eff_31);
  //   // printf("in func get : k1 = %.3e k2 = %.3e k3 = %.3e b123 = %.8e\n",k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,b123);
  //   //exit(0);
  // }




  //Evaluation of background quantities @ z:
  double tau;
  int first_index_back = 0;
  double * pvecback;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  class_call(background_tau_of_z(pba,z,&tau),
             ptsz->error_message,
             ptsz->error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &first_index_back,
                               pvecback),
             ptsz->error_message,
             ptsz->error_message);


  double H_over_c_in_h_over_Mpc = pvecback[pba->index_bg_H]/pba->h;
  double galaxy_normalisation = 1.;//H_over_c_in_h_over_Mpc; // here is normally also the bias but we take it out and multiply afterward
  //galaxy_normalisation *= pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,-2.);///(1.+z);
  double sigmaT_over_mp = 8.305907197761162e-17 * pow(pba->h,2)/pba->h; // !this is sigmaT / m_prot in (Mpc/h)**2/(Msun/h)
  double a = 1. / (1. + z);
  double tau_fac = a*sigmaT_over_mp*ptsz->f_b_gas/ptsz->mu_e*ptsz->f_free;
  double rho_crit  = (3./(8.*_PI_*_G_*_M_sun_))
                                  *pow(_Mpc_over_m_,1)
                                  *pow(_c_,2)
                                  *pvecback[pba->index_bg_rho_crit]
                                  /pow(pba->h,2);
  double tau_normalisation =tau_fac
                            *pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,-2.)
                            *pvecback[pba->index_bg_Omega_m]*rho_crit;
free(pvecback);
// double bg = 1.;
// if (ptsz->use_bg_at_z_in_ksz2g_eff==1)
//   bg = get_mean_galaxy_bias_at_z(z,ptsz);
double vrms2 = get_vrms2_at_z(z,ptsz);
// bg = 1.;
// vrms2 = 1.;

// if (isnan(tau_normalisation) || isinf(tau_normalisation) || tau_normalisation== 0.){
//   printf("in func get : k1 = %.3e k2 = %.3e k3 = %.3e \n",k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc);
// }
result = tau_normalisation*tau_normalisation*vrms2*b123;
// if (k1_in_h_over_Mpc>1e-2 || k1_in_h_over_Mpc>1e-2){ result = 0.;}
// else {result = bg*tau_normalisation*tau_normalisation*vrms2*b123;}
return result;

}


double get_matter_bispectrum_at_z_effective_approach_SC(double k1_in_h_over_Mpc,
                                                     double lambda2,
                                                     double lambda3,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm){

// // get background quantities at z:
double k2_in_h_over_Mpc = lambda2;//*k1_in_h_over_Mpc;
double k3_in_h_over_Mpc = lambda3;//*k1_in_h_over_Mpc;

double * pk_ic = NULL;
double pk1;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
class_call(nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_nonlinear,
                                    k1_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk1, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  ),
                                  pnl->error_message,
                                  pnl->error_message);
//now compute P(k) in units of h^-3 Mpc^3
pk1 *= pow(pba->h,3.);
double pk2;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
class_call(nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_nonlinear,
                                    k2_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk2, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  ),
                                  pnl->error_message,
                                  pnl->error_message);
//now compute P(k) in units of h^-3 Mpc^3
pk2 *= pow(pba->h,3.);
double pk3;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
class_call(nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_nonlinear,
                                    k3_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk3, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  ),
                                  pnl->error_message,
                                  pnl->error_message);
//now compute P(k) in units of h^-3 Mpc^3
pk3 *= pow(pba->h,3.);



  double sigma8_at_z =  get_sigma8_at_z(z,ptsz,pba);

  double knl = get_knl_at_z(z,ptsz);

  double n1 = get_nl_index_at_z_and_k(z,k1_in_h_over_Mpc,ptsz,pnl);
  double n2 = get_nl_index_at_z_and_k(z,k2_in_h_over_Mpc,ptsz,pnl);
  double n3 = get_nl_index_at_z_and_k(z,k3_in_h_over_Mpc,ptsz,pnl);

  // printf("n1 = %.3e \t n2 = %.3e \t n3 = %.3e\n",n1,n2,n3);

  double f2_eff_12 =  bispectrum_f2_kernel_eff_SC(k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,n1,n2,sigma8_at_z,knl);
  double f2_eff_23 =  bispectrum_f2_kernel_eff_SC(k2_in_h_over_Mpc,k3_in_h_over_Mpc,k1_in_h_over_Mpc,n2,n3,sigma8_at_z,knl);
  double f2_eff_31 =  bispectrum_f2_kernel_eff_SC(k3_in_h_over_Mpc,k1_in_h_over_Mpc,k2_in_h_over_Mpc,n3,n1,sigma8_at_z,knl);

  // printf("f1 = %.3e \t f2 = %.3e \t f3 = %.3e\n",f2_eff_12,f2_eff_23,f2_eff_31);

  double b123 = 2.*pk1*pk2*f2_eff_12 + 2.*pk2*pk3*f2_eff_23 + 2.*pk3*pk1*f2_eff_31;
  return b123;

}




double get_matter_bispectrum_at_z_tree_level_PT(double k1_in_h_over_Mpc,
                                                 double lambda2,
                                                 double lambda3,
                                                 double z,
                                                 struct tszspectrum * ptsz,
                                                 struct background * pba,
                                                 struct nonlinear * pnl,
                                                 struct primordial * ppm){

// // get background quantities at z:
double k2_in_h_over_Mpc = lambda2;//*k1_in_h_over_Mpc;
double k3_in_h_over_Mpc = lambda3;//*k1_in_h_over_Mpc;


double * pk_ic = NULL;
double pk1;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
class_call(nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_linear,
                                    k1_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk1, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  ),
                                  pnl->error_message,
                                  pnl->error_message);
//now compute P(k) in units of h^-3 Mpc^3
pk1 *= pow(pba->h,3.);
double pk2;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
class_call(nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_linear,
                                    k2_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk2, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  ),
                                  pnl->error_message,
                                  pnl->error_message);
//now compute P(k) in units of h^-3 Mpc^3
pk2 *= pow(pba->h,3.);
double pk3;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
class_call(nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_linear,
                                    k3_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk3, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  ),
                                  pnl->error_message,
                                  pnl->error_message);
//now compute P(k) in units of h^-3 Mpc^3
pk3 *= pow(pba->h,3.);





  // double sigma8_at_z =  get_sigma8_at_z(z,ptsz,pba);

  // double knl = get_knl_at_z(z,ptsz);
  //
  // double n1 = get_nl_index_at_z_and_k(z,k1_in_h_over_Mpc,ptsz,pnl);
  // double n2 = get_nl_index_at_z_and_k(z,k2_in_h_over_Mpc,ptsz,pnl);
  // double n3 = get_nl_index_at_z_and_k(z,k3_in_h_over_Mpc,ptsz,pnl);

  // printf("n1 = %.3e \t n2 = %.3e \t n3 = %.3e\n",n1,n2,n3);

  double f2_eff_12 =  bispectrum_f2_kernel(k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc);
  double f2_eff_23 =  bispectrum_f2_kernel(k2_in_h_over_Mpc,k3_in_h_over_Mpc,k1_in_h_over_Mpc);
  double f2_eff_31 =  bispectrum_f2_kernel(k3_in_h_over_Mpc,k1_in_h_over_Mpc,k2_in_h_over_Mpc);

  // printf("f1 = %.3e \t f2 = %.3e \t f3 = %.3e\n",f2_eff_12,f2_eff_23,f2_eff_31);

  double b123 = 2.*pk1*pk2*f2_eff_12 + 2.*pk2*pk3*f2_eff_23 + 2.*pk3*pk1*f2_eff_31;
  return b123;

}


double get_matter_bispectrum_at_z_effective_approach_smoothed(double k1_in_h_over_Mpc,
                                                     double lambda2,
                                                     double lambda3,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm){

double k2_in_h_over_Mpc = lambda2;//*k1_in_h_over_Mpc;
double k3_in_h_over_Mpc = lambda3;//*k1_in_h_over_Mpc;
double * pk_ic = NULL;
double pk1;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
class_call(nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_nonlinear,
                                    k1_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk1, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  ),
                                  pnl->error_message,
                                  pnl->error_message);
//now compute P(k) in units of h^-3 Mpc^3
pk1 *= pow(pba->h,3.);
double pk2;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
class_call(nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_nonlinear,
                                    k2_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk2, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  ),
                                  pnl->error_message,
                                  pnl->error_message);
//now compute P(k) in units of h^-3 Mpc^3
pk2 *= pow(pba->h,3.);
double pk3;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
class_call(nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_nonlinear,
                                    k3_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk3, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  ),
                                  pnl->error_message,
                                  pnl->error_message);
//now compute P(k) in units of h^-3 Mpc^3
pk3 *= pow(pba->h,3.);


  double z_asked = log(1.+z);
  double R_asked = log(8./pba->h); //log(R) in Mpc
  double sigma8_at_z =  exp(pwl_interp_2d(ptsz->n_arraySZ,
                             ptsz->ndimSZ,
                             ptsz->array_redshift,
                             ptsz->array_radius,
                             ptsz->array_sigma_at_z_and_R,
                             1,
                             &z_asked,
                             &R_asked));

  double knl = get_knl_at_z(z,ptsz);

  double n1 = get_nl_index_at_z_and_k_no_wiggles(z,k1_in_h_over_Mpc,ptsz,pnl);
  double n2 = get_nl_index_at_z_and_k_no_wiggles(z,k2_in_h_over_Mpc,ptsz,pnl);
  double n3 = get_nl_index_at_z_and_k_no_wiggles(z,k3_in_h_over_Mpc,ptsz,pnl);

  // printf("n1 = %.3e \t n2 = %.3e \t n3 = %.3e\n",n1,n2,n3);

  double f2_eff_12 =  bispectrum_f2_kernel_eff(k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,n1,n2,sigma8_at_z,knl);
  double f2_eff_23 =  bispectrum_f2_kernel_eff(k2_in_h_over_Mpc,k3_in_h_over_Mpc,k1_in_h_over_Mpc,n2,n3,sigma8_at_z,knl);
  double f2_eff_31 =  bispectrum_f2_kernel_eff(k3_in_h_over_Mpc,k1_in_h_over_Mpc,k2_in_h_over_Mpc,n3,n1,sigma8_at_z,knl);

  // printf("f1 = %.3e \t f2 = %.3e \t f3 = %.3e\n",f2_eff_12,f2_eff_23,f2_eff_31);

  double b123 = 2.*pk1*pk2*f2_eff_12 + 2.*pk2*pk3*f2_eff_23 + 2.*pk3*pk1*f2_eff_31;
  return b123;

}




double get_matter_bispectrum_at_z_effective_approach(double k1_in_h_over_Mpc,
                                                     double lambda2,
                                                     double lambda3,
                                                     double z,
                                                     struct tszspectrum * ptsz,
                                                     struct background * pba,
                                                     struct nonlinear * pnl,
                                                     struct primordial * ppm){


double k2_in_h_over_Mpc = lambda2;//*k1_in_h_over_Mpc;
double k3_in_h_over_Mpc = lambda3;//*k1_in_h_over_Mpc;

double * pk_ic = NULL;
double pk1;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
class_call(nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_nonlinear,
                                    k1_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk1, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  ),
                                  pnl->error_message,
                                  pnl->error_message);
//now compute P(k) in units of h^-3 Mpc^3
pk1 *= pow(pba->h,3.);
double pk2;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
class_call(nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_nonlinear,
                                    k2_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk2, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  ),
                                  pnl->error_message,
                                  pnl->error_message);
//now compute P(k) in units of h^-3 Mpc^3
pk2 *= pow(pba->h,3.);
double pk3;
//Input: wavenumber in 1/Mpc
//Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
// printf("z=%.3e\n",z);
class_call(nonlinear_pk_at_k_and_z(
                                    pba,
                                    ppm,
                                    pnl,
                                    pk_nonlinear,
                                    k3_in_h_over_Mpc*pba->h,
                                    z,
                                    pnl->index_pk_m,
                                    &pk3, // number *out_pk_l
                                    pk_ic // array out_pk_ic_l[index_ic_ic]
                                  ),
                                  pnl->error_message,
                                  pnl->error_message);
//now compute P(k) in units of h^-3 Mpc^3
pk3 *= pow(pba->h,3.);


  double z_asked = log(1.+z);
  double R_asked = log(8./pba->h); //log(R) in Mpc
  double sigma8_at_z =  exp(pwl_interp_2d(ptsz->n_arraySZ,
                             ptsz->ndimSZ,
                             ptsz->array_redshift,
                             ptsz->array_radius,
                             ptsz->array_sigma_at_z_and_R,
                             1,
                             &z_asked,
                             &R_asked));

  double knl = get_knl_at_z(z,ptsz);

  double n1 = get_nl_index_at_z_and_k(z,k1_in_h_over_Mpc,ptsz,pnl);
  double n2 = get_nl_index_at_z_and_k(z,k2_in_h_over_Mpc,ptsz,pnl);
  double n3 = get_nl_index_at_z_and_k(z,k3_in_h_over_Mpc,ptsz,pnl);

  // printf("n1 = %.3e \t n2 = %.3e \t n3 = %.3e\n",n1,n2,n3);

  double f2_eff_12 =  bispectrum_f2_kernel_eff(k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,n1,n2,sigma8_at_z,knl);
  double f2_eff_23 =  bispectrum_f2_kernel_eff(k2_in_h_over_Mpc,k3_in_h_over_Mpc,k1_in_h_over_Mpc,n2,n3,sigma8_at_z,knl);
  double f2_eff_31 =  bispectrum_f2_kernel_eff(k3_in_h_over_Mpc,k1_in_h_over_Mpc,k2_in_h_over_Mpc,n3,n1,sigma8_at_z,knl);

  // printf("f1 = %.3e \t f2 = %.3e \t f3 = %.3e\n",f2_eff_12,f2_eff_23,f2_eff_31);

  double b123 = 2.*pk1*pk2*f2_eff_12 + 2.*pk2*pk3*f2_eff_23 + 2.*pk3*pk1*f2_eff_31;
  return b123;

}

// this is in dimensionless unites:
// vrms2/3/c2
double get_vrms2_at_z(double z,
                      struct tszspectrum * ptsz)
  {

// double z = pvectsz[ptsz->index_z];
double z_asked = log(1.+z);

if (z<exp(ptsz->array_redshift[0])-1.)
  z_asked = ptsz->array_redshift[0];
if (z>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
  z_asked =  ptsz->array_redshift[ptsz->n_arraySZ-1];


double vrms2 =  exp(pwl_value_1d(ptsz->n_arraySZ,
                                 ptsz->array_redshift,
                                 ptsz->array_vrms2_at_z,
                                 z_asked));


return vrms2/3./pow(_c_*1e-3,2.);
}




// evaluate normalisation of galaxy terms, i.e., ng_bar
double evaluate_mean_galaxy_number_density_at_z_ngal(
                   double z,
                   int index_g,
                   struct tszspectrum * ptsz)
  {

   // double z = pvectsz[ptsz->index_z];
   double z_asked = log(1.+z);

   if (z<exp(ptsz->array_redshift[0])-1.)
      z_asked = ptsz->array_redshift[0];
   if (z>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
      z_asked =  ptsz->array_redshift[ptsz->n_arraySZ-1];

// int i = 0;
// for (i = 0 ; i<10; i++){
//   printf("i = %d array_z = %.3e array_nbar = %.5e\n",
//           i, ptsz->array_redshift[i],ptsz->array_mean_galaxy_number_density_ngal[index_g][i]);
// }

    return exp(pwl_value_1d(ptsz->n_arraySZ,
                            ptsz->array_redshift,
                            ptsz->array_mean_galaxy_number_density_ngal[index_g],
                            z_asked));
   //pvectsz[ptsz->index_mean_galaxy_number_density] = 1.; // debug BB

// return _SUCCESS_;
}


// evaluate normalisation of galaxy terms, i.e., ng_bar
double evaluate_mean_galaxy_number_density_at_z(
                   double z,
                   struct tszspectrum * ptsz)
  {

   // double z = pvectsz[ptsz->index_z];
   double z_asked = log(1.+z);

   if (z<exp(ptsz->array_redshift[0])-1.)
      z_asked = ptsz->array_redshift[0];
   if (z>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
      z_asked =  ptsz->array_redshift[ptsz->n_arraySZ-1];

// int i = 0;
// for (i = 0 ; i<10; i++){
//   printf("i = %d array_z = %.3e array_nbar = %.5e\n",
//           i, ptsz->array_redshift[i],ptsz->array_mean_galaxy_number_density[i]);
// }

    return exp(pwl_value_1d(ptsz->n_arraySZ,
                            ptsz->array_redshift,
                            ptsz->array_mean_galaxy_number_density,
                            z_asked));
   //pvectsz[ptsz->index_mean_galaxy_number_density] = 1.; // debug BB

// return _SUCCESS_;
}


double get_mean_galaxy_bias_at_z(
                   double z,
                   struct tszspectrum * ptsz)
  {


   double z_asked = log(1.+z);

   if (z<exp(ptsz->array_redshift[0])-1.)
      return 0.;
   if (z>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
      return 0.;

    return exp(pwl_value_1d(ptsz->n_arraySZ,
                            ptsz->array_redshift,
                            ptsz->array_mean_galaxy_bias,
                            z_asked));

}


int evaluate_sigma2_hsv(double * pvecback,
                   double * pvectsz,
                   struct background * pba,
                   struct nonlinear * pnl,
                   struct tszspectrum * ptsz)
  {

   double z = pvectsz[ptsz->index_z];
   double z_asked = log(1.+z);

   if (z<exp(ptsz->array_redshift[0])-1.)
      z_asked = ptsz->array_redshift[0];
   if (z>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
      z_asked =  ptsz->array_redshift[ptsz->n_arraySZ-1];


   pvectsz[ptsz->index_sigma2_hsv] =  exp(pwl_value_1d(ptsz->n_arraySZ,
                                                        ptsz->array_redshift,
                                                        ptsz->array_sigma2_hsv_at_z,
                                                        z_asked));

return _SUCCESS_;
}

int write_output_to_files_ell_indep_ints(struct nonlinear * pnl,
                                         struct background * pba,
                                         struct tszspectrum * ptsz){

  //int index_l;
  char Filepath[_ARGUMENT_LENGTH_MAX_];

   if (ptsz->sz_verbose > 0)
   {
      FILE *fp;

if (ptsz->has_mean_y){
      sprintf(Filepath,"%s%s%s",ptsz->root,"mean_y",".txt");

      printf("-> Writing output files in %s\n",Filepath);
      fp=fopen(Filepath, "w");
      fprintf(fp,"#Input mass bias b = %e\n",
                  1.-1./ptsz->HSEbias);
      fprintf(fp,"#sigma8 = %e\n",
                  pnl->sigma8);
      fprintf(fp,"#Omega_m = %e\n",
                  ptsz->Omega_m_0);
      fprintf(fp,"#h = %e\n",
                  pba->h);

      fprintf(fp,"#Average Compton y parameter\n");
      fprintf(fp,"%e\n",ptsz->y_monopole);
      printf("->Output written in %s\n",Filepath);

      fclose(fp);

    }


      //FILE *fp;
if (ptsz->has_hmf){
      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "hmf_int",
                  ".txt");


      printf("-> Writing output files in %s\n",Filepath);
      fp=fopen(Filepath, "w");
      fprintf(fp,"#Input mass bias b = %e\n",
                  1.-1./ptsz->HSEbias);
      fprintf(fp,"#sigma8 = %e\n",
                  pnl->sigma8[pnl->index_pk_m]);
      fprintf(fp,"#Omega_m = %e\n",
                  ptsz->Omega_m_0);
      fprintf(fp,"#h = %e\n",
                  pba->h);

      fprintf(fp,"#N_tot per steradian\n");
      fprintf(fp,"%.10e\n",ptsz->hmf_int);
      fprintf(fp,"#N_tot full sky (*= 4*pi)\n");
      fprintf(fp,"%.10e\n",ptsz->hmf_int*4.*_PI_); //3.046174198e-4*41253.0= 4*pi // full-sky (4 x pi x 57.3^2=41253 square degrees where 57.3  = 360 / (2 x pi)) // conversion deg2 to steradian 3.046174198e-4
      fprintf(fp,"#N_tot survey (*= 4*pi*f_sky)\n");
      fprintf(fp,"%.10e\n",ptsz->hmf_int*ptsz->Omega_survey);
      printf("->Output written in %s\n",Filepath);

      fclose(fp);

    }



   }

return _SUCCESS_;
}


int write_redshift_dependent_quantities(struct background * pba,
                                        struct tszspectrum * ptsz){

  char Filepath[_ARGUMENT_LENGTH_MAX_];
  FILE *fp;
  sprintf(Filepath,"%s%s%s",ptsz->root,"","redshift_dependent_functions.txt");
  if (ptsz->sz_verbose>=1)
    printf("-> Writing redshift dependent functions in %s\n",Filepath);
  fp=fopen(Filepath, "w");
  fprintf(fp,"# Column 1: redshift\n");
  fprintf(fp,"# Column 2: scale factor\n");
  fprintf(fp,"# Column 3: Hubble parameter [km/s/Mpc]\n");
  fprintf(fp,"# Column 4: Hubble parameter [Mpc^-1]\n"); //pba->H0 = pba->h * 1.e5 / _c_;
  fprintf(fp,"# Column 5: sigma8\n");
  fprintf(fp,"# Column 6: velocity growth rate f = dlnD/dlna\n");
  fprintf(fp,"# Column 7: vrms2 [km/s]\n");
  fprintf(fp,"# Column 8: Delta Chi sigma2_hsv [Mpc/h,see Eq. 33 of Takada and Spergel 2013]\n");
  fprintf(fp,"# Column 9: mvir  [Msun/h] (typical <M>)\n");
  fprintf(fp,"# Column 10: mean galaxy number density ng_bar\n");
  fprintf(fp,"# Column 11: mdel\n");
  fprintf(fp,"# Column 12: mdel_prime\n");
  fprintf(fp,"# Column 13: delta_prime\n");
  fprintf(fp,"# Column 14: delta_prime\n");
  int index_z;
  int n_z = 1e3;
  double z_min = ptsz->z1SZ;
  double z_max = ptsz->z2SZ;
  // double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,ptsz->z_for_pk_hm);
  // double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,ptsz->z_for_pk_hm);
  //background quantities @ z:
  double tau;
  int first_index_back = 0;
  double * pvecback;

  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              pba->error_message);


  for (index_z=0;index_z<n_z;index_z++){

  double ln1pz =  log(1.+z_min)
                  +index_z*(log(1.+z_max)-log(1.+z_min))
                  /(n_z-1.); // log(1+z)

  double z = exp(ln1pz)-1.;



  class_call(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             pba->error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &first_index_back,
                               pvecback),
             pba->error_message,
             pba->error_message);

  double f = pvecback[pba->index_bg_f];
  double a = pvecback[pba->index_bg_a];
  double H = pvecback[pba->index_bg_H]; //in Mpc^-1
  double H_cosmo = pvecback[pba->index_bg_H]*_c_/1e3; // in km/s/Mpc


  double z_asked = log(1.+z);
  double R_asked = log(8./pba->h); //log(R) in Mpc
  // [sigma8 is defined at at R = 8 Mpc/h]
  double sigma8;
  if (ptsz->need_sigma == 1){

   if (z<exp(ptsz->array_redshift[0])-1.)
      z_asked = ptsz->array_redshift[0];
   if (z>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
      z_asked =  ptsz->array_redshift[ptsz->n_arraySZ-1];
   if (R_asked<ptsz->array_radius[0])
      R_asked = ptsz->array_radius[0];
   if (R_asked>ptsz->array_radius[ptsz->ndimSZ-1])
      R_asked =  ptsz->array_radius[ptsz->ndimSZ-1];

    sigma8 =  exp(pwl_interp_2d(ptsz->n_arraySZ,
                       ptsz->ndimSZ,
                       ptsz->array_redshift,
                       ptsz->array_radius,
                       ptsz->array_sigma_at_z_and_R,
                       1,
                       &z_asked,
                       &R_asked));
                     }

  double vrms2 = 0.;
    if (ptsz->has_vrms2)
   vrms2 =  exp(pwl_value_1d(ptsz->n_arraySZ,
                                    ptsz->array_redshift,
                                    ptsz->array_vrms2_at_z,
                                    z_asked));
  double sigma2_hsv = 0.;
    if (ptsz->has_sigma2_hsv)
  sigma2_hsv = exp(pwl_value_1d(ptsz->n_arraySZ,
                                   ptsz->array_redshift,
                                   ptsz->array_sigma2_hsv_at_z,
                                   z_asked));


  double ng = 0.;
  if (ptsz->has_tSZ_gal_1h
     || ptsz->has_gal_gal_1h
     || ptsz->has_gal_gal_2h
     || ptsz->has_gal_lens_1h
     || ptsz->has_gal_lens_2h
     || ptsz->has_gal_lensmag_1h
     || ptsz->has_gal_lensmag_2h
     || ptsz->has_kSZ_kSZ_gal_1h)
  ng  = exp(pwl_value_1d(ptsz->n_arraySZ,
                          ptsz->array_redshift,
                          ptsz->array_mean_galaxy_number_density,
                          z_asked));


  double mvir_over_m200d = 0.;

  double rhoc =  (3./(8.*_PI_*_G_*_M_sun_))
                *pow(_Mpc_over_m_,1)
                *pow(_c_,2)
                *pvecback[pba->index_bg_rho_crit]
                /pow(pba->h,2);
  double Eh = pvecback[pba->index_bg_H]/pba->H0;
  double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
  double delc = Delta_c_of_Omega_m(omega);



  double mvir = pow(10.,14); //m200m for M=10^{13.5} Msun/h
  double rvir = evaluate_rvir_of_mvir(mvir,delc,rhoc,ptsz);

  double cvir;
    if (ptsz->need_hmf) // some relations require sigma
      cvir = evaluate_cvir_of_mvir(mvir,z,ptsz,pba);// 5.72*pow(mvir/1e14,-0.081)/pow(1.+z_asked,0.71);
    else
      cvir = -1;
  double mdel = 1.;

  ///double rs = rvir/cvir;

  double delrho = 200.*omega*rhoc; // 200m
  double delrho_prime = 500.*rhoc; //500c


  // class_call(mVIR_to_mDEL(mvir,
  //                        rvir,
  //                        cvir,
  //                        delrho,
  //                        &mdel,
  //                        ptsz),
  //                 ptsz->error_message,
  //                 ptsz->error_message);

  mvir_over_m200d = mvir/mdel;

  double mvir_recovered = 1.;
  double mvir_over_mvir_recovered;



  // class_call(mDEL_to_mVIR(mdel,
  //                        delrho,
  //                        delc,
  //                        rhoc,
  //                        z,
  //                        &mvir_recovered,
  //                        ptsz,
  //                        pba),
  //                 ptsz->error_message,
  //                 ptsz->error_message);


  double mdel_prime = 1.;
  // class_call(mDEL_to_mDELprime(mdel,
  //                        delrho,
  //                        delrho_prime,
  //                        delc,
  //                        rhoc,
  //                        z,
  //                        &mdel_prime,
  //                        ptsz,
  //                        pba),
  //                 ptsz->error_message,
  //                 ptsz->error_message);
  mvir_over_mvir_recovered = mvir/mvir_recovered;

  double cvir_fac = 1.;
  double cvir_prime = cvir_fac*cvir;

  double delta= 2.5;
  double delta_prime;

  delta_prime = 1.;//delta_to_delta_prime_nfw(delta,cvir,cvir_prime,ptsz);


  fprintf(fp,"%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e\t %.5e\t %.5e\t %.5e\t %.5e\t %.5e\n",
          z,a,H_cosmo,H,sigma8,f,vrms2,sigma2_hsv,mvir,ng,mdel,mdel_prime,delta,delta_prime);

 }

 fclose(fp);

 //exit(0);
 if (ptsz->has_cib_cib_1h){
 sprintf(Filepath,"%s%s%s",ptsz->root,"","cib.txt");
 printf("-> Writing sed cib in %s\n",Filepath);
 fp=fopen(Filepath, "w");

 int index_nu;
 int n_nu = 100;
 double nu_min = 1e1;
 double nu_max = 1e5;
 for (index_nu=0;index_nu<n_nu;index_nu++){

   double ln1pnu =  log(1.+nu_min)
                   +index_nu*(log(1.+nu_max)-log(1.+nu_min))
                   /(n_nu-1.); // log(1+z)

   double nu = exp(ln1pnu)-1.;
   double z = 5.;
   double theta = evaluate_sed_cib(z,  nu, ptsz);

   fprintf(fp,"%.5e \t %.5e\n",nu,theta);
 }

  fclose(fp);
}



  free(pvecback);



                                        }

int write_output_to_files_cl(struct nonlinear * pnl,
                             struct background * pba,
                             struct primordial * ppm,
                             struct tszspectrum * ptsz){


   int index_l;
   char Filepath[_ARGUMENT_LENGTH_MAX_];
   //
   // if (ptsz->sz_verbose > 0)
   // {

/// fill arrays:


if (ptsz->has_sz_ps
  + ptsz->has_sz_trispec
  + ptsz->has_sz_2halo
  + ptsz->has_sz_m_y_y_2h
  + ptsz->has_sz_m_y_y_1h
  + ptsz->has_sz_te_y_y
  + ptsz->has_kSZ_kSZ_gal_1h
  // + ptsz->has_kSZ_kSZ_gal_2h
  + ptsz->has_tSZ_gal_1h
  + ptsz->has_tSZ_gal_2h
  + ptsz->has_tSZ_lens_1h
  + ptsz->has_tSZ_lens_2h
  + ptsz->has_isw_lens
  + ptsz->has_isw_tsz
  + ptsz->has_isw_auto){

double bin_ell_low[ptsz->nlSZ];
double bin_ell_up[ptsz->nlSZ];

for (index_l=0;index_l<ptsz->nlSZ;index_l++){

   double sig_cl_squared;
   double sig_cl_squared_1h;
   double sig_cl_squared_2h;
   double ell = ptsz->ell[index_l];
   sig_cl_squared = 2.*pow((ptsz->cl_sz_1h[index_l]+ptsz->cl_sz_2h[index_l])/(ell*(ell+1.))*2.*_PI_,2.)/(2.*ell+1.);
   sig_cl_squared_1h = 2.*pow((ptsz->cl_sz_1h[index_l])/(ell*(ell+1.))*2.*_PI_,2.)/(2.*ell+1.);
   sig_cl_squared_2h = 2.*pow((ptsz->cl_sz_2h[index_l])/(ell*(ell+1.))*2.*_PI_,2.)/(2.*ell+1.);


   //binned gaussian variance
   double sig_cl_squared_binned;
   double sig_cl_squared_binned_1h;
   double sig_cl_squared_binned_2h;
   double ln_ell_min;
   double ln_ell_max;
   double ln_ell_down;
   double ln_ell_up;
   double n_modes;
   if (index_l == 0){
      ln_ell_up = log(ptsz->ell[index_l+1]);
      ln_ell_max = log(ell) + 0.5*(ln_ell_up-log(ell));
      ln_ell_min = log(ell) - 0.5*(ln_ell_up-log(ell));
      bin_ell_low[index_l] = exp(ln_ell_min);
      bin_ell_up[index_l] = exp(ln_ell_max);
      n_modes = exp(ln_ell_max)-exp(ln_ell_min);
      sig_cl_squared_binned = sig_cl_squared/n_modes;
      sig_cl_squared_binned_1h = sig_cl_squared_1h/n_modes;
      sig_cl_squared_binned_2h = sig_cl_squared_2h/n_modes;
   }
   else if (index_l == ptsz->nlSZ -1){
      ln_ell_down = log(ptsz->ell[index_l-1]);
      ln_ell_min = log(ell) - 0.5*(log(ell)-ln_ell_down);
      ln_ell_max = log(ell) + 0.5*(log(ell)-ln_ell_down);
      bin_ell_low[index_l] = exp(ln_ell_min);
      bin_ell_up[index_l] = exp(ln_ell_max);
      n_modes = exp(ln_ell_max)-exp(ln_ell_min);
      sig_cl_squared_binned = sig_cl_squared/n_modes;
      sig_cl_squared_binned_1h = sig_cl_squared_1h/n_modes;
      sig_cl_squared_binned_2h = sig_cl_squared_2h/n_modes;
   }
   else {
      ln_ell_down = log(ptsz->ell[index_l-1]);
      ln_ell_up = log(ptsz->ell[index_l+1]);
      ln_ell_min = log(ell) - 0.5*(log(ell)-ln_ell_down);
      ln_ell_max = log(ell) + 0.5*(ln_ell_up-log(ell));
      bin_ell_low[index_l] = exp(ln_ell_min);
      bin_ell_up[index_l] = exp(ln_ell_max);
      n_modes = exp(ln_ell_max)-exp(ln_ell_min);
      sig_cl_squared_binned = sig_cl_squared/n_modes;
      sig_cl_squared_binned_1h = sig_cl_squared_1h/n_modes;
      sig_cl_squared_binned_2h = sig_cl_squared_2h/n_modes;
   }

   // normalised cov:
   // see e.g., A 28 of Hill and Pajer 2013
   ptsz->sig_cl_squared_binned[index_l] = sig_cl_squared_binned/ptsz->f_sky;
   ptsz->cov_cl_cl[index_l] = ptsz->sig_cl_squared_binned[index_l] +  ptsz->tllprime_sz[index_l][index_l]/ptsz->Omega_survey;

   sig_cl_squared_binned_1h = sig_cl_squared_binned_1h/ptsz->f_sky;
   sig_cl_squared_binned_2h = sig_cl_squared_binned_2h/ptsz->f_sky;

   double cov_cl_cl_1h = sig_cl_squared_binned_1h + ptsz->tllprime_sz[index_l][index_l]/ptsz->Omega_survey;
   double cov_cl_cl_2h = sig_cl_squared_binned_2h + ptsz->tllprime_sz[index_l][index_l]/ptsz->Omega_survey;



}

if (ptsz->create_ref_trispectrum_for_cobaya){

  if (ptsz->sz_verbose>=1)
    printf("creating reference trispectrum for cobaya sz likelihood\n");

  char Filepath[_ARGUMENT_LENGTH_MAX_];
  FILE *fp;

  double ell;
  double ell_prime;

  int index_l;
  int index_l_prime;

  for (index_l=0;index_l<ptsz->nlSZ;index_l++)
     for (index_l_prime=0;index_l_prime<index_l+1;index_l_prime++) {

       ell_prime = ptsz->ell[index_l_prime];
       ell = ptsz->ell[index_l];


       ptsz->trispectrum_ref[index_l][index_l_prime] = ell*(ell+1.)/(2.*_PI_)*ell_prime*(ell_prime+1.)/(2.*_PI_)*ptsz->tllprime_sz[index_l][index_l_prime];

       if(ptsz->sz_verbose>2) printf("%e\t%e\t%e\n",ell,ell_prime,ptsz->trispectrum_ref[index_l][index_l_prime]);

       ptsz->trispectrum_ref[index_l_prime][index_l] = ptsz->trispectrum_ref[index_l][index_l_prime];
     };

     sprintf(Filepath,
             "%s%s%s%s",
             ptsz->path_to_ref_trispectrum_for_cobaya,
             "/tSZ_trispectrum_ref_",
             ptsz->append_name_trispectrum_ref,
             ".txt");

     fp=fopen(Filepath, "w");

     for (index_l=0;index_l<ptsz->nlSZ;index_l++){
      for (index_l_prime=0;index_l_prime<ptsz->nlSZ;index_l_prime++) {
// Here we save:
// ell*(ell+1.)/(2.*_PI_)*ell_prime*(ell_prime+1.)/(2.*_PI_)*ptsz->tllprime_sz[index_l][index_l_prime];
// where tllprime_sz is the dimensionless trispectrum (set: units for tSZ spectrum = dimensionless)
// and scaled so that T_ll' ~ (10^6 y)^4, so sqrt(T_ll') ~ 10^12 y^2 ~ C_l in the same convention
// This does not include the factor 1/Omega_survey with Omega_survey = 4*pi*f_sky,
// which needs to be taken into account in the computation of the covmat
// this is done in the likelihood codes (in cobaya and montepython).
           fprintf(fp,"%e\t",ptsz->trispectrum_ref[index_l][index_l_prime]);
        }
        fprintf(fp,"\n");
     }
     fclose(fp);


     sprintf(Filepath,
             "%s%s%s%s",
             ptsz->path_to_ref_trispectrum_for_cobaya,
             "/tSZ_c_ell_ref_",
             ptsz->append_name_trispectrum_ref,
             ".txt");

     fp=fopen(Filepath, "w");
     for (index_l=0;index_l<ptsz->nlSZ;index_l++)
           fprintf(fp,"%e\t %e\t %e \t %e \t %e \t %e\n",
           ptsz->ell[index_l],
           ptsz->cl_sz_1h[index_l],
           ptsz->cl_sz_2h[index_l],
           pow(ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)/2./_PI_,2.)*ptsz->sig_cl_squared_binned[index_l], // this includes the fsky factor
           bin_ell_low[index_l],
           bin_ell_up[index_l]);
     fclose(fp);
}


}

// write arrays to output file:
      FILE *fp;

    if ((ptsz->has_sz_ps
      + ptsz->has_sz_trispec
      + ptsz->has_sz_2halo
      + ptsz->has_sz_m_y_y_2h
      + ptsz->has_sz_m_y_y_1h
      + ptsz->has_sz_te_y_y
      + ptsz->has_tSZ_tSZ_tSZ_1halo
      + ptsz->has_kSZ_kSZ_gal_1h
      + ptsz->has_kSZ_kSZ_gal_1h_fft
      + ptsz->has_kSZ_kSZ_gal_2h_fft
      + ptsz->has_kSZ_kSZ_gal_3h_fft
      + ptsz->has_kSZ_kSZ_gallens_1h_fft
      + ptsz->has_kSZ_kSZ_gallens_2h_fft
      + ptsz->has_kSZ_kSZ_gallens_3h_fft
      + ptsz->has_kSZ_kSZ_gallens_hf
      + ptsz->has_kSZ_kSZ_lens_1h_fft
      + ptsz->has_kSZ_kSZ_lens_2h_fft
      + ptsz->has_kSZ_kSZ_lens_3h_fft
      + ptsz->has_kSZ_kSZ_lens_hf
      + ptsz->has_kSZ_kSZ_tSZ_1h
      + ptsz->has_kSZ_kSZ_tSZ_2h
      + ptsz->has_kSZ_kSZ_1h
      + ptsz->has_kSZ_kSZ_2h
      + ptsz->has_kSZ_kSZ_tSZ_3h
      + ptsz->has_kSZ_kSZ_gal_2h
      + ptsz->has_kSZ_kSZ_gal_3h
      + ptsz->has_kSZ_kSZ_gal_hf
      + ptsz->has_kSZ_kSZ_lensmag_1halo
      + ptsz->has_gal_gal_1h
      + ptsz->has_gal_gal_2h
      + ptsz->has_gal_gal_hf
      + ptsz->has_gal_lens_1h
      + ptsz->has_gal_lens_2h
      + ptsz->has_gal_lens_hf
      + ptsz->has_gal_lensmag_1h
      + ptsz->has_gal_lensmag_2h
      + ptsz->has_gal_gallens_1h
      + ptsz->has_gal_gallens_2h
      + ptsz->has_gallens_gallens_1h
      + ptsz->has_gallens_gallens_2h
      + ptsz->has_gallens_lens_1h
      + ptsz->has_gallens_lens_2h
      + ptsz->has_gal_lensmag_hf
      + ptsz->has_lensmag_lensmag_1h
      + ptsz->has_lensmag_lensmag_2h
      + ptsz->has_lensmag_lensmag_hf
      + ptsz->has_lens_lensmag_1h
      + ptsz->has_lens_lensmag_2h
      + ptsz->has_lens_lensmag_hf
      + ptsz->has_tSZ_gal_1h
      + ptsz->has_tSZ_gal_2h
      + ptsz->has_tSZ_cib_1h
      + ptsz->has_tSZ_cib_2h
      + ptsz->has_lens_cib_1h
      + ptsz->has_lens_cib_2h
      + ptsz->has_gal_cib_1h
      + ptsz->has_gal_cib_2h
      + ptsz->has_cib_cib_1h
      + ptsz->has_cib_cib_2h
      + ptsz->has_lens_lens_1h
      + ptsz->has_lens_lens_2h
      + ptsz->has_tSZ_lens_1h
      + ptsz->has_tSZ_lens_2h
      + ptsz->has_isw_lens
      + ptsz->has_isw_tsz
      + ptsz->has_isw_auto > 0) && (ptsz->write_sz>0)){

      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "szpowerspectrum",
                  ".txt");


      //printf("Writing output files in %s\n",Filepath);

      fp=fopen(Filepath, "w");


      fprintf(fp,"# Input mass bias b = %e\n",
                  1.-1./ptsz->HSEbias);
      fprintf(fp,"# sigma8 = %e\n",
                  pnl->sigma8[pnl->index_pk_m]);
      fprintf(fp,"# Omega_m = %e\n",
                  ptsz->Omega_m_0);
      fprintf(fp,"# h = %e\n",
                  pba->h);
      if (ptsz->has_sz_te_y_y){
        char str[100];
        if(ptsz->temperature_mass_relation == 0)
        strcpy(str,"standard");
        if(ptsz->temperature_mass_relation == 0)
        strcpy(str,"lee et al 2019");

      fprintf(fp,"# T-M relation for SZ temperature : %s\n",str);
                            }

       //Halo mass function
       if (ptsz->MF==3) fprintf(fp,"# HMF: Jenkins et al 2001 @ M180m\n");
       if (ptsz->MF==1) fprintf(fp,"# HMF: Tinker et al 2010 @ M200m\n");
       if (ptsz->MF==4) fprintf(fp,"# HMF: Tinker et al 2008 @ M200m\n");
       if (ptsz->MF==5) fprintf(fp,"# HMF: Tinker et al 2008 @ M500c\n");
       if (ptsz->MF==6) fprintf(fp,"# HMF: Tinker et al 2008 @ M1600m\n");
       if (ptsz->MF==2) fprintf(fp,"# HMF: Bocquet et al 2015 @ M200m\n");
       if (ptsz->MF==7) fprintf(fp,"# HMF: Bocquet et al 2015 @ M500c\n\n");

       //Pressure profile
       if (ptsz->pressure_profile == 0)
          fprintf(fp,"# Pressure Profile:  Planck 2013\n");
       if (ptsz->pressure_profile == 2) {
          fprintf(fp,"# Pressure Profile:  Arnaud et al 2010\n");
       }

       if (ptsz->pressure_profile == 3){
          fprintf(fp,"# Pressure Profile:  Custom. GNFW\n");
          fprintf(fp,"# P0GNFW = %e\n",ptsz->P0GNFW);
          fprintf(fp,"# c500 = %e\n",ptsz->c500);
          fprintf(fp,"# gammaGNFW = %e\n",ptsz->gammaGNFW);
          fprintf(fp,"# alphaGNFW = %e\n",ptsz->alphaGNFW);
          fprintf(fp,"# betaGNFW = %e\n",ptsz->betaGNFW);
       }
       if (ptsz->pressure_profile == 4)
          fprintf(fp,"# Pressure Profile:  Battaglia et al 2012\n");


       if (ptsz->exponent_unit == 2.)
          fprintf(fp,"# Dimensions for y power: 'dimensionless'\n");

       if (ptsz->exponent_unit == 0.)
         fprintf(fp,"# Dimensions for y-distortion: 'muK' (micro Kelvin)\n");

       fprintf(fp,"# Omega_survey = %.5e deg2\n", ptsz->Omega_survey*pow(57.3,2));
       fprintf(fp,"# f_sky = %f\n", ptsz->f_sky);
       fprintf(fp,"\n");
       fprintf(fp,"# Columns:\n");
       fprintf(fp,"# 1:multipole\n");
       if (ptsz->exponent_unit == 2.)
          fprintf(fp,"# 2:10^12*ell*(ell+1)/(2*pi)*C_l^tSZ (1-halo term)\n");
       if (ptsz->exponent_unit == 0.)
          fprintf(fp,"# 2:ell*(ell+1)/(2*pi)*C_l^tSZ (1-halo term, [muK^2])\n");
       fprintf(fp,"# 3:Unbinned Gaussian sampling variance (sigma_g_C_l^2), does NOT include f_sky factor\n");
       fprintf(fp,"# 4:Diagonal elements of non-Gaussian sampling variance (trispectrum, T_ll/Omega_survey, i.e.,includes f_sky factor) \n");
       fprintf(fp,"# 5:Binned Gaussian sampling variance (sigma_g_C_l^2_binned, includes f_sky factor)\n");
       fprintf(fp,"# 6:Binned total std dev (Gaussian + non-Gaussian, includes f_sky factor)\n");
        if (ptsz->exponent_unit == 2.)
          fprintf(fp,"# 7:2-halo term 10^12*ell*(ell+1)/(2*pi)*C_l^tSZ (2-halo term)\n");
        if (ptsz->exponent_unit == 0.)
          fprintf(fp,"# 7:2-halo term ell*(ell+1)/(2*pi)*C_l^tSZ (2-halo term, [muK^2])\n");
       fprintf(fp,"# 8:SZ temperature, Te [in keV]\n");
       fprintf(fp,"# 9:cl_kSZ_kSZ_gal_1h [TBC]\n");
       fprintf(fp,"# 10:ell^2*(ell+1)/(2*pi)*C_l^y-phi (1-halo term) in 10^-6 y-units (dimensionless)\n");
       fprintf(fp,"# 11:ell^2*(ell+1)/(2*pi)*C_l^y-phi (2-halo term) in 10^-6 y-units (dimensionless)\n");
       fprintf(fp,"# 12:ell*(ell+1)/(2*pi)*C_l^isw-phi [dimensionless]\n");
       fprintf(fp,"# 13:ell*(ell+1)/(2*pi)*C_l^isw-y [1e-6  dimensionless] \n");
       fprintf(fp,"# 14:ell*(ell+1)/(2*pi)*C_l^isw-isw [dimensionless]\n");
       fprintf(fp,"# 15: signal-to-noise for cl_yy_tot\n");
       fprintf(fp,"# 16: mean mass m_y_y_1h\n");
       fprintf(fp,"# 17: mean mass m_y_y_2h\n");
       fprintf(fp,"# 18: signal-to-noise for cl_yy_1h\n");
       fprintf(fp,"# 19: signal-to-noise for cl_yy_2h\n");
       fprintf(fp,"# 20: sqrt[cov(Y,Y)^ssc/Trispectrum]\n");
       fprintf(fp,"# 21: ell*(ell+1)/(2*pi)*C_l^y-gal (1-halo term) in 10^-6 y-units (dimensionless)\n");
       fprintf(fp,"# 22: b_tSZ_tSZ_tSZ_1halo [TBC]\n");
       fprintf(fp,"# 23: ell*(ell+1)/(2*pi)*C_l^gal-gal (1-halo term) (dimensionless)\n");
       fprintf(fp,"# 24: ell*(ell+1)/(2*pi)*C_l^gal-gal (2-halo term) (dimensionless)\n");
       fprintf(fp,"# 25: ell*(ell+1)/(2*pi)*C_l^gal-lens (1-halo term) (dimensionless)\n");
       fprintf(fp,"# 26: ell*(ell+1)/(2*pi)*C_l^gal-lens (2-halo term) (dimensionless)\n");
       fprintf(fp,"# 27: ell*(ell+1)/(2*pi)*C_l^y-gal (2-halo term) in 10^-6 y-units (dimensionless)\n");
       fprintf(fp,"# 28: ell*(ell+1)/(2*pi)*C_l^lens-lens (1-halo term) (dimensionless)\n");
       fprintf(fp,"# 29: ell*(ell+1)/(2*pi)*C_l^lens-lens (2-halo term) (dimensionless)\n");
       fprintf(fp,"# 30: ell*(ell+1)/(2*pi)*C_l^tSZ-cib (1-halo term) (dim tbc)\n");
       fprintf(fp,"# 31: ell*(ell+1)/(2*pi)*C_l^tSZ-cib (2-halo term) (dim tbc)\n");
       fprintf(fp,"# 32: ell*(ell+1)/(2*pi)*C_l^cib-cib (1-halo term) (dim tbc)\n");
       fprintf(fp,"# 33: ell*(ell+1)/(2*pi)*C_l^cib-cib (2-halo term) (dim tbc)\n");
       fprintf(fp,"# 34: cl_kSZ_kSZ_lensmag_1h (1-halo term) (dim tbc)\n");
       fprintf(fp,"# 35: ell*(ell+1)/(2*pi)*C_l^lens-cib (1-halo term) (dim tbc)\n");
       fprintf(fp,"# 36: ell*(ell+1)/(2*pi)*C_l^lens-cib (2-halo term) (dim tbc)\n");
       fprintf(fp,"# 37: ell*(ell+1)/(2*pi)*C_l^gal-lensmag (1-halo term) (dimensionless)\n");
       fprintf(fp,"# 38: ell*(ell+1)/(2*pi)*C_l^gal-lensmag (2-halo term) (dimensionless)\n");
       fprintf(fp,"# 39: ell*(ell+1)/(2*pi)*C_l^lensmag-lensmag (1-halo term) (dimensionless)\n");
       fprintf(fp,"# 40: ell*(ell+1)/(2*pi)*C_l^lensmag-lensmag (2-halo term) (dimensionless)\n");
       fprintf(fp,"# 41: ell*(ell+1)/(2*pi)*C_l^lens-lensmag (1-halo term) (dimensionless)\n");
       fprintf(fp,"# 42: ell*(ell+1)/(2*pi)*C_l^lens-lensmag (2-halo term) (dimensionless)\n");
       fprintf(fp,"# 43: cl_kSZ_kSZ_gal_2h [TBC]\n");
       fprintf(fp,"# 44: cl_kSZ_kSZ_gal_3h [TBC]\n");
       fprintf(fp,"# 45: C_l^tsz-lensmag 1h [TBC]\n");
       fprintf(fp,"# 46: C_l^tsz-lensmag 2h [TBC]\n");
       fprintf(fp,"# 47: cl_kSZ_kSZ_gal_hf [TBC]\n");
       fprintf(fp,"# 48: ell*(ell+1)/(2*pi)*C_l^gal-gal (effective approach)\n");
       fprintf(fp,"# 49: b_kSZ_kSZ_tSZ_1h \n");
       fprintf(fp,"# 50: b_kSZ_kSZ_tSZ_2h \n");
       fprintf(fp,"# 51: b_kSZ_kSZ_tSZ_3h \n");
       fprintf(fp,"\n");

      for (index_l=0;index_l<ptsz->nlSZ;index_l++){

         double sig_cl_squared;
         double sig_cl_squared_1h;
         double sig_cl_squared_2h;
         double ell = ptsz->ell[index_l];
         sig_cl_squared = 2.*pow((ptsz->cl_sz_1h[index_l]+ptsz->cl_sz_2h[index_l])/(ell*(ell+1.))*2.*_PI_,2.)/(2.*ell+1.);
         sig_cl_squared_1h = 2.*pow((ptsz->cl_sz_1h[index_l])/(ell*(ell+1.))*2.*_PI_,2.)/(2.*ell+1.);
         sig_cl_squared_2h = 2.*pow((ptsz->cl_sz_2h[index_l])/(ell*(ell+1.))*2.*_PI_,2.)/(2.*ell+1.);

         double sig_nl_yy_squared;
         if (ptsz->include_noise_cov_y_y == 1){
         //sig_nl_yy_squared = 0.;
         double nl_yy  =  pwl_value_1d(ptsz->unbinned_nl_yy_size,
                                            ptsz->unbinned_nl_yy_ell,
                                            ptsz->unbinned_nl_yy_n_ell,
                                            ell);

         if (ptsz->nl_yy_is_binned == 1){
           sig_nl_yy_squared = pow(nl_yy/(ell*(ell+1.))*2.*_PI_,2.);
           printf("%d %.3e %.3e\n",ptsz->unbinned_nl_yy_size,ell,nl_yy);
         }
         else {
         sig_nl_yy_squared = 2.*pow(nl_yy/(ell*(ell+1.))*2.*_PI_,2.)/(2.*ell+1.);

          }

         }
         else{
         sig_nl_yy_squared = 0.;
         }

         //binned gaussian variance
         double sig_nl_yy_squared_binned;

         double sig_cl_squared_binned;
         double sig_cl_squared_binned_1h;
         double sig_cl_squared_binned_2h;
         double ln_ell_min;
         double ln_ell_max;
         double ln_ell_down;
         double ln_ell_up;
         double n_modes;




         if (ptsz->nl_yy_is_binned == 0) {
         if (index_l == 0){
            ln_ell_up = log(ptsz->ell[index_l+1]);
            ln_ell_max = log(ell) + 0.5*(ln_ell_up-log(ell));
            ln_ell_min = log(ell) - 0.5*(ln_ell_up-log(ell));
            n_modes = exp(ln_ell_max)-exp(ln_ell_min);
            sig_cl_squared_binned = sig_cl_squared/n_modes;
            sig_cl_squared_binned_1h = sig_cl_squared_1h/n_modes;
            sig_cl_squared_binned_2h = sig_cl_squared_2h/n_modes;
            sig_nl_yy_squared_binned = sig_nl_yy_squared/n_modes;
         }
         else if (index_l == ptsz->nlSZ -1){
            ln_ell_down = log(ptsz->ell[index_l-1]);
            ln_ell_min = log(ell) - 0.5*(log(ell)-ln_ell_down);
            ln_ell_max = log(ell) + 0.5*(log(ell)-ln_ell_down);
            n_modes = exp(ln_ell_max)-exp(ln_ell_min);
            sig_cl_squared_binned = sig_cl_squared/n_modes;
            sig_cl_squared_binned_1h = sig_cl_squared_1h/n_modes;
            sig_cl_squared_binned_2h = sig_cl_squared_2h/n_modes;
            sig_nl_yy_squared_binned = sig_nl_yy_squared/n_modes;
         }
         else {
            ln_ell_down = log(ptsz->ell[index_l-1]);
            ln_ell_up = log(ptsz->ell[index_l+1]);
            ln_ell_min = log(ell) - 0.5*(log(ell)-ln_ell_down);
            ln_ell_max = log(ell) + 0.5*(ln_ell_up-log(ell));
            n_modes = exp(ln_ell_max)-exp(ln_ell_min);
            sig_cl_squared_binned = sig_cl_squared/n_modes;
            sig_cl_squared_binned_1h = sig_cl_squared_1h/n_modes;
            sig_cl_squared_binned_2h = sig_cl_squared_2h/n_modes;
            sig_nl_yy_squared_binned = sig_nl_yy_squared/n_modes;
         }

         //normalised cov:
         // see e.g., A 28 of Hill and Pajer 2013
         ptsz->sig_cl_squared_binned[index_l] = (sig_cl_squared_binned
                                                +sig_nl_yy_squared_binned)
                                                +2.*sqrt(sig_cl_squared_binned*sig_nl_yy_squared_binned)
                                                /ptsz->f_sky;

         }
         else {
           ptsz->sig_cl_squared_binned[index_l] = sig_nl_yy_squared;

         }
         ptsz->cov_cl_cl[index_l] = ptsz->sig_cl_squared_binned[index_l] +  ptsz->tllprime_sz[index_l][index_l]/ptsz->Omega_survey;

         sig_cl_squared_binned_1h = sig_cl_squared_binned_1h/ptsz->f_sky;
         sig_cl_squared_binned_2h = sig_cl_squared_binned_2h/ptsz->f_sky;

         double cov_cl_cl_1h = sig_cl_squared_binned_1h + ptsz->tllprime_sz[index_l][index_l]/ptsz->Omega_survey;
         double cov_cl_cl_2h = sig_cl_squared_binned_2h + ptsz->tllprime_sz[index_l][index_l]/ptsz->Omega_survey;


            fprintf(fp,
                    "%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n",
                    ptsz->ell[index_l],
                    ptsz->cl_sz_1h[index_l],
                    sig_cl_squared,
                    ptsz->tllprime_sz[index_l][index_l]/ptsz->Omega_survey,
                    ptsz->sig_cl_squared_binned[index_l],
                    ell*(ell+1.)/(2.*_PI_)*sqrt(ptsz->cov_cl_cl[index_l]),
                    ptsz->cl_sz_2h[index_l],
                    ptsz->cl_te_y_y[index_l]/ptsz->cl_sz_1h[index_l],
                    ptsz->cl_kSZ_kSZ_gal_1h[index_l],
                    ptsz->cl_tSZ_lens_1h[index_l],
                    ptsz->cl_tSZ_lens_2h[index_l],
                    ptsz->cl_isw_lens[index_l],
                    ptsz->cl_isw_tsz[index_l],
                    ptsz->cl_isw_auto[index_l],
                    1./(ell*(ell+1.)/(2.*_PI_)*sqrt(ptsz->cov_cl_cl[index_l])/(ptsz->cl_sz_1h[index_l]+ptsz->cl_sz_2h[index_l])),
                    ptsz->m_y_y_1h[index_l]/ptsz->cl_sz_1h[index_l],
                    sqrt(ptsz->m_y_y_2h[index_l]/ptsz->cl_sz_2h[index_l]),
                    1./(ell*(ell+1.)/(2.*_PI_)*sqrt(cov_cl_cl_1h)/(ptsz->cl_sz_1h[index_l])),
                    1./(ell*(ell+1.)/(2.*_PI_)*sqrt(cov_cl_cl_2h)/(ptsz->cl_sz_2h[index_l])),
                    sqrt(ptsz->cov_Y_Y_ssc[index_l][index_l]/(ptsz->tllprime_sz[index_l][index_l]/ptsz->Omega_survey)),
                    ptsz->cl_tSZ_gal_1h[index_l],
                    ptsz->b_tSZ_tSZ_tSZ_1halo[index_l],
                    ptsz->cl_gal_gal_1h[index_l],
                    ptsz->cl_gal_gal_2h[index_l],
                    ptsz->cl_gal_lens_1h[index_l],
                    ptsz->cl_gal_lens_2h[index_l],
                    ptsz->cl_tSZ_gal_2h[index_l],
                    ptsz->cl_lens_lens_1h[index_l],
                    ptsz->cl_lens_lens_2h[index_l],
                    ptsz->cl_tSZ_cib_1h[ptsz->id_nu_cib_to_save][index_l],
                    ptsz->cl_tSZ_cib_2h[ptsz->id_nu_cib_to_save][index_l],
                    ptsz->cl_cib_cib_1h[ptsz->id_nu_cib_to_save][ptsz->id_nu_prime_cib_to_save][index_l],
                    ptsz->cl_cib_cib_2h[ptsz->id_nu_cib_to_save][ptsz->id_nu_prime_cib_to_save][index_l],
                    ptsz->cl_kSZ_kSZ_lensmag_1h[index_l],
                    ptsz->cl_lens_cib_1h[ptsz->id_nu_cib_to_save][index_l],
                    ptsz->cl_lens_cib_2h[ptsz->id_nu_cib_to_save][index_l],
                    ptsz->cl_gal_lensmag_1h[index_l],
                    ptsz->cl_gal_lensmag_2h[index_l],
                    ptsz->cl_lensmag_lensmag_1h[index_l],
                    ptsz->cl_lensmag_lensmag_2h[index_l],
                    ptsz->cl_lens_lensmag_1h[index_l],
                    ptsz->cl_lens_lensmag_2h[index_l],
                    ptsz->cl_kSZ_kSZ_gal_2h[index_l],
                    ptsz->cl_kSZ_kSZ_gal_3h[index_l],
                    ptsz->cl_tSZ_lensmag_1h[index_l],
                    ptsz->cl_tSZ_lensmag_2h[index_l],
                    ptsz->cl_kSZ_kSZ_gal_hf[index_l],
                    ptsz->cl_gal_gal_hf[index_l],
                    ptsz->cl_gal_lens_hf[index_l],
                    ptsz->cl_gal_lensmag_hf[index_l],
                    ptsz->cl_lens_lensmag_hf[index_l],
                    ptsz->cl_lensmag_lensmag_hf[index_l],
                    ptsz->cl_gal_cib_1h[ptsz->id_nu_cib_to_save][index_l],
                    ptsz->cl_gal_cib_2h[ptsz->id_nu_cib_to_save][index_l],
                    ptsz->cl_kSZ_kSZ_gal_1h_fft[index_l],
                    ptsz->cl_kSZ_kSZ_gal_2h_fft[index_l],
                    ptsz->cl_kSZ_kSZ_gal_3h_fft[index_l],
                    ptsz->cl_gal_gallens_1h[index_l],
                    ptsz->cl_gal_gallens_2h[index_l],
                    ptsz->b_kSZ_kSZ_tSZ_1h[index_l],
                    ptsz->b_kSZ_kSZ_tSZ_2h[index_l],
                    ptsz->b_kSZ_kSZ_tSZ_3h[index_l],
                    ptsz->cov_ll_kSZ_kSZ_gal[index_l],
                    ptsz->cl_kSZ_kSZ_1h[index_l],
                    ptsz->cl_kSZ_kSZ_2h[index_l],
                    ptsz->cl_kSZ_kSZ_gal_lensing_term[index_l],
                    ptsz->cl_kSZ_kSZ_gallens_1h_fft[index_l],
                    ptsz->cl_kSZ_kSZ_gallens_2h_fft[index_l],
                    ptsz->cl_kSZ_kSZ_gallens_3h_fft[index_l],
                    ptsz->cl_kSZ_kSZ_gallens_lensing_term[index_l],
                    ptsz->cl_kSZ_kSZ_gallens_hf[index_l],
                    ptsz->cl_gallens_gallens_1h[index_l],
                    ptsz->cl_gallens_gallens_2h[index_l],
                    ptsz->cl_gallens_lens_1h[index_l],
                    ptsz->cl_gallens_lens_2h[index_l],
                    ptsz->cov_ll_kSZ_kSZ_gallens[index_l],
                    ptsz->cl_kSZ_kSZ_lens_1h_fft[index_l],
                    ptsz->cl_kSZ_kSZ_lens_2h_fft[index_l],
                    ptsz->cl_kSZ_kSZ_lens_3h_fft[index_l],
                    ptsz->cl_kSZ_kSZ_lens_lensing_term[index_l],
                    ptsz->cl_kSZ_kSZ_lens_hf[index_l],
                    ptsz->cov_ll_kSZ_kSZ_lens[index_l],
                    ptsz->cib_shotnoise[ptsz->id_nu_cib_to_save]
                    );

      }
      if (ptsz->sz_verbose>=1)
        printf("->Output written in %s\n",Filepath);
      fclose(fp);



    }
      int index_l_prime;
      double ell_prime;
      double ell;
      for (index_l=0;index_l<ptsz->nlSZ;index_l++)
         for (index_l_prime=0;index_l_prime<index_l+1;index_l_prime++) {
           ptsz->r_cl_clp[index_l][index_l_prime] = (ptsz->tllprime_sz[index_l][index_l_prime]/ptsz->Omega_survey + ptsz->cov_Y_Y_ssc[index_l][index_l_prime])
                                                    /sqrt(ptsz->cov_cl_cl[index_l] + ptsz->cov_Y_Y_ssc[index_l][index_l])
                                                    /sqrt(ptsz->cov_cl_cl[index_l_prime]+ ptsz->cov_Y_Y_ssc[index_l_prime][index_l_prime]);

           ptsz->r_cl_clp[index_l_prime][index_l] = ptsz->r_cl_clp[index_l][index_l_prime];

         };


if ((ptsz->has_cib_monopole
  >= _TRUE_) && ptsz->write_sz>0){
      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "cib_monopole",
                  ".txt");

      fp=fopen(Filepath, "w");

      int index_k;
      for (index_k=0;index_k<ptsz->n_frequencies_for_cib;index_k++){
      fprintf(fp,"%.4e\t%.4e\n",
      ptsz->frequencies_for_cib[index_k],
      ptsz->cib_monopole[index_k]
    );
    }
      fclose(fp);
  }


if ((ptsz->has_pk_at_z_1h
  +  ptsz->has_pk_at_z_2h
  +  ptsz->has_pk_gg_at_z_1h
  +  ptsz->has_pk_gg_at_z_2h
  +  ptsz->has_pk_bb_at_z_1h
  +  ptsz->has_pk_bb_at_z_2h
  // +  ptsz->has_pk_bb_at_z_2h
  +  ptsz->has_pk_em_at_z_1h
  +  ptsz->has_pk_em_at_z_2h
  +  ptsz->has_pk_HI_at_z_1h
  +  ptsz->has_pk_HI_at_z_2h
  >= _TRUE_) && ptsz->write_sz>0){
      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "halo_model_pk_at_z_k_pk1h_pk2h",
                  ".txt");

      fp=fopen(Filepath, "w");

      int index_k;
      for (index_k=0;index_k<ptsz->n_k_for_pk_hm;index_k++){
      fprintf(fp,"%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n",
      ptsz->k_for_pk_hm[index_k],
      ptsz->pk_at_z_1h[index_k],
      ptsz->pk_at_z_2h[index_k],
      ptsz->pk_gg_at_z_1h[index_k],
      ptsz->pk_gg_at_z_2h[index_k],
      ptsz->pk_bb_at_z_1h[index_k],
      ptsz->pk_bb_at_z_2h[index_k],
      ptsz->pk_HI_at_z_1h[index_k],
      ptsz->pk_HI_at_z_2h[index_k],
      ptsz->pk_em_at_z_1h[index_k],
      ptsz->pk_em_at_z_2h[index_k]
    );
    }
      fclose(fp);
  }


if ((ptsz->has_bk_at_z_1h + ptsz->has_bk_at_z_2h + ptsz->has_bk_at_z_3h >= _TRUE_) && ptsz->write_sz>0){
      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "halo_model_bk_at_z_k_bk1h_bk2h_bk3h",
                  ".txt");

      fp=fopen(Filepath, "w");

      int index_k;
      for (index_k=0;index_k<ptsz->n_k_for_pk_hm;index_k++){
      fprintf(fp,"%.4e\t%.4e\t%.4e\t%.4e\n",ptsz->k_for_pk_hm[index_k],
                                            ptsz->bk_at_z_1h[index_k],
                                            ptsz->bk_at_z_2h[index_k],
                                            ptsz->bk_at_z_3h[index_k]);
  // double k = ptsz->k_for_pk_hm[index_k];
  // double z = ptsz->z_for_pk_hm;
  // double b_hm = ptsz->bk_at_z_3h[index_k];
  // double b_tree = get_matter_bispectrum_at_z_tree_level_PT(k,k,k,z,ptsz,pba,pnl,ppm);
  // printf("bispectrum fields z = %.3e k = %.8e b_hm = %.8e b_tree = %.8e\n",z,k,b_hm,b_tree);

    }
      fclose(fp);
  }


if ((ptsz->has_bk_ttg_at_z_1h + ptsz->has_bk_ttg_at_z_2h + ptsz->has_bk_ttg_at_z_3h >= _TRUE_) && ptsz->write_sz>0){
      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "halo_model_bk_ttg_at_z_k_bk1h_bk2h_bk3h",
                  ".txt");

      fp=fopen(Filepath, "w");

      int index_k;
      for (index_k=0;index_k<ptsz->n_k_for_pk_hm;index_k++){
      fprintf(fp,"%.4e\t%.4e\t%.4e\t%.4e\n",ptsz->k_for_pk_hm[index_k],
                                            ptsz->bk_ttg_at_z_1h[index_k],
                                            ptsz->bk_ttg_at_z_2h[index_k],
                                            ptsz->bk_ttg_at_z_3h[index_k]);
  // double k = ptsz->k_for_pk_hm[index_k];
  // double z = ptsz->z_for_pk_hm;
  // double b_hm = ptsz->bk_at_z_3h[index_k];
  // double b_tree = get_matter_bispectrum_at_z_tree_level_PT(k,k,k,z,ptsz,pba,pnl,ppm);
  // printf("bispectrum fields z = %.3e k = %.8e b_hm = %.8e b_tree = %.8e\n",z,k,b_hm,b_tree);

    }
      fclose(fp);
  }






if (ptsz->has_sz_cov_Y_N && ptsz->write_sz>0){
      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "szpowerspectrum_cov_Y_N",
                  ".txt");

      fp=fopen(Filepath, "w");

      int index_M_bins;
      for (index_l=0;index_l<ptsz->nlSZ;index_l++){
         for (index_M_bins=0;index_M_bins<ptsz->nbins_M;index_M_bins++){
            fprintf(fp,"%e\t",ptsz->cov_Y_N[index_l][index_M_bins]);
      }
         fprintf(fp,"\n");
      }
      fclose(fp);

      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "szpowerspectrum_r_Y_N",
                  ".txt");

      fp=fopen(Filepath, "w");


      for (index_l=0;index_l<ptsz->nlSZ;index_l++){
         for (index_M_bins=0;index_M_bins<ptsz->nbins_M;index_M_bins++){

           ptsz->r_Y_N[index_l][index_M_bins] = (ptsz->cov_Y_N[index_l][index_M_bins]+ptsz->cov_Y_N_next_order[index_l][index_M_bins])
                                                /sqrt(ptsz->cov_N_N[index_M_bins]+ptsz->cov_N_N_hsv[index_M_bins][index_M_bins])
                                                /sqrt(ptsz->cov_cl_cl[index_l]+ptsz->cov_Y_Y_ssc[index_l][index_l]);

            fprintf(fp,"%e\t",ptsz->r_Y_N[index_l][index_M_bins]);
         }
         fprintf(fp,"\n");
      }
      fclose(fp);
  }

  if (ptsz->has_sz_cov_N_N && ptsz->write_sz>0){

    sprintf(Filepath,
              "%s%s%s",
              ptsz->root,
              "szpowerspectrum_cov_N_N_diagonal",
              ".txt");

    fp=fopen(Filepath, "w");

    int index_M_bins;
    for (index_M_bins=0;index_M_bins<ptsz->nbins_M;index_M_bins++)
    fprintf(fp,"%e\t%e\t%e\t%e\n",ptsz->cov_Y_N_mass_bin_edges[index_M_bins],
    ptsz->cov_Y_N_mass_bin_edges[index_M_bins+1],
    ptsz->cov_N_N[index_M_bins],
    ptsz->cov_N_N_hsv[index_M_bins][index_M_bins]);

    fclose(fp);

    sprintf(Filepath,
              "%s%s%s",
              ptsz->root,
              "szpowerspectrum_r_N_N",
              ".txt");

    fp=fopen(Filepath, "w");

    int index_M_bins_1;
    int index_M_bins_2;
    for (index_M_bins_1=0;index_M_bins_1<ptsz->nbins_M;index_M_bins_1++){
    for (index_M_bins_2=0;index_M_bins_2<ptsz->nbins_M;index_M_bins_2++) {
      if (index_M_bins_1 == index_M_bins_2)
      fprintf(fp,"%e\t",1.);
      else
      fprintf(fp,"%e\t",ptsz->cov_N_N_hsv[index_M_bins_1][index_M_bins_2]/sqrt(ptsz->cov_N_N[index_M_bins_1]+ptsz->cov_N_N_hsv[index_M_bins_1][index_M_bins_1])/sqrt(ptsz->cov_N_N[index_M_bins_2]+ptsz->cov_N_N_hsv[index_M_bins_2][index_M_bins_2]));
      }
      fprintf(fp,"\n");
   }
    fclose(fp);


    }
    if (ptsz->has_dndlnM && ptsz->write_sz>0){

      sprintf(Filepath,
                "%s%s%s",
                ptsz->root,
                "szpowerspectrum_dndlnM_masses",
                ".txt");

      fp=fopen(Filepath, "w");

      int index_M_bins;
      for (index_M_bins=0;index_M_bins<ptsz->N_mass_dndlnM;index_M_bins++)
      fprintf(fp,"%e\n",ptsz->dndlnM_array_m[index_M_bins]);

      fclose(fp);

      sprintf(Filepath,
                "%s%s%s",
                ptsz->root,
                "szpowerspectrum_dndlnM_redshifts",
                ".txt");

      fp=fopen(Filepath, "w");

      int index_z_bins;
      for (index_z_bins=0;index_z_bins<ptsz->N_redshift_dndlnM;index_z_bins++)
      fprintf(fp,"%e\n",ptsz->dndlnM_array_z[index_z_bins]);

      fclose(fp);

      sprintf(Filepath,
                "%s%s%s",
                ptsz->root,
                "szpowerspectrum_dndlnM",
                ".txt");

      fp=fopen(Filepath, "w");

      int index_masses;
      int index_redshifts;
      for (index_masses=0;index_masses<ptsz->N_mass_dndlnM;index_masses++) {
      for (index_redshifts=0;index_redshifts<ptsz->N_redshift_dndlnM;index_redshifts++){


        fprintf(fp,"%e\t",ptsz->dndlnM_at_z_and_M[index_redshifts][index_masses]);
        }
        fprintf(fp,"\n");
     }
      fclose(fp);


      }

 if (ptsz->has_sz_trispec && ptsz->write_sz>0){
      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "szpowerspectrum_r_cl_clp",
                  ".txt");

      fp=fopen(Filepath, "w");

      int index_l_prime;
      for (index_l=0;index_l<ptsz->nlSZ;index_l++){
       for (index_l_prime=0;index_l_prime<ptsz->nlSZ;index_l_prime++) {
            fprintf(fp,"%e\t",ptsz->r_cl_clp[index_l][index_l_prime]/ptsz->r_cl_clp[index_l][index_l]);
         }
         fprintf(fp,"\n");
      }
      fclose(fp);


      sprintf(Filepath,
            "%s%s%s",
            ptsz->root,
            "szpowerspectrum_covariance_matrix_1e12_Dl_yxy_with_fsky",
            ".txt");

      fp=fopen(Filepath, "w");
      double ** cov_y_y;
      class_alloc(cov_y_y,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
      //int index_l_prime;
      for (index_l=0;index_l<ptsz->nlSZ;index_l++){
        class_alloc(cov_y_y[index_l],ptsz->nlSZ*sizeof(double),ptsz->error_message);
       for (index_l_prime=0;index_l_prime<index_l+1;index_l_prime++) {
            double ell = ptsz->ell[index_l];
            double ell_prime= ptsz->ell[index_l_prime];

            double dl_fac = ell*(ell+1.)/(2.*_PI_)*ell_prime*(ell_prime+1.)/(2.*_PI_);
            double m_l_lp;
            if (index_l == index_l_prime)
                  m_l_lp = ptsz->cov_cl_cl[index_l];
            else
                  m_l_lp = ptsz->tllprime_sz[index_l][index_l_prime]/ptsz->Omega_survey;

            cov_y_y[index_l][index_l_prime] = m_l_lp*dl_fac;
            cov_y_y[index_l_prime][index_l] = cov_y_y[index_l][index_l_prime];

            //fprintf(fp,"%e\t",m_l_lp*dl_fac);
         }
         //fprintf(fp,"\n");
      }

      for (index_l=0;index_l<ptsz->nlSZ;index_l++){
        for (index_l_prime=0;index_l_prime<ptsz->nlSZ;index_l_prime++) {
          fprintf(fp,"%e\t",cov_y_y[index_l][index_l_prime]);
        }
        fprintf(fp,"\n");
      }
      fclose(fp);


  for (index_l=0;index_l<ptsz->nlSZ;index_l++){
    free(cov_y_y[index_l]);
  }
  free(cov_y_y);
}




   //}

   return _SUCCESS_;
}




int show_preamble_messages(struct background * pba,
                          struct thermo * pth,
                          struct nonlinear * pnl,
                          struct primordial * ppm,
                          struct tszspectrum * ptsz){


   double tau;
   int first_index_back = 0;
   double * pvecback;
   double OmegaM;

   class_alloc(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);


      class_call(background_tau_of_z(pba,
                                     0.0, //TBC: z1SZ?
                                     &tau
                                     ),
                      ptsz->error_message,
                      ptsz->error_message
                      );

     class_call(background_at_tau(pba,
                                  tau,
                                  pba->long_info,
                                  pba->inter_normal,
                                  &first_index_back,
                                  pvecback),
                      ptsz->error_message,
                      ptsz->error_message
                      );

      // ptsz->Rho_crit_0 =
      // (3./(8.*_PI_*_G_*_M_sun_))
      // *pow(_Mpc_over_m_,1)
      // *pow(_c_,2)
      // *pvecback[pba->index_bg_rho_crit]
      // /pow(pba->h,2);
      //
      //
      // ptsz->Omega_m_0 = pvecback[pba->index_bg_Omega_m];
      // ptsz->Omega_r_0 = pvecback[pba->index_bg_Omega_r];
      //
      // ptsz->Omega_ncdm_0 = ptsz->Omega_m_0
      // -pba->Omega0_b
      // -pba->Omega0_cdm;
      //
      // ptsz->Omega0_b = pba->Omega0_b;
      //
      // if (ptsz->f_b_gas == -1.){
      //   ptsz->f_b_gas = pba->Omega0_b/ptsz->Omega_m_0;
      // }

      if (pba->Omega0_lambda != 0.) OmegaM = 1-pba->Omega0_lambda;
      else OmegaM = 1-pba->Omega0_fld;


   if ( pnl->has_pk_m == _TRUE_)
      ptsz->Sigma8OmegaM_SZ = pnl->sigma8[pnl->index_pk_m]*pow(OmegaM/0.28,3./8.);
      //quantities at surface of last scattering
      // conformal distance to redshift of last scattering in Mpc/h


      if (ptsz->sz_verbose > 0)
      {
         printf("Class_sz computations\n");

         //   double z_star;  /**< redshift at which photon optical depth crosses one */
         printf("->lss at z_star = %e\n",pth->z_star);


         //Halo mass function
         if (ptsz->MF==3) printf("->HMF: Jenkins et al 2001 @ M180m\n");
         if (ptsz->MF==1) printf("->HMF: Tinker et al 2010 @ M200m\n");
         if (ptsz->MF==4) printf("->HMF: Tinker et al 2008 @ M200m\n");
         if (ptsz->MF==5) printf("->HMF: Tinker et al 2008 @ M500c\n");
         if (ptsz->MF==6) printf("->HMF: Tinker et al 2008 @ M1600m\n");
         if (ptsz->MF==2) printf("->HMF: Bocquet et al 2015 @ M200m\n");
         if (ptsz->MF==7) printf("->HMF: Bocquet et al 2015 @ M500c\n");
         if (ptsz->MF==8) printf("->HMF: Tinker et al 2008 @ M200c\n");

         if (ptsz->has_sz_ps
           + ptsz->has_sz_2halo
          != _FALSE_){
          printf("->Computing Compton-y quantities\n");
          if (ptsz->has_completeness_for_ps_SZ == 1){
          char sz_cc_type[_ARGUMENT_LENGTH_MAX_];
          if (ptsz->which_ps_sz == 0)
             sprintf(sz_cc_type,"%s","total");
          else if (ptsz->which_ps_sz == 1) //ps resolved
            sprintf(sz_cc_type,"%s","resolved");
          else if (ptsz->which_ps_sz == 2) //ps unresolved
            sprintf(sz_cc_type,"%s","unresolved");

          printf("->Using completeness formalism corresponding to: %s\n", sz_cc_type);
          printf("  with signal-to-noise cut-off = %.3e\n",ptsz->sn_cutoff);
          }

          }
         //Pressure profile
if   (ptsz->has_sz_ps
    + ptsz->has_sz_trispec
    + ptsz->has_sz_2halo
    + ptsz->has_sz_m_y_y_2h
    + ptsz->has_sz_m_y_y_1h
    + ptsz->has_sz_te_y_y
    + ptsz->has_tSZ_tSZ_tSZ_1halo
    + ptsz->has_tSZ_tSZ_tSZ_2h
    + ptsz->has_tSZ_tSZ_tSZ_3h
    + ptsz->has_kSZ_kSZ_tSZ_1h
    + ptsz->has_kSZ_kSZ_tSZ_2h
    + ptsz->has_kSZ_kSZ_tSZ_3h
    + ptsz->has_tSZ_gal_1h
    + ptsz->has_tSZ_gal_2h
    + ptsz->has_tSZ_cib_1h
    + ptsz->has_tSZ_cib_2h
    + ptsz->has_tSZ_lens_1h
    + ptsz->has_tSZ_lens_2h
    + ptsz->has_isw_tsz
    + ptsz->has_mean_y
    + ptsz->has_dydz
    > 0
){
         if (ptsz->pressure_profile == 0)
            printf("->Pressure Profile:  Planck 2013\n");
         if (ptsz->pressure_profile == 2)
            printf("->Pressure Profile:  Arnaud et al 2010\n");
         if (ptsz->pressure_profile == 3){
            printf("->Pressure Profile:  Custom. GNFW\n");
            printf("P0GNFW = %e\n",ptsz->P0GNFW);
            printf("c500 = %e\n",ptsz->c500);
            printf("gammaGNFW = %e\n",ptsz->gammaGNFW);
            printf("alphaGNFW = %e\n",ptsz->alphaGNFW);
            printf("betaGNFW = %e\n",ptsz->betaGNFW);
         }
         if (ptsz->pressure_profile == 4)
            printf("->Pressure Profile:  Battaglia et al 2012\n");
}
         //Concentration-Mass relations
         // if (ptsz->MF!=5 && ptsz->MF!=7)
         // {
         printf("The following concentration-mass relation is set.\n");
         printf("(Only used when NFW profiles are involved, e.g., for conversions between different mass definitions.)\n");
            if (ptsz->concentration_parameter==0)
               printf("->C-M relation:  Duffy et al 2008\n");
            if (ptsz->concentration_parameter==1)
               printf("->C-M relation:  Seljak 2000\n");
            if (ptsz->concentration_parameter==2)
               printf("->C-M relation:  Klypin 2010\n");
            if (ptsz->concentration_parameter==3)
               printf("->C-M relation:  Sanchez-Conde 2014\n");
            if (ptsz->concentration_parameter==4)
               printf("->C-M relation:  Zhao 2009\n");
            if (ptsz->concentration_parameter==5)
               printf("->C-M relation:  DM14\n");
            if (ptsz->concentration_parameter==6)
               printf("->C-M relation:  Bhattacharya et al 2013\n");


         // }
         printf("->h = %e\n",pba->h);
         printf("->OmegaM (all except DE/Lambda) = %e\n",OmegaM);
         printf("->OmegaL = %e\n",1.-OmegaM);
         if ( pnl->has_pk_m == _TRUE_)
            printf("->sigma8 = %e\n",pnl->sigma8[pnl->index_pk_m]);
         printf("->Bias B = %e\n",ptsz->HSEbias);printf("->Bias b = %e\n",1.-1./ptsz->HSEbias);
      }

   if ( pnl->has_pk_m == _TRUE_){
  ptsz->sigma8_Pcb = pnl->sigma8[pnl->index_pk_cb];
   if (ptsz->sz_verbose > 0)
      printf("->sigma8_cb= %e\n",ptsz->sigma8_Pcb);
    }



    free(pvecback);

return _SUCCESS_;
}


double gnu_tsz_of_nu_in_ghz(double nu_in_ghz,
                            double Tcmb){

      double frequency_in_Hz = 1e9*nu_in_ghz;

      return ((_h_P_*frequency_in_Hz/(_k_B_*Tcmb))*(1./tanh((_h_P_*frequency_in_Hz
              /(_k_B_*Tcmb))
              /2.))-4.);
}

int show_results(struct background * pba,
                         struct nonlinear * pnl,
                         struct primordial * ppm,
                         struct tszspectrum * ptsz){


  if (ptsz->has_sz_ps){
printf("\n\n");
printf("########################################\n");
printf("tSZ power spectrum 1-halo term:\n");
printf("########################################\n");
printf("\n");

   int index_l;
   for (index_l=0;index_l<ptsz->nlSZ;index_l++){

   printf("ell = %e\t\t cl_y_y (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_sz_1h[index_l]);


 if (ptsz->has_sz_trispec){

    printf("\n");
      int index_l_prime;
      for (index_l_prime=0;index_l_prime<index_l+1;index_l_prime++)
      {

        if (index_l_prime == index_l)
            printf("ell = %e \t\t ell_prime = %e\t\t trispectrum/Omega_survey = %e\t\t sigma_g_cl_squared_binned = %e\n",
                  ptsz->ell[index_l],
                  ptsz->ell[index_l_prime],
                  ptsz->tllprime_sz[index_l][index_l_prime]/ptsz->Omega_survey,
                  ptsz->sig_cl_squared_binned[index_l]);
        else
            printf("ell = %e \t\t ell_prime = %e\t\t trispectrum/Omega_survey = %e \n",
                      ptsz->ell[index_l],
                      ptsz->ell[index_l_prime],
                      ptsz->tllprime_sz[index_l][index_l_prime]/ptsz->Omega_survey);
      }
printf("\n");
}

 if (ptsz->has_sz_cov_Y_N){
      int index_M_bins;
      for (index_M_bins=0;index_M_bins<ptsz->nbins_M;index_M_bins++)
      {


         printf("ell = %.3e\t M_min = %.3e\t M_max = %.3e\t cov_Y_N = %.3e\t cov_Y_N [ssc] = %.3e\t cov_N_N = %.3e\t cov_N_N [ssc] = %.3e\t cov_Y_Y = %.3e\t cov_Y_Y [ssc] = %.3e\n",
                   ptsz->ell[index_l],
                   ptsz->cov_Y_N_mass_bin_edges[index_M_bins],
                   ptsz->cov_Y_N_mass_bin_edges[index_M_bins+1],
                   ptsz->cov_Y_N[index_l][index_M_bins],
                   ptsz->cov_Y_N_next_order[index_l][index_M_bins],
                   ptsz->cov_N_N[index_M_bins],
                   ptsz->cov_N_N_hsv[index_M_bins][index_M_bins],
                   ptsz->cov_cl_cl[index_l],
                   ptsz->cov_Y_Y_ssc[index_l][index_l]);
      }

        printf("\n");
}
 }
}

 if (ptsz->has_sz_2halo){


printf("\n\n");
printf("########################################\n");
printf("tSZ power spectrum 2-halo term:\n");
printf("########################################\n");
printf("\n");

int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){
 printf("ell = %e\t\t cl_y_y (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_sz_2h[index_l]);
}
 }



 if (ptsz->has_mean_y)
   printf("mean_y =  %e \n",ptsz->y_monopole);
if (ptsz->has_hmf){
  printf("\n\n");
  printf("######################################################\n");
  printf("total number of halos (accounting for completeness):\n");
  printf("######################################################\n");
  printf("\n");

   printf("N_tot =  %e (per steradian)\n",ptsz->hmf_int);
   printf("N_tot =  %e (full-sky)\n",ptsz->hmf_int*4.*_PI_); // 1deg2 = 3.046174198e-4 sr
   printf("N_tot =  %e (fsky_from_skyfracs = %.3e)\n",ptsz->hmf_int*4.*_PI_*ptsz->fsky_from_skyfracs,ptsz->fsky_from_skyfracs);
 }

if (ptsz->has_kSZ_kSZ_gal_1h){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_gal_1h  = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gal_1h[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gal_covmat){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_t2t2f = %e  cov_ll_kSZ_kSZ_gal = %e \n",ptsz->ell[index_l],ptsz->cl_t2t2f[index_l],ptsz->cov_ll_kSZ_kSZ_gal[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gal_lensing_term){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t kSZ_kSZ_gal (lensing) = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gal_lensing_term[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gal_1h_fft){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_gal_1h_fft = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gal_1h_fft[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gal_2h_fft){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_gal_2h_fft = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gal_2h_fft[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gal_3h_fft){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_gal_3h_fft = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gal_3h_fft[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gal_2h){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_gal_2h = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gal_2h[index_l]);
}
}
if (ptsz->has_kSZ_kSZ_gal_3h){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_gal_3h = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gal_3h[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gal_hf){
printf("\n\n");
printf("###################################################\n");
printf("kSZ x kSZ x galaxy power spectrum halofit approach:\n");
printf("###################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_gal_hf (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gal_hf[index_l]);
}
}


if (ptsz->has_tSZ_tSZ_tSZ_1halo){
printf("\n\n");
printf("###################################################\n");
printf("1-halo tsz bispectrum:\n");
printf("###################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t bl_tSZ_tSZ_tSZ (1h) = %e \n",ptsz->ell[index_l],ptsz->b_tSZ_tSZ_tSZ_1halo[index_l]);
}
}

if (ptsz->has_tSZ_tSZ_tSZ_2h){
printf("\n\n");
printf("###################################################\n");
printf("2-halo tsz bispectrum:\n");
printf("###################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t bl_tSZ_tSZ_tSZ (2h) = %e \n",ptsz->ell[index_l],ptsz->b_tSZ_tSZ_tSZ_2h[index_l]);
}
}

if (ptsz->has_tSZ_tSZ_tSZ_3h){
printf("\n\n");
printf("###################################################\n");
printf("2-halo tsz bispectrum:\n");
printf("###################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t bl_tSZ_tSZ_tSZ (3h) = %e \n",ptsz->ell[index_l],ptsz->b_tSZ_tSZ_tSZ_3h[index_l]);
}
}


if (ptsz->has_kSZ_kSZ_lens_covmat){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cov_ll_kSZ_kSZ_lens = %e \n",ptsz->ell[index_l],ptsz->cov_ll_kSZ_kSZ_lens[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_lens_lensing_term){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t kSZ_kSZ_lens (lensing) = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_lens_lensing_term[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_lens_1h_fft){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_lens_1h_fft = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_lens_1h_fft[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_lens_2h_fft){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_lens_2h_fft = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_lens_2h_fft[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_lens_3h_fft){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_lens_3h_fft = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_lens_3h_fft[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_lens_hf){
printf("\n\n");
printf("###################################################\n");
printf("kSZ x kSZ x cmb lensing power spectrum halofit approach:\n");
printf("###################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_lens_hf (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_lens_hf[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gallens_covmat){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cov_ll_kSZ_kSZ_gallens = %e \n",ptsz->ell[index_l],ptsz->cov_ll_kSZ_kSZ_gallens[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gallens_lensing_term){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t kSZ_kSZ_gallens (lensing) = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gallens_lensing_term[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gallens_1h_fft){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_gallens_1h_fft = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gallens_1h_fft[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gallens_2h_fft){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_gallens_2h_fft = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gallens_2h_fft[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gallens_3h_fft){
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_gallens_3h_fft = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gallens_3h_fft[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_gallens_hf){
printf("\n\n");
printf("###################################################\n");
printf("kSZ x kSZ x galaxy lensing power spectrum halofit approach:\n");
printf("###################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_kSZ_kSZ_gallens_hf (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gallens_hf[index_l]);
}
}

if (ptsz->has_gallens_gallens_1h){
printf("\n\n");
printf("##########################################################\n");
printf("galaxy lensing x galaxy lensing power spectrum 1-halo term:\n");
printf("##########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gallens_gallens (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_gallens_gallens_1h[index_l]);
}
}
if (ptsz->has_gallens_gallens_2h){
printf("\n\n");
printf("##########################################################\n");
printf("galaxy lensing x galaxy lensing power spectrum 2-halo term:\n");
printf("##########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gallens_gallens (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_gallens_gallens_2h[index_l]);
}
}
if (ptsz->has_gallens_lens_1h){
printf("\n\n");
printf("##########################################################\n");
printf("galaxy lensing x cmb lensing power spectrum 1-halo term:\n");
printf("##########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gallens_lens (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_gallens_lens_1h[index_l]);
}
}
if (ptsz->has_gallens_lens_2h){
printf("\n\n");
printf("##########################################################\n");
printf("galaxy lensing x cmb lensing power spectrum 2-halo term:\n");
printf("##########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gallens_lens (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_gallens_lens_2h[index_l]);
}
}
if (ptsz->has_tSZ_gal_1h){
printf("\n\n");
printf("########################################\n");
printf("tSZ x galaxy power spectrum 1-halo term:\n");
printf("########################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_y_gal (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_1h[index_l]);
}
}
if (ptsz->has_tSZ_gal_2h){
printf("\n\n");
printf("########################################\n");
printf("tSZ x galaxy power spectrum 2-halo term:\n");
printf("########################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_y_gal (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_tSZ_gal_2h[index_l]);
}
}

if (ptsz->has_tSZ_lens_1h){
printf("\n\n");
printf("########################################\n");
printf("tSZ x lensing power spectrum 1-halo term:\n");
printf("########################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_y_kappa (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_tSZ_lens_1h[index_l]);
}
}
if (ptsz->has_tSZ_lens_2h){
printf("\n\n");
printf("########################################\n");
printf("tSZ x lensing power spectrum 2-halo term:\n");
printf("########################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_y_kappa (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_tSZ_lens_2h[index_l]);
}
}

if (ptsz->has_sz_rates){
printf("\n\n");
printf("########################################\n");
printf("sz rates:\n");
printf("########################################\n");
printf("\n");
// int index_rate;
// for (index_rate=0;index_rate<ptsz->szcat_size;index_rate++){
//
// printf("cluster id = %d\trate = %.4e \n",index_rate,ptsz->szrate[index_rate]);
// }

printf("cluster id = %d\trate = %.4e \n",0,ptsz->szrate[0]);
printf("cluster id = %d\trate = %.4e \n",ptsz->szcat_size-1,ptsz->szrate[ptsz->szcat_size-1]);

}

// if (ptsz->has_gal_cib_1h){
// printf("\n\n");
// printf("########################################\n");
// printf("tSZ x galaxy power spectrum 1-halo term:\n");
// printf("########################################\n");
// printf("\n");
// int index_l;
// for (index_l=0;index_l<ptsz->nlSZ;index_l++){
//
// printf("ell = %e\t\t cl_gal_cib (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_cib_1h[index_l]);
// }
// }
// if (ptsz->has_gal_cib_2h){
// printf("\n\n");
// printf("########################################\n");
// printf("tSZ x galaxy power spectrum 2-halo term:\n");
// printf("########################################\n");
// printf("\n");
// int index_l;
// for (index_l=0;index_l<ptsz->nlSZ;index_l++){
//
// printf("ell = %e\t\t cl_gal_cib (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_cib_2h[index_l]);
// }
// }


if (ptsz->has_gal_gal_1h){
printf("\n\n");
printf("########################################\n");
printf("galaxy x galaxy power spectrum 1-halo term:\n");
printf("########################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gal_gal (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_gal_1h[index_l]);
}
}
if (ptsz->has_gal_gal_2h){
printf("\n\n");
printf("###########################################\n");
printf("galaxy x galaxy power spectrum 2-halo term:\n");
printf("###########################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gal_gal (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_gal_2h[index_l]);
}
}

if (ptsz->has_gal_gal_hf){
printf("\n\n");
printf("###########################################\n");
printf("galaxy x galaxy effective approach :\n");
printf("###########################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gal_gal (effective approach) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_gal_hf[index_l]);
}
}

if (ptsz->has_tSZ_lensmag_1h){
printf("\n\n");
printf("#######################################################\n");
printf("tSZ x lensing magnification power spectrum 1-halo term:\n");
printf("#######################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_y_lensmag (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_tSZ_lensmag_1h[index_l]);
}
}
if (ptsz->has_tSZ_lensmag_2h){
printf("\n\n");
printf("#######################################################\n");
printf("tSZ x lensing magnification power spectrum 2-halo term:\n");
printf("#######################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_y_lensmag (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_tSZ_lensmag_2h[index_l]);
}
}
if (ptsz->has_gal_lens_1h){
printf("\n\n");
printf("##########################################################\n");
printf("galaxy x lensing power spectrum 1-halo term:\n");
printf("##########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gal_lens (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_lens_1h[index_l]);
}
}
if (ptsz->has_gal_lens_2h){
printf("\n\n");
printf("##########################################################\n");
printf("galaxy x lensing power spectrum 2-halo term:\n");
printf("##########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gal_lens (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_lens_2h[index_l]);
}
}

if (ptsz->has_gal_lens_hf){
printf("\n\n");
printf("##########################################################\n");
printf("galaxy x lensing power spectrum (effective approach):\n");
printf("##########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gal_lens (hf) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_lens_hf[index_l]);
}
}
if (ptsz->has_gal_lensmag_1h){
printf("\n\n");
printf("##########################################################\n");
printf("galaxy x lensing magnification power spectrum 1-halo term:\n");
printf("##########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gal_lensmag (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_lensmag_1h[index_l]);
}
}
if (ptsz->has_gal_lensmag_2h){
printf("\n\n");
printf("##########################################################\n");
printf("galaxy x lensing magnification power spectrum 2-halo term:\n");
printf("##########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gal_lensmag (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_lensmag_2h[index_l]);
}
}
if (ptsz->has_gal_gallens_1h){
printf("\n\n");
printf("##########################################################\n");
printf("galaxy x galaxy lensing power spectrum 1-halo term:\n");
printf("##########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gal_gallens (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_gallens_1h[index_l]);
}
}
if (ptsz->has_gal_gallens_2h){
printf("\n\n");
printf("##########################################################\n");
printf("galaxy x galaxy lensing power spectrum 2-halo term:\n");
printf("##########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gal_gallens (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_gallens_2h[index_l]);
}
}
if (ptsz->has_gal_lensmag_hf){
printf("\n\n");
printf("########################################################################\n");
printf("galaxy x lensing magnification power spectrum (effective approach):\n");
printf("########################################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_gal_lensmag (hf) = %e \n",ptsz->ell[index_l],ptsz->cl_gal_lensmag_hf[index_l]);
}
}
if (ptsz->has_lensmag_lensmag_1h){
printf("\n\n");
printf("#########################################################################\n");
printf("lensing magnification x lensing magnification power spectrum 1-halo term:\n");
printf("#########################################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_lensmag_lensmag (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_lensmag_lensmag_1h[index_l]);
}
}
if (ptsz->has_lensmag_lensmag_2h){
printf("\n\n");
printf("#########################################################################\n");
printf("lensing magnification x lensing magnification power spectrum 2-halo term:\n");
printf("#########################################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_lensmag_lensmag (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_lensmag_lensmag_2h[index_l]);
}
}
if (ptsz->has_lensmag_lensmag_hf){
printf("\n\n");
printf("#######################################################################################\n");
printf("lensing magnification x lensing magnification power spectrum (effective approach):\n");
printf("#######################################################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_lensmag_lensmag (hf) = %e \n",ptsz->ell[index_l],ptsz->cl_lensmag_lensmag_hf[index_l]);
}
}

if (ptsz->has_lens_lensmag_1h){
printf("\n\n");
printf("###########################################################\n");
printf("lensing x lensing magnification power spectrum 1-halo term:\n");
printf("###########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_lens_lensmag (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_lens_lensmag_1h[index_l]);
}
}
if (ptsz->has_lens_lensmag_2h){
printf("\n\n");
printf("###########################################################\n");
printf("lensing x lensing magnification power spectrum 2-halo term:\n");
printf("###########################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_lens_lensmag (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_lens_lensmag_2h[index_l]);
}
}

if (ptsz->has_lens_lensmag_hf){
printf("\n\n");
printf("########################################################################\n");
printf("lensing x lensing magnification power spectrum (effective approach):\n");
printf("########################################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_lens_lensmag (hf) = %e \n",ptsz->ell[index_l],ptsz->cl_lens_lensmag_hf[index_l]);
}
}

if (ptsz->has_pk_HI_at_z_1h){
printf("\n\n");
printf("#######################################\n");
printf("HI x HI power spectrum (1h):\n");
printf("#######################################\n");
printf("\n");
int index_k;
printf("ptsz->n_k_for_pk_hm = %d\n",ptsz->n_k_for_pk_hm);
for (index_k=0;index_k<ptsz->n_k_for_pk_hm;index_k++){

printf("k = %e [h/Mpc]\t\t pk_HI (1h) = %e  [Mpc/h]^3\n",ptsz->k_for_pk_hm[index_k],ptsz->pk_HI_at_z_1h[index_k]);
}
}

if (ptsz->has_pk_HI_at_z_2h){
printf("\n\n");
printf("#######################################\n");
printf("HI x HI power spectrum (2h):\n");
printf("#######################################\n");
printf("\n");
int index_k;
for (index_k=0;index_k<ptsz->n_k_for_pk_hm;index_k++){

printf("k = %e [h/Mpc]\t\t pk_HI (2h) = %e  [Mpc/h]^3\n",ptsz->k_for_pk_hm[index_k],ptsz->pk_HI_at_z_2h[index_k]);
}
}



if (ptsz->has_kSZ_kSZ_1h){
printf("\n\n");
printf("########################################################################\n");
printf("kSZ x kSZ power spectrum 1-halo term:\n");
printf("########################################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_ksz_ksz (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_1h[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_2h){
printf("\n\n");
printf("########################################################################\n");
printf("kSZ x kSZ power spectrum 2-halo term:\n");
printf("########################################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_ksz_ksz (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_2h[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_tSZ_1h){
printf("\n\n");
printf("########################################################################\n");
printf("kSZ x kSZ x tSZ bispectrum 1-halo term:\n");
printf("########################################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t b (1h) = %e \n",ptsz->ell[index_l],ptsz->b_kSZ_kSZ_tSZ_1h[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_tSZ_2h){
printf("\n\n");
printf("########################################################################\n");
printf("kSZ x kSZ x tSZ bispectrum 2-halo term:\n");
printf("########################################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t b (2h) = %e \n",ptsz->ell[index_l],ptsz->b_kSZ_kSZ_tSZ_2h[index_l]);
}
}

if (ptsz->has_kSZ_kSZ_tSZ_3h){
printf("\n\n");
printf("########################################################################\n");
printf("kSZ x kSZ x tSZ bispectrum 3-halo term:\n");
printf("########################################################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t b (3h) = %e \n",ptsz->ell[index_l],ptsz->b_kSZ_kSZ_tSZ_3h[index_l]);
}
}

if (ptsz->has_lens_lens_1h){
printf("\n\n");
printf("#######################################\n");
printf("lens x lens power spectrum 1-halo term:\n");
printf("#######################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_lens_lens (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_lens_lens_1h[index_l]);
}
}
if (ptsz->has_lens_lens_2h){
printf("\n\n");
printf("#######################################\n");
printf("lens x lens power spectrum 2-halo term:\n");
printf("#######################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_lens_lens (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_lens_lens_2h[index_l]);
}
}
if (ptsz->has_lens_lens_hf){
printf("\n\n");
printf("#######################################\n");
printf("lens x lens power spectrum (hf):\n");
printf("#######################################\n");
printf("\n");
int index_l;
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_lens_lens (hf) = %e \n",ptsz->ell[index_l],ptsz->cl_lens_lens_hf[index_l]);
}
}


if (ptsz->has_lens_cib_2h){
printf("\n\n");
printf("#######################################\n");
printf("lensing x cib power spectrum w-halo term:\n");
printf("#######################################\n");
printf("\n");
int index_l;
int index_nu;
// int index_nu_prime;
for (index_nu=0;index_nu<ptsz->cib_frequency_list_num;index_nu++){
    // for (index_nu_prime=0;index_nu_prime<ptsz->cib_frequency_list_num;index_nu_prime++){
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_lens_cib (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_lens_cib_2h[index_nu][index_l]);
}
// }
}
}


if (ptsz->has_cib_cib_1h){
printf("\n\n");
printf("#######################################\n");
printf("cib x cib power spectrum 1-halo term:\n");
printf("#######################################\n");
printf("\n");
int index_l;
int index_nu;
int index_nu_prime;
for (index_nu=0;index_nu<ptsz->cib_frequency_list_num;index_nu++){
    for (index_nu_prime=0;index_nu_prime<ptsz->cib_frequency_list_num;index_nu_prime++){
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_cib_cib (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_cib_cib_1h[index_nu][index_nu_prime][index_l]);
}
}
}
}

if (ptsz->has_ngal_ngal_1h){
printf("\n\n");
printf("#######################################\n");
printf("gal[n] x gal[n] power spectrum 1-halo term:\n");
printf("#######################################\n");
printf("\n");
int index_l;
int index_g;
int index_g_prime;
for (index_g=0;index_g<ptsz->galaxy_samples_list_num;index_g++){
    for (index_g_prime=0;index_g_prime<ptsz->galaxy_samples_list_num;index_g_prime++){
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_ngal[%d]_ngal[%d] (1h) = %e \n",ptsz->ell[index_l],index_g,index_g_prime,ptsz->cl_ngal_ngal_1h[index_g][index_g_prime][index_l]);
}
}
}
}

if (ptsz->has_ngal_ngal_2h){
printf("\n\n");
printf("#######################################\n");
printf("gal[n] x gal[n] power spectrum 2-halo term:\n");
printf("#######################################\n");
printf("\n");
int index_l;
int index_g;
int index_g_prime;
for (index_g=0;index_g<ptsz->galaxy_samples_list_num;index_g++){
    for (index_g_prime=0;index_g_prime<ptsz->galaxy_samples_list_num;index_g_prime++){
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_ngal[%d]_ngal[%d] (2h) = %e \n",ptsz->ell[index_l],index_g,index_g_prime,ptsz->cl_ngal_ngal_2h[index_g][index_g_prime][index_l]);
}
}
}
}

if (ptsz->has_ngal_ngal_hf){
printf("\n\n");
printf("#######################################\n");
printf("gal[n] x gal[n] power spectrum halofit method:\n");
printf("#######################################\n");
printf("\n");
int index_l;
int index_g;
int index_g_prime;
for (index_g=0;index_g<ptsz->galaxy_samples_list_num;index_g++){
    for (index_g_prime=0;index_g_prime<ptsz->galaxy_samples_list_num;index_g_prime++){
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_ngal[%d]_ngal[%d] (hf) = %e \n",ptsz->ell[index_l],index_g,index_g_prime,ptsz->cl_ngal_ngal_hf[index_g][index_g_prime][index_l]);
}
}
}
}

if (ptsz->has_ngal_lens_1h){
printf("\n\n");
printf("#######################################\n");
printf("gal[n] x lens power spectrum 1-halo term:\n");
printf("#######################################\n");
printf("\n");
int index_l;
int index_g;

for (index_g=0;index_g<ptsz->galaxy_samples_list_num;index_g++){

for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_ngal[%d]_lens (1h) = %e \n",ptsz->ell[index_l],index_g,ptsz->cl_ngal_lens_1h[index_g][index_l]);
}
}
}


if (ptsz->has_ngal_lens_2h){
printf("\n\n");
printf("#######################################\n");
printf("gal[n] x lens power spectrum 2-halo term:\n");
printf("#######################################\n");
printf("\n");
int index_l;
int index_g;

for (index_g=0;index_g<ptsz->galaxy_samples_list_num;index_g++){

for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_ngal[%d]_lens (2h) = %e \n",ptsz->ell[index_l],index_g,ptsz->cl_ngal_lens_2h[index_g][index_l]);
}
}
}


if (ptsz->has_ngal_lens_hf){
printf("\n\n");
printf("#######################################\n");
printf("gal[n] x lens power spectrum halofit method:\n");
printf("#######################################\n");
printf("\n");
int index_l;
int index_g;

for (index_g=0;index_g<ptsz->galaxy_samples_list_num;index_g++){

for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_ngal[%d]_lens (hf) = %e \n",ptsz->ell[index_l],index_g,ptsz->cl_ngal_lens_hf[index_g][index_l]);
}
}
}





if (ptsz->has_cib_cib_2h){
printf("\n\n");
printf("#######################################\n");
printf("cib x cib power spectrum 2-halo term:\n");
printf("#######################################\n");
printf("\n");
int index_l;
int index_nu;
int index_nu_prime;
for (index_nu=0;index_nu<ptsz->cib_frequency_list_num;index_nu++){
    for (index_nu_prime=0;index_nu_prime<ptsz->cib_frequency_list_num;index_nu_prime++){
for (index_l=0;index_l<ptsz->nlSZ;index_l++){

printf("ell = %e\t\t cl_cib_cib (2h) = %e \n",ptsz->ell[index_l],ptsz->cl_cib_cib_2h[index_nu][index_nu_prime][index_l]);
}
}
}
}

if (ptsz->has_cib_monopole){
printf("\n\n");
printf("#######################################\n");
printf("cib monopole:\n");
printf("#######################################\n");
printf("\n");
int index_nu;
for (index_nu=0;index_nu<ptsz->n_frequencies_for_cib;index_nu++){


printf("nu = %e\t\t cib intensity = %e \n",ptsz->frequencies_for_cib[index_nu],ptsz->cib_monopole[index_nu]);
}

}

if (ptsz->has_cib_shotnoise){
printf("\n\n");
printf("#######################################\n");
printf("cib shotnoise:\n");
printf("#######################################\n");
printf("\n");
int index_nu;
for (index_nu=0;index_nu<ptsz->cib_frequency_list_num;index_nu++){

double nu = ptsz->cib_frequency_list[index_nu];
printf("nu = %e\t\t cib shotnoise = %e \n",nu,ptsz->cib_shotnoise[index_nu]);
}

}

   return _SUCCESS_;
}




int select_multipole_array(struct tszspectrum * ptsz)
{

   class_alloc(ptsz->ell_plc,
                     26*sizeof(double),
                     ptsz->error_message);

   ptsz->ell_plc[0] = 10.;
   ptsz->ell_plc[1] = 13.5;
   ptsz->ell_plc[2] = 18.;
   ptsz->ell_plc[3] = 23.5;
   ptsz->ell_plc[4] = 30.5;
   ptsz->ell_plc[5] = 40.;
   ptsz->ell_plc[6] = 52.5;
   ptsz->ell_plc[7] = 68.5;
   ptsz->ell_plc[8] = 89.5;
   ptsz->ell_plc[9] = 117.;
   ptsz->ell_plc[10] = 152.5;
   ptsz->ell_plc[11] = 198.;
   ptsz->ell_plc[12] = 257.5;
   ptsz->ell_plc[13] = 335.5;
   ptsz->ell_plc[14] = 436.5;
   ptsz->ell_plc[15] = 567.5;
   ptsz->ell_plc[16] = 738.;
   ptsz->ell_plc[17] = 959.5;
   ptsz->ell_plc[18] = 1247.5;
   ptsz->ell_plc[19] = 1622.;
   ptsz->ell_plc[20] = 2109.;
   ptsz->ell_plc[21] = 2742.;
   ptsz->ell_plc[22] = 3000.;
   ptsz->ell_plc[23] = 7000.;
   ptsz->ell_plc[24] = 10000.;
   ptsz->ell_plc[25] = 20000.;

   class_alloc(ptsz->ell_plc_no_low_ell,
                     13*sizeof(double),
                     ptsz->error_message);
   ptsz->ell_plc_no_low_ell[0] = 40.;
   ptsz->ell_plc_no_low_ell[1] = 52.5;
   ptsz->ell_plc_no_low_ell[2] = 68.5;
   ptsz->ell_plc_no_low_ell[3] = 89.5;
   ptsz->ell_plc_no_low_ell[4] = 117.;
   ptsz->ell_plc_no_low_ell[5] = 152.5;
   ptsz->ell_plc_no_low_ell[6] = 198.;
   ptsz->ell_plc_no_low_ell[7] = 257.5;
   ptsz->ell_plc_no_low_ell[8] = 335.5;
   ptsz->ell_plc_no_low_ell[9] = 436.5;
   ptsz->ell_plc_no_low_ell[10] = 567.5;
   ptsz->ell_plc_no_low_ell[11] = 738.;
   ptsz->ell_plc_no_low_ell[12] = 959.5;

   class_alloc(ptsz->ell_plc_low,
                     29*sizeof(double),
                     ptsz->error_message);

   ptsz->ell_plc_low[0] = 2.;
   ptsz->ell_plc_low[1] = 5.;
   ptsz->ell_plc_low[2] = 8.;
   ptsz->ell_plc_low[3] = 10.;
   ptsz->ell_plc_low[4] = 13.5;
   ptsz->ell_plc_low[5] = 18.;
   ptsz->ell_plc_low[6] = 23.5;
   ptsz->ell_plc_low[7] = 30.5;
   ptsz->ell_plc_low[8] = 40.;
   ptsz->ell_plc_low[9] = 52.5;
   ptsz->ell_plc_low[10] = 68.5;
   ptsz->ell_plc_low[11] = 89.5;
   ptsz->ell_plc_low[12] = 117.;
   ptsz->ell_plc_low[13] = 152.5;
   ptsz->ell_plc_low[14] = 198.;
   ptsz->ell_plc_low[15] = 257.5;
   ptsz->ell_plc_low[16] = 335.5;
   ptsz->ell_plc_low[17] = 436.5;
   ptsz->ell_plc_low[18] = 567.5;
   ptsz->ell_plc_low[19] = 738.;
   ptsz->ell_plc_low[20] = 959.5;
   ptsz->ell_plc_low[21] = 1247.5;
   ptsz->ell_plc_low[22] = 1622.;
   ptsz->ell_plc_low[23] = 2109.;
   ptsz->ell_plc_low[24] = 2742.;
   ptsz->ell_plc_low[25] = 3000.;
   ptsz->ell_plc_low[26] = 7000.;
   ptsz->ell_plc_low[27] = 10000.;
   ptsz->ell_plc_low[28] = 20000.;



   class_alloc(ptsz->ell_trispectrum,
                     10*sizeof(double),
                     ptsz->error_message);

   ptsz->ell_trispectrum[0] = 100.;
   ptsz->ell_trispectrum[1] = 200.;
   ptsz->ell_trispectrum[2] = 300.;
   ptsz->ell_trispectrum[3] = 500.;
   ptsz->ell_trispectrum[4] = 1000.;


   ptsz->ell_trispectrum[5] = 2000.;
   ptsz->ell_trispectrum[6] = 3000.;
   ptsz->ell_trispectrum[7] = 5000.;
   ptsz->ell_trispectrum[8] = 10000.;
   ptsz->ell_trispectrum[9] = 20000.;

   class_alloc(ptsz->frequencies_for_cib,
                     10*sizeof(double),
                     ptsz->error_message);
    if (ptsz->dfreq == 0.){
      ptsz->n_frequencies_for_cib = (int)((log(ptsz->freq_max) - log(ptsz->freq_min))/ptsz->dlogfreq) + 1;
      class_realloc(ptsz->frequencies_for_cib,
                    ptsz->frequencies_for_cib,
                    ptsz->n_frequencies_for_cib*sizeof(double),
                    ptsz->error_message);

      int i;
      for (i=0;i<ptsz->n_frequencies_for_cib;i++)
       ptsz->frequencies_for_cib[i] = exp(log(ptsz->freq_min)+i*ptsz->dlogfreq);
     }
     else{
       ptsz->n_frequencies_for_cib = (int)(ptsz->freq_max -ptsz->freq_min)/ptsz->dfreq + 1;
       class_realloc(ptsz->frequencies_for_cib,
                     ptsz->frequencies_for_cib,
                     ptsz->n_frequencies_for_cib*sizeof(double),
                     ptsz->error_message);
     int i;
     for (i=0;i<ptsz->n_frequencies_for_cib;i++)
      ptsz->frequencies_for_cib[i] = ptsz->freq_min+i*ptsz->dfreq;
     }


   class_alloc(ptsz->ell_mock,
                     10*sizeof(double),
                     ptsz->error_message);




   if(ptsz->ell_sz==4){
    if (ptsz->dell == 0.){
      ptsz->nlSZ = (int)((log(ptsz->ell_max_mock) - log(ptsz->ell_min_mock))/ptsz->dlogell) + 1;
      class_realloc(ptsz->ell_mock,
                    ptsz->ell_mock,
                    ptsz->nlSZ*sizeof(double),
                    ptsz->error_message);

      int i;
      for (i=0;i<ptsz->nlSZ;i++)
       ptsz->ell_mock[i] = exp(log(ptsz->ell_min_mock)+i*ptsz->dlogell);
     }
     else{
       ptsz->nlSZ = (int)(ptsz->ell_max_mock -ptsz->ell_min_mock)/ptsz->dell + 1;
       class_realloc(ptsz->ell_mock,
                     ptsz->ell_mock,
                     ptsz->nlSZ*sizeof(double),
                     ptsz->error_message);
     int i;
     for (i=0;i<ptsz->nlSZ;i++)
      ptsz->ell_mock[i] = ptsz->ell_min_mock+i*ptsz->dell;
     }
   }

  if (ptsz->has_kSZ_kSZ_gal_1h
   || ptsz->has_kSZ_kSZ_gal_2h
   || ptsz->has_kSZ_kSZ_gal_3h
   || ptsz->has_kSZ_kSZ_gal_hf
   || ptsz->has_kSZ_kSZ_gallens_hf
   || ptsz->has_kSZ_kSZ_lens_hf
   // || ptsz->has_kSZ_kSZ_gal_covmat // not needed for this one...
   || ptsz->has_kSZ_kSZ_gal_lensing_term
   || ptsz->has_kSZ_kSZ_gallens_lensing_term
   || ptsz->has_kSZ_kSZ_lens_lensing_term
   || ptsz->has_kSZ_kSZ_lensmag_1halo){
  class_alloc(ptsz->ell_kSZ2_gal_multipole_grid,
              ptsz->N_kSZ2_gal_multipole_grid*sizeof(double),
              ptsz->error_message);
  class_alloc(ptsz->theta_kSZ2_gal_theta_grid,
              ptsz->N_kSZ2_gal_theta_grid*sizeof(double),
              ptsz->error_message);

  int index_l;
  for (index_l=0;index_l<ptsz->N_kSZ2_gal_multipole_grid;index_l++){
  // ptsz->ell_kSZ2_gal_multipole_grid[index_l] = exp(log(2.) + index_l*(log(5.*ptsz->ell_max_mock)
  //                                                  - log(2.))/(ptsz->N_kSZ2_gal_multipole_grid-1.));

   ptsz->ell_kSZ2_gal_multipole_grid[index_l] = log(ptsz->ell_min_kSZ2_gal_multipole_grid) + index_l*(log(ptsz->ell_max_kSZ2_gal_multipole_grid)
                                                    - log(ptsz->ell_min_kSZ2_gal_multipole_grid))/(ptsz->N_kSZ2_gal_multipole_grid-1.);


   // ptsz->ell_kSZ2_gal_multipole_grid[index_l] = 2. + index_l*(1e4
   //                                                  - 2.)/(ptsz->N_kSZ2_gal_multipole_grid-1.);


  }

  double theta_min = 0.;
  double theta_max = 2.*_PI_;
  // double theta_min = -1.;
  // double theta_max = 1.;

  for (index_l=0;index_l<ptsz->N_kSZ2_gal_theta_grid;index_l++){

   ptsz->theta_kSZ2_gal_theta_grid[index_l] = theta_min + index_l*(theta_max
                                                    - theta_min)/(ptsz->N_kSZ2_gal_theta_grid-1.);

  }


  }

   class_alloc(ptsz->ell,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   int index_l;
   for (index_l=0;index_l<ptsz->nlSZ;index_l++)
   {
      if (ptsz->ell_sz == 0)
         ptsz->ell[index_l] = pow(10.,1.+1.*index_l*0.2);
      else if (ptsz->ell_sz == 1){
         ptsz->ell[index_l] = ptsz->ell_plc[index_l];
         //ptsz->nlSZ = 18;
       }
      else if (ptsz->ell_sz == 2)
         ptsz->ell[index_l] = ptsz->ell_trispectrum[index_l];
      else if (ptsz->ell_sz == 3)
         ptsz->ell[index_l] = ptsz->ell_plc_low[index_l];
      else if (ptsz->ell_sz == 4)
         ptsz->ell[index_l] = ptsz->ell_mock[index_l];
     else if (ptsz->ell_sz == 5)
        ptsz->ell[index_l] = ptsz->ell_plc_no_low_ell[index_l];
   }

   if (ptsz->sz_verbose>0){
     printf("\n");
     printf("-> If angular power spectra requested:\n");
     printf("-> Number of multipoles computed by class_sz = %d\n", ptsz->nlSZ);
     printf("\n");
   }


   free(ptsz->ell_trispectrum);
   free(ptsz->ell_plc);
   free(ptsz->ell_plc_low);
   free(ptsz->ell_plc_no_low_ell);
   free(ptsz->ell_mock);

   return _SUCCESS_;

}



int initialise_and_allocate_memory(struct tszspectrum * ptsz){


  if (ptsz->use_redshift_dependent_M_min){
    load_M_min_of_z(ptsz);
  }

   ptsz->Omega_survey = 4.*_PI_*ptsz->f_sky;

// here make sure we bypass all the hm consistency in case its not asked.
   if (ptsz->hm_consistency == 1){
   ptsz->m_min_counter_terms =  ptsz->M1SZ;
 }
 else{
   ptsz->m_min_counter_terms = -1;
 }

   if (ptsz->has_sz_ps
      +ptsz->has_sz_2halo
      +ptsz->has_sz_te_y_y
      // +ptsz->has_sz_cov_N_Cl
      +ptsz->has_sz_cov_Y_N
      +ptsz->has_sz_cov_Y_N_next_order
      +ptsz->has_sz_cov_Y_Y_ssc
      +ptsz->has_mean_y
      +ptsz->has_dydz
      +ptsz->has_sz_m_y_y_1h
      +ptsz->has_sz_m_y_y_2h
      +ptsz->has_tSZ_tSZ_tSZ_1halo
      +ptsz->has_tSZ_tSZ_tSZ_2h
      +ptsz->has_tSZ_tSZ_tSZ_3h
      +ptsz->has_kSZ_kSZ_tSZ_1h
      +ptsz->has_kSZ_kSZ_tSZ_2h
      +ptsz->has_kSZ_kSZ_tSZ_3h
      +ptsz->has_tSZ_lens_1h
      +ptsz->has_tSZ_lens_2h
      +ptsz->has_tSZ_lensmag_2h
      +ptsz->has_tSZ_lensmag_1h
      +ptsz->has_tSZ_gal_1h
      +ptsz->has_tSZ_gal_2h
      +ptsz->has_tSZ_cib_1h
      +ptsz->has_tSZ_cib_2h
      +ptsz->has_sz_trispec != _FALSE_){
      ptsz->has_electron_pressure = 1;


      // if battaglia pressure profile need 200c
      if (ptsz->pressure_profile == 4){
        ptsz->has_200c = 1;
        ptsz->delta_def_electron_pressure = 1;
      }
      else if (ptsz->pressure_profile == 0){ // Planck 2013
        ptsz->has_500c = 1;
        ptsz->delta_def_electron_pressure = 2;
      }
      else if (ptsz->pressure_profile == 2){ // Arnaud et al 2010
        ptsz->has_500c = 1;
        ptsz->delta_def_electron_pressure = 2;
      }
      else if (ptsz->pressure_profile == 3){ // Custom. GNFW
        // ptsz->delta_def_electron_pressure  = 2;
        if (ptsz->delta_def_electron_pressure == 2){
          ptsz->has_500c = 1;
        }
        else if (ptsz->delta_def_electron_pressure == 1){
          ptsz->has_200c = 1;
        }
        else if (ptsz->delta_def_electron_pressure == 0){
          ptsz->has_200m = 1;
        }
      }
    }



   // if (
   //    +ptsz->has_kSZ_kSZ_gal_1h
   //    +ptsz->has_kSZ_kSZ_gal_1h_fft
   //    +ptsz->has_kSZ_kSZ_gal_2h_fft
   //    +ptsz->has_kSZ_kSZ_gal_3h_fft
   //    +ptsz->has_kSZ_kSZ_gal_2h
   //    +ptsz->has_kSZ_kSZ_gal_3h
   //    +ptsz->has_kSZ_kSZ_lensmag_1halo
   //    != _FALSE_)
   //    ptsz->has_200c = 1;

   if (ptsz->has_kSZ_kSZ_gal_1h
      +ptsz->has_kSZ_kSZ_gal_1h_fft
      +ptsz->has_kSZ_kSZ_gal_2h_fft
      +ptsz->has_kSZ_kSZ_gal_3h_fft
      +ptsz->has_kSZ_kSZ_gallens_1h_fft
      +ptsz->has_kSZ_kSZ_gallens_2h_fft
      +ptsz->has_kSZ_kSZ_gallens_3h_fft
      +ptsz->has_kSZ_kSZ_lens_1h_fft
      +ptsz->has_kSZ_kSZ_lens_2h_fft
      +ptsz->has_kSZ_kSZ_lens_3h_fft
      +ptsz->has_kSZ_kSZ_gal_2h
      +ptsz->has_kSZ_kSZ_gal_3h
      +ptsz->has_kSZ_kSZ_tSZ_1h
      +ptsz->has_kSZ_kSZ_tSZ_2h
      +ptsz->has_gas_density_profile_2h
      +ptsz->has_kSZ_kSZ_1h
      +ptsz->has_kSZ_kSZ_2h
      +ptsz->has_kSZ_kSZ_tSZ_3h
      +ptsz->has_kSZ_kSZ_lensmag_1halo != _FALSE_){

      ptsz->has_electron_density = 1;

      if (ptsz->tau_profile == 1){ // battaglia tau profile, need m200c
          ptsz->has_200c = 1;
          ptsz->delta_def_electron_density = 1;
          }
      else if (ptsz->tau_profile == 0){ // nfw profile
        if (ptsz->delta_def_electron_density == 0){
          ptsz->has_200m = 1;
        }
        else if (ptsz->delta_def_electron_density == 1){
          ptsz->has_200c = 1;
        }
        else if (ptsz->delta_def_electron_density == 2){
          ptsz->has_500c = 1;
        }
      }
    }

  if (ptsz->has_kSZ_kSZ_gal_1h_fft
  ||  ptsz->has_kSZ_kSZ_gal_2h_fft
  ||  ptsz->has_kSZ_kSZ_gal_3h_fft
  ||  ptsz->has_kSZ_kSZ_gal_covmat
  ||  ptsz->has_kSZ_kSZ_gallens_1h_fft
  ||  ptsz->has_kSZ_kSZ_gallens_2h_fft
  ||  ptsz->has_kSZ_kSZ_gallens_3h_fft
  ||  ptsz->has_kSZ_kSZ_gallens_covmat
  ||  ptsz->has_kSZ_kSZ_lens_1h_fft
  ||  ptsz->has_kSZ_kSZ_lens_2h_fft
  ||  ptsz->has_kSZ_kSZ_lens_3h_fft
  ||  ptsz->has_kSZ_kSZ_lens_covmat
  ||  ptsz->convert_cls_to_gamma
  ||  ptsz->has_pk_b_at_z_2h
  ||  ptsz->has_sz_counts_fft
  ||  ptsz->has_gas_density_profile_2h
  ||  ptsz->has_gas_pressure_profile_2h
  ||  ptsz->use_fft_for_profiles_transform
  ||  ptsz->has_n5k){
    if(ptsz->sz_verbose>1) printf("constructing fftw plan\n");
    // ptsz->N_samp_fftw = 100;
    fftw_complex* a_tmp;
    fftw_complex* b_tmp;
    a_tmp = fftw_alloc_complex(ptsz->N_samp_fftw);
    b_tmp = fftw_alloc_complex(ptsz->N_samp_fftw);
    ptsz->forward_plan = fftw_plan_dft_1d(ptsz->N_samp_fftw,
                                    (fftw_complex*) a_tmp,
                                    (fftw_complex*) b_tmp,
                                    -1, FFTW_ESTIMATE);
    ptsz->reverse_plan = fftw_plan_dft_1d(ptsz->N_samp_fftw,
                                    (fftw_complex*) b_tmp,
                                    (fftw_complex*) b_tmp,
                                    +1, FFTW_ESTIMATE);

    fftw_free(a_tmp);
    fftw_free(b_tmp);
  }


   if (ptsz->has_tSZ_gal_1h
      +ptsz->has_tSZ_gal_2h
      +ptsz->has_kSZ_kSZ_gal_1h
      +ptsz->has_kSZ_kSZ_gal_1h_fft
      +ptsz->has_kSZ_kSZ_gal_2h_fft
      +ptsz->has_kSZ_kSZ_gal_3h_fft
      +ptsz->has_kSZ_kSZ_gal_2h
      +ptsz->has_kSZ_kSZ_gal_3h
      //+ptsz->has_kSZ_kSZ_gal_hf
      +ptsz->has_kSZ_kSZ_lensmag_1halo
      +ptsz->has_gal_gal_1h
      +ptsz->has_ngal_ngal_1h
      +ptsz->has_ngal_ngal_2h
      +ptsz->has_ngal_lens_1h
      +ptsz->has_ngal_lens_2h
      // +ptsz->has_gal_gal_hf
      +ptsz->has_gal_lens_1h
      +ptsz->has_gal_lens_2h
      +ptsz->has_gal_cib_1h
      +ptsz->has_gal_cib_2h
      +ptsz->has_gal_lensmag_1h
      +ptsz->has_gal_lensmag_2h
      +ptsz->has_gal_gallens_1h
      +ptsz->has_gal_gallens_2h
      +ptsz->has_tSZ_lensmag_1h
      +ptsz->has_tSZ_lensmag_2h
      // +ptsz->has_lensmag_lensmag_1h
      // +ptsz->has_lensmag_lensmag_2h
      // +ptsz->has_lens_lensmag_1h
      // +ptsz->has_lens_lensmag_2h
      +ptsz->has_gal_gal_2h != _FALSE_){
      ptsz->has_galaxy = 1;

        if (ptsz->delta_def_galaxies == 0){
          ptsz->has_200m = 1;
        }
        else if (ptsz->delta_def_galaxies == 1){
          ptsz->has_200c = 1;
        }
        else if (ptsz->delta_def_galaxies == 2){
          ptsz->has_500c = 1;
        }

  }



if (ptsz->has_gal_cib_1h
  +ptsz->has_gal_cib_2h
  +ptsz->has_tSZ_cib_1h
  +ptsz->has_tSZ_cib_2h
  +ptsz->has_lens_cib_1h
  +ptsz->has_lens_cib_2h
  +ptsz->has_cib_cib_1h
  +ptsz->has_cib_cib_2h
  +ptsz->has_cib_monopole
  +ptsz->has_cib_shotnoise
  +ptsz->has_dcib0dz
  )
  {
    ptsz->has_cib = 1;

    if (ptsz->delta_def_cib == 0){
      ptsz->has_200m = 1;
    }
    else if (ptsz->delta_def_cib == 1){
      ptsz->has_200c = 1;
    }
    else if (ptsz->delta_def_cib == 2){
      ptsz->has_500c = 1;
    }

  }

if (
  ptsz->has_pk_at_z_1h
  +ptsz->has_pk_at_z_2h
  +ptsz->has_pk_em_at_z_1h
  +ptsz->has_pk_em_at_z_2h
  +ptsz->has_bk_at_z_1h
  +ptsz->has_bk_at_z_2h
  +ptsz->has_bk_at_z_3h
  // +ptsz->has_bk_at_z_hf
  )
  {
    ptsz->has_matter_density = 1;

    if (ptsz->delta_def_matter_density == 0){
      ptsz->has_200m = 1;
    }
    else if (ptsz->delta_def_matter_density == 1){
      ptsz->has_200c = 1;
    }
    else if (ptsz->delta_def_matter_density == 2){
      ptsz->has_500c = 1;
    }

  }
if (
  ptsz->has_bk_ttg_at_z_1h
  +ptsz->has_bk_ttg_at_z_2h
  +ptsz->has_bk_ttg_at_z_3h
  +ptsz->has_pk_bb_at_z_1h
  +ptsz->has_pk_bb_at_z_2h
  +ptsz->has_pk_b_at_z_2h
  +ptsz->has_pk_em_at_z_1h
  +ptsz->has_pk_em_at_z_2h
  // +ptsz->has_bk_at_z_hf
  )
  {
    ptsz->has_electron_density = 1;

    if (ptsz->delta_def_electron_density == 0){
      ptsz->has_200m = 1;
    }
    else if (ptsz->delta_def_electron_density == 1){
      ptsz->has_200c = 1;
    }
    else if (ptsz->delta_def_electron_density == 2){
      ptsz->has_500c = 1;
    }
    ptsz->has_galaxy = 1;

    if (ptsz->delta_def_galaxies == 0){
      ptsz->has_200m = 1;
    }
    else if (ptsz->delta_def_galaxies == 1){
      ptsz->has_200c = 1;
    }
    else if (ptsz->delta_def_galaxies == 2){
      ptsz->has_500c = 1;
    }

  }

if (ptsz->has_pk_HI_at_z_1h
   +ptsz->has_pk_HI_at_z_2h){
     ptsz->has_HI_density = 1;

     if (ptsz->delta_def_HI_density == 0){
       ptsz->has_200m = 1;
     }
     else if (ptsz->delta_def_HI_density == 1){
       ptsz->has_200c = 1;
     }
     else if (ptsz->delta_def_HI_density == 2){
       ptsz->has_500c = 1;
     }


   }


if (ptsz->has_kSZ_kSZ_lensmag_1halo
  +ptsz->has_gal_lens_1h
  +ptsz->has_gal_lens_2h
  +ptsz->has_ngal_lens_1h
  +ptsz->has_ngal_lens_2h
  +ptsz->has_gal_lensmag_1h
  +ptsz->has_gal_lensmag_2h
  +ptsz->has_gal_gallens_1h
  +ptsz->has_gal_gallens_2h
  +ptsz->has_gallens_gallens_1h
  +ptsz->has_gallens_gallens_2h
  +ptsz->has_gallens_lens_1h
  +ptsz->has_gallens_lens_2h
  +ptsz->has_kSZ_kSZ_gallens_1h_fft
  +ptsz->has_kSZ_kSZ_gallens_2h_fft
  +ptsz->has_kSZ_kSZ_gallens_3h_fft
  +ptsz->has_kSZ_kSZ_lens_1h_fft
  +ptsz->has_kSZ_kSZ_lens_2h_fft
  +ptsz->has_kSZ_kSZ_lens_3h_fft
  +ptsz->has_lens_lens_1h
  +ptsz->has_lens_lens_2h
  +ptsz->has_lens_lensmag_1h
  +ptsz->has_lens_lensmag_2h
  +ptsz->has_lensmag_lensmag_1h
  +ptsz->has_lensmag_lensmag_2h
  +ptsz->has_tSZ_lens_1h
  +ptsz->has_tSZ_lens_2h
  +ptsz->has_tSZ_lensmag_1h
  +ptsz->has_tSZ_lensmag_2h
  +ptsz->has_lens_cib_1h
  +ptsz->has_lens_cib_2h
  )
  {
    ptsz->has_lensing = 1;

    if (ptsz->delta_def_matter_density == 0){
      ptsz->has_200m = 1;
    }
    else if (ptsz->delta_def_matter_density == 1){
      ptsz->has_200c = 1;
    }
    else if (ptsz->delta_def_matter_density == 2){
      ptsz->has_500c = 1;
    }

  }


 if (ptsz->has_sz_counts == _TRUE_){
   ptsz->has_500c = 1;
 }

 if (ptsz->has_sz_rates == _TRUE_){
   ptsz->has_500c = 1;
 }

 // if (ptsz->has_cib
 //  +  ptsz->has_cib_monopole
 //  +  ptsz->has_dcib0dz
 //  +  ptsz->has_gal_gal_1h
 //  +  ptsz->has_gal_gal_2h
 //  +  ptsz->has_galaxy
 //   != _FALSE_){
 //   ptsz->has_200m = 1;
 // }

if (ptsz->need_hmf){
  if (ptsz->integrate_wrt_m200m == 1 && ptsz->has_500c == 1)
    ptsz->need_m200m_to_m500c = 1;

  if (ptsz->integrate_wrt_m200m == 1 && ptsz->has_200c == 1)
    ptsz->need_m200m_to_m200c = 1;

  if (ptsz->integrate_wrt_m200c == 1 && ptsz->has_200m == 1)
    ptsz->need_m200c_to_m200m = 1;

  if (ptsz->integrate_wrt_m200c == 1 && ptsz->has_500c == 1){
    ptsz->need_m200c_to_m500c = 1;
  }

  if (ptsz->integrate_wrt_m500c == 1 && ptsz->has_200c == 1)
    ptsz->need_m500c_to_m200c = 1;
}

   class_alloc(ptsz->ln_x_for_pp, ptsz->ln_x_size_for_pp*sizeof(double),ptsz->error_message);
   int i;
   for (i=0;i<ptsz->ln_x_size_for_pp;i++){
      ptsz->ln_x_for_pp[i] = log(ptsz->x_inSZ)+i*(log(ptsz->x_outSZ)-log(ptsz->x_inSZ))/(ptsz->ln_x_size_for_pp-1.);
      //printf("M = %e \t i=%d\n",ptsz->M_bins[i],i);
   }

   class_alloc(ptsz->x_for_pp, ptsz->x_size_for_pp*sizeof(double),ptsz->error_message);
   //int i;
   for (i=0;i<ptsz->x_size_for_pp;i++){
      ptsz->x_for_pp[i] = ptsz->x_inSZ+i*(ptsz->x_outSZ-ptsz->x_inSZ)/(ptsz->x_size_for_pp-1.);
      //printf("M = %e \t i=%d\n",ptsz->M_bins[i],i);
   }





   //mass bins for covariance between cluster counts and power spectrum
   //ptsz->dlogM = 1.;
   //ptsz->nbins_M = (int)((log(ptsz->M2SZ) - log(ptsz->M1SZ))/ptsz->dlogM) + 1;
   //ptsz->nbins_M = 20;
   //printf("mass bins for N-C_l cov, nbins = %d\n", ptsz->nbins_M);
   //M_bins is the array of bin edges from m_min to m_max
   class_alloc(ptsz->M_bins,
                        (ptsz->nbins_M+1)*sizeof(double),
                        ptsz->error_message);
   for (i=0;i<ptsz->nbins_M;i++){
      ptsz->M_bins[i] = exp(log(ptsz->M1SZ)+i*(log(ptsz->M2SZ)-log(ptsz->M1SZ))/(ptsz->nbins_M));
      //printf("M = %e \t i=%d\n",ptsz->M_bins[i],i);
   }

   class_alloc(ptsz->cov_Y_N_mass_bin_edges,
                        (ptsz->nbins_M+1)*sizeof(double),
                        ptsz->error_message);
   for (i=0;i<ptsz->nbins_M + 1;i++){
      ptsz->cov_Y_N_mass_bin_edges[i] = exp(log(ptsz->M1SZ)+i*(log(ptsz->M2SZ)-log(ptsz->M1SZ))/(ptsz->nbins_M));
   }

   class_alloc(ptsz->dndlnM_array_z,
              ptsz->N_redshift_dndlnM*sizeof(double),
              ptsz->error_message);
   for (i=0;i<ptsz->N_redshift_dndlnM;i++){
     if (ptsz->N_redshift_dndlnM == 1)
     ptsz->dndlnM_array_z[i] =  ptsz->z1SZ_dndlnM;
     else
      ptsz->dndlnM_array_z[i] =  ptsz->z1SZ_dndlnM +
                                  +i*(ptsz->z2SZ_dndlnM-ptsz->z1SZ_dndlnM)
                                  /(ptsz->N_redshift_dndlnM-1.);
 }



   class_alloc(ptsz->dndlnM_array_m,
              ptsz->N_mass_dndlnM*sizeof(double),
              ptsz->error_message);
   for (i=0;i<ptsz->N_mass_dndlnM;i++){
     if (ptsz->N_mass_dndlnM == 1)
     ptsz->dndlnM_array_m[i] = ptsz->M1SZ_dndlnM;
     else
      ptsz->dndlnM_array_m[i] =  exp(log(ptsz->M1SZ_dndlnM)+i*(log(ptsz->M2SZ_dndlnM)-log(ptsz->M1SZ_dndlnM))/(ptsz->N_mass_dndlnM-1.));
  }




   //quantities that are specified for each integrand
   ptsz->index_integrand_id = 0;
   ptsz->index_multipole_for_lensing_profile = ptsz->index_integrand_id + 1;
   ptsz->index_characteristic_multipole_for_nfw_profile = ptsz->index_multipole_for_lensing_profile + 1;
   ptsz->index_multipole_for_nfw_profile = ptsz->index_characteristic_multipole_for_nfw_profile + 1;
   ptsz->index_multipole_for_tau_profile = ptsz->index_multipole_for_nfw_profile + 1;
   ptsz->index_dlnMdeltadlnM = ptsz->index_multipole_for_tau_profile + 1;
   ptsz->index_te_of_m = ptsz->index_dlnMdeltadlnM + 1;
   ptsz->index_multipole_for_pressure_profile = ptsz->index_te_of_m  + 1;
   ptsz->index_multipole_prime = ptsz->index_multipole_for_pressure_profile +1;
   ptsz->index_md = ptsz->index_multipole_prime +1;
   ptsz->index_Rho_crit = ptsz->index_md +1;
   ptsz->index_Delta_c  = ptsz->index_Rho_crit +1;
   ptsz->index_rVIR = ptsz->index_Delta_c +1;
   ptsz->index_cVIR = ptsz->index_rVIR +1;
   ptsz->index_mVIR = ptsz->index_cVIR +1;
   ptsz->index_m500 = ptsz->index_mVIR +1;
   ptsz->index_r500 = ptsz->index_m500 +1;
   ptsz->index_l500 = ptsz->index_r500 +1;
   ptsz->index_m200 = ptsz->index_l500 +1;
   ptsz->index_m200m = ptsz->index_m200 +1;
   ptsz->index_m1600m  = ptsz->index_m200m +1;
   ptsz->index_m180m  = ptsz->index_m1600m +1;
   ptsz->index_mass_for_hmf  = ptsz->index_m180m +1;
   ptsz->index_m500c = ptsz->index_mass_for_hmf +1;
   ptsz->index_r500c = ptsz->index_m500c +1;
   ptsz->index_Rh = ptsz->index_r500c +1;
   ptsz->index_mf = ptsz->index_Rh +1;
   ptsz->index_dlognudlogRh  = ptsz->index_mf +1;
   ptsz->index_lognu  = ptsz->index_dlognudlogRh +1;
   ptsz->index_dlogSigma2dlogRh  = ptsz->index_lognu +1;
   ptsz->index_dndlogRh  = ptsz->index_dlogSigma2dlogRh +1;
   ptsz->index_logSigma2  = ptsz->index_dndlogRh +1;
   ptsz->index_z = ptsz->index_logSigma2 +1;
   ptsz->index_rs = ptsz->index_z +1;
   ptsz->index_ls = ptsz->index_rs +1;
   ptsz->index_c200c = ptsz->index_ls +1;
   ptsz->index_m200c = ptsz->index_c200c +1;
   ptsz->index_r200c = ptsz->index_m200c +1;
   ptsz->index_l200c = ptsz->index_r200c +1;
   ptsz->index_multipole = ptsz->index_l200c +1;
   ptsz->index_pressure_profile = ptsz->index_multipole +1;
   ptsz->index_tau_profile = ptsz->index_pressure_profile +1;
   ptsz->index_lensing_profile = ptsz->index_tau_profile +1;
   ptsz->index_completeness = ptsz->index_lensing_profile +1;
   ptsz->index_volume = ptsz->index_completeness + 1;
   ptsz->index_chi2 = ptsz->index_volume + 1;
   ptsz->index_vrms2 = ptsz->index_chi2 + 1;
   ptsz->index_dgdz = ptsz->index_vrms2 + 1;
   ptsz->index_lensing_Sigma_crit = ptsz->index_dgdz + 1;
   ptsz->index_hmf = ptsz->index_lensing_Sigma_crit + 1;
   ptsz->index_halo_bias_b2 = ptsz->index_hmf + 1;
   ptsz->index_halo_bias = ptsz->index_halo_bias_b2 + 1;
   ptsz->index_k_value_for_halo_bias = ptsz->index_halo_bias +1;
   ptsz->index_pk_for_halo_bias = ptsz->index_k_value_for_halo_bias +1;
   ptsz->index_multipole_for_pk = ptsz->index_pk_for_halo_bias + 1;
   ptsz->index_mass_bin_1 =  ptsz->index_multipole_for_pk + 1;
   ptsz->index_mass_bin_2 = ptsz->index_mass_bin_1 + 1;
   ptsz->index_multipole_1 = ptsz->index_mass_bin_2  + 1;
   ptsz->index_multipole_2 = ptsz->index_multipole_1 + 1;
   ptsz->index_sigma2_hsv = ptsz->index_multipole_2 + 1;
   ptsz->index_part_id_cov_hsv = ptsz->index_sigma2_hsv + 1;
   ptsz->index_redshift_for_dndlnM = ptsz->index_part_id_cov_hsv + 1;
   ptsz->index_mass_for_dndlnM = ptsz->index_redshift_for_dndlnM +1;

   ptsz->index_phi_galaxy_counts = ptsz->index_mass_for_dndlnM + 1;
   ptsz->index_mean_galaxy_number_density = ptsz->index_phi_galaxy_counts+1;
   ptsz->index_c500c = ptsz->index_mean_galaxy_number_density+1;
   ptsz->index_multipole_for_galaxy_profile = ptsz->index_c500c+1;
   ptsz->index_multipole_for_truncated_nfw_profile = ptsz->index_multipole_for_galaxy_profile+1;
   ptsz->index_galaxy_profile =  ptsz->index_multipole_for_truncated_nfw_profile + 1;
   ptsz->index_multipole_for_cib_profile = ptsz->index_galaxy_profile + 1;
   ptsz->index_cib_profile = ptsz->index_multipole_for_cib_profile + 1;
   ptsz->index_frequency_for_cib_profile = ptsz->index_cib_profile + 1;
   ptsz->index_frequency_prime_for_cib_profile = ptsz->index_frequency_for_cib_profile + 1;
   ptsz->index_ngal_for_galaxy_profile = ptsz->index_frequency_prime_for_cib_profile + 1;
   ptsz->index_ngal_prime_for_galaxy_profile = ptsz->index_ngal_for_galaxy_profile + 1;
   ptsz->index_multipole_3 =  ptsz->index_ngal_prime_for_galaxy_profile + 1;
   ptsz->index_W_lensmag = ptsz->index_multipole_3 + 1;
   ptsz->index_c200m = ptsz->index_W_lensmag + 1;
   ptsz->index_r200m = ptsz->index_c200m + 1;
   ptsz->index_k_for_pk_hm = ptsz->index_r200m + 1;
   ptsz->index_density_profile =  ptsz->index_k_for_pk_hm + 1;

   ptsz->index_integrate_wrt_mvir = ptsz->index_density_profile + 1; // not needed
   ptsz->index_integrate_wrt_m500c = ptsz->index_integrate_wrt_mvir + 1; // not needed
   ptsz->index_integrate_wrt_m200m = ptsz->index_integrate_wrt_m500c + 1; // not needed

   ptsz->index_has_electron_pressure = ptsz->index_integrate_wrt_m200m + 1;
   ptsz->index_has_electron_density = ptsz->index_has_electron_pressure + 1;
   ptsz->index_has_HI_density = ptsz->index_has_electron_density + 1;
   ptsz->index_has_galaxy = ptsz->index_has_HI_density + 1;
   ptsz->index_has_matter_density = ptsz->index_has_galaxy + 1;
   ptsz->index_has_lensing = ptsz->index_has_matter_density + 1;
   ptsz->index_has_cib = ptsz->index_has_lensing + 1;
   ptsz->index_has_isw = ptsz->index_has_cib + 1;

   ptsz->index_has_vir = ptsz->index_has_isw + 1;
   ptsz->index_has_500c = ptsz->index_has_vir + 1;
   ptsz->index_has_200m = ptsz->index_has_500c + 1;
   ptsz->index_has_200c = ptsz->index_has_200m + 1;

   ptsz->index_ell_1 = ptsz->index_has_200c + 1;
   ptsz->index_ell_2 = ptsz->index_ell_1 + 1;
   ptsz->index_ell_3 = ptsz->index_ell_2 + 1;

   //quantities integrated over redshift
   ptsz->index_integral =  ptsz->index_ell_3 + 1;

   ptsz->index_integral_over_m = ptsz->index_integral+1;




   //integrands at m and z
   ptsz->index_integrand =  ptsz->index_integral_over_m + 1;


   // mass/concentration/radiues:
  ptsz->index_mass_for_galaxies = ptsz->index_integrand + 1;
  ptsz->index_mass_for_cib = ptsz->index_mass_for_galaxies + 1;
  ptsz->index_mass_for_matter_density = ptsz->index_mass_for_cib+ 1;
  ptsz->index_mass_for_electron_pressure = ptsz->index_mass_for_matter_density + 1;
  ptsz->index_mass_for_electron_density = ptsz->index_mass_for_electron_pressure + 1;
  ptsz->index_concentration_for_galaxies = ptsz->index_mass_for_electron_density + 1;
  ptsz->index_concentration_for_cib = ptsz->index_concentration_for_galaxies + 1;
  ptsz->index_concentration_for_matter_density = ptsz->index_concentration_for_cib  + 1;
  ptsz->index_concentration_for_electron_pressure = ptsz->index_concentration_for_matter_density + 1;
  ptsz->index_concentration_for_electron_density = ptsz->index_concentration_for_electron_pressure + 1;
  ptsz->index_radius_for_galaxies = ptsz->index_concentration_for_electron_density + 1;
  ptsz->index_radius_for_cib = ptsz->index_radius_for_galaxies + 1;
  ptsz->index_radius_for_matter_density = ptsz->index_radius_for_cib + 1;
  ptsz->index_radius_for_electron_pressure = ptsz->index_radius_for_matter_density + 1;
  ptsz->index_radius_for_electron_density = ptsz->index_radius_for_electron_pressure + 1;

  ptsz->index_W_gallens_sources = ptsz->index_radius_for_electron_density + 1;
  ptsz->index_szrate = ptsz->index_W_gallens_sources +1;

  ptsz->index_mass_for_HI_density = ptsz->index_szrate +1;
  ptsz->index_radius_for_HI_density = ptsz->index_mass_for_HI_density +1;
  ptsz->index_concentration_for_HI_density = ptsz->index_radius_for_HI_density+1;
   //final size of pvecsz vector
  ptsz->tsz_size  = ptsz->index_szrate + 1;


//
   //printf("cib_dim = %d\n",1000);

   if (ptsz->has_cib_cib_1h
     + ptsz->has_cib_cib_2h
     + ptsz->has_cib_shotnoise
     + ptsz->has_tSZ_cib_1h
     + ptsz->has_tSZ_cib_2h
     + ptsz->has_lens_cib_1h
     + ptsz->has_lens_cib_2h
     + ptsz->has_gal_cib_1h
     + ptsz->has_gal_cib_2h
     != _FALSE_){
   if (ptsz->sz_verbose>=1)
   printf("-> [cib] number of cib frequencies = %d\n",ptsz->cib_frequency_list_num);
   int index_cib_freq;
   if (ptsz->sz_verbose>=1){
   for (index_cib_freq = 0; index_cib_freq < ptsz->cib_frequency_list_num; index_cib_freq++)
   {
     printf("-> [cib] frequency %d = %.3e\n",index_cib_freq, ptsz->cib_frequency_list[index_cib_freq]);
   }
   }

   ptsz->cib_dim = ptsz->cib_frequency_list_num*(ptsz->cib_frequency_list_num+1)/2*ptsz->nlSZ;
   if (ptsz->sz_verbose>=1)
   printf("-> [cib] cib_dim = %d\n",ptsz->cib_dim);
   //exit(0);
   }


   if ((ptsz->has_ngal_ngal_1h
       +ptsz->has_ngal_ngal_2h
       +ptsz->has_ngal_ngal_hf
       +ptsz->has_ngal_lens_1h
       +ptsz->has_ngal_lens_2h
       +ptsz->has_ngal_lens_hf
     )
     != _FALSE_){
   if (ptsz->sz_verbose>=1)
   printf("-> [ngal] number of galaxy samples = %d\n",ptsz->galaxy_samples_list_num);
   int index_ngal;
   if (ptsz->sz_verbose>=1){
   for (index_ngal = 0; index_ngal < ptsz->galaxy_samples_list_num; index_ngal++)
   {
     printf("-> [ngal] sample %d has id #%d\n",index_ngal, ptsz->galaxy_samples_list[index_ngal]);
   }
   }

   ptsz->ngal_dim = ptsz->galaxy_samples_list_num*(ptsz->galaxy_samples_list_num+1)/2*ptsz->nlSZ;
   if (ptsz->sz_verbose>=1)
   printf("-> [ngal] ngal_dim = %d\n",ptsz->ngal_dim);
   //exit(0);
   }



   if (ptsz->has_pk_at_z_1h
     + ptsz->has_pk_at_z_2h
     + ptsz->has_pk_gg_at_z_1h
     + ptsz->has_pk_gg_at_z_2h
     + ptsz->has_pk_bb_at_z_1h
     + ptsz->has_pk_bb_at_z_2h
     + ptsz->has_pk_b_at_z_2h
     + ptsz->has_pk_em_at_z_1h
     + ptsz->has_pk_em_at_z_2h
     + ptsz->has_pk_HI_at_z_1h
     + ptsz->has_pk_HI_at_z_2h
     + ptsz->has_bk_at_z_1h
     + ptsz->has_bk_at_z_2h
     + ptsz->has_bk_at_z_3h
     + ptsz->has_bk_ttg_at_z_1h
     + ptsz->has_bk_ttg_at_z_2h
     + ptsz->has_bk_ttg_at_z_3h
     >= _TRUE_){

     class_alloc(ptsz->k_for_pk_hm,
                       10*sizeof(double),
                       ptsz->error_message);
    // printf("%.5e %.5e %.5e\n",ptsz->k_max_for_pk_hm,ptsz->k_min_for_pk_hm,ptsz->dlnk_for_pk_hm);
     ptsz->n_k_for_pk_hm = (int)((log(ptsz->k_max_for_pk_hm) - log(ptsz->k_min_for_pk_hm))/ptsz->dlnk_for_pk_hm) + 1;

     class_realloc(ptsz->k_for_pk_hm,
                   ptsz->k_for_pk_hm,
                   ptsz->n_k_for_pk_hm*sizeof(double),
                   ptsz->error_message);
                   int i;

   class_alloc(ptsz->pk_at_z_1h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->pk_at_z_2h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->pk_gg_at_z_1h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->pk_gg_at_z_2h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->pk_bb_at_z_1h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->pk_bb_at_z_2h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->pk_b_at_z_2h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->pk_em_at_z_1h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->pk_em_at_z_2h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->pk_HI_at_z_1h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->pk_HI_at_z_2h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   for (i=0;i<ptsz->n_k_for_pk_hm;i++){
    ptsz->k_for_pk_hm[i] = exp(log(ptsz->k_min_for_pk_hm)+i*ptsz->dlnk_for_pk_hm);
    ptsz->pk_at_z_1h[i] = 0.;
    ptsz->pk_at_z_2h[i] = 0.;
    ptsz->pk_gg_at_z_1h[i] = 0.;
    ptsz->pk_gg_at_z_2h[i] = 0.;
    ptsz->pk_bb_at_z_1h[i] = 0.;
    ptsz->pk_bb_at_z_2h[i] = 0.;
    ptsz->pk_b_at_z_2h[i] = 0.;
    ptsz->pk_em_at_z_1h[i] = 0.;
    ptsz->pk_em_at_z_2h[i] = 0.;
    ptsz->pk_HI_at_z_1h[i] = 0.;
    ptsz->pk_HI_at_z_2h[i] = 0.;
   }

   class_alloc(ptsz->bk_at_z_1h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->bk_at_z_2h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->bk_at_z_3h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);


   for (i=0;i<ptsz->n_k_for_pk_hm;i++){
    ptsz->k_for_pk_hm[i] = exp(log(ptsz->k_min_for_pk_hm)+i*ptsz->dlnk_for_pk_hm);
    ptsz->bk_at_z_1h[i] = 0.;
    ptsz->bk_at_z_2h[i] = 0.;
    ptsz->bk_at_z_3h[i] = 0.;
   }

   class_alloc(ptsz->bk_ttg_at_z_1h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->bk_ttg_at_z_2h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->bk_ttg_at_z_3h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);


   for (i=0;i<ptsz->n_k_for_pk_hm;i++){
    ptsz->k_for_pk_hm[i] = exp(log(ptsz->k_min_for_pk_hm)+i*ptsz->dlnk_for_pk_hm);
    ptsz->bk_ttg_at_z_1h[i] = 0.;
    ptsz->bk_ttg_at_z_2h[i] = 0.;
    ptsz->bk_ttg_at_z_3h[i] = 0.;
   }

   }

//exit(0);


   ptsz->hmf_int = 0.;
   ptsz->y_monopole = 0.;


   class_alloc(ptsz->cl_sz_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cib_monopole,sizeof(double *)*ptsz->n_frequencies_for_cib,ptsz->error_message);
   class_alloc(ptsz->cib_shotnoise,sizeof(double *)*ptsz->cib_frequency_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_isw_lens,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_isw_tsz,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_isw_auto,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_tSZ_gal_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_tSZ_gal_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_tSZ_lensmag_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_tSZ_lensmag_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   //class_alloc(ptsz->cl_tSZ_cib_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   //class_alloc(ptsz->cl_tSZ_cib_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);

   class_alloc(ptsz->cl_gal_gal_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_gal_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_gal_hf,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_lens_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_lens_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_lens_hf,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_lensmag_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_lensmag_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_gallens_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_gallens_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gallens_gallens_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gallens_gallens_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gallens_lens_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gallens_lens_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_lensmag_hf,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_lensmag_lensmag_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_lensmag_lensmag_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_lensmag_lensmag_hf,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_lens_lensmag_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_lens_lensmag_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_lens_lensmag_hf,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_lens_lens_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_lens_lens_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_lens_lens_hf,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_tSZ_lens_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_tSZ_lens_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   // if (ptsz->has_sz_rates)
   //  class_alloc(ptsz->szrate,sizeof(double *)*ptsz->szcat_size,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gal_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cov_ll_kSZ_kSZ_gal,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_t2t2f,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gal_lensing_term,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gal_1h_fft,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gal_2h_fft,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gal_3h_fft,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cov_ll_kSZ_kSZ_gallens,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gallens_lensing_term,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gallens_1h_fft,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gallens_2h_fft,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gallens_3h_fft,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gallens_hf,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cov_ll_kSZ_kSZ_lens,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_lens_lensing_term,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_lens_1h_fft,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_lens_2h_fft,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_lens_3h_fft,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_lens_hf,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->b_kSZ_kSZ_tSZ_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->b_kSZ_kSZ_tSZ_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->b_kSZ_kSZ_tSZ_3h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gal_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gal_3h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gal_hf,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_lensmag_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->b_tSZ_tSZ_tSZ_1halo,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->b_tSZ_tSZ_tSZ_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->b_tSZ_tSZ_tSZ_3h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_te_y_y,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->m_y_y_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->m_y_y_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cov_cl_cl,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->sig_cl_squared_binned,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);

   class_alloc(ptsz->cl_sz_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);

   class_alloc(ptsz->dndlnM_at_z_and_M,ptsz->N_redshift_dndlnM*sizeof(double *),ptsz->error_message);
   int index_z;
   for (index_z = 0; index_z<ptsz->N_redshift_dndlnM;index_z ++)
   {
     class_alloc(ptsz->dndlnM_at_z_and_M[index_z],(ptsz->N_mass_dndlnM)*sizeof(double),ptsz->error_message);
   }

   if (ptsz->has_sz_rates+ptsz->has_sz_counts_fft){
   class_alloc(ptsz->szrate,ptsz->szcat_size*sizeof(double *),ptsz->error_message);
   int irate;
   for (irate = 0; irate < ptsz->szcat_size ; irate++){
     ptsz->szrate[irate] = 0.;
   }
   // exit(0);
 }

   class_alloc(ptsz->tllprime_sz,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
   class_alloc(ptsz->trispectrum_ref,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
   class_alloc(ptsz->r_cl_clp,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
   class_alloc(ptsz->cov_N_N,ptsz->nbins_M*sizeof(double *),ptsz->error_message);
   class_alloc(ptsz->cov_N_N_hsv,ptsz->nbins_M*sizeof(double *),ptsz->error_message);
   class_alloc(ptsz->cov_Y_N,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
   class_alloc(ptsz->cov_Y_N_next_order,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
   class_alloc(ptsz->cov_Y_Y_ssc,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
   class_alloc(ptsz->r_Y_N,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
   int index_l,index_l_prime;


   for (index_l=0;index_l<ptsz->nlSZ;index_l++){
      ptsz->cl_sz_1h[index_l] = 0.;
      ptsz->cl_isw_lens[index_l] = 0.;
      ptsz->cl_isw_tsz[index_l] = 0.;
      ptsz->cl_isw_auto[index_l] = 0.;
      ptsz->cl_gal_gal_1h[index_l] = 0.;
      ptsz->cl_gal_gal_2h[index_l] = 0.;
      ptsz->cl_gal_gal_hf[index_l] = 0.;
      //ptsz->cl_tSZ_cib_1h[index_l] = 0.;
      //ptsz->cl_tSZ_cib_2h[index_l] = 0.;
      //ptsz->cl_cib_cib_1h[index_l] = 0.;
      //ptsz->cl_cib_cib_2h[index_l] = 0.;
      ptsz->cl_gal_lens_1h[index_l] = 0.;
      ptsz->cl_gal_lens_2h[index_l] = 0.;
      ptsz->cl_gal_lens_hf[index_l] = 0.;
      ptsz->cl_gal_lensmag_1h[index_l] = 0.;
      ptsz->cl_gal_lensmag_2h[index_l] = 0.;
      ptsz->cl_gal_gallens_1h[index_l] = 0.;
      ptsz->cl_gal_gallens_2h[index_l] = 0.;
      ptsz->cl_gallens_gallens_1h[index_l] = 0.;
      ptsz->cl_gallens_gallens_2h[index_l] = 0.;
      ptsz->cl_gallens_lens_1h[index_l] = 0.;
      ptsz->cl_gallens_lens_2h[index_l] = 0.;
      ptsz->cl_gal_lensmag_hf[index_l] = 0.;
      ptsz->cl_lensmag_lensmag_1h[index_l] = 0.;
      ptsz->cl_lensmag_lensmag_2h[index_l] = 0.;
      ptsz->cl_lensmag_lensmag_hf[index_l] = 0.;
      ptsz->cl_lens_lensmag_1h[index_l] = 0.;
      ptsz->cl_lens_lensmag_2h[index_l] = 0.;
      ptsz->cl_lens_lensmag_hf[index_l] = 0.;
      ptsz->cl_tSZ_gal_1h[index_l] = 0.;
      ptsz->cl_tSZ_gal_2h[index_l] = 0.;
      ptsz->cl_tSZ_lensmag_1h[index_l] = 0.;
      ptsz->cl_tSZ_lensmag_2h[index_l] = 0.;
      ptsz->cl_tSZ_lens_1h[index_l] = 0.;
      ptsz->cl_lens_lens_1h[index_l] = 0.;
      ptsz->cl_lens_lens_2h[index_l] = 0.;
      ptsz->cl_lens_lens_hf[index_l] = 0.;
      ptsz->cl_tSZ_lens_2h[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gal_1h[index_l] = 0.;
      ptsz->cov_ll_kSZ_kSZ_gal[index_l] = 0.;
      ptsz->cl_t2t2f[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gal_lensing_term[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gal_1h_fft[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gal_2h_fft[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gal_3h_fft[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gallens_1h_fft[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gallens_2h_fft[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gallens_3h_fft[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gallens_lensing_term[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gallens_hf[index_l] = 0.;
      ptsz->cov_ll_kSZ_kSZ_gallens[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_lens_1h_fft[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_lens_2h_fft[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_lens_3h_fft[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_lens_lensing_term[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_lens_hf[index_l] = 0.;
      ptsz->cov_ll_kSZ_kSZ_lens[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gal_2h[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gal_3h[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gal_hf[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_1h[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_2h[index_l] = 0.;
      ptsz->b_kSZ_kSZ_tSZ_1h[index_l] = 0.;
      ptsz->b_kSZ_kSZ_tSZ_2h[index_l] = 0.;
      ptsz->b_kSZ_kSZ_tSZ_3h[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_lensmag_1h[index_l] = 0.;
      ptsz->b_tSZ_tSZ_tSZ_1halo[index_l] = 0.;
      ptsz->b_tSZ_tSZ_tSZ_2h[index_l] = 0.;
      ptsz->b_tSZ_tSZ_tSZ_3h[index_l] = 0.;
      ptsz->cl_te_y_y[index_l] = 0.;
      ptsz->m_y_y_1h[index_l] = 0.;
      ptsz->m_y_y_2h[index_l] = 0.;
      ptsz->cl_sz_2h[index_l] = 0.;
      ptsz->cov_cl_cl[index_l] = 0.;
      ptsz->sig_cl_squared_binned[index_l] = 0.;

      class_alloc(ptsz->tllprime_sz[index_l],(index_l+1)*sizeof(double),ptsz->error_message);
      class_alloc(ptsz->trispectrum_ref[index_l],ptsz->nlSZ*sizeof(double),ptsz->error_message);
      class_alloc(ptsz->r_cl_clp[index_l],ptsz->nlSZ*sizeof(double),ptsz->error_message);
      class_alloc(ptsz->cov_Y_N[index_l],(ptsz->nbins_M)*sizeof(double),ptsz->error_message);
      class_alloc(ptsz->cov_Y_N_next_order[index_l],(ptsz->nbins_M)*sizeof(double),ptsz->error_message);
      class_alloc(ptsz->r_Y_N[index_l],(ptsz->nbins_M)*sizeof(double),ptsz->error_message);

      for (index_l_prime = 0; index_l_prime<index_l+1; index_l_prime ++){
         ptsz->tllprime_sz[index_l][index_l_prime] = 0.;
      }

      for (index_l_prime = 0; index_l_prime<ptsz->nlSZ; index_l_prime ++){
        ptsz->r_cl_clp[index_l][index_l_prime] = 0.;
        ptsz->trispectrum_ref[index_l][index_l_prime] = 0.;
          }



   int index_cov_Y_N_M_bins;
   for (index_cov_Y_N_M_bins = 0; index_cov_Y_N_M_bins<ptsz->nbins_M; index_cov_Y_N_M_bins ++){
      ptsz->cov_Y_N[index_l][index_cov_Y_N_M_bins] = 0.;
      ptsz->cov_Y_N_next_order[index_l][index_cov_Y_N_M_bins] = 0.;
      ptsz->r_Y_N[index_l][index_cov_Y_N_M_bins] = 0.;
   }
}

   int index_M_bins_1;
   for (index_M_bins_1 = 0; index_M_bins_1<ptsz->nbins_M; index_M_bins_1 ++){
      ptsz->cov_N_N[index_M_bins_1] = 0.;
      class_alloc(ptsz->cov_N_N_hsv[index_M_bins_1],(ptsz->nbins_M)*sizeof(double),ptsz->error_message);
      int index_M_bins_2;
      for (index_M_bins_2 = 0; index_M_bins_2<ptsz->nbins_M; index_M_bins_2 ++)
      ptsz->cov_N_N_hsv[index_M_bins_1][index_M_bins_2] = 0.;

   }

   int index_multipole_1;
   for (index_multipole_1 = 0; index_multipole_1<ptsz->nlSZ; index_multipole_1 ++){
      class_alloc(ptsz->cov_Y_Y_ssc[index_multipole_1],(ptsz->nlSZ)*sizeof(double),ptsz->error_message);
      int index_multipole_2;
      for (index_multipole_2 = 0; index_multipole_2<ptsz->nlSZ; index_multipole_2 ++)
      ptsz->cov_Y_Y_ssc[index_multipole_1][index_multipole_2] = 0.;

   }

   class_alloc(ptsz->cl_cib_cib_1h,sizeof(double ***)*ptsz->cib_frequency_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_cib_cib_2h,sizeof(double ***)*ptsz->cib_frequency_list_num,ptsz->error_message);

   int index_nu;
   int index_nu_prime;
   for (index_nu=0;index_nu<ptsz->n_frequencies_for_cib;index_nu++){
     ptsz->cib_monopole[index_nu] = 0.;
   }

   for (index_nu=0;index_nu<ptsz->cib_frequency_list_num;index_nu++){
     class_alloc(ptsz->cl_cib_cib_1h[index_nu],sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
     class_alloc(ptsz->cl_cib_cib_2h[index_nu],sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
     ptsz->cib_shotnoise[index_nu] = 0.;
        for (index_nu_prime=0;index_nu_prime<ptsz->cib_frequency_list_num;index_nu_prime++){
          class_alloc(ptsz->cl_cib_cib_1h[index_nu][index_nu_prime],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
          class_alloc(ptsz->cl_cib_cib_2h[index_nu][index_nu_prime],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
                  for (index_l=0;index_l<ptsz->nlSZ;index_l++){
     ptsz->cl_cib_cib_1h[index_nu][index_nu_prime][index_l] = 0.;
     ptsz->cl_cib_cib_2h[index_nu][index_nu_prime][index_l] = 0.;
   }
   }
   }

   if (ptsz->has_ngal_ngal_1h
     + ptsz->has_ngal_ngal_2h
     + ptsz->has_ngal_ngal_hf
     + ptsz->has_ngal_lens_1h
     + ptsz->has_ngal_lens_2h
     + ptsz->has_ngal_lens_hf
   ){

   class_alloc(ptsz->cl_ngal_ngal_1h,sizeof(double ***)*ptsz->galaxy_samples_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_ngal_ngal_2h,sizeof(double ***)*ptsz->galaxy_samples_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_ngal_ngal_hf,sizeof(double ***)*ptsz->galaxy_samples_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_ngal_lens_1h,sizeof(double **)*ptsz->galaxy_samples_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_ngal_lens_2h,sizeof(double **)*ptsz->galaxy_samples_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_ngal_lens_hf,sizeof(double **)*ptsz->galaxy_samples_list_num,ptsz->error_message);
   int index_g;
   int index_g_prime;
   for (index_g=0;index_g<ptsz->galaxy_samples_list_num;index_g++){
     class_alloc(ptsz->cl_ngal_ngal_1h[index_g],sizeof(double **)*ptsz->galaxy_samples_list_num,ptsz->error_message);
     class_alloc(ptsz->cl_ngal_ngal_2h[index_g],sizeof(double **)*ptsz->galaxy_samples_list_num,ptsz->error_message);
     class_alloc(ptsz->cl_ngal_ngal_hf[index_g],sizeof(double **)*ptsz->galaxy_samples_list_num,ptsz->error_message);

     class_alloc(ptsz->cl_ngal_lens_1h[index_g],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
     class_alloc(ptsz->cl_ngal_lens_2h[index_g],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
     class_alloc(ptsz->cl_ngal_lens_hf[index_g],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);


                  for (index_l=0;index_l<ptsz->nlSZ;index_l++){
                       ptsz->cl_ngal_lens_1h[index_g][index_l] = 0.;
                       ptsz->cl_ngal_lens_2h[index_g][index_l] = 0.;
                       ptsz->cl_ngal_lens_hf[index_g][index_l] = 0.;
                         // ptsz->cl_cib_cib_2h[index_g][index_g_prime][index_l] = 0.;
                       }
     // class_alloc(ptsz->cl_cib_cib_2h[index_g],sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
     // ptsz->cib_shotnoise[index_g] = 0.;
        for (index_g_prime=0;index_g_prime<ptsz->galaxy_samples_list_num;index_g_prime++){
          class_alloc(ptsz->cl_ngal_ngal_1h[index_g][index_g_prime],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
          class_alloc(ptsz->cl_ngal_ngal_2h[index_g][index_g_prime],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
          class_alloc(ptsz->cl_ngal_ngal_hf[index_g][index_g_prime],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
          // class_alloc(ptsz->cl_cib_cib_2h[index_g][index_g_prime],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
                  for (index_l=0;index_l<ptsz->nlSZ;index_l++){
                     ptsz->cl_ngal_ngal_1h[index_g][index_g_prime][index_l] = 0.;
                     ptsz->cl_ngal_ngal_2h[index_g][index_g_prime][index_l] = 0.;
                     ptsz->cl_ngal_ngal_hf[index_g][index_g_prime][index_l] = 0.;
                     // ptsz->cl_cib_cib_2h[index_g][index_g_prime][index_l] = 0.;
                   }
   }
   }
 }



   class_alloc(ptsz->cl_tSZ_cib_1h,sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_tSZ_cib_2h,sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_lens_cib_1h,sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_lens_cib_2h,sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_gal_cib_1h,sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_gal_cib_2h,sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
   for (index_nu=0;index_nu<ptsz->cib_frequency_list_num;index_nu++){
     class_alloc(ptsz->cl_tSZ_cib_1h[index_nu],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
     class_alloc(ptsz->cl_tSZ_cib_2h[index_nu],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
     class_alloc(ptsz->cl_lens_cib_1h[index_nu],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
     class_alloc(ptsz->cl_lens_cib_2h[index_nu],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
     class_alloc(ptsz->cl_gal_cib_1h[index_nu],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
     class_alloc(ptsz->cl_gal_cib_2h[index_nu],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);

for (index_l=0;index_l<ptsz->nlSZ;index_l++){
  ptsz->cl_tSZ_cib_1h[index_nu][index_l] = 0.;
  ptsz->cl_tSZ_cib_2h[index_nu][index_l] = 0.;
  ptsz->cl_lens_cib_1h[index_nu][index_l] = 0.;
  ptsz->cl_lens_cib_2h[index_nu][index_l] = 0.;
  ptsz->cl_gal_cib_1h[index_nu][index_l] = 0.;
  ptsz->cl_gal_cib_2h[index_nu][index_l] = 0.;

}
   }


 // printf("counting integrands\n");
   int last_index_integrand_id = 0;

   ptsz->index_integrand_id_dndlnM_first  = 0;
   ptsz->index_integrand_id_dndlnM_last = ptsz->index_integrand_id_dndlnM_first + ptsz->N_redshift_dndlnM*ptsz->N_mass_dndlnM - 1;
   ptsz->index_integrand_id_hmf = ptsz->index_integrand_id_dndlnM_last + 1;
   ptsz->index_integrand_id_mean_y = ptsz->index_integrand_id_hmf + 1;
   ptsz->index_integrand_id_cib_monopole_first = ptsz->index_integrand_id_mean_y + 1;
   ptsz->index_integrand_id_cib_monopole_last = ptsz->index_integrand_id_cib_monopole_first + ptsz->n_frequencies_for_cib - 1;
   ptsz->index_integrand_id_cib_shotnoise_first = ptsz->index_integrand_id_cib_monopole_last + 1;
   ptsz->index_integrand_id_cib_shotnoise_last = ptsz->index_integrand_id_cib_shotnoise_first + ptsz->cib_frequency_list_num - 1;
   ptsz->index_integrand_id_sz_ps_first = ptsz->index_integrand_id_cib_shotnoise_last + 1;
   ptsz->index_integrand_id_sz_ps_last = ptsz->index_integrand_id_sz_ps_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_trispectrum_first = ptsz->index_integrand_id_sz_ps_last + 1;
   ptsz->index_integrand_id_trispectrum_last = ptsz->index_integrand_id_trispectrum_first + ptsz->nlSZ*(ptsz->nlSZ+1)/2 - 1;
   ptsz->index_integrand_id_sz_ps_2halo_first = ptsz->index_integrand_id_trispectrum_last +1;
   ptsz->index_integrand_id_sz_ps_2halo_last = ptsz->index_integrand_id_sz_ps_2halo_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_sz_ps_te_y_y_first = ptsz->index_integrand_id_sz_ps_2halo_last +1;
   ptsz->index_integrand_id_sz_ps_te_y_y_last = ptsz->index_integrand_id_sz_ps_te_y_y_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_sz_ps_m_y_y_1h_first = ptsz->index_integrand_id_sz_ps_te_y_y_last +1;
   ptsz->index_integrand_id_sz_ps_m_y_y_1h_last = ptsz->index_integrand_id_sz_ps_m_y_y_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_sz_ps_m_y_y_2h_first = ptsz->index_integrand_id_sz_ps_m_y_y_1h_last +1;
   ptsz->index_integrand_id_sz_ps_m_y_y_2h_last = ptsz->index_integrand_id_sz_ps_m_y_y_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_cov_Y_N_first = ptsz->index_integrand_id_sz_ps_m_y_y_2h_last + 1;
   ptsz->index_integrand_id_cov_Y_N_last = ptsz->index_integrand_id_cov_Y_N_first + ptsz->nlSZ*ptsz->nbins_M - 1;
   ptsz->index_integrand_id_cov_N_N_hsv_first = ptsz->index_integrand_id_cov_Y_N_last + 1;
   ptsz->index_integrand_id_cov_N_N_hsv_last = ptsz->index_integrand_id_cov_N_N_hsv_first +ptsz->nbins_M*(ptsz->nbins_M+1)/2 - 1;
   ptsz->index_integrand_id_cov_Y_Y_ssc_first = ptsz->index_integrand_id_cov_N_N_hsv_last + 1;
   ptsz->index_integrand_id_cov_Y_Y_ssc_last = ptsz->index_integrand_id_cov_Y_Y_ssc_first +ptsz->nlSZ*(ptsz->nlSZ+1)/2 - 1;
   ptsz->index_integrand_id_cov_N_N_first = ptsz->index_integrand_id_cov_Y_Y_ssc_last + 1;
   ptsz->index_integrand_id_cov_N_N_last = ptsz->index_integrand_id_cov_N_N_first + ptsz->nbins_M - 1;
   ptsz->index_integrand_id_cov_Y_N_next_order_first = ptsz->index_integrand_id_cov_N_N_last + 1;
   ptsz->index_integrand_id_cov_Y_N_next_order_last = ptsz->index_integrand_id_cov_Y_N_next_order_first + ptsz->nlSZ*ptsz->nbins_M - 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_1h_first = ptsz->index_integrand_id_cov_Y_N_next_order_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_1h_last = ptsz->index_integrand_id_kSZ_kSZ_gal_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_1h_fft_first = ptsz->index_integrand_id_kSZ_kSZ_gal_1h_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_1h_fft_last = ptsz->index_integrand_id_kSZ_kSZ_gal_1h_fft_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_2h_fft_first = ptsz->index_integrand_id_kSZ_kSZ_gal_1h_fft_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_2h_fft_last = ptsz->index_integrand_id_kSZ_kSZ_gal_2h_fft_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_3h_fft_first = ptsz->index_integrand_id_kSZ_kSZ_gal_2h_fft_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_3h_fft_last = ptsz->index_integrand_id_kSZ_kSZ_gal_3h_fft_first + ptsz->nlSZ - 1;

   ptsz->index_integrand_id_kSZ_kSZ_gallens_1h_fft_first = ptsz->index_integrand_id_kSZ_kSZ_gal_3h_fft_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_gallens_1h_fft_last = ptsz->index_integrand_id_kSZ_kSZ_gallens_1h_fft_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_kSZ_kSZ_gallens_2h_fft_first = ptsz->index_integrand_id_kSZ_kSZ_gallens_1h_fft_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_gallens_2h_fft_last = ptsz->index_integrand_id_kSZ_kSZ_gallens_2h_fft_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_kSZ_kSZ_gallens_3h_fft_first = ptsz->index_integrand_id_kSZ_kSZ_gallens_2h_fft_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_gallens_3h_fft_last = ptsz->index_integrand_id_kSZ_kSZ_gallens_3h_fft_first + ptsz->nlSZ - 1;

   ptsz->index_integrand_id_kSZ_kSZ_lens_1h_fft_first = ptsz->index_integrand_id_kSZ_kSZ_gallens_3h_fft_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_lens_1h_fft_last = ptsz->index_integrand_id_kSZ_kSZ_lens_1h_fft_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_kSZ_kSZ_lens_2h_fft_first = ptsz->index_integrand_id_kSZ_kSZ_lens_1h_fft_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_lens_2h_fft_last = ptsz->index_integrand_id_kSZ_kSZ_lens_2h_fft_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_kSZ_kSZ_lens_3h_fft_first = ptsz->index_integrand_id_kSZ_kSZ_lens_2h_fft_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_lens_3h_fft_last = ptsz->index_integrand_id_kSZ_kSZ_lens_3h_fft_first + ptsz->nlSZ - 1;


   ptsz->index_integrand_id_kSZ_kSZ_gal_2h_first = ptsz->index_integrand_id_kSZ_kSZ_lens_3h_fft_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_2h_last = ptsz->index_integrand_id_kSZ_kSZ_gal_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_3h_first = ptsz->index_integrand_id_kSZ_kSZ_gal_2h_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_3h_last = ptsz->index_integrand_id_kSZ_kSZ_gal_3h_first + ptsz->nlSZ - 1;

   ptsz->index_integrand_id_kSZ_kSZ_gallens_hf_first = ptsz->index_integrand_id_kSZ_kSZ_gal_3h_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_gallens_hf_last = ptsz->index_integrand_id_kSZ_kSZ_gallens_hf_first + ptsz->nlSZ - 1;

   ptsz->index_integrand_id_kSZ_kSZ_lens_hf_first =  ptsz->index_integrand_id_kSZ_kSZ_gallens_hf_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_lens_hf_last =  ptsz->index_integrand_id_kSZ_kSZ_lens_hf_first +  ptsz->nlSZ - 1;

   ptsz->index_integrand_id_kSZ_kSZ_gal_hf_first = ptsz->index_integrand_id_kSZ_kSZ_lens_hf_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_hf_last = ptsz->index_integrand_id_kSZ_kSZ_gal_hf_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_first =  ptsz->index_integrand_id_kSZ_kSZ_gal_hf_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_last = ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_first + ptsz->nlSZ - 1;

   ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_first = ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_last + 1;
   ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_last = ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_tSZ_tSZ_tSZ_2h_first = ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_last + 1;
   ptsz->index_integrand_id_tSZ_tSZ_tSZ_2h_last = ptsz->index_integrand_id_tSZ_tSZ_tSZ_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_tSZ_tSZ_tSZ_3h_first = ptsz->index_integrand_id_tSZ_tSZ_tSZ_2h_last + 1;
   ptsz->index_integrand_id_tSZ_tSZ_tSZ_3h_last = ptsz->index_integrand_id_tSZ_tSZ_tSZ_3h_first + ptsz->nlSZ - 1;

   ptsz->index_integrand_id_tSZ_gal_1h_first = ptsz->index_integrand_id_tSZ_tSZ_tSZ_3h_last + 1;
   ptsz->index_integrand_id_tSZ_gal_1h_last = ptsz->index_integrand_id_tSZ_gal_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_tSZ_gal_2h_first = ptsz->index_integrand_id_tSZ_gal_1h_last + 1;
   ptsz->index_integrand_id_tSZ_gal_2h_last = ptsz->index_integrand_id_tSZ_gal_2h_first + ptsz->nlSZ - 1;

   last_index_integrand_id = ptsz->index_integrand_id_tSZ_gal_2h_last;
   if(ptsz->has_pk_at_z_1h+ptsz->has_pk_at_z_2h >= _TRUE_){
   ptsz->index_integrand_id_pk_at_z_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_pk_at_z_1h_last = ptsz->index_integrand_id_pk_at_z_1h_first + ptsz->n_k_for_pk_hm - 1;
   ptsz->index_integrand_id_pk_at_z_2h_first = ptsz->index_integrand_id_pk_at_z_1h_last + 1;
   ptsz->index_integrand_id_pk_at_z_2h_last = ptsz->index_integrand_id_pk_at_z_2h_first + ptsz->n_k_for_pk_hm - 1;
   last_index_integrand_id =  ptsz->index_integrand_id_pk_at_z_2h_last;
 }

   if(ptsz->has_sz_rates >= _TRUE_){
   ptsz->index_integrand_id_szrates_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_szrates_last = ptsz->index_integrand_id_szrates_first + ptsz->szcat_size - 1;
   last_index_integrand_id =  ptsz->index_integrand_id_szrates_last;
 }

   if(ptsz->has_pk_gg_at_z_1h+ptsz->has_pk_gg_at_z_2h >= _TRUE_){
   ptsz->index_integrand_id_pk_gg_at_z_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_pk_gg_at_z_1h_last = ptsz->index_integrand_id_pk_gg_at_z_1h_first + ptsz->n_k_for_pk_hm - 1;
   ptsz->index_integrand_id_pk_gg_at_z_2h_first = ptsz->index_integrand_id_pk_gg_at_z_1h_last + 1;
   ptsz->index_integrand_id_pk_gg_at_z_2h_last = ptsz->index_integrand_id_pk_gg_at_z_2h_first + ptsz->n_k_for_pk_hm - 1;
   last_index_integrand_id =  ptsz->index_integrand_id_pk_gg_at_z_2h_last;
 }

   if(ptsz->has_pk_bb_at_z_1h+ptsz->has_pk_bb_at_z_2h >= _TRUE_){
   ptsz->index_integrand_id_pk_bb_at_z_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_pk_bb_at_z_1h_last = ptsz->index_integrand_id_pk_bb_at_z_1h_first + ptsz->n_k_for_pk_hm - 1;
   ptsz->index_integrand_id_pk_bb_at_z_2h_first = ptsz->index_integrand_id_pk_bb_at_z_1h_last + 1;
   ptsz->index_integrand_id_pk_bb_at_z_2h_last = ptsz->index_integrand_id_pk_bb_at_z_2h_first + ptsz->n_k_for_pk_hm - 1;
   last_index_integrand_id =  ptsz->index_integrand_id_pk_bb_at_z_2h_last;
 }

 if(ptsz->has_pk_b_at_z_2h >= _TRUE_){
 ptsz->index_integrand_id_pk_b_at_z_2h_first = last_index_integrand_id + 1;
 ptsz->index_integrand_id_pk_b_at_z_2h_last = ptsz->index_integrand_id_pk_b_at_z_2h_first + ptsz->n_k_for_pk_hm - 1;
 last_index_integrand_id =  ptsz->index_integrand_id_pk_b_at_z_2h_last;
}


   if(ptsz->has_pk_em_at_z_1h+ptsz->has_pk_em_at_z_2h >= _TRUE_){
   ptsz->index_integrand_id_pk_em_at_z_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_pk_em_at_z_1h_last = ptsz->index_integrand_id_pk_em_at_z_1h_first + ptsz->n_k_for_pk_hm - 1;
   ptsz->index_integrand_id_pk_em_at_z_2h_first = ptsz->index_integrand_id_pk_em_at_z_1h_last + 1;
   ptsz->index_integrand_id_pk_em_at_z_2h_last = ptsz->index_integrand_id_pk_em_at_z_2h_first + ptsz->n_k_for_pk_hm - 1;
   last_index_integrand_id =  ptsz->index_integrand_id_pk_em_at_z_2h_last;
 }

   if(ptsz->has_pk_HI_at_z_1h+ptsz->has_pk_HI_at_z_2h >= _TRUE_){
   ptsz->index_integrand_id_pk_HI_at_z_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_pk_HI_at_z_1h_last = ptsz->index_integrand_id_pk_HI_at_z_1h_first + ptsz->n_k_for_pk_hm - 1;
   ptsz->index_integrand_id_pk_HI_at_z_2h_first = ptsz->index_integrand_id_pk_HI_at_z_1h_last + 1;
   ptsz->index_integrand_id_pk_HI_at_z_2h_last = ptsz->index_integrand_id_pk_HI_at_z_2h_first + ptsz->n_k_for_pk_hm - 1;
   last_index_integrand_id =  ptsz->index_integrand_id_pk_HI_at_z_2h_last;
 }

   if(ptsz->has_bk_at_z_1h+ptsz->has_bk_at_z_2h +ptsz->has_bk_at_z_3h >= _TRUE_){
   ptsz->index_integrand_id_bk_at_z_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_bk_at_z_1h_last = ptsz->index_integrand_id_bk_at_z_1h_first + ptsz->n_k_for_pk_hm - 1;
   ptsz->index_integrand_id_bk_at_z_2h_first = ptsz->index_integrand_id_bk_at_z_1h_last + 1;
   ptsz->index_integrand_id_bk_at_z_2h_last = ptsz->index_integrand_id_bk_at_z_2h_first + ptsz->n_k_for_pk_hm - 1;
   ptsz->index_integrand_id_bk_at_z_3h_first = ptsz->index_integrand_id_bk_at_z_2h_last + 1;
   ptsz->index_integrand_id_bk_at_z_3h_last = ptsz->index_integrand_id_bk_at_z_3h_first + ptsz->n_k_for_pk_hm - 1;
   last_index_integrand_id =  ptsz->index_integrand_id_bk_at_z_3h_last;
 }

   if(ptsz->has_kSZ_kSZ_tSZ_1h+ptsz->has_kSZ_kSZ_tSZ_2h +ptsz->has_kSZ_kSZ_tSZ_3h >= _TRUE_){
   ptsz->index_integrand_id_kSZ_kSZ_tSZ_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_kSZ_kSZ_tSZ_1h_last = ptsz->index_integrand_id_kSZ_kSZ_tSZ_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_kSZ_kSZ_tSZ_2h_first = ptsz->index_integrand_id_kSZ_kSZ_tSZ_1h_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_tSZ_2h_last = ptsz->index_integrand_id_kSZ_kSZ_tSZ_2h_first + ptsz->nlSZ  - 1;
   ptsz->index_integrand_id_kSZ_kSZ_tSZ_3h_first = ptsz->index_integrand_id_kSZ_kSZ_tSZ_2h_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_tSZ_3h_last = ptsz->index_integrand_id_kSZ_kSZ_tSZ_3h_first + ptsz->nlSZ  - 1;
   last_index_integrand_id =  ptsz->index_integrand_id_kSZ_kSZ_tSZ_3h_last;
 }

   if(ptsz->has_kSZ_kSZ_1h+ptsz->has_kSZ_kSZ_2h >= _TRUE_){
   ptsz->index_integrand_id_kSZ_kSZ_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_kSZ_kSZ_1h_last = ptsz->index_integrand_id_kSZ_kSZ_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_kSZ_kSZ_2h_first = ptsz->index_integrand_id_kSZ_kSZ_1h_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_2h_last = ptsz->index_integrand_id_kSZ_kSZ_2h_first + ptsz->nlSZ  - 1;
   last_index_integrand_id =  ptsz->index_integrand_id_kSZ_kSZ_2h_last;
 }


   if(ptsz->has_bk_ttg_at_z_1h+ptsz->has_bk_ttg_at_z_2h +ptsz->has_bk_ttg_at_z_3h >= _TRUE_){
   ptsz->index_integrand_id_bk_ttg_at_z_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_bk_ttg_at_z_1h_last = ptsz->index_integrand_id_bk_ttg_at_z_1h_first + ptsz->n_k_for_pk_hm - 1;
   ptsz->index_integrand_id_bk_ttg_at_z_2h_first = ptsz->index_integrand_id_bk_ttg_at_z_1h_last + 1;
   ptsz->index_integrand_id_bk_ttg_at_z_2h_last = ptsz->index_integrand_id_bk_ttg_at_z_2h_first + ptsz->n_k_for_pk_hm - 1;
   ptsz->index_integrand_id_bk_ttg_at_z_3h_first = ptsz->index_integrand_id_bk_ttg_at_z_2h_last + 1;
   ptsz->index_integrand_id_bk_ttg_at_z_3h_last = ptsz->index_integrand_id_bk_ttg_at_z_3h_first + ptsz->n_k_for_pk_hm - 1;
   last_index_integrand_id =  ptsz->index_integrand_id_bk_ttg_at_z_3h_last;
 }

   ptsz->index_integrand_id_gal_gal_hf_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_gal_gal_hf_last = ptsz->index_integrand_id_gal_gal_hf_first + ptsz->nlSZ - 1;
   last_index_integrand_id = ptsz->index_integrand_id_gal_gal_hf_last;

   ptsz->index_integrand_id_gal_lens_hf_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_gal_lens_hf_last = ptsz->index_integrand_id_gal_lens_hf_first + ptsz->nlSZ - 1;
   last_index_integrand_id = ptsz->index_integrand_id_gal_lens_hf_last;

   ptsz->index_integrand_id_gal_lensmag_hf_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_gal_lensmag_hf_last = ptsz->index_integrand_id_gal_lensmag_hf_first + ptsz->nlSZ - 1;
   last_index_integrand_id = ptsz->index_integrand_id_gal_lensmag_hf_last;

   ptsz->index_integrand_id_gal_gallens_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_gal_gallens_1h_last = ptsz->index_integrand_id_gal_gallens_1h_first + ptsz->nlSZ - 1;
   last_index_integrand_id = ptsz->index_integrand_id_gal_gallens_1h_last;

   ptsz->index_integrand_id_gal_gallens_2h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_gal_gallens_2h_last = ptsz->index_integrand_id_gal_gallens_2h_first + ptsz->nlSZ - 1;
   last_index_integrand_id = ptsz->index_integrand_id_gal_gallens_2h_last;

   ptsz->index_integrand_id_lens_lensmag_hf_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_lens_lensmag_hf_last = ptsz->index_integrand_id_lens_lensmag_hf_first + ptsz->nlSZ - 1;
   last_index_integrand_id = ptsz->index_integrand_id_lens_lensmag_hf_last;

   ptsz->index_integrand_id_lensmag_lensmag_hf_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_lensmag_lensmag_hf_last = ptsz->index_integrand_id_lensmag_lensmag_hf_first + ptsz->nlSZ - 1;
   last_index_integrand_id = ptsz->index_integrand_id_lensmag_lensmag_hf_last;

   ptsz->index_integrand_id_lens_lens_hf_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_lens_lens_hf_last = ptsz->index_integrand_id_lens_lens_hf_first + ptsz->nlSZ - 1;
   last_index_integrand_id = ptsz->index_integrand_id_lens_lens_hf_last;


   ptsz->index_integrand_id_gal_gal_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_gal_gal_1h_last = ptsz->index_integrand_id_gal_gal_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_gal_gal_2h_first = ptsz->index_integrand_id_gal_gal_1h_last + 1;
   ptsz->index_integrand_id_gal_gal_2h_last = ptsz->index_integrand_id_gal_gal_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_gal_lens_1h_first = ptsz->index_integrand_id_gal_gal_2h_last + 1;
   ptsz->index_integrand_id_gal_lens_1h_last = ptsz->index_integrand_id_gal_lens_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_gal_lensmag_1h_first = ptsz->index_integrand_id_gal_lens_1h_last + 1;
   ptsz->index_integrand_id_gal_lensmag_1h_last = ptsz->index_integrand_id_gal_lensmag_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_gal_gallens_1h_first = ptsz->index_integrand_id_gal_lensmag_1h_last + 1;
   ptsz->index_integrand_id_gal_gallens_1h_last = ptsz->index_integrand_id_gal_gallens_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_gal_gallens_2h_first = ptsz->index_integrand_id_gal_gallens_1h_last + 1;
   ptsz->index_integrand_id_gal_gallens_2h_last = ptsz->index_integrand_id_gal_gallens_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_gallens_gallens_1h_first = ptsz->index_integrand_id_gal_gallens_2h_last + 1;
   ptsz->index_integrand_id_gallens_gallens_1h_last = ptsz->index_integrand_id_gallens_gallens_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_gallens_gallens_2h_first = ptsz->index_integrand_id_gallens_gallens_1h_last + 1;
   ptsz->index_integrand_id_gallens_gallens_2h_last = ptsz->index_integrand_id_gallens_gallens_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_gallens_lens_1h_first = ptsz->index_integrand_id_gallens_gallens_2h_last + 1;
   ptsz->index_integrand_id_gallens_lens_1h_last = ptsz->index_integrand_id_gallens_lens_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_gallens_lens_2h_first = ptsz->index_integrand_id_gallens_lens_1h_last + 1;
   ptsz->index_integrand_id_gallens_lens_2h_last = ptsz->index_integrand_id_gallens_lens_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_tSZ_lensmag_1h_first = ptsz->index_integrand_id_gallens_lens_2h_last + 1;
   ptsz->index_integrand_id_tSZ_lensmag_1h_last = ptsz->index_integrand_id_tSZ_lensmag_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_lensmag_lensmag_1h_first = ptsz->index_integrand_id_tSZ_lensmag_1h_last + 1;
   ptsz->index_integrand_id_lensmag_lensmag_1h_last = ptsz->index_integrand_id_lensmag_lensmag_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_lens_lensmag_1h_first = ptsz->index_integrand_id_lensmag_lensmag_1h_last + 1;
   ptsz->index_integrand_id_lens_lensmag_1h_last = ptsz->index_integrand_id_lens_lensmag_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_cib_cib_2h_first = ptsz->index_integrand_id_lens_lensmag_1h_last + 1;
   ptsz->index_integrand_id_cib_cib_2h_last = ptsz->index_integrand_id_cib_cib_2h_first + ptsz->cib_dim - 1;
   ptsz->index_integrand_id_tSZ_cib_1h_first = ptsz->index_integrand_id_cib_cib_2h_last + 1;
   ptsz->index_integrand_id_tSZ_cib_1h_last = ptsz->index_integrand_id_tSZ_cib_1h_first + ptsz->nlSZ*ptsz->cib_frequency_list_num - 1;
   ptsz->index_integrand_id_lens_cib_1h_first = ptsz->index_integrand_id_tSZ_cib_1h_last + 1;
   ptsz->index_integrand_id_lens_cib_1h_last = ptsz->index_integrand_id_lens_cib_1h_first + ptsz->nlSZ*ptsz->cib_frequency_list_num - 1;
   ptsz->index_integrand_id_gal_cib_1h_first = ptsz->index_integrand_id_lens_cib_1h_last + 1;
   ptsz->index_integrand_id_gal_cib_1h_last = ptsz->index_integrand_id_gal_cib_1h_first + ptsz->nlSZ*ptsz->cib_frequency_list_num - 1;
   ptsz->index_integrand_id_cib_cib_1h_first = ptsz->index_integrand_id_gal_cib_1h_last + 1;
   ptsz->index_integrand_id_cib_cib_1h_last = ptsz->index_integrand_id_cib_cib_1h_first + ptsz->cib_dim - 1;

   last_index_integrand_id = ptsz->index_integrand_id_cib_cib_1h_last;
   if (ptsz->has_ngal_ngal_1h
     + ptsz->has_ngal_ngal_2h
     + ptsz->has_ngal_ngal_hf
     + ptsz->has_ngal_lens_1h
     + ptsz->has_ngal_lens_2h
     + ptsz->has_ngal_lens_hf
   ){
   ptsz->index_integrand_id_ngal_ngal_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_ngal_ngal_1h_last = ptsz->index_integrand_id_ngal_ngal_1h_first + ptsz->ngal_dim - 1;
   ptsz->index_integrand_id_ngal_ngal_2h_first = ptsz->index_integrand_id_ngal_ngal_1h_last+ 1;
   ptsz->index_integrand_id_ngal_ngal_2h_last = ptsz->index_integrand_id_ngal_ngal_2h_first + ptsz->ngal_dim - 1;
   ptsz->index_integrand_id_ngal_ngal_hf_first = ptsz->index_integrand_id_ngal_ngal_2h_last + 1;
   ptsz->index_integrand_id_ngal_ngal_hf_last = ptsz->index_integrand_id_ngal_ngal_hf_first + ptsz->ngal_dim - 1;

   ptsz->index_integrand_id_ngal_lens_1h_first = ptsz->index_integrand_id_ngal_ngal_hf_first + 1;
   ptsz->index_integrand_id_ngal_lens_1h_last = ptsz->index_integrand_id_ngal_lens_1h_first + ptsz->nlSZ*ptsz->galaxy_samples_list_num - 1;
   ptsz->index_integrand_id_ngal_lens_2h_first = ptsz->index_integrand_id_ngal_lens_1h_last+ 1;
   ptsz->index_integrand_id_ngal_lens_2h_last = ptsz->index_integrand_id_ngal_lens_2h_first + ptsz->nlSZ*ptsz->galaxy_samples_list_num - 1;
   ptsz->index_integrand_id_ngal_lens_hf_first = ptsz->index_integrand_id_ngal_lens_2h_last + 1;
   ptsz->index_integrand_id_ngal_lens_hf_last = ptsz->index_integrand_id_ngal_lens_hf_first + ptsz->nlSZ*ptsz->galaxy_samples_list_num - 1;

   last_index_integrand_id = ptsz->index_integrand_id_ngal_lens_hf_last;
}
   ptsz->index_integrand_id_tSZ_cib_2h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_tSZ_cib_2h_last = ptsz->index_integrand_id_tSZ_cib_2h_first + ptsz->nlSZ*ptsz->cib_frequency_list_num - 1;
   ptsz->index_integrand_id_lens_cib_2h_first = ptsz->index_integrand_id_tSZ_cib_2h_last + 1;
   ptsz->index_integrand_id_lens_cib_2h_last = ptsz->index_integrand_id_lens_cib_2h_first + ptsz->nlSZ*ptsz->cib_frequency_list_num - 1;
   ptsz->index_integrand_id_gal_cib_2h_first = ptsz->index_integrand_id_lens_cib_2h_last + 1;
   ptsz->index_integrand_id_gal_cib_2h_last = ptsz->index_integrand_id_gal_cib_2h_first + ptsz->nlSZ*ptsz->cib_frequency_list_num - 1;
   ptsz->index_integrand_id_gal_lens_2h_first = ptsz->index_integrand_id_gal_cib_2h_last + 1;
   ptsz->index_integrand_id_gal_lens_2h_last = ptsz->index_integrand_id_gal_lens_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_gal_lensmag_2h_first = ptsz->index_integrand_id_gal_lens_2h_last + 1;
   ptsz->index_integrand_id_gal_lensmag_2h_last = ptsz->index_integrand_id_gal_lensmag_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_tSZ_lensmag_2h_first = ptsz->index_integrand_id_gal_lensmag_2h_last + 1;
   ptsz->index_integrand_id_tSZ_lensmag_2h_last = ptsz->index_integrand_id_tSZ_lensmag_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_lensmag_lensmag_2h_first = ptsz->index_integrand_id_tSZ_lensmag_2h_last + 1;
   ptsz->index_integrand_id_lensmag_lensmag_2h_last = ptsz->index_integrand_id_lensmag_lensmag_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_lens_lensmag_2h_first = ptsz->index_integrand_id_lensmag_lensmag_2h_last + 1;
   ptsz->index_integrand_id_lens_lensmag_2h_last = ptsz->index_integrand_id_lens_lensmag_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_lens_lens_1h_first = ptsz->index_integrand_id_lens_lensmag_2h_last + 1;
   ptsz->index_integrand_id_lens_lens_1h_last = ptsz->index_integrand_id_lens_lens_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_lens_lens_2h_first = ptsz->index_integrand_id_lens_lens_1h_last + 1;
   ptsz->index_integrand_id_lens_lens_2h_last = ptsz->index_integrand_id_lens_lens_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_tSZ_lens_1h_first = ptsz->index_integrand_id_lens_lens_2h_last + 1;
   ptsz->index_integrand_id_tSZ_lens_1h_last = ptsz->index_integrand_id_tSZ_lens_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_tSZ_lens_2h_first = ptsz->index_integrand_id_tSZ_lens_1h_last + 1;
   ptsz->index_integrand_id_tSZ_lens_2h_last = ptsz->index_integrand_id_tSZ_lens_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_isw_lens_first = ptsz->index_integrand_id_tSZ_lens_2h_last + 1;
   ptsz->index_integrand_id_isw_lens_last = ptsz->index_integrand_id_isw_lens_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_isw_tsz_first = ptsz->index_integrand_id_isw_lens_last + 1;
   ptsz->index_integrand_id_isw_tsz_last = ptsz->index_integrand_id_isw_tsz_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_isw_auto_first = ptsz->index_integrand_id_isw_tsz_last + 1;
   ptsz->index_integrand_id_isw_auto_last = ptsz->index_integrand_id_isw_auto_first + ptsz->nlSZ - 1;

   ptsz->number_of_integrands =  ptsz->index_integrand_id_isw_auto_last + 1;
 // printf("counting integrands %d\n",ptsz->number_of_integrands );




}

double evaluate_dlnMdeltadlnM(double logM,
                             double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct nonlinear * pnl,
                             struct tszspectrum * ptsz)
                         {
  // printf("dlnM\n");
//! JCH edit: I think Komatsu has forgotten the Jacobian factor dlnMdel/dlnM
//! as discussed in Eq. (5) of Komatsu-Seljak (2002)
//! Approximate via standard three-point finite difference
//! (checked w/ Mathematica implementation -- agrees very well)
double result;
double tol=1.e-6;
double delc,rhoc,omega;
//double mvir;
double mvir1,mvir2,rvir1,rvir2;
double cvir1,cvir2,rs1,rs2,m200d2,m200d1;

delc = pvectsz[ptsz->index_Delta_c];
rhoc = pvectsz[ptsz->index_Rho_crit];
omega = pvecback[pba->index_bg_Omega_m]; // rho_m / rho_crit; with rho_m including all non-rel components (e.g., rho_ncdm - 3.* p_ncdm)

double z = pvectsz[ptsz->index_z];

  mvir1=exp(logM-tol);
  mvir2=exp(logM+tol);
  rvir1=evaluate_rvir_of_mvir(mvir1,delc,rhoc,ptsz);//pow(3.*mvir1/4./_PI_/delc/rhoc,1./3.); //! virial radius, h^-1 Mpc
  rvir2=evaluate_rvir_of_mvir(mvir2,delc,rhoc,ptsz);//pow(3.*mvir2/4./_PI_/delc/rhoc,1./3.); //! virial radius, h^-1 Mpc

  //! JCH edit: from personal communication with Battaglia: in their paper,
  //! they changed Duffy's M_pivot from 2d12 Msun/h to 1d14 Msun/h, so
  //! normalization becomes 5.72 instead of 7.85
  //! N.B.: this is the same as Eq. (D17) of Komatsu et al., arXiv:1001.4538

  cvir1=evaluate_cvir_of_mvir(mvir1,z,ptsz,pba);//5.72*pow(mvir1/1e14,-0.081)/pow(1.+z,0.71);
  cvir2=evaluate_cvir_of_mvir(mvir2,z,ptsz,pba);//5.72*pow(mvir2/1e14,-0.081)/pow(1.+z,0.71);
  rs1=rvir1/cvir1; //! NFW scale radius, h^-1 Mpc
  rs2=rvir2/cvir2; //! NFW scale radius, h^-1 Mpc

  class_call(mVIR_to_mDEL(mvir1,
                       rs1,
                       cvir1,
                       200.*omega*rhoc,
                       &m200d1,
                       ptsz),
                  ptsz->error_message,
                  ptsz->error_message);
  class_call(mVIR_to_mDEL(mvir2,
                       rs2,
                       cvir2,
                       200.*omega*rhoc,
                       &m200d2,
                       ptsz),
                  ptsz->error_message,
                  ptsz->error_message);



  double dlnm200ddlnm = (log(m200d2)-log(m200d1))/2./tol;
  result = dlnm200ddlnm;


return result;
}


double radial_kernel_W_galaxy_lensing_at_z(double z,
                                           double * pvectsz,
                                           struct background * pba,
                                           struct tszspectrum * ptsz){

    evaluate_redshift_int_gallens_sources(pvectsz,ptsz);
    double redshift_int_sources = pvectsz[ptsz->index_W_gallens_sources];
    // double chi = sqrt(V->pvectsz[V->ptsz->index_chi2]);
    double chi = sqrt(pvectsz[ptsz->index_chi2]);
    // pvectsz[ptsz->index_lensing_Sigma_crit] = pow(3.*pow(pba->H0/pba->h,2)/2./ptsz->Rho_crit_0,-1)*pow((1.+z),1.)/(chi*redshift_int_sources);

    pvectsz[ptsz->index_lensing_Sigma_crit] = 1./(redshift_int_sources);

    if (isnan(pvectsz[ptsz->index_lensing_Sigma_crit])||isinf(pvectsz[ptsz->index_lensing_Sigma_crit])){
      printf("%.3e\n",redshift_int_sources);
      printf("nan or inf in sigmacrit\n");
      exit(0);
    }

    return 1./pvectsz[ptsz->index_lensing_Sigma_crit];


                                           }

double radial_kernel_W_galaxy_lensing_magnification_at_z(
                                           double z,
                                           double * pvectsz,
                                           struct background * pba,
                                           struct tszspectrum * ptsz){


    double chi = sqrt(pvectsz[ptsz->index_chi2]);
    // compute kernel for lensing magnification
    // lensing of galaxies

    evaluate_redshift_int_lensmag(pvectsz,ptsz);
    double redshift_int_lensmag = pvectsz[ptsz->index_W_lensmag];
    // double chi = sqrt(V->pvectsz[V->ptsz->index_chi2]);
    // pvectsz[ptsz->index_lensing_Sigma_crit] = pow(3.*pow(pba->H0/pba->h,2)/2./ptsz->Rho_crit_0,-1)*pow((1.+z),1.)/(chi*redshift_int_lensmag);
    pvectsz[ptsz->index_lensing_Sigma_crit] = 1./(redshift_int_lensmag);

    if (isnan(pvectsz[ptsz->index_lensing_Sigma_crit])||isinf(pvectsz[ptsz->index_lensing_Sigma_crit])){
      printf("%.3e\n",redshift_int_lensmag);
      printf("nan or inf in sigmacrit\n");
      exit(0);
    }

    return 1./pvectsz[ptsz->index_lensing_Sigma_crit];


                                           }

double radial_kernel_W_cmb_lensing_at_z(double z,
                                           double * pvectsz,
                                           struct background * pba,
                                           struct tszspectrum * ptsz){

    double chi = sqrt(pvectsz[ptsz->index_chi2]);
    double chi_star =  ptsz->chi_star;  // in Mpc/h
    // double sigma_crit_kappa =  pow(3.*pow(pba->H0/pba->h,2)/2./ptsz->Rho_crit_0,-1)*pow((1.+z),1.)*chi_star/chi/(chi_star-chi);
    double sigma_crit_kappa =  chi_star/(chi_star-chi);

    pvectsz[ptsz->index_lensing_Sigma_crit] = sigma_crit_kappa;

    if (isnan(pvectsz[ptsz->index_lensing_Sigma_crit])||isinf(pvectsz[ptsz->index_lensing_Sigma_crit])){
      printf("%.3e\n",pvectsz[ptsz->index_lensing_Sigma_crit]);
      printf("nan or inf in sigmacrit\n");
      exit(0);
    }

    return 1./pvectsz[ptsz->index_lensing_Sigma_crit];


                                           }

double radial_kernel_W_lensing_at_z(  double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct primordial * ppm,
                                    struct nonlinear * pnl,
                                    struct tszspectrum * ptsz){
double chi = sqrt(pvectsz[ptsz->index_chi2]);
double chi_star =  ptsz->chi_star;  // in Mpc/h
double Omega_m = ptsz->Omega_m_0;
double z = pvectsz[ptsz->index_z];
double H_over_c_in_h_over_Mpc = pvecback[pba->index_bg_H]/pba->h;
int index_l = (int)  pvectsz[ptsz->index_multipole];
double ell = ptsz->ell[index_l];

double result = 3.*pow(Omega_m,1.)*pow(pba->H0/pba->h,2)/2.*(chi_star-chi)*pow(1.+z,1.)/chi/chi_star; // homogeneous to (H0/h)^2 * h/Mpc

return result;
}

double radial_kernel_W_lensing_magnification_at_z(  double * pvecback,
                                                    double * pvectsz,
                                                    struct background * pba,
                                                    struct primordial * ppm,
                                                    struct nonlinear * pnl,
                                                    struct tszspectrum * ptsz){
double chi = sqrt(pvectsz[ptsz->index_chi2]);
double chi_star =  ptsz->chi_star;  // in Mpc/h
double Omega_m = ptsz->Omega_m_0;
double z = pvectsz[ptsz->index_z];
double H_over_c_in_h_over_Mpc = pvecback[pba->index_bg_H]/pba->h;
int index_l = (int)  pvectsz[ptsz->index_multipole];
double ell = ptsz->ell[index_l];

evaluate_redshift_int_lensmag(pvectsz,ptsz);

double redshift_int_lensmag = pvectsz[ptsz->index_W_lensmag];
double result = 3.*pow(Omega_m,1.)*pow(pba->H0/pba->h,2)/2.*pow(1.+z,1.)/chi*redshift_int_lensmag; // homogeneous to (H0/h)^2 * h/Mpc

return result;
}

double delta_ell_lens_at_ell_and_z(  double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct primordial * ppm,
                                    struct nonlinear * pnl,
                                    struct tszspectrum * ptsz){

double chi = sqrt(pvectsz[ptsz->index_chi2]);
double chi_star =  ptsz->chi_star;  // in Mpc/h
double Omega_m = ptsz->Omega_m_0;
double z = pvectsz[ptsz->index_z];

double H_over_c_in_h_over_Mpc = pvecback[pba->index_bg_H]/pba->h;
int index_l = (int)  pvectsz[ptsz->index_multipole];
double ell = ptsz->ell[index_l];

double D = pvecback[pba->index_bg_D];
// note the factor 2
double result = 2.*3.*pow(Omega_m,1.)*pow(pba->H0/pba->h,2)/2.*(chi_star-chi)*pow(1.+z,1.)/chi/chi_star/(pow(ell+0.5,2.))*D; // homogeneous to (H0/h)^2 * h/Mpc

return result;
}

double delta_ell_isw_at_ell_and_z( double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct primordial * ppm,
                                    struct nonlinear * pnl,
                                    struct tszspectrum * ptsz){

// double result;
//
double chi = sqrt(pvectsz[ptsz->index_chi2]);
double chi_star =  ptsz->chi_star;
double Omega_m = ptsz->Omega_m_0;
double c = _c_;
double H0 = 100.*pba->h;
double z = pvectsz[ptsz->index_z];
double dgdz = pvectsz[ptsz->index_dgdz]; // d(D/a)/dz = D(1-f)
double D = pvecback[pba->index_bg_D];
double H_over_c_in_h_over_Mpc = pvecback[pba->index_bg_H]/pba->h;
int index_l = (int)  pvectsz[ptsz->index_multipole];
double ell = ptsz->ell[index_l];



double result = 3.*pow(Omega_m,1.)*pow(pba->H0/pba->h,2)/(pow(ell+0.5,2.))*dgdz*H_over_c_in_h_over_Mpc; // homogeneous to (H0/h)^2 * h/Mpc

return result;
                                }


double radial_kernel_W_galaxy_at_z( double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct tszspectrum * ptsz){

double H_over_c_in_h_over_Mpc = pvecback[pba->index_bg_H]/pba->h;
//unwise & halofit model
// if (ptsz->galaxy_sample==1 && ptsz->use_hod==0){
int index_md = (int) pvectsz[ptsz->index_md];
//unwise case
//f_dndz is fphi, which is f/If or bdNdz/ ∫bdNdz
if ((ptsz->galaxy_sample==1 && index_md == ptsz->index_md_gal_gal_hf)
 || (ptsz->galaxy_sample==1 && index_md == ptsz->index_md_gal_lens_hf)
 || (ptsz->galaxy_sample==1 && index_md == ptsz->index_md_gal_lensmag_hf)
 || (ptsz->galaxy_sample==1 && index_md == ptsz->index_md_lens_lensmag_hf)
 || (ptsz->galaxy_sample==1 && index_md == ptsz->index_md_lensmag_lensmag_hf)
 || (ptsz->galaxy_sample==1 && index_md == ptsz->index_md_kSZ_kSZ_gal_hf  && ptsz->use_fdndz_for_ksz2g_eff == 1)
)
{
evaluate_galaxy_number_counts_fdndz(pvecback,pvectsz,pba,ptsz);
}
else{
evaluate_galaxy_number_counts(pvecback,pvectsz,pba,ptsz);
// printf("ng = %.5e\n",pvectsz[ptsz->index_phi_galaxy_counts]);
}
//evaluate_galaxy_number_counts_fdndz(V->pvecback,V->pvectsz,V->pba,V->ptsz);

//  normalized galaxy redshift distribution in the bin:
double phi_galaxy_at_z = pvectsz[ptsz->index_phi_galaxy_counts];
//double phi_galaxy_at_z = 1.; // BB debug
// H_over_c_in_h_over_Mpc = dz/dChi
// phi_galaxy_at_z = dng/dz normalized
double result = H_over_c_in_h_over_Mpc*phi_galaxy_at_z;






return result;
}


double radial_kernel_W_galaxy_ngal_at_z( int index_g,
                                         double * pvecback,
                                         double * pvectsz,
                                         struct background * pba,
                                         struct tszspectrum * ptsz){

double H_over_c_in_h_over_Mpc = pvecback[pba->index_bg_H]/pba->h;

double z_asked  = pvectsz[ptsz->index_z];
double phig = 0.;

if(z_asked<ptsz->normalized_dndz_ngal_z[index_g][0])
   phig = 1e-100;
else if (z_asked>ptsz->normalized_dndz_ngal_z[index_g][ptsz->normalized_dndz_ngal_size[index_g]-1])
   phig = 1e-100;
else  phig =  pwl_value_1d(ptsz->normalized_dndz_ngal_size[index_g],
                         ptsz->normalized_dndz_ngal_z[index_g],
                         ptsz->normalized_dndz_ngal_phig[index_g],
                         z_asked);

//  normalized galaxy redshift distribution in the bin:
// double phi_galaxy_at_z = pvectsz[ptsz->index_phi_galaxy_counts];
//double phi_galaxy_at_z = 1.; // BB debug
// H_over_c_in_h_over_Mpc = dz/dChi
// phi_galaxy_at_z = dng/dz normalized
double result = H_over_c_in_h_over_Mpc*phig;


return result;
}



double evaluate_galaxy_number_counts( double * pvecback,
                                      double * pvectsz,
                                      struct background * pba,
                                      struct tszspectrum * ptsz){


    double z_asked  = pvectsz[ptsz->index_z];
//     double phig = 0.;
//
//
//   if(z_asked<ptsz->normalized_dndz_z[0])
//      phig = 1e-100;
//   else if (z_asked>ptsz->normalized_dndz_z[ptsz->normalized_dndz_size-1])
//      phig = 1e-100;
// else  phig =  pwl_value_1d(ptsz->normalized_dndz_size,
//                            ptsz->normalized_dndz_z,
//                            ptsz->normalized_dndz_phig,
//                            z_asked);

pvectsz[ptsz->index_phi_galaxy_counts] = get_galaxy_number_counts(z_asked,ptsz);

                                    }




double get_source_galaxy_number_counts(double z,
                                struct tszspectrum * ptsz){


    double z_asked  = z-ptsz->Delta_z_source;
    double phig = 0.;

  if(z_asked<ptsz->normalized_source_dndz_z[0])
     phig = 1e-100;
  else if (z_asked>ptsz->normalized_source_dndz_z[ptsz->normalized_source_dndz_size-1])
     phig = 1e-100;
  else  phig =  pwl_value_1d(ptsz->normalized_source_dndz_size,
                           ptsz->normalized_source_dndz_z,
                           ptsz->normalized_source_dndz_phig,
                           z_asked);

return phig;

                                    }


double get_galaxy_number_counts(double z,
                                struct tszspectrum * ptsz){


    double z_asked  = z-ptsz->Delta_z_lens;
    double phig = 0.;
  //
  // if (ptsz->galaxy_sample == 1)
  // if(z_asked<ptsz->normalized_cosmos_dndz_z[0])
  //    phig = 1e-100;
  // else if (z_asked>ptsz->normalized_cosmos_dndz_z[ptsz->normalized_cosmos_dndz_size-1])
  //    phig = 1e-100;
  // else  phig =  pwl_value_1d(ptsz->normalized_cosmos_dndz_size,
  //                              ptsz->normalized_cosmos_dndz_z,
  //                              ptsz->normalized_cosmos_dndz_phig,
  //                              z_asked);

  if(z_asked<ptsz->normalized_dndz_z[0])
     phig = 1e-100;
  else if (z_asked>ptsz->normalized_dndz_z[ptsz->normalized_dndz_size-1])
     phig = 1e-100;
else  phig =  pwl_value_1d(ptsz->normalized_dndz_size,
                           ptsz->normalized_dndz_z,
                           ptsz->normalized_dndz_phig,
                           z_asked);

return phig;

                                    }



double evaluate_galaxy_number_counts_fdndz( double * pvecback,
                                            double * pvectsz,
                                            struct background * pba,
                                            struct tszspectrum * ptsz){


    double z_asked  = pvectsz[ptsz->index_z];
    double phig = 0.;

  if(z_asked<ptsz->normalized_fdndz_z[0])
     phig = 1e-100;
  else if (z_asked>ptsz->normalized_fdndz_z[ptsz->normalized_fdndz_size-1])
     phig = 1e-100;
  else  phig =  pwl_value_1d(ptsz->normalized_fdndz_size,
                               ptsz->normalized_fdndz_z,
                               ptsz->normalized_fdndz_phig,
                               z_asked);

pvectsz[ptsz->index_phi_galaxy_counts] = phig;

                                    }

double HOD_mean_number_of_central_galaxies(double z,
                                           double M_halo,
                                           double M_min,
                                           double sigma_log10M,
                                           double f_cen,
                                           struct tszspectrum * ptsz,
                                           struct background * pba){
 double result = 0.;
 result = f_cen*0.5*(1.+gsl_sf_erf((log10(M_halo/M_min)/sigma_log10M)));
 return result;
}


double HOD_mean_number_of_satellite_galaxies(double z,
                                             double M_halo,
                                             double Nc_mean,
                                             double M_min,
                                             double alpha_s,
                                             double M1_prime,
                                             struct tszspectrum * ptsz,
                                             struct background * pba){


double result;
if (M_halo>M_min){

result = Nc_mean*pow((M_halo-M_min)/M1_prime,alpha_s);

 }
else {
result = 0.;
}
return result;
}


// int evaluate_cib_profile(double m_delta,
//                          double r_delta,
//                          double c_delta,
//                          double * pvecback,
//                          double * pvectsz,
//                          struct background * pba,
//                          struct tszspectrum * ptsz){
//
// // this implementation uses 200m
// double M_halo = m_delta/pba->h; // convet to Msun
//
// double Lc_nu;
// double Lc_nu_prime;
// double Ls_nu;
// double Ls_nu_prime;
//
// double z = pvectsz[ptsz->index_z];
//
// pvectsz[ptsz->index_multipole_for_truncated_nfw_profile] = pvectsz[ptsz->index_multipole_for_cib_profile];
// double xout = 1.;//ptsz->x_out_truncated_nfw_profile_satellite_galaxies;
// double l = pvectsz[ptsz->index_multipole_for_truncated_nfw_profile];
// double chi = sqrt(pvectsz[ptsz->index_chi2]);
// double k = (l+0.5)/chi;
// // printf("r = %.8e m = %.8e z = %.8e\n",r_delta,m_delta,z);
// double us = evaluate_truncated_nfw_profile(k,r_delta,c_delta,xout,pvectsz,pba,ptsz);
//
// double ug_at_ell;
// double nu;
// double nu_prime;
//
// int index_md = (int) pvectsz[ptsz->index_md];
// int index_nu = (int) pvectsz[ptsz->index_frequency_for_cib_profile];
// int index_nu_prime = (int) pvectsz[ptsz->index_frequency_prime_for_cib_profile];
//
//
// // 2-halo terms
// if (_tSZ_cib_2h_
//   ||_cib_monopole_
//   ||_dcib0dz_
//   ||_lens_cib_2h_
//   ||_gal_cib_2h_
//   ||_tSZ_cib_1h_
//   ||_lens_cib_1h_
//   ||_gal_cib_1h_
//   ||_cib_cib_2h_
//    ){
//
// if (_cib_cib_2h_){
// if (index_nu != index_nu_prime){
// if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1){
// nu = ptsz->cib_frequency_list[index_nu];
// Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz,pba);
// Ls_nu = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu,pba,ptsz);
// }
// else if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2){
// nu = ptsz->cib_frequency_list[index_nu_prime];
// Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz,pba);
// Ls_nu = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu_prime,pba,ptsz);
// }
// }
// else{
// nu = ptsz->cib_frequency_list[index_nu];
// Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz,pba);
// Ls_nu = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu,pba,ptsz);
// }
// }
// else if(_cib_monopole_){
// nu = ptsz->frequencies_for_cib[index_nu];
// Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz,pba);
// Ls_nu = get_L_sat_at_z_M_nu(z,M_halo,nu,ptsz);
// us = 1.;
//
// if (isnan(Ls_nu)){
//   printf("nan at %.5e %.5e %.5e\n",z,M_halo,nu);
//   exit(0);
// }
//
// }
// else if(_dcib0dz_){
// nu =  exp(ptsz->array_dcib0dz_nu[index_nu]);
// Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz,pba);
// Ls_nu = get_L_sat_at_z_M_nu(z,M_halo,nu,ptsz);
// us = 1.;
// }
// // cross terms
// else {
// nu = ptsz->frequencies_for_cib[index_nu];
// Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz,pba);
// Ls_nu = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu,pba,ptsz);
// }
//
// // eq. 13 of MM20
// ug_at_ell  = 1./(4.*_PI_)*(Lc_nu+Ls_nu*us);
//
// }// end 2halo terms and monopole
//
// // 1-halo terms
// else if( _cib_cib_1h_){
//
// nu = ptsz->cib_frequency_list[index_nu];
//
// Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz,pba);
// Ls_nu = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu,pba,ptsz);
//
// nu_prime = ptsz->cib_frequency_list[index_nu_prime];
//
// Lc_nu_prime = Luminosity_of_central_galaxies(z,M_halo,nu_prime,pvectsz,ptsz,pba);
// Ls_nu_prime = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu_prime,pba,ptsz);
// // eq. 15 of MM20
// ug_at_ell  = 1./(4.*_PI_)*sqrt(Ls_nu*Ls_nu_prime*us*us+Lc_nu*Ls_nu_prime*us+Lc_nu_prime*Ls_nu*us);
//
// }
// // need to fix units:
// // from Mpc^2 to (Mpc/h)^2
// ug_at_ell *= pow(pba->h,2.);
//
//
// pvectsz[ptsz->index_cib_profile] = ug_at_ell;
//
// }

int evaluate_cib_profile(double m_delta,
                         double r_delta,
                         double c_delta,
                         double * pvecback,
                         double * pvectsz,
                         struct background * pba,
                         struct tszspectrum * ptsz){


double M_halo = m_delta/pba->h; // convert to Msun

double Lc_nu;
double Lc_nu_prime;
double Ls_nu;
double Ls_nu_prime;

double z = pvectsz[ptsz->index_z];

pvectsz[ptsz->index_multipole_for_truncated_nfw_profile] = pvectsz[ptsz->index_multipole_for_cib_profile];
double xout = 1.;//ptsz->x_out_truncated_nfw_profile_satellite_galaxies;
double l = pvectsz[ptsz->index_multipole_for_truncated_nfw_profile];
double chi = sqrt(pvectsz[ptsz->index_chi2]);
double chi_in_Mpc = chi/pba->h;
double k = (l+0.5)/chi;
// double k = (l)/chi;
// printf("r = %.8e m = %.8e z = %.8e\n",r_delta,m_delta,z);
double us = evaluate_truncated_nfw_profile(z,k,r_delta,c_delta,xout);


double ug_at_ell;
double nu;
double nu_prime;

int index_md = (int) pvectsz[ptsz->index_md];
int index_nu = (int) pvectsz[ptsz->index_frequency_for_cib_profile];
int index_nu_prime = (int) pvectsz[ptsz->index_frequency_prime_for_cib_profile];

double L_nu;
double S_nu;
// 2-halo terms
if (_tSZ_cib_2h_
  ||_cib_monopole_
  ||_cib_shotnoise_
  ||_dcib0dz_
  ||_lens_cib_2h_
  ||_gal_cib_2h_
  ||_tSZ_cib_1h_
  ||_lens_cib_1h_
  ||_gal_cib_1h_
  ||_cib_cib_2h_
   ){

if (_cib_cib_2h_){
if (index_nu != index_nu_prime){
if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1){
nu = ptsz->cib_frequency_list[index_nu];
Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz,pba);
Ls_nu = get_L_sat_at_z_M_nu(z,M_halo,nu,ptsz);

if (ptsz->has_cib_flux_cut == 1){
L_nu = Lc_nu+Ls_nu;
S_nu = L_nu/4./_PI_/(1.+z)/chi_in_Mpc/chi_in_Mpc; // eq. 20 of shang et al
if (S_nu*1e3 > ptsz->cib_Snu_cutoff_list_in_mJy[index_nu]){
  Lc_nu = 0.;
  Ls_nu = 0.;
}
}

}
else if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2){
nu = ptsz->cib_frequency_list[index_nu_prime];
Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz,pba);
Ls_nu = get_L_sat_at_z_M_nu(z,M_halo,nu,ptsz);
if (ptsz->has_cib_flux_cut == 1){
L_nu = Lc_nu+Ls_nu;
S_nu = L_nu/4./_PI_/(1.+z)/chi_in_Mpc/chi_in_Mpc;
if (S_nu*1e3 > ptsz->cib_Snu_cutoff_list_in_mJy[index_nu_prime]){
  Lc_nu = 0.;
  Ls_nu = 0.;
}
}
}
}
else{
nu = ptsz->cib_frequency_list[index_nu];
Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz,pba);
Ls_nu = get_L_sat_at_z_M_nu(z,M_halo,nu,ptsz);
// double Ls_nu_old = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu,pba,ptsz);
// printf("%.8e %.8e\n",Ls_nu,Ls_nu_old);
if (ptsz->has_cib_flux_cut == 1){
L_nu = Lc_nu+Ls_nu;
S_nu = L_nu/4./_PI_/(1.+z)/chi_in_Mpc/chi_in_Mpc;
if (S_nu*1e3 > ptsz->cib_Snu_cutoff_list_in_mJy[index_nu]){
  Lc_nu = 0.;
  Ls_nu = 0.;
}
}

}
}
else if(_cib_monopole_){
nu = ptsz->frequencies_for_cib[index_nu];
// printf("nu = %.8e\n",nu);
// exit(0);
Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu/(1.+ptsz->z_obs_cib),pvectsz,ptsz,pba);
Ls_nu = get_L_sat_at_z_M_nu(z,M_halo,nu/(1.+ptsz->z_obs_cib),ptsz);
us = 1.;

// cant apply fluxcut to monopole in the current setup
// careful with frequencies_for_cib vs cib_frequency_list
// if (ptsz->has_cib_flux_cut == 1){
// L_nu = Lc_nu+Ls_nu;
// S_nu = L_nu/4./_PI_/(1.+z)/chi_in_Mpc/chi_in_Mpc;
// if (S_nu*1e3 > ptsz->cib_Snu_cutoff_list_in_mJy[index_nu]){
//   Lc_nu = 0.;
//   Ls_nu = 0.;
// }
// }

if (isnan(Ls_nu)){
  printf("nan at %.5e %.5e %.5e\n",z,M_halo,nu);
  exit(0);
}

}
else if(_cib_shotnoise_){
nu = ptsz->cib_frequency_list[index_nu];
// nu = ptsz->cib_frequency_list[index_nu];
// printf("nu = %.8e\n",nu);
// exit(0);
Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu/(1.+ptsz->z_obs_cib),pvectsz,ptsz,pba);
Ls_nu = get_L_sat_at_z_M_nu(z,M_halo,nu/(1.+ptsz->z_obs_cib),ptsz);
us = 1.;

double Snu_cent = Lc_nu/4./_PI_/(1.+z)/chi_in_Mpc/chi_in_Mpc;
if (ptsz->has_cib_flux_cut == 1){
if (Snu_cent*1e3 > ptsz->cib_Snu_cutoff_list_in_mJy[index_nu]){
  Snu_cent = 0.;
  Lc_nu = 0.;
  Ls_nu = 0.;
}
}
double Snu_sat = Ls_nu/4./_PI_/(1.+z)/chi_in_Mpc/chi_in_Mpc;
if (ptsz->has_cib_flux_cut == 1){
if (Snu_sat*1e3 > ptsz->cib_Snu_cutoff_list_in_mJy[index_nu]){
  Snu_sat = 0.;
  Lc_nu = 0.;
  Ls_nu = 0.;
}
}

ug_at_ell = pow(1./4./_PI_*Lc_nu,2.)+ pow(1./4./_PI_*Ls_nu,2.);
// need to fix units:
// from Mpc^4 to (Mpc/h)^4
ug_at_ell *= pow(pba->h,4.);


pvectsz[ptsz->index_cib_profile] = ug_at_ell;
return _SUCCESS_;

// if (ptsz->has_cib_flux_cut == 1){
// L_nu = Lc_nu+Ls_nu;
// S_nu = L_nu/4./_PI_/(1.+z)/chi_in_Mpc/chi_in_Mpc;
// if (S_nu*1e3 > ptsz->cib_Snu_cutoff_list_in_mJy[index_nu]){
//   Lc_nu = 0.;
//   Ls_nu = 0.;
// }
// }
//
// if (isnan(Ls_nu)){
//   printf("nan at %.5e %.5e %.5e\n",z,M_halo,nu);
//   exit(0);
// }

// ug_at_ell = S_nu;
// // need to fix units:
// // from Mpc^2 to (Mpc/h)^2
// ug_at_ell *= pow(pba->h,2.);
//
//
// pvectsz[ptsz->index_cib_profile] = ug_at_ell;
// return _SUCCESS_;


}
else if(_dcib0dz_){
nu =  exp(ptsz->array_dcib0dz_nu[index_nu]);
Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu/(1.+ptsz->z_obs_cib),pvectsz,ptsz,pba);
Ls_nu = get_L_sat_at_z_M_nu(z,M_halo,nu/(1.+ptsz->z_obs_cib),ptsz);
us = 1.;
if (ptsz->has_cib_flux_cut == 1){
L_nu = Lc_nu+Ls_nu;
S_nu = L_nu/4./_PI_/(1.+z)/chi_in_Mpc/chi_in_Mpc;
if (S_nu*1e3 > ptsz->cib_Snu_cutoff_list_in_mJy[index_nu]){
  Lc_nu = 0.;
  Ls_nu = 0.;
}
}
}
// cross terms
else {
// nu = ptsz->frequencies_for_cib[index_nu];
nu = ptsz->cib_frequency_list[index_nu];
Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz,pba);
Ls_nu = get_L_sat_at_z_M_nu(z,M_halo,nu,ptsz);
if (ptsz->has_cib_flux_cut == 1){
L_nu = Lc_nu+Ls_nu;
S_nu = L_nu/4./_PI_/(1.+z)/chi_in_Mpc/chi_in_Mpc;
if (S_nu*1e3 > ptsz->cib_Snu_cutoff_list_in_mJy[index_nu]){
  Lc_nu = 0.;
  Ls_nu = 0.;
}
}
}

// eq. 13 of MM20
ug_at_ell  = 1./(4.*_PI_)*(Lc_nu+Ls_nu*us);

}// end 2halo terms and monopole

// 1-halo terms
else if( _cib_cib_1h_){

nu = ptsz->cib_frequency_list[index_nu];

Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz,pba);
Ls_nu = get_L_sat_at_z_M_nu(z,M_halo,nu,ptsz);

if (ptsz->has_cib_flux_cut == 1){
L_nu = Lc_nu+Ls_nu;
S_nu = L_nu/4./_PI_/(1.+z)/chi_in_Mpc/chi_in_Mpc;
if (S_nu*1e3 > ptsz->cib_Snu_cutoff_list_in_mJy[index_nu]){
  Lc_nu = 0.;
  Ls_nu = 0.;
}
}

nu_prime = ptsz->cib_frequency_list[index_nu_prime];

Lc_nu_prime = Luminosity_of_central_galaxies(z,M_halo,nu_prime,pvectsz,ptsz,pba);
Ls_nu_prime = get_L_sat_at_z_M_nu(z,M_halo,nu_prime,ptsz);
if (ptsz->has_cib_flux_cut == 1){
L_nu = Lc_nu+Ls_nu;
S_nu = L_nu/4./_PI_/(1.+z)/chi_in_Mpc/chi_in_Mpc;
if (S_nu*1e3 > ptsz->cib_Snu_cutoff_list_in_mJy[index_nu_prime]){
  Lc_nu = 0.;
  Ls_nu = 0.;
}
}

// eq. 15 of MM20
ug_at_ell  = 1./(4.*_PI_)*sqrt(Ls_nu*Ls_nu_prime*us*us+Lc_nu*Ls_nu_prime*us+Lc_nu_prime*Ls_nu*us);

}

if (ptsz->use_maniyar_cib_model == 0){
// need to fix units:
// from Mpc^2 to (Mpc/h)^2
ug_at_ell *= pow(pba->h,2.);
}


pvectsz[ptsz->index_cib_profile] = ug_at_ell;
return _SUCCESS_;
}



double Luminosity_of_central_galaxies(double z,
                                      double  M_halo,
                                      double nu,
                                      double * pvectsz,
                                      struct tszspectrum * ptsz,
                                      struct background * pba){
double result = 0.;

if (ptsz->use_maniyar_cib_model){
  M_halo = M_halo*(1.-ptsz->maniyar_cib_fsub);
}

double L_gal = evaluate_galaxy_luminosity(z, M_halo, nu, ptsz);
double nc;
// nc = HOD_mean_number_of_central_galaxies(z,M_halo,ptsz->M_min_HOD,ptsz->sigma_log10M_HOD,pvectsz,ptsz,pba);
if (M_halo>=ptsz->M_min_HOD) nc = 1.;
else nc = 0.;

if (ptsz->use_nc_1_for_all_halos_cib_HOD == 1){
  nc = 1.;
}

return result =  nc*L_gal;
                                      }

// Old stuff not used
// double Luminosity_of_satellite_galaxies(double z,
//                                         double M_halo,
//                                         double nu,
//                                         struct tszspectrum * ptsz,
//                                         struct background * pba){
//
// double result = 0.;
// double lnMs_min = log(1e6);
// double lnMs_max = log(1e11);
// double dlnM = (lnMs_max - lnMs_min)/10.;
//
// double L_sat = 0.;
// double L_gal;
// double lnMs = lnMs_min;
// double M_sub;
// double M_host =  M_halo;
// double dNdlnMs;
// while (lnMs<lnMs_max){
// M_sub = exp(lnMs);
// L_gal = evaluate_galaxy_luminosity(z, M_sub, nu, ptsz);
// //printf("Lgal = %.3e\n",L_gal);
//
//
// dNdlnMs = subhalo_hmf_dndlnMs(M_host,M_sub);
// L_sat += L_gal*dNdlnMs;
// lnMs += dlnM;
// }
// result = L_sat;
// //printf("%.3e\n",result);
// return result;
//                                       }


double subhalo_hmf_dndlnMs(double M_host,double M_sub,struct tszspectrum * ptsz){
if (ptsz->SHMF==1){
// Subhalo mass function: Equation 12 of https://iopscience.iop.org/article/10.1088/0004-637X/719/1/88/pdf
return 0.30*pow(M_sub/M_host,-0.7)*exp(-9.9*pow(M_sub/M_host,2.5));
}
else if (ptsz->SHMF==2){
// eq 3.9 in 2001.08787
// F. Jiang and F. C. van den Bosch, Generating merger trees for dark matter haloes: a comparison of
// methods, MNRAS 440 (2014) 193 [1311.5225].
double g1,a1,g2,a2,b,x;
g1 = 0.13;
a1 = -0.83;
g2 =  1.33;
a2 = -0.02;
b = 5.67;
x = 1.19;
double mr = M_sub/M_host;
return (g1*pow(mr,a1)+g2*pow(mr,a2))*exp(-b*pow(mr,x));
}
}

// see https://github.com/abhimaniyar/halomodel_cib_tsz_cibxtsz/blob/master/Cell_cib.py
double maniyar_cib_Mdot(double M, double z, struct tszspectrum * ptsz){
  double om0 = ptsz->Omega_m_0;
  double ol0 = 1.-ptsz->Omega_m_0;
  // M in msun
  //result in Msun/yr
  return 46.1*(1.+1.11*z)*sqrt(om0*pow(1.+z,3.)+ ol0)*pow(M/1e12,1.1);
}
// see https://github.com/abhimaniyar/halomodel_cib_tsz_cibxtsz/blob/master/Cell_cib.py
double maniyar_cib_sfr(double M, double z, struct tszspectrum * ptsz){

double mhdot = maniyar_cib_Mdot(M,z,ptsz);
double f_b = ptsz->Omega0_b/ptsz->Omega_m_0;
double logM = log(M);
double logMeff = log(ptsz->m_eff_cib);
double sigM2;
if (M<ptsz->m_eff_cib){
  sigM2 = ptsz->sigma2_LM_cib;
}
else{
  sigM2 = pow(sqrt(ptsz->sigma2_LM_cib)-ptsz->maniyar_cib_tau*r8_max(0.,ptsz->maniyar_cib_zc-z),2.);
}
double sfrmhdot = ptsz->maniyar_cib_etamax*exp(-pow(logM-logMeff,2.)/(2.*sigM2));
return mhdot*f_b*sfrmhdot/1e-10; // 1e-10 is the Kennicutt constant for Chabrier IMF see  https://arxiv.org/pdf/2006.16329.pdf units: Msun/yr/Lsun
// result is in Lsun
}

double evaluate_Sigma_cib(double M, struct tszspectrum * ptsz){
// eq. 33 MM20
double result;
double log10M = log10(M);
double log10Meff = log10(ptsz->m_eff_cib);
result = M/sqrt(2.*_PI_*ptsz->sigma2_LM_cib)*exp(-pow(log10M-log10Meff,2.)/(2.*ptsz->sigma2_LM_cib));

return result;
                                      }

double evaluate_phi_cib(double z, struct tszspectrum * ptsz){
// eq. 31 of MM20
double result = pow(1.+z,ptsz->delta_cib);
if (z>ptsz->z_plateau_cib)
  result = pow(1.+z,0.);
return result;
                                      }


/// this agrees with FMC but may be different from the Websky model,
/// in websky the normalization seeem to not depend on z.
double evaluate_sed_cib(double z, double nu, struct tszspectrum * ptsz){
double result = 0.;
double Td = evaluate_dust_temperature(z,ptsz);

//nu = nu*(1.+z);

double nu0;
double x = -(3.+ptsz->beta_cib+ptsz->gamma_cib)*exp(-(3.+ptsz->beta_cib+ptsz->gamma_cib));
double hplanck=6.62607004e-34; //#m2 kg / s
double kb = 1.38064852e-23; //#m2 kg s-2 K-1
double clight = 299792458.;
nu0 = 1e-9*kb*Td/hplanck*(3.+ptsz->beta_cib+ptsz->gamma_cib+gsl_sf_lambert_W0(x));

if (ptsz->cib_nu0_norm == 1){
//printf("nu=%.3e, nu0 = %.3e\n",nu,nu0);
if (nu>=nu0){
result = pow(nu/nu0,-ptsz->gamma_cib);
}
else{
double B_nu_at_Td = (2.*hplanck/clight/clight)*pow(nu*1e9,3.)
                    /(exp(hplanck*nu*1e9/kb/Td)-1.);
double B_nu0_at_Td = (2.*hplanck/clight/clight)*pow(nu0*1e9,3.)
                    /(exp(hplanck*nu0*1e9/kb/Td)-1.);
result = pow(nu/nu0,ptsz->beta_cib)*B_nu_at_Td/B_nu0_at_Td;

}
}
else{
  // double z_fid = 0.
  // double Td_fid = evaluate_dust_temperature(z_fid,ptsz);
  // double nu0_fid;
  // // double x = -(3.+ptsz->beta_cib+ptsz->gamma_cib)*exp(-(3.+ptsz->beta_cib+ptsz->gamma_cib));
  // // double hplanck=6.62607004e-34; //#m2 kg / s
  // // double kb = 1.38064852e-23; //#m2 kg s-2 K-1
  // // double clight = 299792458.;
  // nu0 = 1e-9*kb*Td_fid/hplanck*(3.+ptsz->beta_cib+ptsz->gamma_cib+gsl_sf_lambert_W0(x));
  if (nu>=nu0){

  result = pow(nu,-ptsz->gamma_cib);
  }
  else{
    double B_nu_at_Td = (2.*hplanck/clight/clight)*pow(nu*1e9,3.)
                        /(exp(hplanck*nu*1e9/kb/Td)-1.);
    // double B_nu0_at_Td = (2.*hplanck/clight/clight)*pow(nu0*1e9,3.)
    //                     /(exp(hplanck*nu0*1e9/kb/Td)-1.);
    result = pow(nu,ptsz->beta_cib)*B_nu_at_Td/pow(Td,4.+ptsz->beta_cib);

  }

}

//printf("nu0=%.6e\n",nu0);
//exit(0);
//
return result;
                                      }

double evaluate_dust_temperature(double z, struct tszspectrum * ptsz){
double result = 0.;
double T0 = ptsz->T0_cib;
double alpha = ptsz->alpha_cib;

result = T0*pow(1.+z,alpha);

return result;
                                      }

double evaluate_galaxy_luminosity(double z, double M, double nu, struct tszspectrum * ptsz){
double result = 0.;


if (ptsz->use_maniyar_cib_model){
  // printf("M = %.3e z = %.3e nu = %.3e Lgal = %.3e\n",M,z,nu,result);
 result = 4.*_PI_*maniyar_cib_sfr(M,z,ptsz)*get_cib_Snu_z_and_nu(z,nu,ptsz);
 // printf("M = %.3e z = %.3e nu = %.3e Lgal = %.3e\n",M,z,nu,result);
}
else{
double L0 = ptsz->L0_cib;
double Phi = evaluate_phi_cib(z,ptsz);
double Theta =  evaluate_sed_cib(z,nu*(1.+z),ptsz);
double Sigma = evaluate_Sigma_cib(M,ptsz);
result = L0*Phi*Sigma*Theta;
//printf("Phi =  %.3e, Theta = %.3e, Sigma = %.3e\n",Phi, Theta, Sigma);
}




return result;
                                      }


// double get_galaxy_profile_at_z_m_l_1h(double z,
//                                       double m,
//                                       double r_delta,
//                                       double c_delta,
//                                       double l,
//                                       struct tszspectrum * ptsz,
//                                       struct background * pba){
//   double tau;
//   int first_index_back = 0;
//   double * pvecback;
//   class_alloc(pvecback,
//         pba->bg_size*sizeof(double),
//         ptsz->error_message);
//
//   class_call(background_tau_of_z(pba,z,&tau),
//        ptsz->error_message,
//        ptsz->error_message);
//
//   class_call(background_at_tau(pba,
//                          tau,
//                          pba->long_info,
//                          pba->inter_normal,
//                          &first_index_back,
//                          pvecback),
//        ptsz->error_message,
//        ptsz->error_message);
//
//
// double chi = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi
// // double k1 = (l+0.5)/chi;
// // double l = chi*k-0.5;
// double M_halo = m; // Msun_over_h
// double ng_bar = evaluate_mean_galaxy_number_density_at_z(z,ptsz);
// // printf("%.3e\n",ng_bar);
// // double ng_bar = pow(pba->h,-3.);//pvectsz[ptsz->index_mean_galaxy_number_density];
// double nc;
// double ns;
// double us;
// // double z = pvectsz[ptsz->index_z];
//
//
// double M_min;
// double M0;
// double M1_prime;
// double sigma_log10M;
//
//
// M_min = ptsz->M_min_HOD;
// M1_prime = ptsz->M1_prime_HOD;
// sigma_log10M = ptsz->sigma_log10M_HOD;
// M0 = ptsz->M0_HOD;
//
// nc = HOD_mean_number_of_central_galaxies(z,M_halo,M_min,sigma_log10M,ptsz->f_cen_HOD,ptsz,pba);
//
// ns = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,M0,ptsz->alpha_s_HOD,M1_prime,ptsz,pba);
// double xout = ptsz->x_out_truncated_nfw_profile_satellite_galaxies;
// // pvectsz[ptsz->index_multipole_for_truncated_nfw_profile] = pvectsz[ptsz->index_multipole_for_galaxy_profile];
// // printf("nc ns %.3e  %.3e\n",nc,ns);
// double k = (l+0.5)/chi;
// // double delta = 200.*pvecback[pba->index_bg_Omega_m];
// us = get_truncated_nfw_profile_at_z_m_k_xout(z,m,r_delta,c_delta,k,xout,pba,ptsz);
// // printf("us %.3e\n",us);
// double ug_at_ell;
// ug_at_ell  =(1./ng_bar)*sqrt(ns*ns*us*us+2.*ns*us);
//
//
// free(pvecback);
// return ug_at_ell;
//
// }
//

int evaluate_galaxy_profile_2h(
                            double kl,
                            double m_delta,
                            double r_delta,
                            double c_delta,
                            double * pvecback,
                            double * pvectsz,
                            struct background * pba,
                            struct tszspectrum * ptsz){

double M_halo;

M_halo = m_delta; // Msun_over_h
double z = pvectsz[ptsz->index_z];
double ng_bar =  evaluate_mean_galaxy_number_density_at_z(z,ptsz);//pvectsz[ptsz->index_mean_galaxy_number_density];
double nc;
double ns;
double us;



double M_min;
double M0;
double M1_prime;
double sigma_log10M;

M_min = ptsz->M_min_HOD;
M1_prime = ptsz->M1_prime_HOD;
sigma_log10M = ptsz->sigma_log10M_HOD;
M0 = ptsz->M0_HOD;

nc = HOD_mean_number_of_central_galaxies(z,M_halo,M_min,sigma_log10M,ptsz->f_cen_HOD,ptsz,pba);

ns = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,M0,ptsz->alpha_s_HOD,M1_prime,ptsz,pba);
double xout = ptsz->x_out_truncated_nfw_profile_satellite_galaxies;
// pvectsz[ptsz->index_multipole_for_truncated_nfw_profile] = pvectsz[ptsz->index_multipole_for_galaxy_profile];
// double l = pvectsz[ptsz->index_multipole_for_truncated_nfw_profile];
// double chi = sqrt(pvectsz[ptsz->index_chi2]);
// double k = (l+0.5)/chi;
us = evaluate_truncated_nfw_profile(z,kl,r_delta,c_delta,xout);

double ug_at_ell;


ug_at_ell  = (1./ng_bar)*(nc+ns*us);
if (isinf(ug_at_ell) || isnan(ug_at_ell)){
printf("ng_bar = %.3e nc = %.3e ns = %.3e us = %.3e\n",ng_bar,nc,ns,us);
exit(0);
}
// printf("ug_at_ell = %.3e ngb = %.3e nc = %.5e ns = %.5e us = %.5e\n",
//       ug_at_ell,ng_bar,nc,ns,us);

pvectsz[ptsz->index_galaxy_profile] = ug_at_ell;

}



int evaluate_galaxy_profile_1h(
                            double kl,
                            double m_delta,
                            double r_delta,
                            double c_delta,
                            double * pvecback,
                            double * pvectsz,
                            struct background * pba,
                            struct tszspectrum * ptsz){




double M_halo;

M_halo = m_delta; // Msun_over_h
double ng_bar = pvectsz[ptsz->index_mean_galaxy_number_density];
double nc;
double ns;
double us;
double z = pvectsz[ptsz->index_z];


double M_min;
double M0;
double M1_prime;
double sigma_log10M;

M_min = ptsz->M_min_HOD;
M1_prime = ptsz->M1_prime_HOD;
sigma_log10M = ptsz->sigma_log10M_HOD;
M0 = ptsz->M0_HOD;

nc = HOD_mean_number_of_central_galaxies(z,M_halo,M_min,sigma_log10M,ptsz->f_cen_HOD,ptsz,pba);

ns = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,M0,ptsz->alpha_s_HOD,M1_prime,ptsz,pba);
double xout = ptsz->x_out_truncated_nfw_profile_satellite_galaxies;
// pvectsz[ptsz->index_multipole_for_truncated_nfw_profile] = pvectsz[ptsz->index_multipole_for_galaxy_profile];
// double l = pvectsz[ptsz->index_multipole_for_truncated_nfw_profile];
// double chi = sqrt(pvectsz[ptsz->index_chi2]);
// double k = (l+0.5)/chi;
us = evaluate_truncated_nfw_profile(z,kl,r_delta,c_delta,xout);


double ug_at_ell;

ug_at_ell  = (1./ng_bar)*sqrt(ns*ns*us*us+2.*ns*us);

pvectsz[ptsz->index_galaxy_profile] = ug_at_ell;

}





int evaluate_galaxy_profile_ngal(
                            double kl,
                            double m_delta,
                            double r_delta,
                            double c_delta,
                            double * pvecback,
                            double * pvectsz,
                            struct background * pba,
                            struct tszspectrum * ptsz){




double M_halo;

M_halo = m_delta; // Msun_over_h
int index_md = (int) pvectsz[ptsz->index_md];
int index_g = (int) pvectsz[ptsz->index_ngal_for_galaxy_profile];
int index_g_prime;

// printf("_ngal_lens_1h_ = %d\n",_ngal_lens_1h_);
// exit(0);

if (_ngal_lens_1h_
  ||_ngal_lens_2h_){
  index_g_prime = index_g;
}
 else{
   index_g_prime = (int) pvectsz[ptsz->index_ngal_prime_for_galaxy_profile];
 }

double z = pvectsz[ptsz->index_z];
double ng_bar;
double nc;
double ns;
double us;

double M_min;
double M0;
double M1_prime;
double sigma_log10M;
double alpha_s_HOD = ptsz->alpha_s_HOD_ngal[index_g];
double f_cen_HOD = ptsz->f_cen_HOD_ngal[index_g];

double xout = ptsz->x_out_truncated_nfw_profile_satellite_galaxies_ngal[index_g];

M_min = ptsz->M_min_HOD_ngal[index_g];
M1_prime = ptsz->M1_prime_HOD_ngal[index_g];
sigma_log10M = ptsz->sigma_log10M_HOD_ngal[index_g];
M0 = ptsz->M0_HOD_ngal[index_g];

double ng_bar_galprime;
double nc_galprime;
double ns_galprime;
double us_galprime;

double M_min_galprime;
double M0_galprime;
double M1_prime_galprime;
double sigma_log10M_galprime;

double alpha_s_HOD_galprime = ptsz->alpha_s_HOD_ngal[index_g_prime];
double f_cen_HOD_galprime = ptsz->f_cen_HOD_ngal[index_g_prime];

double xout_galprime = ptsz->x_out_truncated_nfw_profile_satellite_galaxies_ngal[index_g_prime];

M_min_galprime = ptsz->M_min_HOD_ngal[index_g_prime];
M1_prime_galprime = ptsz->M1_prime_HOD_ngal[index_g_prime];
sigma_log10M_galprime = ptsz->sigma_log10M_HOD_ngal[index_g_prime];
M0_galprime = ptsz->M0_HOD_ngal[index_g_prime];


// if (index_g_prime == index_g){



// }



double ug_at_ell;

// 1-halo auto terms
if (_ngal_ngal_1h_){
  nc = HOD_mean_number_of_central_galaxies(z,M_halo,M_min,sigma_log10M,f_cen_HOD,ptsz,pba);
  ns = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,M0,alpha_s_HOD,M1_prime,ptsz,pba);
  us = evaluate_truncated_nfw_profile(z,kl,r_delta,c_delta,xout);
  ng_bar = evaluate_mean_galaxy_number_density_at_z_ngal(z,index_g,ptsz);
    if (index_g_prime == index_g){
      ug_at_ell  = (1./ng_bar)*sqrt(ns*ns*us*us+2.*ns*us);
    }
    else{
      // printf("doing cross\n");
      // exit(0);
      ng_bar_galprime = evaluate_mean_galaxy_number_density_at_z_ngal(z,index_g_prime,ptsz);
      nc_galprime = HOD_mean_number_of_central_galaxies(z,M_halo,M_min_galprime,sigma_log10M_galprime,f_cen_HOD_galprime,ptsz,pba);
      ns_galprime = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,M0_galprime,alpha_s_HOD_galprime,M1_prime_galprime,ptsz,pba);
      us_galprime = evaluate_truncated_nfw_profile(z,kl,r_delta,c_delta,xout_galprime);

      ug_at_ell  = sqrt((1./ng_bar)*(nc+ns*us)*(1./ng_bar_galprime)*(nc_galprime+ns_galprime*us_galprime));
    }
    }
if (_ngal_ngal_2h_
  ||_ngal_lens_1h_
  ||_ngal_lens_2h_
  ){
  // if (index_g_prime == index_g){
  //   nc = HOD_mean_number_of_central_galaxies(z,M_halo,M_min,sigma_log10M,f_cen_HOD,ptsz,pba);
  //   ns = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,M0,alpha_s_HOD,M1_prime,ptsz,pba);
  //   us = evaluate_truncated_nfw_profile(z,kl,r_delta,c_delta,xout);
  //   ng_bar = evaluate_mean_galaxy_number_density_at_z_ngal(z,index_g,ptsz);
  //   ug_at_ell  = (1./ng_bar)*(nc+ns*us);
  // }
  // else{
    if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1
       || (int) pvectsz[ptsz->index_part_id_cov_hsv] ==  0) {
    nc = HOD_mean_number_of_central_galaxies(z,M_halo,M_min,sigma_log10M,f_cen_HOD,ptsz,pba);
    ns = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,M0,alpha_s_HOD,M1_prime,ptsz,pba);
    us = evaluate_truncated_nfw_profile(z,kl,r_delta,c_delta,xout);
    ng_bar = evaluate_mean_galaxy_number_density_at_z_ngal(z,index_g,ptsz);
    ug_at_ell  = (1./ng_bar)*(nc+ns*us);

    // printf("ug_at_ell = %.3e ngb = %.3e nc = %.5e ns = %.5e us = %.5e\n",
    //       ug_at_ell,ng_bar,nc,ns,us);

    }
    if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) { // gal_gal_2h for index_g_prime
    ng_bar_galprime = evaluate_mean_galaxy_number_density_at_z_ngal(z,index_g_prime,ptsz);
    nc_galprime = HOD_mean_number_of_central_galaxies(z,M_halo,M_min_galprime,sigma_log10M_galprime,f_cen_HOD_galprime,ptsz,pba);
    ns_galprime = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,M0_galprime,alpha_s_HOD_galprime,M1_prime_galprime,ptsz,pba);
    us_galprime = evaluate_truncated_nfw_profile(z,kl,r_delta,c_delta,xout_galprime);
    ng_bar_galprime = evaluate_mean_galaxy_number_density_at_z_ngal(z,index_g_prime,ptsz);
    ug_at_ell  = (1./ng_bar_galprime)*(nc_galprime+ns_galprime*us_galprime);
    }
  // }
  }

if (_ngal_lens_1h_
  ||_ngal_lens_2h_
  ){
  // if (index_g_prime == index_g){
  //   nc = HOD_mean_number_of_central_galaxies(z,M_halo,M_min,sigma_log10M,f_cen_HOD,ptsz,pba);
  //   ns = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,M0,alpha_s_HOD,M1_prime,ptsz,pba);
  //   us = evaluate_truncated_nfw_profile(z,kl,r_delta,c_delta,xout);
  //   ng_bar = evaluate_mean_galaxy_number_density_at_z_ngal(z,index_g,ptsz);
  //   ug_at_ell  = (1./ng_bar)*(nc+ns*us);
  // }
  // else{
    // if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1
    //    || (int) pvectsz[ptsz->index_part_id_cov_hsv] ==  0) {
    nc = HOD_mean_number_of_central_galaxies(z,M_halo,M_min,sigma_log10M,f_cen_HOD,ptsz,pba);
    ns = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,M0,alpha_s_HOD,M1_prime,ptsz,pba);
    us = evaluate_truncated_nfw_profile(z,kl,r_delta,c_delta,xout);
    ng_bar = evaluate_mean_galaxy_number_density_at_z_ngal(z,index_g,ptsz);
    ug_at_ell  = (1./ng_bar)*(nc+ns*us);

    // printf("ug_at_ell = %.3e ngb = %.3e nc = %.5e ns = %.5e us = %.5e\n",
    //       ug_at_ell,ng_bar,nc,ns,us);

    // }
    // if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) { // gal_gal_2h for index_g_prime
    // ng_bar_galprime = evaluate_mean_galaxy_number_density_at_z_ngal(z,index_g_prime,ptsz);
    // nc_galprime = HOD_mean_number_of_central_galaxies(z,M_halo,M_min_galprime,sigma_log10M_galprime,f_cen_HOD_galprime,ptsz,pba);
    // ns_galprime = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,M0_galprime,alpha_s_HOD_galprime,M1_prime_galprime,ptsz,pba);
    // us_galprime = evaluate_truncated_nfw_profile(z,kl,r_delta,c_delta,xout_galprime);
    // ng_bar_galprime = evaluate_mean_galaxy_number_density_at_z_ngal(z,index_g_prime,ptsz);
    // ug_at_ell  = (1./ng_bar_galprime)*(nc_galprime+ns_galprime*us_galprime);
    // }
  // }
  }

pvectsz[ptsz->index_galaxy_profile] = ug_at_ell;

}




//
// // analytical truncated NFW profile
// double get_truncated_nfw_profile_at_z_m_k_xout(//double * pvecback,
//                                       double z,
//                                       double m,
//                                       double r_delta,
//                                       double c_delta,
//                                       double k,
//                                       double xout,
//                                       // double delta,
//                                       struct background * pba,
//                                       struct tszspectrum * ptsz)
// {
//   // double c_delta, r_delta;
//   double tau;
//   int first_index_back = 0;
//   double * pvecback;
//   class_alloc(pvecback,
//         pba->bg_size*sizeof(double),
//         ptsz->error_message);
//
//   class_call(background_tau_of_z(pba,z,&tau),
//        ptsz->error_message,
//        ptsz->error_message);
//
//   class_call(background_at_tau(pba,
//                          tau,
//                          pba->long_info,
//                          pba->inter_normal,
//                          &first_index_back,
//                          pvecback),
//        ptsz->error_message,
//        ptsz->error_message);
//
//
// double chi = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi
// double rho_crit = (3./(8.*_PI_*_G_*_M_sun_))
//                                 *pow(_Mpc_over_m_,1)
//                                 *pow(_c_,2)
//                                 *pvecback[pba->index_bg_rho_crit]
//                                 /pow(pba->h,2);
//
// double q = k*r_delta/c_delta*(1.+z);
//
// double denominator = m_nfw(c_delta); //normalization
//
//
// double numerator = cos(q)*(gsl_sf_Ci((1.+xout*c_delta)*q)-gsl_sf_Ci(q))
//                    +sin(q)*(gsl_sf_Si((1.+xout*c_delta)*q)-gsl_sf_Si(q))
//                    -sin(xout*c_delta*q)/((1.+xout*c_delta)*q);
//                    // printf("%.3e %.3e\n",numerator, denominator);
//
//
// free(pvecback);
// return numerator/denominator;
// //return 1.; // BB debug
// }



// analytical truncated NFW profile
// truncated at r_out = xout*r_delta
double evaluate_truncated_nfw_profile(//double * pvecback,
                                      double z,
                                      double k,
                                      double r_delta,
                                      double c_delta,
                                      double xout)//, // so: r_out = xout*r_delta
                                      //double * pvectsz,
                                      // struct background * pba,
                                      // struct tszspectrum * ptsz)
{

//double z = pvectsz[ptsz->index_z];
// c_delta = 5.;
double q = k*r_delta/c_delta*(1.+z); // uk -> 1 when q->0
double denominator = m_nfw(xout*c_delta); //normalization



double numerator = cos(q)*(gsl_sf_Ci((1.+xout*c_delta)*q)-gsl_sf_Ci(q))
                   +sin(q)*(gsl_sf_Si((1.+xout*c_delta)*q)-gsl_sf_Si(q))
                   -sin(xout*c_delta*q)/((1.+xout*c_delta)*q);
// printf("%d %.8e\n",ptsz->delta_def_matter_density,numerator/denominator);
if (isnan(numerator/denominator) || isinf(numerator/denominator)){
  printf("evaluate_truncated_nfw_profile: r %.3e c %.3e  k %.3e z %.3e\n",r_delta,c_delta,k,z);
  printf("evaluate_truncated_nfw_profile: num %.3e den %.3e q %.3e x %.3e  k %.3e z %.3e\n",numerator, denominator,q,xout,k,z);
  exit(0);
}

return numerator/denominator;

}

// from https://arxiv.org/pdf/1112.5479.pdf
double get_c200c_at_m_and_z_B13(double M,
                                double z,
                                struct background * pba,
                                struct tszspectrum * ptsz){
  double * pvecback;
  double tau;
  int first_index_back = 0;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);


  class_call(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             pba->error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &first_index_back,
                               pvecback),
             pba->error_message,
             pba->error_message);

double D = pvecback[pba->index_bg_D];
double nu = 1./D*(1.12*pow(M/5e13,0.3)+0.53); //use the nu as defined in the B13 paper and pivot mass in Msun/h
// double c200m  = pow(D,0.9)*7.7*pow(nu,-0.29); // vir
// double nu = sqrt(get_nu_at_z_and_m(z,M,ptsz,pba));
double c200c  = pow(D,0.54)*5.9*pow(nu,-0.35); // 200c
free(pvecback);
return c200c;
}

double get_c200c_at_m_and_z(double M,
                            double z,
                            struct background * pba,
                            struct tszspectrum * ptsz)
                            {
double c;
if (ptsz->concentration_parameter==0){
c = get_c200c_at_m_and_z_D08(M,z);
}
else if (ptsz->concentration_parameter==6){
c = get_c200c_at_m_and_z_B13(M,z,pba,ptsz);
}
else if (ptsz->concentration_parameter==7){
c = 5.;
}
return  c;
                            }


double get_c200m_at_m_and_z(double M,
                            double z,
                            struct background * pba,
                            struct tszspectrum * ptsz)
                            {
double c;
if (ptsz->concentration_parameter==0){
c = get_c200m_at_m_and_z_D08(M,z);
}
else if (ptsz->concentration_parameter==6){
c = get_c200m_at_m_and_z_B13(M,z,pba,ptsz);
}
return  c;
                            }

// from https://arxiv.org/pdf/1112.5479.pdf
double get_c200m_at_m_and_z_B13(double M,
                                double z,
                                struct background * pba,
                                struct tszspectrum * ptsz){
  double * pvecback;
  double tau;
  int first_index_back = 0;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);


  class_call(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             pba->error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &first_index_back,
                               pvecback),
             pba->error_message,
             pba->error_message);

double D = pvecback[pba->index_bg_D];
double nu = 1./D*(1.12*pow(M/5e13,0.3)+0.53); // pivot mass in Msun/h
// double nu = sqrt(get_nu_at_z_and_m(z,M,ptsz,pba));


// double c200m  = pow(D,0.9)*7.7*pow(nu,-0.29); // vir
double c200m  = pow(D,1.15)*9.*pow(nu,-0.29); // 200m
free(pvecback);
return c200m;
}


double get_c200m_at_m_and_z_D08(double M,
                                double z){
// if (ptsz->concentration_parameter==0){
  // double M = pvectsz[ptsz->index_m200m];// mass in  Msun/h
  double A = 10.14;
  double B = -0.081;
  double C = -1.01;
  double M_pivot = 2.e12; // pivot mass in Msun/h

  // double z = pvectsz[ptsz->index_z];
  double c200m =A*pow(M/M_pivot,B)*pow(1.+z,C);
  return c200m;
}




double get_c200c_at_m_and_z_D08(double M,
                             double z)
{
double A = 5.71;
double B = -0.084;
double C = -0.47;
double M_pivot = 2.e12; // pivot mass in Msun/h
double c200c =A*pow(M/M_pivot,B)*pow(1.+z,C);
return c200c;
}




double  get_c500c_at_m_and_z(//double * pvecback,
                        double m,
                        double z,
                        struct background * pba,
                        struct tszspectrum * ptsz)
{
// implemented only KA20 so far...

double A = 3.67;
double B = -0.0903;
double C = -0.51;
double M_pivot = 2.78164e12*pba->h; // pivot mass in Msun/h

double c500 =A*pow(m/M_pivot,B)*pow(1.+z,C);
return c500;
}




double evaluate_unwise_m_min_cut(double z,
                                 int sample_id,
                                 struct tszspectrum * ptsz)
{
//z = 2.1; // BB debug
double m_cut;

double zall_red[4] = {0.75,  1.00,  1.5,  2.0};
double mcut_all_red[4] = {12.00, 12.00, 12.6, 13.6};

double zall_green[11] = {0.00, 0.25, 0.4, 0.5, 0.65, 0.75, 1.00, 1.25, 1.50, 2.00, 2.50};
double mcut_all_green[11] = {11.9, 12.0, 12.15, 12.15, 11.75, 11.75, 12.4, 12.6, 12.75, 13.25, 13.25};

double zall_green_shallow[6] = {0.25, 0.4, 0.5, 0.65, 0.75, 1.00};
double mcut_all_green_shallow[6] = {11.5, 12, 12, 11, 11, 12.72};

// 0 is red
if (sample_id == 0){

if (z<=0.75)
m_cut = 12.;
else if (z>0.75 && z<=2.)
m_cut = pwl_value_1d(4,zall_red,mcut_all_red,z);
else if (z>=2.)
 m_cut = 13.6;

}

// 1 is green
if (sample_id == 1){

if (z<=2.5)
m_cut = pwl_value_1d(11,zall_green,mcut_all_green,z);
else
m_cut = 13.55;

}

// 2 is green_shallow
if (sample_id == 2){

if (z<=1.)
m_cut = pwl_value_1d(6,zall_green_shallow,mcut_all_green_shallow,z);
else if (z < 2.5 )
m_cut = -0.5161*pow(z,4)+2.919*pow(z,3)-5.384*pow(z,2)+3.842*z+12.01;
else
m_cut = 13.42;



}


// 3 is blue
if (sample_id == 3){
m_cut =  11.65 + z;
}

// adjust the mass
double mass_fac = ptsz->M_min_HOD_mass_factor_unwise;

m_cut = pow(10.,m_cut*mass_fac); // as in hmvec 'hod_A_log10mthresh'



// if (sample_id == 0) // red
//   mass_fac = 1.7;
// else if (sample_id == 1) // green
//   mass_fac = 1.2;
// else if (sample_id == 3) // blue
//   mass_fac = 1.;

return m_cut; // all in m_sun/h
}






double integrand_kSZ2_X_at_theta(double ln_ell_prime, void *p){
double ell_prime  = exp(ln_ell_prime);
  //double ell_prime  = ln_ell_prime;

double integrand_cl_kSZ2_X_at_theta = 0.;
struct Parameters_for_integrand_kSZ2_X_at_theta *V = ((struct Parameters_for_integrand_kSZ2_X_at_theta *) p);



     double ell = V->ptsz->ell[V->index_ell_3];
     double abs_ell_minus_ell_prime = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos(V->theta));
     // double abs_ell_minus_ell_prime = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*V->theta);

     double ell_1 = abs_ell_minus_ell_prime;
     double ell_2 = ell_prime;
     double ell_3 = ell;

     // check bispectrum condition
     int bispec_cd;
     bispec_cd = 1.;//bispectrum_condition(ell_1,ell_2,ell_3);
//
     // if (bispec_cd == 1){
//        //
//        // double ln_ell1 = log(ell_1);
//        // double ln_ell2 = log(ell_2);
//        // // double db =  exp(pwl_interp_2d(V->ptsz->N_kSZ2_gal_multipole_grid,
//        // //                                V->ptsz->N_kSZ2_gal_multipole_grid,
//        // //                                V->ln_ell,
//        // //                                V->ln_ell,
//        // //                                V->b_l1_l2_l_1d,
//        // //                                1,
//        // //                                &ln_ell1,
//        // //                                &ln_ell2));
//        // double db =  pwl_interp_2d(V->ptsz->N_kSZ2_gal_multipole_grid,
//        //                                V->ptsz->N_kSZ2_gal_multipole_grid,
//        //                                V->ln_ell,
//        //                                V->ln_ell,
//        //                                V->b_l1_l2_l_1d,
//        //                                1,
//        //                                &ln_ell1,
//        //                                &ln_ell2);
//
       // double theta_1 = cos(V->theta);
       double theta_1 = V->theta;
       double ln_ell2 = log(ell_2);
       // printf("interpolating @ l = %.3e with l_min = %.3e l_max = %.3e\n",
       // ell_2,
       // exp(V->ln_ell[0]),
       // exp(V->ln_ell[V->ptsz->N_kSZ2_gal_multipole_grid-1]));
       // printf("interpolating @ theta = %.3e with theta_min = %.3e theta_max = %.3e\n",
       // theta_1,
       // V->ptsz->theta_kSZ2_gal_theta_grid[0],
       // V->ptsz->theta_kSZ2_gal_theta_grid[V->ptsz->N_kSZ2_gal_theta_grid-1]);


       double db =  pwl_interp_2d(V->ptsz->N_kSZ2_gal_theta_grid,
                                  V->ptsz->N_kSZ2_gal_multipole_grid,
                                  V->ptsz->theta_kSZ2_gal_theta_grid,
                                  V->ln_ell,
                                  V->b_l1_l2_l_1d,
                                  1,
                                  &theta_1,
                                  &ln_ell2);
      // printf("theta = %.8e db = %.8e\n",V->theta,db);
if (isnan(db) || isinf(db)){
  // db = 0.;
if (isnan(db))
printf("found nan in interpolation of b_l1_l2_l_1d\n");
if (isinf(db))
printf("found inf in interpolation of b_l1_l2_l_1d\n");
printf("theta = %.3e \t l2 = %.3e \t l = %.3e\n",theta_1,exp(ln_ell2),ell_3);

printf("\n\n");
exit(0);
}
       // printf("end interpolation db = %.3e\n",db);
      // Alternatively, compute the bispectrum here instaead of interpolation:
      //
      //     V->Pvectsz[V->ptsz->index_ell_1] = ell_1;
      //     V->Pvectsz[V->ptsz->index_ell_2] = ell_2;
      //     V->Pvectsz[V->ptsz->index_ell_3] = ell_3;
      //     V->Pvectsz[V->ptsz->index_md] = V->ptsz->index_md_kSZ_kSZ_gal_hf;
      //
      //     class_call(integrate_over_redshift(V->pba,
      //                                        V->pnl,
      //                                        V->ppm,
      //                                        V->ptsz,
      //                                        V->Pvecback,
      //                                        V->Pvectsz),
      //                     V->ptsz->error_message,
      //                     V->ptsz->error_message);
      // //
      //  double db = V->Pvectsz[V->ptsz->index_integral];
      // printf("db = %.3e\n",db);
      /////////////

      double fl_prime = 1.;
      if  (ell_prime <= V->ptsz->l_unwise_filter[0] || ell_prime >= V->ptsz->l_unwise_filter[V->ptsz->unwise_filter_size-1])
        fl_prime = 0.;
      else
        fl_prime = pwl_value_1d(V->ptsz->unwise_filter_size,
                                V->ptsz->l_unwise_filter,
                                V->ptsz->f_unwise_filter,
                                ell_prime);

      double fl_minus_l_prime = 1.;
      if  (abs_ell_minus_ell_prime <= V->ptsz->l_unwise_filter[0] || abs_ell_minus_ell_prime >= V->ptsz->l_unwise_filter[V->ptsz->unwise_filter_size-1])
        fl_minus_l_prime = 0.;
      else
        fl_minus_l_prime = pwl_value_1d(V->ptsz->unwise_filter_size,
                                        V->ptsz->l_unwise_filter,
                                        V->ptsz->f_unwise_filter,
                                        abs_ell_minus_ell_prime);



       integrand_cl_kSZ2_X_at_theta = fl_minus_l_prime*fl_prime*ell_prime*ell_prime*db/(2.*_PI_)/(2.*_PI_);
       //integrand_cl_kSZ2_X_at_theta = fl_minus_l_prime*fl_prime*ell_prime*db/(2.*_PI_)/(2.*_PI_);


     // }
    //  if (ell>5100.){
    //   printf("ell = %.3e l-l' = %.3e ell_prime = %.3e theta = %.3e \t integrand = %.3e fl = %.3e fl' = %.3e\n",
    //   ell,
    //   abs_ell_minus_ell_prime,
    //   ell_prime,
    //   V->theta,
    //   integrand_cl_kSZ2_X_at_theta,
    //   fl_minus_l_prime,
    //   fl_prime
    // );
  //}
       return integrand_cl_kSZ2_X_at_theta;

}




double integrand_kSZ2_X_lensing_term_at_theta(double ln_ell_prime, void *p){

struct Parameters_for_integrand_kSZ2_X_lensing_term_at_theta *V = ((struct Parameters_for_integrand_kSZ2_X_lensing_term_at_theta *) p);



     double ell = V->ptsz->ell[V->index_ell];
     double ell_prime  = exp(ln_ell_prime);
     double integrand_cl_kSZ2_X_at_theta = 0.;


       double theta = V->theta;

       double db =  pwl_interp_2d(V->ptsz->N_kSZ2_gal_theta_grid,
                                  V->ptsz->N_kSZ2_gal_multipole_grid,
                                  V->ptsz->theta_kSZ2_gal_theta_grid,
                                  V->ln_ellprime,
                                  V->integrand_l_lprime_phi,
                                  1,
                                  &theta,
                                  &ln_ell_prime);
      // printf("theta = %.8e db = %.8e\n",V->theta,db);
if (isnan(db) || isinf(db)){
  // db = 0.;
if (isnan(db))
printf("found nan in interpolation of b_l1_l2_l_1d\n");
if (isinf(db))
printf("found inf in interpolation of b_l1_l2_l_1d\n");
printf("theta = %.3e \t l2 = %.3e \t l = %.3e\n",theta,exp(ln_ell_prime),ell);

printf("\n\n");
exit(0);
}



       integrand_cl_kSZ2_X_at_theta = db;

       return integrand_cl_kSZ2_X_at_theta;

}





double integrand_kSZ2_X(double theta, void *p){
// double ell_prime  = exp(ln_ell_prime);
  //double ell_prime  = ln_ell_prime;

double integrand_cl_kSZ2_X = 0.;
struct Parameters_for_integrand_kSZ2_X *W = ((struct Parameters_for_integrand_kSZ2_X *) p);


  // adaptative integration
   struct Parameters_for_integrand_kSZ2_X_at_theta V;
   V.pnl = W->pnl;
   V.ppm = W->ppm;
   V.ptsz = W->ptsz;
   V.pba = W->pba;
   V.Pvecback = W->Pvecback;
   V.Pvectsz = W->Pvectsz;
   V.ln_ell = W->ln_ell;
   V.index_ell_3 = W->index_ell_3;
   V.b_l1_l2_l_1d = W->b_l1_l2_l_1d;

   //printf("index l3 = %d\n",V.index_ell_3);

   void * params;

   double r; //result of the integral
   double epsrel= 1.e-6;//ptsz->redshift_epsrel;//ptsz->patterson_epsrel;
   double epsabs= 1.e-50;//ptsz->redshift_epsabs;//ptsz->patterson_epsabs;
   int show_neval = 0;//ptsz->patterson_show_neval;

  // while(theta < 2.*_PI_){

    V.theta = theta;
    params = &V;

    double ell_min = exp(V.ptsz->ell_kSZ2_gal_multipole_grid[0]);
    // double ell_min = V.ptsz->l_unwise_filter[0];
    double ell_max = 2.*V.ptsz->l_unwise_filter[V.ptsz->unwise_filter_size-1];

    // printf("ell_min = %.3e \t ell_max = %.3e\n", ell_min,ell_max);
    // exit(0);
    r=Integrate_using_Patterson_adaptive(log(ell_min), log(ell_max),
                                        epsrel, epsabs,
                                        integrand_kSZ2_X_at_theta,
                                        params,show_neval);

// printf("r at theta = %.8e, r = %.8e \n",theta,r);
    integrand_cl_kSZ2_X = r;

    return integrand_cl_kSZ2_X;

}




double integrand_kSZ2_X_lensing_term(double theta, void *p){
// double ell_prime  = exp(ln_ell_prime);
  //double ell_prime  = ln_ell_prime;

double integrand_cl_kSZ2_X = 0.;
struct Parameters_for_integrand_kSZ2_X_lensing_term *W = ((struct Parameters_for_integrand_kSZ2_X_lensing_term *) p);


  // adaptative integration
   struct Parameters_for_integrand_kSZ2_X_lensing_term_at_theta V;
   V.pnl = W->pnl;
   V.ppm = W->ppm;
   V.ptsz = W->ptsz;
   V.pba = W->pba;
   // V.Pvecback = W->Pvecback;
   // V.Pvectsz = W->Pvectsz;
   V.ln_ellprime = W->ln_ellprime;
   V.index_ell = W->index_ell;
   V.integrand_l_lprime_phi = W->integrand_l_lprime_phi;

   //printf("index l3 = %d\n",V.index_ell_3);

   void * params;

   double r; //result of the integral
   double epsrel= 1.e-6;//ptsz->redshift_epsrel;//ptsz->patterson_epsrel;
   double epsabs= 1.e-50;//ptsz->redshift_epsabs;//ptsz->patterson_epsabs;
   int show_neval = 0;//ptsz->patterson_show_neval;

  // while(theta < 2.*_PI_){

    V.theta = theta;
    params = &V;

    double ell_min = exp(V.ptsz->ell_kSZ2_gal_multipole_grid[0]);
    // double ell_min = V.ptsz->l_unwise_filter[0];
    double ell_max = 2.*V.ptsz->l_unwise_filter[V.ptsz->unwise_filter_size-1];

    // printf("ell_min = %.3e \t ell_max = %.3e\n", ell_min,ell_max);
    // exit(0);
    r=Integrate_using_Patterson_adaptive(log(ell_min), log(ell_max),
                                        epsrel, epsabs,
                                        integrand_kSZ2_X_lensing_term_at_theta,
                                        params,show_neval);

// printf("r at theta = %.8e, r = %.8e \n",theta,r);
    integrand_cl_kSZ2_X = r;

    return integrand_cl_kSZ2_X;

}



// eq. 16 of https://arxiv.org/pdf/1510.04075.pdf
double bispectrum_f2_kernel(double k1, double k2, double k3){
// double cos_theta = 1.;
// double term1 = 5./7.;
// double term2 = 1./2.*cos_theta*(k1/k2+k2/k1);
// double term3 = 2./7*pow(cos_theta,2.);
// return term1+term2+term3;


// // Eq. 13 of Cooray & Hu - https://iopscience.iop.org/article/10.1086/318660/fulltext/51716.text.html
// double omega_m = 1.; // for simplicity here, not sure whether it should be omega_m(z)
//
// double term1 = 1.-2./7.*pow(omega_m,-2./63);
// double term2a = pow((k3*k3-k1*k1-k2*k2)/(2.*k1*k2),2.);
// double term2b = (k1*k1+k2*k2)/(k3*k3-k1*k1-k2*k2)+2./7.*pow(omega_m,-2./63);
//
// return term1+term2a*term2b;

// double a1 = 1.;//bispectrum_f2_kernel_eff_a(k1,n1,sig8_at_z,knl);
// double a2 = 1.;//bispectrum_f2_kernel_eff_a(k2,n2,sig8_at_z,knl);
//
// double b1 = 1.;//bispectrum_f2_kernel_eff_b(k1,n1,knl);
// double b2 = 1.;//bispectrum_f2_kernel_eff_b(k2,n2,knl);
//
// double c1 = 1.;//bispectrum_f2_kernel_eff_c(k1,n1,knl);
// double c2 = 1.;//bispectrum_f2_kernel_eff_c(k2,n2,knl);
//
// // Eq. 13 of Cooray & Hu - https://iopscience.iop.org/article/10.1086/318660/fulltext/51716.text.html
// double omega_m = 1.; // for simplicity here, not sure whether it should be omega_m(z)

double term1 = 5./7.;//1.-2./7.*pow(omega_m,-2./63);
double cos_theta_12 = (k3*k3-k1*k1-k2*k2)/(2.*k1*k2);
double term2a = 1./2.*cos_theta_12*(k2/k1+k1/k2);//*b1*b2;
double term2b = cos_theta_12*cos_theta_12*2./7.;//*pow(omega_m,-2./63)*c1*c2;

return term1 + term2a + term2b;//;//;//+term2a+term2b;
// +k3*k3/(4.*k1*k1)
// +k3*k3/(4.*k2*k2)
// -1./4.
// -(k1*k1)/(4.*k2*k2)
// -1./4.
// -(k2*k2)/(4.*k1*k1);
// +2./7.*(k3*k3*k3*k3)/(4.*k1*k1*k2*k2)
// -2.*2./7.*k3*k3/(4.*k2*k2)
// -2.*2./7.*k3*k3/(4.*k1*k1)
// +2./7.*k1*k1/(4.*k2*k2)
// +2./7.*k2*k2/(4.*k1*k1)
// +1./7.;
//
// 10./28.*pk3*pk1 --> ok
//
// +3./7.*k2*k2/(4.*k3*k3)*pk3*pk1 --> ook
//
// +3./7.*k2*k2/(4.*k1*k1)*pk3*pk1 --> ok
//
// -5./7.*(k3*k3)/(4.*k1*k1)*pk3*pk1 --> ok
//
// -5./7.*(k1*k1)/(4.*k3*k3)*pk3*pk1 --> ok
//
// +2./7.*(k2*k2*k2*k2)/(4.*k3*k3*k1*k1)*pk3*pk1

// 10./28.*pk2*pk3 --> ok
//
// +3./7.*k1*k1/(4.*k2*k2)*pk2*pk3 --> ok
//
// +3./7.*k1*k1/(4.*k3*k3)*pk2*pk3 --> ok
//
// -5./7.*(k2*k2)/(4.*k3*k3)*pk2*pk3
//
// -5./7.*(k3*k3)/(4.*k2*k2)*pk2*pk3
//
// +2./7.*(k1*k1*k1*k1)/(4.*k2*k2*k3*k3)*pk2*pk3

//
//
//
// ;
}




// eq. 2.6 of https://arxiv.org/pdf/1111.4477.pdf
// as used in Hill, Ferraro et al 2016 projected field papers
double bispectrum_f2_kernel_eff(double k1, double k2, double k3,
                                double n1, double n2, double sig8_at_z, double knl){
// double cos_theta = 1.;
// double term1 = 5./7.;
// double term2 = 1./2.*cos_theta*(k1/k2+k2/k1);
// double term3 = 2./7*pow(cos_theta,2.);
// return term1+term2+term3;

double a1 = bispectrum_f2_kernel_eff_a(k1,n1,sig8_at_z,knl);
double a2 = bispectrum_f2_kernel_eff_a(k2,n2,sig8_at_z,knl);

double b1 = bispectrum_f2_kernel_eff_b(k1,n1,knl);
double b2 = bispectrum_f2_kernel_eff_b(k2,n2,knl);

double c1 = bispectrum_f2_kernel_eff_c(k1,n1,knl);
double c2 = bispectrum_f2_kernel_eff_c(k2,n2,knl);

// Eq. 13 of Cooray & Hu - https://iopscience.iop.org/article/10.1086/318660/fulltext/51716.text.html
double omega_m = 1.; // for simplicity here, not sure whether it should be omega_m(z)

double term1 = 1.-2./7.*pow(omega_m,-2./63);
double cos_theta_12 = (k3*k3-k1*k1-k2*k2)/(2.*k1*k2);
double term2a = 1./2.*cos_theta_12*(k2/k1+k1/k2)*b1*b2;
double term2b = cos_theta_12*cos_theta_12*2./7.*pow(omega_m,-2./63)*c1*c2;

return term1*a1*a2+term2a+term2b;
}
double bispectrum_f2_kernel_eff_SC(double k1, double k2, double k3,
                                double n1, double n2, double sig8_at_z, double knl){
// double cos_theta = 1.;
// double term1 = 5./7.;
// double term2 = 1./2.*cos_theta*(k1/k2+k2/k1);
// double term3 = 2./7*pow(cos_theta,2.);
// return term1+term2+term3;

double a1 = bispectrum_f2_kernel_eff_a_SC(k1,n1,sig8_at_z,knl);
double a2 = bispectrum_f2_kernel_eff_a_SC(k2,n2,sig8_at_z,knl);

double b1 = bispectrum_f2_kernel_eff_b_SC(k1,n1,knl);
double b2 = bispectrum_f2_kernel_eff_b_SC(k2,n2,knl);

double c1 = bispectrum_f2_kernel_eff_c_SC(k1,n1,knl);
double c2 = bispectrum_f2_kernel_eff_c_SC(k2,n2,knl);

// Eq. 13 of Cooray & Hu - https://iopscience.iop.org/article/10.1086/318660/fulltext/51716.text.html
double omega_m = 1.; // for simplicity here, not sure whether it should be omega_m(z)

double term1 = 5./7.;//1.-2./7.*pow(omega_m,-2./63);
double cos_theta_12 = (k3*k3-k1*k1-k2*k2)/(2.*k1*k2);
double term2a = 1./2.*cos_theta_12*(k2/k1+k1/k2)*b1*b2;
double term2b = cos_theta_12*cos_theta_12*2./7.*c1*c2;//*pow(omega_m,-2./63)*c1*c2;

return term1*a1*a2+term2a+term2b;
}
double bispectrum_f2_kernel_eff_a_SC(double k1,double n1,double sig8_at_z,double knl){
double result;
double a6 = -0.2;
double a1 = 0.25;
double a2 = 3.5;
double q = k1/knl;
double Q3 = bispectrum_f2_kernel_eff_Q3(n1);
double term1 = 1. + pow(sig8_at_z,a6)*pow(0.7*Q3,1./2.)*pow(q*a1,n1+a2);
double term2 = 1.+pow(q*a1,n1+a2);
result = term1/term2;
return result;
}
double bispectrum_f2_kernel_eff_b_SC(double k1,double n1,double knl){
double result;
double a3 = 2.;
double q = k1/knl;
double term1 = 1.+0.2*a3*(n1+3.)*pow(q,n1+3.);
double term2 = 1.+pow(q,n1+3.5);
result = term1/term2;
return result;
}
double bispectrum_f2_kernel_eff_c_SC(double k1,double n1,double knl){
double result;
double a4 = 1.;
double a5 = 2.;
double q = k1/knl;
double term1 = 1. + 4.5*a4/(1.5+pow(n1+3.,4))*pow(q*a5,n1+3.);
double term2 = 1.+pow(q*a5,n1+3.5);
result = term1/term2;
return result;
}
double bispectrum_f2_kernel_eff_a(double k1,double n1,double sig8_at_z,double knl){
double result;
double a6 = -0.575;
double a1 = 0.484;
double a2 = 3.740;
double q = k1/knl;
double Q3 = bispectrum_f2_kernel_eff_Q3(n1);
double term1 = 1. + pow(sig8_at_z,a6)*pow(0.7*Q3,1./2.)*pow(q*a1,n1+a2);
double term2 = 1.+pow(q*a1,n1+a2);
result = term1/term2;
return result;
}
double bispectrum_f2_kernel_eff_b(double k1,double n1,double knl){
double result;
double a3 = -0.849;
double a7 = 0.128;
double a8 = -0.722;
double q = k1/knl;
double term1 = 1.+0.2*a3*(n1+3.)*pow(q*a7,n1+3.+a8);
double term2 = 1.+pow(q*a7,n1+3.5+a8);
result = term1/term2;
return result;
}
double bispectrum_f2_kernel_eff_c(double k1,double n1,double knl){
double result;
double a4 = 0.392;
double a5 = 1.013;
double a9 = -0.926;
double q = k1/knl;
double term1 = 1. + 4.5*a4/(1.5+pow(n1+3.,4))*pow(q*a5,n1+3.+a9);
double term2 = 1.+pow(q*a5,n1+3.5+a9);
result = term1/term2;
return result;
}
double bispectrum_f2_kernel_eff_Q3(double n1){
double result;
result = (4.-pow(2.,n1))/(1.+pow(2.,n1+1.));
return result;
}
