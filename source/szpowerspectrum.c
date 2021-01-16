/** @file input.c Documented SZ module.
 *
 * Boris Bolliet, 2020
 */


#define _DEBUG


#include "szpowerspectrum.h"
#include "sz_tools.h"
#include "Patterson.h"





int szpowerspectrum_init(
                          struct background * pba,
                          struct thermo * pth,
                          struct nonlinear * pnl,
                          struct primordial * ppm,
                          struct tszspectrum * ptsz
			                    )
{

   // Skip the module if no SZ/halo-model computations are requested:
    if (ptsz->has_sz_ps
      + ptsz->has_hmf
      + ptsz->has_pk_at_z_1h
      + ptsz->has_pk_at_z_2h
      + ptsz->has_mean_y
      + ptsz->has_sz_2halo
      + ptsz->has_sz_trispec
      + ptsz->has_sz_m_y_y_1h
      + ptsz->has_sz_m_y_y_2h
      + ptsz->has_sz_te_y_y
      + ptsz->has_sz_cov_N_N
      + ptsz->has_tSZ_tSZ_tSZ_1halo
      + ptsz->has_kSZ_kSZ_gal_1halo
      + ptsz->has_kSZ_kSZ_lensmag_1halo
      + ptsz->has_tSZ_gal_1h
      + ptsz->has_tSZ_gal_2h
      + ptsz->has_tSZ_cib_1h
      + ptsz->has_tSZ_cib_2h
      + ptsz->has_lens_cib_1h
      + ptsz->has_lens_cib_2h
      + ptsz->has_cib_cib_1h
      + ptsz->has_cib_cib_2h
      + ptsz->has_gal_gal_1h
      + ptsz->has_gal_gal_2h
      + ptsz->has_gal_lens_1h
      + ptsz->has_gal_lens_2h
      + ptsz->has_lens_lens_1h
      + ptsz->has_lens_lens_2h
      + ptsz->has_tSZ_lens_1h
      + ptsz->has_tSZ_lens_2h
      + ptsz->has_isw_lens
      + ptsz->has_isw_tsz
      + ptsz->has_isw_auto
      + ptsz->has_dndlnM == _FALSE_)
   {
      if (ptsz->sz_verbose > 0)
         printf("->No SZ-y or N quantities requested. SZ ps module skipped.\n");
   }

   else
   {


   select_multipole_array(ptsz);
   show_preamble_messages(pba,pth,pnl,ppm,ptsz);

   if ((ptsz->experiment == 0 && ptsz->has_completeness_for_ps_SZ == 1) || (ptsz->experiment == 0 && ptsz->has_sz_counts  == 1))
      read_Planck_noise_map(ptsz);


   external_pressure_profile_init(ptsz);

   // rho_nfw loaded only if lens/kSZ requested, see func. def.
   load_rho_nfw_profile(ptsz);

   if (ptsz->concentration_parameter == 4)
      read_Zhao_CM_init(ptsz);

   tabulate_sigma_and_dsigma_from_pk(pba,pnl,ppm,ptsz);


   tabulate_L_sat_at_nu_and_nu_prime(pba,ptsz);

   if (ptsz->has_sigma2_hsv)
   tabulate_sigma2_hsv_from_pk(pba,pnl,ppm,ptsz);



   initialise_and_allocate_memory(ptsz);

   // load alpha(z) normalisation
   // of Tinker et al 2010 HMF
   if (ptsz->MF==1)
   load_T10_alpha_norm(ptsz);


   if (ptsz->has_dndlnM
     // || ptsz->has_tSZ_gal_1h
     // || ptsz->has_tSZ_gal_2h
     // || ptsz->has_kSZ_kSZ_gal_1halo
     // || ptsz->has_gal_gal_1h
     // || ptsz->has_gal_gal_2h
     // || ptsz->has_gal_lens_1h
     // || ptsz->has_gal_lens_2h
   )
   tabulate_dndlnM(pba,pnl,ppm,ptsz);

   if (ptsz->has_vrms2)
   tabulate_vrms2_from_pk(pba,pnl,ppm,ptsz);


      if (ptsz->has_tSZ_gal_1h
       || ptsz->has_tSZ_gal_2h
       || ptsz->has_kSZ_kSZ_gal_1halo
       || ptsz->has_kSZ_kSZ_lensmag_1halo
       || ptsz->has_gal_gal_1h
       || ptsz->has_gal_gal_2h
       || ptsz->has_gal_lens_1h
       || ptsz->has_gal_lens_2h){
   tabulate_mean_galaxy_number_density(pba,pnl,ppm,ptsz);

   // only performed if requested:
   load_normalized_dndz(ptsz);

   if (ptsz->has_kSZ_kSZ_gal_1halo
    || ptsz->has_kSZ_kSZ_lensmag_1halo )
   load_unwise_filter(ptsz);



 }



  // tabulate lensing magnificaion integral, only used when requested
  tabulate_redshift_int_lensmag(ptsz,pba);


   // only performed when requested:
   load_unbinned_nl_yy(ptsz);

   if (ptsz->write_sz>0)
   write_redshift_dependent_quantities(pba,ptsz);

   //SO data and Functions

   if (ptsz->experiment == 1){
      read_SO_Qfit(ptsz);
      read_SO_noise(ptsz);}





   double * Pvecback;
   double * Pvectsz;
   int index_integrand;

   int abort;

#ifdef _OPENMP
   double tstart, tstop;
#endif

   abort = _FALSE_;
    //printf("number_of_integrands=%d\n",ptsz->number_of_integrands);


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
//number_of_threads= 8;
#pragma omp parallel \
   shared(abort,pba,ptsz,ppm,pnl)\
   private(tstart,tstop,Pvectsz,Pvecback,index_integrand)\
   num_threads(number_of_threads)
	 {

#ifdef _OPENMP
	   tstart = omp_get_wtime();
#endif

	   class_alloc_parallel(Pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
       int i;
       for(i = 0; i<ptsz->tsz_size;i++) Pvectsz[i] = 0.;

	   class_alloc_parallel(Pvecback,pba->bg_size*sizeof(double),ptsz->error_message);





//Loop over integrands
//the computation is parallelized with respect to the integrands
#pragma omp for schedule (dynamic)
for (index_integrand=0;index_integrand<ptsz->number_of_integrands;index_integrand++)
	     {
#pragma omp flush(abort)

       Pvectsz[ptsz->index_integrand_id] = index_integrand;
       class_call_parallel(compute_sz(pba,
                                      pnl,
                                      ppm,
                                      ptsz,
                                      Pvecback,
                                      Pvectsz),
                                     ptsz->error_message,
                                     ptsz->error_message);

          }
#ifdef _OPENMP
      tstop = omp_get_wtime();
      if (ptsz->sz_verbose > 0)
         printf("In %s: time spent in parallel region (loop over X's) = %e s for thread %d\n",
                   __func__,tstop-tstart,omp_get_thread_num());


#endif
   free(Pvecback);
   free(Pvectsz);
	} //end of parallel region

   if (abort == _TRUE_) return _FAILURE_;

   ////////////////end - cl



   if (ptsz->write_sz>0 || ptsz->create_ref_trispectrum_for_cobaya){

   write_output_to_files_cl(pnl,pba,ptsz);
   write_output_to_files_ell_indep_ints(pnl,pba,ptsz);
   if (ptsz->sz_verbose>1) show_results(pba,pnl,ppm,ptsz);
 }
   }

   return _SUCCESS_;
}

int szpowerspectrum_free(struct tszspectrum *ptsz)
{
   free(ptsz->ell);
   free(ptsz->cl_sz_1h);
   free(ptsz->cl_isw_lens);
   free(ptsz->cl_isw_tsz);
   free(ptsz->cl_isw_auto);
   free(ptsz->cl_gal_gal_1h);
   free(ptsz->cl_gal_gal_2h);
   free(ptsz->cl_gal_lens_1h);
   free(ptsz->cl_gal_lens_2h);
   free(ptsz->cl_tSZ_gal_1h);
   free(ptsz->cl_tSZ_gal_2h);
   free(ptsz->cl_tSZ_cib_1h);
   free(ptsz->cl_tSZ_cib_2h);
   free(ptsz->cl_lens_cib_1h);
   free(ptsz->cl_lens_cib_2h);
   free(ptsz->cl_cib_cib_1h);
   free(ptsz->cl_cib_cib_2h);
   free(ptsz->cl_lens_lens_1h);
   free(ptsz->cl_lens_lens_2h);
   free(ptsz->cl_tSZ_lens_1h);
   free(ptsz->cl_tSZ_lens_2h);
   free(ptsz->cl_kSZ_kSZ_gal_1h);
   free(ptsz->cl_kSZ_kSZ_lensmag_1h);
   free(ptsz->b_tSZ_tSZ_tSZ_1halo);
   free(ptsz->cl_te_y_y);
   free(ptsz->m_y_y_1h);
   free(ptsz->m_y_y_2h);
   free(ptsz->cov_cl_cl);
   free(ptsz->sig_cl_squared_binned);
   free(ptsz->cl_sz_2h);
   free(ptsz->tllprime_sz);
   free(ptsz->cov_Y_N);
   free(ptsz->cov_Y_N_next_order);
   free(ptsz->cov_Y_Y_ssc);
   free(ptsz->cov_N_N);
   free(ptsz->cov_N_N_hsv);
   free(ptsz->r_Y_N);
   free(ptsz->r_cl_clp);
   free(ptsz->trispectrum_ref);

   free(ptsz->thetas);
   free(ptsz->skyfracs);
   free(ptsz->ylims);
   free(ptsz->sky_averaged_ylims);



if (ptsz->experiment == 1){
    free(ptsz->SO_Qfit);
    free(ptsz->SO_thetas);
    free(ptsz->SO_RMS);
    free(ptsz->SO_skyfrac);
}
   free(ptsz->w_gauss);
   free(ptsz->x_gauss);

if(ptsz->has_kSZ_kSZ_gal_1halo
|| ptsz->has_kSZ_kSZ_lensmag_1halo){
  free(ptsz->ell_kSZ2_gal_multipole_grid);
  free(ptsz->l_unwise_filter);
  free(ptsz->f_unwise_filter);
}

if (ptsz->has_kSZ_kSZ_lensmag_1halo){
  free(ptsz->array_W_lensmag);
  free(ptsz->array_z_W_lensmag);
}

if (ptsz->has_cib_cib_1h
  ||ptsz->has_cib_cib_2h
  ||ptsz->has_tSZ_cib_1h
  ||ptsz->has_tSZ_cib_2h
  ||ptsz->has_lens_cib_1h
  ||ptsz->has_lens_cib_2h
  ){

free(ptsz->array_m_L_sat);
free(ptsz->array_z_L_sat);
free(ptsz->array_L_sat_at_z_and_M_at_nu);
//free(ptsz->array_L_sat_at_z_and_M_at_nu_prime);

  }

if (ptsz->has_dndlnM
  // || ptsz->has_gal_gal_1h
  // || ptsz->has_gal_gal_2h
  // || ptsz->has_gal_lens_1h
  // || ptsz->has_gal_lens_2h
  // || ptsz->has_tSZ_gal_1h
  // || ptsz->has_tSZ_gal_2h
  // || ptsz->has_kSZ_kSZ_gal_1halo
){
   free(ptsz->array_m_dndlnM);
   free(ptsz->array_z_dndlnM);
   free(ptsz->array_dndlnM_at_z_and_M);}

   free(ptsz->dndlnM_at_z_and_M);
   free(ptsz->dndlnM_array_z);
   free(ptsz->dndlnM_array_m);

   free(ptsz->array_radius);
   free(ptsz->array_redshift);
   free(ptsz->array_sigma_at_z_and_R);
   free(ptsz->array_dsigma2dR_at_z_and_R);

   free(ptsz->array_vrms2_at_z);

if (ptsz->has_tSZ_gal_1h
   || ptsz->has_tSZ_gal_2h
   || ptsz->has_kSZ_kSZ_gal_1halo
   || ptsz->has_kSZ_kSZ_lensmag_1halo
   || ptsz->has_gal_gal_1h
   || ptsz->has_gal_gal_2h
   || ptsz->has_gal_lens_1h
   || ptsz->has_gal_lens_2h){
   free(ptsz->array_mean_galaxy_number_density);
   free(ptsz->normalized_dndz_z);
   free(ptsz->normalized_dndz_phig);
}

if (ptsz->include_noise_cov_y_y==1){
   free(ptsz->unbinned_nl_yy_ell);
   free(ptsz->unbinned_nl_yy_n_ell);
}

   free(ptsz->array_sigma2_hsv_at_z);


   free(ptsz->PP_lnx);
   free(ptsz->PP_lnI);
   free(ptsz->PP_d2lnI);

if(ptsz->has_pk_at_z_1h + ptsz->has_pk_at_z_2h >= _TRUE_){

free(ptsz->k_for_pk_hm);
free(ptsz->pk_at_z_1h);
free(ptsz->pk_at_z_2h);
}


 if (ptsz->has_tSZ_lens_1h
  || ptsz->has_tSZ_lens_2h
  || ptsz->has_lens_lens_1h
  || ptsz->has_lens_lens_2h
  || ptsz->has_gal_lens_1h
  || ptsz->has_gal_lens_2h
  || ptsz->has_lens_cib_1h
  || ptsz->has_lens_cib_2h
  || ptsz->has_kSZ_kSZ_gal_1halo){
   free(ptsz->RNFW_lnx);
   free(ptsz->RNFW_lnI);
 }

  if (ptsz->MF==1){
    free(ptsz->T10_ln1pz);
    free(ptsz->T10_lnalpha);
  }

   free(ptsz->CM_redshift);
   free(ptsz->CM_logM);
   free(ptsz->CM_logC);

   free(ptsz->M_bins);
   free(ptsz->cov_Y_N_mass_bin_edges);

   free(ptsz->ln_x_for_pp);
   free(ptsz->x_for_pp);


return _SUCCESS_;
}




int compute_sz(struct background * pba,
                struct nonlinear * pnl,
                struct primordial * ppm,
                struct tszspectrum * ptsz,
                double * Pvecback,
                double * Pvectsz){




   int index_integrand = (int) Pvectsz[ptsz->index_integrand_id];

   //printf("index_integrand = %d\n",index_integrand);
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
   else if (index_integrand == ptsz->index_integrand_id_mean_y && ptsz->has_mean_y) {
      Pvectsz[ptsz->index_md] = ptsz->index_md_mean_y;
      if (ptsz->sz_verbose > 0) printf("computing mean y\n");
   }
   else if (index_integrand>=ptsz->index_integrand_id_sz_ps_first && index_integrand <= ptsz->index_integrand_id_sz_ps_last && ptsz->has_sz_ps){
      Pvectsz[ptsz->index_md] = ptsz->index_md_sz_ps;
      Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_sz_ps_first);
      if (ptsz->sz_verbose == 1 && index_integrand ==ptsz->index_integrand_id_sz_ps_first ) printf("computing cl^yy @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      if (ptsz->sz_verbose == 1 && index_integrand ==ptsz->index_integrand_id_sz_ps_last ) printf("computing cl^yy @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      if (ptsz->sz_verbose >1) printf("computing cl^yy @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
    }

    else if (index_integrand>=ptsz->index_integrand_id_trispectrum_first && index_integrand <= ptsz->index_integrand_id_trispectrum_last && ptsz->has_sz_trispec){
       Pvectsz[ptsz->index_md] = ptsz->index_md_trispectrum;
       int index_ell_ell_prime = (int) (index_integrand - ptsz->index_integrand_id_trispectrum_first);
       int n = (-1.+sqrt(1. + 4.*2.*index_ell_ell_prime))/2.;
       int index_ell = floor(n);
       int index_ell_prime = index_ell_ell_prime -index_ell*(index_ell+1)/2;
       Pvectsz[ptsz->index_multipole] = (double) (index_ell);
       Pvectsz[ptsz->index_multipole_prime] = (double) (index_ell_prime);
       if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_trispectrum_first) printf("computing trispectrum @ ell_id = %.0f, ell_id_prime = %.0f\n",Pvectsz[ptsz->index_multipole],Pvectsz[ptsz->index_multipole_prime]);
       if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_trispectrum_last) printf("computing trispectrum @ ell_id = %.0f, ell_id_prime = %.0f\n",Pvectsz[ptsz->index_multipole],Pvectsz[ptsz->index_multipole_prime]);
       if (ptsz->sz_verbose > 1) printf("computing trispectrum @ ell_id = %.0f, ell_id_prime = %.0f\n",Pvectsz[ptsz->index_multipole],Pvectsz[ptsz->index_multipole_prime]);
     }
   else if (index_integrand>=ptsz->index_integrand_id_sz_ps_2halo_first && index_integrand <= ptsz->index_integrand_id_sz_ps_2halo_last && ptsz->has_sz_2halo){
      Pvectsz[ptsz->index_md] = ptsz->index_md_2halo;
      Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_sz_ps_2halo_first);
      if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_sz_ps_2halo_first) printf("computing cl^yy 2-halo term @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_sz_ps_2halo_last) printf("computing cl^yy 2-halo term @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      if (ptsz->sz_verbose > 1) printf("computing cl^yy 2-halo term @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
    }


   else if (index_integrand>=ptsz->index_integrand_id_sz_ps_te_y_y_first && index_integrand <= ptsz->index_integrand_id_sz_ps_te_y_y_last && ptsz->has_sz_te_y_y){
      Pvectsz[ptsz->index_md] = ptsz->index_md_te_y_y;
      Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_sz_ps_te_y_y_first);
      if (ptsz->sz_verbose > 0) printf("computing cl^Teyy @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
    }

    else if (index_integrand>=ptsz->index_integrand_id_sz_ps_m_y_y_1h_first && index_integrand <= ptsz->index_integrand_id_sz_ps_m_y_y_1h_last && ptsz->has_sz_m_y_y_1h){
       Pvectsz[ptsz->index_md] = ptsz->index_md_m_y_y_1h;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_sz_ps_m_y_y_1h_first);
       if (ptsz->sz_verbose > 0) printf("computing m_y_y_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
     else if (index_integrand>=ptsz->index_integrand_id_sz_ps_m_y_y_2h_first && index_integrand <= ptsz->index_integrand_id_sz_ps_m_y_y_2h_last && ptsz->has_sz_m_y_y_2h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_m_y_y_2h;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_sz_ps_m_y_y_2h_first);
        if (ptsz->sz_verbose > 0) printf("computing m_y_y_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
    else if (index_integrand>=ptsz->index_integrand_id_cov_Y_N_first && index_integrand <= ptsz->index_integrand_id_cov_Y_N_last && ptsz->has_sz_cov_Y_N){
       Pvectsz[ptsz->index_md] = ptsz->index_md_cov_Y_N;
       int index_ell_mass = (int) (index_integrand - ptsz->index_integrand_id_cov_Y_N_first);
       int index_ell = index_ell_mass / ptsz->nbins_M;
       int index_mass = index_ell_mass % ptsz->nbins_M;
       Pvectsz[ptsz->index_multipole] = (double) (index_ell);
       Pvectsz[ptsz->index_mass_bin_1] = (double) (index_mass);
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
         if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_cov_Y_Y_ssc_first) printf("computing cov(Y,Y) [ssc] @ ell_id = %.0f, ell_id = %.0f\n",Pvectsz[ptsz->index_multipole_1],Pvectsz[ptsz->index_multipole_2]);
         if (ptsz->sz_verbose == 1 && index_integrand==ptsz->index_integrand_id_cov_Y_Y_ssc_last) printf("computing cov(Y,Y) [ssc] @ ell_id = %.0f, ell_id = %.0f\n",Pvectsz[ptsz->index_multipole_1],Pvectsz[ptsz->index_multipole_2]);
         if (ptsz->sz_verbose > 1) printf("computing cov(Y,Y) [ssc] @ ell_id = %.0f, ell_id = %.0f\n",Pvectsz[ptsz->index_multipole_1],Pvectsz[ptsz->index_multipole_2]);
       }

    else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_gal_1halo_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_gal_1halo_last && ptsz->has_kSZ_kSZ_gal_1halo){
       Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gal_1halo;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_gal_1halo_first);
       if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-gal @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
     else if (index_integrand>=ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_first && index_integrand <= ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_last && ptsz->has_kSZ_kSZ_lensmag_1halo){
        Pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_lensmag_1halo;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^kSZ-kSZ-lensmag @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
    else if (index_integrand>=ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_first && index_integrand <= ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_last && ptsz->has_tSZ_tSZ_tSZ_1halo){
       Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_tSZ_tSZ_1halo;
       Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_first);
       if (ptsz->sz_verbose > 0) printf("computing b^y-y-y @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
     }
     else if (index_integrand>=ptsz->index_integrand_id_pk_at_z_1h_first && index_integrand <= ptsz->index_integrand_id_pk_at_z_1h_last && ptsz->has_pk_at_z_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_pk_at_z_1h;
        Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_at_z_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing pk^1h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_pk_at_z_2h_first && index_integrand <= ptsz->index_integrand_id_pk_at_z_2h_last && ptsz->has_pk_at_z_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_pk_at_z_2h;
         Pvectsz[ptsz->index_k_for_pk_hm] = (double) (index_integrand - ptsz->index_integrand_id_pk_at_z_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing pk^2h @ k_id = %.0f\n",Pvectsz[ptsz->index_k_for_pk_hm]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_gal_gal_1h_first && index_integrand <= ptsz->index_integrand_id_gal_gal_1h_last && ptsz->has_gal_gal_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_gal_gal_1h;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_gal_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^gal-gal_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
     else if (index_integrand>=ptsz->index_integrand_id_gal_gal_2h_first && index_integrand <= ptsz->index_integrand_id_gal_gal_2h_last && ptsz->has_gal_gal_2h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_gal_gal_2h;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_gal_2h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^gal-gal_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_gal_lens_1h_first && index_integrand <= ptsz->index_integrand_id_gal_lens_1h_last && ptsz->has_gal_lens_1h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_gal_lens_1h;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_lens_1h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^gal-lens_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_gal_lens_2h_first && index_integrand <= ptsz->index_integrand_id_gal_lens_2h_last && ptsz->has_gal_lens_2h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_gal_lens_2h;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_gal_lens_2h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^gal-lens_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_tSZ_cib_1h_first && index_integrand <= ptsz->index_integrand_id_tSZ_cib_1h_last && ptsz->has_tSZ_cib_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_cib_1h;
        int index_multipole_cib1 = (double) (index_integrand - ptsz->index_integrand_id_tSZ_cib_1h_first);
        int index_cib1 = index_multipole_cib1 / ptsz->nlSZ;
        int index_multipole = index_multipole_cib1 % ptsz->nlSZ;
        Pvectsz[ptsz->index_multipole] = (double) index_multipole ;
         Pvectsz[ptsz->index_frequency_for_cib_profile] = (double) index_cib1;
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
          if (ptsz->sz_verbose > 0) printf("computing cl^lens-cib_2h @ frequency_id = %.0f, ell_id = %.0f\n",
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
         if (ptsz->sz_verbose > 0) printf("computing cl^cib-cib_1h @ frequency_id = %.0f, frequency_prime_id = %.0f, ell_id = %.0f\n",
                                          Pvectsz[ptsz->index_frequency_for_cib_profile],
                                          Pvectsz[ptsz->index_frequency_prime_for_cib_profile],
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
         //  if (ptsz->sz_verbose > 0) printf("computing cl^cib-cib_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
         if (ptsz->sz_verbose > 0) printf("computing cl^cib-cib_2h @ frequency_id = %.0f, frequency_prime_id = %.0f, ell_id = %.0f\n",
                                Pvectsz[ptsz->index_frequency_for_cib_profile],
                                Pvectsz[ptsz->index_frequency_prime_for_cib_profile],
                                Pvectsz[ptsz->index_multipole]);
        }
     else if (index_integrand>=ptsz->index_integrand_id_tSZ_gal_1h_first && index_integrand <= ptsz->index_integrand_id_tSZ_gal_1h_last && ptsz->has_tSZ_gal_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_gal_1h;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_gal_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^y-gal_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_tSZ_gal_2h_first && index_integrand <= ptsz->index_integrand_id_tSZ_gal_2h_last && ptsz->has_tSZ_gal_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_gal_2h;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_gal_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^y-gal_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_lens_lens_1h_first && index_integrand <= ptsz->index_integrand_id_lens_lens_1h_last && ptsz->has_lens_lens_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_lens_lens_1h;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_lens_lens_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^lens-lens_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_lens_lens_2h_first && index_integrand <= ptsz->index_integrand_id_lens_lens_2h_last && ptsz->has_lens_lens_2h){
         Pvectsz[ptsz->index_md] = ptsz->index_md_lens_lens_2h;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_lens_lens_2h_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^lens-lens_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_tSZ_lens_1h_first && index_integrand <= ptsz->index_integrand_id_tSZ_lens_1h_last && ptsz->has_tSZ_lens_1h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_lens_1h;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_lens_1h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^y-phi_1h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
     else if (index_integrand>=ptsz->index_integrand_id_tSZ_lens_2h_first && index_integrand <= ptsz->index_integrand_id_tSZ_lens_2h_last && ptsz->has_tSZ_lens_2h){
        Pvectsz[ptsz->index_md] = ptsz->index_md_tSZ_lens_2h;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_tSZ_lens_2h_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^y-phi_2h @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
      else if (index_integrand>=ptsz->index_integrand_id_isw_lens_first && index_integrand <= ptsz->index_integrand_id_isw_lens_last && ptsz->has_isw_lens){
         Pvectsz[ptsz->index_md] = ptsz->index_md_isw_lens;
         Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_isw_lens_first);
         if (ptsz->sz_verbose > 0) printf("computing cl^isw-phi @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
       }
     else if (index_integrand>=ptsz->index_integrand_id_isw_tsz_first && index_integrand <= ptsz->index_integrand_id_isw_tsz_last && ptsz->has_isw_tsz){
        Pvectsz[ptsz->index_md] = ptsz->index_md_isw_tsz;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_isw_tsz_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^isw-y @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }
     else if (index_integrand>=ptsz->index_integrand_id_isw_auto_first && index_integrand <= ptsz->index_integrand_id_isw_auto_last && ptsz->has_isw_auto){
        Pvectsz[ptsz->index_md] = ptsz->index_md_isw_auto;
        Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_isw_auto_first);
        if (ptsz->sz_verbose > 0) printf("computing cl^isw-isw @ ell_id = %.0f\n",Pvectsz[ptsz->index_multipole]);
      }


     else
     {
       //printf("id not found. index_integrand = %d \n",index_integrand);
     return 0;
}

 //return 0;

 int index_md = (int) Pvectsz[ptsz->index_md];


 if (_dndlnM_){
   int index_z = (int) Pvectsz[ptsz->index_redshift_for_dndlnM];
   int index_m = (int) Pvectsz[ptsz->index_mass_for_dndlnM];

   double z_asked = ptsz->dndlnM_array_z[index_z];
   double m_asked = ptsz->dndlnM_array_m[index_m];

  ptsz->dndlnM_at_z_and_M[index_z][index_m] = get_dndlnM_at_z_and_M(z_asked,m_asked,ptsz);

  //printf("index_z = %d index_m =%d\n",index_z,index_m);
  //printf("z = %.4e m = %.4e dndlnM = %.4e\n",ptsz->dndlnM_array_z[index_z],ptsz->dndlnM_array_m[index_m],ptsz->dndlnM_at_z_and_M[index_z][index_m]);


 }
 else {

   if (_kSZ_kSZ_gal_1halo_
    || _kSZ_kSZ_lensmag_1halo_){
     // loop over l1,l2 for each ell
     int index_ell_1;
     int index_ell_2;
     int index_ell_3 = (int) Pvectsz[ptsz->index_multipole];


    double **b_l1_l2_l;
    class_alloc(b_l1_l2_l,
                ptsz->N_kSZ2_gal_multipole_grid*sizeof(double *),
                ptsz->error_message);

                for (index_ell_1=0;
                     index_ell_1<ptsz->N_kSZ2_gal_multipole_grid;
                     index_ell_1++)
                {
                  class_alloc(b_l1_l2_l[index_ell_1],
                              ptsz->N_kSZ2_gal_multipole_grid*sizeof(double),
                              ptsz->error_message);
                }



     for (index_ell_1=0;index_ell_1<ptsz->N_kSZ2_gal_multipole_grid;index_ell_1++){
       for (index_ell_2=0;index_ell_2<=index_ell_1;index_ell_2++){

         if (ptsz->sz_verbose > 0)
          printf("computing b_kSZ_kSZ_X_1h @ l3_id = %d and (l1,l2) = (%d,%d)\n",
                  index_ell_3, index_ell_1, index_ell_2);
          Pvectsz[ptsz->index_multipole_1] = index_ell_1;
          Pvectsz[ptsz->index_multipole_2] = index_ell_2;
          Pvectsz[ptsz->index_multipole_3] = index_ell_3;

          class_call(integrate_over_redshift(pba,
                                             pnl,
                                             ppm,
                                             ptsz,
                                             Pvecback,
                                             Pvectsz),
                          ptsz->error_message,
                          ptsz->error_message);
          //Pvectsz[ptsz->index_integral] = 0.;

       b_l1_l2_l[index_ell_1][index_ell_2] = Pvectsz[ptsz->index_integral];
       if (index_ell_1 != index_ell_2)
       b_l1_l2_l[index_ell_2][index_ell_1] = Pvectsz[ptsz->index_integral];
       }
     }



   // put bispectrum in 1d format for 2d interpolation
   double * b_l1_l2_l_1d;
   double * ln_ell;
   class_alloc(b_l1_l2_l_1d,
               sizeof(double *)*ptsz->N_kSZ2_gal_multipole_grid*ptsz->N_kSZ2_gal_multipole_grid,
               ptsz->error_message);
   class_alloc(ln_ell,
               sizeof(double *)*ptsz->N_kSZ2_gal_multipole_grid,
               ptsz->error_message);

   int index_l1_l2 = 0;
   for (index_ell_1=0;index_ell_1<ptsz->N_kSZ2_gal_multipole_grid;index_ell_1++){
     ln_ell[index_ell_1] = log(ptsz->ell_kSZ2_gal_multipole_grid[index_ell_1]);
     for (index_ell_2=0;index_ell_2<ptsz->N_kSZ2_gal_multipole_grid;index_ell_2++){
     b_l1_l2_l_1d[index_l1_l2] = log(b_l1_l2_l[index_ell_1][index_ell_2]);
     index_l1_l2 += 1;
     }
   }


   free(b_l1_l2_l);

   // now we integrate the bispectrum to compute power spectrum
   double cl_kSZ2_gal = 0.;
   //int ell_prime;
   int abs_ell_minus_ell_prime;
   double ell_min = ptsz->ell_kSZ2_gal_multipole_grid[0];
   double ell_max = ptsz->l_unwise_filter[ptsz->unwise_filter_size-1];//ptsz->ell_kSZ2_gal_multipole_grid[ptsz->N_kSZ2_gal_multipole_grid-1];
   //double ell_max = ptsz->ell_kSZ2_gal_multipole_grid[ptsz->N_kSZ2_gal_multipole_grid-1];

   //integrate over theta
   double theta = 0.;
   double dtheta = 0.01;

   double dell_prime =  1.;



  // adaptative integration
   struct Parameters_for_integrand_kSZ2_X_at_theta V;

   V.ptsz = ptsz;
   V.ln_ell = ln_ell;
   V.index_ell_3 = index_ell_3;
   V.b_l1_l2_l_1d = b_l1_l2_l_1d;
   void * params;

   double r; //result of the integral
   double epsrel= 1.e-4;//ptsz->redshift_epsrel;//ptsz->patterson_epsrel;
   double epsabs= 1.e-35;//ptsz->redshift_epsabs;//ptsz->patterson_epsabs;
   int show_neval = 0;//ptsz->patterson_show_neval;

   while(theta < 2.*_PI_){

    V.theta = theta;
    params = &V;
    //printf("ell_min = %.3e \t ell_max = %.3e\n", ell_min,ell_max);
     r=Integrate_using_Patterson_adaptive(log(ell_min), log(ell_max),
                                          epsrel, epsabs,
                                          integrand_kSZ2_X_at_theta,
                                          params,show_neval);
   //printf("integral result = %.3e\n",r);
    //
    // double lnell_prime = log(ell_min);
    // //double lnell_prime = ell_min;
    // double dlnell_prime =  .1;
    // double r=0.;
    //  while (lnell_prime <= log(ell_max)){
    //     //while (lnell_prime <= ell_max){
    //   r += dlnell_prime*integrand_kSZ2_X_at_theta(lnell_prime,&V);
    //   lnell_prime += dlnell_prime;
    //  }


     cl_kSZ2_gal += dtheta*r;
     theta += dtheta;
   }





   free(b_l1_l2_l_1d);
   free(ln_ell);






   int index_l = index_ell_3;

   if (_kSZ_kSZ_gal_1halo_)
    ptsz->cl_kSZ_kSZ_gal_1h[index_l] = cl_kSZ2_gal;
   else if(_kSZ_kSZ_lensmag_1halo_)
    ptsz->cl_kSZ_kSZ_lensmag_1h[index_l] = cl_kSZ2_gal;


   }

   else {
   class_call(integrate_over_redshift(pba,
                                      pnl,
                                      ppm,
                                      ptsz,
                                      Pvecback,
                                      Pvectsz),
                   ptsz->error_message,
                   ptsz->error_message);
          }

 }

   if (_hmf_){

      ptsz->hmf_int = Pvectsz[ptsz->index_integral];

   }

   if (_mean_y_){
      ptsz->y_monopole = Pvectsz[ptsz->index_integral]/pow(ptsz->Tcmb_gNU,1)/1.e6; //1e6 to convert Tcmb in micro Kelvins

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

 //  if (_kSZ_kSZ_gal_1halo_){
 //    int index_l = (int) Pvectsz[ptsz->index_multipole];
 //    //int index_m = (int) Pvectsz[ptsz->index_mass_bin_1];
 //   ptsz->cl_kSZ_kSZ_gal_1h[index_l] = Pvectsz[ptsz->index_integral];///pow(ptsz->Tcmb_gNU,1)/1.e6;
 //
 //
 // }

  if (_tSZ_tSZ_tSZ_1halo_){
    int index_l = (int) Pvectsz[ptsz->index_multipole];
    //int index_m = (int) Pvectsz[ptsz->index_mass_bin_1];
   ptsz->b_tSZ_tSZ_tSZ_1halo[index_l] = Pvectsz[ptsz->index_integral]/pow(ptsz->Tcmb_gNU,3)/1.e18; // dimensionless


 }

 if (_pk_at_z_1h_){
  int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
  ptsz->pk_at_z_1h[index_k] = Pvectsz[ptsz->index_integral];

}

if (_pk_at_z_2h_){
 int index_k = (int) Pvectsz[ptsz->index_k_for_pk_hm];
 ptsz->pk_at_z_2h[index_k] = Pvectsz[ptsz->index_integral];

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

return _SUCCESS_;
}





double integrand_at_m_and_z(double logM,
                             double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct primordial * ppm,
                             struct nonlinear * pnl,
                             struct tszspectrum * ptsz)
    {


   evaluate_HMF(logM,pvecback,pvectsz,pba,pnl,ptsz);
   // double z_asked = pvectsz[ptsz->index_z];
   // double m_asked = exp(logM);
   // // double calc_dn = get_dndlnM_at_z_and_M(z_asked,m_asked,ptsz);
   // // printf("z = %.2e M = %.4e \t eHMF = %.4e \t get_dn = %.4e\n",pvectsz[ptsz->index_z],m_asked,pvectsz[ptsz->index_hmf],
   // // calc_dn );
   //
   // pvectsz[ptsz->index_hmf] = get_dndlnM_at_z_and_M(z_asked,m_asked,ptsz);

   evaluate_completeness(pvecback,pvectsz,pba,ptsz);
   int index_l = (int) pvectsz[ptsz->index_multipole];
   pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];
   evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];




      int index_md = (int) pvectsz[ptsz->index_md];

   if (_hmf_){

      pvectsz[ptsz->index_integrand] = //pvectsz[ptsz->index_chi2]
                                       pvectsz[ptsz->index_hmf]
                                       *pvectsz[ptsz->index_completeness];
   }




   else if (_mean_y_){
     pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                       pvectsz[ptsz->index_hmf]
                                       *pvectsz[ptsz->index_completeness]
                                       *pow(pressure_profile_at_ell,1.);

   }


   else if (_tSZ_power_spectrum_){

      //int index_l = (int) pvectsz[ptsz->index_multipole];

         pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                           pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_dlnMdeltadlnM]
                                           *pvectsz[ptsz->index_completeness]
                                           *pow(pressure_profile_at_ell,2.);

    }


   else if (_trispectrum_){



     int index_l_prime = (int) pvectsz[ptsz->index_multipole_prime];
     pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l_prime];
     evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);
     double pressure_profile_at_ell_prime = pvectsz[ptsz->index_pressure_profile];

     pvectsz[ptsz->index_integrand] = //pvectsz[ptsz->index_chi2]
                                     pvectsz[ptsz->index_hmf]
                                     *pvectsz[ptsz->index_completeness]
                                     *pow(pressure_profile_at_ell,2.)
                                     *pow(pressure_profile_at_ell_prime,2.);
 }

 else if (_2halo_){

     //int index_l = (int) pvectsz[ptsz->index_multipole];
     //double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];
     evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);


     //this integrand is squared afterward
     pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                        *pvectsz[ptsz->index_dlnMdeltadlnM]
                                        *pvectsz[ptsz->index_halo_bias]
                                        *pvectsz[ptsz->index_completeness]
                                        *pow(pressure_profile_at_ell,1.);
                                        //*pow(pvectsz[ptsz->index_chi2]
                                        //*pvectsz[ptsz->index_pk_for_halo_bias],0.5);

 }

 else if  (_te_y_y_){

          //int index_l = (int) pvectsz[ptsz->index_multipole];
           evaluate_temperature_mass_relation(pvecback,pvectsz,pba,ptsz);

          pvectsz[ptsz->index_integrand] = pvectsz[ptsz->index_te_of_m]
                                           //*pvectsz[ptsz->index_chi2]
                                           *pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_completeness]
                                           *pow(pvectsz[ptsz->index_pressure_profile],2.);


        }
 else if  (_m_y_y_1h_){

          //int index_l = (int) pvectsz[ptsz->index_multipole];
           //evaluate_temperature_mass_relation(pvecback,pvectsz,pba,ptsz);

           pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                             pvectsz[ptsz->index_mass_for_hmf]
                                             *pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_dlnMdeltadlnM]
                                             *pvectsz[ptsz->index_completeness]
                                             *pow(pressure_profile_at_ell,2.);
        }
 else if  (_m_y_y_2h_){

           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);
           //evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);

           //this integrand is squared afterward
           pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                              *pvectsz[ptsz->index_mass_for_hmf]
                                              *pvectsz[ptsz->index_dlnMdeltadlnM]
                                              *pvectsz[ptsz->index_halo_bias]
                                              *pvectsz[ptsz->index_completeness]
                                              *pow(pressure_profile_at_ell,1.);
                                              //*pow(pvectsz[ptsz->index_pk_for_halo_bias],0.5);

        }
 else if  (_cov_Y_N_){


           pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                             pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness]
                                             *pow(pvectsz[ptsz->index_pressure_profile],2.);


        }

 else if  (_cov_N_N_){

           pvectsz[ptsz->index_integrand] =  ptsz->Omega_survey
                                             //*pvectsz[ptsz->index_chi2]
                                             *pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness];
                      }

 else if  (_cov_N_N_hsv_){


           if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

           evaluate_sigma2_hsv(pvecback,pvectsz,pba,pnl,ptsz);
           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);
           //evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);

           //printf("%.4e \t %.4e \t %4.e \n",pvectsz[ptsz->index_sigma2_hsv],pvectsz[ptsz->index_halo_bias],pvectsz[ptsz->index_completeness]);

           pvectsz[ptsz->index_integrand] =  ptsz->Omega_survey
                                             //*pvectsz[ptsz->index_chi2]
                                             *pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness]
                                             *pvectsz[ptsz->index_halo_bias]
                                             *pvectsz[ptsz->index_sigma2_hsv];
                                           }

           if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

           pvectsz[ptsz->index_integrand] =  ptsz->Omega_survey
                                             //*pvectsz[ptsz->index_chi2]
                                             *pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness]
                                             *pvectsz[ptsz->index_halo_bias];

                                           }

        }


 else if  (_cov_Y_N_next_order_){

           if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

           evaluate_sigma2_hsv(pvecback,pvectsz,pba,pnl,ptsz);
           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

           pvectsz[ptsz->index_integrand] =  ptsz->Omega_survey
                                             //*pvectsz[ptsz->index_chi2]
                                             *pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness]
                                             *pvectsz[ptsz->index_halo_bias]
                                             *pvectsz[ptsz->index_sigma2_hsv];
                                           }

           if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

           pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                             pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness]
                                             *pvectsz[ptsz->index_halo_bias]
                                             *pow(pressure_profile_at_ell,2.);

                                           }
        }
 else if  (_cov_Y_Y_ssc_){

           if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

           evaluate_sigma2_hsv(pvecback,pvectsz,pba,pnl,ptsz);
           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

           pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                             pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness]
                                             *pvectsz[ptsz->index_halo_bias]
                                             *pow(pressure_profile_at_ell,2.)
                                             *pvectsz[ptsz->index_sigma2_hsv];
                                           }

           if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

           evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

           pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                             pvectsz[ptsz->index_hmf]
                                             *pvectsz[ptsz->index_completeness]
                                             *pvectsz[ptsz->index_halo_bias]
                                             *pow(pressure_profile_at_ell,2.);

                                           }
        }
 else if (_tSZ_tSZ_tSZ_1halo_){

    //int index_l = (int) pvectsz[ptsz->index_multipole];
   pvectsz[ptsz->index_multipole_for_pressure_profile] = 10.;//the actual multipole
   evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_1 = pvectsz[ptsz->index_pressure_profile];
   int index_l = (int) pvectsz[ptsz->index_multipole];
   pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l];
   evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_2 = pvectsz[ptsz->index_pressure_profile];
   pvectsz[ptsz->index_multipole_for_pressure_profile] = ptsz->ell[index_l]+10.;
   evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);
   double pressure_profile_at_ell_3 = pvectsz[ptsz->index_pressure_profile];

       pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                         pvectsz[ptsz->index_hmf]
                                         *pvectsz[ptsz->index_dlnMdeltadlnM]
                                         //*pvectsz[ptsz->index_completeness]
                                         *pressure_profile_at_ell_1
                                         *pressure_profile_at_ell_2
                                         *pressure_profile_at_ell_3;
     // printf("l1 = %.3e\n",pressure_profile_at_ell_1);
     // printf("l2 = %.3e\n",pressure_profile_at_ell_2);
     // printf("l3 = %.3e\n",pressure_profile_at_ell_3);


  }

 else if (_kSZ_kSZ_gal_1halo_){

   int index_l_1 = (int) pvectsz[ptsz->index_multipole_1];
   pvectsz[ptsz->index_multipole_for_tau_profile] = ptsz->ell_kSZ2_gal_multipole_grid[index_l_1];//the actual multipole
   evaluate_tau_profile(pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];

   int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
   pvectsz[ptsz->index_multipole_for_tau_profile] = ptsz->ell_kSZ2_gal_multipole_grid[index_l_2];
   evaluate_tau_profile(pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];

   int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
   pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l_3];
   evaluate_galaxy_profile(pvecback,pvectsz,pba,ptsz);
   double galaxy_profile_at_ell_3 = pvectsz[ptsz->index_galaxy_profile];

   // //velocity dispersion (kSZ)
   // evaluate_vrms2(pvecback,pvectsz,pba,pnl,ptsz);

       pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                         *tau_profile_at_ell_1
                                         *tau_profile_at_ell_2
                                         *galaxy_profile_at_ell_3;


  }


 else if (_kSZ_kSZ_lensmag_1halo_){

   int index_l_1 = (int) pvectsz[ptsz->index_multipole_1];
   pvectsz[ptsz->index_multipole_for_tau_profile] = ptsz->ell_kSZ2_gal_multipole_grid[index_l_1];//the actual multipole
   evaluate_tau_profile(pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_1 = pvectsz[ptsz->index_tau_profile];

   int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
   pvectsz[ptsz->index_multipole_for_tau_profile] = ptsz->ell_kSZ2_gal_multipole_grid[index_l_2];
   evaluate_tau_profile(pvecback,pvectsz,pba,ptsz);
   double tau_profile_at_ell_2 = pvectsz[ptsz->index_tau_profile];

   int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
   pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l_3];
   evaluate_lensing_profile(pvecback,pvectsz,pba,ptsz);
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



  else if (_pk_at_z_1h_){

    evaluate_density_profile(pvecback,pvectsz,pba,ptsz);
    double density_profile_at_k_1 = pvectsz[ptsz->index_density_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *pvectsz[ptsz->index_mass_for_hmf]
                                      *pvectsz[ptsz->index_mass_for_hmf]
                                      //
                                      //*pvectsz[ptsz->index_mVIR]
                                      //*pvectsz[ptsz->index_mVIR]
                                      *density_profile_at_k_1
                                      *density_profile_at_k_1;
   }



  else if (_gal_gal_1h_){

    int index_l = (int) pvectsz[ptsz->index_multipole];
    pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
    evaluate_galaxy_profile(pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_ell_1 = pvectsz[ptsz->index_galaxy_profile];

        pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                          *galaxy_profile_at_ell_1
                                          *galaxy_profile_at_ell_1;
   }

   else if (_gal_gal_2h_){

     int index_l = (int) pvectsz[ptsz->index_multipole];
     pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
     evaluate_galaxy_profile(pvecback,pvectsz,pba,ptsz);
     double galaxy_profile_at_ell_1 = pvectsz[ptsz->index_galaxy_profile];

     evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias]
                                           *galaxy_profile_at_ell_1;

    }

   else if (_gal_lens_1h_){
              // Here we follow only KFW20


              // // WIxSC (KA20) use Tinker et al 2010 linear bias (halo bias):
              // if (ptsz->galaxy_sample==0)
              // evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);
              //  // unWISE use an effective linear bias Eq. 6.2 of KFW20 (1909.07412).
              //  else  if (ptsz->galaxy_sample==1)
              //  evaluate_effective_galaxy_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

              int index_l = (int) pvectsz[ptsz->index_multipole];
              pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
              evaluate_lensing_profile(pvecback,pvectsz,pba,ptsz);
              //double lensing_profile_at_ell_1 =

             //int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_galaxy_profile] =  ptsz->ell[index_l];
             evaluate_galaxy_profile(pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               //*pvectsz[ptsz->index_dlnMdeltadlnM]
                                               //*pvectsz[ptsz->index_completeness]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_galaxy_profile];
                                               //*pvectsz[ptsz->index_halo_bias];

          }
   else if (_gal_lens_2h_){

              // WIxSC (KA20) use Tinker et al 2010 linear bias (halo bias):
              if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

              evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_galaxy_profile] =  ptsz->ell[index_l];
             evaluate_galaxy_profile(pvecback,pvectsz,pba,ptsz);



             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_galaxy_profile]
                                               *pvectsz[ptsz->index_halo_bias];


            }

          if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

          evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

         int index_l = (int) pvectsz[ptsz->index_multipole];


         pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
         evaluate_lensing_profile(pvecback,pvectsz,pba,ptsz);

         pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                           //*pvectsz[ptsz->index_dlnMdeltadlnM]
                                           //*pvectsz[ptsz->index_completeness]
                                           *pvectsz[ptsz->index_lensing_profile]
                                          // *pvectsz[ptsz->index_galaxy_profile]
                                           *pvectsz[ptsz->index_halo_bias];
      }
        }


  else if (_cib_cib_1h_){

    int index_l = (int) pvectsz[ptsz->index_multipole];
    pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];
    evaluate_cib_profile(pvecback,pvectsz,pba,ptsz);
    double cib_profile_at_ell_1 = pvectsz[ptsz->index_cib_profile];

        pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                          *cib_profile_at_ell_1
                                          *cib_profile_at_ell_1;
   }

   else if (_cib_cib_2h_){
             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];
             //pvectsz[ptsz->index_frequency_for_cib_profile] = ptsz->nu_cib_GHz;
             evaluate_cib_profile(pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_cib_profile]
                                               *pvectsz[ptsz->index_halo_bias];

    }

  else if (_tSZ_cib_1h_){

    int index_l = (int) pvectsz[ptsz->index_multipole];
    pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];
    evaluate_cib_profile(pvecback,pvectsz,pba,ptsz);
    double cib_profile_at_ell_1 = pvectsz[ptsz->index_cib_profile];
    pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
    evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);
    double pressure_profile_at_ell_2 = pvectsz[ptsz->index_pressure_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *cib_profile_at_ell_1
                                      *pressure_profile_at_ell_2;
   }
   else if  (_tSZ_cib_2h_){


             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
             evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_pressure_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];
             evaluate_cib_profile(pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_cib_profile]
                                               *pvectsz[ptsz->index_halo_bias];

                                             }

          }

  else if (_lens_cib_1h_){

    int index_l = (int) pvectsz[ptsz->index_multipole];
    pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];
    evaluate_cib_profile(pvecback,pvectsz,pba,ptsz);
    double cib_profile_at_ell_1 = pvectsz[ptsz->index_cib_profile];
    pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
    evaluate_lensing_profile(pvecback,pvectsz,pba,ptsz);
    double lensing_profile_at_ell_2 = pvectsz[ptsz->index_lensing_profile];

    pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                      *cib_profile_at_ell_1
                                      *lensing_profile_at_ell_2;
   }
   else if  (_lens_cib_2h_){


             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_lensing_profile] =  ptsz->ell[index_l];
             evaluate_lensing_profile(pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_cib_profile] = ptsz->ell[index_l];
             evaluate_cib_profile(pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_cib_profile]
                                               *pvectsz[ptsz->index_halo_bias];

                                             }

          }


  else if (_tSZ_gal_1h_){

    int index_l = (int) pvectsz[ptsz->index_multipole];
    pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
    evaluate_galaxy_profile(pvecback,pvectsz,pba,ptsz);
    double galaxy_profile_at_ell_1 = pvectsz[ptsz->index_galaxy_profile];
    pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
    evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);
    double pressure_profile_at_ell_2 = pvectsz[ptsz->index_pressure_profile];

        pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                          *galaxy_profile_at_ell_1
                                          *pressure_profile_at_ell_2;
   }

   else if  (_tSZ_gal_2h_){


             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {
             // WIxSC (KA20) use Tinker et al 2010 linear bias (halo bias):
             //if (ptsz->galaxy_sample==0)
             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);
             //  // unWISE use an effective linear bias Eq. 6.2 of KFW20 (1909.07412).
             // else  if (ptsz->galaxy_sample==1)
             // evaluate_effective_galaxy_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

             //evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);
             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
             evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_dlnMdeltadlnM]
                                               *pvectsz[ptsz->index_completeness]
                                               *pvectsz[ptsz->index_pressure_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

             // WIxSC (KA20) use Tinker et al 2010 linear bias (halo bias):
             //if (ptsz->galaxy_sample==0)
             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);
              // unWISE use an effective linear bias Eq. 6.2 of KFW20 (1909.07412).
             // else  if (ptsz->galaxy_sample==1)
             // evaluate_effective_galaxy_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

             //evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);
             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_galaxy_profile] = ptsz->ell[index_l];
             evaluate_galaxy_profile(pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_dlnMdeltadlnM]
                                               *pvectsz[ptsz->index_completeness]
                                               *pvectsz[ptsz->index_galaxy_profile]
                                               *pvectsz[ptsz->index_halo_bias];

                                             }

          }
  else if (_lens_lens_1h_){

    int index_l = (int) pvectsz[ptsz->index_multipole];
    pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
    evaluate_lensing_profile(pvecback,pvectsz,pba,ptsz);
    double lensing_profile_at_ell_1 = pvectsz[ptsz->index_lensing_profile];

        pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                          pvectsz[ptsz->index_hmf]
                                          *pvectsz[ptsz->index_dlnMdeltadlnM]
                                          *pvectsz[ptsz->index_completeness]
                                          *lensing_profile_at_ell_1
                                          *lensing_profile_at_ell_1;
   }

   else if (_lens_lens_2h_){

     int index_l = (int) pvectsz[ptsz->index_multipole];
     pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
     evaluate_lensing_profile(pvecback,pvectsz,pba,ptsz);
     double lensing_profile_at_ell_1 = pvectsz[ptsz->index_lensing_profile];

     evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_halo_bias]
                                           *lensing_profile_at_ell_1;
    // printf("bias = %.3e\n",pvectsz[ptsz->index_halo_bias]);
    // printf("hmf = %.3e\n",pvectsz[ptsz->index_hmf]);
    // printf("ug = %.3e\n",galaxy_profile_at_ell_1);
    }


  else if (_tSZ_lens_1h_){

    int index_l = (int) pvectsz[ptsz->index_multipole];
    pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
    evaluate_lensing_profile(pvecback,pvectsz,pba,ptsz);
    double lensing_profile_at_ell_1 = pvectsz[ptsz->index_lensing_profile];
    pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
    evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);
    double pressure_profile_at_ell_2 = pvectsz[ptsz->index_pressure_profile];

        pvectsz[ptsz->index_integrand] =  //pvectsz[ptsz->index_chi2]
                                          pvectsz[ptsz->index_hmf]
                                          *pvectsz[ptsz->index_dlnMdeltadlnM]
                                          *pvectsz[ptsz->index_completeness]
                                          *lensing_profile_at_ell_1
                                          *pressure_profile_at_ell_2;
   }

   else if  (_tSZ_lens_2h_){


             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);
             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_pressure_profile] =  ptsz->ell[index_l];
             evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_dlnMdeltadlnM]
                                               *pvectsz[ptsz->index_completeness]
                                               *pvectsz[ptsz->index_pressure_profile]
                                               *pvectsz[ptsz->index_halo_bias];
                                             }

             if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2) {

             evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);
             int index_l = (int) pvectsz[ptsz->index_multipole];
             pvectsz[ptsz->index_multipole_for_lensing_profile] = ptsz->ell[index_l];
             evaluate_lensing_profile(pvecback,pvectsz,pba,ptsz);

             pvectsz[ptsz->index_integrand] =   pvectsz[ptsz->index_hmf]
                                               *pvectsz[ptsz->index_dlnMdeltadlnM]
                                               *pvectsz[ptsz->index_completeness]
                                               *pvectsz[ptsz->index_lensing_profile]
                                               *pvectsz[ptsz->index_halo_bias];

                                             }

          }



   else if (_isw_tsz_){

     evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);
     //evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
     pvectsz[ptsz->index_integrand] =   //pvectsz[ptsz->index_chi2]
                                       pvectsz[ptsz->index_hmf]
                                       *pvecback[pba->index_bg_D]
                                       *pvectsz[ptsz->index_dlnMdeltadlnM]
                                       *pvectsz[ptsz->index_completeness]
                                       *pvectsz[ptsz->index_halo_bias]
                                       *pressure_profile_at_ell;


   }

   return pvectsz[ptsz->index_integrand];
}


int evaluate_temperature_mass_relation( double * pvecback,
                                        double * pvectsz,
                                        struct background * pba,
                                        struct tszspectrum * ptsz)
  {


   double mass = pvectsz[ptsz->index_m500]; //biased mass = M/B (X-ray mass)
   //double mass = pvectsz[ptsz->index_m500]*ptsz->HSEbias; //biased mass = M/B (X-ray mass)

   double Eh = pvecback[pba->index_bg_H]/ptsz->H0_in_class_units;

if (ptsz->temperature_mass_relation == 0){
   pvectsz[ptsz->index_te_of_m] = 5.*pow(Eh*mass //biased mass
                                         /3.e14,2./3.); //kB*Te in keV

      }
   //lee et al [1912.07924v1]
   else if (ptsz->temperature_mass_relation == 1){


      double A[3] = {4.763,4.353,3.997};
      double B[3] = {0.581,0.571,0.593};
      double C[3] = {0.013,0.008,0.009};
      double z_interp[3] = {0.0,0.5,1.0};

      double Ap,Bp,Cp,zp;

      zp = pvectsz[ptsz->index_z];
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
   double M = mass*ptsz->HSEbias; //true mass

   pvectsz[ptsz->index_te_of_m] = Ap*pow(Eh,2./3.)*pow(M/Mfid,Bp+Cp*log(M/Mfid)); //kB*Te in keV

   }


   return _SUCCESS_;
}


int evaluate_tau_profile(double * pvecback,
                        double * pvectsz,
                        struct background * pba,
                        struct tszspectrum * ptsz)
{

   double result;
   double characteristic_radius;
   double characteristic_multipole;
   double tau_normalisation;


   characteristic_radius = pvectsz[ptsz->index_rs];///ptsz->cvir_tau_profile_factor;// in Mpc/h
   characteristic_multipole = pvectsz[ptsz->index_ls];//*ptsz->cvir_tau_profile_factor;
   //printf("l1tau = %f\n", pvectsz[ptsz->index_multipole_for_tau_profile]);
   pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile] = ptsz->cvir_tau_profile_factor*characteristic_multipole;
   //printf("l1char = %f\n", pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile]);

   pvectsz[ptsz->index_multipole_for_nfw_profile] = pvectsz[ptsz->index_multipole_for_tau_profile];

   class_call(two_dim_ft_nfw_profile(ptsz,pba,pvectsz,&result,1),
                                     ptsz->error_message,
                                     ptsz->error_message);


   // double lnx_asked = log((pvectsz[ptsz->index_multipole_for_nfw_profile]+0.5)/characteristic_multipole);
   //
   // if(lnx_asked<ptsz->RNFW_lnx[0])
   //    result = ptsz->RNFW_lnI[0];
   // else if (lnx_asked>ptsz->RNFW_lnx[ptsz->RNFW_lnx_size-1])
   //    result = -100.;
   // else  result =  pwl_value_1d(ptsz->RNFW_lnx_size,
   //                              ptsz->RNFW_lnx,
   //                              ptsz->RNFW_lnI,
   //                              lnx_asked);
   // result = exp(result);
   double rs = pvectsz[ptsz->index_rs];
   double rvir = pvectsz[ptsz->index_rVIR];
   double cvir_prime = ptsz->cvir_tau_profile_factor*pvectsz[ptsz->index_cVIR];
   double cvir = pvectsz[ptsz->index_cVIR];
   double xout = ptsz->x_out_nfw_profile*cvir;
   double xout_prime = ptsz->x_out_nfw_profile*cvir_prime;

  //printf("nfw_result=%.4e m(xout)=%.4e\n", result,m_nfw(xout_prime));
  pvectsz[ptsz->index_tau_profile] = result;

  // JCH
  // ls=da(z)/rs
  // rho0=mvir/(4d0*pi*rs**3d0*(dlog(1d0+cvir)-cvir/(1d0+cvir))) !Eq. (2.66) of Binney & Tremaine
  // tau2d=(obh2/h0**2d0)/om0/mu_e*f_free/h0 * & !h0 here converts Msun/h to Msun
  //      8.30702e-17 * h0**2d0 & !this is sigmaT / m_prot in (Mpc/h)**2/Msun
  //      * 4d0*pi*rho0*rs/(ls**2d0)*&
  //      rombint2(rhoNFW,xin,xoutrho*cvir,tol,(ell+0.5d0)/ls,cvir)


  tau_normalisation = pba->Omega0_b/ptsz->Omega_m_0/ptsz->mu_e*ptsz->f_free/pba->h;

  double sigmaT_over_mpc2 = 8.30702e-17 * pow(pba->h,2)/pba->h; // !this is sigmaT / m_prot in (Mpc/h)**2/(Msun/h)

  // adjust normalisation so that total mass is conserved
  // within same limiting radius r_max = xout*r_vir

  double factor_norm = 1;

  // rho0=mvir/(4d0*pi*rs**3d0*(dlog(1d0+cvir)-cvir/(1d0+cvir))) !Eq. (2.66) of Binney & Tremaine
  // double rho0 = pvectsz[ptsz->index_mVIR]/(4.*_PI_*pow(pvectsz[ptsz->index_rs],3.)
  //               *(log(1.+ptsz->cvir_tau_profile_factor*pvectsz[ptsz->index_cVIR])-ptsz->cvir_tau_profile_factor*pvectsz[ptsz->index_cVIR]
  //               /(1.+ptsz->cvir_tau_profile_factor*pvectsz[ptsz->index_cVIR])));

  double rho0 = pvectsz[ptsz->index_mVIR]/m_nfw(cvir_prime);
   //rho0 = 1.;

   rho0 *= m_nfw(xout)/m_nfw(xout_prime);
   rho0 *= m_nfw(cvir_prime)/m_nfw(cvir);

   pvectsz[ptsz->index_tau_profile] =  sigmaT_over_mpc2
                                       *tau_normalisation
                                       *rho0 // in (Msun/h)/(Mpc/h)**3
                                       *pvectsz[ptsz->index_tau_profile]
                                       *pow(pvecback[pba->index_bg_ang_distance]*pba->h,-2.); //(rs*ls)^2 in [Mpc/h]^2


   return _SUCCESS_;
}

int evaluate_density_profile(double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct tszspectrum * ptsz)
{

  int index_k = (int) pvectsz[ptsz->index_k_for_pk_hm];
  double k = ptsz->k_for_pk_hm[index_k];

   //int index_md = (int) pvectsz[ptsz->index_md];

   //int index_l = (int) pvectsz[ptsz->index_multipole_for_tau_profile];
   //printf("ell pp=%e\n",ptsz->ell[index_l]);

   double result;
   double characteristic_radius;
   //double characteristic_multipole;
   double density_normalisation;


   characteristic_radius = pvectsz[ptsz->index_rs]; // in Mpc/h

    class_call(two_dim_ft_nfw_profile(ptsz,pba,pvectsz,&result,0),
                                      ptsz->error_message,
                                      ptsz->error_message);



   pvectsz[ptsz->index_density_profile] = result;



  density_normalisation = 1.;

  // rho0=mvir/(4d0*pi*rs**3d0*(dlog(1d0+cvir)-cvir/(1d0+cvir))) !Eq. (2.66) of Binney & Tremaine
  double rho0 = pvectsz[ptsz->index_mVIR]/(4.*_PI_*pow(pvectsz[ptsz->index_rs],3.)*(log(1.+pvectsz[ptsz->index_cVIR])-pvectsz[ptsz->index_cVIR]/(1.+pvectsz[ptsz->index_cVIR])));

  //
  // printf("lens rho0 = %.3e\n",rho0);
  // printf("lens norm = %.3e\n",lensing_normalisation);
  // printf("Sigma_crit = %.3e\n",pvectsz[ptsz->index_lensing_Sigma_crit]);


   pvectsz[ptsz->index_density_profile] =  density_normalisation
                                           *rho0
                                           *pvectsz[ptsz->index_density_profile]
                                           /pvectsz[ptsz->index_mVIR]
                                           *(4*_PI_)
                                           *pow(characteristic_radius,3.); //rs in Mpc/h

  //printf("lens prof = %.3e\n",pvectsz[ptsz->index_lensing_profile]);

   return _SUCCESS_;
}



int evaluate_lensing_profile(double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct tszspectrum * ptsz)
{

   //int index_md = (int) pvectsz[ptsz->index_md];

   //int index_l = (int) pvectsz[ptsz->index_multipole_for_tau_profile];
   //printf("ell pp=%e\n",ptsz->ell[index_l]);

   double result;
   double characteristic_radius;
   double characteristic_multipole;
   double lensing_normalisation;


   characteristic_radius = pvectsz[ptsz->index_rs]; // in Mpc/h
   characteristic_multipole = pvectsz[ptsz->index_ls];
   pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile] = characteristic_multipole;
   pvectsz[ptsz->index_multipole_for_nfw_profile] = pvectsz[ptsz->index_multipole_for_lensing_profile];

    class_call(two_dim_ft_nfw_profile(ptsz,pba,pvectsz,&result,0),
                                      ptsz->error_message,
                                      ptsz->error_message);


      // double lnx_asked = log((pvectsz[ptsz->index_multipole_for_nfw_profile]+0.5)/characteristic_multipole);
      //
      // if(lnx_asked<ptsz->RNFW_lnx[0])
      //    result = ptsz->RNFW_lnI[0];
      // else if (lnx_asked>ptsz->RNFW_lnx[ptsz->RNFW_lnx_size-1])
      //    result = -100.;
      // else  result =  pwl_value_1d(ptsz->RNFW_lnx_size,
      //                              ptsz->RNFW_lnx,
      //                              ptsz->RNFW_lnI,
      //                              lnx_asked);
      // result = exp(result);

   // this is the transform of the nfw, dimensionless
   // x**(-1d0)*(1d0+x)**(-2d0)*x**2d0*dsin(y*x)/(y*x)
   pvectsz[ptsz->index_lensing_profile] = result;


  // lensing normalisation is dimensionless
  // we write down things in terms of phi here
  // (\kappa = l(l+1)\phi_l/2, see Hill & Spergel 1312.4525)

  int index_md = (int) pvectsz[ptsz->index_md];
  if ( _gal_lens_2h_
    || _gal_lens_1h_
    //|| _cib_lens_1h_
    //|| _cib_lens_1h_
    || _kSZ_kSZ_lensmag_1halo_
    || _lens_lens_1h_
    || _lens_lens_2h_)
  lensing_normalisation = 1.; // kappa
  else // phi
  lensing_normalisation = 2./(pvectsz[ptsz->index_multipole_for_lensing_profile]*(pvectsz[ptsz->index_multipole_for_lensing_profile]+1.));


  // rho0=mvir/(4d0*pi*rs**3d0*(dlog(1d0+cvir)-cvir/(1d0+cvir))) !Eq. (2.66) of Binney & Tremaine
  double rho0 = pvectsz[ptsz->index_mVIR]/(4.*_PI_*pow(pvectsz[ptsz->index_rs],3.)*(log(1.+pvectsz[ptsz->index_cVIR])-pvectsz[ptsz->index_cVIR]/(1.+pvectsz[ptsz->index_cVIR])));

  //
  // printf("lens rho0 = %.3e\n",rho0);
  // printf("lens norm = %.3e\n",lensing_normalisation);
  // printf("Sigma_crit = %.3e\n",pvectsz[ptsz->index_lensing_Sigma_crit]);


   pvectsz[ptsz->index_lensing_profile] =  lensing_normalisation
                                           *rho0
                                           *pvectsz[ptsz->index_lensing_profile]
                                           /pvectsz[ptsz->index_lensing_Sigma_crit]
                                           *(4*_PI_)
                                           *pow(characteristic_multipole,-2)
                                           *characteristic_radius; //rs in Mpc/h

  //printf("lens prof = %.3e\n",pvectsz[ptsz->index_lensing_profile]);

   return _SUCCESS_;
}

int evaluate_pressure_profile(double * pvecback,
                              double * pvectsz,
                              struct background * pba,
                              struct tszspectrum * ptsz)
{
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

   //int index_l = (int) pvectsz[ptsz->index_multipole_for_pressure_profile];
   //printf("ell pp=%e\n",ptsz->ell[index_l]);

   //custom gNFW pressure profile or Battaglia et al 2012
   if (ptsz->pressure_profile == 3 || ptsz->pressure_profile == 4 ){

      class_call(two_dim_ft_pressure_profile(ptsz,pba,pvectsz,&result),
                                              ptsz->error_message,
                                              ptsz->error_message);

   }

  else {
      //printf("ell pp=%e\n",result);

      lnx_asked = log((pvectsz[ptsz->index_multipole_for_pressure_profile]+0.5)/pvectsz[ptsz->index_l500]);

      if(lnx_asked<ptsz->PP_lnx[0] || _mean_y_)
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
      double f_b = pba->Omega0_b/ptsz->Omega_m_0;

      double Eh = pvecback[pba->index_bg_H]/ptsz->H0_in_class_units;

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

      // //Custom. GNFW
      // else if (ptsz->pressure_profile == 3)
      // pressure_normalisation = C_pressure
      //                          *ptsz->P0GNFW
      //                          *pow(0.7/pba->h, 1.5); // assuming X-ray data based pressure profile
     //Custom. GNFW
     else if (ptsz->pressure_profile == 3)
     pressure_normalisation = C_pressure
                              *ptsz->P0GNFW
                              *pow(0.7/pba->h, 1.); // assuming SZ data based pressure profile


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

int evaluate_completeness(double * pvecback,
                          double * pvectsz,
                          struct background * pba,
                          struct tszspectrum * ptsz){

    double comp_at_M_and_z = 0.;


    if (ptsz->has_completeness_for_ps_SZ == 1){
    comp_at_M_and_z = 0.;

    double mp_bias = pvectsz[ptsz->index_m500]; //biased mass = M/B
    //double redshift = pvectsz[ptsz->index_z];

    //printf("mass m500 = %e\n",mp_bias);
    //printf("mass m200 = %e\n",pvectsz[ptsz->index_m200]); //true mass
    //printf("bias = %e\n",ptsz->HSEbias);
    //printf("redshift = %e\n",redshift);
    double Eh = pvecback[pba->index_bg_H]/ptsz->H0_in_class_units;
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

    }

   if (ptsz->which_ps_sz == 0)
      pvectsz[ptsz->index_completeness] = 1.;
   else if (ptsz->which_ps_sz == 1) //ps resolved
      pvectsz[ptsz->index_completeness] = comp_at_M_and_z;
   else if (ptsz->which_ps_sz == 2) //ps unresolved
      pvectsz[ptsz->index_completeness] = (1.-comp_at_M_and_z);

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

  // blue sample
  if (ptsz->unwise_galaxy_sample_id==3)
    b = 0.8+1.2*z;
  // green (two greens, depending on the min halo mass cut)
  else if (ptsz->unwise_galaxy_sample_id==2 || ptsz->unwise_galaxy_sample_id==1)
    b = gsl_max(1.6*z*z, 1.);
  // red
  else if (ptsz->unwise_galaxy_sample_id==0)
    b = gsl_max(2.*pow(z,1.5),1.);


   pvectsz[ptsz->index_halo_bias] = b;
   return _SUCCESS_;
}





int evaluate_halo_bias(double * pvecback,
                       double * pvectsz,
                       struct background * pba,
                       struct primordial * ppm,
                       struct nonlinear * pnl,
                       struct tszspectrum * ptsz)
{

   double nu = exp(pvectsz[ptsz->index_lognu]);
   double nuTink = sqrt(nu); //in T10 paper: nu=delta/sigma while here nu=(delta/sigma)^2

   double Delta;

   //Tinker et al 2008 @ M1600-mean
   if (ptsz->MF==6)
   Delta = 1600.;

   //Jenkins et al 2001
   else if (ptsz->MF==3)
   Delta = 180.;

   //Tinker et al 2008 @ m500
   //Boquet et al 2015
   else if (ptsz->MF==5 || ptsz->MF==7)
   Delta = 500./pvecback[pba->index_bg_Omega_m];

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

   pvectsz[ptsz->index_halo_bias] = 1.-ATT*(pow(nuTink,aTTT)/(pow(nuTink,aTTT)+pow(ptsz->delta_cSZ,aTTT)))
                                                      +BTT*pow(nuTink,bTTT)+CTT*pow(nuTink,cTTT);
  //
  //  int index_l = (int)  pvectsz[ptsz->index_multipole];
  //  double z = pvectsz[ptsz->index_z];
  //  //identical to sqrt(pvectsz[index_chi2])
  //  double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi
  //
  //  pvectsz[ptsz->index_k_value_for_halo_bias] = (ptsz->ell[index_l]+0.5)/d_A; //units h/Mpc
  //
  //
  //  double k = pvectsz[ptsz->index_k_value_for_halo_bias]; //in h/Mpc
  //
  //  double pk;
  //  double * pk_ic = NULL;
  //
  //
  //
  // //Input: wavenumber in 1/Mpc
  // //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
  //    class_call(nonlinear_pk_at_k_and_z(
  //                                      pba,
  //                                      ppm,
  //                                      pnl,
  //                                      pk_linear,
  //                                      k*pba->h,
  //                                      z,
  //                                      pnl->index_pk_m,
  //                                      &pk, // number *out_pk_l
  //                                      pk_ic // array out_pk_ic_l[index_ic_ic]
  //                                    ),
  //                                    pnl->error_message,
  //                                    pnl->error_message);
  //
  //
  //  //now compute P(k) in units of h^-3 Mpc^3
  //  pvectsz[ptsz->index_pk_for_halo_bias] = pk*pow(pba->h,3.); //in units Mpc^3/h^3
   return _SUCCESS_;
}



//Compute P(k) in units of h^-3 Mpc^3
int evaluate_pk_at_ell_plus_one_half_over_chi(double * pvecback,
                                             double * pvectsz,
                                             struct background * pba,
                                             struct primordial * ppm,
                                             struct nonlinear * pnl,
                                             struct tszspectrum * ptsz)
{

   // double nu = exp(pvectsz[ptsz->index_lognu]);
   // double nuTink = sqrt(nu); //in T10 paper: nu=delta/sigma while here nu=(delta/sigma)^2
   //
   // double Delta;
   //
   // //Tinker et al 2008 @ M1600-mean
   // if (ptsz->MF==6)
   // Delta = 1600.;
   //
   // //Jenkins et al 2001
   // else if (ptsz->MF==3)
   // Delta = 180.;
   //
   // //Tinker et al 2008 @ m500
   // //Boquet et al 2015
   // else if (ptsz->MF==5 || ptsz->MF==7)
   // Delta = 500./pvecback[pba->index_bg_Omega_m];
   //
   // else
   // Delta = 200.;
   //
   //
   // // Table 2 of Tinker et al 2010:
   // double y = log10(Delta);
   // double ATT = 1.+0.24*y*exp(-pow(4./y,4.));
   // double aTTT = 0.44*y-0.88;
   // double BTT = 0.183;
   // double bTTT = 1.5;
   // double CTT = 0.019+0.107*y+0.19*exp(-pow(4./y,4.));
   // double cTTT = 2.4;
   //
   // pvectsz[ptsz->index_halo_bias] = 1.-ATT*(pow(nuTink,aTTT)/(pow(nuTink,aTTT)+pow(ptsz->delta_cSZ,aTTT)))
   //                                                    +BTT*pow(nuTink,bTTT)+CTT*pow(nuTink,cTTT);

  double k;
  int index_md = (int) pvectsz[ptsz->index_md];
  double z = pvectsz[ptsz->index_z];

if (_pk_at_z_2h_){
  int index_k = (int) pvectsz[ptsz->index_k_for_pk_hm];
  k = ptsz->k_for_pk_hm[index_k];

}
else {
   int index_l = (int)  pvectsz[ptsz->index_multipole];
   //identical to sqrt(pvectsz[index_chi2])
   double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi

   pvectsz[ptsz->index_k_value_for_halo_bias] = (ptsz->ell[index_l]+0.5)/d_A; //units h/Mpc


   k = pvectsz[ptsz->index_k_value_for_halo_bias]; //in h/Mpc
 }

   double pk;
   double * pk_ic = NULL;

   //printf("evaluating pk at k=%.3e h/Mpc\n",k);



   //int index_md = (int) pvectsz[ptsz->index_md];
   if (((ptsz->has_gal_gal_2h == _TRUE_) && (index_md == ptsz->index_md_gal_gal_2h) && (ptsz->galaxy_sample==1) && (ptsz->use_hod ==0)) // unWISE
       //|| ((ptsz->has_gal_gal_2h == _TRUE_) && (index_md == ptsz->index_md_gal_gal_2h) && (ptsz->galaxy_sample==0)) // WIxSC
       || ((ptsz->has_lens_lens_2h == _TRUE_) && (index_md == ptsz->index_md_lens_lens_2h) && (ptsz->use_hod ==0))
       || ((ptsz->has_gal_lens_2h == _TRUE_) && (index_md == ptsz->index_md_gal_lens_2h) && (ptsz->galaxy_sample==1) && (ptsz->use_hod ==0))) // unWISE
   {

       //Input: wavenumber in 1/Mpc
       //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$

          // halofit approach: uses non-linear pk
          class_call(nonlinear_pk_at_k_and_z(
                                            pba,
                                            ppm,
                                            pnl,
                                            pk_nonlinear,
                                            //pk_linear,
                                            k*pba->h,
                                            z,
                                            pnl->index_pk_m,
                                            &pk, // number *out_pk_l
                                            pk_ic // array out_pk_ic_l[index_ic_ic]
                                          ),
                                          pnl->error_message,
                                          pnl->error_message);

   }
   else{

       //Input: wavenumber in 1/Mpc
       //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
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
                                          pnl->error_message,
                                          pnl->error_message);
}



   //now compute P(k) in units of h^-3 Mpc^3
   pvectsz[ptsz->index_pk_for_halo_bias] = pk*pow(pba->h,3.); //in units Mpc^3/h^3
   return _SUCCESS_;
}






double evaluate_pk_halofit_over_pk_linear_at_ell_plus_one_half_over_chi(double * pvecback,
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

   //printf("evaluating pk at k=%.3e h/Mpc\n",k);



   int index_md = (int) pvectsz[ptsz->index_md];


          class_call(nonlinear_pk_at_k_and_z(
                                            pba,
                                            ppm,
                                            pnl,
                                            //pk_nonlinear,
                                            pk_nonlinear,
                                            k*pba->h,
                                            z,
                                            pnl->index_pk_m,
                                            &pk, // number *out_pk_l
                                            pk_ic // array out_pk_ic_l[index_ic_ic]
                                          ),
                                          pnl->error_message,
                                          pnl->error_message);

   double pk_halofit = pk;



       //Input: wavenumber in 1/Mpc
       //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
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
                                          pnl->error_message,
                                          pnl->error_message);

    double pk_lin  = pk;




   double  result = pk_halofit/pk_lin;
   return result;
}





int evaluate_pk_at_ell_plus_one_half_over_chi_today(double * pvecback,
                                                   double * pvectsz,
                                                   struct background * pba,
                                                   struct primordial * ppm,
                                                   struct nonlinear * pnl,
                                                   struct tszspectrum * ptsz)
{

   // double nu = exp(pvectsz[ptsz->index_lognu]);
   // double nuTink = sqrt(nu); //in T10 paper: nu=delta/sigma while here nu=(delta/sigma)^2
   //
   // double Delta;
   //
   // //Tinker et al 2008 @ M1600-mean
   // if (ptsz->MF==6)
   // Delta = 1600.;
   //
   // //Jenkins et al 2001
   // else if (ptsz->MF==3)
   // Delta = 180.;
   //
   // //Tinker et al 2008 @ m500
   // //Boquet et al 2015
   // else if (ptsz->MF==5 || ptsz->MF==7)
   // Delta = 500./pvecback[pba->index_bg_Omega_m];
   //
   // else
   // Delta = 200.;
   //
   //
   // // Table 2 of Tinker et al 2010:
   // double y = log10(Delta);
   // double ATT = 1.+0.24*y*exp(-pow(4./y,4.));
   // double aTTT = 0.44*y-0.88;
   // double BTT = 0.183;
   // double bTTT = 1.5;
   // double CTT = 0.019+0.107*y+0.19*exp(-pow(4./y,4.));
   // double cTTT = 2.4;
   //
   // pvectsz[ptsz->index_halo_bias] = 1.-ATT*(pow(nuTink,aTTT)/(pow(nuTink,aTTT)+pow(ptsz->delta_cSZ,aTTT)))
   //                                                    +BTT*pow(nuTink,bTTT)+CTT*pow(nuTink,cTTT);

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


int evaluate_HMF(double logM,
                 double * pvecback,
                 double * pvectsz,
                 struct background * pba,
                 struct nonlinear * pnl,
                 struct tszspectrum * ptsz)
{

   double z = pvectsz[ptsz->index_z];
   pvectsz[ptsz->index_dlnMdeltadlnM] = 1.;

  //mass conversions

  // Tinker et al 2008 @ M500c : the mass integral runs over m500c -> logM = logM500c
  if (ptsz->MF==5){

  pvectsz[ptsz->index_m500c] = exp(logM);
  pvectsz[ptsz->index_r500c] = pow(3.*pvectsz[ptsz->index_m500c]/(4.*_PI_*500.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc

  // //convert from m500c to mVIR:
  class_call(mDEL_to_mVIR(pvectsz[ptsz->index_m500c],
                           500.*(pvectsz[ptsz->index_Rho_crit]),
                           pvectsz[ptsz->index_Delta_c],
                           pvectsz[ptsz->index_Rho_crit],
                           z,
                           &pvectsz[ptsz->index_mVIR],
                           ptsz),
                  ptsz->error_message,
                  ptsz->error_message);
   // pvectsz[ptsz->index_mVIR] = 0.;
  }
  // for the other HMF the mass integral runs over mVIR
  else {
   //beginning of mass function
   pvectsz[ptsz->index_mVIR] = exp(logM);
}

   pvectsz[ptsz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[ptsz->index_mVIR],pvectsz[ptsz->index_Delta_c],pvectsz[ptsz->index_Rho_crit],ptsz);
   //compute concentration_parameter using mVIR

   pvectsz[ ptsz->index_cVIR] = evaluate_cvir_of_mvir(pvectsz[ptsz->index_mVIR],z,ptsz);

   //Scale radius:
   //(note that rs is actually r200c/c200 for SC14 concentration-mass relation ---> TBC)
   pvectsz[ptsz->index_rs] = pvectsz[ptsz->index_rVIR]/pvectsz[ptsz->index_cVIR];
   pvectsz[ptsz->index_ls] = pvecback[pba->index_bg_ang_distance]*pba->h/pvectsz[ptsz->index_rs];


   //Compute m500c
   //for the pressure profile
   //(except when HMF is at m500c, i.e., T08@M500 and B16M500c)
   if (ptsz->MF!=5 && ptsz->MF!=7){

      class_call(mVIR_to_mDEL(pvectsz[ptsz->index_mVIR],
                           pvectsz[ptsz->index_rVIR],
                           pvectsz[ptsz->index_cVIR],
                           500.*(pvectsz[ptsz->index_Rho_crit]),
                           &pvectsz[ptsz->index_m500c],
                           ptsz),
                      ptsz->error_message,
                      ptsz->error_message);

   }


if (ptsz->pressure_profile == 4){
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
    pvectsz[ptsz->index_l200c] = pvecback[pba->index_bg_ang_distance]*pba->h/pvectsz[ptsz->index_r200c];
}

   // //if HMF=T08@m500c -> m200=m500=M
   // else pvectsz[ptsz->index_m500] = exp(logM);

   if (ptsz->mass_dependent_bias == 1)
      ptsz->HSEbias = 1./(0.8/(1.+ ptsz->Ap*pow(pvectsz[ptsz->index_m500]/3.e14,ptsz->alpha_b)));


  //m500 X-ray for the pressure profiles A10 and P13
   pvectsz[ptsz->index_m500] = pvectsz[ptsz->index_m500c]/ptsz->HSEbias;

   //r500 X-ray for the pressure profiles A10 and P13
   pvectsz[ptsz->index_r500] = pow(3.*pvectsz[ptsz->index_m500]/(4.*_PI_*500.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc
   pvectsz[ptsz->index_l500] = pvecback[pba->index_bg_ang_distance]*pba->h/pvectsz[ptsz->index_r500];




  // compute m180m
  //Jenkins et al 2001 @ m180-mean
  //Jenkins et al 2001
  if (ptsz->MF==3){
     class_call(mVIR_to_mDEL(pvectsz[ptsz->index_mVIR],
                          pvectsz[ptsz->index_rVIR],
                          pvectsz[ptsz->index_cVIR],
                          180.*( pvecback[pba->index_bg_Omega_m])
                          *pvectsz[ptsz->index_Rho_crit],
                          &pvectsz[ptsz->index_m180m],
                          ptsz),
                     ptsz->error_message,
                     ptsz->error_message);}


  // compute m1600m
  //Tinker et al 2008 @ m1600-mean
  else if (ptsz->MF==6){
   class_call(mVIR_to_mDEL(pvectsz[ptsz->index_mVIR],
                        pvectsz[ptsz->index_rVIR],
                        pvectsz[ptsz->index_cVIR],
                        1600.*( pvecback[pba->index_bg_Omega_m])
                        *pvectsz[ptsz->index_Rho_crit],
                        &pvectsz[ptsz->index_m1600m],
                        ptsz),
                   ptsz->error_message,
                   ptsz->error_message);}



  //compute m200m
  //Tinker et al 2008 @ m200-mean
  //Tinker et al 2010 @ m200-mean
  //Boquet et al 2015 @ m200-mean
  //convert mVIR to m200m
  //used in all cases except T08@M500 [5] and B15@M500 [7]
  else if (ptsz->MF!=5 && ptsz->MF!=7) {
     class_call(mVIR_to_mDEL(pvectsz[ptsz->index_mVIR],
                          pvectsz[ptsz->index_rVIR],
                          pvectsz[ptsz->index_cVIR],
                          200.*pvecback[pba->index_bg_Omega_m]
                          *pvectsz[ptsz->index_Rho_crit],
                          &pvectsz[ptsz->index_m200m],
                          ptsz),
                     ptsz->error_message,
                     ptsz->error_message);
    // compute r200m:
    pvectsz[ptsz->index_r200m]  = pow(3.*pvectsz[ptsz->index_m200m]
                                      /(4.*_PI_*200.*pvecback[pba->index_bg_Omega_m]
                                      *pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc

   // Duffy et al 08
   if (ptsz->concentration_parameter==0){
    evaluate_c200m_D08(pvecback,pvectsz,pba,ptsz);
   }
   else {
     if (ptsz->sz_verbose>=3)
     printf("c200m not implemented yet for this concentration mass relation\n");
   }


   //dlnm200ddlnm;
   pvectsz[ptsz->index_dlnMdeltadlnM]= evaluate_dlnMdeltadlnM(log(pvectsz[ptsz->index_mVIR]),
                                                              pvecback,
                                                              pvectsz,
                                                              pba,
                                                              pnl,
                                                              ptsz);
                                                            }



  double m_for_hmf;

  //Tinker et al 2010
  if (ptsz->MF==1)
    m_for_hmf = pvectsz[ptsz->index_m200m];

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


  pvectsz[ptsz->index_mass_for_hmf] = m_for_hmf;

  if (ptsz->HMF_prescription_NCDM == 0) //Matter
    pvectsz[ptsz->index_Rh] = pow(3.*pvectsz[ptsz->index_mass_for_hmf]/(4*_PI_*(pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0),1./3.);

  else if (ptsz->HMF_prescription_NCDM == 1) //CDM
    pvectsz[ptsz->index_Rh] = pow(3.*pvectsz[ptsz->index_mass_for_hmf]/(4*_PI_*(pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0),1./3.);

  else if (ptsz->HMF_prescription_NCDM == 2) //No-pres
    pvectsz[ptsz->index_Rh] = pow(3.*pvectsz[ptsz->index_mass_for_hmf]/(4*_PI_*ptsz->Omega_m_0*ptsz->Rho_crit_0),1./3.);



   double z_asked = log(1.+z);
   double R_asked = log(exp(log(pvectsz[ptsz->index_Rh]))/pba->h);

   if (z<exp(ptsz->array_redshift[0])-1.)
      z_asked = ptsz->array_redshift[0];
   if (z>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
      z_asked =  ptsz->array_redshift[ptsz->n_arraySZ-1];
   if (log(exp(log(pvectsz[ptsz->index_Rh]))/pba->h)<ptsz->array_radius[0])
      R_asked = ptsz->array_radius[0];
   if (log(exp(log(pvectsz[ptsz->index_Rh]))/pba->h)>ptsz->array_radius[ptsz->ndimSZ-1])
      R_asked =  ptsz->array_radius[ptsz->ndimSZ-1];


   pvectsz[ptsz->index_logSigma2] =
   exp(pwl_interp_2d(ptsz->n_arraySZ,
                     ptsz->ndimSZ,
                     ptsz->array_redshift,
                     ptsz->array_radius,
                     ptsz->array_sigma_at_z_and_R,
                     1,
                     &z_asked,
                     &R_asked));


   pvectsz[ptsz->index_logSigma2] *= pvectsz[ptsz->index_logSigma2];
   pvectsz[ptsz->index_logSigma2] = log(pvectsz[ptsz->index_logSigma2]);


   pvectsz[ptsz->index_lognu] = 2.*log(ptsz->delta_cSZ) - pvectsz[ptsz->index_logSigma2];



   pvectsz[ptsz->index_dlogSigma2dlogRh] =
   pwl_interp_2d(
                  ptsz->n_arraySZ,
                  ptsz->ndimSZ,
                  ptsz->array_redshift,
                  ptsz->array_radius,
                  ptsz->array_dsigma2dR_at_z_and_R,
                  1,
                  &z_asked,
                  &R_asked
                  );


   pvectsz[ptsz->index_dlogSigma2dlogRh] *=
   exp(log(pvectsz[ptsz->index_Rh]))/pba->h
   /exp(pvectsz[ptsz->index_logSigma2]);


   pvectsz[ptsz->index_dlognudlogRh] =
   -pvectsz[ptsz->index_dlogSigma2dlogRh];

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
      class_call(
                      MF_T08(
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
   //Tinker et al 2008 @ m500
   else if (ptsz->MF==5) {
      class_call(
                      MF_T08_m500(
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


   pvectsz[ptsz->index_dndlogRh] = 3./(4.*_PI_*pow(pvectsz[ptsz->index_Rh],3))
                                                   *pvectsz[ptsz->index_dlognudlogRh]
                                                   *pvectsz[ptsz->index_mf];

   //Return the HMF - dn/dlogM in units of h^3 Mpc^-3
   pvectsz[ptsz->index_hmf] = pvectsz[ptsz->index_dndlogRh]/3.;

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


// evaluate normalisation of galaxy terms, i.e., ng_bar
int evaluate_mean_galaxy_number_density_at_z(
                   double * pvectsz,
                   struct tszspectrum * ptsz)
  {

   double z = pvectsz[ptsz->index_z];
   double z_asked = log(1.+z);

   if (z<exp(ptsz->array_redshift[0])-1.)
      z_asked = ptsz->array_redshift[0];
   if (z>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
      z_asked =  ptsz->array_redshift[ptsz->n_arraySZ-1];


   pvectsz[ptsz->index_mean_galaxy_number_density] =  exp(pwl_value_1d(ptsz->n_arraySZ,
                                                          ptsz->array_redshift,
                                                          ptsz->array_mean_galaxy_number_density,
                                                          z_asked));

return _SUCCESS_;
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

      printf("Writing output files in %s\n",Filepath);
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


      printf("Writing output files in %s\n",Filepath);
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

   if (z<exp(ptsz->array_redshift[0])-1.)
      z_asked = ptsz->array_redshift[0];
   if (z>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
      z_asked =  ptsz->array_redshift[ptsz->n_arraySZ-1];
   if (R_asked<ptsz->array_radius[0])
      R_asked = ptsz->array_radius[0];
   if (R_asked>ptsz->array_radius[ptsz->ndimSZ-1])
      R_asked =  ptsz->array_radius[ptsz->ndimSZ-1];

    double sigma8 =  exp(pwl_interp_2d(ptsz->n_arraySZ,
                       ptsz->ndimSZ,
                       ptsz->array_redshift,
                       ptsz->array_radius,
                       ptsz->array_sigma_at_z_and_R,
                       1,
                       &z_asked,
                       &R_asked));

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
     || ptsz->has_kSZ_kSZ_gal_1halo)
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
  double Eh = pvecback[pba->index_bg_H]/ptsz->H0_in_class_units;
  double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
  double delc = Delta_c_of_Omega_m(omega);

  double mvir = pow(10.,14); //m200m for M=10^{13.5} Msun/h
  double rvir = evaluate_rvir_of_mvir(mvir,delc,rhoc,ptsz);
  double cvir = evaluate_cvir_of_mvir(mvir,z,ptsz);// 5.72*pow(mvir/1e14,-0.081)/pow(1.+z_asked,0.71);
  double mdel;
  ///double rs = rvir/cvir;

  double delrho = 200.*omega*rhoc; // 200m
  double delrho_prime = 500.*rhoc; //500c


  class_call(mVIR_to_mDEL(mvir,
                         rvir,
                         cvir,
                         delrho,
                         &mdel,
                         ptsz),
                  ptsz->error_message,
                  ptsz->error_message);

  mvir_over_m200d = mvir/mdel;

  double mvir_recovered = 1.;
  double mvir_over_mvir_recovered;



  class_call(mDEL_to_mVIR(mdel,
                         delrho,
                         delc,
                         rhoc,
                         z,
                         &mvir_recovered,
                         ptsz),
                  ptsz->error_message,
                  ptsz->error_message);


  double mdel_prime;
  class_call(mDEL_to_mDELprime(mdel,
                         delrho,
                         delrho_prime,
                         delc,
                         rhoc,
                         z,
                         &mdel_prime,
                         ptsz),
                  ptsz->error_message,
                  ptsz->error_message);
  mvir_over_mvir_recovered = mvir/mvir_recovered;

  double cvir_fac = 1.;
  double cvir_prime = cvir_fac*cvir;

  double delta= 2.5;
  double delta_prime;

  delta_prime = delta_to_delta_prime_nfw(delta,cvir,cvir_prime,ptsz);

  //
  // printf("z = %.3e mvir = %.4e mdel = %.4e mvir_rec = %.4e mdel_prime = %.4e cvir = %.4e delta = %.4e delta_prime = %.4e\n",
  //         z,mvir,mdel,mvir_recovered,mdel_prime, cvir,delta,delta_prime);


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
   double z = 0.;
   double theta = evaluate_sed_cib(z,  nu, ptsz);

   fprintf(fp,"%.5e \t %.5e\n",nu,theta);
 }

  fclose(fp);
}



  free(pvecback);


                                        }

int write_output_to_files_cl(struct nonlinear * pnl,
                             struct background * pba,
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
  + ptsz->has_kSZ_kSZ_gal_1halo
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
             ptsz->append_name_cobaya_ref,
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
             ptsz->append_name_cobaya_ref,
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
      + ptsz->has_kSZ_kSZ_gal_1halo
      + ptsz->has_kSZ_kSZ_lensmag_1halo
      + ptsz->has_gal_gal_1h
      + ptsz->has_gal_gal_2h
      + ptsz->has_gal_lens_1h
      + ptsz->has_gal_lens_2h
      + ptsz->has_tSZ_gal_1h
      + ptsz->has_tSZ_gal_2h
      + ptsz->has_tSZ_cib_1h
      + ptsz->has_tSZ_cib_2h
      + ptsz->has_lens_cib_1h
      + ptsz->has_lens_cib_2h
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
       if (ptsz->pressure_profile == 2)
          fprintf(fp,"# Pressure Profile:  Arnaud et al 2010\n");
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
                    "%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n",
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
                    ptsz->cl_lens_cib_2h[ptsz->id_nu_cib_to_save][index_l]
                    );

      }
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
          // ptsz->r_cl_clp[index_l][index_l_prime] = (ptsz->cov_Y_Y_ssc[index_l][index_l_prime])
          //                                          /sqrt(ptsz->cov_cl_cl[index_l] + ptsz->cov_Y_Y_ssc[index_l][index_l])
          //                                          /sqrt(ptsz->cov_cl_cl[index_l_prime]+ ptsz->cov_Y_Y_ssc[index_l_prime][index_l_prime]);



           ptsz->r_cl_clp[index_l_prime][index_l] = ptsz->r_cl_clp[index_l][index_l_prime];

           // if (index_l==index_l_prime)
           //   ptsz->r_cl_clp[index_l_prime][index_l] = 1.;

           // ell_prime = ptsz->ell[index_l_prime];
           // ell = ptsz->ell[index_l];
           // ptsz->trispectrum_ref[index_l][index_l_prime] = ell*(ell+1.)/(2.*_PI_)*ell_prime*(ell_prime+1.)/(2.*_PI_)*ptsz->tllprime_sz[index_l][index_l_prime];
           //
           // ptsz->trispectrum_ref[index_l_prime][index_l] = ptsz->trispectrum_ref[index_l][index_l_prime];

         };




if ((ptsz->has_pk_at_z_1h + ptsz->has_pk_at_z_2h >= _TRUE_) && ptsz->write_sz>0){
      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "halo_model_pk_at_z_k_pk1h_pk2h",
                  ".txt");

      fp=fopen(Filepath, "w");

      int index_k;
      for (index_k=0;index_k<ptsz->n_k_for_pk_hm;index_k++){
      fprintf(fp,"%.4e\t%.4e\t%.4e\n",ptsz->k_for_pk_hm[index_k],ptsz->pk_at_z_1h[index_k],ptsz->pk_at_z_2h[index_k]);
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
    fprintf(fp,"%e\t%e\t%e\t%e\n",ptsz->cov_Y_N_mass_bin_edges[index_M_bins],ptsz->cov_Y_N_mass_bin_edges[index_M_bins+1],ptsz->cov_N_N[index_M_bins],ptsz->cov_N_N_hsv[index_M_bins][index_M_bins]);

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

      ptsz->Rho_crit_0 =
      (3./(8.*_PI_*_G_*_M_sun_))
      *pow(_Mpc_over_m_,1)
      *pow(_c_,2)
      *pvecback[pba->index_bg_rho_crit]
      /pow(pba->h,2);

      //ptsz->D_0 =  pvecback[pba->index_bg_D];
      //printf("D0 = %e\n",ptsz->D_0);

      ptsz->Omega_m_0 = pvecback[pba->index_bg_Omega_m];
      ptsz->Omega_ncdm_0 = ptsz->Omega_m_0
      -pba->Omega0_b
      -pba->Omega0_cdm;

      ptsz->H0_in_class_units = pvecback[pba->index_bg_H];

      if (pba->Omega0_lambda != 0.) OmegaM = 1-pba->Omega0_lambda;
      else OmegaM = 1-pba->Omega0_fld;


      ptsz->Sigma8OmegaM_SZ = pnl->sigma8[pnl->index_pk_m]*pow(OmegaM/0.28,3./8.);


      //quantities at surface of last scattering
      // conformal distance to redshift of last scattering in Mpc/h
      ptsz->chi_star = pth->ra_star*pba->h;





      if (ptsz->sz_verbose > 0)
      {
         printf("Class_sz computations\n");

         //   double z_star;  /**< redshift at which photon optical depth crosses one */
         printf("->lss at z_star = %e\n",pth->z_star);
         //Cosmological parameters
         //printf("pba->H0 = %e, %e\n",pba->H0,ptsz->H0_in_class_units);

         // printf("Omega0_b = %e\n",pba->Omega0_b);
         // printf("Omega0_cdm = %e\n",pba->Omega0_cdm);
         // printf("Omega0_m = %e\n",ptsz->Omega_m_0 );
         // printf("Omega0_ncdm = %e\n",
         //           ptsz->Omega_m_0
         //           -pba->Omega0_b
         //           -pba->Omega0_cdm);
         // printf("Rho_crit_0 = %e in units of h^2 M_sun/Mpc^3\n",ptsz->Rho_crit_0 );


         //Halo mass function
         if (ptsz->MF==3) printf("->HMF: Jenkins et al 2001 @ M180m\n");
         if (ptsz->MF==1) printf("->HMF: Tinker et al 2010 @ M200m\n");
         if (ptsz->MF==4) printf("->HMF: Tinker et al 2008 @ M200m\n");
         if (ptsz->MF==5) printf("->HMF: Tinker et al 2008 @ M500c\n");
         if (ptsz->MF==6) printf("->HMF: Tinker et al 2008 @ M1600m\n");
         if (ptsz->MF==2) printf("->HMF: Bocquet et al 2015 @ M200m\n");
         if (ptsz->MF==7) printf("->HMF: Bocquet et al 2015 @ M500c\n\n");

         //Pressure profile
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
         // }
         printf("->h = %e\n",pba->h);
         printf("->OmegaM (all except DE/Lambda) = %e\n",OmegaM);
         printf("->OmegaL = %e\n",1.-OmegaM);
         printf("->sigma8 = %e\n",pnl->sigma8[pnl->index_pk_m]);
         printf("->Bias B = %e\n",ptsz->HSEbias);printf("->Bias b = %e\n",1.-1./ptsz->HSEbias);
      }

 // }

   //class_call(spectra_sigma_ncdm(pba,
  //                                              ppm,
  //                                              pnl,
  //                                              8./pba->h,
  //                                              0.,&(ptsz->sigma8_Pcb)),
  //                 pnl->error_message,
  //                 pnl->error_message);


  ptsz->sigma8_Pcb = pnl->sigma8[pnl->index_pk_cb];
   if (ptsz->sz_verbose > 0)
      printf("->sigma8_cb= %e\n",ptsz->sigma8_Pcb);

   //compute T_cmb*gNU at 150GHz
   double frequency_in_Hz = 150.0e9;
   //double frequency_in_Hz = 0.0;
   //double frequency_in_Hz = 100.0e9;
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
                      /2.))
       -4.);
   if (ptsz->sz_verbose > 0){
       printf("->Tcmb_gNU at 150GHz= %e Kelvins\n",ptsz->Tcmb_gNU_at_150GHz);
   //    printf("->gNU at 150GHz= %e\n",ptsz->Tcmb_gNU/pba->T_cmb);
   //    printf("->Tcmb = %e K\n",pba->T_cmb);
   //
    }


    ptsz->Omega_survey = 4.*_PI_*ptsz->f_sky;
return _SUCCESS_;
}




int show_results(struct background * pba,
                         struct nonlinear * pnl,
                         struct primordial * ppm,
                         struct tszspectrum * ptsz){

  //printf("-> SZ results:\n");
  if (ptsz->has_sz_ps){
   int index_l;
   for (index_l=0;index_l<ptsz->nlSZ;index_l++){

   printf("ell = %e\t\t C_ell (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_sz_1h[index_l]);


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

 if (ptsz->has_mean_y)
   printf("mean_y =  %e \n",ptsz->y_monopole);
if (ptsz->has_hmf){
   printf("N_tot =  %e (per steradian)\n",ptsz->hmf_int);
   printf("N_tot =  %e (full-sky)\n",ptsz->hmf_int*3.046174198e-4*41253.0);
 }

 if (ptsz->has_kSZ_kSZ_gal_1halo){
  int index_l;
  for (index_l=0;index_l<ptsz->nlSZ;index_l++){

  printf("ell = %e\t\t cl_kSZ_kSZ_gal_1h (1h) = %e \n",ptsz->ell[index_l],ptsz->cl_kSZ_kSZ_gal_1h[index_l]);
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



   class_alloc(ptsz->ell_mock,
                     10*sizeof(double),
                     ptsz->error_message);




   if(ptsz->ell_sz==4){
      ptsz->nlSZ = (int)((log(ptsz->ell_max_mock) - log(ptsz->ell_min_mock))/ptsz->dlogell) + 1;
      class_realloc(ptsz->ell_mock,
                    ptsz->ell_mock,
                    ptsz->nlSZ*sizeof(double),
                    ptsz->error_message);

      int i;
      for (i=0;i<ptsz->nlSZ;i++)
       ptsz->ell_mock[i] = exp(log(ptsz->ell_min_mock)+i*ptsz->dlogell);
   }

  if (ptsz->has_kSZ_kSZ_gal_1halo
   || ptsz->has_kSZ_kSZ_lensmag_1halo){
  class_alloc(ptsz->ell_kSZ2_gal_multipole_grid,
              ptsz->N_kSZ2_gal_multipole_grid*sizeof(double),
              ptsz->error_message);

  int index_l;
  for (index_l=0;index_l<ptsz->N_kSZ2_gal_multipole_grid;index_l++){
  // ptsz->ell_kSZ2_gal_multipole_grid[index_l] = exp(log(2.) + index_l*(log(5.*ptsz->ell_max_mock)
  //                                                  - log(2.))/(ptsz->N_kSZ2_gal_multipole_grid-1.));

   ptsz->ell_kSZ2_gal_multipole_grid[index_l] = exp(log(2.) + index_l*(log(4e3)
                                                    - log(2.))/(ptsz->N_kSZ2_gal_multipole_grid-1.));

  }

  }







   class_alloc(ptsz->ell,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   int index_l;
   for (index_l=0;index_l<ptsz->nlSZ;index_l++)
   {
      if (ptsz->ell_sz == 0)
         ptsz->ell[index_l] = pow(10.,1.+1.*index_l*0.2);
      else if (ptsz->ell_sz == 1)
         ptsz->ell[index_l] = ptsz->ell_plc[index_l];
      else if (ptsz->ell_sz == 2)
         ptsz->ell[index_l] = ptsz->ell_trispectrum[index_l];
      else if (ptsz->ell_sz == 3)
         ptsz->ell[index_l] = ptsz->ell_plc_low[index_l];
      else if (ptsz->ell_sz == 4)
         ptsz->ell[index_l] = ptsz->ell_mock[index_l];
     else if (ptsz->ell_sz == 5)
        ptsz->ell[index_l] = ptsz->ell_plc_no_low_ell[index_l];
   }


   free(ptsz->ell_trispectrum);
   free(ptsz->ell_plc);
   free(ptsz->ell_plc_low);
   free(ptsz->ell_plc_no_low_ell);
   free(ptsz->ell_mock);

   //free(ptsz->cib_frequency_list);


   return _SUCCESS_;

}



int initialise_and_allocate_memory(struct tszspectrum * ptsz){
   //printf("cib_dim = %d\n",1000);


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



   class_alloc(ptsz->x_gauss, 6*sizeof(double),ptsz->error_message);
   class_alloc(ptsz->w_gauss, 6*sizeof(double),ptsz->error_message);

   ptsz->x_gauss[0]=0.0;
   ptsz->x_gauss[1]=0.1488743389;
   ptsz->x_gauss[2]=0.4333953941;
   ptsz->x_gauss[3]=0.6794095682;
   ptsz->x_gauss[4]=0.8650633666;
   ptsz->x_gauss[5]=0.9739065285;

   ptsz->w_gauss[0]=0.0;
   ptsz->w_gauss[1]=0.2955242247;
   ptsz->w_gauss[2]=0.2692667193;
   ptsz->w_gauss[3]=0.2190863625;
   ptsz->w_gauss[4]=0.1494513491;
   ptsz->w_gauss[5]=0.0666713443;


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
   ptsz->index_multipole_for_lensing_profile = 0;
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
   ptsz->index_rVIR  = ptsz->index_Delta_c +1;
   ptsz->index_cVIR  = ptsz->index_rVIR +1;
   ptsz->index_mVIR  = ptsz->index_cVIR +1;
   ptsz->index_m500  = ptsz->index_mVIR +1;
   ptsz->index_r500  = ptsz->index_m500 +1;
   ptsz->index_l500  = ptsz->index_r500 +1;
   ptsz->index_m200  = ptsz->index_l500 +1;
   ptsz->index_m200m  = ptsz->index_m200 +1;
   ptsz->index_m1600m  = ptsz->index_m200m +1;
   ptsz->index_m180m  = ptsz->index_m1600m +1;
   ptsz->index_mass_for_hmf  = ptsz->index_m180m +1;
   ptsz->index_m500c  = ptsz->index_mass_for_hmf +1;
   ptsz->index_r500c  = ptsz->index_m500c +1;
   ptsz->index_Rh  = ptsz->index_r500c +1;
   ptsz->index_mf  = ptsz->index_Rh +1;
   ptsz->index_dlognudlogRh  = ptsz->index_mf +1;
   ptsz->index_lognu  = ptsz->index_dlognudlogRh +1;
   ptsz->index_dlogSigma2dlogRh  = ptsz->index_lognu +1;
   ptsz->index_dndlogRh  = ptsz->index_dlogSigma2dlogRh +1;
   ptsz->index_logSigma2  = ptsz->index_dndlogRh +1;
   ptsz->index_z = ptsz->index_logSigma2 +1;
   ptsz->index_rs = ptsz->index_z +1;
   ptsz->index_ls = ptsz->index_rs +1;
   ptsz->index_m200c = ptsz->index_ls +1;
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
   ptsz->index_halo_bias = ptsz->index_hmf + 1;
   ptsz->index_k_value_for_halo_bias = ptsz->index_halo_bias +1;
   ptsz->index_pk_for_halo_bias = ptsz->index_k_value_for_halo_bias +1;
   ptsz->index_mass_bin_1 = ptsz->index_pk_for_halo_bias + 1;
   ptsz->index_mass_bin_2 = ptsz->index_mass_bin_1 + 1;
   ptsz->index_multipole_1 = ptsz->index_mass_bin_2  + 1;
   ptsz->index_multipole_2 = ptsz->index_multipole_1 + 1;
   ptsz->index_sigma2_hsv = ptsz->index_multipole_2 + 1;
   ptsz->index_part_id_cov_hsv = ptsz->index_sigma2_hsv + 1;
   ptsz->index_redshift_for_dndlnM = ptsz->index_part_id_cov_hsv + 1;
   ptsz->index_mass_for_dndlnM = ptsz->index_redshift_for_dndlnM +1;

   ptsz->index_phi_galaxy_counts = ptsz->index_mass_for_dndlnM + 1;
   ptsz->index_mean_galaxy_number_density = ptsz->index_phi_galaxy_counts+1;
   ptsz->index_c500c_KA20 = ptsz->index_mean_galaxy_number_density+1;
   ptsz->index_multipole_for_galaxy_profile = ptsz->index_c500c_KA20+1;
   ptsz->index_galaxy_profile =  ptsz->index_multipole_for_galaxy_profile + 1;
   ptsz->index_multipole_for_cib_profile = ptsz->index_galaxy_profile + 1;
   ptsz->index_cib_profile = ptsz->index_multipole_for_cib_profile + 1;
   ptsz->index_frequency_for_cib_profile = ptsz->index_cib_profile + 1;
   ptsz->index_frequency_prime_for_cib_profile = ptsz->index_frequency_for_cib_profile + 1;
   ptsz->index_multipole_3 =  ptsz->index_frequency_prime_for_cib_profile + 1;
   ptsz->index_W_lensmag = ptsz->index_multipole_3 + 1;
   ptsz->index_c200m = ptsz->index_W_lensmag + 1;
   ptsz->index_r200m = ptsz->index_c200m + 1;
   ptsz->index_k_for_pk_hm = ptsz->index_r200m + 1;
   ptsz->index_density_profile =  ptsz->index_k_for_pk_hm + 1;
   //quantities integrated over redshift
   ptsz->index_integral =  ptsz->index_density_profile + 1;
   ptsz->index_integral_over_m = ptsz->index_integral+1;




   //integrands at m and z
   ptsz->index_integrand =  ptsz->index_integral + 1;

   //final size of pvecsz vector
   ptsz->tsz_size  = ptsz->index_integrand + 1;


//
   //printf("cib_dim = %d\n",1000);

   if (ptsz->has_cib_cib_1h
     + ptsz->has_cib_cib_2h
     + ptsz->has_tSZ_cib_1h
     + ptsz->has_tSZ_cib_2h
     + ptsz->has_lens_cib_1h
     + ptsz->has_lens_cib_2h
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


   if (ptsz->has_pk_at_z_1h + ptsz->has_pk_at_z_2h >= _TRUE_){
     class_alloc(ptsz->k_for_pk_hm,
                       10*sizeof(double),
                       ptsz->error_message);
     ptsz->n_k_for_pk_hm = (int)((log(ptsz->k_max_for_pk_hm) - log(ptsz->k_min_for_pk_hm))/ptsz->dlnk_for_pk_hm) + 1;
     class_realloc(ptsz->k_for_pk_hm,
                   ptsz->k_for_pk_hm,
                   ptsz->n_k_for_pk_hm*sizeof(double),
                   ptsz->error_message);
                   int i;

   class_alloc(ptsz->pk_at_z_1h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);
   class_alloc(ptsz->pk_at_z_2h,sizeof(double *)*ptsz->n_k_for_pk_hm,ptsz->error_message);

   for (i=0;i<ptsz->n_k_for_pk_hm;i++){
    ptsz->k_for_pk_hm[i] = exp(log(ptsz->k_min_for_pk_hm)+i*ptsz->dlnk_for_pk_hm);
    ptsz->pk_at_z_1h[i] = 0.;
    ptsz->pk_at_z_2h[i] = 0.;
   }

   }

//exit(0);


   ptsz->hmf_int = 0.;
   ptsz->y_monopole = 0.;


   class_alloc(ptsz->cl_sz_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_isw_lens,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_isw_tsz,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_isw_auto,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_tSZ_gal_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_tSZ_gal_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   //class_alloc(ptsz->cl_tSZ_cib_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   //class_alloc(ptsz->cl_tSZ_cib_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);

   class_alloc(ptsz->cl_gal_gal_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_gal_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_lens_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_gal_lens_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_lens_lens_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_lens_lens_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_tSZ_lens_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_tSZ_lens_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_gal_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_kSZ_kSZ_lensmag_1h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->b_tSZ_tSZ_tSZ_1halo,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
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
      //ptsz->cl_tSZ_cib_1h[index_l] = 0.;
      //ptsz->cl_tSZ_cib_2h[index_l] = 0.;
      //ptsz->cl_cib_cib_1h[index_l] = 0.;
      //ptsz->cl_cib_cib_2h[index_l] = 0.;
      ptsz->cl_gal_lens_1h[index_l] = 0.;
      ptsz->cl_gal_lens_2h[index_l] = 0.;
      ptsz->cl_tSZ_gal_1h[index_l] = 0.;
      ptsz->cl_tSZ_gal_2h[index_l] = 0.;
      ptsz->cl_tSZ_lens_1h[index_l] = 0.;
      ptsz->cl_lens_lens_1h[index_l] = 0.;
      ptsz->cl_lens_lens_2h[index_l] = 0.;
      ptsz->cl_tSZ_lens_2h[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_gal_1h[index_l] = 0.;
      ptsz->cl_kSZ_kSZ_lensmag_1h[index_l] = 0.;
      ptsz->b_tSZ_tSZ_tSZ_1halo[index_l] = 0.;
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
   for (index_nu=0;index_nu<ptsz->cib_frequency_list_num;index_nu++){
     class_alloc(ptsz->cl_cib_cib_1h[index_nu],sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
     class_alloc(ptsz->cl_cib_cib_2h[index_nu],sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
        for (index_nu_prime=0;index_nu_prime<ptsz->cib_frequency_list_num;index_nu_prime++){
          class_alloc(ptsz->cl_cib_cib_1h[index_nu][index_nu_prime],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
          class_alloc(ptsz->cl_cib_cib_2h[index_nu][index_nu_prime],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
                  for (index_l=0;index_l<ptsz->cib_frequency_list_num;index_l++){
     ptsz->cl_cib_cib_1h[index_nu][index_nu_prime][index_l] = 0.;
     ptsz->cl_cib_cib_2h[index_nu][index_nu_prime][index_l] = 0.;
   }
   }
   }


   class_alloc(ptsz->cl_tSZ_cib_1h,sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_tSZ_cib_2h,sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_lens_cib_1h,sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
   class_alloc(ptsz->cl_lens_cib_2h,sizeof(double **)*ptsz->cib_frequency_list_num,ptsz->error_message);
   for (index_nu=0;index_nu<ptsz->cib_frequency_list_num;index_nu++){
     class_alloc(ptsz->cl_tSZ_cib_1h[index_nu],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
     class_alloc(ptsz->cl_tSZ_cib_2h[index_nu],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
     class_alloc(ptsz->cl_lens_cib_1h[index_nu],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
     class_alloc(ptsz->cl_lens_cib_2h[index_nu],sizeof(double *)*ptsz->nlSZ,ptsz->error_message);

for (index_l=0;index_l<ptsz->cib_frequency_list_num;index_l++){
  ptsz->cl_tSZ_cib_1h[index_nu][index_l] = 0.;
  ptsz->cl_tSZ_cib_2h[index_nu][index_l] = 0.;
  ptsz->cl_lens_cib_1h[index_nu][index_l] = 0.;
  ptsz->cl_lens_cib_2h[index_nu][index_l] = 0.;

}
   }



   int last_index_integrand_id = 0;

   ptsz->index_integrand_id_dndlnM_first  = 0;
   ptsz->index_integrand_id_dndlnM_last = ptsz->index_integrand_id_dndlnM_first + ptsz->N_redshift_dndlnM*ptsz->N_mass_dndlnM - 1;
   ptsz->index_integrand_id_hmf = ptsz->index_integrand_id_dndlnM_last + 1;
   ptsz->index_integrand_id_mean_y = ptsz->index_integrand_id_hmf + 1;
   ptsz->index_integrand_id_sz_ps_first = ptsz->index_integrand_id_mean_y + 1;
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
   ptsz->index_integrand_id_kSZ_kSZ_gal_1halo_first = ptsz->index_integrand_id_cov_Y_N_next_order_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_gal_1halo_last = ptsz->index_integrand_id_kSZ_kSZ_gal_1halo_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_first =  ptsz->index_integrand_id_kSZ_kSZ_gal_1halo_last + 1;
   ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_last = ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_first = ptsz->index_integrand_id_kSZ_kSZ_lensmag_1halo_last + 1;
   ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_last = ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_tSZ_gal_1h_first = ptsz->index_integrand_id_tSZ_tSZ_tSZ_1halo_last + 1;
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

   ptsz->index_integrand_id_gal_gal_1h_first = last_index_integrand_id + 1;
   ptsz->index_integrand_id_gal_gal_1h_last = ptsz->index_integrand_id_gal_gal_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_gal_gal_2h_first = ptsz->index_integrand_id_gal_gal_1h_last + 1;
   ptsz->index_integrand_id_gal_gal_2h_last = ptsz->index_integrand_id_gal_gal_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_gal_lens_1h_first = ptsz->index_integrand_id_gal_gal_2h_last + 1;
   ptsz->index_integrand_id_gal_lens_1h_last = ptsz->index_integrand_id_gal_lens_1h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_cib_cib_2h_first = ptsz->index_integrand_id_gal_lens_1h_last + 1;
   ptsz->index_integrand_id_cib_cib_2h_last = ptsz->index_integrand_id_cib_cib_2h_first + ptsz->cib_dim - 1;
   ptsz->index_integrand_id_tSZ_cib_1h_first = ptsz->index_integrand_id_cib_cib_2h_last + 1;
   ptsz->index_integrand_id_tSZ_cib_1h_last = ptsz->index_integrand_id_tSZ_cib_1h_first + ptsz->nlSZ*ptsz->cib_frequency_list_num - 1;
   ptsz->index_integrand_id_lens_cib_1h_first = ptsz->index_integrand_id_tSZ_cib_1h_last + 1;
   ptsz->index_integrand_id_lens_cib_1h_last = ptsz->index_integrand_id_lens_cib_1h_first + ptsz->nlSZ*ptsz->cib_frequency_list_num - 1;
   ptsz->index_integrand_id_cib_cib_1h_first = ptsz->index_integrand_id_lens_cib_1h_last + 1;
   ptsz->index_integrand_id_cib_cib_1h_last = ptsz->index_integrand_id_cib_cib_1h_first + ptsz->cib_dim - 1;
   ptsz->index_integrand_id_tSZ_cib_2h_first = ptsz->index_integrand_id_cib_cib_1h_last + 1;
   ptsz->index_integrand_id_tSZ_cib_2h_last = ptsz->index_integrand_id_tSZ_cib_2h_first + ptsz->nlSZ*ptsz->cib_frequency_list_num - 1;
   ptsz->index_integrand_id_lens_cib_2h_first = ptsz->index_integrand_id_tSZ_cib_2h_last + 1;
   ptsz->index_integrand_id_lens_cib_2h_last = ptsz->index_integrand_id_lens_cib_2h_first + ptsz->nlSZ*ptsz->cib_frequency_list_num - 1;
   ptsz->index_integrand_id_gal_lens_2h_first = ptsz->index_integrand_id_lens_cib_2h_last + 1;
   ptsz->index_integrand_id_gal_lens_2h_last = ptsz->index_integrand_id_gal_lens_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integrand_id_lens_lens_1h_first = ptsz->index_integrand_id_gal_lens_2h_last + 1;
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

}

double evaluate_dlnMdeltadlnM(double logM,
                             double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct nonlinear * pnl,
                             struct tszspectrum * ptsz)
                         {
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

  cvir1=evaluate_cvir_of_mvir(mvir1,z,ptsz);//5.72*pow(mvir1/1e14,-0.081)/pow(1.+z,0.71);
  cvir2=evaluate_cvir_of_mvir(mvir2,z,ptsz);//5.72*pow(mvir2/1e14,-0.081)/pow(1.+z,0.71);
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
evaluate_galaxy_number_counts(pvecback,pvectsz,pba,ptsz);
//  normalized galaxy redshift distribution in the bin:
double phi_galaxy_at_z = pvectsz[ptsz->index_phi_galaxy_counts];
// H_over_c_in_h_over_Mpc = dz/dChi
// phi_galaxy_at_z = dng/dz normalized
double result = H_over_c_in_h_over_Mpc*phi_galaxy_at_z;




return result;
}

double evaluate_galaxy_number_counts( double * pvecback,
                                    double * pvectsz,
                                    struct background * pba,
                                    struct tszspectrum * ptsz){


    double z_asked  = pvectsz[ptsz->index_z];
    double phig = 0.;

  if(z_asked<ptsz->normalized_dndz_z[0])
     phig = 1e-100;
  else if (z_asked>ptsz->normalized_dndz_z[ptsz->normalized_dndz_size-1])
     phig = 1e-100;
  else  phig =  pwl_value_1d(ptsz->normalized_dndz_size,
                               ptsz->normalized_dndz_z,
                               ptsz->normalized_dndz_phig,
                               z_asked);

pvectsz[ptsz->index_phi_galaxy_counts] = phig;

                                    }


double HOD_mean_number_of_central_galaxies(double z,
                                           double M_halo,
                                           double M_min,
                                           double sigma_lnM,
                                           double * pvectsz,
                                           struct tszspectrum * ptsz){
 double result = 0.;

 int index_md = (int) pvectsz[ptsz->index_md];
 if (_cib_cib_1h_
   || _cib_cib_2h_
   || _tSZ_cib_1h_
   || _tSZ_cib_2h_
   || _lens_cib_1h_
   || _lens_cib_2h_){
   //printf("cib\n");
 if (M_halo>=M_min) result = 1.;
 else result = 0.;
 }
 else {
 if (ptsz->galaxy_sample == 1){ // KFSW20
 M_min = evaluate_unwise_m_min_cut(z,ptsz->unwise_galaxy_sample_id,ptsz);
 result = 0.5*(1.+gsl_sf_erf((log10(M_halo/M_min)/(sqrt(2.)*ptsz->sigma_lnM_HOD))));
 }

 else {
 result = 0.5*(1.+gsl_sf_erf((log10(M_halo/M_min)/sigma_lnM)));
 }

}


 return result;
}


double HOD_mean_number_of_satellite_galaxies(double z,
                                             double M_halo,
                                             double Nc_mean,
                                             double M_min,
                                             double alpha_s,
                                             double M1_prime,
                                             struct tszspectrum * ptsz){
double result =  0.;

   if (ptsz->galaxy_sample == 1){ // KFSW20
   M_min = evaluate_unwise_m_min_cut(z,ptsz->unwise_galaxy_sample_id,ptsz);
   if (M_halo>ptsz->M_min_HOD_satellite_mass_factor_unwise*M_min){
   result = pow((M_halo-ptsz->M_min_HOD_satellite_mass_factor_unwise*M_min)/(ptsz->M1_prime_HOD_factor*M_min),alpha_s);
    }
 else result = 0.;
}// end KFSW20
   else {
  if (M_halo>M_min){
   result = Nc_mean*pow((M_halo-M_min)/M1_prime,alpha_s);
 }
 else result = 0.;
  }

return result;
}


int evaluate_cib_profile(double * pvecback,
                         double * pvectsz,
                         struct background * pba,
                         struct tszspectrum * ptsz){


double M_halo = pvectsz[ptsz->index_mass_for_hmf]/pba->h; // convet to Msun
//printf("%.3e %.3e\n",M_halo,pvectsz[ptsz->index_mVIR]);

//double frequency_for_cib_profile = pvectsz[ptsz->index_frequency_for_cib_profile]; // in GHz


double Lc_nu;
double Lc_nu_prime;
double Ls_nu;
double Ls_nu_prime;

double z = pvectsz[ptsz->index_z];

pvectsz[ptsz->index_multipole_for_galaxy_profile] = pvectsz[ptsz->index_multipole_for_cib_profile];;
double us = evaluate_truncated_nfw_profile(pvecback,pvectsz,pba,ptsz);

double ug_at_ell;
double nu;
double nu_prime;

int index_md = (int) pvectsz[ptsz->index_md];
int index_nu = (int) pvectsz[ptsz->index_frequency_for_cib_profile];
int index_nu_prime = (int) pvectsz[ptsz->index_frequency_prime_for_cib_profile];

//printf("%.3e %.3e\n",ptsz->cib_frequency_list[index_nu],ptsz->cib_frequency_list[index_nu_prime]);

// index_nu = 0;
// index_nu_prime = 0;


// 2-halo terms
if (_tSZ_cib_2h_
  ||_lens_cib_2h_
  ||_cib_cib_2h_
   ){

if (_cib_cib_2h_){
if (index_nu != index_nu_prime){
if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  1){
nu = ptsz->cib_frequency_list[index_nu];
Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz);
Ls_nu = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu,pba,ptsz);
}
else if ((int) pvectsz[ptsz->index_part_id_cov_hsv] ==  2){
nu = ptsz->cib_frequency_list[index_nu_prime];
Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz);
Ls_nu = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu_prime,pba,ptsz);
}
}
else{
nu = ptsz->cib_frequency_list[index_nu];
Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz);
Ls_nu = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu,pba,ptsz);
}
}
// cross terms
else {
nu = ptsz->cib_frequency_list[index_nu];
Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz);
Ls_nu = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu,pba,ptsz);
}

// eq. 13 of MM20
ug_at_ell  = 1./(4.*_PI_)*(Lc_nu+Ls_nu*us);

}

// 1-halo terms
else if(_tSZ_cib_1h_
     || _lens_cib_1h_
     || _cib_cib_1h_
       ){
if (_cib_cib_1h_){
nu = ptsz->cib_frequency_list[index_nu];
Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz);
//Ls_nu = Luminosity_of_satellite_galaxies(z,M_halo,nu,ptsz);
Ls_nu = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu,pba,ptsz);

nu_prime = ptsz->cib_frequency_list[index_nu_prime];
Lc_nu_prime = Luminosity_of_central_galaxies(z,M_halo,nu_prime,pvectsz,ptsz);
//Ls_nu_prime = Luminosity_of_satellite_galaxies(z,M_halo,nu,ptsz);
Ls_nu_prime = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu_prime,pba,ptsz);
// eq. 15 of MM20
ug_at_ell  = 1./(4.*_PI_)*sqrt(Ls_nu*Ls_nu_prime*us*us+Lc_nu*Ls_nu_prime*us+Lc_nu_prime*Ls_nu*us);
}

else {
  nu = ptsz->cib_frequency_list[index_nu];
  Lc_nu = Luminosity_of_central_galaxies(z,M_halo,nu,pvectsz,ptsz);
  Ls_nu = get_L_sat_at_z_and_M_at_nu(z,M_halo,index_nu,pba,ptsz);
  Ls_nu_prime = Ls_nu;
  // eq 28 of MM20
  // printf("Lc_nu = %.3e\n",Lc_nu);
  // printf("Ls_nu = %.3e\n",Ls_nu);
  // printf("us = %.3e\n",us);
  ug_at_ell  = 1./(4.*_PI_)*(Lc_nu+Ls_nu*us);
  //ug_at_ell  = 1./(4.*_PI_)*sqrt(Ls_nu*Ls_nu_prime*us*us+Lc_nu*Ls_nu_prime*us+Lc_nu_prime*Ls_nu*us);

}
}
// need to fix units:
// from Mpc^2 to (Mpc/h)^2
ug_at_ell *= pow(pba->h,2.);
//printf("cib = %.3e\n",ug_at_ell);

pvectsz[ptsz->index_cib_profile] = ug_at_ell;

}


double Luminosity_of_central_galaxies(double z,
                                      double  M_halo,
                                      double nu,
                                      double * pvectsz,
                                      struct tszspectrum * ptsz){
double result = 0.;

double L_gal = evaluate_galaxy_luminosity(z, M_halo, nu, ptsz);
double nc = HOD_mean_number_of_central_galaxies(z,M_halo,ptsz->M_min_HOD,ptsz->sigma_lnM_HOD,pvectsz,ptsz);
return result =  nc*L_gal;
                                      }

double Luminosity_of_satellite_galaxies(double z,
                                        double M_halo,
                                        double nu,
                                        struct tszspectrum * ptsz){

double result = 0.;
double lnMs_min = log(1e6);
double lnMs_max = log(1e11);
double dlnM = (lnMs_max - lnMs_min)/10.;

double L_sat = 0.;
double L_gal;
double lnMs = lnMs_min;
double M_sub;
double M_host =  M_halo;
double dNdlnMs;
while (lnMs<lnMs_max){
M_sub = exp(lnMs);
L_gal = evaluate_galaxy_luminosity(z, M_sub, nu, ptsz);
//printf("Lgal = %.3e\n",L_gal);


dNdlnMs = subhalo_hmf_dndlnMs(M_host,M_sub);
L_sat += L_gal*dNdlnMs;
lnMs += dlnM;
}
result = L_sat;
//printf("%.3e\n",result);
return result;
                                      }


double subhalo_hmf_dndlnMs(double M_host,double M_sub){
// Subhalo mass function: Equation 12 of https://iopscience.iop.org/article/10.1088/0004-637X/719/1/88/pdf
  return 0.30*pow(M_sub/M_host,-0.7)*exp(-9.9*pow(M_sub/M_host,2.5));
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
return result;
                                      }

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

double L0 = ptsz->L0_cib;
double Phi = evaluate_phi_cib(z,ptsz);
double Theta =  evaluate_sed_cib(z,nu*(1.+z),ptsz);
double Sigma = evaluate_Sigma_cib(M,ptsz);
result = L0*Phi*Sigma*Theta;
//printf("Phi =  %.3e, Theta = %.3e, Sigma = %.3e\n",Phi, Theta, Sigma);
return result;
                                      }


int evaluate_galaxy_profile(double * pvecback,
                            double * pvectsz,
                            struct background * pba,
                            struct tszspectrum * ptsz){


double M_halo = pvectsz[ptsz->index_mass_for_hmf]; // Msun_over_h
double ng_bar = pvectsz[ptsz->index_mean_galaxy_number_density];
double nc;
double ns;
double us;
double z = pvectsz[ptsz->index_z];

nc = HOD_mean_number_of_central_galaxies(z,M_halo,ptsz->M_min_HOD,ptsz->sigma_lnM_HOD,pvectsz,ptsz);
ns = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,ptsz->M_min_HOD,ptsz->alpha_s_HOD,ptsz->M1_prime_HOD,ptsz);
us = evaluate_truncated_nfw_profile(pvecback,pvectsz,pba,ptsz);

double ug_at_ell;

int index_md = (int) pvectsz[ptsz->index_md];


// 2-halo terms or 1-halo cross
if (_gal_gal_2h_
  || _gal_lens_2h_
  || _tSZ_gal_2h_
   ){

ug_at_ell  = (1./ng_bar)*(nc+ns*us);

}

// 1-halo terms
else if(_gal_gal_1h_
     || _gal_lens_1h_
     || _tSZ_gal_1h_
     || _tSZ_lens_1h_
     || _kSZ_kSZ_gal_1halo_) {

  ug_at_ell  =(1./ng_bar)*sqrt(ns*ns*us*us+2.*ns*us);
  }

pvectsz[ptsz->index_galaxy_profile] = ug_at_ell;

}

double evaluate_truncated_nfw_profile(double * pvecback,
                                      double * pvectsz,
                                      struct background * pba,
                                      struct tszspectrum * ptsz)
{
double c_delta, r_delta;
int index_md = (int) pvectsz[ptsz->index_md];

if (_cib_cib_1h_
  || _cib_cib_2h_
  || _tSZ_cib_1h_
  || _tSZ_cib_2h_
  || _lens_cib_1h_
  || _lens_cib_2h_
){
  //printf("cib\n");
  r_delta = pvectsz[ptsz->index_r200m];
  c_delta = pvectsz[ptsz->index_c200m];
}
else{
// WIxSC KA20 with T08@M500c
if (ptsz->galaxy_sample == 0 && ptsz->MF == 5) {
evaluate_c500c_KA20(pvecback,pvectsz,pba,ptsz);
c_delta = pvectsz[ptsz->index_c500c_KA20]; //Eq. 27 of KA20
r_delta = pvectsz[ptsz->index_r500c];
// printf("c500c_KA20\n");
// exit(0);
}
else if (ptsz->galaxy_sample == 0 && ptsz->MF == 1){
  r_delta = pvectsz[ptsz->index_r200m];
  c_delta = pvectsz[ptsz->index_c200m];
}
// unWISE
else if (ptsz->galaxy_sample == 1) {
// mail from alex:
// The code that we use is here:
// https://github.com/bccp/simplehod/blob/master/simplehod/simplehod.pyx
// in particular, see the get_nfw_r function within the RNGAdapter class (line 499).  I believe this is NFW truncated at the virial radius.  If the concentration is > 1, it is set to 1 (I *think* this is what the linked code is doing in lines 507-515).
 // The concentration is coming from the fitting formula of Dutton & Maccio 2014:
// https://nbodykit.readthedocs.io/en/latest/api/_autosummary/nbodykit.transform.html#nbodykit.transform.HaloConcentration
// In particular, concentration, virial radius and velocity dispersion are all generated from the mass as described in the README here:
// https://github.com/bccp/simplehod
// the factor 2.5 is to agree with the cutoff of the lensing profile

  r_delta = ptsz->x_out_truncated_nfw_profile*pvectsz[ptsz->index_rVIR];
  c_delta = pvectsz[ptsz->index_cVIR];
  if (c_delta>1.) c_delta = 1.;
}

else {
  // r_delta =1.17*pvectsz[ptsz->index_rs];
  // c_delta = pvectsz[ptsz->index_cVIR];
  //default like KA20, needs to run with M500 and P13
  evaluate_c500c_KA20(pvecback,pvectsz,pba,ptsz);
  c_delta = pvectsz[ptsz->index_c500c_KA20]; //Eq. 27 of KA20
  r_delta = pvectsz[ptsz->index_r500c];
}
}




double z = pvectsz[ptsz->index_z];

double ell = pvectsz[ptsz->index_multipole_for_galaxy_profile];
double chi = sqrt(pvectsz[ptsz->index_chi2]);
double k = (ell+0.5)/chi;

double q = k*r_delta/c_delta*(1.+z);//TBC: (1+z) needs to be there to match KA20,  but is it consistent?
double denominator = (log(1.+c_delta)-c_delta/(1.+c_delta));


double numerator = cos(q)*(gsl_sf_Ci((1.+c_delta)*q)-gsl_sf_Ci(q))
                   +sin(q)*(gsl_sf_Si((1.+c_delta)*q)-gsl_sf_Si(q))
                   -sin(c_delta*q)/((1.+c_delta)*q);


return numerator/denominator;
}


int evaluate_c200m_D08(double * pvecback,
                        double * pvectsz,
                        struct background * pba,
                        struct tszspectrum * ptsz)
{
double M = pvectsz[ptsz->index_m200m]; // mass in  Msun/h
double A = 10.14;
double B = -0.081;
double C = -1.01;
double M_pivot = 2.e12; // pivot mass in Msun/h

double z = pvectsz[ptsz->index_z];
double c200m =A*pow(M/M_pivot,B)*pow(1.+z,C);
pvectsz[ptsz->index_c200m] = c200m;
}



int evaluate_c500c_KA20(double * pvecback,
                        double * pvectsz,
                        struct background * pba,
                        struct tszspectrum * ptsz)
{
double M = pvectsz[ptsz->index_m500c];
double A = 3.67;
double B = -0.0903;
double C = -0.51;
double M_pivot = 2.78164e12*pba->h; // pivot mass in Msun/h

double z = pvectsz[ptsz->index_z];
double c500 =A*pow(M/M_pivot,B)*pow(1.+z,C);
pvectsz[ptsz->index_c500c_KA20] = c500;
}


// Mass cuts for unwise galaxy
// This functions returns an m_cut value depending on the sample (red, blue or gree)
// We typically use it for M200m, but this is an arbitrary choice
// We restrict the mass integral to go from m_cut to the max set in the param file ('M2SZ', typically 1e15Msun/h)

// def(mcut, zz, isamp, green_option='default'):
//     '''Gives mcut for 5-param Zheng HOD for wise sample isamp at redshift zz.
//     Satellite fractions 5-10% for red and green, and 25% for blue.
//     green_option = 'default' or 'shallower'; 'default' is the default
//     bias evolution, and shallower is a somewhat shallower bias evolution,
//     but with higher masses.'''
//     if isamp=='red':
//         zall = [ 0.75,  1.00,  1.5,  2.0]
//         mcut_all = [12.00, 12.00, 12.6, 13.6]
//         if zz <= 0.75:
//             mcut = 12.00
//         elif (zz > 0.75) & (zz <= 2.00):
//             mcut = np.interp(zz,zall,mcut_all)
//         elif zz >= 2.00:
//             mcut = 13.6
//     elif isamp=='green':
//         if green_option == 'default':
//             zall = [0.00, 0.25, 0.4, 0.5, 0.65, 0.75, 1.00, 1.25, 1.50, 2.00, 2.50]
//             mcut_all = [11.9, 12.0, 12.15, 12.15, 11.75, 11.75, 12.4, 12.6, 12.75, 13.25, 13.25]
//             if zz <= 2.5:
//                mcut = np.interp(zz, zall,mcut_all)
//             else:
//                 mcut = 13.55
//         elif green_option == 'shallower':
//             zall = [0.25, 0.4, 0.5, 0.65, 0.75, 1.00]
//             mcut_all = [11.5, 12, 12, 11, 11, 12.72]
//             if zz <= 1.0:
//                 mcut = np.interp(zz, zall,mcut_all)
//             else:
//                 mcut = -0.5161*zz**4+2.919*zz**3-5.384*zz**2+\
//                            3.842*zz+12.01 if zz < 2.5 else 13.42
//     elif isamp=='blue':
//         mcut = 11.65 + zz
//     return mcut

double evaluate_unwise_m_min_cut(double z,
                                 int sample_id,
                                 struct tszspectrum * ptsz)
{
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

m_cut = pow(10.,m_cut);


// adjust the mass
double mass_fac = ptsz->M_min_HOD_mass_factor_unwise;
// if (sample_id == 0) // red
//   mass_fac = 1.7;
// else if (sample_id == 1) // green
//   mass_fac = 1.2;
// else if (sample_id == 3) // blue
//   mass_fac = 1.;

return mass_fac*m_cut;
}






double integrand_kSZ2_X_at_theta(double ln_ell_prime, void *p){
  double ell_prime  = exp(ln_ell_prime);
  //double ell_prime  = ln_ell_prime;

double integrand_cl_kSZ2_X_at_theta = 0.;
struct Parameters_for_integrand_kSZ2_X_at_theta *V = ((struct Parameters_for_integrand_kSZ2_X_at_theta *) p);



     double ell = (int) V->ptsz->ell[V->index_ell_3];
     double abs_ell_minus_ell_prime = sqrt(ell*ell+ell_prime*ell_prime-2.*ell*ell_prime*cos(V->theta));

     double ell_1 = abs_ell_minus_ell_prime;
     double ell_2 = ell_prime;
     double ell_3 = ell;

     // check bispectrum condition
     int bispec_cd;
     bispec_cd = bispectrum_condition(ell_1,ell_2,ell_3);

     if (bispec_cd == 1){

       double ln_ell1 = log(ell_1);
       double ln_ell2 = log(ell_2);
       double db =  exp(pwl_interp_2d(V->ptsz->N_kSZ2_gal_multipole_grid,
                                      V->ptsz->N_kSZ2_gal_multipole_grid,
                                      V->ln_ell,
                                      V->ln_ell,
                                      V->b_l1_l2_l_1d,
                                      1,
                                      &ln_ell1,
                                      &ln_ell2));

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


     }
      //printf("ell = %.3e ell_prime = %.3e theta = %.3e \t integrand = %.3e\n",ell,ell_prime,V->theta, integrand_cl_kSZ2_X_at_theta);
       return integrand_cl_kSZ2_X_at_theta;

}
