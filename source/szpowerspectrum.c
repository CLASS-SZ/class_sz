/** @file szpowerspectrum.c Documented SZ module.
 *
 * Boris Bolliet, 10.2017
 *
 *This module is dedicated to the computation of
 *the SZ power spectrum, following Komatsu integrand.f90
 */
#define _DEBUG


#include "szpowerspectrum.h"
#include "sz_tools.h"





int szpowerspectrum_init(
                          struct background * pba,
                          struct nonlinear * pnl,
                          struct primordial * ppm,
                          struct tszspectrum * ptsz
			                    )
{
   select_multipole_array(ptsz);
   show_preamble_messages(pba,pnl,ppm,ptsz);

   if (ptsz->experiment == 0)
      read_Planck_noise_map(ptsz);

   //obsolete:
   external_pressure_profile_init(ptsz);

   if (ptsz->concentration_parameter == 4)
      read_Zhao_CM_init(ptsz);

   tabulate_sigma_and_dsigma_from_pk(pba,pnl,ppm,ptsz);

   initialise_and_allocate_memory(ptsz);

   //SO data and Functions

   if (ptsz->experiment == 1){
      read_SO_Qfit(ptsz);
      read_SO_noise(ptsz);}


   if (ptsz->has_sz_ps + ptsz->has_hmf + ptsz->has_mean_y == _FALSE_)
   {
      if (ptsz->sz_verbose > 0)
         printf("->No SZ-y quantities requested. SZ ps module skipped.\n");
   }

   else
   {
   double * Pvecback;
   double * Pvectsz;
   int index_integrand;

   int abort;
#ifdef _OPENMP
   double tstart, tstop;
#endif
   abort = _FALSE_;
    //printf("number_of_integrands=%d\n",ptsz->number_of_integrands);

#pragma omp parallel \
   shared(abort,pba,ptsz,ppm,pnl)\
   private(tstart,tstop,Pvectsz,Pvecback,index_integrand)
	 {

#ifdef _OPENMP
	   tstart = omp_get_wtime();
#endif

	   class_alloc_parallel(Pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
       int i;
       for(i = 0; i<ptsz->tsz_size;i++) Pvectsz[i] = 0.;

	   class_alloc_parallel(Pvecback,pba->bg_size*sizeof(double),ptsz->error_message);

#pragma omp for schedule (dynamic)



//Loop over integrands
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
                     printf("In %s: time spent in parallel region (loop over l's) = %e s for thread %d\n",
                               __func__,tstop-tstart,omp_get_thread_num());


#endif
   free(Pvecback);
   free(Pvectsz);


	} //end of parallel region

   if (abort == _TRUE_) return _FAILURE_;

   ////////////////end - cl

   if (ptsz->create_ref_trispectrum_for_cobaya){

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
                    "%s%s%s",
                    ptsz->path_to_class,
                    "/sz_auxiliary_files/cobaya_class_sz_likelihoods/cobaya_reference_trispectrum/",
                    "tSZ_trispectrum_ref.txt");

        fp=fopen(Filepath, "w");

        for (index_l=0;index_l<ptsz->nlSZ;index_l++){
         for (index_l_prime=0;index_l_prime<ptsz->nlSZ;index_l_prime++) {
              fprintf(fp,"%e\t",ptsz->trispectrum_ref[index_l][index_l_prime]);
           }
           fprintf(fp,"\n");
        }
        fclose(fp);


        sprintf(Filepath,
            "%s%s%s",
            ptsz->path_to_class,
            "/sz_auxiliary_files/cobaya_class_sz_likelihoods/cobaya_reference_trispectrum/",
            "tSZ_c_ell_ref.txt");

        fp=fopen(Filepath, "w");
        for (index_l=0;index_l<ptsz->nlSZ;index_l++)
              fprintf(fp,"%e\t %e\n",ptsz->ell[index_l],ptsz->cl_sz[index_l]);
        fclose(fp);
  }

   if (ptsz->sz_verbose>0){
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
   free(ptsz->cl_sz);
   free(ptsz->cl_te_y_y);
   free(ptsz->cov_cl_cl);
   free(ptsz->cl_sz_2h);
   free(ptsz->tllprime_sz);
   free(ptsz->cov_N_cl);
   free(ptsz->r_N_cl);
   free(ptsz->r_cl_clp);
   free(ptsz->trispectrum_ref);

   free(ptsz->thetas);
   free(ptsz->skyfracs);
   free(ptsz->ylims);

   free(ptsz->SO_Qfit);
   free(ptsz->SO_thetas);


   free(ptsz->SO_RMS);
   free(ptsz->SO_skyfrac);

   free(ptsz->w_gauss);
   free(ptsz->x_gauss);

   free(ptsz->array_radius);
   free(ptsz->array_redshift);
   free(ptsz->array_sigma_at_z_and_R);
   free(ptsz->array_dsigma2dR_at_z_and_R);


   free(ptsz->PP_lnx);
   free(ptsz->PP_lnI);
   free(ptsz->PP_d2lnI);


   free(ptsz->CM_redshift);
   free(ptsz->CM_logM);
   free(ptsz->CM_logC);

   free(ptsz->M_bins);


   free(ptsz->ln_x_for_pp);


return _SUCCESS_;
}




int compute_sz(struct background * pba,
                      struct nonlinear * pnl,
                      struct primordial * ppm,
                      struct tszspectrum * ptsz,
                      double * Pvecback,
                      double * Pvectsz){


   int index_integrand = (int) Pvectsz[ptsz->index_integrand_id];
   //Pvectsz[ptsz->index_multipole] = 0.;

   if(index_integrand == ptsz->index_integrand_id_hmf){
      Pvectsz[ptsz->index_md] = ptsz->index_md_hmf;
      //printf("computing hmf int\n");
   }
   else if (index_integrand == ptsz->index_integrand_id_mean_y) {
      Pvectsz[ptsz->index_md] = ptsz->index_md_mean_y;
      //printf("computing mean y\n");
   }
   else if (index_integrand>=ptsz->index_integrand_id_sz_ps_first && index_integrand <= ptsz->index_integrand_id_sz_ps_last){
      Pvectsz[ptsz->index_md] = ptsz->index_md_sz_ps;
      Pvectsz[ptsz->index_multipole] = (double) (index_integrand - ptsz->index_integrand_id_sz_ps_first);
    }



   class_call(integrate_over_redshift_at_each_ell(pba,
                                                  pnl,
                                                  ppm,
                                                  ptsz,
                                                  Pvecback,
                                                  Pvectsz),
                   ptsz->error_message,
                   ptsz->error_message);


   int index_md = (int) Pvectsz[ptsz->index_md];
   //printf("index_md = %d\n",index_md );

   if (_hmf_){

      ptsz->hmf_int = Pvectsz[ptsz->index_integral];
      //printf("hmf_int = %e\n",ptsz->hmf_int );
   }

   if (_mean_y_){
      ptsz->y_monopole = Pvectsz[ptsz->index_integral];
      //if (ptsz->sz_verbose>0)
      //printf("mean_y = %e\n",ptsz->y_monopole );
   }


   if (_tSZ_power_spectrum_){
      int index_l = (int) Pvectsz[ptsz->index_multipole];
       ptsz->cl_sz[index_l] = Pvectsz[ptsz->index_integral];
       ptsz->cl_te_y_y[index_l] = Pvectsz[ptsz->index_integral_te_y_y];

       //if (ptsz->sz_verbose>0)
       //printf("ell = %e\tcl_sz = %e\n",ptsz->ell[index_l],ptsz->cl_sz[index_l] );

      if (_cov_N_Cl_){
         int index_l = (int) Pvectsz[ptsz->index_multipole];
         int index_M_bins;
         for (index_M_bins=0;index_M_bins<ptsz->nbins_M-1;index_M_bins++){

            ptsz->cov_N_cl[index_l][index_M_bins] = Pvectsz[ptsz->index_integral_cov_N_cl_first+index_M_bins];

            //normalised cov:
            //will be divided by cov_cl_cl at the end
            ptsz->r_N_cl[index_l][index_M_bins] = ptsz->cov_N_cl[index_l][index_M_bins]
                                                  /sqrt(Pvectsz[ptsz->index_integral_N_for_cov_N_cl_first+index_M_bins]);

            //if (ptsz->sz_verbose>0)
            //printf("ell = %e\t M = %e\t cov = %e\n",
          //            ptsz->ell[index_l],
        //              ptsz->M_bins[index_M_bins],
      //                ptsz->cov_N_cl[index_l][index_M_bins]);

         }
      }

   if (_tSZ_trispectrum_){
      int index_l = (int) Pvectsz[ptsz->index_multipole];
      int index_l_prime;
      for (index_l_prime=0; index_l_prime<index_l+1;index_l_prime++)
      {
         ptsz->tllprime_sz[index_l][index_l_prime] = Pvectsz[ptsz->index_integral_2h_first+index_l_prime];
         //if (ptsz->sz_verbose>0)
         //printf("ell = %e\tell_prime = %e\t t_llp = %e\n",
        //        ptsz->ell[index_l],ptsz->ell[index_l_prime],
        //        ptsz->tllprime_sz[index_l][index_l_prime]);

      }


   }
      if (_tSZ_2halo_){
         int index_l = (int) Pvectsz[ptsz->index_multipole];
         ptsz->cl_sz_2h[index_l] = Pvectsz[ptsz->index_integral_2halo_term];
         //if (ptsz->sz_verbose>0)
         //printf("ell = %e\tcl_sz_2h = %e\n",ptsz->ell[index_l],ptsz->cl_sz_2h[index_l] );
      }

      if (_tSZ_te_y_y_){
         int index_l = (int) Pvectsz[ptsz->index_multipole];
         ptsz->cl_te_y_y[index_l] = Pvectsz[ptsz->index_integral_te_y_y];
         //if (ptsz->sz_verbose>0)
         //printf("ell = %e\tcl_te_y_y = %e\n",ptsz->ell[index_l],ptsz->cl_te_y_y[index_l] );
      }
   }


return _SUCCESS_;
}





int integrand_at_m_and_z(double logM,
                                     double * pvecback,
                                     double * pvectsz,
                                     struct background * pba,
                                     struct primordial * ppm,
                                     struct nonlinear * pnl,
                                     struct tszspectrum * ptsz)
{

pvectsz[ptsz->index_dlnMdeltadlnM] = 1.;
//printf("before eval HMF = %e\n",pvectsz[ptsz->index_dlnMdeltadlnM] );
   evaluate_HMF(logM,pvecback,pvectsz,pba,pnl,ptsz);
//printf("after eval HMF = %e\n",pvectsz[ptsz->index_dlnMdeltadlnM] );

   double z = pvectsz[ptsz->index_z];


   //volume element in units h^-3 Mpc^3
   pvectsz[ptsz->index_volume] = pow(1.+z,2)
                                 *pow(pvecback[pba->index_bg_ang_distance]*pba->h,2)
                                 *_c_*1.e-5
                                 /(pvecback[pba->index_bg_H]/pba->H0);

   evaluate_completeness(pvecback,pvectsz,pba,ptsz);

   pvectsz[ptsz->index_multipole_for_pressure_profile] = pvectsz[ptsz->index_multipole];

   evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);




   //Return the HMF - dn/dlogM in units of h^3 Mpc^-3
   pvectsz[ptsz->index_hmf] = pvectsz[ptsz->index_dndlogRh]/3.;

      int index_md = (int) pvectsz[ptsz->index_md];

   if (_hmf_){

      pvectsz[ptsz->index_integrand] = pvectsz[ptsz->index_volume]
                                       *pvectsz[ptsz->index_hmf]
                                       *pvectsz[ptsz->index_completeness];
   }




   else if (_mean_y_){
      pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_volume]
                                                         *pvectsz[ptsz->index_hmf]
                                                         *pvectsz[ptsz->index_completeness]
                                                         *pow(pvectsz[ptsz->index_pressure_profile],1.)
                                                         /pow(ptsz->Tcmb_gNU,1)/1.e6;
   }


   else if (_tSZ_power_spectrum_){

      int index_l = (int) pvectsz[ptsz->index_multipole];
      int flag_cov_N_cl = (int) pvectsz[ptsz->index_flag_cov_N_cl];

      double pressure_profile_at_ell = pvectsz[ptsz->index_pressure_profile];

      if (_cov_N_Cl_ && flag_cov_N_cl == _TRUE_){

      for (int index_bins_M = 0; index_bins_M <ptsz->nbins_M-1;index_bins_M++){

         pvectsz[ptsz->index_integrand_cov_N_cl_first+index_bins_M] =  pvectsz[ptsz->index_volume]
                                                                       *pvectsz[ptsz->index_hmf]
                                                                       *pvectsz[ptsz->index_completeness]
                                                                       *pow(pvectsz[ptsz->index_pressure_profile],2.)
                                                                       /pow(ptsz->Tcmb_gNU,ptsz->exponent_unit);


         pvectsz[ptsz->index_integrand_N_for_cov_N_cl_first+index_bins_M] =  pvectsz[ptsz->index_volume]
                                                                             *pvectsz[ptsz->index_hmf]
                                                                             *pvectsz[ptsz->index_completeness];



                  }

      }

      else {
         pvectsz[ptsz->index_integrand] =  pvectsz[ptsz->index_volume]
                                           *pvectsz[ptsz->index_hmf]
                                           *pvectsz[ptsz->index_dlnMdeltadlnM]
                                           *pvectsz[ptsz->index_completeness]
                                           *pow(pvectsz[ptsz->index_pressure_profile],2.)
                                           *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                           /(2*_PI_*pow(ptsz->Tcmb_gNU,ptsz->exponent_unit));






            if (_tSZ_trispectrum_) {

      int index_l = (int) pvectsz[ptsz->index_multipole];
         for (int index_l_prime = 0; index_l_prime<index_l+1;index_l_prime++){
      pvectsz[ptsz->index_multipole_for_pressure_profile] = index_l_prime;
      evaluate_pressure_profile(pvecback,pvectsz,pba,ptsz);
      double pressure_profile_at_ell_prime = pvectsz[ptsz->index_pressure_profile];

      pvectsz[ptsz->index_integrand_2h_first+index_l_prime] = pvectsz[ptsz->index_volume]
                                                              *pvectsz[ptsz->index_hmf]
                                                              *pvectsz[ptsz->index_completeness]
                                                              *pow(pressure_profile_at_ell,2.)
                                                              *pow(pressure_profile_at_ell_prime,2.)
                                                              /(4*_PI_*pow(ptsz->Tcmb_gNU,2.*ptsz->exponent_unit));
         }
   }

      if (_tSZ_2halo_){

         evaluate_halo_bias(pvecback,pvectsz,pba,ppm,pnl,ptsz);

         pvectsz[ptsz->index_integrand_2halo_term] =  pvectsz[ptsz->index_hmf]
                                                      *pvectsz[ptsz->index_dlnMdeltadlnM]
                                                      *pvectsz[ptsz->index_halo_bias]
                                                      *pvectsz[ptsz->index_completeness]
                                                      *pow(pressure_profile_at_ell,1.)
                                                      *pow(pvectsz[ptsz->index_volume]
                                                      *pvectsz[ptsz->index_pk_for_halo_bias]
                                                      *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                                      /(2*_PI_*pow(ptsz->Tcmb_gNU,ptsz->exponent_unit)),0.5);

      }

         if (_tSZ_te_y_y_){

            evaluate_temperature_mass_relation(pvecback,pvectsz,pba,ptsz);

            pvectsz[ptsz->index_integrand_te_y_y] = pvectsz[ptsz->index_te_of_m]
                                                    *pvectsz[ptsz->index_volume]
                                                    *pvectsz[ptsz->index_hmf]
                                                    *pvectsz[ptsz->index_completeness]
                                                    *pow(pvectsz[ptsz->index_pressure_profile],2.)
                                                    *ptsz->ell[index_l]*(ptsz->ell[index_l]+1.)
                                                    /(2*_PI_*pow(ptsz->Tcmb_gNU,ptsz->exponent_unit));


         }



    }

   }
   else {
      for (int i = 0; i<ptsz->number_of_integrals_per_thread; i++) pvectsz[ptsz->index_integrand + i] = 0.;
   }





   return _SUCCESS_;
}


int evaluate_temperature_mass_relation(double * pvecback,
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

      double A,B,C;
   //at z=0.0
   A = 4.763;
   B = 0.581;
   C = 0.013;

   //at z=0.5
   A = 4.353;
   B = 0.571;
   C = 0.008;

   //at z=1.0
   //A = 3.997;
   //B = 0.593;
   //C = 0.009;


   double Mfid = 3.e14 ; //Msun/h
   double M = mass*ptsz->HSEbias; //true mass

   pvectsz[ptsz->index_te_of_m] = A*pow(Eh,2./3.)*pow(M/Mfid,B+C*log(M/Mfid)); //kB*Te in keV

   }


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

   int index_l = (int) pvectsz[ptsz->index_multipole_for_pressure_profile];
   //printf("ell pp=%e\n",ptsz->ell[index_l]);

   //custom gNFW pressure profile or Battaglia et al 2012
   //if (ptsz->pressure_profile == 3 || ptsz->pressure_profile == 4 ){

      /*class_call(two_dim_ft_pressure_profile(ptsz,pba,pvectsz,&result),
                                              ptsz->error_message,
                                              ptsz->error_message);
*/
   //}

  // else {
      //printf("ell pp=%e\n",result);

      lnx_asked = log(ptsz->ell[index_l]/pvectsz[ptsz->index_l500]);

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


  // }



   pvectsz[ptsz->index_pressure_profile] = result;

   //in units of Mpc^-1*micro Kelvins
   double sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb = 283./0.5176; //1./0.5176=1.932=(5Xh+3)/2(Xh+1) with Xh = 0.76 and Pth=1.932Pe

   double characteristic_radius;
   double characteristic_multipole;
   double pressure_normalisation;


   //Battaglia et al 2012 pressure profile
   if (ptsz->pressure_profile == 4) {

      double P_200; //in units of eV/cm^3, corresponds to Pth

      double rho_crit_at_z = pvectsz[ptsz->index_Rho_crit]; //in units of h^2 M_sun/Mpc^3
      double R_200crit = pvectsz[ptsz->index_r200c]; //in units of h^-1 Mpc
      double f_b = pba->Omega0_b/ptsz->Omega_m_0;
      double _G_in_eV_Mpc_over_Msun2 = _G_/(_eV_ *_Mpc_over_m_ /_M_sun_/_M_sun_);
      double Eh = pvecback[pba->index_bg_H]/ptsz->H0_in_class_units;


      //double P_200_boris = _G_in_eV_Mpc_over_Msun2*pvectsz[ptsz->index_m200c]
      //            /pba->h*200.*rho_crit_at_z*pba->h*pba->h
      //            *f_b/2./(R_200crit/pba->h)/pow(_Mpc_over_m_,3.)/1e6;

      P_200 = pvectsz[ptsz->index_m200c]/(R_200crit)*f_b
              *2.61051e-18*pow(100.*pba->h*Eh,2.);




      pressure_normalisation = P_200;///1.932; //1./0.5176=1.932=(5Xh+3)/2(Xh+1) with Xh = 0.76
      //divide by 1.932 to obtain Pe
      characteristic_radius = pvectsz[ptsz->index_r200c]/pba->h; //in Mpc
      characteristic_multipole = pvectsz[ptsz->index_l200c];





   }

   else {


      double C_pressure = 1.65  // formula D1 of WMAP 7 year paper (normalisation of pressure profile)
                                    *pow(pba->h/0.7,2)
                                    *pow(pvecback[pba->index_bg_H]/pba->H0,8./3.)
                                    *pow(pvectsz[ptsz->index_m500]/(3.e14*0.7),2./3.+ptsz->alpha_p)
                                    *pow(pvectsz[ptsz->index_m500]*ptsz->HSEbias/3.e14, ptsz->delta_alpha);

      pressure_normalisation = C_pressure
                               *ptsz->P0GNFW
                               *pow(0.7/pba->h, 1.5); // formula D3 of WMAP 7 year paper (normalisation of pressure profile)
                               //last term should not be there when computing P13 pressure profile? need to ask planck people
                               //this factor actually is in komatsu's integrand.f90...

      characteristic_radius = pvectsz[ptsz->index_r500]/pba->h; // in Mpc
      characteristic_multipole = pvectsz[ptsz->index_l500];


   }

   pvectsz[ptsz->index_pressure_profile] = -0.953652      //gnu at 150 GHz, should be -0.9533281807274405
                                           *sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb
                                           *pressure_normalisation
                                           *pvectsz[ptsz->index_pressure_profile]
                                           *(4*_PI_)
                                           *pow(characteristic_multipole,-2)
                                           *characteristic_radius //rs in Mpc
                                           /50.; //pressure normalised to 50eV/cm^3






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
    double redshift = pvectsz[ptsz->index_z];

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
    double sn_cutoff = 6.;
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

    int index_patches;
    double sum_skyfracs = 0.;

    for (index_patches =0;
    index_patches<ptsz->nskyfracs;
    index_patches++){

    double y1 = ptsz->ylims[index_patches][l1];
    double y2 = ptsz->ylims[index_patches][l2];
    y = y1 + (y2-y1)/(th2-th1)*(thp-th1);

    double c2 = erf_compl_ps(yp,y,sn_cutoff);

    comp_at_M_and_z += c2*ptsz->skyfracs[index_patches];
    sum_skyfracs += ptsz->skyfracs[index_patches];
    }
    //Now divide by total sky fraction
    comp_at_M_and_z = comp_at_M_and_z/sum_skyfracs;
    }
   if (ptsz->which_ps_sz == 0)
      pvectsz[ptsz->index_completeness] = 1.;
   else if (ptsz->which_ps_sz == 1) //ps resolved
      pvectsz[ptsz->index_completeness] = comp_at_M_and_z;
   else if (ptsz->which_ps_sz == 2) //ps unresolved
      pvectsz[ptsz->index_completeness] = (1.-comp_at_M_and_z);
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
   double y = log10(200.);
   double ATT = 1.+0.24*y*exp(-pow(4./y,4.));
   double aTTT = 0.44*y-0.88;
   double BTT = 0.183;
   double bTTT = 1.5;
   double CTT = 0.019+0.107*y+0.19*exp(-pow(4./y,4.));
   double cTTT = 2.4;

   pvectsz[ptsz->index_halo_bias] = 1.-ATT*(pow(nuTink,aTTT)/(pow(nuTink,aTTT)+pow(ptsz->delta_cSZ,aTTT)))
                                                      +BTT*pow(nuTink,bTTT)+CTT*pow(nuTink,cTTT);

   int index_l = (int)  pvectsz[ptsz->index_multipole];
   double z = pvectsz[ptsz->index_z];
   double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z); //multiply by h to get in Mpc/h

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
                                       z,
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

   //beginning of mass function
   pvectsz[ptsz->index_rVIR] =
   pow(3.*exp(logM)/(4*_PI_*pvectsz[ptsz->index_Delta_c]*pvectsz[ptsz->index_Rho_crit]),1./3.);

   //For SC14 C-M relation
   //we need r200crit in place of rVIR:
   if (ptsz->concentration_parameter==3)
      pvectsz[ptsz->index_rVIR] =
      pow(3.*exp(logM)/(4*_PI_*200.*pvectsz[ptsz->index_Rho_crit]),1./3.);



   //D08 c-m relation
   if (ptsz->concentration_parameter==0){
      pvectsz[ptsz->index_cVIR] =7.85*pow(exp(logM)/2.e12,-0.081)*pow(1.+z,-0.71);

      //pvectsz[ptsz->index_cVIR] =5.72*pow(exp(logM)/1.e14,-0.081)*pow(1.+z,-0.71);

   }

   //S00 c-m relation
   else if (ptsz->concentration_parameter==1){
      pvectsz[ptsz->index_cVIR] =10.*pow(exp(logM)/3.42e12,-0.2)/(1.+z);
   }

   //K10 c-m relation
   else if (ptsz->concentration_parameter==2){
      class_call(CvirMvirKLYPIN(&pvectsz[ptsz->index_cVIR],
                                             logM,
                                             z,
                                             ptsz),
                      ptsz->error_message,
                      ptsz->error_message);
   }

   //SC14 c-m relation
   else if (ptsz->concentration_parameter==3){
      class_call(C200M200SC14(&pvectsz[ptsz->index_cVIR],
                                          logM,
                                          z,
                                          ptsz),
                      ptsz->error_message,
                      ptsz->error_message);
   }


   //Z09 interpolated c-m relation
   else if (ptsz->concentration_parameter==4){
      class_call(CvirMvirZHAO(&pvectsz[ptsz->index_cVIR],
                                          logM,
                                          log(z),
                                          ptsz),
                      ptsz->error_message,
                      ptsz->error_message);
   }


   //Scale radius:
   //(note that rs is actually r200c/c200 for SC14 concentration-mass relation)
   pvectsz[ptsz->index_rs] = pvectsz[ptsz->index_rVIR]/pvectsz[ptsz->index_cVIR];


   //ell_s is *not* used
   pvectsz[ptsz->index_ls] = pvecback[pba->index_bg_ang_distance]*pba->h/pvectsz[ptsz->index_rs];


   //Convert m to m500
   //for the pressure profile
   //(except when HMF is at m500c, i.e., T08@M500 and B16M500c)
   if (ptsz->MF!=5 && ptsz->MF!=7){

      class_call(m_to_mDEL(exp(logM),
                                     pvectsz[ptsz->index_rs],
                                     pvectsz[ptsz->index_cVIR],
                                     500.*(pvectsz[ptsz->index_Rho_crit]),
                                     &pvectsz[ptsz->index_m500],
                                     ptsz),
                      ptsz->error_message,
                      ptsz->error_message);


      class_call(m_to_mDEL(exp(logM),
                                     pvectsz[ptsz->index_rs],
                                     pvectsz[ptsz->index_cVIR],
                                     200.*(pvectsz[ptsz->index_Rho_crit]),
                                     &pvectsz[ptsz->index_m200c],
                                     ptsz),
                      ptsz->error_message,
                      ptsz->error_message);

   }

   //if HMF=T08@m500c -> m200=m500=M
   else pvectsz[ptsz->index_m500] = exp(logM);

   if (ptsz->mass_dependent_bias == 1)
      ptsz->HSEbias = 1./(0.8/(1.+ ptsz->Ap*pow(pvectsz[ptsz->index_m500]/3.e14,ptsz->alpha_b)));



   pvectsz[ptsz->index_m500] = pvectsz[ptsz->index_m500]/ptsz->HSEbias;


   pvectsz[ptsz->index_r500] = pow(3.*pvectsz[ptsz->index_m500]/(4.*_PI_*500.*pvectsz[ptsz->index_Rho_crit]),1./3.);
   //in units of h^-1 Mpc


   pvectsz[ptsz->index_r200c] = pow(3.*pvectsz[ptsz->index_m200c]/(4.*_PI_*200.*pvectsz[ptsz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc


   pvectsz[ptsz->index_l500] = pvecback[pba->index_bg_ang_distance]*pba->h/pvectsz[ptsz->index_r500];

   pvectsz[ptsz->index_l200c] = pvecback[pba->index_bg_ang_distance]*pba->h/pvectsz[ptsz->index_r200c];

   //m200-mean or m1600-mean for HMF:
   //m180-mean for Jenkins 2001 HMF.
   //Bypassed for T08@m500 B16@m500.

   if (ptsz->MF!=5 && ptsz->MF!=7){

      //Jenkins et al 2001 @ m180-mean
      if (ptsz->MF==3){
         class_call(m_to_mDEL(exp(logM),
                                        pvectsz[ptsz->index_rs],
                                        pvectsz[ptsz->index_cVIR],
                                        180.*( pvecback[pba->index_bg_Omega_m])
                                        *pvectsz[ptsz->index_Rho_crit],
                                        &pvectsz[ptsz->index_m200],
                                        ptsz),
                         ptsz->error_message,
                         ptsz->error_message);
      }


      //Tinker et al 2008 @ m1600-mean
      else if (ptsz->MF==6){
         class_call(m_to_mDEL(exp(logM),
                                        pvectsz[ptsz->index_rs],
                                        pvectsz[ptsz->index_cVIR],
                                        1600.*( pvecback[pba->index_bg_Omega_m])
                                        *pvectsz[ptsz->index_Rho_crit],
                                        &pvectsz[ptsz->index_m200],
                                        ptsz),
                         ptsz->error_message,
                         ptsz->error_message);
      }

      //Tinker et al 2008 @ m200-mean
      //Tinker et al 2010 @ m200-mean
      //Boquet et al 2015 @ m200-mean

      else{
         class_call(m_to_mDEL(exp(logM),
                                        pvectsz[ptsz->index_rs],
                                        pvectsz[ptsz->index_cVIR],
                                        200.*pvecback[pba->index_bg_Omega_m]
                                        *pvectsz[ptsz->index_Rho_crit],
                                        &pvectsz[ptsz->index_m200],
                                        ptsz),
                         ptsz->error_message,
                         ptsz->error_message);

   pvectsz[ptsz->index_dlnMdeltadlnM]= evaluate_dlnMdeltadlnM(logM,
                                                    pvecback,
                                                    pvectsz,
                                                    pba,
                                                    pnl,
                                                    ptsz);

   }
   }

   else pvectsz[ptsz->index_m200] = exp(logM);




   //No-pres h^-1 Mpc
   pvectsz[ptsz->index_Rh] =
   pow(3.*pvectsz[ptsz->index_m200]/
         (4*_PI_*ptsz->Omega_m_0
          *ptsz->Rho_crit_0),1./3.);

   if (ptsz->HMF_prescription_NCDM == 0) //Matter
      pvectsz[ptsz->index_Rh] =
      pow(3.*pvectsz[ptsz->index_m200]/
            (4*_PI_*(pba->Omega0_cdm+pba->Omega0_b)
             *ptsz->Rho_crit_0),1./3.);

   else if (ptsz->HMF_prescription_NCDM == 1) //CDM
      pvectsz[ptsz->index_Rh] = pow(3.*pvectsz[ptsz->index_m200]/
                                                      (4*_PI_*(pba->Omega0_cdm+pba->Omega0_b)
                                                       *ptsz->Rho_crit_0),1./3.);




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

   //pvectsz[ptsz->index_dndlogRh] = 3.;//pvectsz[ptsz->index_mf];


   //printf("mf = %e\t%e\t%e\n",pvectsz[ptsz->index_Rh],pvectsz[ptsz->index_dlognudlogRh],pvectsz[ptsz->index_mf]);
   ///end of mass function evaluations



return _SUCCESS_;
}


int write_output_to_files_ell_indep_ints(struct nonlinear * pnl,
                                                             struct background * pba,
                                                             struct tszspectrum * ptsz){
   //This block has to be commented
   //when using the Python wrapper
   //or alternatively set sz_verbose=0
   int index_l;

    char Filepath[_ARGUMENT_LENGTH_MAX_];

   if (ptsz->sz_verbose > 0)
   {
      FILE *fp;


if (ptsz->has_mean_y){
         sprintf(Filepath,
                     "%s%s%s",
                     ptsz->root,
                     "mean_y",
                     ".txt");


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
      //fclose(fp);

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

      fprintf(fp,"#N_tot\n");
      fprintf(fp,"%e\n",ptsz->hmf_int);
      printf("->Output written in %s\n",Filepath);

      fclose(fp);

    }



   }

return _SUCCESS_;
}

int write_output_to_files_cl(struct nonlinear * pnl,
                             struct background * pba,
                             struct tszspectrum * ptsz){


   int index_l;
   char Filepath[_ARGUMENT_LENGTH_MAX_];

   if (ptsz->sz_verbose > 0)
   {
      FILE *fp;

      if (ptsz->has_sz_ps){

      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "szpowerspectrum",
                  ".txt");


      //printf("Writing output files in %s\n",Filepath);

      fp=fopen(Filepath, "w");


      fprintf(fp,"#Input mass bias b = %e\n",
                  1.-1./ptsz->HSEbias);
      fprintf(fp,"#sigma8 = %e\n",
                  pnl->sigma8[pnl->index_pk_m]);
      fprintf(fp,"#Omega_m = %e\n",
                  ptsz->Omega_m_0);
      fprintf(fp,"#h = %e\n",
                  pba->h);
      fprintf(fp,"#temperature_mass_relation = %d\n",
                              ptsz->temperature_mass_relation);


         fprintf(fp,"#(3-5 need to be divided by f_sky, 6 by sqrt(f_sky)\n");
         if (ptsz->exponent_unit == 2.)
            fprintf(fp,"#1:l\t 2:10^12*D_l^tSZ\t 3:sigma_g_C_l^2\t 4:T_ll \t 5:sigma_g_C_l^2_binned \t 6:sig_D_l_binned\t 7:D_l_2h\t 8:Te \n");
         if (ptsz->exponent_unit == 0.)
            fprintf(fp,"#1:l\t 2:D_l^tSZ [muK^2]\t 3:sigma_g_C_l^2\t 4:T_ll \t 5:sigma_g_C_l^2_binned \t 6:sig_D_l_binned\t 7:D_l_2h\t 8:Te \n");


      for (index_l=0;index_l<ptsz->nlSZ;index_l++){

         double sig_cl_squared;
         double ell = ptsz->ell[index_l];
         sig_cl_squared = 2.*pow(ptsz->cl_sz[index_l]/(ell*(ell+1.))*2.*_PI_,2.)/(2.*ell+1.);


         //binned gaussian variance
         double sig_cl_squared_binned;
         double ln_ell_min;
         double ln_ell_max;
         double ln_ell_down;
         double ln_ell_up;
         double n_modes;
         if (index_l == 0){
            ln_ell_up = log(ptsz->ell[index_l+1]);
            ln_ell_max = log(ell) + 0.5*(ln_ell_up-log(ell));
            ln_ell_min = log(ell) - 0.5*(ln_ell_up-log(ell));
            n_modes = exp(ln_ell_max)-exp(ln_ell_min);
            sig_cl_squared_binned = sig_cl_squared/n_modes;
         }
         else if (index_l == ptsz->nlSZ -1){
            ln_ell_down = log(ptsz->ell[index_l-1]);
            ln_ell_min = log(ell) - 0.5*(log(ell)-ln_ell_down);
            ln_ell_max = log(ell) + 0.5*(log(ell)-ln_ell_down);
            n_modes = exp(ln_ell_max)-exp(ln_ell_min);
            sig_cl_squared_binned = sig_cl_squared/n_modes;
         }
         else {
            ln_ell_down = log(ptsz->ell[index_l-1]);
            ln_ell_up = log(ptsz->ell[index_l+1]);
            ln_ell_min = log(ell) - 0.5*(log(ell)-ln_ell_down);
            ln_ell_max = log(ell) + 0.5*(ln_ell_up-log(ell));
            n_modes = exp(ln_ell_max)-exp(ln_ell_min);
            sig_cl_squared_binned = sig_cl_squared/n_modes;
         }

         //normalised cov:

         ptsz->cov_cl_cl[index_l] = sig_cl_squared_binned +  ptsz->tllprime_sz[index_l][index_l];

         int index_M_bins;
         for (index_M_bins=0;index_M_bins<ptsz->nbins_M-1;index_M_bins++){
            ptsz->r_N_cl[index_l][index_M_bins] = ptsz->r_N_cl[index_l][index_M_bins]/sqrt(ptsz->cov_cl_cl[index_l]);

         }

            fprintf(fp,
                    "%e\t\t %e\t\t %e\t\t %e\t\t %e\t\t%e\t\t%e\t\t%e\n",
                    ptsz->ell[index_l],
                    ptsz->cl_sz[index_l],
                    sig_cl_squared,
                    ptsz->tllprime_sz[index_l][index_l],
                    sig_cl_squared_binned,
                    ell*(ell+1.)/(2.*_PI_)
                    *sqrt(sig_cl_squared_binned+ptsz->tllprime_sz[index_l][index_l]),
                    ptsz->cl_sz_2h[index_l],
                    ptsz->cl_te_y_y[index_l]/ptsz->cl_sz[index_l]
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
           ptsz->r_cl_clp[index_l][index_l_prime] = ptsz->tllprime_sz[index_l][index_l_prime]
                                                                        /sqrt(ptsz->cov_cl_cl[index_l])
                                                                        /sqrt(ptsz->cov_cl_cl[index_l_prime]);

           ptsz->r_cl_clp[index_l_prime][index_l] = ptsz->r_cl_clp[index_l][index_l_prime];

           ell_prime = ptsz->ell[index_l_prime];
           ell = ptsz->ell[index_l];
           ptsz->trispectrum_ref[index_l][index_l_prime] = ell*(ell+1.)/(2.*_PI_)*ell_prime*(ell_prime+1.)/(2.*_PI_)*ptsz->tllprime_sz[index_l][index_l_prime];

           ptsz->trispectrum_ref[index_l_prime][index_l] = ptsz->trispectrum_ref[index_l][index_l_prime];



         };



      //FILE *fp;

if (ptsz->has_sz_cov_N_Cl){
      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "szpowerspectrum_cov_N_cl",
                  ".txt");

      fp=fopen(Filepath, "w");

      int index_M_bins;
      for (index_l=0;index_l<ptsz->nlSZ;index_l++){
         for (index_M_bins=0;index_M_bins<ptsz->nbins_M-1;index_M_bins++){
            fprintf(fp,"%e\t",ptsz->cov_N_cl[index_l][index_M_bins]);
      }
         fprintf(fp,"\n");
      }
      fclose(fp);

      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "szpowerspectrum_r_N_cl",
                  ".txt");

      fp=fopen(Filepath, "w");


      for (index_l=0;index_l<ptsz->nlSZ;index_l++){
         for (index_M_bins=0;index_M_bins<ptsz->nbins_M-1;index_M_bins++){
            fprintf(fp,"%e\t",ptsz->r_N_cl[index_l][index_M_bins]);
         }
         fprintf(fp,"\n");
      }
      fclose(fp);

    }



 if (ptsz->has_sz_trispec){
      sprintf(Filepath,
                  "%s%s%s",
                  ptsz->root,
                  "szpowerspectrum_r_cl_clp",
                  ".txt");

      fp=fopen(Filepath, "w");

      int index_l_prime;
      for (index_l=0;index_l<ptsz->nlSZ;index_l++){
       for (index_l_prime=0;index_l_prime<ptsz->nlSZ;index_l_prime++) {
            fprintf(fp,"%e\t",ptsz->r_cl_clp[index_l][index_l_prime]);
         }
         fprintf(fp,"\n");
      }
      fclose(fp);
}




   }

   return _SUCCESS_;
}




int show_preamble_messages(struct background * pba,
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

      ptsz->Omega_m_0 = pvecback[pba->index_bg_Omega_m];
      ptsz->Omega_ncdm_0 = ptsz->Omega_m_0
      -pba->Omega0_b
      -pba->Omega0_cdm;

      ptsz->H0_in_class_units = pvecback[pba->index_bg_H];

      if (pba->Omega0_lambda != 0.) OmegaM = 1-pba->Omega0_lambda;
      else OmegaM = 1-pba->Omega0_fld;


      ptsz->Sigma8OmegaM_SZ = pnl->sigma8[pnl->index_pk_m]*pow(OmegaM/0.28,3./8.);

      if (ptsz->sz_verbose > 0)
      {
         printf("Computing SZ angular power spectrum\n");

         //Cosmological parameters
         //printf("pba->H0 = %e, %e\n",pba->H0,ptsz->H0_in_class_units);
         printf("pba->h = %e\n",pba->h);
         printf("Omega0_b = %e\n",pba->Omega0_b);
         printf("Omega0_cdm = %e\n",pba->Omega0_cdm);
         printf("Omega0_m = %e\n",ptsz->Omega_m_0 );
         printf("Omega0_ncdm = %e\n",
                   ptsz->Omega_m_0
                   -pba->Omega0_b
                   -pba->Omega0_cdm);
         printf("Rho_crit_0 = %e in units of h^2 M_sun/Mpc^3\n",ptsz->Rho_crit_0 );


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
         if (ptsz->MF!=5 && ptsz->MF!=7)
         {
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
         }

         printf("->sigma8(OmegaM/0.28)^3/8 = %e\n",ptsz->Sigma8OmegaM_SZ);
         printf("->sigma8*(OmegaM/B)^3/8*h^-1/5 = %e\n",
                   pnl->sigma8[pnl->index_pk_m]*pow(OmegaM/ptsz->HSEbias,3./8.)*pow(pba->h,-1./5.));
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
   ptsz->Tcmb_gNU =
   pba->T_cmb
   *((_h_P_*frequency_in_Hz
       /(_k_B_*pba->T_cmb))
      *(1./tanh((_h_P_*frequency_in_Hz
                      /(_k_B_*pba->T_cmb))
                     /2.))
      -4.);
   if (ptsz->sz_verbose > 0){
      printf("->Tcmb_gNU at 150GHz= %e\n",ptsz->Tcmb_gNU);
      printf("->gNU at 150GHz= %e\n",ptsz->Tcmb_gNU/pba->T_cmb);
      printf("->Tcmb = %e K\n",pba->T_cmb);

   }

return _SUCCESS_;
}




int show_results(struct background * pba,
                         struct nonlinear * pnl,
                         struct primordial * ppm,
                         struct tszspectrum * ptsz){

  printf("#1:ell\t\t\t 2:y^2 (tSZ)\n");
   int index_l;
   for (index_l=0;index_l<ptsz->nlSZ;index_l++){

   //if (ptsz->ell[index_l]==0) ptsz->y_monopole = cl/pow(ptsz->Tcmb_gNU,1)/1.e6;
   //divide by gNU at 150GHz cause it was in Komatsu's formula
   //divide by 1e6


   printf("%e\t\t %e \n",ptsz->ell[index_l],ptsz->cl_sz[index_l]);



 if (ptsz->has_sz_trispec){

    printf("\n");
      int index_l_prime;
      for (index_l_prime=0;index_l_prime<index_l+1;index_l_prime++)
      {


            printf("%e\t\t %e \n",
                      ptsz->ell[index_l_prime],
                      ptsz->tllprime_sz[index_l][index_l_prime]);
      }

printf("\n");

}
 if (ptsz->has_sz_cov_N_Cl){
      int index_M_bins;
      for (index_M_bins=0;index_M_bins<ptsz->nbins_M-1;index_M_bins++)
      {


         printf("%e\t\t %e \n",
                   ptsz->M_bins[index_M_bins],
                   ptsz->r_N_cl[index_l][index_M_bins]);
      }

        printf("\n");
}



   }

 if (ptsz->has_mean_y)
   printf("mean_y =  %e \n",ptsz->y_monopole);
if (ptsz->has_hmf)
   printf("N_tot =  %e \n",ptsz->hmf_int);

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
   }


   free(ptsz->ell_trispectrum);
   free(ptsz->ell_plc);
   free(ptsz->ell_plc_low);
   free(ptsz->ell_mock);

   return _SUCCESS_;

}



int initialise_and_allocate_memory(struct tszspectrum * ptsz){


   class_alloc(ptsz->ln_x_for_pp, ptsz->ln_x_size_for_pp*sizeof(double),ptsz->error_message);
   int i;
   for (i=0;i<ptsz->ln_x_size_for_pp;i++){
      ptsz->ln_x_for_pp[i] = log(ptsz->x_inSZ)+i*(log(ptsz->x_outSZ)-log(ptsz->x_inSZ))/(ptsz->ln_x_size_for_pp-1.);
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
   ptsz->nbins_M = 50;
   //printf("mass bins for N-C_l cov, nbins = %d\n", ptsz->nbins_M);
   class_alloc(ptsz->M_bins,
                        ptsz->nbins_M*sizeof(double),
                        ptsz->error_message);

   for (int i=0;i<ptsz->nbins_M;i++){
      ptsz->M_bins[i] = exp(log(ptsz->M1SZ)+i*(log(ptsz->M2SZ)-log(ptsz->M1SZ))/(ptsz->nbins_M-1.));
      //printf("M = %e \t i=%d\n",ptsz->M_bins[i],i);
   }


   //function of redshift
   ptsz->index_dlnMdeltadlnM = 0;
   ptsz->index_te_of_m = ptsz->index_dlnMdeltadlnM + 1;
   ptsz->index_flag_cov_N_cl = ptsz->index_te_of_m + 1;
   ptsz->index_multipole_for_pressure_profile = ptsz->index_flag_cov_N_cl + 1;
   ptsz->index_multipole_prime = ptsz->index_multipole_for_pressure_profile +1;
   ptsz->index_md = ptsz->index_multipole_prime +1;
   ptsz->index_Rho_crit = ptsz->index_md +1;
   ptsz->index_Delta_c  = ptsz->index_Rho_crit +1;
   ptsz->index_rVIR  = ptsz->index_Delta_c +1;
   ptsz->index_cVIR  = ptsz->index_rVIR +1;
   ptsz->index_m500  = ptsz->index_cVIR +1;
   ptsz->index_r500  = ptsz->index_m500 +1;
   ptsz->index_l500  = ptsz->index_r500 +1;
   ptsz->index_m200  = ptsz->index_l500 +1;
   ptsz->index_Rh  = ptsz->index_m200 +1;
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
   ptsz->index_completeness = ptsz->index_pressure_profile +1;
   ptsz->index_volume = ptsz->index_completeness + 1;
   ptsz->index_hmf = ptsz->index_volume + 1;
   ptsz->index_halo_bias = ptsz->index_hmf + 1;
   ptsz->index_k_value_for_halo_bias = ptsz->index_halo_bias +1;
   ptsz->index_pk_for_halo_bias = ptsz->index_k_value_for_halo_bias +1;



   //quantities integrated over redshift
   ptsz->index_integral = ptsz->index_pk_for_halo_bias+1;
   ptsz->index_integral_te_y_y = ptsz->index_integral + 1;
   ptsz->index_integral_2halo_term = ptsz->index_integral_te_y_y + 1;
   ptsz->index_integral_2h_first =  ptsz->index_integral_2halo_term +1;
   ptsz->index_integral_2h_last = ptsz->index_integral_2h_first + ptsz->nlSZ - 1;
   ptsz->index_integral_cov_N_cl_first =  ptsz->index_integral_2h_last +1;
   ptsz->index_integral_cov_N_cl_last =  ptsz->index_integral_cov_N_cl_first + ptsz->nbins_M - 2;

   ptsz->index_integral_N_for_cov_N_cl_first =  ptsz->index_integral_cov_N_cl_last +1;
   ptsz->index_integral_N_for_cov_N_cl_last =  ptsz->index_integral_N_for_cov_N_cl_first + ptsz->nbins_M - 2;


   //quantities integrated over m
   ptsz->index_integral_over_m = ptsz->index_integral_N_for_cov_N_cl_last+1;
   ptsz->index_integral_te_y_y_over_m = ptsz->index_integral_over_m + 1;
   ptsz->index_integral_2halo_term_over_m = ptsz->index_integral_te_y_y_over_m + 1;
   ptsz->index_integral_2h_first_over_m =  ptsz->index_integral_2halo_term_over_m +1;
   ptsz->index_integral_2h_last_over_m = ptsz->index_integral_2h_first_over_m + ptsz->nlSZ - 1;
   ptsz->index_integral_cov_N_cl_first_over_m =  ptsz->index_integral_2h_last_over_m +1;
   ptsz->index_integral_cov_N_cl_last_over_m =  ptsz->index_integral_cov_N_cl_first_over_m + ptsz->nbins_M - 2;

   ptsz->index_integral_N_for_cov_N_cl_first_over_m =  ptsz->index_integral_cov_N_cl_last_over_m +1;
   ptsz->index_integral_N_for_cov_N_cl_last_over_m =  ptsz->index_integral_N_for_cov_N_cl_first_over_m + ptsz->nbins_M - 2;





   //integrands at m and z
   ptsz->index_integrand =  ptsz->index_integral_N_for_cov_N_cl_last_over_m + 1;
   ptsz->index_integrand_te_y_y = ptsz->index_integrand + 1;

   ptsz->index_integrand_2halo_term = ptsz->index_integrand_te_y_y +1;

   ptsz->index_integrand_2h_first = ptsz->index_integrand_2halo_term + 1;
   ptsz->index_integrand_2h_last = ptsz->index_integrand_2h_first + ptsz->nlSZ - 1;

   ptsz->index_integrand_cov_N_cl_first = ptsz->index_integrand_2h_last + 1;
   ptsz->index_integrand_cov_N_cl_last = ptsz->index_integrand_cov_N_cl_first + ptsz->nbins_M - 2;

   ptsz->index_integrand_N_for_cov_N_cl_first = ptsz->index_integrand_cov_N_cl_last + 1;
   ptsz->index_integrand_N_for_cov_N_cl_last =  ptsz->index_integrand_N_for_cov_N_cl_first + ptsz->nbins_M - 2;



   //final size of pvecsz vector
   ptsz->tsz_size  = ptsz->index_integrand_N_for_cov_N_cl_last + 1;



   ptsz->index_integrands_first = ptsz->index_integrand;
   ptsz->index_integrands_last = ptsz->index_integrand_N_for_cov_N_cl_last;


   ptsz->index_integrals_over_m_first = ptsz->index_integral_over_m;
   ptsz->index_integrals_over_m_last = ptsz->index_integral_N_for_cov_N_cl_last_over_m;

   ptsz->index_integrals_over_z_first = ptsz->index_integral;
   ptsz->index_integrals_over_z_last = ptsz->index_integral_N_for_cov_N_cl_last;


   ptsz->number_of_integrals_per_thread = ptsz->index_integrals_over_z_last-ptsz->index_integrals_over_z_first+1;



   ptsz->hmf_int = 0.;
   ptsz->y_monopole = 0.;


   class_alloc(ptsz->cl_sz,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_te_y_y,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cov_cl_cl,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);
   class_alloc(ptsz->cl_sz_2h,sizeof(double *)*ptsz->nlSZ,ptsz->error_message);

   class_alloc(ptsz->tllprime_sz,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
   class_alloc(ptsz->trispectrum_ref,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
   class_alloc(ptsz->r_cl_clp,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
   class_alloc(ptsz->cov_N_cl,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
   class_alloc(ptsz->r_N_cl,ptsz->nlSZ*sizeof(double *),ptsz->error_message);
   int index_l,index_l_prime;
   for (index_l=0;index_l<ptsz->nlSZ;index_l++){
      ptsz->cl_sz[index_l] = 0.;
      ptsz->cl_te_y_y[index_l] = 0.;
      ptsz->cl_sz_2h[index_l] = 0.;
      ptsz->cov_cl_cl[index_l] = 0.;

      class_alloc(ptsz->tllprime_sz[index_l],(index_l+1)*sizeof(double),ptsz->error_message);
      class_alloc(ptsz->trispectrum_ref[index_l],ptsz->nlSZ*sizeof(double),ptsz->error_message);
      class_alloc(ptsz->r_cl_clp[index_l],ptsz->nlSZ*sizeof(double),ptsz->error_message);
      class_alloc(ptsz->cov_N_cl[index_l],(ptsz->nbins_M-1)*sizeof(double),ptsz->error_message);
      class_alloc(ptsz->r_N_cl[index_l],(ptsz->nbins_M-1)*sizeof(double),ptsz->error_message);

      for (index_l_prime = 0; index_l_prime<index_l+1; index_l_prime ++){
         ptsz->tllprime_sz[index_l][index_l_prime] = 0.;
      }

      for (index_l_prime = 0; index_l_prime<ptsz->nlSZ; index_l_prime ++){
        ptsz->r_cl_clp[index_l][index_l_prime] = 0.;
        ptsz->trispectrum_ref[index_l][index_l_prime] = 0.;
          }

      int index_M_bins;
      for (index_M_bins = 0; index_M_bins<ptsz->nbins_M-1; index_M_bins ++){
         ptsz->cov_N_cl[index_l][index_M_bins] = 0.;
         ptsz->r_N_cl[index_l][index_M_bins] = 0.;
      }
   }





   ptsz->index_integrand_id_hmf = 0;
   ptsz->index_integrand_id_mean_y = ptsz->index_integrand_id_hmf + 1;
   ptsz->index_integrand_id_sz_ps_first = ptsz->index_integrand_id_mean_y + 1;
   ptsz->index_integrand_id_sz_ps_last = ptsz->index_integrand_id_sz_ps_first + ptsz->nlSZ - 1;
   //ptsz->index_integrand_id_sz_trispec_first = ptsz->index_integrand_id_sz_ps_last + 1;
   //ptsz->index_integrand_id_sz_trispec_last = ptsz->index_integrand_id_sz_trispec_first + ptsz->nlSZ*(ptsz->nlSZ+1)/2 - 1;
   //ptsz->number_of_integrands = ptsz->index_integrand_id_sz_trispec_last + 1;
   ptsz->number_of_integrands = ptsz->index_integrand_id_sz_ps_last + 1;
   //ptsz->number_of_integrands = 1;

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
double mvir;
double mvir1,mvir2,rvir1,rvir2;
double cvir1,cvir2,rs1,rs2,dlnM200ddlnm,m200d2,m200d1;

delc = pvectsz[ptsz->index_Delta_c];
rhoc = pvectsz[ptsz->index_Rho_crit];
omega = pvecback[pba->index_bg_Omega_m];

double z = pvectsz[ptsz->index_z];

  mvir1=exp(logM-tol);
  mvir2=exp(logM+tol);
  rvir1=pow(3.*mvir1/4./_PI_/delc/rhoc,1./3.); //! virial radius, h^-1 Mpc
  rvir2=pow(3.*mvir2/4./_PI_/delc/rhoc,1./3.); //! virial radius, h^-1 Mpc

  //! JCH edit: from personal communication with Battaglia: in their paper,
  //! they changed Duffy's M_pivot from 2d12 Msun/h to 1d14 Msun/h, so
  //! normalization becomes 5.72 instead of 7.85
  //! N.B.: this is the same as Eq. (D17) of Komatsu et al., arXiv:1001.4538

  cvir1=5.72*pow(mvir1/1e14,-0.081)/pow(1.+z,0.71);
  cvir2=5.72*pow(mvir2/1e14,-0.081)/pow(1.+z,0.71);
  rs1=rvir1/cvir1; //! NFW scale radius, h^-1 Mpc
  rs2=rvir2/cvir2; //! NFW scale radius, h^-1 Mpc

  class_call(m_to_mDEL(mvir1,
                       rs1,
                       cvir1,
                       200.*omega*rhoc,
                       &m200d1,
                       ptsz),
                  ptsz->error_message,
                  ptsz->error_message);
  class_call(m_to_mDEL(mvir2,
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
