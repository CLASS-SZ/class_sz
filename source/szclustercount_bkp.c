/** @file szpowerspectrum.c Documented SZ module.
 *
 * Boris Bolliet and Florian Ruppin, 11.2017
 *
 *This module is dedicated to the computation of
 *the number counts from Halo Mass Functions (HMF)
 */

#include "szclustercount.h"

int szcount_init(struct background * pba,
                 struct nonlinear * pnl,
                 struct primordial * ppm,
                 struct tszspectrum * ptsz,
                 struct szcount * pcsz)
{
  pcsz->sz_verbose = ptsz->sz_verbose;
  if (pcsz->has_sz_counts == _FALSE_)
  {
    if (pcsz->sz_verbose > 0)
      printf("->No SZ cluster counts requested. SZ cluster count module skipped.\n");
  }

  else
  {

    if (pcsz->sz_verbose > 0)
      printf("->Computing SZ cluster counts.\n");

   // // if ((ptsz->experiment == 0 && ptsz->has_completeness_for_ps_SZ == 1)
   // //  || (ptsz->experiment == 0 && ptsz->has_sz_counts  == 1))
   //    read_Planck_noise_map(ptsz);


    double * Pvecback;
    double * Pvectsz;

    class_alloc(Pvectsz,
                ptsz->tsz_size*sizeof(double),
                pcsz->error_message);

    class_alloc(Pvecback,
                pba->bg_size*sizeof(double),
                pcsz->error_message);
    initialise_and_allocate_memory_cc(ptsz,pcsz);

    class_call(compute_count_sz(pba,
                                pnl,
                                ppm,
                                ptsz,
                                pcsz,
                                Pvecback,
                                Pvectsz),
               pcsz->error_message,
               pcsz->error_message);


    free(Pvecback);
    free(Pvectsz);

  }

  return _SUCCESS_;

}


int szcount_free(struct szcount * pcsz)
{
  if (pcsz->has_sz_counts == _TRUE_){
  free(pcsz->redshift);
  free(pcsz->dNdzdy_theoretical);
  free(pcsz->z_center);
  free(pcsz->steps_m);
  free(pcsz->steps_z);
  free(pcsz->logy);
  }

  return _SUCCESS_;
}


int compute_count_sz(struct background * pba,
                     struct nonlinear * pnl,
                     struct primordial * ppm,
                     struct tszspectrum * ptsz,
                     struct szcount * pcsz,
                     double * Pvecback,
                     double * Pvectsz)
{
  //clock_t begin = clock();


  const int dim_y = pcsz->Nbins_y+1;

  int index_m, index_z, index_y;




  //clock_t end = clock();
  //double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  //printf("time spent in allocating = %e\n",time_spent);
  if (pcsz->sz_verbose > 0)
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

#pragma omp parallel \
shared(abort,pba,ppm,pnl,ptsz,pcsz)\
private(tstart, tstop,pvecsz)
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
    for (index_y=0; index_y<pcsz->Nbins_y+1; index_y ++){
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
    if (pcsz->sz_verbose > 0)
      printf("In %s: time spent in parallel region (loop over s/n's) = %e s for thread %d\n",
             __func__,tstop-tstart,omp_get_thread_num());
#endif
    free(pvecsz);
  }
  if (abort == _TRUE_) return _FAILURE_;


 //
 //  double H0_class_units =  pba->H0;
 //  int first_index_back = 0;
 //  double tau;
 //  double * pvecback;
 //  class_alloc(pvecback,
 //              pba->bg_size*sizeof(double),
 //              ptsz->error_message);
 //
 //
 // for (index_y=0;index_y <(pcsz->Nbins_y+1);index_y++){
 //    for (index_z=0;index_z<pcsz->Nbins_z;index_z++){
 //      double z_bin_min = pcsz->z_center[index_z]-0.5*pcsz->dz;
 //      double z_bin_max = pcsz->z_center[index_z]+0.5*pcsz->dz;
 //
 //
 //      int j1=0;
 //      int j2=0;
 //
 //      int c;
 //
 //      double dif_test;
 //
 //      dif_test = fabs(pcsz->steps_z[0]-z_bin_min);
 //      for (c = 1; c < pcsz->nsteps_z; c++)
 //      {
 //        if (fabs(pcsz->steps_z[c] -z_bin_min)< dif_test)
 //        {
 //          dif_test = fabs(pcsz->steps_z[c]-z_bin_min);
 //          j1 = c;
 //        }
 //      }
 //
 //
 //      dif_test = fabs(pcsz->steps_z[0]-z_bin_max);
 //      for (c = 1; c < pcsz->nsteps_z; c++)
 //      {
 //        if (fabs(pcsz->steps_z[c] -z_bin_max)< dif_test)
 //        {
 //          dif_test = fabs(pcsz->steps_z[c]-z_bin_max);
 //          j2 = c;
 //        }
 //      }
 //
 //      //printf("index_y=%d\t,index_z=%d\tj1=%d\tj2=%d\n",index_y,index_z,j1,j2);
 //
 //      //Here we compute the number counts per SNR and redshift bin:
 //      int grid_index_z,grid_index_m;
 //      double SUM2 = 0.;
 //
 //
 //      for (grid_index_z=j1;grid_index_z<j2;grid_index_z++){
 //        for (grid_index_m=0;grid_index_m<pcsz->nsteps_m;grid_index_m++){
 //          double z1 = pcsz->steps_z[grid_index_z];
 //          double z2 = pcsz->steps_z[grid_index_z+1];
 //
 //
 //
 //          double m_asked =  exp(pcsz->steps_m[grid_index_m]);
 //          double z_asked,Eh,d_A,rz,volume;
 //          //hmf 1
 //          z_asked = z1;
 //
 //          class_call(background_tau_of_z(pba,
 //                                         z_asked,
 //                                         &tau),
 //                     ptsz->error_message,
 //                     ptsz->error_message);
 //
 //          class_call(background_at_tau(pba,
 //                                       tau,
 //                                       pba->long_info,
 //                                       pba->inter_normal,
 //                                       &first_index_back,
 //                                       pvecback),
 //                     ptsz->error_message,
 //                     ptsz->error_message);
 //
 //          Eh = pvecback[pba->index_bg_H]/H0_class_units;
 //          d_A = pvecback[pba->index_bg_ang_distance]*pba->h; //class dA multiply by h to get it in units Mpc/h
 //          rz = d_A*(1.+z_asked);
 //          volume = 3.0e8/1.0e5*rz*rz/Eh;
 //
 //          double f1 = volume*get_dndlnM_at_z_and_M(z_asked,m_asked,ptsz);// = grid[grid_index_m][grid_index_z+1];
 //
 //          //hmf 2
 //          z_asked = z2;
 //
 //          class_call(background_tau_of_z(pba,
 //                                         z_asked,
 //                                         &tau),
 //                     ptsz->error_message,
 //                     ptsz->error_message);
 //
 //          class_call(background_at_tau(pba,
 //                                       tau,
 //                                       pba->long_info,
 //                                       pba->inter_normal,
 //                                       &first_index_back,
 //                                       pvecback),
 //                     ptsz->error_message,
 //                     ptsz->error_message);
 //
 //          Eh = pvecback[pba->index_bg_H]/H0_class_units;
 //          d_A = pvecback[pba->index_bg_ang_distance]*pba->h; //class dA multiply by h to get it in units Mpc/h
 //          rz = d_A*(1.+z_asked);
 //          volume = 3.0e8/1.0e5*rz*rz/Eh;
 //
 //          double f2 = volume*get_dndlnM_at_z_and_M(z_asked,m_asked,ptsz);// = grid[grid_index_m][grid_index_z+1];
 //
 //
 //          double c1 = completeness_2d[grid_index_m][grid_index_z][index_y];
 //          double c2 = completeness_2d[grid_index_m][grid_index_z+1][index_y];
 //
 //
 //
 //          if (pcsz->has_completeness == 0){
 //            c1 = 1.;
 //            c2 = 1.;
 //          }
 //
 //          SUM2 = SUM2 +0.5*(f1*c1+f2*c2)*(z2-z1)*pcsz->dlnM;
 //        } // integrate over full mass range
 //      } // integrate over redshift inside each bin
 //
 //      //survey area
 //      double deg2= 3.046174198e-4; // conversion deg2 to steradian
 //      deg2 *= 41253.0; //Planck full-sky (4 x pi x 57.3^2=41253 square degrees where 57.3  = 360 / (2 x pi))
 //
 //      pcsz->dNdzdy_theoretical[index_z][index_y]=deg2*SUM2;
 //
 //    }//end loop z bins for lkl
 //  }//end loop y bins for lkl

  //end bin in mass
  if (pcsz->sz_verbose > 0)
    printf("->SZ_counts computations done.\n");

  write_output_cluster_counts(pcsz);




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

  const int dim_1 = pcsz->Ny;
  const int dim_2 = ptsz->nthetas;


  int index1,index2;

  double ** erfs = NULL;

  class_alloc(erfs,
              dim_1*sizeof(double *),
              pcsz->error_message);


  for (index1=0;index1<dim_1;index1++)
  {
    class_alloc(erfs[index1],dim_2*sizeof(double*),pcsz->error_message);


    for (index2=0;index2<dim_2;index2++){

        erfs[index1][index2]=0.;


    }
  }

  int l_array[3];
  double theta_array[3];

  int index_z, index_m, index_y;

  const int dim_mass = pcsz->nsteps_m;
  const int dim_redshift = pcsz->nsteps_z;

  double ** completeness_2d = NULL;
  class_alloc(completeness_2d,dim_mass*sizeof(double *),pcsz->error_message);


  for (index_m=0;index_m<dim_mass;index_m++)
  {
    class_alloc(completeness_2d[index_m],dim_redshift*sizeof(double*),pcsz->error_message);

    for (index_z=0;index_z<dim_redshift;index_z++){
      completeness_2d[index_m][index_z]=0.;

    }
  }



  index_y = (int) pvecsz[pcsz->index_y];

  double y_min = pow(10., pcsz->logy[index_y] - pcsz->dlogy/2.);
  double y_max = pow(10., pcsz->logy[index_y] + pcsz->dlogy/2.);

  if (pcsz->sz_verbose > 3){
    printf("->SZ_counts grid_C_2d.\n");
    //printf("->In signal-to-noise bin:\n");
    // printf("->bin id = %d y_min = %.3e y_max = %.3e\n",index_y,y_min,y_max);
    }


  if (pcsz->sigmaM == 0.){
    for (index_z=0;index_z<pcsz->nsteps_z;index_z++){
      double zp = pcsz->steps_z[index_z];

      for (index_m=0;index_m<pcsz->nsteps_m;index_m++){

        //compute_theta_and_y_at_z_and_m
        double mp= exp(pcsz->steps_m[index_m]);
        double yp = get_y_at_m_and_z(mp,zp,ptsz,pba);
        double thp = get_theta_at_m_and_z(mp,zp,ptsz,pba);
        //Planck


        find_theta_bin(ptsz,thp,l_array,theta_array);
        int l1 = l_array[1];
        int l2 = l_array[2];
        double th1 = theta_array[1];
        double th2 = theta_array[2];


        int index_patches;
        for (index_patches =0;index_patches<ptsz->nskyfracs;index_patches++){

          double y1 = ptsz->ylims[index_patches][l1];
          double y2 = ptsz->ylims[index_patches][l2];
          double y = y1 + (y2-y1)/(th2-th1)*(thp-th1);

          double c2 = erf_compl(yp,y,pcsz->sn_cutoff);
          c2 *= erf_compl(yp,y,y_min);
          c2 *= (1.-erf_compl(yp,y,y_max));


          if (index_y == 0){
            c2 = erf_compl(yp,y,pcsz->sn_cutoff);
            c2 *= (1.-erf_compl(yp,y,y_max));


          }

          if (index_y == pcsz->Nbins_y){
            c2 = erf_compl(yp,y,pcsz->sn_cutoff) ;
            c2 *= erf_compl(yp,y,y_min);


          }

          completeness_2d[index_m][index_z] += c2*ptsz->skyfracs[index_patches];
        } // end loop patches

      }//end m loop
    }//end z loop
  }//end if sigmaM=0

  else {
    double fac =1./sqrt(2.*_PI_*pow(pcsz->sigmaM,2));

    if (pcsz->sz_verbose > 3)
      printf("->SZ_counts grid_C_2d debug 1.\n");
    double fsky = 0.;
    int index_patches;
    int index2,index1;

    ////// tabulate erfs as a function of theta and y in each s/n bin
    for (index2=0;index2<ptsz->nthetas;index2++){
      double lny = pcsz->lnymin;

      for (index1=0;index1<pcsz->Ny;index1++)
      {
        double y0=exp(lny);
        lny=lny+pcsz->dlny;
        for (index_patches=0;index_patches<ptsz->nskyfracs;index_patches++){
          fsky += ptsz->skyfracs[index_patches];
          double y1 = ptsz->ylims[index_patches][index2];
          int k = index_y;
          double qmin=pcsz->logy[k]-pcsz->dlogy/2.;
          double qmax=pcsz->logy[k]+pcsz->dlogy/2.;
          qmin=pow(10.,qmin);
          qmax=pow(10.,qmax);

          double c2=erf_compl(y0,y1,pcsz->sn_cutoff)*erf_compl(y0,y1,qmin)*(1.-erf_compl(y0,y1,qmax));

          if (k==0)  c2=erf_compl(y0,y1,pcsz->sn_cutoff)*(1.-erf_compl(y0,y1,qmax));
          if (k==pcsz->Nbins_y) c2=erf_compl(y0,y1,qmin)*erf_compl(y0,y1,pcsz->sn_cutoff);
          erfs[index1][index2]=erfs[index1][index2]+c2*ptsz->skyfracs[index_patches];


        } //end loop patches
      } //end loop y
    } //end loop thetas

    ////// end tabulate erfs as a function of theta and y


    // tabulate completeness as a function of z and m
    // integrate erfs wrt y at all (z,M) to get completeness in a z,m grid
  if (pcsz->sz_verbose > 3)
      printf("->SZ_counts grid_C_2d debug 2.\n");
    for (index_z=0;index_z<pcsz->nsteps_z;index_z++){

      double zp = pcsz->steps_z[index_z];
      // if (pcsz->sz_verbose > 3)
      //   printf("->SZ_counts grid_C_2d debug 3, z = %.4e.\n",zp);

      for (index_m=0;index_m<pcsz->nsteps_m;index_m++){

        double mp= exp(pcsz->steps_m[index_m]);
        double yp = get_y_at_m_and_z(mp,zp,ptsz,pba);
        double thp = get_theta_at_m_and_z(mp,zp,ptsz,pba);

        find_theta_bin(ptsz,thp,l_array,theta_array);

        // if (pcsz->sz_verbose > 3)
        //   printf("->SZ_counts grid_C_2d debug 4, z = %.4e.\n",zp);
        int l1 = l_array[1];
        int l2 = l_array[2];
        // printf("l1 = %d, l2 = %d\n",l1,l2);
        // exit(0);

        double th1 = theta_array[1];
        double th2 = theta_array[2];

        double y = yp;
        double mu = log(y);


        double int_comp =0.;
        double lny=pcsz->lnymin;
        int k;

        // if (pcsz->sz_verbose > 3)
        //   printf("->SZ_counts grid_C_2d debug 5, z = %.4e.\n",zp);

        // at a fixed theta(z,m)
        // integrate over y, erf(theta,y)*fac/y*exp(-arg(y))

        for (k=0;k<pcsz->Ny-1;k++){
          // printf("k = %d int_comp1 = %e\n",k,int_comp);
          double y0=exp(lny);
          // printf("k = %d int_comp2 = %e\n",k,int_comp);
          y=exp(lny+pcsz->dlny);
          // printf("k = %d int_comp3 = %e\n",k,int_comp);
          double dy=y-y0;
          // printf("k = %d int_comp4 = %e\n",k,int_comp);
          double arg0=((lny-mu)/(sqrt(2.)*ptsz->sigmaM_ym));

          double win0=erfs[k][l1]+(erfs[k][l2]-erfs[k][l1])/(th2-th1)*(thp-th1);
          double win=erfs[k+1][l1]+(erfs[k+1][l2]-erfs[k+1][l1])/(th2-th1)*(thp-th1);

          lny=lny+pcsz->dlny;
          double arg=((lny-mu)/(sqrt(2.)*pcsz->sigmaM));
          double py=(win0*fac/y0*exp(-arg0*arg0)+win*fac/y*exp(-arg*arg))*0.5;

          int_comp=int_comp+py*dy;
          // printf("k = %d int_comp15 = %e\n",k,int_comp);
        }
        if (int_comp > fsky) int_comp=fsky;
        if (int_comp < 0.) int_comp=0.;


        completeness_2d[index_m][index_z] =int_comp;

      }//end m loop

    }//end z loop

      // end tabulate completeness as a function of z and m
  //  exit(0);
      if (pcsz->sz_verbose > 3)
      printf("->SZ_counts grid_C_2d debug 3.\n");
  }// end else sigmaM != 0

  free(erfs);


    for (index_z=0;index_z<pcsz->Nbins_z;index_z++){
      double z_bin_min = pcsz->z_center[index_z]-0.5*pcsz->dz;
      double z_bin_max = pcsz->z_center[index_z]+0.5*pcsz->dz;


      int j1=0;
      int j2=0;

      int c;

      double dif_test;

      dif_test = fabs(pcsz->steps_z[0]-z_bin_min);
      for (c = 1; c < pcsz->nsteps_z; c++)
      {
        if (fabs(pcsz->steps_z[c] -z_bin_min)< dif_test)
        {
          dif_test = fabs(pcsz->steps_z[c]-z_bin_min);
          j1 = c;
        }
      }


      dif_test = fabs(pcsz->steps_z[0]-z_bin_max);
      for (c = 1; c < pcsz->nsteps_z; c++)
      {
        if (fabs(pcsz->steps_z[c] -z_bin_max)< dif_test)
        {
          dif_test = fabs(pcsz->steps_z[c]-z_bin_max);
          j2 = c;
        }
      }

      //printf("index_y=%d\t,index_z=%d\tj1=%d\tj2=%d\n",index_y,index_z,j1,j2);

      //Here we compute the number counts per SNR and redshift bin:
      int grid_index_z,grid_index_m;
      double SUM2 = 0.;


      for (grid_index_z=j1;grid_index_z<j2;grid_index_z++){
        for (grid_index_m=0;grid_index_m<pcsz->nsteps_m;grid_index_m++){
          double z1 = pcsz->steps_z[grid_index_z];
          double z2 = pcsz->steps_z[grid_index_z+1];



          double m_asked =  exp(pcsz->steps_m[grid_index_m]);

          double f1 = get_volume_at_z(z1,pba)*get_dndlnM_at_z_and_M(z1,m_asked,ptsz);// = grid[grid_index_m][grid_index_z+1];
          double f2 = get_volume_at_z(z2,pba)*get_dndlnM_at_z_and_M(z2,m_asked,ptsz);// = grid[grid_index_m][grid_index_z+1];

          double c1 = completeness_2d[grid_index_m][grid_index_z];
          double c2 = completeness_2d[grid_index_m][grid_index_z+1];



          if (pcsz->has_completeness == 0){
            c1 = 1.;
            c2 = 1.;
          }

          SUM2 = SUM2 +0.5*(f1*c1+f2*c2)*(z2-z1)*pcsz->dlnM;
        } // integrate over full mass range
      } // integrate over redshift inside each redshift bin

      pcsz->dNdzdy_theoretical[index_z][index_y]=4.*_PI_*SUM2;

    }//end loop z bins for lkl

  free(completeness_2d);



  return _SUCCESS_;
}





int write_output_cluster_counts(struct szcount * pcsz){

  char Filepath[_ARGUMENT_LENGTH_MAX_];
  int i,index_m,index_z;

if (pcsz->sz_verbose > 0)
{
  FILE *fp;
  sprintf(Filepath,"%s%s%s",pcsz->root,"dndzdy",".txt");
  fp=fopen(Filepath, "w");

  if(fp == NULL)
    exit(-1);
    int j;
  for (j=0;j<pcsz->Nbins_y+1;j++){
    for (i=0;i<pcsz->Nbins_z;i++){

      fprintf(fp,"%e\n",pcsz->dNdzdy_theoretical[i][j]);
    }}

  double total_counts = 0.;
  for (j=0;j<pcsz->Nbins_z;j++){
    if (j== 0) {
        printf("snr (mid)\t");
        for (i=0;i<pcsz->Nbins_y+1;i++){
          printf("y=%.3e\t",pow(10,pcsz->logy[i]));
        }
        printf(" ------ \n");
    }
      printf("z=%.3e\t",pcsz->z_center[j]);
    for (i=0;i<pcsz->Nbins_y+1;i++){

      total_counts += pcsz->dNdzdy_theoretical[j][i];
      printf("%e\t",pcsz->dNdzdy_theoretical[j][i]);
    }
    printf(" ------ \n");
  }
  if (pcsz->has_completeness == 0)
    total_counts = total_counts/(pcsz->Nbins_y+1.);

    printf("total counts = %e\n", total_counts);


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

    for (i=0;i<pcsz->Nbins_y+1;i++){
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

  //pcsz->nzSZ = 20.; //cosmomc settings
  //pcsz->size_logM = 105; //cosmomc settings

  //pcsz->rho_m_at_z = ptsz->Omega_m_0*ptsz->Rho_crit_0*pow((1.+pcsz->redshift_for_dndm),3);


  class_alloc(pcsz->redshift,sizeof(double)*pcsz->nzSZ,pcsz->error_message);




//Planck cut_off = 6.;
//SO cut_off = 5.;
// if(ptsz->experiment == 0) pcsz->sn_cutoff = 6.;
// if(ptsz->experiment == 1) pcsz->sn_cutoff = 5.;
  // pcsz->sn_cutoff = 5.;
  // pcsz->alpha;
  // //pcsz->ystar = pow(10.,pcsz->ystar)/pow(2., pcsz->alpha)*0.00472724;//8.9138435358806980e-004;
  // pcsz->beta = 0.66;
  // pcsz->thetastar = 6.997;
  // pcsz->alpha_theta = 1./3.;


  //grid for mass
  if (pcsz->mass_range == 0){
    pcsz->lnM_max = log(ptsz->M2SZ);
    pcsz->lnM_min = log(ptsz->M1SZ);
  }

  else {
    pcsz->lnM_max = 37.; //cosmomc/szount.f90 range
    pcsz->lnM_min = 31.54;

  }

  pcsz->dlnM = 0.05; //0.05 ref value in szcounts.f90


  pcsz->nsteps_m = floor((pcsz->lnM_max - pcsz->lnM_min) /pcsz->dlnM);

  double lnM = pcsz->lnM_min;
  int index_m;

  class_alloc(
              pcsz->steps_m,
              pcsz->nsteps_m*sizeof(double),
              pcsz->error_message
              );

  for (index_m=0; index_m<pcsz->nsteps_m; index_m++){
    pcsz->steps_m[index_m] = lnM + pcsz->dlnM/2.;
    lnM += pcsz->dlnM;
  }


  //grid for redshift
  //# Redshift bin parameters
  pcsz->z_0 = 0.;
  //pcsz->z_max : 1. for the Planck lkl (default), ste in param file.
  if (ptsz->experiment == 0){
    pcsz->z_max = 1.;
  }
  else if (ptsz->experiment == 1){
    pcsz->z_max = 2.8;
  }
  pcsz->dz = 0.1;

  pcsz->Nbins_z =floor((pcsz->z_max - pcsz->z_0)/pcsz->dz) + 1;

// printf("%d\n",pcsz->Nbins_z);
// exit(0);
  class_alloc(pcsz->z_center,pcsz->Nbins_z*sizeof(double),pcsz->error_message);
  int index_z;
  for (index_z = 0; index_z<pcsz->Nbins_z; index_z ++){
    pcsz->z_center[index_z] = pcsz->z_0 + 0.5*pcsz->dz + index_z*pcsz->dz;
    //printf("index_z=%d, z_center=%e\n",index_z,z_center[index_z]);
  }

  if(pcsz->z_0==0.) pcsz->z_center[0] += 1.e-8;

  double binz=pcsz->z_center[1]-pcsz->z_center[0];


  double z_max = pcsz->z_center[pcsz->Nbins_z-1] + 0.5*pcsz->dz;
  double z_i = pcsz->z_0;
  pcsz->nsteps_z = 0;
  while (z_i <= z_max) {
    z_i = next_z(z_i,binz);
    pcsz->nsteps_z += 1;
  }

//printf("nsteps_z=%d\n",pcsz->nsteps_z);

  class_alloc(pcsz->steps_z,
              pcsz->nsteps_z*sizeof(double),
              pcsz->error_message);

  z_i = pcsz->z_0;

  for(index_z = 0; index_z<pcsz->nsteps_z; index_z++){
    pcsz->steps_z[index_z] = z_i;
    z_i = next_z(z_i,binz);
  }

  if (pcsz->steps_z[0]==0) pcsz->steps_z[0] = 1.e-5;


  //grid for s/n

  //grid for y
  //# y bin parameters
  //# Logymin corresponds to log10 of S/N_min (5 or 6)
  //# Logymax corresponds to log10 of S/N_max  (~32 for planck) (but the s/n is higher in the next to last bin)
  if (ptsz->experiment==0){
  pcsz->logy_min = 0.7;
  pcsz->logy_max = 1.5;
  pcsz->dlogy = 0.25;
}
else if (ptsz->experiment==1){
  pcsz->logy_min = 0.6989700043360189;
  pcsz->logy_max = 1.8124259665302023;
  pcsz->dlogy = 0.1;
}
  pcsz->Nbins_y = floor((pcsz->logy_max - pcsz->logy_min)/pcsz->dlogy)+1;
  // printf("%d\n",pcsz->Nbins_y);
  //exit(0);
  class_alloc(pcsz->logy,(pcsz->Nbins_y+1)*sizeof(double),pcsz->error_message);
  int index_y;
  double y_i = pcsz->logy_min + pcsz->dlogy/2.;
  for (index_y = 0; index_y<pcsz->Nbins_y+1; index_y ++){
    pcsz->logy[index_y] = y_i;
    y_i += pcsz->dlogy;
    //printf("index_y=%d, logy=%e\n",index_y,logy[index_y]);
  }


  //y_500 grid
  pcsz->lnymin = -11.5;
  pcsz->lnymax = 10.;
  pcsz->dlny = 0.05;

  pcsz->Ny = (pcsz->lnymax-pcsz->lnymin)/pcsz->dlny;


  class_alloc(pcsz->dNdzdy_theoretical,
              pcsz->Nbins_z*sizeof(double *),
              pcsz->error_message);


  for (index_z=0;
       index_z<pcsz->Nbins_z;
       index_z++)
  {
    class_alloc(pcsz->dNdzdy_theoretical[index_z],
                (pcsz->Nbins_y+1)*sizeof(double),
                pcsz->error_message);
  }


  pcsz->index_y = 0;
  pcsz->pvecsz_size = pcsz->index_y + 1;

    return _SUCCESS_;
}


int find_theta_bin(struct tszspectrum * ptsz, double thp, int * l_array, double * theta_array){
  int l1,l2;
  double th1,th2;

  if (thp > ptsz->theta_bin_max){
    l1 = ptsz->nthetas - 1;
    l2 = ptsz->nthetas - 2;
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
