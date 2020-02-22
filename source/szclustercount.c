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
  free(pcsz->dndz);
  free(pcsz->dndmdz);
  free(pcsz->logM_at_z);
  free(pcsz->dNdzdy_theoretical);
  free(pcsz->dNdzdm_theoretical);
  free(pcsz->temp_0_theoretical);
  free(pcsz->temp_1_theoretical);
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

  const int dim1 = pcsz->nsteps_m;
  const int dim2 = pcsz->nsteps_z;
  const int dim3 = pcsz->Nbins_y+1;

  int index_m, index_z, index_y;

  double *** completeness_2d = NULL;
  double *** d_completeness_2d_dq = NULL;
  class_alloc(completeness_2d,dim1*sizeof(double **),pcsz->error_message);
  class_alloc(d_completeness_2d_dq,dim1*sizeof(double **),pcsz->error_message);

  for (index_m=0;index_m<dim1;index_m++)
  {
    class_alloc(completeness_2d[index_m],dim2*sizeof(double*),pcsz->error_message);
    class_alloc(d_completeness_2d_dq[index_m],dim2*sizeof(double*),pcsz->error_message);

    for (index_z=0;index_z<dim2;index_z++){
      class_alloc(completeness_2d[index_m][index_z],
                  dim3*sizeof(double),
                  pcsz->error_message);
      class_alloc(d_completeness_2d_dq[index_m][index_z],
                  dim3*sizeof(double),
                  pcsz->error_message);

      for (index_y=0;index_y<dim3;index_y++)
      { completeness_2d[index_m][index_z][index_y]=0.;
        d_completeness_2d_dq[index_m][index_z][index_y]=0.;
      }
    }
  }


  const int dim_1 = pcsz->Ny;
  const int dim_2 = ptsz->nthetas;
  const int dim_3 = pcsz->Nbins_y+1;

  int index1,index2,index3;

  double *** erfs = NULL;
  double *** d_erfs_dq = NULL;
  class_alloc(erfs,
              dim_1*sizeof(double **),
              pcsz->error_message);
  class_alloc(d_erfs_dq,
              dim_1*sizeof(double **),
              pcsz->error_message);

  for (index1=0;index1<dim_1;index1++)
  {
    class_alloc(erfs[index1],dim_2*sizeof(double*),pcsz->error_message);
    class_alloc(d_erfs_dq[index1],dim_2*sizeof(double*),pcsz->error_message);

    for (index2=0;index2<dim_2;index2++){
      class_alloc(erfs[index1][index2],dim_3*sizeof(double),pcsz->error_message);
      class_alloc(d_erfs_dq[index1][index2],dim_3*sizeof(double),pcsz->error_message);

      for (index3=0;index3<dim_3;index3++)
      {
        erfs[index1][index2][index3]=0.;
        d_erfs_dq[index1][index2][index3]=0.;
      }
    }
  }
  //clock_t end = clock();
  //double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  //printf("time spent in allocating = %e\n",time_spent);

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
shared(abort,completeness_2d,d_completeness_2d_dq,erfs,d_erfs_dq,pba,ppm,pnl,ptsz,pcsz)\
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
      class_call_parallel(grid_C_2d(d_completeness_2d_dq,
                                    completeness_2d,
                                    d_erfs_dq,
                                    erfs,
                                    pvecsz,
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



    double ** grid = NULL;
    double ** grid_temp_0 = NULL;
    double ** grid_temp_1 = NULL;



  class_alloc(grid,
              pcsz->nsteps_m*sizeof(double *),
              pcsz->error_message);

  class_alloc(grid_temp_0,
              pcsz->nsteps_m*sizeof(double *),
              pcsz->error_message);

  class_alloc(grid_temp_1,
              pcsz->nsteps_m*sizeof(double *),
              pcsz->error_message);

  for (index_m=0;
       index_m<pcsz->nsteps_m;
       index_m++)
  {
    class_alloc(grid[index_m],
                pcsz->nsteps_z*sizeof(double),
                pcsz->error_message);
    class_alloc(grid_temp_0[index_m],
                pcsz->nsteps_z*sizeof(double),
                pcsz->error_message);
    class_alloc(grid_temp_1[index_m],
                pcsz->nsteps_z*sizeof(double),
                pcsz->error_message);
  }




  /////start get_grid

  double * log_X;
  double * log_Y;
  double * log_Z;

  class_alloc(log_X,sizeof(double)*pcsz->nzSZ,pcsz->error_message);
  class_alloc(log_Y,sizeof(double)*pcsz->size_logM,pcsz->error_message);
  class_alloc(log_Z,sizeof(double *)*pcsz->nzSZ*pcsz->size_logM,ptsz->error_message);


  double z1SZ = ptsz->z1SZ;
  double z2SZ = ptsz->z2SZ;

  Pvectsz[ptsz->index_md] = ptsz->index_md_hmf;

  int i;
  for (i=0;i<pcsz->nzSZ;i++)
  {
    pcsz->redshift[i] =
    exp(log(1.+ z1SZ)+i*(log(1.+ z2SZ)-log(1.+ z1SZ))/(1.*pcsz->nzSZ-1.))-1.;
    log_X[i]=log(1.+pcsz->redshift[i]);
    Pvectsz[ptsz->index_z] = pcsz->redshift[i];
    Pvectsz[ptsz->index_multipole] = 0;

    /*
    class_call(integrate_over_m_at_z_qgaus_sz(Pvecback,
                                              Pvectsz,
                                              pba,
                                              pnl,
                                              ppm,
                                              ptsz),
               pcsz->error_message,
               pcsz->error_message);

    pcsz->dndz[i] = Pvectsz[ptsz->index_integrals_over_m_first];
    */
    int j;
    for (j=0;j<pcsz->size_logM;j++)
    {
      if (i==0) log_Y[j]=pcsz->logM_at_z[j]*log(10.);

      class_call(integrand_at_m_and_z(pcsz->logM_at_z[j]*log(10.),
                                      Pvecback,
                                      Pvectsz,
                                      pba,
                                      ppm,
                                      pnl,
                                      ptsz),
                 pcsz->error_message,
                 pcsz->error_message);

      // dn/dlogM in units of h^3 Mpc^-3
      pcsz->dndmdz[j][i] = Pvectsz[ptsz->index_hmf];


      if (pcsz->dndmdz[j][i] == 0.)
        log_Z[i*pcsz->size_logM +j]=log(1e-100);
      else
        log_Z[i*pcsz->size_logM +j]=log(pcsz->dndmdz[j][i]);
    }
  }

  /////end get_grid
/*
/////start output some quantities to file
  FILE *fp;
  fp=fopen("output/data_sigmaR_class_sz.dat", "w");
  double z_out = 2.;

  Pvectsz[ptsz->index_md] = ptsz->index_md_hmf;
  Pvectsz[ptsz->index_z] = z_out;
  Pvectsz[ptsz->index_multipole] = 0;
  double lnM;


    for (int j=0;j<100;j++)
    {

      lnM =  log(1.e13) + j*(log(1.e16)-log(1.e13))/99.;

      class_call(integrand_at_m_and_z(lnM,
                                      Pvecback,
                                      Pvectsz,
                                      pba,
                                      ppm,
                                      pnl,
                                      ptsz),
                 pcsz->error_message,
                 pcsz->error_message);

      fprintf(fp,"%d\t%e\t%e\t%e\t%e\n",j+1,exp(lnM),Pvectsz[ptsz->index_hmf],exp(Pvectsz[ptsz->index_logSigma2]/2.),1./2.*Pvectsz[ptsz->index_dlogSigma2dlogRh]);

  }
  fclose(fp);
  /////end output to file
*/

  int first_index_back = 0;
  double tau;
  double * pvecback;
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              ptsz->error_message);


  class_call(background_tau_of_z(pba,0.,&tau),
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
  double H0_class_units = pvecback[pba->index_bg_H];

  double m_asked;
  double z_asked;

  for (index_z=0;index_z<pcsz->nsteps_z;index_z++){

    class_call(background_tau_of_z(pba,
                                   pcsz->steps_z[index_z],
                                   &tau),
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

    double Eh = pvecback[pba->index_bg_H]/H0_class_units;
    double d_A = pvecback[pba->index_bg_ang_distance]*pba->h; //class dA multiply by h to get it in units Mpc/h
    double rz = d_A*(1.+pcsz->steps_z[index_z]);
    double volume = 3.0e8/1.0e5*rz*rz/Eh;

    //survey area
    double deg2= 3.046174198e-4;
    if (ptsz->experiment == 0) deg2 *= 41253.0; //Planck full-sky
    if (ptsz->experiment == 1) deg2 *= 599.353; //SO

    double HMF;

    for (index_m=0;index_m<pcsz->nsteps_m;index_m++){
      m_asked = pcsz->steps_m[index_m];

      if ((index_z == 0) && (pcsz->steps_z[0]<pcsz->redshift[0]))
        z_asked =  log(1.+pcsz->redshift[0]);
      else
        z_asked = log(1.+pcsz->steps_z[index_z]);


      HMF =  exp(pwl_interp_2d(pcsz->size_logM,
                               pcsz->nzSZ,
                               log_Y,
                               log_X,
                               log_Z,
                               1,
                               &m_asked,
                               &z_asked));


      grid[index_m][index_z] = HMF * deg2 * volume;
      //grid[index_m][index_z] = 1.;

      //printf("%d \t %d\n",index_m,index_z);

      double mp= exp(pcsz->steps_m[index_m]);
      double Tp = 5.*pow(Eh*mp/3.e14* (1./pba->h),2./3.);
      double mp_bias = mp/ptsz->HSEbias;

      double ystar2 = pcsz->ystar;
      ystar2 *=  pow(100.*pba->h/70.,-2.+pcsz->alpha);
      double y500_for_mp_at_zp =  ystar2 * pow(mp_bias/3.e14* (1./pba->h),pcsz->alpha);
      y500_for_mp_at_zp *=   pow(Eh,pcsz->beta) *pow(100.*d_A/(500.0*100.*pba->h),-2.);

      grid_temp_0[index_m][index_z] = y500_for_mp_at_zp*HMF * deg2 * volume;
      grid_temp_1[index_m][index_z] = Tp*y500_for_mp_at_zp*HMF * deg2 * volume;

    }//end loop over masses - steps_m
  }//end loop over redshifts - steps_z


 //integrate_m_zq in szcounts.f90
 //printf("Nbins_y=%d\tNbins_z=%d\n",pcsz->Nbins_y,pcsz->Nbins_z);
 for (index_y=0;index_y <(pcsz->Nbins_y+1);index_y++){
    for (index_z=0;index_z<pcsz->Nbins_z;index_z++){
//printf("z=%e\n",pcsz->z_center[index_z]);
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


      int jj,ii;
      double SUM2 = 0.;
      double SUM2_temp_0 = 0.;
      double SUM2_temp_1 = 0.;

      for (jj=j1;jj<j2;jj++){
        for (ii=0;ii<pcsz->nsteps_m;ii++){
          double x1 = pcsz->steps_z[jj];
          double x2 = pcsz->steps_z[jj+1];

          double f1 = grid[ii][jj];
          double f2 = grid[ii][jj+1];

          double f1_temp_0 = grid_temp_0[ii][jj];
          double f2_temp_0 = grid_temp_0[ii][jj+1];

          double f1_temp_1 = grid_temp_1[ii][jj];
          double f2_temp_1 = grid_temp_1[ii][jj+1];


          double c1 =  completeness_2d[ii][jj][index_y];
          double c2 = completeness_2d[ii][jj+1][index_y];

          double d_c1_dq =  d_completeness_2d_dq[ii][jj][index_y];
          double d_c2_dq = d_completeness_2d_dq[ii][jj+1][index_y];


          if (pcsz->has_completeness == 0){
            c1 = 1.;
            c2 = 1.;
            d_c1_dq = 0.;
            d_c2_dq = 0.;
          }

          SUM2 = SUM2 +0.5*(f1*c1+f2*c2)*(x2-x1)*pcsz->dlnM;
          //SUM2 = SUM2 +0.5*(c1+c2)*(x2-x1)*pcsz->dlnM;
          SUM2_temp_0 = SUM2_temp_0 +0.5*(f1_temp_0*d_c1_dq+f2_temp_0*d_c2_dq)*(x2-x1)*pcsz->dlnM;
          SUM2_temp_1 = SUM2_temp_1 +0.5*(f1_temp_1*d_c1_dq+f2_temp_1*d_c2_dq)*(x2-x1)*pcsz->dlnM;
        }
      }

      pcsz->dNdzdy_theoretical[index_z][index_y]=SUM2;
      pcsz->temp_0_theoretical[index_z][index_y]=SUM2_temp_0;
      if (SUM2_temp_0 != 0.){
        pcsz->temp_1_theoretical[index_z][index_y]=SUM2_temp_1/SUM2_temp_0;
      }
      else {
        pcsz->temp_0_theoretical[index_z][index_y]=0.;
        pcsz->temp_1_theoretical[index_z][index_y]=0.;
      }
    }//end loop z bins for lkl
  }//end loop y bins for lkl



  for (index_m=0;index_m <pcsz->nsteps_m;index_m++){ //loop over massess
    for (index_z=0;index_z<pcsz->Nbins_z;index_z++){

      double z_bin_min = pcsz->z_center[index_z]-0.5*pcsz->dz;
      double z_bin_max = pcsz->z_center[index_z]+0.5*pcsz->dz;

      double dif_test = fabs(pcsz->steps_z[0]-z_bin_min);
      int j1=0;
      int j2=0;

      int c;

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

      double SUM2 = 0.;

      int jj;
      for (jj=j1;jj<j2-1;jj++){

        double x1 = pcsz->steps_z[jj];
        double x2 = pcsz->steps_z[jj+1];

        double f1 = grid[index_m][jj];
        double f2 = grid[index_m][jj+1];

        double c1 = 1.;
        double c2 = 1.;

        SUM2 = SUM2 +0.5*(f1*c1+f2*c2)*(x2-x1)*pcsz->dlnM;


      }

      pcsz->dNdzdm_theoretical[index_z][index_m]=SUM2;

    }//end loop z
  }//end loop m

  //end bin in mass


  write_output_cluster_counts(pcsz);




  free(log_X);
  free(log_Y);
  free(log_Z);

  free(pvecback);


  for (index_m=0;index_m<pcsz->nsteps_m;index_m++)
  {
    free(grid[index_m]);
    free(grid_temp_0[index_m]);
    free(grid_temp_1[index_m]);
  }

  free(grid);
  free(grid_temp_0);
  free(grid_temp_1);

  for (index_m=0;
       index_m<dim1;
       index_m++)
  {

    for (index_z=0;
         index_z<dim2;
         index_z++){
      free(completeness_2d[index_m][index_z]);
      free(d_completeness_2d_dq[index_m][index_z]);
    }

    free(completeness_2d[index_m]);


    free(d_completeness_2d_dq[index_m]);
  }


  free(completeness_2d);
  free(d_completeness_2d_dq);



  for (index1=0;
       index1<dim_1;
       index1++)
  {

    for (index2=0;
         index2<dim_2;
         index2++){
      free(erfs[index1][index2]);
      free(d_erfs_dq[index1][index2]);
    }

    free(erfs[index1]);
    free(d_erfs_dq[index1]);
  }


  free(erfs);
  free(d_erfs_dq);



  return _SUCCESS_;
}




int grid_C_2d(
              double *** d_completeness_2d_dq,
              double *** completeness_2d,
              double *** d_erfs_dq,
              double *** erfs,
              double * pvecsz,
              struct background *pba,
              struct primordial * ppm,
              struct nonlinear * pnl,
              struct tszspectrum * ptsz,
              struct szcount * pcsz
              ){
  int i;

  int l_array[3];
  double theta_array[3];

  int index_z, index_m, index_y;

  index_y = (int) pvecsz[pcsz->index_y];

  double y_min = pow(10., pcsz->logy[index_y] - pcsz->dlogy/2.);
  double y_max = pow(10., pcsz->logy[index_y] + pcsz->dlogy/2.);



  double tau;
  int first_index_back = 0;
  double * pvecback;
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              ptsz->error_message);

  class_call(background_tau_of_z(pba,0.,&tau),
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
  double H0_class_units = pvecback[pba->index_bg_H];



  if (pcsz->sigmaM == 0.){
    for (index_z=0;index_z<pcsz->nsteps_z;index_z++){
      double zp = pcsz->steps_z[index_z];
      class_call(background_tau_of_z(pba,zp,&tau),
                 pcsz->error_message,
                 pcsz->error_message);

      class_call(background_at_tau(pba,
                                   tau,
                                   pba->long_info,
                                   pba->inter_normal,
                                   &first_index_back,
                                   pvecback),
                 pcsz->error_message,
                 pcsz->error_message);


      double Eh = pvecback[pba->index_bg_H]/H0_class_units;
      double d_A = pvecback[pba->index_bg_ang_distance]*pba->h;

      for (index_m=0;index_m<pcsz->nsteps_m;index_m++){

        //compute_theta_and_y_at_z_and_m
        double mp= exp(pcsz->steps_m[index_m]);
        double mp_bias = mp/ptsz->HSEbias;
        double H0 = pba->h*100.;
        double thetastar2 = pcsz->thetastar * pow(H0/70.,-2./3.);
        double theta500_for_mp_at_zp =  thetastar2 * pow(mp_bias/3.e14* (100./H0),pcsz->alpha_theta);
        theta500_for_mp_at_zp *=    pow(Eh,-2./3) *pow(100.*d_A/(500.0*H0),-1.);
        double thp = theta500_for_mp_at_zp;

        double ystar2 = pcsz->ystar;
        ystar2 *=  pow(H0/70.,-2.+pcsz->alpha);
        double y500_for_mp_at_zp =  ystar2 * pow(mp_bias/3.e14* (100./H0),pcsz->alpha);
        y500_for_mp_at_zp *=   pow(Eh,pcsz->beta) *pow(100.*d_A/(500.0*H0),-2.);
        double yp = y500_for_mp_at_zp;


        //Planck
        if(ptsz->experiment == 0){

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

          double d_c2_dq =  d_erf_compl_dq(yp,y,pcsz->sn_cutoff)
                            *erf_compl(yp,y,y_min)
                            *(1.-erf_compl(yp,y,y_max))
                            +erf_compl(yp,y,pcsz->sn_cutoff)
                            *d_erf_compl_dq(yp,y,y_min)
                            *(1.-erf_compl(yp,y,y_max))
                            +erf_compl(yp,y,pcsz->sn_cutoff)
                            *erf_compl(yp,y,y_min)
                            *(-d_erf_compl_dq(yp,y,y_max));

          if (index_y == 0){
            c2 = erf_compl(yp,y,pcsz->sn_cutoff);
            c2 *= (1.-erf_compl(yp,y,y_max));

            d_c2_dq = d_erf_compl_dq(yp,y,pcsz->sn_cutoff)
                      *(1.-erf_compl(yp,y,y_max))
                      +erf_compl(yp,y,pcsz->sn_cutoff)
                      *(-d_erf_compl_dq(yp,y,y_max));

          }

          if (index_y == pcsz->Nbins_y){
            c2 = erf_compl(yp,y,pcsz->sn_cutoff) ;
            c2 *= erf_compl(yp,y,y_min);

            d_c2_dq = d_erf_compl_dq(yp,y,pcsz->sn_cutoff)
                      *erf_compl(yp,y,y_min)
                      +erf_compl(yp,y,pcsz->sn_cutoff)
                      *d_erf_compl_dq(yp,y,y_min);

          }

          completeness_2d[index_m][index_z][index_y] += c2*ptsz->skyfracs[index_patches];
          d_completeness_2d_dq[index_m][index_z][index_y] += d_c2_dq*ptsz->skyfracs[index_patches];
        } // end loop patches
      }//end ptsz->experiment = Planck

      //The Simons Observatory
      if(ptsz->experiment == 1){

        double total_area = 0.;
        //SO scaling relations
        double Qp = pwl_value_1d ( ptsz->SO_Q_size, ptsz->SO_thetas, ptsz->SO_Qfit, thp );
        //printf("Qp = %e\n", Qp);
        double A = 4.95e-5;
        double B = 0.08;
        double t = -0.00848*pow(mp/(3.e14*70./(pba->h*100.))*Eh,-0.585);
        double f_rel = 1. + 3.79*t -28.2*t*t;
        double yp = A*pow(Eh,2.)*pow(mp/(3.e14*70./(pba->h*100.)),1.+B)*Qp*f_rel;
        //y0_SO = A_a*(Eh(z)**2.)*((m2/(3.e14*msun*70./cosmopar%H0))**(1. + B_a))*Qfunc_SO(m*msun, z, theta, Q0)*relfn(m*msun, z)


      int index_patches;
      for (index_patches =0;index_patches<ptsz->SO_RMS_size;index_patches++){

        double y = ptsz->SO_RMS[index_patches];

        double c2 = erf_compl(yp,y,pcsz->sn_cutoff);
        c2 *= erf_compl(yp,y,y_min);
        c2 *= (1.-erf_compl(yp,y,y_max));

        double d_c2_dq =  d_erf_compl_dq(yp,y,pcsz->sn_cutoff)
                          *erf_compl(yp,y,y_min)
                          *(1.-erf_compl(yp,y,y_max))
                          +erf_compl(yp,y,pcsz->sn_cutoff)
                          *d_erf_compl_dq(yp,y,y_min)
                          *(1.-erf_compl(yp,y,y_max))
                          +erf_compl(yp,y,pcsz->sn_cutoff)
                          *erf_compl(yp,y,y_min)
                          *(-d_erf_compl_dq(yp,y,y_max));

        if (index_y == 0){
          c2 = erf_compl(yp,y,pcsz->sn_cutoff);
          c2 *= (1.-erf_compl(yp,y,y_max));

          d_c2_dq = d_erf_compl_dq(yp,y,pcsz->sn_cutoff)
                    *(1.-erf_compl(yp,y,y_max))
                    +erf_compl(yp,y,pcsz->sn_cutoff)
                    *(-d_erf_compl_dq(yp,y,y_max));

        }

        if (index_y == pcsz->Nbins_y){
          c2 = erf_compl(yp,y,pcsz->sn_cutoff) ;
          c2 *= erf_compl(yp,y,y_min);

          d_c2_dq = d_erf_compl_dq(yp,y,pcsz->sn_cutoff)
                    *erf_compl(yp,y,y_min)
                    +erf_compl(yp,y,pcsz->sn_cutoff)
                    *d_erf_compl_dq(yp,y,y_min);

        }
        //test 1d
        y_min = 5.;
        y_max = 1e4;
        c2 = erf_compl(yp,y,pcsz->sn_cutoff) ;
        //c2 *= erf_compl(yp,y,y_min);


        //c2 =1.; //no selfn
        completeness_2d[index_m][index_z][index_y] += c2*ptsz->SO_skyfrac[index_patches];
        total_area += ptsz->SO_skyfrac[index_patches];
        d_completeness_2d_dq[index_m][index_z][index_y] += d_c2_dq*ptsz->SO_skyfrac[index_patches];
      } // end loop patches

      //printf("total_area = %e\n",total_area);
      completeness_2d[index_m][index_z][index_y] = completeness_2d[index_m][index_z][index_y]/total_area;
    }//end ptsz->experiment = SO



      }//end m loop
    }//end z loop
  }//end if sigmaM=0

  else {
    double fac =1./sqrt(2.*_PI_*pow(pcsz->sigmaM,2));


    double fsky = 0.;
    int index_patches;
    int index2,index1;

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

          double d_c2_dq=  d_erf_compl_dq(y0,y1,pcsz->sn_cutoff)*erf_compl(y0,y1,qmin)*(1.-erf_compl(y0,y1,qmax))
                          +erf_compl(y0,y1,pcsz->sn_cutoff)*d_erf_compl_dq(y0,y1,qmin)*(1.-erf_compl(y0,y1,qmax))
                          +erf_compl(y0,y1,pcsz->sn_cutoff)*erf_compl(y0,y1,qmin)*(-d_erf_compl_dq(y0,y1,qmax));

            if (k==0)  d_c2_dq=d_erf_compl_dq(y0,y1,pcsz->sn_cutoff)*(1.-erf_compl(y0,y1,qmax))
            +erf_compl(y0,y1,pcsz->sn_cutoff)*(-d_erf_compl_dq(y0,y1,qmax));
          if (k==pcsz->Nbins_y) d_c2_dq=erf_compl(y0,y1,qmin)*d_erf_compl_dq(y0,y1,pcsz->sn_cutoff)
            +d_erf_compl_dq(y0,y1,qmin)*erf_compl(y0,y1,pcsz->sn_cutoff);

          erfs[index1][index2][k]=erfs[index1][index2][k]+c2*ptsz->skyfracs[index_patches];
          d_erfs_dq[index1][index2][k]=d_erfs_dq[index1][index2][k]+d_c2_dq*ptsz->skyfracs[index_patches];

        } //end loop patches
      } //end loop y
    } //end loop thetas


    for (index_z=0;index_z<pcsz->nsteps_z;index_z++){

      double zp = pcsz->steps_z[index_z];
      class_call(background_tau_of_z(pba,
                                     zp,
                                     &tau),
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


      double Eh = pvecback[pba->index_bg_H]/H0_class_units;
      double d_A = pvecback[pba->index_bg_ang_distance]*pba->h;

      for (index_m=0;index_m<pcsz->nsteps_m;index_m++){

        //compute_theta_and_y_at_z_and_m
        double mp= exp(pcsz->steps_m[index_m]);
        double mp_bias = mp/ptsz->HSEbias;
        double H0 = pba->h*100.;
        double thetastar2 = pcsz->thetastar * pow(H0/70.,-2./3.);
        double theta500_for_mp_at_zp =  thetastar2 * pow(mp_bias/3.e14* (100./H0),pcsz->alpha_theta);
        theta500_for_mp_at_zp *=    pow(Eh,-2./3) *pow(100.*d_A/(500.0*H0),-1.);
        double thp = theta500_for_mp_at_zp;
        double ystar2 = pcsz->ystar;
        ystar2 *=  pow(H0/70.,-2.+pcsz->alpha);
        double y500_for_mp_at_zp =  ystar2 * pow(mp_bias/3.e14* (100./H0),pcsz->alpha);
        y500_for_mp_at_zp *=   pow(Eh,pcsz->beta) *pow(100.*d_A/(500.0*H0),-2.);
        double yp = y500_for_mp_at_zp;



        find_theta_bin(ptsz,thp,l_array,theta_array);
        int l1 = l_array[1];
        int l2 = l_array[2];

        double th1 = theta_array[1];
        double th2 = theta_array[2];

        double y = yp;
        double mu = log(y);


        double int_comp =0.;
        double d_int_comp_dq =0.;
        double lny=pcsz->lnymin;
        int k;
        for (k=0;k<pcsz->Ny-1;k++){
          double y0=exp(lny);
          y=exp(lny+pcsz->dlny);
          double dy=y-y0;
          double arg0=((lny-mu)/(sqrt(2.)*pcsz->sigmaM));
          double win0=erfs[k][l1][index_y]+(erfs[k][l2][index_y]-erfs[k][l1][index_y])/(th2-th1)*(thp-th1);
          double win=erfs[k+1][l1][index_y]+(erfs[k+1][l2][index_y]-erfs[k+1][l1][index_y]
                                             )/(th2-th1)*(thp-th1);

          double d_win0_dq=d_erfs_dq[k][l1][index_y]+(d_erfs_dq[k][l2][index_y]-d_erfs_dq[k][l1][index_y])/(th2-th1)*(thp-th1);
          double d_win_dq=d_erfs_dq[k+1][l1][index_y]+(d_erfs_dq[k+1][l2][index_y]-d_erfs_dq[k+1][l1][index_y]
                                                       )/(th2-th1)*(thp-th1);
          lny=lny+pcsz->dlny;
          double arg=((lny-mu)/(sqrt(2.)*pcsz->sigmaM));
          double py=(win0*fac/y0*exp(-arg0*arg0)+win*fac/y*exp(-arg*arg))*0.5;
          int_comp=int_comp+py*dy;

          double d_py_dq=(d_win0_dq*fac/y0*exp(-arg0*arg0)+d_win_dq*fac/y*exp(-arg*arg))*0.5;
          d_int_comp_dq=d_int_comp_dq+d_py_dq*dy;
        }
        if (int_comp > fsky) int_comp=fsky;
        if (int_comp < 0.) int_comp=0.;
        if (d_int_comp_dq > fsky) d_int_comp_dq=fsky; //TBC
        if (d_int_comp_dq < 0.) d_int_comp_dq=0.;
        //printf("int_comp = %e\n",int_comp);

        completeness_2d[index_m][index_z][index_y] =int_comp;
        d_completeness_2d_dq[index_m][index_z][index_y] =d_int_comp_dq;

      }//end m loop
    }//end z loop
  }// end else sigmaM != 0

  free(pvecback);


  return _SUCCESS_;
}


int write_output_cluster_counts(struct szcount * pcsz){

  int i,j,index_m,index_z;

if (pcsz->sz_verbose > 0)
{
  FILE *fp;
  fp=fopen("output/dndz_SZ_CLASS.txt", "w");
  if(fp == NULL)
    exit(-1);
    int j;
  for (j=0;j<pcsz->Nbins_y+1;j++){

    for (i=0;i<pcsz->Nbins_z;i++){

      fprintf(fp,"%e\n",pcsz->dNdzdy_theoretical[i][j]);
    }}

  double total_counts = 0.;
  for (j=0;j<pcsz->Nbins_z;j++){
    for (i=0;i<pcsz->Nbins_y+1;i++){
      total_counts += pcsz->dNdzdy_theoretical[j][i];
      printf("%e\t",pcsz->dNdzdy_theoretical[j][i]);
    }
    printf(" ------ \n");
  }
  printf("total counts = %e\n", total_counts);


  fp=fopen("output/dndz_SZ_bins_z_center.txt", "w");
  if(fp == NULL)
    exit(-1);

    for (i=0;
         i<pcsz->Nbins_z;
         i++){

      fprintf(fp,
              "%e\n",
              pcsz->z_center[i]);
    }




  fp=fopen("output/dndz_SZ_bins_y_center.txt", "w");
  if(fp == NULL)
    exit(-1);

    for (i=0;
         i<pcsz->Nbins_y+1;
         i++){

      fprintf(fp,
              "%e\n",
              pow(10.,pcsz->logy[i]));
    }


  fp=fopen("output/dndz_SZ_bins_m_center.txt", "w");
  if(fp == NULL)
    exit(-1);

    for (i=0;
         i<pcsz->nsteps_m;
         i++){

      fprintf(fp,
              "%e\n",
              pcsz->steps_m[i]);
    }





  fp=fopen("output/dNdzdm.txt", "w");
  for (index_z=0;index_z<pcsz->Nbins_z;index_z++){
    for (index_m=0;index_m<pcsz->nsteps_m;index_m++){

      fprintf(fp,
              "%e\t %e\t %e\n",
              pcsz->dNdzdm_theoretical[index_z][index_m],
              log10(exp(pcsz->steps_m[index_m])),
              pcsz->z_center[index_z]);
    }
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

  pcsz->rho_m_at_z = ptsz->Omega_m_0*ptsz->Rho_crit_0*pow((1.+pcsz->redshift_for_dndm),3);


  class_alloc(pcsz->redshift,sizeof(double)*pcsz->nzSZ,pcsz->error_message);
  class_alloc(pcsz->dndz,sizeof(double)*pcsz->nzSZ,pcsz->error_message);
  class_alloc(pcsz->logM_at_z,sizeof(double)*pcsz->size_logM,pcsz->error_message);
  class_alloc(pcsz->dndmdz,pcsz->size_logM*sizeof(double *),pcsz->error_message);

  int i;
  for (i=0;i<pcsz->size_logM;i++){
    pcsz->logM_at_z[i] = 10.+i*7./(pcsz->size_logM-1);
    class_alloc(pcsz->dndmdz[i],pcsz->nzSZ*sizeof(double),pcsz->error_message);
  }



//Planck cut_off = 6.;
//SO cut_off = 5.;
if(ptsz->experiment == 0) pcsz->sn_cutoff = 6.;
if(ptsz->experiment == 1) pcsz->sn_cutoff = 5.;

  pcsz->alpha;
  pcsz->ystar = pow(10.,pcsz->ystar)/pow(2., pcsz->alpha)*0.00472724;//8.9138435358806980e-004;
  pcsz->beta = 0.66;
  pcsz->thetastar = 6.997;
  pcsz->alpha_theta = 1./3.;


  //grid for mass
  if (pcsz->mass_range == 0){
    pcsz->lnM_max = log(ptsz->M2SZ);
    pcsz->lnM_min = log(ptsz->M1SZ);
  }

  else {
    pcsz->lnM_max = 37.;
    pcsz->lnM_min = 31.54; //[TBC]

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
  pcsz->z_max = 2.; //z_max = 1. for the Planck lkl
  pcsz->dz = 0.1;

  pcsz->Nbins_z =floor((pcsz->z_max - pcsz->z_0)/pcsz->dz) + 1;

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
  //# Logymin corresponds to S/N = 6
  //# Logymax corresponds to S/N = 32 (but the s/n is higher in the next to last bin)
  pcsz->logy_min = 0.7;
  pcsz->logy_max = 1.5;
  pcsz->dlogy = 0.25;
  pcsz->Nbins_y = floor((pcsz->logy_max - pcsz->logy_min)/pcsz->dlogy)+1;
  class_alloc(pcsz->logy,(pcsz->Nbins_y+1)*sizeof(double),pcsz->error_message);
  int index_y;
  double y_i = pcsz->logy_min + pcsz->dlogy/2.;
  for (index_y = 0; index_y<pcsz->Nbins_y+1; index_y ++){
    pcsz->logy[index_y] = y_i;
    y_i += pcsz->dlogy;
    //printf("index_y=%d, logy=%e\n",index_y,logy[index_y]);
  }


  //noise??
  pcsz->lnymin = -11.5;
  pcsz->lnymax = 10.;
  pcsz->dlny = 0.05;

  pcsz->Ny = (pcsz->lnymax-pcsz->lnymin)/pcsz->dlny;


  class_alloc(pcsz->dNdzdy_theoretical,
              pcsz->Nbins_z*sizeof(double *),
              pcsz->error_message);
  class_alloc(pcsz->dNdzdm_theoretical,
              pcsz->Nbins_z*sizeof(double *),
              pcsz->error_message);
  class_alloc(pcsz->temp_0_theoretical,
              pcsz->Nbins_z*sizeof(double *),
              pcsz->error_message);
  class_alloc(pcsz->temp_1_theoretical,
              pcsz->Nbins_z*sizeof(double *),
              pcsz->error_message);


  for (index_z=0;
       index_z<pcsz->Nbins_z;
       index_z++)
  {
    class_alloc(pcsz->dNdzdy_theoretical[index_z],
                (pcsz->Nbins_y+1)*sizeof(double),
                pcsz->error_message);
    class_alloc(pcsz->temp_0_theoretical[index_z],
                (pcsz->Nbins_y+1)*sizeof(double),
                pcsz->error_message);
    class_alloc(pcsz->temp_1_theoretical[index_z],
                (pcsz->Nbins_y+1)*sizeof(double),
                pcsz->error_message);
    class_alloc(pcsz->dNdzdm_theoretical[index_z],
                pcsz->nsteps_m*sizeof(double),
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
