# include "class_sz.h"
# include "class_sz_tools.h"
# include "Patterson.h"
# include "r8lib.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "fft.h"

/////////////////////////////////SZ-TOOLS//////////

int zbrent_sz_delta_to_delta_prime_nfw(double x1,
                                        double x2,
                                        double tol,
                                        double cvir,
                                        double cvir_prime,
                                        double delta,
                                        double fa,
                                        double fb,
                                        double * delta_prime,
                                       struct tszspectrum * ptsz){
  int iter;
  int ITMAX = 100;

  double a;
  double b;
  double c;
  double d;
  double e;
  double min1;
  double min2;
  double fc;
  double p;
  double q;
  double r;
  double tol1;
  double s;
  double xm;
  double EPS2;


  EPS2=3.e-8;
  a =x1;
  b =x2;

  class_call(
               dtod_prime_nfw(
                      a,
                      delta,
                      cvir,
                      cvir_prime,
                      &fa
                      ),
             ptsz->error_message,
             ptsz->error_message);

  class_call(
               dtod_prime_nfw(
                      b,
                      delta,
                      cvir,
                      cvir_prime,
                      &fb
                      ),
             ptsz->error_message,
             ptsz->error_message);


  if ((fb)*(fa) > 0.0)  {
    printf("Root must be bracketed in ZBRENT\n");
    return _FAILURE_;
  }

  fc=fb;

  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb)*(fc) > 0.0) {
      c=a;
      fc=fa;
      e=d=b-a;
    }

    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*(EPS2)*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0)  {
      *delta_prime = b;


      return _SUCCESS_;
    }

    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/(fa);
      if (a == c) {
        p=2.0*(xm)*(s);
        q=1.0-s;
      }
      else {
        q=fa/(fc);
        r=fb/(fc);
        p=s*(2.0*(xm)*(q)*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0)  q = -q;
      p=fabs(p);
      min1=3.0*(xm)*(q)-fabs(tol1*(q));
      min2=fabs(e*(q));
      if (2.0*(p) < (min1 < min2 ? min1 : min2))
      {
        e=d;
        d=p/(q);
      }
      else {
        d=xm;
        e=d;
      }
    }
    else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
    class_call(
               dtod_prime_nfw(
                      b,
                      delta,
                      cvir,
                      cvir_prime,
                      &fb
                      ),
               ptsz->error_message,
               ptsz->error_message);
  }

  printf("Max. num. of ite. exceeded in ZBRENT\n");

  return _FAILURE_;



                                        }

 int dtod_prime_nfw( double delta_prime,
                     double delta,
                     double cvir,
                     double cvir_prime,
                     double * dRES
                   ){


*dRES = m_nfw(delta_prime*cvir_prime)/m_nfw(cvir_prime) - m_nfw(delta*cvir)/m_nfw(cvir);

return _SUCCESS_;
                   }


//Routine used for
//the conversion between delta's
double delta_to_delta_prime_nfw(
  double delta,
  double cvir,
  double cvir_prime,
  struct tszspectrum * ptsz
              )
{
  double delta_prime;

  double  var;

  double  lTEST;

  double  fa;
  double  fb;
  double  m1;
  double  m2;
  double  mLO;
  double  mUP;
  double  logMDEL;



  int  i;
  int iMAX = 50;

  double * mTEST;
  class_alloc(mTEST,
              iMAX*sizeof( double ),
              ptsz->error_message);



  mTEST[0] = delta;


  class_call(
             dtod_prime_nfw(
                    mTEST[0],
                    delta,
                    cvir,
                    cvir_prime,
                    &lTEST
                    ),
             ptsz->error_message,
             ptsz->error_message
             );
//printf("lTEST = %.3e delta = %.3e\n", lTEST, delta);

  if (lTEST <= 0.) {
    for (i=1;i<iMAX;i++ ) {

      mTEST[i] = 2.*mTEST[i-1];

      class_call(
        dtod_prime_nfw(
                     mTEST[i],
                     delta,
                     cvir,
                     cvir_prime,
                     &lTEST
                     ),
                 ptsz->error_message,
                 ptsz->error_message
                 );

      if (lTEST > 0.)
      {
        m1 = mTEST[i];
        m2 = mTEST[i-1];
        break;
      }
    }
  }
  else
  {
    for (i=1;i<iMAX;i++ )
    {
      mTEST[i] = mTEST[i-1]/2.;

      class_call(
        dtod_prime_nfw(
               mTEST[i],
               delta,
               cvir,
               cvir_prime,
               &lTEST
               ),
                 ptsz->error_message,
                 ptsz->error_message);

    //printf("lTEST = %.3e i = %d\n", lTEST, i);

      if(lTEST < 0.)
      {
        m1 = mTEST[i];
        m2 = mTEST[i-1];
        break;
      }
    }
  }

  mLO=MIN(m1,m2);
  mUP=MAX(m1,m2);

  //printf("mLO = %.3e mUP = %.3e\n", mLO,mUP);

  class_call(zbrent_sz_delta_to_delta_prime_nfw(
                       mLO,
                       mUP,
                       1.e-4,
                       cvir,
                       cvir_prime,
                       delta,
                       fa,
                       fb,
                       &delta_prime,
                      ptsz),
             ptsz->error_message,
             ptsz->error_message);





  free(mTEST);
  return delta_prime;
}


//Root finding algorithm
//for the virial mass mVIR to
//overdensity mass mDEL (e.g. m200)
int zbrent_sz(
              double x1,
              double x2,
              double tol,
              double VAR1,
              double VAR2,
              double VAR3,
              double VAR4,
              double fa,
              double fb,
              double * logMDEL,
              struct tszspectrum * ptsz
              )
{
  int iter;
  int ITMAX = 100;

  double a;
  double b;
  double c;
  double d;
  double e;
  double min1;
  double min2;
  double fc;
  double p;
  double q;
  double r;
  double tol1;
  double s;
  double xm;
  double EPS2;


  EPS2=3.e-8;
  a =x1;
  b =x2;

  class_call(
             mVtomD(
                    a,
                    VAR1,
                    VAR2,
                    VAR3,
                    VAR4,
                    &fa,
                    ptsz
                    ),
             ptsz->error_message,
             ptsz->error_message);

  class_call(
             mVtomD(
                    b,
                    VAR1,
                    VAR2,
                    VAR3,
                    VAR4,
                    &fb,
                    ptsz
                    ),
             ptsz->error_message,
             ptsz->error_message);


  if ((fb)*(fa) > 0.0)  {
    printf("Root must be bracketed in ZBRENT\n");
    return _FAILURE_;
  }

  fc=fb;

  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb)*(fc) > 0.0) {
      c=a;
      fc=fa;
      e=d=b-a;
    }

    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*(EPS2)*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0)  {
      *logMDEL = b;


      return _SUCCESS_;
    }

    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/(fa);
      if (a == c) {
        p=2.0*(xm)*(s);
        q=1.0-s;
      }
      else {
        q=fa/(fc);
        r=fb/(fc);
        p=s*(2.0*(xm)*(q)*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0)  q = -q;
      p=fabs(p);
      min1=3.0*(xm)*(q)-fabs(tol1*(q));
      min2=fabs(e*(q));
      if (2.0*(p) < (min1 < min2 ? min1 : min2))
      {
        e=d;
        d=p/(q);
      }
      else {
        d=xm;
        e=d;
      }
    }
    else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
    class_call(
               mVtomD(
                      b,
                      VAR1,
                      VAR2,
                      VAR3,
                      VAR4,
                      &fb,
                      ptsz
                      ),
               ptsz->error_message,
               ptsz->error_message);
  }

  printf("Max. num. of ite. exceeded in ZBRENT\n");

  return _FAILURE_;
}

//Root finding algorithm
//for the virial mass mDEL to
//overdensity mass mVIR (e.g. m200)
int zbrent_D_to_V_sz(
              double x1,
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
              struct tszspectrum * ptsz,
              struct background * pba
              )
{
  int iter;
  int ITMAX = 100;

  double a;
  double b;
  double c;
  double d;
  double e;
  double min1;
  double min2;
  double fc;
  double p;
  double q;
  double r;
  double tol1;
  double s;
  double xm;
  double EPS2;

  double mvir_test;
  double cvir_test;
  double rvir_test;



  EPS2=3.e-8;
  a =x1;
  b =x2;


  mvir_test = exp(a);
  cvir_test = evaluate_cvir_of_mvir(mvir_test,z,ptsz,pba);
  rvir_test = evaluate_rvir_of_mvir(mvir_test,delc,rhoc,ptsz);



  class_call(
             mDtomV(
                    a,
                    mDEL,
                    rvir_test,
                    cvir_test,
                    delrho,
                    &fa,
                    ptsz
                    ),
             ptsz->error_message,
             ptsz->error_message);

  mvir_test = exp(b);
  cvir_test = evaluate_cvir_of_mvir(mvir_test,z,ptsz,pba);
  rvir_test = evaluate_rvir_of_mvir(mvir_test,delc,rhoc,ptsz);


  class_call(
             mDtomV(
                    b,
                    mDEL,
                    rvir_test,
                    cvir_test,
                    delrho,
                    &fb,
                    ptsz
                    ),
             ptsz->error_message,
             ptsz->error_message);


  if ((fb)*(fa) > 0.0)  {
    printf("Root must be bracketed in ZBRENT\n");
    return _FAILURE_;
  }

  fc=fb;

  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb)*(fc) > 0.0) {
      c=a;
      fc=fa;
      e=d=b-a;
    }

    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*(EPS2)*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0)  {
      *logMVIR = b;


      return _SUCCESS_;
    }

    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/(fa);
      if (a == c) {
        p=2.0*(xm)*(s);
        q=1.0-s;
      }
      else {
        q=fa/(fc);
        r=fb/(fc);
        p=s*(2.0*(xm)*(q)*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0)  q = -q;
      p=fabs(p);
      min1=3.0*(xm)*(q)-fabs(tol1*(q));
      min2=fabs(e*(q));
      if (2.0*(p) < (min1 < min2 ? min1 : min2))
      {
        e=d;
        d=p/(q);
      }
      else {
        d=xm;
        e=d;
      }
    }
    else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));

  mvir_test = exp(b);
  cvir_test = evaluate_cvir_of_mvir(mvir_test,z,ptsz,pba);
  rvir_test = evaluate_rvir_of_mvir(mvir_test,delc,rhoc,ptsz);

    class_call(
               mDtomV(
                      b,
                      mDEL,
                      rvir_test,
                      cvir_test,
                      delrho,
                      &fb,
                      ptsz
                      ),
               ptsz->error_message,
               ptsz->error_message);
  }

  printf("Max. num. of ite. exceeded in ZBRENT\n");

  return _FAILURE_;
}


//Root finding algorithm
//for the nonlinear scale
int zbrent_pkl_to_knl(
              double x1,
              double x2,
              double tol,
              double fa,
              double fb,
              double * knl,
              double z,
              struct tszspectrum * ptsz,
              struct background * pba,
              struct nonlinear * pnl,
              struct primordial * ppm
              )
{
  int iter;
  int ITMAX = 100;

  double a;
  double b;
  double c;
  double d;
  double e;
  double min1;
  double min2;
  double fc;
  double p;
  double q;
  double r;
  double tol1;
  double s;
  double xm;
  double EPS2;

  double knl_test;




  EPS2=3.e-8;
  a =x1;
  b =x2;


  // knl_test = exp(a);


  class_call(
             pkl_to_knl(
                        a,
                        &fa,
                        z,
                        ptsz,
                        pba,
                        pnl,
                        ppm
                        ),
             ptsz->error_message,
             ptsz->error_message);

  // knl_test = exp(b);


  class_call(
             pkl_to_knl(
                    b,
                    &fb,
                    z,
                    ptsz,
                    pba,
                    pnl,
                    ppm
                    ),
             ptsz->error_message,
             ptsz->error_message);


  if ((fb)*(fa) > 0.0)  {
    printf("Root must be bracketed in ZBRENT\n");
    return _FAILURE_;
  }

  fc=fb;

  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb)*(fc) > 0.0) {
      c=a;
      fc=fa;
      e=d=b-a;
    }

    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*(EPS2)*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0)  {
      *knl = b;


      return _SUCCESS_;
    }

    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/(fa);
      if (a == c) {
        p=2.0*(xm)*(s);
        q=1.0-s;
      }
      else {
        q=fa/(fc);
        r=fb/(fc);
        p=s*(2.0*(xm)*(q)*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0)  q = -q;
      p=fabs(p);
      min1=3.0*(xm)*(q)-fabs(tol1*(q));
      min2=fabs(e*(q));
      if (2.0*(p) < (min1 < min2 ? min1 : min2))
      {
        e=d;
        d=p/(q);
      }
      else {
        d=xm;
        e=d;
      }
    }
    else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));



    class_call(
               pkl_to_knl(
                      b,
                      &fb,
                      z,
                      ptsz,
                      pba,
                      pnl,
                      ppm
                      ),
               ptsz->error_message,
               ptsz->error_message);
  }

  printf("Max. num. of ite. exceeded in ZBRENT\n");

  return _FAILURE_;
}




//This routine reads the tabulated
//C-M relation Zhao2009,
//and stores the tabulated values.
//C(redshift[i],lnmass[j])

int read_Zhao_CM_init(
                      struct tszspectrum * ptsz
                      )
{
  //read the redshift and ln mass tables
  char line[_LINE_LENGTH_MAX_];
  FILE *process;
  int n_data_guess, n_data = 0;
  double *lnx = NULL, *tmp = NULL, **logC = NULL;
  double this_lnx;
  int status;
  int index_x;
  int index_redshift;
  int index_mass;


  class_alloc(ptsz->CM_redshift,sizeof(double *)*100,ptsz->error_message);
  class_alloc(ptsz->CM_logM,sizeof(double *)*100,ptsz->error_message);



  n_data = 0;
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));

  class_open(process,"sz_auxiliary_files/C-M_Zhao09/lnconcentration_vs_z_and_lnm-redshits.txt", "r",ptsz->error_message);

  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 ptsz->error_message,
                 "Error allocating memory to read the C-M relation Zhao et al 2009.\n");
      lnx = tmp;
    };


    /* Store */
    lnx[n_data]   = this_lnx;
    n_data++;
  }

  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  ptsz->CM_redshift_size = n_data;

  class_realloc(ptsz->CM_redshift,
                ptsz->CM_redshift,
                ptsz->CM_redshift_size*sizeof(double),
                ptsz->error_message);


  /** Store them */
  for (index_x=0; index_x<ptsz->CM_redshift_size; index_x++) {
    ptsz->CM_redshift[index_x] = lnx[index_x];
  };


  //Masses

  n_data = 0;
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));


  class_open(process,"sz_auxiliary_files/C-M_Zhao09/lnconcentration_vs_z_and_lnm-masses.txt", "r",ptsz->error_message);



  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 ptsz->error_message,
                 "Error allocating memory to read the C-M relation Zhao et al 2009.\n");
      lnx = tmp;
    };


    /* Store */
    lnx[n_data]   = this_lnx;
    n_data++;
  }

  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  ptsz->CM_logM_size = n_data;

  class_realloc(ptsz->CM_logM,
                ptsz->CM_logM,
                ptsz->CM_logM_size*sizeof(double),
                ptsz->error_message);


  /** Store them */
  for (index_x=0; index_x<ptsz->CM_logM_size; index_x++) {
    ptsz->CM_logM[index_x] = lnx[index_x];
  };


  /** Release the memory used locally */
  free(lnx);

  //Read concentration (lnC)

  class_alloc(ptsz->CM_logC,
              sizeof(double *)*ptsz->CM_redshift_size*ptsz->CM_logM_size,
              ptsz->error_message);

  class_alloc(logC,
              ptsz->CM_redshift_size*sizeof(double *),
              ptsz->error_message);


  for (index_redshift=0;
       index_redshift<ptsz->CM_redshift_size;
       index_redshift++)
  {
    class_alloc(logC[index_redshift],
                ptsz->CM_logM_size*sizeof(double),
                ptsz->error_message);
  }


  class_open(process,"sz_auxiliary_files/C-M_Zhao09/lnconcentration_vs_z_and_lnm.txt", "r",ptsz->error_message);



  int z =0;
  while (fgets(line, sizeof(line)-1, process) != NULL) {
    int i=0;
    char *err, *p = line;
    double val;
    while (*p) {
      val = strtod(p, &err);
      logC[z][i] = val;
      p = err + 1;
      i+=1;
    }
    z+=1;
  }


  int index = 0;
  for (index_redshift=0;
       index_redshift<ptsz->CM_redshift_size;
       index_redshift++){
    for (index_mass=0;
         index_mass<ptsz->CM_logM_size;
         index_mass++){

      ptsz->CM_logC[index] = logC[index_redshift][index_mass];
      index += 1;
    }
  }

  status = fclose(process);


  for (index_redshift=0;
       index_redshift<ptsz->CM_redshift_size;
       index_redshift++){
         free(logC[index_redshift]);
       }
  free(logC);

  return _SUCCESS_;
}

//Zhao et al 2009
//concentration mass relation
//cVIR-mVIR computed with mandc-1.03main
//for PL15 BF cosmo.
//Read tabulated values and interpolate
int  CvirMvirZHAO(
                  double * result,
                  double logM ,
                  double logz,
                  struct tszspectrum * ptsz
                  )
{


  double logz_asked = logz;
  double logM_asked = logM;

  if (logz<ptsz->CM_redshift[0])
    logz_asked = ptsz->CM_redshift[0];
  if (logz>ptsz->CM_redshift[ptsz->CM_redshift_size-1])
    logz_asked =  ptsz->CM_redshift[ptsz->CM_redshift_size-1];
  if (logM<ptsz->CM_logM[0])
    logM_asked = ptsz->CM_logM[0];
  if (logM>ptsz->CM_logM[ptsz->CM_logM_size-1])
    logM_asked =  ptsz->CM_logM[ptsz->CM_logM_size-1];

  *result = exp(pwl_interp_2d(
                              ptsz->CM_redshift_size,
                              ptsz->CM_logM_size,
                              ptsz->CM_redshift,
                              ptsz->CM_logM,
                              ptsz->CM_logC,
                              1,
                              &logz_asked,
                              &logM_asked
                              ));

  return _SUCCESS_;
}


//Sanchez-Conde & Prada 2014
//concentration mass relation
//c200-m200 crit
int  C200M200SC14(
                  double * result,
                  double logM ,
                  double z,
                  struct tszspectrum * ptsz
                  )
{
  double c_array[6] =
  {
    37.5153,
    -1.5093,
    1.636e-2,
    3.66e-4,
    -2.89237e-5,
    5.32e-7
  };
  *result =
  c_array[0]
  *pow(logM,0)
  *pow(1.+z,-1.)
  +
  c_array[1]
  *pow(logM,1)
  *pow(1.+z,-1.)
  +
  c_array[2]
  *pow(logM,2)
  *pow(1.+z,-1.)
  +
  c_array[3]
  *pow(logM,3)
  *pow(1.+z,-1.)
  +
  c_array[4]
  *pow(logM,4)
  *pow(1.+z,-1.)
  +
  c_array[5]
  *pow(logM,5)
  *pow(1.+z,-1.);

  return _SUCCESS_;
}

int  CvirMvirKLYPIN (
                     double * result,
                     double logM ,
                     double z,
                     struct tszspectrum * ptsz
                     )
{

  double z_tab[20] =
  {
    0.,
    0.31578947,
    0.63157895,
    0.94736842,
    1.26315789,
    1.57894737,
    1.89473684,
    2.21052632,
    2.52631579,
    2.84210526,
    3.15789474,
    3.47368421,
    3.78947368,
    4.10526316,
    4.42105263,
    4.73684211,
    5.05263158,
    5.36842105,
    5.68421053,
    6.
  };


  double c0_tab[20] =
  {
    9.6,
    7.89848895,
    6.57388797,
    5.59198421,
    4.82413741,
    4.2543651,
    3.80201899,
    3.4341066,
    3.15047911,
    2.92643281,
    2.74396076,
    2.60306296,
    2.50373941,
    2.44412709,
    2.40000661,
    2.36392585,
    2.33588481,
    2.31588348,
    2.30392188,
    2.3
  };

  double d2c0_tab[20] =
  {
    0.08062053,
    0.08062053,
    0.08062053,
    0.08062053,
    0.08062053,
    0.08062053,
    0.41689737,
    0.41689737,
    0.41689737,
    0.41689737,
    0.41689737,
    0.41689737,
    0.41689737,
    0.84668258,
    0.84668258,
    0.84668258,
    0.84668258,
    0.84668258,
    0.84668258,
    0.84668258
  };

  double lnM0_tab[20] =
  {
    45.46291469,
    41.51644832,
    38.29554435,
    35.80033605,
    34.03132449,
    32.96668373,
    32.12518764,
    31.309971,
    30.52126833,
    29.79801323,
    29.16858499,
    28.6329836,
    28.19120908,
    27.84158345,
    27.56229323,
    27.34662656,
    27.19458346,
    27.10616391,
    27.08136792,
    27.12019549
  };


  double d2lnM0_tab[20] =
  {
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291,
    0.63800291
  };

  double * c0z,* M0z;

  class_alloc(c0z,
              1*sizeof(double),
              ptsz->error_message);
  class_alloc(M0z,
              1*sizeof(double),
              ptsz->error_message);

  splint(z_tab,c0_tab,d2c0_tab,20,z,c0z);
  splint(z_tab,lnM0_tab,d2lnM0_tab,20,z,M0z);


  double Mvir = exp(logM);


  //Eq. 12 of 1002.3660v4
  *result =
  (*c0z)
  *pow(Mvir/1.e12,-0.075)
  *(1.+pow(Mvir/(exp(*M0z)),0.26));

  free(c0z);
  free(M0z);

  return _SUCCESS_;

}

double evaluate_rvir_of_mvir(double mvir,
                            double delc,
                            double rhoc,
                            struct tszspectrum * ptsz){

return pow(3.*mvir/(4*_PI_*delc*rhoc),1./3.);
                            }


double evaluate_cvir_of_mvir(double mvir,
                            double z,
                            struct tszspectrum * ptsz,
                            struct background * pba){
double cvir;
//D08 c-m relation
if (ptsz->concentration_parameter==0){
cvir = 7.85*pow(mvir/2.e12,-0.081)*pow(1.+z,-0.71);
}

//S00 c-m relation
else if (ptsz->concentration_parameter==1){
  cvir =10.*pow(mvir/3.42e12,-0.2)/(1.+z);
}

//K10 c-m relation
else if (ptsz->concentration_parameter==2){
   class_call(CvirMvirKLYPIN(&cvir,log(mvir),z,ptsz),
                   ptsz->error_message,
                   ptsz->error_message);
}


//SC14 c-m relation
// TBD m200c, r200c
else if (ptsz->concentration_parameter==3){
  printf("Warning: implementation of this concentration needs check.\n");
  exit(0);
   class_call(C200M200SC14(&cvir,
                           log(mvir),
                           z,
                           ptsz),
                   ptsz->error_message,
                   ptsz->error_message);
}


 //Z09 interpolated c-m relation
 else if (ptsz->concentration_parameter==4){
    printf("Warning: implementation of this concentration needs check.\n");
    exit(0);
    class_call(CvirMvirZHAO(&cvir,log(mvir),log(z),ptsz),
                    ptsz->error_message,
                    ptsz->error_message);
 }

// Dutton and Maccio 2014 (https://arxiv.org/pdf/1402.7073.pdf)
else if (ptsz->concentration_parameter==5){
  // here for virial mass in Msun/h:
  // see ea. 7 of 1402.7073
  double a =  0.537 + (1.025-0.537)*exp(-0.718*pow(z,1.08));
  double b = -0.097 + 0.024*z;
  double log10cvir = a + b*log10(mvir/1.e12);
  cvir = pow(10.,log10cvir);
}

// else if (ptsz->HMF==1 && ptsz->tau_profile == 1){
//
// }

else if (ptsz->concentration_parameter==6){

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
// double nu = 1./D*(1.12*pow(mvir/5e13,0.3)+0.53);

// # Compute the spherical collapse threshold of Nakamura-Suto, 1997.
// Om_mz = self.cosmology._Omega_m()
// dc0 = (3./20.)*pow(12.*np.pi,2./3.);
// self.delta_c = dc0*(1.+0.012299*np.log10(Om_mz));
// nu = delta_c / sig

  // double sig = get_sigma_at_z_and_m(z,mvir,ptsz,pba);

double nu = sqrt(get_nu_at_z_and_m(z,mvir,ptsz,pba));

cvir  = pow(D,0.9)*7.7*pow(nu,-0.29); // vir
free(pvecback);
}


return cvir;
                            }


//Routine used for
//the conversion between
//the viral mass and the overdensity mass
int mVtomD (
            double logMD ,
            double mVIR,
            double rvir,
            double c,
            double delrho,
            double * mRES,
            struct tszspectrum * ptsz
            )
{
  double  C;
  double rs = rvir/c;

  // here C is r_delta/r_s
  C = pow(3.*exp(logMD)/(4*_PI_*delrho),1./3.)/rs;

  *mRES =
  exp(logMD)/mVIR
  -(log(1.+C)
    -C/(1.+C))
  /(log(1.+c)
    -c/(1.+c));


  return _SUCCESS_;
}



//Routine used for
//finding the non-linear scale
int pkl_to_knl (
            double knl,
            double * mRES,
            double z,
            struct tszspectrum * ptsz,
            struct background * pba,
            struct nonlinear * pnl,
            struct primordial * ppm
            )
{
  double  knl_mpc,pkl_mpc;
  knl_mpc = knl*pba->h;

    //   double tau;
    //   int first_index_back = 0;
    //
    //
    //   class_call_(background_tau_of_z(pba,z,&tau),
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
    // free(pvecback);
    // free(pvectsz);

  enum pk_outputs pk_for_knl;
  pk_for_knl = pk_linear;
  double * pk_ic = NULL;
  double pk;
  double k;

  k = knl_mpc;
  // printf("knl=%.3e k=%.3e z=%.3e\n",knl,k,z);
    //Input: wavenumber in 1/Mpc
    //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
   class_call(nonlinear_pk_at_k_and_z(
                                     pba,
                                     ppm,
                                     pnl,
                                     pk_linear,
                                     k,
                                     z,
                                     pnl->index_pk_cb,
                                     &pk, // number *out_pk_l
                                     pk_ic // array out_pk_ic_l[index_ic_ic]
                                   ),
                                   pnl->error_message,
                                   pnl->error_message);

// printf("pk=%.3e\n",pk);

  pkl_mpc = pk;


  *mRES =
  pow(knl_mpc,3.)*pkl_mpc
  -2.*_PI_*_PI_;



  return _SUCCESS_;
}


//Routine used for
//the conversion between
//the viral mass and the overdensity mass
int mDtomV (
            double logMVIR ,
            double mD,
            double rvir,
            double c,
            double delrho,
            double * mRES,
            struct tszspectrum * ptsz
            )
{
  double  C;

  double mvir = exp(logMVIR);
  double rs = rvir/c;

  C = pow(3.*mD/(4*_PI_*delrho),1./3.)/rs;


    *mRES =
    mvir/mD
    -1./((log(1.+C)
      -C/(1.+C))
    /(log(1.+c)
      -c/(1.+c)));



  return _SUCCESS_;
}


 int mDEL_to_mDELprime(
               double mDEL ,
               double delrho,
               double delrho_prime,
               double delc,
               double rhoc,
               double z,
               double * mDELprime,
               struct tszspectrum * ptsz,
               struct background * pba
             ){

double mvir;
//first go from mDEL to mVIR:
class_call(mDEL_to_mVIR(mDEL,
                        delrho,
                        delc,
                        rhoc,
                        z,
                        &mvir,
                        ptsz,
                        pba),
                ptsz->error_message,
                ptsz->error_message);

//then go from mvir to mdel_prime
double rvir = evaluate_rvir_of_mvir(mvir,delc,rhoc,ptsz);
double cvir = evaluate_cvir_of_mvir(mvir,z,ptsz,pba);

class_call(mVIR_to_mDEL(mvir,
                     rvir,
                     cvir,
                     delrho_prime,
                     mDELprime,
                     ptsz),
                ptsz->error_message,
                ptsz->error_message);

return _SUCCESS_;
             }




//Routine used for
//the non linear scale
int solve_pkl_to_knl(
              double * result,
              double z,
              struct tszspectrum * ptsz,
              struct background * pba,
              struct nonlinear * pnl,
              struct primordial * ppm
              )
{

  double  mDEL;
  double  var;

  double  lTEST;

  double  fa;
  double  fb;
  double  m1;
  double  m2;
  double  mLO;
  double  mUP;
  double  logMDEL;



  int  i;
  int iMAX = 50;

  double * mTEST;
  class_alloc(mTEST,
              iMAX*sizeof( double ),
              ptsz->error_message);



  mTEST[0] = 1.;

 // printf("res 0 ini : %.3e\n",lTEST);
  class_call(
             pkl_to_knl(
                    mTEST[0],
                    &lTEST,
                    z,
                    ptsz,
                    pba,
                    pnl,
                    ppm
                    ),
             ptsz->error_message,
             ptsz->error_message
             );
 // printf("res 0 : %.3e\n",lTEST);
 //exit(0);

  if (lTEST <= 0.) {
    for (i=1;i<iMAX;i++ ) {

      mTEST[i] = 2.*mTEST[i-1];

      class_call(
                 pkl_to_knl(
                        mTEST[i],
                        &lTEST,
                        z,
                        ptsz,
                        pba,
                        pnl,
                        ppm
                        ),
                 ptsz->error_message,
                 ptsz->error_message
                 );

      if (lTEST > 0.)
      {
        m1 = mTEST[i];
        m2 = mTEST[i-1];
        break;
      }
    }
  }
  else
  {
    for (i=1;i<iMAX;i++ )
    {
      mTEST[i] = mTEST[i-1]/2.;

      class_call(
                 pkl_to_knl(
                        mTEST[i],
                        &lTEST,
                        z,
                        ptsz,
                        pba,
                        pnl,
                        ppm
                        ),
                 ptsz->error_message,
                 ptsz->error_message);

      if(lTEST < 0.)
      {
        m1 = mTEST[i];
        m2 = mTEST[i-1];
        break;
      }
    }
  }

  mLO=MIN(m1,m2);
  mUP=MAX(m1,m2);

  class_call(zbrent_pkl_to_knl(
                               mLO,
                               mUP,
                               1.e-4,
                               fa,
                               fb,
                               &logMDEL,
                               z,
                               ptsz,
                               pba,
                               pnl,
                               ppm
                               ),
             ptsz->error_message,
             ptsz->error_message);

  mDEL = logMDEL;
  *result = mDEL;


  free(mTEST);


  return _SUCCESS_;
}



//Routine used for
//the conversion between masses
int mVIR_to_mDEL(
              double mVIR ,
              double rvir,
              double c ,
              double delrho,
              double * result,
              struct tszspectrum * ptsz
              )
{
  double  mDEL;
  double  var;

  double  lTEST;

  double  fa;
  double  fb;
  double  m1;
  double  m2;
  double  mLO;
  double  mUP;
  double  logMDEL;



  int  i;
  int iMAX = 50;

  double * mTEST;
  class_alloc(mTEST,
              iMAX*sizeof( double ),
              ptsz->error_message);



  mTEST[0] = mVIR;


  class_call(
             mVtomD(
                    log(mTEST[0]),
                    mVIR,
                    rvir,
                    c,
                    delrho,
                    &lTEST,
                    ptsz
                    ),
             ptsz->error_message,
             ptsz->error_message
             );

  if (lTEST <= 0.) {
    for (i=1;i<iMAX;i++ ) {

      mTEST[i] = 2.*mTEST[i-1];

      class_call(
                 mVtomD(
                        log(mTEST[i]),
                        mVIR,
                        rvir,
                        c,
                        delrho,
                        &lTEST,
                        ptsz
                        ),
                 ptsz->error_message,
                 ptsz->error_message
                 );

      if (lTEST > 0.)
      {
        m1 = log(mTEST[i]);
        m2 = log(mTEST[i-1]);
        break;
      }
    }
  }
  else
  {
    for (i=1;i<iMAX;i++ )
    {
      mTEST[i] = mTEST[i-1]/2.;

      class_call(
                 mVtomD(
                        log(mTEST[i]),
                        mVIR,
                        rvir,
                        c,
                        delrho,
                        &lTEST,
                        ptsz
                        ),
                 ptsz->error_message,
                 ptsz->error_message);

      if(lTEST < 0.)
      {
        m1 = log(mTEST[i]);
        m2 = log(mTEST[i-1]);
        break;
      }
    }
  }

  mLO=MIN(m1,m2);
  mUP=MAX(m1,m2);

  class_call(zbrent_sz(
                       mLO,
                       mUP,
                       1.e-4,
                       mVIR,
                       rvir,
                       c,
                       delrho,
                       fa,
                       fb,
                       &logMDEL,
                       ptsz
                       ),
             ptsz->error_message,
             ptsz->error_message);

  mDEL = exp(logMDEL);
  *result = mDEL;


  free(mTEST);


  return _SUCCESS_;
}

//Routine used for
//the conversion between masses
int mDEL_to_mVIR(
              double mDEL ,
              double delrho,
              double delc,
              double rhoc,
              double z,
              double * result,
              struct tszspectrum * ptsz,
              struct background * pba
              )
{
  double  mVIR;
  double  * mTEST;
  double  lTEST;

  double  fa;
  double  fb;
  double  m1;
  double  m2;
  double  mLO;
  double  mUP;
  double  logMVIR;

  int  i;
  int iMAX = 50;


  class_alloc(mTEST,
              iMAX*sizeof( double ),
              ptsz->error_message);


  // var[0] = mDEL;
  // // var[1] = rs;
  // // var[2] = c;
  // var[1] = delrho;
  double mvir_test,cvir_test,rvir_test;

  mTEST[0] = mDEL;

  mvir_test = mTEST[0];
  cvir_test = evaluate_cvir_of_mvir(mvir_test,z,ptsz,pba);
  rvir_test = evaluate_rvir_of_mvir(mvir_test,delc,rhoc,ptsz);




  class_call(
             mDtomV(
                    log(mTEST[0]),
                    mDEL,
                    rvir_test,
                    cvir_test,
                    delrho,
                    &lTEST,
                    ptsz
                    ),
             ptsz->error_message,
             ptsz->error_message
             );

  if (lTEST <= 0.) {
    for (i=1;i<iMAX;i++ ) {

      mTEST[i] = 2.*mTEST[i-1];

      mvir_test = mTEST[i];
      cvir_test = evaluate_cvir_of_mvir(mvir_test,z,ptsz,pba);
      rvir_test = evaluate_rvir_of_mvir(mvir_test,delc,rhoc,ptsz);


      class_call(
                 mDtomV(
                        log(mTEST[i]),
                        mDEL,
                        rvir_test,
                        cvir_test,
                        delrho,
                        &lTEST,
                        ptsz
                        ),
                 ptsz->error_message,
                 ptsz->error_message
                 );

      if (lTEST > 0.)
      {
        m1 = log(mTEST[i]);
        m2 = log(mTEST[i-1]);
        break;
      }
    }
  }
  else
  {
    for (i=1;i<iMAX;i++ )
    {
      mTEST[i] = mTEST[i-1]/2.;

      mvir_test = mTEST[i];
      cvir_test = evaluate_cvir_of_mvir(mvir_test,z,ptsz,pba);
      rvir_test = evaluate_rvir_of_mvir(mvir_test,delc,rhoc,ptsz);

      class_call(
                 mDtomV(
                        log(mTEST[i]),
                        mDEL,
                        rvir_test,
                        cvir_test,
                        delrho,
                        &lTEST,
                        ptsz
                        ),
                 ptsz->error_message,
                 ptsz->error_message);

      if(lTEST < 0.)
      {
        m1 = log(mTEST[i]);
        m2 = log(mTEST[i-1]);
        break;
      }
    }
  }

  mLO=MIN(m1,m2);
  mUP=MAX(m1,m2);

  //printf("z= %.5e mLO = %.8e mUP = %.8e\n", z,mLO,mUP);

  class_call(zbrent_D_to_V_sz(
                       mLO,
                       mUP,
                       1.e-4,
                       mDEL,
                       delrho,
                       fa,
                       fb,
                       z,
                       delc,
                       rhoc,
                       &logMVIR,
                       ptsz,
                       pba
                       ),
             ptsz->error_message,
             ptsz->error_message);

  mVIR = exp(logMVIR);
  *result = mVIR;


  free(mTEST);

  return _SUCCESS_;
}


double m_nfw(double x){
return log(1.+x)-x/(1.+x);
}




struct Parameters_for_integrand_sigma2_hsv{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double z;
};



double integrand_sigma2_hsv(double lnk, void *p){

  struct Parameters_for_integrand_sigma2_hsv *V = ((struct Parameters_for_integrand_sigma2_hsv *) p);

  double pk;
  double k = exp(lnk);

  double * pk_ic = NULL;



  double W;


  //background quantities @ z:
  double tau;
  int first_index_back = 0;
  double * pvecback;
  class_alloc(pvecback,
              V->pba->bg_size*sizeof(double),
              V->pba->error_message);

  double z;
  if (V->z == 0.)
    z = V->z + 1e-10; // distance diverges at low-z
  else
    z = V->z;
  class_call(background_tau_of_z(V->pba,z,&tau),
             V->pba->error_message,
             V->pba->error_message);

  class_call(background_at_tau(V->pba,
                               tau,
                               V->pba->long_info,
                               V->pba->inter_normal,
                               &first_index_back,
                               pvecback),
             V->pba->error_message,
             V->pba->error_message);


  double Theta_s = sqrt(V->ptsz->Omega_survey/_PI_); // see below Eq. 45 of Takada and Spergel 2013
  double Chi = pvecback[V->pba->index_bg_ang_distance]*(1.+V->z);  //'Chi' comoving distance in Mpc
  double r_hsv = Chi*Theta_s; // in Mpc


  free(pvecback);


  //here k in 1/Mpc
  double x_hsv = k*r_hsv;

  W = 2.*gsl_sf_bessel_J1(x_hsv)/x_hsv; //see e.g., below Eq. 45 of Takada and Spergel 2013


    //Input: wavenumber in 1/Mpc
    //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
   class_call(nonlinear_pk_at_k_and_z(
                                     V->pba,
                                     V->ppm,
                                     V->pnl,
                                     pk_linear,
                                     k,
                                     V->z,
                                     V->pnl->index_pk_cb,
                                     &pk, // number *out_pk_l
                                     pk_ic // array out_pk_ic_l[index_ic_ic]
                                   ),
                                   V->pnl->error_message,
                                   V->pnl->error_message);
  double result = k*k*pk*W*W;




  return result;

}

int spectra_sigma2_hsv(
                   struct background * pba,
                   struct primordial * ppm,
                   struct nonlinear *pnl,
                   struct tszspectrum * ptsz,
                   double z,
                   double * sigma2_hsv
                   ) {

double k_min = exp(ptsz->ln_k_for_tSZ[0]);
double k_max = exp(ptsz->ln_k_for_tSZ[ptsz->ln_k_size_for_tSZ-1]);


struct Parameters_for_integrand_sigma2_hsv V;
  V.pnl = pnl;
  V.ppm = ppm;
  V.ptsz = ptsz;
  V.pba = pba;
  V.z = z;

  void * params = &V;
  double r; //result of the integral

  double epsrel = 1e-6;
  double epsabs = 1e-30;
  //int show_neval = ptsz->patterson_show_neval;

  r=Integrate_using_Patterson_adaptive(log(k_min),
                                        log(k_max),
                                        epsrel, epsabs,
                                        integrand_sigma2_hsv,
                                        params,0);



  //
  // gsl_function F;
  // F.function = &integrand_sigma2_hsv;
  // F.params = params;
  //
  // int n_subintervals_gsl = 300;
  //
  // // double epsrel=ptsz->mass_epsrel;
  // // double epsabs=ptsz->mass_epsabs;
  //
  //
  // gsl_integration_workspace * w = gsl_integration_workspace_alloc (n_subintervals_gsl);
  //
  // double result_gsl, error;
  // int key = 4;
  // gsl_integration_qag(&F,log(k_min),log(k_max),epsabs,epsrel,n_subintervals_gsl,key,w,&result_gsl,&error);
  // gsl_integration_workspace_free(w);
  //
  // r = result_gsl;

  *sigma2_hsv = r/(2.*_PI_)*pba->h;


                     }



struct Parameters_for_integrand_lensmag{
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
};



double integrand_lensmag(double ln1pzs, void *p){

  struct Parameters_for_integrand_lensmag *V = ((struct Parameters_for_integrand_lensmag *) p);

  double integrand;
  double  zs = exp(ln1pzs)-1.;




  double W;


  //background quantities @ zs:
  double tau;
  int first_index_back = 0;
  double * pvecback;
  class_alloc(pvecback,
              V->pba->bg_size*sizeof(double),
              V->pba->error_message);

  class_call(background_tau_of_z(V->pba,zs,&tau),
             V->pba->error_message,
             V->pba->error_message);

  class_call(background_at_tau(V->pba,
                               tau,
                               V->pba->long_info,
                               V->pba->inter_normal,
                               &first_index_back,
                               pvecback),
             V->pba->error_message,
             V->pba->error_message);



  double Chi_at_zs = pvecback[V->pba->index_bg_ang_distance]*(1.+zs);  //'Chi' comoving distance in Mpc
  double Chi_at_z = sqrt(V->pvectsz[V->ptsz->index_chi2])/V->pba->h;  //'Chi' comoving distance in Mpc



  free(pvecback);



  W = (Chi_at_zs-Chi_at_z)/Chi_at_zs;

  double dndzs = 0.;

/////////////////////////////////
  double z_asked  = zs;
  double phig = 0.;
//unwise: use Cosmos cross-match dndz
if (V->ptsz->galaxy_sample==1){
if(z_asked<V->ptsz->normalized_cosmos_dndz_z[0])
   phig = 1e-100;
else if (z_asked>V->ptsz->normalized_cosmos_dndz_z[V->ptsz->normalized_cosmos_dndz_size-1])
   phig = 1e-100;
else  phig =  pwl_value_1d(V->ptsz->normalized_cosmos_dndz_size,
                             V->ptsz->normalized_cosmos_dndz_z,
                             V->ptsz->normalized_cosmos_dndz_phig,
                             z_asked);
// printf("integrand ok phig\n");
}
else{
  if(z_asked<V->ptsz->normalized_dndz_z[0])
     phig = 1e-100;
  else if (z_asked>V->ptsz->normalized_dndz_z[V->ptsz->normalized_dndz_size-1])
     phig = 1e-100;
  else  phig =  pwl_value_1d(V->ptsz->normalized_dndz_size,
                               V->ptsz->normalized_dndz_z,
                               V->ptsz->normalized_dndz_phig,
                               z_asked);

}
 dndzs = phig;
////////////////////////////////

  integrand = dndzs*W;

  integrand *= (1.+zs);

  // printf("-> integrand z = %.3e phig = =%.3e\n",z_asked,phig);

  return integrand;

}

int redshift_int_lensmag(
                  struct tszspectrum * ptsz,
                  struct background * pba,
                  double * pvectsz,
                  double * result
                   ) {

double z =  pvectsz[ptsz->index_z];

double zs_min = z;
double zs_max = ptsz->z2SZ;


struct Parameters_for_integrand_lensmag V;
  V.pvectsz = pvectsz;
  V.ptsz = ptsz;
  V.pba = pba;


  void * params = &V;
  double r; //result of the integral

  double epsrel = 1e-6;
  double epsabs = 1e-30;
  //int show_neval = ptsz->patterson_show_neval;

  r=Integrate_using_Patterson_adaptive(log(1.+zs_min),
                                        log(1.+zs_max),
                                        epsrel, epsabs,
                                        integrand_lensmag,
                                        params,0);


  *result = r;
                     }



// velocity dispersion for kSZ quantities
// see e.g., eq. 29 of 1807.07324
// also Appendix B of 1711.07879 for a different approach

int spectra_vrms2(
                   struct background * pba,
                   struct primordial * ppm,
                   struct nonlinear *pnl,
                   struct tszspectrum * ptsz,
                   double z,
                   double * vrms2
                   //double * sigma_prime
                   ) {

  double pk;
  double * pk_ic = NULL;
  //double * tk = NULL; //transfer

  double * array_for_sigma;
  //double tk_cdm,tk_b,tk_m,tk_ncdm,Omega_cdm,Omega_b,Omega_ncdm;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double k,W;

  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              ptsz->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);

    //background quantities @ z:
    double tau;
    int first_index_back = 0;
    double * pvecback;
    class_alloc(pvecback,
                pba->bg_size*sizeof(double),
                pba->error_message);

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
    double aH = pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]; //in Mpc^-1
    aH *= _c_/1e5*1e2; //in km/s/Mpc

    W = f*aH ;
    // printf("ok z = %e\n",W);

    free(pvecback);



      for (i=0;i<ptsz->ln_k_size_for_tSZ;i++) {
        k=exp(ptsz->ln_k_for_tSZ[i]);
        if (i == (ptsz->ln_k_size_for_tSZ-1)) k *= 0.9999999;
// printf("ok k = %e I = %e\n",k,pk*W*W);
    // //Input: wavenumber in 1/Mpc
    // //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
  enum pk_outputs pk_for_vrms2;
  if (ptsz->pk_nonlinear_for_vrms2 == 1){
    pk_for_vrms2 = pk_nonlinear;
  }
  else {
    pk_for_vrms2 = pk_linear;
  }
  // printf("ok k2 = %e I = %e\n",k,pk*W*W);

   class_call(nonlinear_pk_at_k_and_z(
                                     pba,
                                     ppm,
                                     pnl,
                                     pk_for_vrms2,
                                     k,
                                     z,
                                     pnl->index_pk_cb,
                                     &pk, // number *out_pk_l
                                     pk_ic // array out_pk_ic_l[index_ic_ic]
                                   ),
                                   pnl->error_message,
                                   pnl->error_message);

  // printf("ok k3 = %e I = %e\n",k,pk*W*W);


    array_for_sigma[i*index_num+index_k]=k;
    array_for_sigma[i*index_num+index_y]=pk*W*W;
    // printf("ok k = %e I = %e\n",k,pk*W*W);
  }
// printf("ok z = %e\n",W);
  class_call(array_spline(array_for_sigma,
                          index_num,
                          ptsz->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);
//printf("ok z = %e\n",W);
  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        ptsz->ln_k_size_for_tSZ,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        vrms2,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);
//printf("ok z = %e\n",W);
  free(array_for_sigma);
  *vrms2 = *vrms2/(2.*_PI_*_PI_);
// printf("ok z = %e\n",*vrms2);
  return _SUCCESS_;

}


/**
 * This routine computes sigma(R) given P(k) and ncdm species (does not check that k_max is large
 * enough)
 *
 * @param pba   Input: pointer to background structure
 * @param ppm   Input: pointer to primordial structure
 * @param pnl   Input: pointer to spectra structure
 * @param z     Input: redshift
 * @param R     Input: radius in Mpc
 * @param sigma Output: variance in a sphere of radius R (dimensionless)
 */

int spectra_sigma_ncdm(
                       struct background * pba,
                       struct primordial * ppm,
                       struct nonlinear *pnl,
                       struct tszspectrum * ptsz,
                       double R,
                       double z,
                       double * sigma
                       ) {

    // printf("entering \n");

  double pk;
  double * pk_ic = NULL;
  //double * tk = NULL; //transfer

  double * array_for_sigma;
  //double tk_cdm,tk_b,tk_m,tk_ncdm,Omega_cdm,Omega_b,Omega_ncdm;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double k,W,x;





  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              ptsz->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);


  // printf("entering2 \n");

  for (i=0;i<ptsz->ln_k_size_for_tSZ;i++) {
    k=exp(ptsz->ln_k_for_tSZ[i]);
    if (i == (ptsz->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
    W=3./x/x/x*(sin(x)-x*cos(x));

  // printf("entering4 %.3e \n",R);

   class_call(nonlinear_pk_at_k_and_z(
                                     pba,
                                     ppm,
                                     pnl,
                                     pk_linear,
                                     k, //Input: wavenumber in 1/Mpc
                                     z,
                                     pnl->index_pk_cb,
                                     &pk, // number *out_pk_l
                                     pk_ic // array out_pk_ic_l[index_ic_ic]
                                   ),
                                   pnl->error_message,
                                   pnl->error_message);
  // printf("entering5 jsjsjsjs %d \n",i);
  //
  //   printf("pk sig =%.3e\n",pk);
  //   exit(0);

    array_for_sigma[i*index_num+index_k]=k;
    array_for_sigma[i*index_num+index_y]=k*k*pk*W*W;
  }

  class_call(array_spline(array_for_sigma,
                          index_num,
                          ptsz->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        ptsz->ln_k_size_for_tSZ,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  free(array_for_sigma);
  *sigma = sqrt(*sigma/(2.*_PI_*_PI_));

  return _SUCCESS_;

}





//This routine computes dSigma2/dR
//at R and z for ncdm species

int spectra_sigma_ncdm_prime(
                             struct background * pba,
                             struct primordial * ppm,
                             struct nonlinear *pnl,
                             struct tszspectrum * ptsz,
                             double R,
                             double z,
                             double * sigma_prime
                             ) {

  double pk;
  double * pk_ic = NULL;
  //double * tk = NULL; //transfer
  //double tk_cdm,tk_b,tk_m,tk_ncdm,Omega_cdm,Omega_b,Omega_ncdm;
  //double Omega_cdm,Omega_b,Omega_ncdm;
  double * array_for_sigma;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double k,W,x,W_prime;




  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              ptsz->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);

  for (i=0;i<ptsz->ln_k_size_for_tSZ;i++) {
    k=exp(ptsz->ln_k_for_tSZ[i]);
    if (i == (ptsz->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
    W=3./x/x/x*(sin(x)-x*cos(x));
    W_prime=3./x/x*sin(x)-9./x/x/x/x*(sin(x)-x*cos(x));


class_call(nonlinear_pk_at_k_and_z(
                                  pba,
                                  ppm,
                                  pnl,
                                  pk_linear,
                                  k,
                                  z,
                                  pnl->index_pk_cb,
                                  &pk, // number *out_pk_l
                                  pk_ic // array out_pk_ic_l[index_ic_ic]
                                ),
                                pnl->error_message,
                                pnl->error_message);



    array_for_sigma[i*index_num+index_k]=k;
    array_for_sigma[i*index_num+index_y]=k*k*pk*k*2.*W*W_prime;
  }

  class_call(array_spline(array_for_sigma,
                          index_num,
                          ptsz->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        ptsz->ln_k_size_for_tSZ,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma_prime,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  free(array_for_sigma);


  *sigma_prime = *sigma_prime/(2.*_PI_*_PI_);

  return _SUCCESS_;

}


//Spline interpolation routine
//for interpolating T08 HMF at m500
//and Klypin 2010 c-m relation
int splint(
           double xa[],
           double ya[],
           double y2a[],
           int npoints,
           double x,
           double *y
           )
{
  int klo,khi,k;
  float h,b,a;

  klo=0;
  khi = npoints-1;
  while (khi-klo > 1)
  {
    k = (khi+klo) >> 1;
    if (xa[k] > x)
      khi = k;
    else
      klo = k;
  }

  h = xa[khi] - xa[klo];
  if (h == 0.0) return 0; /* bad input */

  a = (xa[khi] - x)/h;
  b = (x-xa[klo])/h;

  *y =
  a*ya[klo]
  +b*ya[khi]
  +((a*a*a-a)*y2a[klo]
    +(b*b*b-b)*y2a[khi])
  *(h*h)
  /6.0;

  return 1;
}

//spline integration for mass integral
//
//
// int integrate_over_m_at_z_spline(struct tszspectrum * ptsz,
//                                  struct background * pba,
//                                  double * pvectsz,
//                                  double * result) {
//
//
//   double * array_for_integral;
//   int index_num;
//   int index_x;
//   int index_y;
//   int index_ddy;
//   int i;
//
//   double x,W;
//
//
//   i=0;
//   index_x=i;
//   i++;
//   index_y=i;
//   i++;
//   index_ddy=i;
//   i++;
//   index_num=i;
//
//   double integrand_value = 0.;
//
//   class_alloc(array_for_integral,
//               ptsz->ln_M_size*index_num*sizeof(double),
//               ptsz->error_message);
//
//   for (i=0;i<ptsz->ln_M_size;i++) {
//     x=exp(ptsz->ln_x_for_pp[i]);
//
//     plc_gnfw(&plc_gnfw_at_x,x,pvectsz,pba,ptsz);
//
//
//
//     double pp_at_x_and_ell_over_ell_char = x*plc_gnfw_at_x;
//     array_for_integral[i*index_num+index_x]= log(x);
//     array_for_integral[i*index_num+index_y]= pp_at_x_and_ell_over_ell_char;
//   }
//
//
//
//
//   class_call(array_spline(array_for_integral,
//                           index_num,
//                           ptsz->ln_x_size_for_pp,
//                           index_x,
//                           index_y,
//                           index_ddy,
//                           _SPLINE_EST_DERIV_,
//                           ptsz->error_message),
//              ptsz->error_message,
//              ptsz->error_message);
//
//   class_call(array_integrate_all_spline(array_for_integral,
//                                         index_num,
//                                         ptsz->ln_x_size_for_pp,
//                                         index_x,
//                                         index_y,
//                                         index_ddy,
//                                         result,
//                                         ptsz->error_message),
//              ptsz->error_message,
//              ptsz->error_message);
//
//   free(array_for_integral);
//
//   return _SUCCESS_;
//
// }
//


struct Parameters_for_integrand_patterson_pp{
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
};


struct Parameters_for_integrand_nfw_profile{
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
};


double integrand_nfw_profile(double x, void *p){

  struct Parameters_for_integrand_nfw_profile *V = ((struct Parameters_for_integrand_nfw_profile *) p);

    double nfw_profile_at_x = 0.;
    rho_gnfw(&nfw_profile_at_x,x,V->pvectsz,V->pba,V->ptsz);
    double result = nfw_profile_at_x;
  return result;

}

double integrand_patterson_test_pp(double x, void *p){

  struct Parameters_for_integrand_patterson_pp *V = ((struct Parameters_for_integrand_patterson_pp *) p);

    //double x=exp(ln_x);

    double plc_gnfw_at_x = 0.;
    plc_gnfw(&plc_gnfw_at_x,x,V->pvectsz,V->pba,V->ptsz);

    double result = plc_gnfw_at_x;

  return result;

}

int two_dim_ft_nfw_profile(struct tszspectrum * ptsz,
                          struct background * pba,
                          double * pvectsz,
                          double * result
                          ) {

  struct Parameters_for_integrand_nfw_profile V;
  V.ptsz = ptsz;
  V.pba = pba;
  V.pvectsz = pvectsz;


  void * params = &V;

  gsl_function F;
  F.function = &integrand_nfw_profile;
  F.params = params;

  double eps_abs = ptsz->nfw_profile_epsabs;
  double eps_rel = ptsz->nfw_profile_epsrel;

  double result_gsl, error;

  double xin = 1.e-5;
  double c_nfw;

  //Battaglia 16 case:

  double rvir = pvectsz[ptsz->index_rVIR]; //in Mpc/h
  // double r200c = pvectsz[ptsz->index_r200c]; //in Mpc/h
  double rs = pvectsz[ptsz->index_rs]; //in Mpc/h
  // xout = 50.*rvir/rs; // as in hmvec (default 20, but set to 50 in example file)
  double xout = 50.; // as in hmvec (default 20, but set to 50 in example file)

  c_nfw = 1.;


// QAWO

  double delta_l = xout - xin;

  gsl_integration_workspace * w;
  gsl_integration_qawo_table * wf;

  int size_w = 500;
  w = gsl_integration_workspace_alloc(size_w);


  double w0;

  int index_md = (int) pvectsz[ptsz->index_md];
  double y_eff;
  y_eff = (pvectsz[ptsz->index_multipole_for_nfw_profile]+0.5)
             /pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile];

  w0 = y_eff;


  wf = gsl_integration_qawo_table_alloc(w0, delta_l,GSL_INTEG_SINE,100); // default 30


  int limit = size_w; //number of sub interval
  gsl_integration_qawo(&F,xin,eps_abs,eps_rel,limit,w,wf,&result_gsl,&error);

  *result = result_gsl;

  gsl_integration_qawo_table_free(wf);
  gsl_integration_workspace_free(w);



}



int rho_gnfw(double * rho_nfw_x,
            double x ,
            double * pvectsz,
            struct background * pba,
            struct tszspectrum * ptsz)
{

 int index_md = (int) pvectsz[ptsz->index_md];

 double z = pvectsz[ptsz->index_z];

 double y_eff;
   y_eff = (pvectsz[ptsz->index_multipole_for_nfw_profile]+0.5)
            /pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile];


    double A_rho0 = ptsz->A_rho0;
    double A_alpha = ptsz->A_alpha;
    double A_beta = ptsz->A_beta;

    double alpha_m_rho0 = ptsz->alpha_m_rho0;
    double alpha_m_alpha = ptsz->alpha_m_alpha;
    double alpha_m_beta = ptsz->alpha_m_beta;

    double alpha_z_rho0 = ptsz->alpha_z_rho0;
    double alpha_z_alpha = ptsz->alpha_z_alpha;
    double alpha_z_beta = ptsz->alpha_z_beta;

  // Eq. A1 and A2:
  double m200_over_msol = pvectsz[ptsz->index_m200c]/pba->h; // convert to Msun
  double rho0 = 1.;
  double alpha = A_alpha*pow(m200_over_msol/1e14,alpha_m_alpha)*pow(1.+z,alpha_z_alpha);
  double beta = A_beta*pow(m200_over_msol/1e14,alpha_m_beta)*pow(1.+z,alpha_z_beta);

  double gamma = ptsz->gamma_B16;
  double xc = 0.5;

  *rho_nfw_x = get_gas_profile_at_x_M_z_b16_200c(x,
                                                 pvectsz[ptsz->index_m200c],
                                                 z,
                                                 A_rho0,
                                                 A_alpha,
                                                 A_beta,
                                                 alpha_m_rho0,
                                                 alpha_m_alpha,
                                                 alpha_m_beta,
                                                 alpha_z_rho0,
                                                 alpha_z_alpha,
                                                 alpha_z_beta,
                                                 gamma,
                                                 pba,
                                                 ptsz)*pow(x,2)/(x*y_eff);

}






double get_rho_crit_at_z(double z_asked,
                         struct background * pba,
                         struct tszspectrum * ptsz){
double result;
double rho_crit;




double * pvecback;
double * pvectsz;



double tau;
double z = z_asked;
int first_index_back = 0;

class_alloc(pvecback,
            pba->bg_size*sizeof(double),
            ptsz->error_message);

class_alloc(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
 int i;
 for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

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




pvectsz[ptsz->index_z] = z;
pvectsz[ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);

rho_crit = pvectsz[ptsz->index_Rho_crit];


free(pvecback);
free(pvectsz);



result = rho_crit;

return result;
}


double get_gas_profile_at_x_M_z_nfw_200m(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct tszspectrum * ptsz){
double result;
double r_asked = 0.;
double rho_s;
double delta;
double c_delta;
double r_delta;
double r_s;
double rho_crit;
double f_b = pba->Omega0_b/ptsz->Omega_m_0;
double x;
double p_x;



double * pvecback;
double * pvectsz;



double tau;
double z = z_asked;
int first_index_back = 0;

class_alloc(pvecback,
            pba->bg_size*sizeof(double),
            ptsz->error_message);

class_alloc(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
 int i;
 for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

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




pvectsz[ptsz->index_z] = z;
pvectsz[ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);

rho_crit = pvectsz[ptsz->index_Rho_crit];
delta = 200.*pvecback[pba->index_bg_Omega_m];
c_delta = get_c200m_at_m_and_z(m_asked,z,pba,ptsz);
r_delta = pow(3.*m_asked/(4.*_PI_*delta*rho_crit),1./3.); //in units of h^-1 Mpc

// rho_s = pow(c_delta,3.)*delta*rho_crit/3./m_nfw(c_delta);



r_s = r_delta/c_delta;
x = r_asked/r_s;
x = x_asked;
p_x = 1./x*1./pow(1.+x,2);
rho_s = m_asked/m_nfw(c_delta)/4./_PI_/pow(r_s,3.);

free(pvecback);
free(pvectsz);



// result = rho_s*f_b*p_x/rho_crit/f_b;
result = rho_s*p_x*f_b;
return result;
}



double get_gas_profile_at_x_M_z_b16_200c(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         double A_rho0,
                                         double A_alpha,
                                         double A_beta,
                                         double alpha_m_rho0,
                                         double alpha_m_alpha,
                                         double alpha_m_beta,
                                         double alpha_z_rho0,
                                         double alpha_z_alpha,
                                         double alpha_z_beta,
                                         double gamma,
                                         struct background * pba,
                                         struct tszspectrum * ptsz){
double result;
double r_asked = 0.;
double rho_s;
double delta;
double c_delta;
double r_delta;
double r_s;
double rho_crit;
double f_b = pba->Omega0_b/ptsz->Omega_m_0;
double x;
double p_x;



double * pvecback;
double * pvectsz;



double tau;
double z = z_asked;
int first_index_back = 0;

class_alloc(pvecback,
            pba->bg_size*sizeof(double),
            ptsz->error_message);

class_alloc(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
 int i;
 for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

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




pvectsz[ptsz->index_z] = z;
pvectsz[ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);

rho_crit = pvectsz[ptsz->index_Rho_crit];
delta = 200.;
//pvectsz[ptsz->index_m200c] = m_asked;

c_delta = 1.;
r_delta = pow(3.*m_asked/(4.*_PI_*delta*rho_crit),1./3.); //in units of h^-1 Mpc

// rho_s = pow(c_delta,3.)*delta*rho_crit/3./m_nfw(c_delta);

r_s = r_delta/c_delta;
x = r_asked/r_s;
x = x_asked;
// p_x = 1./x*1./pow(1.+x,2);

free(pvecback);
free(pvectsz);


    // double A_rho0;
    // double A_alpha;
    // double A_beta;
    //
    // double alpha_m_rho0;
    // double alpha_m_alpha;
    // double alpha_m_beta;
    //
    // double alpha_z_rho0;
    // double alpha_z_alpha;
    // double alpha_z_beta;

  // // Battaglia 16 -- https://arxiv.org/pdf/1607.02442.pdf
  // // Table 2
  // if (ptsz->tau_profile_mode == 0){
  //   // agn feedback
  //   // A_rho0 = 4.e3;
  //   A_alpha = 0.88;
  //   A_beta = 3.83;
  //
  //   alpha_m_rho0 = 0.29;
  //   alpha_m_alpha = -0.03;
  //   alpha_m_beta = 0.04;
  //
  //   alpha_z_rho0 = -0.66;
  //   alpha_z_alpha = 0.19;
  //   alpha_z_beta = -0.025;
  //   }
  // else if (ptsz->tau_profile_mode == 1){
  //   // shock heating
  //   // A_rho0 = 1.9e4;
  //   A_alpha = 0.70;
  //   A_beta = 4.43;
  //
  //   alpha_m_rho0 = 0.09;
  //   alpha_m_alpha = -0.017;
  //   alpha_m_beta = 0.005;
  //
  //   alpha_z_rho0 = -0.95;
  //   alpha_z_alpha = 0.27;
  //   alpha_z_beta = 0.037;
  // }

  // Eq. A1 and A2:
  double m200_over_msol = m_asked/pba->h; // convert to Msun
  // double rho0  = 1.;
  double rho0 = A_rho0*pow(m200_over_msol/1e14,alpha_m_rho0)*pow(1.+z,alpha_z_rho0);
  double alpha = A_alpha*pow(m200_over_msol/1e14,alpha_m_alpha)*pow(1.+z,alpha_z_alpha);
  double beta = A_beta*pow(m200_over_msol/1e14,alpha_m_beta)*pow(1.+z,alpha_z_beta);

  // double gamma = -0.2;
  double xc = 0.5;

  p_x = pow(x/xc,gamma)*pow(1.+ pow(x/xc,alpha),-(beta+gamma)/alpha);
  // p_x = m200_over_msol;

  //
  // double rho0 = A_rho0*pow(m200_over_msol/1e14,alpha_m_rho0)*pow(1.+z,alpha_z_rho0);
  // tau_normalisation = rho0*4.*_PI_*pow(pvectsz[ptsz->index_r200c],3)
  //                     *pvectsz[ptsz->index_Rho_crit];


// result = rho0*rho_crit*f_b*p_x/rho_crit/f_b;
result = rho0*rho_crit*p_x*f_b;
return result;
}




double get_gas_density_profile_at_l_M_z(double l_asked, double m_asked, double z_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
  double l = log(l_asked);


  // if (ptsz->tau_profile == 1){
  // find the closest l's in the grid:
  int id_l_low;
  int id_l_up;
  int n_ell = ptsz->n_ell_density_profile;
  int n_m = ptsz->n_m_density_profile;
  int n_z = ptsz->n_z_density_profile;
  r8vec_bracket(n_ell,ptsz->array_profile_ln_l,l,&id_l_low,&id_l_up);

  // interpolate 2d at l_low:

 double ln_rho_low = pwl_interp_2d(n_m,
                                n_z,
                                ptsz->array_profile_ln_m,
                                ptsz->array_profile_ln_1pz,
                                ptsz->array_profile_ln_rho_at_lnl_lnM_z[id_l_low-1],
                                1,
                                &m,
                                &z);

 double ln_rho_up = pwl_interp_2d(n_m,
                                n_z,
                                ptsz->array_profile_ln_m,
                                ptsz->array_profile_ln_1pz,
                                ptsz->array_profile_ln_rho_at_lnl_lnM_z[id_l_up-1],
                                1,
                                &m,
                                &z);
 double ln_l_low = ptsz->array_profile_ln_l[id_l_low-1];
 double ln_l_up = ptsz->array_profile_ln_l[id_l_up-1];


 return ln_rho_low + ((l - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low);


}




// Tabulate 2D Fourier transform of density profile on a [z - ln_M - ln_ell] grid
// this is the tau profile for kSZ
// here we adopt m200c !
int tabulate_gas_density_profile(struct background * pba,
                             struct tszspectrum * ptsz){

if (ptsz->has_kSZ_kSZ_lensmag_1halo
+ ptsz->has_kSZ_kSZ_gal_1h
+ ptsz->has_kSZ_kSZ_gal_1h_fft
+ ptsz->has_kSZ_kSZ_gal_2h_fft
+ ptsz->has_kSZ_kSZ_gal_3h_fft
+ ptsz->has_kSZ_kSZ_gal_2h
+ ptsz->has_kSZ_kSZ_gal_3h
== _FALSE_
)
  return 0;


 // array of multipoles:
 double ln_ell_min = log(1e-3);
 double ln_ell_max = log(5e4);
 int n_ell = ptsz->n_ell_density_profile;
 int n_m = ptsz->n_m_density_profile;
 int n_z = ptsz->n_z_density_profile;

 class_alloc(ptsz->array_profile_ln_l,sizeof(double *)*n_ell,ptsz->error_message);

 // array of masses:
 double ln_m_min = log(1e6);
 double ln_m_max = log(1e18);


 class_alloc(ptsz->array_profile_ln_m,sizeof(double *)*n_m,ptsz->error_message);


 // array of redshifts:
 double ln_1pz_min = log(1.+ptsz->z1SZ);
 double ln_1pz_max = log(1.+ptsz->z2SZ);


 class_alloc(ptsz->array_profile_ln_1pz,sizeof(double *)*n_z,ptsz->error_message);
int index_m_z;

int index_l;
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
  ptsz->array_profile_ln_l[index_l] = ln_ell_min
                                      +index_l*(ln_ell_max-ln_ell_min)
                                      /(n_ell-1.);
}

int index_m;
for (index_m=0;
     index_m<n_m;
     index_m++)
{
  ptsz->array_profile_ln_m[index_m] = ln_m_min
                                      +index_m*(ln_m_max-ln_m_min)
                                      /(n_m-1.);
}

int index_z;
for (index_z=0;
     index_z<n_z;
     index_z++)
{
  ptsz->array_profile_ln_1pz[index_z] = ln_1pz_min
                                   +index_z*(ln_1pz_max-ln_1pz_min)
                                   /(n_z-1.);
}


class_alloc(ptsz->array_profile_ln_rho_at_lnl_lnM_z,n_ell*sizeof(double *),ptsz->error_message);
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
class_alloc(ptsz->array_profile_ln_rho_at_lnl_lnM_z[index_l],n_m*n_z*sizeof(double *),ptsz->error_message);
index_m_z = 0;
for (index_m=0;
     index_m<n_m;
     index_m++){

for (index_z=0;
     index_z<n_z;
     index_z++)
{
  ptsz->array_profile_ln_rho_at_lnl_lnM_z[index_l][index_m_z] = 1e-100; // initialize with super small number
  index_m_z += 1;
}

     }
}

int has_ksz_bkp = ptsz->has_kSZ_kSZ_gal_1h;
ptsz->has_kSZ_kSZ_gal_1h = _TRUE_; //pretend we need the tau_profile

//Parallelization of profile computation
/* initialize error management flag */


int abort;
double tstart, tstop;
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,\
ptsz,pba)\
private(tstart, tstop,index_l,index_z,index_m,index_m_z) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
#pragma omp flush(abort)
double * pvectsz;
double * pvecback;

class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);
class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
int index_pvectsz;
for (index_pvectsz=0;
     index_pvectsz<ptsz->tsz_size;
     index_pvectsz++){
       pvectsz[index_pvectsz] = 0.; // set everything to 0.
     }
index_m_z = 0;
for (index_z=0;
     index_z<n_z;
     index_z++){
for (index_m=0;
     index_m<n_m;
     index_m++){



  double z = exp(ptsz->array_profile_ln_1pz[index_z])-1.;
  double lnM = ptsz->array_profile_ln_m[index_m];
  double ell = exp(ptsz->array_profile_ln_l[index_l]);


  double tau;
  int first_index_back = 0;


  class_call_parallel(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             pba->error_message);

  class_call_parallel(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &first_index_back,
                               pvecback),
             pba->error_message,
             pba->error_message);


  // fill relevant entries
  pvectsz[ptsz->index_z] = z;

  pvectsz[ptsz->index_multipole_for_nfw_profile] = ell;
  pvectsz[ptsz->index_md] = 0; // avoid the if condition in rho_nfw for the pk mode computation

  pvectsz[ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);
  pvectsz[ptsz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
  double omega = pvecback[pba->index_bg_Omega_m];
  pvectsz[ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);


  double result;


  // pvectsz[ptsz->index_has_electron_density] = 1;
  // do_mass_conversions(lnM,z,pvecback,pvectsz,pba,ptsz);

 // only  do the integration of Battaglia profile
 // nfw has an analytical formula
 if (ptsz->tau_profile == 1){
  pvectsz[ptsz->index_m200c] = exp(lnM);
  class_call_parallel(mDEL_to_mVIR(pvectsz[ptsz->index_m200c],
                                   200.*(pvectsz[ptsz->index_Rho_crit]),
                                   pvectsz[ptsz->index_Delta_c],
                                   pvectsz[ptsz->index_Rho_crit],
                                   z,
                                   &pvectsz[ptsz->index_mVIR],
                                   ptsz,
                                   pba),
                  ptsz->error_message,
                  ptsz->error_message);
 //
 //  // rvir needed to cut off the integral --> e.g., xout = 50.*rvir/r200c
  pvectsz[ptsz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[ptsz->index_mVIR],pvectsz[ptsz->index_Delta_c],pvectsz[ptsz->index_Rho_crit],ptsz);

  pvectsz[ptsz->index_r200c] = pow(3.*pvectsz[ptsz->index_m200c]/(4.*_PI_*200.*pvectsz[ptsz->index_Rho_crit]),1./3.);
  pvectsz[ptsz->index_l200c] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_r200c];
  pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile] = pvectsz[ptsz->index_l200c];
  pvectsz[ptsz->index_rs] = pvectsz[ptsz->index_r200c];

 double result_int;
 class_call_parallel(two_dim_ft_nfw_profile(ptsz,pba,pvectsz,&result_int),
                     ptsz->error_message,
                     ptsz->error_message);


 result = result_int;
 double tau_normalisation = 1.;
 tau_normalisation = 4.*_PI_*pow(pvectsz[ptsz->index_r200c],3);

 result *= tau_normalisation;

}
else if (ptsz->tau_profile == 0){ // truncated nfw profile

   pvectsz[ptsz->index_m200c] = exp(lnM);
   pvectsz[ptsz->index_r200c] = pow(3.*pvectsz[ptsz->index_m200c]/(4.*_PI_*200.*pvectsz[ptsz->index_Rho_crit]),1./3.);
   pvectsz[ptsz->index_c200c] = get_c200c_at_m_and_z(pvectsz[ptsz->index_m200c],z,pba,ptsz);
   // set 1 for matter_type = tau
   double xout = 1.;//ptsz->x_out_truncated_nfw_profile;
   double r_delta = pvectsz[ptsz->index_r200c];
   double c_delta = pvectsz[ptsz->index_c200c];
   double m_delta = pvectsz[ptsz->index_m200c];

   // double r_delta = pvectsz[ptsz->index_radius_for_electron_density];
   // double c_delta = pvectsz[ptsz->index_concentration_for_electron_density];
   // double m_delta = pvectsz[ptsz->index_mass_for_electron_density];
   pvectsz[ptsz->index_rs] = r_delta/c_delta;

  pvectsz[ptsz->index_multipole_for_truncated_nfw_profile] = pvectsz[ptsz->index_multipole_for_nfw_profile];
  double l = pvectsz[ptsz->index_multipole_for_truncated_nfw_profile];
  double chi = sqrt(pvectsz[ptsz->index_chi2]);
  double k = (l+0.5)/chi;
   result =  evaluate_truncated_nfw_profile(k,r_delta,c_delta,xout,pvectsz,pba,ptsz);
   result *= m_delta/(4.*_PI_*pow(pvectsz[ptsz->index_rs],3));
   double f_b = pba->Omega0_b/ptsz->Omega_m_0;
   result *= f_b;

   double tau_normalisation = 4.*_PI_*pow(pvectsz[ptsz->index_rs],3);
   result *= tau_normalisation;
 }

  ptsz->array_profile_ln_rho_at_lnl_lnM_z[index_l][index_m_z] = result;

  index_m_z += 1;
     }


     }

     free(pvectsz);
     free(pvecback);
}

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in tab profile parallel region (loop over ell's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

// restore initial state:
ptsz->has_kSZ_kSZ_gal_1h = has_ksz_bkp;


                                      }





double get_pressure_profile_at_l_M_z(double l_asked, double m_asked, double z_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
  double l = log(l_asked);


  // if (ptsz->tau_profile == 1){
  // find the closest l's in the grid:
  int id_l_low;
  int id_l_up;
  int n_ell = ptsz->n_ell_pressure_profile;
  int n_m = ptsz->n_m_pressure_profile;
  int n_z = ptsz->n_z_pressure_profile;
  r8vec_bracket(n_ell,ptsz->array_pressure_profile_ln_l,l,&id_l_low,&id_l_up);

  // interpolate 2d at l_low:

 double ln_rho_low = pwl_interp_2d(n_m,
                                n_z,
                                ptsz->array_pressure_profile_ln_m,
                                ptsz->array_pressure_profile_ln_1pz,
                                ptsz->array_pressure_profile_ln_p_at_lnl_lnM_z[id_l_low-1],
                                1,
                                &m,
                                &z);

 double ln_rho_up = pwl_interp_2d(n_m,
                                n_z,
                                ptsz->array_pressure_profile_ln_m,
                                ptsz->array_pressure_profile_ln_1pz,
                                ptsz->array_pressure_profile_ln_p_at_lnl_lnM_z[id_l_up-1],
                                1,
                                &m,
                                &z);
 double ln_l_low = ptsz->array_pressure_profile_ln_l[id_l_low-1];
 double ln_l_up = ptsz->array_pressure_profile_ln_l[id_l_up-1];

 return ln_rho_low + ((l - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low);



}


// Tabulate 2D Fourier transform of pressure profile on a [z - ln_M - ln_ell] grid
int tabulate_pressure_profile_B12(struct background * pba,
                                  struct tszspectrum * ptsz){

 // array of multipoles:
 double ln_ell_min = log(0.5);
 double ln_ell_max = log(5e4);
 int n_ell = ptsz->n_ell_pressure_profile;
 int n_m = ptsz->n_m_pressure_profile;
 int n_z = ptsz->n_z_pressure_profile;

 class_alloc(ptsz->array_pressure_profile_ln_l,sizeof(double *)*n_ell,ptsz->error_message);

 // array of masses:
 double ln_m_min = log(1e7);
 double ln_m_max = log(1e17);


 class_alloc(ptsz->array_pressure_profile_ln_m,sizeof(double *)*n_m,ptsz->error_message);


 // array of redshifts:
 double ln_1pz_min = log(1.+ptsz->z1SZ);
 double ln_1pz_max = log(1.+ptsz->z2SZ);


 class_alloc(ptsz->array_pressure_profile_ln_1pz,sizeof(double *)*n_z,ptsz->error_message);
int index_m_z;

int index_l;
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
  ptsz->array_pressure_profile_ln_l[index_l] = ln_ell_min
                                      +index_l*(ln_ell_max-ln_ell_min)
                                      /(n_ell-1.);
}

int index_m;
for (index_m=0;
     index_m<n_m;
     index_m++)
{
  ptsz->array_pressure_profile_ln_m[index_m] = ln_m_min
                                      +index_m*(ln_m_max-ln_m_min)
                                      /(n_m-1.);
}

int index_z;
for (index_z=0;
     index_z<n_z;
     index_z++)
{
  ptsz->array_pressure_profile_ln_1pz[index_z] = ln_1pz_min
                                   +index_z*(ln_1pz_max-ln_1pz_min)
                                   /(n_z-1.);
}


class_alloc(ptsz->array_pressure_profile_ln_p_at_lnl_lnM_z,n_ell*sizeof(double *),ptsz->error_message);
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
class_alloc(ptsz->array_pressure_profile_ln_p_at_lnl_lnM_z[index_l],n_m*n_z*sizeof(double *),ptsz->error_message);
index_m_z = 0;
for (index_m=0;
     index_m<n_m;
     index_m++){
//class_alloc(ptsz->array_profile_ln_rho_at_lnl_lnM_z[index_l][index_m_z],n_z*sizeof(double ),ptsz->error_message);

for (index_z=0;
     index_z<n_z;
     index_z++)
{
  // ptsz->array_profile_ln_rho_at_lnl_lnM_z[index_l][index_m_z] = -100.; // initialize with super small number
  ptsz->array_pressure_profile_ln_p_at_lnl_lnM_z[index_l][index_m_z] = 1e-100; // initialize with super small number
  index_m_z += 1;
}

     }
}


//Parallelization of profile computation
/* initialize error management flag */


int abort;
double tstart, tstop;
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,\
ptsz,pba)\
private(tstart, tstop,index_l,index_z,index_m,index_m_z) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
#pragma omp flush(abort)
double * pvectsz;
double * pvecback;
class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);
class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
int index_pvectsz;
for (index_pvectsz=0;
     index_pvectsz<ptsz->tsz_size;
     index_pvectsz++){
       pvectsz[index_pvectsz] = 0.; // set everything to 0.
     }
index_m_z = 0;
for (index_z=0;
     index_z<n_z;
     index_z++){
for (index_m=0;
     index_m<n_m;
     index_m++){



  double z = exp(ptsz->array_pressure_profile_ln_1pz[index_z])-1.;
  double lnM = ptsz->array_pressure_profile_ln_m[index_m];
  double ell = exp(ptsz->array_pressure_profile_ln_l[index_l]);


  double tau;
  int first_index_back = 0;


  class_call_parallel(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             pba->error_message);

  class_call_parallel(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &first_index_back,
                               pvecback),
             pba->error_message,
             pba->error_message);


  // fill relevant entries
  pvectsz[ptsz->index_z] = z;

  pvectsz[ptsz->index_multipole_for_pressure_profile] = ell;
  pvectsz[ptsz->index_md] = 0; // avoid the if condition in rho_nfw for the pk mode computation

  pvectsz[ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);
  pvectsz[ptsz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
  double omega = pvecback[pba->index_bg_Omega_m];
  pvectsz[ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);


  double result;
  pvectsz[ptsz->index_m200c] = exp(lnM);
  class_call_parallel(mDEL_to_mVIR(pvectsz[ptsz->index_m200c],
                                   200.*(pvectsz[ptsz->index_Rho_crit]),
                                   pvectsz[ptsz->index_Delta_c],
                                   pvectsz[ptsz->index_Rho_crit],
                                   z,
                                   &pvectsz[ptsz->index_mVIR],
                                   ptsz,
                                   pba),
                  ptsz->error_message,
                  ptsz->error_message);
 //
 //  // rvir needed to cut off the integral --> e.g., xout = 50.*rvir/r200c
  pvectsz[ptsz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[ptsz->index_mVIR],pvectsz[ptsz->index_Delta_c],pvectsz[ptsz->index_Rho_crit],ptsz);
  pvectsz[ptsz->index_r200c] = pow(3.*pvectsz[ptsz->index_m200c]/(4.*_PI_*200.*pvectsz[ptsz->index_Rho_crit]),1./3.);
  pvectsz[ptsz->index_l200c] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_r200c];
  // pvectsz[ptsz->index_characteristic_multipole_for_profile] = pvectsz[ptsz->index_l200c];
  pvectsz[ptsz->index_rs] = pvectsz[ptsz->index_r200c];///pvectsz[ptsz->index_c200c];
 double result_int;
 class_call_parallel(two_dim_ft_pressure_profile(ptsz,pba,pvectsz,&result_int),
                                    ptsz->error_message,
                                    ptsz->error_message);
 result = result_int;



 // double tau_normalisation = 1.;//pvectsz[ptsz->index_m200m];///(4.*_PI_*pow(pvectsz[ptsz->index_rs],3.));
 // tau_normalisation = 4.*_PI_*pow(pvectsz[ptsz->index_r200c],3);//*pvectsz[ptsz->index_Rho_crit];
 // result *= tau_normalisation;


  // ptsz->array_profile_ln_rho_at_lnl_lnM_z[index_l][index_m_z] = log(result);
  ptsz->array_pressure_profile_ln_p_at_lnl_lnM_z[index_l][index_m_z] = result;
  // printf("ell = %.3e z = %.3e m = %.3e lnrho = %.3e\n",ell,z,exp(lnM),log(result));
  index_m_z += 1;
     }


     }

free(pvectsz);
free(pvecback);

}

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in tab profile parallel region (loop over ell's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;

                                      }



// Tabulate 2D Fourier transform of density profile on a [ln_ell_over_ell_char] grid
int tabulate_pressure_profile_gNFW(struct background * pba,
                                   struct tszspectrum * ptsz){

// if (ptsz->has_kSZ_kSZ_gal_1h + ptsz->has_kSZ_kSZ_gal_2h + ptsz->has_kSZ_kSZ_gal_3h == _FALSE_)
//   return 0;


 // array of multipoles:
 // this is (l+0.5)/ls
 double ln_ell_min = log(1e-2);
 double ln_ell_max = log(50.);
 int n_ell = ptsz->array_profile_ln_PgNFW_at_lnl_over_ls_size;


 class_alloc(ptsz->array_profile_ln_l_over_ls,sizeof(double *)*n_ell,ptsz->error_message);

 // // array of masses:
 // double ln_m_min = log(1e8);
 // double ln_m_max = log(1e18);

 //
 // class_alloc(ptsz->array_profile_ln_m,sizeof(double *)*n_m,ptsz->error_message);
 //

 // // array of redshifts:
 // double ln_1pz_min = log(1.+ptsz->z1SZ);
 // double ln_1pz_max = log(1.+ptsz->z2SZ);


//  class_alloc(ptsz->array_profile_ln_1pz,sizeof(double *)*n_z,ptsz->error_message);
// int index_m_z;

int index_l;
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
  ptsz->array_profile_ln_l_over_ls[index_l] = ln_ell_min
                                              +index_l*(ln_ell_max-ln_ell_min)
                                              /(n_ell-1.);
}


class_alloc(ptsz->array_profile_ln_PgNFW_at_lnl_over_ls,n_ell*sizeof(double *),ptsz->error_message);
for (index_l=0;
     index_l<n_ell;
     index_l++)
{

  ptsz->array_profile_ln_PgNFW_at_lnl_over_ls[index_l] = 1e-100; // initialize with super small number

}

// int has_ksz_bkp = ptsz->has_kSZ_kSZ_gal_1h;
// ptsz->has_kSZ_kSZ_gal_1h = _TRUE_; //pretend we need the tau_profile

//Parallelization of profile computation
/* initialize error management flag */


int abort;
double tstart, tstop;
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,\
ptsz,pba)\
private(tstart, tstop,index_l) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
#pragma omp flush(abort)
double * pvectsz;
double * pvecback;
class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);
class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
int index_pvectsz;
for (index_pvectsz=0;
     index_pvectsz<ptsz->tsz_size;
     index_pvectsz++){
       pvectsz[index_pvectsz] = 0.; // set everything to 0.
     }





  pvectsz[ptsz->index_l500] = 1.;
  pvectsz[ptsz->index_multipole_for_pressure_profile] = exp(ptsz->array_profile_ln_l_over_ls[index_l])-0.5;
  double result;


  class_call_parallel(two_dim_ft_pressure_profile(ptsz,pba,pvectsz,&result),
                                     ptsz->error_message,
                                     ptsz->error_message);



  // ptsz->array_profile_ln_rho_at_lnl_lnM_z[index_l][index_m_z] = log(result);
  ptsz->array_profile_ln_PgNFW_at_lnl_over_ls[index_l] = log(result);
  // printf("ell/ells = %.3e ln_pgnfw = %.3e\n",exp(ptsz->array_profile_ln_l_over_ls[index_l]),log(result));
  // printf("ell/ells = %.3e ln_pgnfw = %.3e\n",exp(ptsz->array_profile_ln_l_over_ls[index_l]),log(result));



  // printf("freeing pp pvectsz, pvecback\n");
  // free(pvectsz);
  // free(pvecback);
  // printf("freed\n");


}

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in tab profile parallel region (loop over ell's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

return _SUCCESS_;

                                      }




/**
 * This routine computes 2d ft of pressure profile at ell/ell_characteristic
 *
 * @param pba   Input: pointer to background structure
 * @param ppm   Input: pointer to primordial structure
 * @param psp   Input: pointer to spectra structure
 * @param z     Input: redshift
 * @param R     Input: radius in Mpc
 * @param sigma Output: variance in a sphere of radius R (dimensionless)
 */

int two_dim_ft_pressure_profile(struct tszspectrum * ptsz,
                                struct background * pba,
                                double * pvectsz,
                                double * result
                          ) {


  //GSL
  // QAWO

  double xin = ptsz->x_inSZ;
  double xout = 0.;
if (ptsz->pressure_profile == 4) { //for Battaglia et al 2012 pressure profile
    double rvir = pvectsz[ptsz->index_rVIR]; //in Mpc/h
    double r200c = pvectsz[ptsz->index_r200c]; //in Mpc/h
    xout = 1.5*rvir/r200c; // the truncation radius is in multiples of rvir
    }
else{
    xout = ptsz->x_outSZ; // in all other cases the truncation radius is in multiples of rs=r_delta/c_delta
    }

  double delta_l = xout - xin;

  gsl_integration_workspace * w;
  gsl_integration_qawo_table * wf;

  int size_w = 3000;
  w = gsl_integration_workspace_alloc(size_w);


  double w0;
  if (ptsz->pressure_profile == 4) //for Battaglia et al 2012 pressure profile
  w0 = (pvectsz[ptsz->index_multipole_for_pressure_profile]+0.5)/pvectsz[ptsz->index_l200c];
  else
  w0 =  (pvectsz[ptsz->index_multipole_for_pressure_profile]+0.5)/pvectsz[ptsz->index_l500];

  wf = gsl_integration_qawo_table_alloc(w0, delta_l,GSL_INTEG_SINE,50);

  struct Parameters_for_integrand_patterson_pp V;
  V.ptsz = ptsz;
  V.pba = pba;
  V.pvectsz = pvectsz;

  void * params = &V;

  gsl_function F;
  F.function = &integrand_patterson_test_pp;
  F.params = params;

  double eps_abs = ptsz->pressure_profile_epsabs;
  double eps_rel = ptsz->pressure_profile_epsrel;

  double result_gsl, error;
  int limit = size_w; //number of sub interval
  gsl_integration_qawo(&F,ptsz->x_for_pp[0],eps_abs,eps_rel,limit,w,wf,&result_gsl,&error);

  *result = result_gsl;

  gsl_integration_qawo_table_free(wf);
  gsl_integration_workspace_free(w);

  return _SUCCESS_;

}





/**
 * This routine computes sigma(R) given P(k) (does not check that k_max is large
 * enough)
 *
 * @param pba   Input: pointer to background structure
 * @param ppm   Input: pointer to primordial structure
 * @param psp   Input: pointer to spectra structure
 * @param z     Input: redshift
 * @param R     Input: radius in Mpc
 * @param sigma Output: variance in a sphere of radius R (dimensionless)
 */

int spectra_sigma_for_tSZ(
                          struct background * pba,
                          struct primordial * ppm,
                          struct nonlinear *pnl,
                          struct tszspectrum * ptsz,
                          double R,
                          double z,
                          double * sigma
                          ) {

  double pk;
  double * pk_ic = NULL;

  double * array_for_sigma;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double k,W,x;


  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              ptsz->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);

  for (i=0;i<ptsz->ln_k_size_for_tSZ;i++) {
    k=exp(ptsz->ln_k_for_tSZ[i]);
    if (i == (ptsz->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
    W=3./x/x/x*(sin(x)-x*cos(x));

    /*
    class_call(spectra_pk_at_k_and_z(pba,ppm,psp,k,z,&pk,pk_ic),
               psp->error_message,
               psp->error_message);*/

     class_call(nonlinear_pk_at_k_and_z(
                                       pba,
                                       ppm,
                                       pnl,
                                       pk_linear,
                                       k,
                                       z,
                                       pnl->index_pk_m,
                                       &pk, // number *out_pk_l
                                       pk_ic // array out_pk_ic_l[index_ic_ic]
                                     ),
                                     pnl->error_message,
                                     pnl->error_message);


    array_for_sigma[i*index_num+index_k]=k;
    array_for_sigma[i*index_num+index_y]=k*k*pk*W*W;
  }




  class_call(array_spline(array_for_sigma,
                          index_num,
                          ptsz->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        ptsz->ln_k_size_for_tSZ,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  free(array_for_sigma);


  *sigma = sqrt(*sigma/(2.*_PI_*_PI_));

  return _SUCCESS_;

}


//This routine computes dSigma2/dR
//at R and z

int spectra_sigma_prime(
                        struct background * pba,
                        struct primordial * ppm,
                        struct nonlinear *pnl,
                        struct tszspectrum * ptsz,
                        double R,
                        double z,
                        double * sigma_prime
                        ) {

  double pk;
  double * pk_ic = NULL;

  double * array_for_sigma;
  int index_num;
  int index_k;
  int index_y;
  int index_ddy;
  int i;

  double k,W,x,W_prime;



  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              ptsz->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);

  for (i=0;i<ptsz->ln_k_size_for_tSZ;i++) {
    k=exp(ptsz->ln_k_for_tSZ[i]);
    if (i == (ptsz->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
    W=3./x/x/x*(sin(x)-x*cos(x));
    W_prime=3./x/x*sin(x)-9./x/x/x/x*(sin(x)-x*cos(x));

    //class_call(spectra_pk_at_k_and_z(pba,ppm,psp,k,z,&pk,pk_ic),
    //           psp->error_message,
    //           psp->error_message);

    class_call(nonlinear_pk_at_k_and_z(
                                      pba,
                                      ppm,
                                      pnl,
                                      pk_linear,
                                      k,
                                      z,
                                      pnl->index_pk_m,
                                      &pk, // number *out_pk_l
                                      pk_ic // array out_pk_ic_l[index_ic_ic]
                                    ),
                                    pnl->error_message,
                                    pnl->error_message);


    array_for_sigma[i*index_num+index_k]=k;
    array_for_sigma[i*index_num+index_y]=k*k*pk*k*2.*W*W_prime;
  }

  class_call(array_spline(array_for_sigma,
                          index_num,
                          ptsz->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        ptsz->ln_k_size_for_tSZ,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma_prime,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  free(array_for_sigma);



  *sigma_prime = *sigma_prime/(2.*_PI_*_PI_);

  return _SUCCESS_;

}



//This routine reads the tabulated
//pressure profiles,
//and stores the tabulated values.

int external_pressure_profile_init(struct precision * ppr, struct tszspectrum * ptsz)
{

if (ptsz->pressure_profile != 0 && ptsz->pressure_profile != 2 )
  return 0;
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
     +ptsz->has_sz_2halo
     +ptsz->has_sz_ps
     +ptsz->has_mean_y
     +ptsz->has_dydz
     == 0)
     return 0;

  if (ptsz->sz_verbose > 0)
    printf("-> Using tabulated pressure profile transform\n");

  class_alloc(ptsz->PP_lnx,sizeof(double *)*100,ptsz->error_message);
  class_alloc(ptsz->PP_lnI,sizeof(double *)*100,ptsz->error_message);
  class_alloc(ptsz->PP_d2lnI,sizeof(double *)*100,ptsz->error_message);

  //char arguments[_ARGUMENT_LENGTH_MAX_];
  char line[_LINE_LENGTH_MAX_];
  //char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  FILE *process;
  int n_data_guess, n_data = 0;
  double *lnx = NULL, *lnI = NULL, *d2lnI = NULL, *tmp = NULL;
  double this_lnx, this_lnI, this_d2lnI;
  int status;
  int index_x;


  /** 1. Initialization */
  /* Prepare the data (with some initial size) */
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));
  lnI = (double *)malloc(n_data_guess*sizeof(double));
  d2lnI = (double *)malloc(n_data_guess*sizeof(double));


  /* Prepare the command */
  /* If the command is just a "cat", no arguments need to be passed */
  // if(strncmp("cat ", ptsz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];
  if (ptsz->pressure_profile==0){
    class_open(process,ptsz->P13_file, "r",ptsz->error_message);
  }
  else if (ptsz->pressure_profile==2){
    if (ptsz->sz_verbose > 0)
      printf("-> Openning the pressure profile file for A10\n");
    //class_open(process,"sz_auxiliary_files/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_A10.txt", "r",ptsz->error_message);
    class_open(process,ptsz->A10_file, "r",ptsz->error_message);
    printf("-> File Name: %s\n",ptsz->A10_file);
    // printf("-> File Name: %s\n",ppr->sBBN_file);

  }

    // sprintf(Filepath,
    //         "%s%s",
    //         // "%s%s%s",
    //         "cat ",
    //         // ptsz->path_to_class,
    //         "/sz_auxiliary_files/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_P13.txt");
  //   sprintf(Filepath,
  //           "%s%s",
  //           // "%s%s%s",
  //           "cat ",
  //           ptsz->path_to_class,
  //           "/sz_auxiliary_files/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_A10.txt");
  // process = popen(Filepath, "r");
  if (ptsz->sz_verbose > 0)
    printf("-> Scanning the pressure profile file\n");
  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf %lf %lf", &this_lnx, &this_lnI, &this_d2lnI);
    // printf("lnx = %e\n",this_lnx);




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
      tmp = (double *)realloc(d2lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 ptsz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      d2lnI = tmp;
    };
    /* Store */
    lnx[n_data]   = this_lnx;
    lnI[n_data]   = this_lnI;
    d2lnI[n_data] = this_d2lnI;
    n_data++;
    /* Check ascending order of the k's */
    if(n_data>1) {
      class_test(lnx[n_data-1] <= lnx[n_data-2],
                 ptsz->error_message,
                 "The ell/ell500's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  // printf("ptsz->PP_lnI[index_x] = %e\n",lnI[0]);

  status = fclose(process);
  // printf("ptsz->PP_lnI[index_x] = %e\n",lnI[0]);

  // fclose(process);
  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  ptsz->PP_lnx_size = n_data;
  /** Make room */
  // printf("ptsz->PP_lnI[index_x] = %d\n",n_data);

  class_realloc(ptsz->PP_lnx,
                ptsz->PP_lnx,
                ptsz->PP_lnx_size*sizeof(double),
                ptsz->error_message);
  class_realloc(ptsz->PP_lnI,
                ptsz->PP_lnI,
                ptsz->PP_lnx_size*sizeof(double),
                ptsz->error_message);
  class_realloc(ptsz->PP_d2lnI,
                ptsz->PP_d2lnI,
                ptsz->PP_lnx_size*sizeof(double),
                ptsz->error_message);


  /** Store them */
  for (index_x=0; index_x<ptsz->PP_lnx_size; index_x++) {
    ptsz->PP_lnx[index_x] = lnx[index_x];
    ptsz->PP_lnI[index_x] = lnI[index_x];
    ptsz->PP_d2lnI[index_x] = d2lnI[index_x];
    // printf("ptsz->PP_lnI[index_x] = %e\n",ptsz->PP_lnI[index_x]);

  };

  /** Release the memory used locally */
  free(lnx);
  free(lnI);
  free(d2lnI);

  return _SUCCESS_;
}



//This routine reads the tabulated
//noise curve for yxy covariance,
//and stores the tabulated values.

int load_unbinned_nl_yy(struct tszspectrum * ptsz)
{


// don't load  if none of the following are required:
if ( (ptsz->include_noise_cov_y_y != 1 )){
  // if (ptsz->sz_verbose>=1)
  //   printf("-> noise curve for yxy covariance not requested\n");
  return 0;
}

if (ptsz->sz_verbose >= 1)
  printf("-> loading the noise curve for yxy covariance\n");


  class_alloc(ptsz->unbinned_nl_yy_ell,sizeof(double *)*100,ptsz->error_message);
  class_alloc(ptsz->unbinned_nl_yy_n_ell,sizeof(double *)*100,ptsz->error_message);
  //class_alloc(ptsz->PP_d2lnI,sizeof(double *)*100,ptsz->error_message);

  //char arguments[_ARGUMENT_LENGTH_MAX_];
  char line[_LINE_LENGTH_MAX_];
  //char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  FILE *process;
  int n_data_guess, n_data = 0;
  double *lnx = NULL, *lnI = NULL,  *tmp = NULL;
  double this_lnx, this_lnI;
  int status;
  int index_x;


  /** 1. Initialization */
  /* Prepare the data (with some initial size) */
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));
  lnI = (double *)malloc(n_data_guess*sizeof(double));



  /* Prepare the command */
  /* If the command is just a "cat", no arguments need to be passed */
  // if(strncmp("cat ", ptsz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];

    sprintf(Filepath,
            "%s%s",
            "cat ",
            ptsz->full_path_to_noise_curve_for_y_y);
  //printf("-> HI2 loading the noise curve for yxy covariance\n");
  process = popen(Filepath, "r");
  printf("-> %s\n",Filepath);

  //int il = 0;
  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {
    //printf("%d\n",il);
    //il++;
    sscanf(line, "%lf %lf ", &this_lnx, &this_lnI);



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
  status = pclose(process);
  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  ptsz->unbinned_nl_yy_size = n_data;
  /** Make room */

  class_realloc(ptsz->unbinned_nl_yy_ell,
                ptsz->unbinned_nl_yy_ell,
                ptsz->unbinned_nl_yy_size*sizeof(double),
                ptsz->error_message);
  class_realloc(ptsz->unbinned_nl_yy_n_ell,
                ptsz->unbinned_nl_yy_n_ell,
                ptsz->unbinned_nl_yy_size*sizeof(double),
                ptsz->error_message);



  /** Store them */
  for (index_x=0; index_x<ptsz->unbinned_nl_yy_size; index_x++) {
    ptsz->unbinned_nl_yy_ell[index_x] = lnx[index_x];
    ptsz->unbinned_nl_yy_n_ell[index_x] = lnI[index_x];
    //printf("z=%.3e phig=%.3e\n",ptsz->unbinned_nl_yy_ell[index_x],ptsz->unbinned_nl_yy_n_ell[index_x]);
  };

  //exit(0);
  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}




//This routine reads the tabulated
//dndz galaxy counts,
//and stores the tabulated values.

int load_normalized_dndz(struct tszspectrum * ptsz)
{

// don't load the unwise  dndz  if none of the following are required:
// all quantities requiring galaxy or lensmag need that:
if (   (ptsz->has_tSZ_gal_1h != _TRUE_ )
    && (ptsz->has_tSZ_gal_2h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_1h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_1h_fft != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_2h_fft != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_3h_fft != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_2h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_3h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_hf != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_lensmag_1halo != _TRUE_ )
    && (ptsz->has_gal_gal_1h != _TRUE_ )
    && (ptsz->has_gal_lens_1h != _TRUE_ )
    && (ptsz->has_gal_lens_2h != _TRUE_ )
    && (ptsz->has_gal_cib_1h != _TRUE_ )
    && (ptsz->has_gal_cib_2h != _TRUE_ )
    && (ptsz->has_gal_lens_hf != _TRUE_ )
    && (ptsz->has_gal_lensmag_1h != _TRUE_ )
    && (ptsz->has_gal_lensmag_2h != _TRUE_ )
    && (ptsz->has_gal_lensmag_hf != _TRUE_ )
    && (ptsz->has_tSZ_lensmag_1h != _TRUE_ )
    && (ptsz->has_tSZ_lensmag_2h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_1h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_2h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_hf != _TRUE_ )
    && (ptsz->has_lens_lensmag_1h != _TRUE_ )
    && (ptsz->has_lens_lensmag_2h != _TRUE_ )
    && (ptsz->has_lens_lensmag_hf != _TRUE_ )
    && (ptsz->has_gal_gal_2h != _TRUE_ )
    && (ptsz->has_gal_gal_hf != _TRUE_ ))
  return 0;

if (ptsz->sz_verbose>=1){
    if (ptsz->galaxy_sample == 0)
    printf("-> Loading dndz WIxSC\n");
    if (ptsz->galaxy_sample == 1)
    printf("-> Loading dndz unwise\n");
    if (ptsz->galaxy_sample == 2)
    printf("-> Loading dndz file\n");
    }

  class_alloc(ptsz->normalized_dndz_z,sizeof(double *)*100,ptsz->error_message);
  class_alloc(ptsz->normalized_dndz_phig,sizeof(double *)*100,ptsz->error_message);

  //class_alloc(ptsz->PP_d2lnI,sizeof(double *)*100,ptsz->error_message);

  //char arguments[_ARGUMENT_LENGTH_MAX_];
  char line[_LINE_LENGTH_MAX_];
  //char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  FILE *process;
  int n_data_guess, n_data = 0;
  double *lnx = NULL, *lnI = NULL,  *tmp = NULL;
  double this_lnx, this_lnI, this_lnJ, this_lnK;
  int status;
  int index_x;


  /** 1. Initialization */
  /* Prepare the data (with some initial size) */
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));
  lnI = (double *)malloc(n_data_guess*sizeof(double));



  /* Prepare the command */
  /* If the command is just a "cat", no arguments need to be passed */
  // if(strncmp("cat ", ptsz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];

  //unwise
  if (ptsz->galaxy_sample == 1){
    if (ptsz->sz_verbose > 0){
      printf("-> Openning the dndz file for unWISE galaxies\n");
      printf("-> File Name: %s\n",ptsz->UNWISE_dndz_file);
      // printf("-> File Name: %s\n",ptsz->UNWISE_fdndz_file);
      // printf("-> File Name: %s\n",ptsz->A10_file);
    }
  class_open(process,ptsz->UNWISE_dndz_file, "r",ptsz->error_message);
  if (ptsz->sz_verbose > 0)
    printf("-> File opened successfully\n");
    // sprintf(Filepath,
    //         "%s%s",
    //         "cat ",
    //         //ptsz->path_to_class,
    //         "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/UNWISE_galaxy_ditributions/normalised_dndz.txt");
            }

  else if (ptsz->galaxy_sample == 0){
    if (ptsz->sz_verbose > 0){
      printf("-> Openning the dndz file for WISC3 galaxies\n");
      printf("-> File Name: %s\n",ptsz->WISC3_dndz_file);
      // printf("-> File Name: %s\n",ptsz->UNWISE_fdndz_file);
      // printf("-> File Name: %s\n",ptsz->A10_file);
    }
  class_open(process,ptsz->WISC3_dndz_file, "r",ptsz->error_message);
    if (ptsz->sz_verbose > 0)
      printf("-> File opened successfully\n");
 //  sprintf(Filepath,
 //          "%s%s",
 //          "cat ",
 //          //ptsz->path_to_class,
 //          "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/run_scripts/yxg/data/dndz/WISC_bin3.txt");
 // process = popen(Filepath, "r");
        }

  else if (ptsz->galaxy_sample == 2){
    if (ptsz->sz_verbose > 0){
      printf("-> Openning the dndz file for unWISE galaxies\n");
      printf("-> File Name: %s\n",ptsz->full_path_to_dndz_gal);
      // printf("-> File Name: %s\n",ptsz->UNWISE_fdndz_file);
      // printf("-> File Name: %s\n",ptsz->A10_file);
    }
  class_open(process,ptsz->full_path_to_dndz_gal, "r",ptsz->error_message);
  if (ptsz->sz_verbose > 0)
    printf("-> File opened successfully\n");
  // sprintf(Filepath,
  //         "%s%s",
  //         "cat ",
  //         ptsz->full_path_to_dndz_gal);
  //         //"/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/run_scripts/yxg/data/dndz/unwise_red.txt");
  // process = popen(Filepath, "r");
        }


  // process = popen(Filepath, "r");

  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {

    // unWISE load and read column depending on the requested color
      if (ptsz->galaxy_sample == 1){
    sscanf(line, "%lf %lf %lf %lf", &this_lnx, &this_lnI, &this_lnJ, &this_lnK);
    //sscanf(line, "%lf %lf", &this_lnx, &this_lnI);

    // red
    if (ptsz->unwise_galaxy_sample_id == 0)
    this_lnI = this_lnK;

    // green
    if (ptsz->unwise_galaxy_sample_id == 1 || ptsz->unwise_galaxy_sample_id == 2)
    this_lnI = this_lnJ;

    // blue
    //if (ptsz->unwise_galaxy_sample_id == 3)
    //this_lnI = this_lnI;

    // printf("lnx = %e\n",this_lnx);
                                    }

  // WIxSC and "other": just two column files
  else if ((ptsz->galaxy_sample == 0) || (ptsz->galaxy_sample == 2)){

    sscanf(line, "%lf %lf ", &this_lnx, &this_lnI);
    // printf("lnx = %e\n",this_lnx);
  }


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
  // if (ptsz->galaxy_sample == 2){
  // // if (ptsz->galaxy_sample == 2 || ptsz->galaxy_sample == 0 ){
  //   status = pclose(process);
  // }
  // else{
    status = fclose(process);
  //}
  // status = pclose(process);

  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  ptsz->normalized_dndz_size = n_data;
  /** Make room */

  class_realloc(ptsz->normalized_dndz_z,
                ptsz->normalized_dndz_z,
                ptsz->normalized_dndz_size*sizeof(double),
                ptsz->error_message);
  class_realloc(ptsz->normalized_dndz_phig,
                ptsz->normalized_dndz_phig,
                ptsz->normalized_dndz_size*sizeof(double),
                ptsz->error_message);



  /** Store them */
  for (index_x=0; index_x<ptsz->normalized_dndz_size; index_x++) {
    ptsz->normalized_dndz_z[index_x] = lnx[index_x];
    ptsz->normalized_dndz_phig[index_x] = lnI[index_x];
    //print("z=%.3e phig=%.3e\n",ptsz->normalized_dndz_z[index_x])
  };

  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}


int load_normalized_fdndz(struct tszspectrum * ptsz)
{

// don't load the unwise  dndz  if none of the following are required:
if (   (ptsz->has_tSZ_gal_1h != _TRUE_ )
    && (ptsz->has_tSZ_gal_2h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_1h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_1h_fft != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_2h_fft != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_3h_fft != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_2h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_3h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_hf != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_lensmag_1halo != _TRUE_ )
    && (ptsz->has_gal_gal_1h != _TRUE_ )
    && (ptsz->has_gal_lens_1h != _TRUE_ )
    && (ptsz->has_gal_lens_2h != _TRUE_ )
    && (ptsz->has_gal_cib_1h != _TRUE_ )
    && (ptsz->has_gal_cib_2h != _TRUE_ )
    && (ptsz->has_gal_lens_hf != _TRUE_ )
    && (ptsz->has_gal_lensmag_1h != _TRUE_ )
    && (ptsz->has_gal_lensmag_2h != _TRUE_ )
    && (ptsz->has_gal_lensmag_hf != _TRUE_ )
    && (ptsz->has_tSZ_lensmag_1h != _TRUE_ )
    && (ptsz->has_tSZ_lensmag_2h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_1h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_2h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_hf != _TRUE_ )
    && (ptsz->has_lens_lensmag_1h != _TRUE_ )
    && (ptsz->has_lens_lensmag_2h != _TRUE_ )
    && (ptsz->has_lens_lensmag_hf != _TRUE_ )
    && (ptsz->has_gal_gal_2h != _TRUE_ )
    && (ptsz->has_gal_gal_hf != _TRUE_ ))
  return 0;

if ((ptsz->galaxy_sample == 0) || (ptsz->galaxy_sample == 2))
  return 0;

if (ptsz->sz_verbose >= 1)
printf("-> Loading fdndz unwise\n");

  class_alloc(ptsz->normalized_fdndz_z,sizeof(double *)*100,ptsz->error_message);
  class_alloc(ptsz->normalized_fdndz_phig,sizeof(double *)*100,ptsz->error_message);

  //class_alloc(ptsz->PP_d2lnI,sizeof(double *)*100,ptsz->error_message);

  //char arguments[_ARGUMENT_LENGTH_MAX_];
  char line[_LINE_LENGTH_MAX_];
  //char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  FILE *process;
  int n_data_guess, n_data = 0;
  double *lnx = NULL, *lnI = NULL,  *tmp = NULL;
  double this_lnx, this_lnI, this_lnJ, this_lnK;
  int status;
  int index_x;


  /** 1. Initialization */
  /* Prepare the data (with some initial size) */
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));
  lnI = (double *)malloc(n_data_guess*sizeof(double));



  /* Prepare the command */
  /* If the command is just a "cat", no arguments need to be passed */
  // if(strncmp("cat ", ptsz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  // char Filepath[_ARGUMENT_LENGTH_MAX_];
  //
  // //unwise
  //
  //   sprintf(Filepath,
  //           "%s%s",
  //           "cat ",
  //           //ptsz->path_to_class,
  //           "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/UNWISE_galaxy_ditributions/normalised_fdndz.txt");

  class_open(process,ptsz->UNWISE_fdndz_file, "r",ptsz->error_message);

  // process = popen(Filepath, "r");

  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {

    // unWISE load and read column depending on the requested color

    sscanf(line, "%lf %lf %lf %lf", &this_lnx, &this_lnI, &this_lnJ, &this_lnK);
    //sscanf(line, "%lf %lf", &this_lnx, &this_lnI);

    // red
    if (ptsz->unwise_galaxy_sample_id == 0)
    this_lnI = this_lnK;

    // green
    if (ptsz->unwise_galaxy_sample_id == 1 || ptsz->unwise_galaxy_sample_id == 2)
    this_lnI = this_lnJ;

    // blue
    //if (ptsz->unwise_galaxy_sample_id == 3)
    //this_lnI = this_lnI;

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

  /** 3. Store the read results into CLASS structures */
  ptsz->normalized_fdndz_size = n_data;
  /** Make room */

  class_realloc(ptsz->normalized_fdndz_z,
                ptsz->normalized_fdndz_z,
                ptsz->normalized_fdndz_size*sizeof(double),
                ptsz->error_message);
  class_realloc(ptsz->normalized_fdndz_phig,
                ptsz->normalized_fdndz_phig,
                ptsz->normalized_fdndz_size*sizeof(double),
                ptsz->error_message);



  /** Store them */
  for (index_x=0; index_x<ptsz->normalized_fdndz_size; index_x++) {
    ptsz->normalized_fdndz_z[index_x] = lnx[index_x];
    ptsz->normalized_fdndz_phig[index_x] = lnI[index_x];
    //print("z=%.3e phig=%.3e\n",ptsz->normalized_dndz_z[index_x])
  };

  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}

// unwise dndz deduced from cross-match with spectroscopic surveys
int load_normalized_cosmos_dndz(struct tszspectrum * ptsz)
{
  // printf("gs = %d\n",ptsz->galaxy_sample);
// don't load the unwise  dndz  if none of the following are required:
if (   (ptsz->has_tSZ_gal_1h != _TRUE_ )
    && (ptsz->has_tSZ_gal_2h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_1h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_1h_fft != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_2h_fft != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_3h_fft != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_2h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_3h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_hf != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_lensmag_1halo != _TRUE_ )
    && (ptsz->has_gal_gal_1h != _TRUE_ )
    && (ptsz->has_gal_lens_1h != _TRUE_ )
    && (ptsz->has_gal_lens_2h != _TRUE_ )
    && (ptsz->has_gal_cib_1h != _TRUE_ )
    && (ptsz->has_gal_cib_2h != _TRUE_ )
    && (ptsz->has_gal_lens_hf != _TRUE_ )
    && (ptsz->has_gal_lensmag_1h != _TRUE_ )
    && (ptsz->has_gal_lensmag_2h != _TRUE_ )
    && (ptsz->has_gal_lensmag_hf != _TRUE_ )
    && (ptsz->has_tSZ_lensmag_1h != _TRUE_ )
    && (ptsz->has_tSZ_lensmag_2h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_1h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_2h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_hf != _TRUE_ )
    && (ptsz->has_lens_lensmag_1h != _TRUE_ )
    && (ptsz->has_lens_lensmag_2h != _TRUE_ )
    && (ptsz->has_lens_lensmag_hf != _TRUE_ )
    && (ptsz->has_gal_gal_2h != _TRUE_ )
    && (ptsz->has_gal_gal_hf != _TRUE_ ))
  return 0;



if ((ptsz->galaxy_sample == 0) || (ptsz->galaxy_sample == 2))
  return 0;

if (ptsz->sz_verbose >= 1)
printf("-> Loading cosmos dndz unwise\n");

  class_alloc(ptsz->normalized_cosmos_dndz_z,sizeof(double *)*100,ptsz->error_message);
  class_alloc(ptsz->normalized_cosmos_dndz_phig,sizeof(double *)*100,ptsz->error_message);

  //class_alloc(ptsz->PP_d2lnI,sizeof(double *)*100,ptsz->error_message);

  //char arguments[_ARGUMENT_LENGTH_MAX_];
  char line[_LINE_LENGTH_MAX_];
  //char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  FILE *process;
  int n_data_guess, n_data = 0;
  double *lnx = NULL, *lnI = NULL,  *tmp = NULL;
  double this_lnx, this_lnI, this_lnJ, this_lnK;
  int status;
  int index_x;


  /** 1. Initialization */
  /* Prepare the data (with some initial size) */
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));
  lnI = (double *)malloc(n_data_guess*sizeof(double));



  /* Prepare the command */
  /* If the command is just a "cat", no arguments need to be passed */
  // if(strncmp("cat ", ptsz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];

  //unwise
    class_open(process,ptsz->UNWISE_cosmos_dndz_file, "r",ptsz->error_message);
    if (ptsz->sz_verbose > 0){
      printf("-> Openning the cosmos dndz file for unWISE galaxies\n");
      printf("-> File Name: %s\n",ptsz->UNWISE_cosmos_dndz_file);
      // printf("-> File Name: %s\n",ptsz->UNWISE_fdndz_file);
      // printf("-> File Name: %s\n",ptsz->A10_file);
    }

    // sprintf(Filepath,
    //         "%s%s",
    //         "cat ",
    //         //ptsz->path_to_class,
    //         "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/UNWISE_galaxy_ditributions/normalised_dndz_cosmos.txt");
    //


  // process = popen(Filepath, "r");

  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {

    // unWISE load and read column depending on the requested color

    sscanf(line, "%lf %lf %lf %lf", &this_lnx, &this_lnI, &this_lnJ, &this_lnK);
    // sscanf(line, "%lf %lf", &this_lnx, &this_lnI);

    // red
    if (ptsz->unwise_galaxy_sample_id == 0)
    this_lnI = this_lnK;

    // green
    if (ptsz->unwise_galaxy_sample_id == 1 || ptsz->unwise_galaxy_sample_id == 2)
    this_lnI = this_lnJ;

    // blue
    //if (ptsz->unwise_galaxy_sample_id == 3)
    //this_lnI = this_lnI;

    // printf("lnx = %e\n",this_lnI);





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

  // printf("closing the process\n");

  /* Close the process */
  status = fclose(process);
  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  ptsz->normalized_cosmos_dndz_size = n_data;
  /** Make room */

  class_realloc(ptsz->normalized_cosmos_dndz_z,
                ptsz->normalized_cosmos_dndz_z,
                ptsz->normalized_cosmos_dndz_size*sizeof(double),
                ptsz->error_message);
  class_realloc(ptsz->normalized_cosmos_dndz_phig,
                ptsz->normalized_cosmos_dndz_phig,
                ptsz->normalized_cosmos_dndz_size*sizeof(double),
                ptsz->error_message);


  // printf("fillling the arrays of dim = %d\n",ptsz->normalized_cosmos_dndz_size);
  /** Store them */
  for (index_x=0; index_x<ptsz->normalized_cosmos_dndz_size; index_x++) {
    ptsz->normalized_cosmos_dndz_z[index_x] = lnx[index_x];
    ptsz->normalized_cosmos_dndz_phig[index_x] = lnI[index_x];
    // printf("z=%.3e phig=%.3e\n",ptsz->normalized_dndz_z[index_x]);
  };
// exit(0);
  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}




int load_ksz_filter(struct tszspectrum * ptsz)
{

  if (ptsz->sz_verbose >= 1)
    printf("-> loading the filter f(l) for cl^kSZ2_gal\n");


  class_alloc(ptsz->l_unwise_filter,sizeof(double *)*100,ptsz->error_message);
  class_alloc(ptsz->f_unwise_filter,sizeof(double *)*100,ptsz->error_message);
  //class_alloc(ptsz->PP_d2lnI,sizeof(double *)*100,ptsz->error_message);

  //char arguments[_ARGUMENT_LENGTH_MAX_];
  char line[_LINE_LENGTH_MAX_];
  //char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  FILE *process;
  int n_data_guess, n_data = 0;
  double *lnx = NULL, *lnI = NULL,  *tmp = NULL;
  double this_lnx, this_lnI;
  int status;
  int index_x;


  /** 1. Initialization */
  /* Prepare the data (with some initial size) */
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));
  lnI = (double *)malloc(n_data_guess*sizeof(double));


  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];

  class_open(process,ptsz->ksz_filter_file, "r",ptsz->error_message);
  printf("-> File Name: %s\n",ptsz->ksz_filter_file);


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

  /** 3. Store the read results into CLASS structures */
  ptsz->unwise_filter_size = n_data;
  /** Make room */

  class_realloc(ptsz->l_unwise_filter,
                ptsz->l_unwise_filter,
                ptsz->unwise_filter_size*sizeof(double),
                ptsz->error_message);
  class_realloc(ptsz->f_unwise_filter,
                ptsz->f_unwise_filter,
                ptsz->unwise_filter_size*sizeof(double),
                ptsz->error_message);



  /** Store them */
  for (index_x=0; index_x<ptsz->unwise_filter_size; index_x++) {
    ptsz->l_unwise_filter[index_x] = lnx[index_x];
    ptsz->f_unwise_filter[index_x] = lnI[index_x];
  };

  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}




//This routine reads the tabulated
//alpha(z) normalisation for Tinker et al 2010 HMF
//and stores the tabulated values.

int load_T10_alpha_norm(struct tszspectrum * ptsz)
{


  class_alloc(ptsz->T10_ln1pz,sizeof(double *)*100,ptsz->error_message);
  class_alloc(ptsz->T10_lnalpha,sizeof(double *)*100,ptsz->error_message);
  //class_alloc(ptsz->PP_d2lnI,sizeof(double *)*100,ptsz->error_message);

  //char arguments[_ARGUMENT_LENGTH_MAX_];
  char line[_LINE_LENGTH_MAX_];
  //char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_];
  FILE *process;
  int n_data_guess, n_data = 0;
  double *lnx = NULL, *lnI = NULL,  *tmp = NULL;
  double this_lnx, this_lnI;
  int status;
  int index_x;


  /** 1. Initialization */
  /* Prepare the data (with some initial size) */
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));
  lnI = (double *)malloc(n_data_guess*sizeof(double));

  // class_open(process,"sz_auxiliary_files/Tinker_et_al_10_alpha_consistency_msyriac.txt", "r",ptsz->error_message);
  class_open(process,ptsz->Tinker_et_al_10_alpha_consistency_msyriac_file, "r",ptsz->error_message);


  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf %lf", &this_lnx, &this_lnI);
    //printf("lnx = %e  lnI = %e \n",this_lnx,this_lnI);




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

  /** 3. Store the read results into CLASS structures */
  ptsz->T10_lnalpha_size = n_data;
  /** Make room */

  class_realloc(ptsz->T10_ln1pz,
                ptsz->T10_ln1pz,
                ptsz->T10_lnalpha_size*sizeof(double),
                ptsz->error_message);
  class_realloc(ptsz->T10_lnalpha,
                ptsz->T10_lnalpha,
                ptsz->T10_lnalpha_size*sizeof(double),
                ptsz->error_message);



  /** Store them */
  for (index_x=0; index_x<ptsz->T10_lnalpha_size; index_x++) {
    ptsz->T10_ln1pz[index_x] = log(1.+lnx[index_x]);
    ptsz->T10_lnalpha[index_x] = log(lnI[index_x]);
  };

  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}



//HMF Tinker et al 2008
//interpolated at mdeltac
int MF_T08_m500(
                double * result,
                double * lognu ,
                double z ,
                double delta_crit,
                struct tszspectrum * ptsz
                )
{
  //T08@m500
  if(z>3.) z=3.;
  double om0 = ptsz->Omega_m_0;
  double ol0 = 1.-ptsz->Omega_m_0;
  double Omega_m_z = om0*pow(1.+z,3.)/(om0*pow(1.+z,3.)+ ol0);
  double  delta_mean = delta_crit/Omega_m_z;

  delta_mean = log10(delta_mean);

  double delta_mean_tab[9]=
  {
    200.,
    300.,
    400.,
    600.,
    800.,
    1200.,
    1600.,
    2400.,
    3200.
  };

  int i;
  for (i=0;i<9;i++)
    delta_mean_tab[i] =
    log10(delta_mean_tab[i]);

  double A_tab[9] =
  {
    0.186,
    0.200,
    0.212,
    0.218,
    0.248,
    0.255,
    0.260,
    0.260,
    0.260
  };

  double aa_tab[9] =
  {
    1.47,
    1.52,
    1.56,
    1.61,
    1.87,
    2.13,
    2.30,
    2.53,
    2.66
  };

  double b_tab[9] =
  {
    2.57,
    2.25,
    2.05,
    1.87,
    1.59,
    1.51,
    1.46,
    1.44,
    1.41
  };

  double c_tab[9] =
  {
    1.19,
    1.27,
    1.34,
    1.45,
    1.58,
    1.80,
    1.97,
    2.24,
    2.44
  };

  //#Table 3 of Tinker 2008 : 2nd derivatives
  double d2_A_tab[9] =
  {
    0.00,
    0.50,
    -1.56,
    3.05,
    -2.95,
    1.07,
    -0.71,
    0.21,
    0.00
  };

  double d2_aa_tab[9] =
  {
    0.00,
    1.19,
    -6.34,
    21.36,
    -10.95,
    2.59,
    -0.85,
    -2.07,
    0.00
  };

  double d2_b_tab[9] =
  {
    0.00,
    -1.08,
    12.61,
    -20.96,
    24.08,
    -6.64,
    3.84,
    -2.09,
    0.00
  };

  double d2_c_tab[9] =
  {
    0.00,
    0.94,
    -0.43,
    4.61,
    0.01,
    1.21,
    1.43,
    0.33,
    0.00
  };

  double * Ap0,* a0,* b0,* c0;

  class_alloc(Ap0,
              1*sizeof(double),
              ptsz->error_message);
  class_alloc(a0,
              1*sizeof(double),
              ptsz->error_message);
  class_alloc(b0,
              1*sizeof(double),
              ptsz->error_message);
  class_alloc(c0,
              1*sizeof(double),
              ptsz->error_message);

  splint(delta_mean_tab,
         A_tab,
         d2_A_tab,
         9,
         delta_mean,
         Ap0);

  splint(delta_mean_tab,
         aa_tab,
         d2_aa_tab,
         9,
         delta_mean,
         a0);

  splint(delta_mean_tab,
         b_tab,
         d2_b_tab,
         9,
         delta_mean,
         b0);

  splint(delta_mean_tab,
         c_tab,
         d2_c_tab,
         9,
         delta_mean,
         c0);

  double alphaT08 =
  pow(10.,-pow(0.75/log10(pow(10.,delta_mean)/75.),1.2));
  double   Ap=*Ap0*pow(1.+z,-0.14);
  double   a=*a0*pow(1.+z,-0.06);
  double   b=*b0*pow(1.+z,-alphaT08);
  double   c=*c0;
  double   nu= exp(*lognu);
  double sigma= 1.6865/sqrt(nu);

  free(Ap0);
  free(a0);
  free(b0);
  free(c0);

  //sigma = 1.;

  *result = 0.5*(Ap*(pow(sigma/b,-a)+1.)*exp(-c/pow(sigma,2.)));
  //*result = 0.5*(Ap*(pow(1.0/b,-a)+1.)*exp(-c/pow(1.0,2.)));

  return _SUCCESS_;
}



int plc_gnfw(double * plc_gnfw_x,
             double x ,
             double * pvectsz,
             struct background * pba,
             struct tszspectrum * ptsz)
{
  //Custom. GNFW pressure profile
  //int index_l = (int) pvectsz[ptsz->index_multipole_for_pressure_profile];
  int index_md = (int) pvectsz[ptsz->index_md];
      *plc_gnfw_x = 0.;

  // Example Arnaud 2010
  // ptsz->P0GNFW = 8.130;
  // ptsz->c500 = 1.156;
  // ptsz->gammaGNFW = 0.3292;
  // ptsz->alphaGNFW = 1.0620;
  // ptsz->betaGNFW = 5.4807;

  //Custom. GNFW
  //if(ptsz->pressure_profile == 3){
      *plc_gnfw_x = (1./(pow(ptsz->c500*x,ptsz->gammaGNFW)
                    *pow(1.+ pow(ptsz->c500*x,ptsz->alphaGNFW),
                        (ptsz->betaGNFW-ptsz->gammaGNFW)/ptsz->alphaGNFW)))
                    *pow(x,2)
                    //*sin(x*(ptsz->ell[index_l]+0.5)/pvectsz[ptsz->index_l500])
                    /(x*(pvectsz[ptsz->index_multipole_for_pressure_profile]+0.5)/pvectsz[ptsz->index_l500]);


    if (_mean_y_ || _dydz_)
      *plc_gnfw_x = (1./(pow(ptsz->c500*x,ptsz->gammaGNFW)
                       *pow(1.+ pow(ptsz->c500*x,ptsz->alphaGNFW),
                           (ptsz->betaGNFW-ptsz->gammaGNFW)/ptsz->alphaGNFW)))
                      *pow(x,2);

  //}



  //Battaglia et al 2012 pressure profile
  //Eq. 10
  if(ptsz->pressure_profile == 4){
    //ptsz->P0_B12 = 18.1;
    //ptsz->xc_B12 = 0.497;
    //ptsz->beta_B12 = 4.35;

    //ptsz->alpha_m_P0_B12 = 0.154;
    //ptsz->alpha_m_xc_B12 = -0.00865;
    //ptsz->alpha_m_beta_B12 = 0.0393;

    //ptsz->alpha_z_P0_B12 = -0.758;
    //ptsz->alpha_z_xc_B12 = 0.731;
    //ptsz->alpha_z_beta_B12 = 0.415;



    double xc;
    double beta;
    double P0;

    double m200_over_msol = pvectsz[ptsz->index_m200c]/pba->h; // convert to Msun
    double z = pvectsz[ptsz->index_z];


    P0 = ptsz->P0_B12*pow(m200_over_msol/1e14,ptsz->alpha_m_P0_B12)*pow(1+z,ptsz->alpha_z_P0_B12);
    xc = ptsz->xc_B12*pow(m200_over_msol/1e14,ptsz->alpha_m_xc_B12)*pow(1+z,ptsz->alpha_z_xc_B12);
    beta = ptsz->beta_B12*pow(m200_over_msol/1e14,ptsz->alpha_m_beta_B12)*pow(1+z,ptsz->alpha_z_beta_B12);

    double gamma = -0.3;
    double alpha = 1.0;

      *plc_gnfw_x = P0*pow(x/xc,gamma)*pow(1.+ pow(x/xc,alpha),-beta)
                    *pow(x,2)
                    /(x*(pvectsz[ptsz->index_multipole_for_pressure_profile]+0.5)/pvectsz[ptsz->index_l200c]);
    if (_mean_y_ || _dydz_)
      *plc_gnfw_x = P0*pow(x/xc,gamma)*pow(1.+ pow(x/xc,alpha),-beta)*pow(x,2);



  }
  return _SUCCESS_;
}



// HMF Tinker 2010
// https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumdks/szfastdks/mf_tinker10.f90
int MF_T10 (
            double * result,
            double * lognu ,
            double z ,
            struct tszspectrum * ptsz
            )
{
  *result = get_f_tinker10_at_nu_and_z(exp(*lognu),z,ptsz->hm_consistency,ptsz);
  return _SUCCESS_;
}


double get_f_tinker10_at_nu_and_z(double nu, double z, int hm_consistency, struct tszspectrum * ptsz){
  if(z>3.) z=3.;

  double alpha;

  // always do alpha(z)

  alpha = get_T10_alpha_at_z(z,ptsz);

  double lognu = log(nu);


  double result = 0.5
            *alpha
            *(1.+pow(pow(ptsz->beta0SZ*pow(1.+z,0.2),2.)
                   *exp(lognu),
                   -ptsz->phi0SZ
                   *pow(1.+z,-0.08)))*pow(
                    exp(lognu),
                    ptsz->eta0SZ
                    *pow(1.+z,0.27))
            *exp(-ptsz->gamma0SZ
                 *pow(1.+z,-0.01)
                 *exp(lognu)/2.)
            *sqrt(exp(lognu));
  return result;


}



//HMF Tinker et al 2008
//@ M200m
int MF_T08(
           double * result,
           double * lognu ,
           double z ,
           struct tszspectrum * ptsz
           )
{
  double alphaT08 = pow(10.,-pow(0.75/log10(200./75.),1.2));

  double   Ap=0.186*pow(1.+z,-0.14);
  double   a=1.47*pow(1.+z,-0.06);
  double   b=2.57*pow(1.+z,-alphaT08);
  double   c=1.19;
  double   nu= exp(*lognu);
  double sigma= 1.6865/sqrt(nu);

  *result = 0.5*(Ap*(pow(sigma/b,-a)+1.)*exp(-c/pow(sigma,2.)));

  return _SUCCESS_;
}



//HMF Tinker et al 2008
//@ M1600m
int MF_T08_M1600m(
                  double * result,
                  double * lognu ,
                  double z ,
                  struct tszspectrum * ptsz
                  )
{
  double alphaT08 = pow(10.,-pow(0.75/log10(1600./75.),1.2));
  double   Ap=0.260*pow(1.+z,-0.14);
  double   a=2.30*pow(1.+z,-0.06);
  double   b=1.46*pow(1.+z,-alphaT08);
  double   c=1.97;
  double   nu= exp(*lognu);
  double sigma= 1.6865/sqrt(nu);

  *result = 0.5*(Ap*(pow(sigma/b,-a)+1.)*exp(-c/pow(sigma,2.)));
  return _SUCCESS_;
}


//HMF Boquet et al 2015
int MF_B15(
           double * result,
           double * lognu ,
           double z ,
           struct tszspectrum * ptsz
           )
{
  //B15
  double   Ap=ptsz->Ap0*pow(1.+z,0.285);
  double   a=ptsz->a0*pow(1.+z,-0.058);
  double   b=ptsz->b0*pow(1.+z,-0.366);
  double   c=ptsz->c0*pow(1.+z,-0.045);
  double   nu= exp(*lognu);
  double sigma= 1.6865/sqrt(nu);

  *result = 0.5*(Ap*(pow(sigma/b,-a)+1.)*exp(-c/pow(sigma,2.)));

  return _SUCCESS_;
}


//HMF Boquet et al 2015 @ M500c
//TBC (11.04.19)
int MF_B15_M500c(double * result,
                 double * lognu ,
                 double z ,
                 struct tszspectrum * ptsz)
{
  //B15

  double   Ap=0.180*pow(1.+z,1.088);
  double   a=2.29*pow(1.+z,0.150);
  double   b=2.44*pow(1.+z,-1.008);
  double   c=1.97*pow(1.+z,-0.322);
  double   nu= exp(*lognu);
  double sigma= 1.6865/sqrt(nu);

  *result = 0.5*(Ap*(pow(sigma/b,-a)+1.)*exp(-c/pow(sigma,2.)));

  return _SUCCESS_;
}



int MF_J01(double * result,
           double * lognu ,
           struct tszspectrum * ptsz)
{
  double   nu= exp(*lognu);
  double sigma= 1.6865/sqrt(nu);

  *result =
  0.5
  *0.301
  *exp(-pow(fabs(log(1./sigma)+0.64),3.82));

  return _SUCCESS_;
}

double erf_compl_ps(double y,
                    double sn,
                    double q){
  //Completeness with error function
  double arg = (y - q * sn)/(sqrt(2.) * sn);
  double erf_compl = (erf(arg) + 1.)/2.;
  return erf_compl;
}

double next_z(double z_i, double dz, struct tszspectrum * ptsz){
  // Compute redshift bins where bins are defined with higher resolution for z<0.2
  double dz_i;
  double highres_z_cutoff = ptsz->cluster_count_completeness_grid_z_cutoff_low;
  if (z_i < highres_z_cutoff)
    dz_i = ptsz->dz_cluster_count_completeness_grid_low_z;
  //else if ((z_i >= highres_z_cutoff) && (z_i <= 1.))
  else if ((z_i >= highres_z_cutoff) && (z_i <= ptsz->cluster_count_completeness_grid_z_cutoff_mid))
    dz_i = ptsz->dz_cluster_count_completeness_grid_mid_z;
  else
    dz_i = ptsz->dz_cluster_count_completeness_grid_high_z;;

  double next_z = z_i + dz_i;
  return next_z;
}


double erf_compl(double y,
                 double sn,
                 double q){
  //Completeness with error function
  double arg = (y - q * sn)/(sqrt(2.) * sn);
  double erf_compl = (erf(arg) + 1.)/2.;
  return erf_compl;
}

double d_erf_compl_dq(double y,
                      double sn,
                      double q){
  //Completeness with error function
  double arg = (y - q * sn)/(sqrt(2.) * sn);
  double erf_compl = (y/sn)*exp(-arg*arg)/sqrt(2.*_PI_);
  return erf_compl;
}


double Delta_c_of_Omega_m(double Omega_m){
double Delta_c = 18.*pow(_PI_,2) + 82.*(Omega_m-1.) - 39.*pow((Omega_m-1.),2);
return Delta_c;
}

struct Parameters_for_integrand_redshift{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvecback;
  double * pvectsz;
  //double * llprime_grid;
};

double integrand_redshift(double ln1pz, void *p){

  struct Parameters_for_integrand_redshift *V = ((struct Parameters_for_integrand_redshift *) p);

   V->pvectsz[V->ptsz->index_z] = exp(ln1pz)-1.;

  double z =  V->pvectsz[V->ptsz->index_z];

  //Evaluation of background quantities @ z:
  double tau;
  int first_index_back = 0;

  class_call(background_tau_of_z(V->pba,z,&tau),
             V->ptsz->error_message,
             V->ptsz->error_message);

  class_call(background_at_tau(V->pba,
                               tau,
                               V->pba->long_info,
                               V->pba->inter_normal,
                               &first_index_back,
                               V->pvecback),
             V->ptsz->error_message,
             V->ptsz->error_message);


  //volume element in units h^-3 Mpc^3
  //volume = dv/(dzdOmega)*(c/H)
  // Chi^2 dChi = dV/(dzdOmega)*(c/H) dz

  V->pvectsz[V->ptsz->index_volume] = pow(1.+z,2)
                                      *pow(V->pvecback[V->pba->index_bg_ang_distance]*V->pba->h,2)
                                      *_c_*1.e-5
                                      /(V->pvecback[V->pba->index_bg_H]/V->pba->H0);


  V->pvectsz[V->ptsz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2); // conformal distance squared in [Mpc/h]^2

  V->pvectsz[V->ptsz->index_dgdz] = V->pvecback[V->pba->index_bg_D]*(1.-V->pvecback[V->pba->index_bg_f]); // d/dz(D/a)


  int index_md = (int) V->pvectsz[V->ptsz->index_md];

  if (((V->ptsz->has_kSZ_kSZ_lensmag_1halo == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_lensmag_1halo))
    ||((V->ptsz->has_tSZ_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_lensmag_1h))
    ||((V->ptsz->has_tSZ_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_lensmag_2h))
    ||((V->ptsz->has_gal_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_1h))
    ||((V->ptsz->has_gal_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_2h))
    ||((V->ptsz->has_lensmag_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_lensmag_lensmag_1h))
    ||((V->ptsz->has_lensmag_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lensmag_lensmag_2h))
){
    // compute kernel for lensing magnification
    // lensing of galaxies

    evaluate_redshift_int_lensmag(V->pvectsz,V->ptsz);
    double redshift_int_lensmag = V->pvectsz[V->ptsz->index_W_lensmag];
    double chi = sqrt(V->pvectsz[V->ptsz->index_chi2]);
    V->pvectsz[V->ptsz->index_lensing_Sigma_crit] = pow(3.*pow(V->pba->H0/V->pba->h,2)/2./V->ptsz->Rho_crit_0,-1)*pow((1.+z),1.)/(chi*redshift_int_lensmag);

  }
  else if (
      ((V->ptsz->has_lens_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_1h))
    ||((V->ptsz->has_lens_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_2h))
  )
  {

double chi = sqrt(V->pvectsz[V->ptsz->index_chi2]);
double chi_star =  V->ptsz->chi_star;  // in Mpc/h
// sigma_crit_lensmag:
evaluate_redshift_int_lensmag(V->pvectsz,V->ptsz);
double redshift_int_lensmag = V->pvectsz[V->ptsz->index_W_lensmag];

double sigma_crit_lensmag = pow(3.*pow(V->pba->H0/V->pba->h,2)/2./V->ptsz->Rho_crit_0,-1)*pow((1.+z),1.)/(chi*redshift_int_lensmag);
// sigma crit kappa:
double sigma_crit_kappa =  pow(3.*pow(V->pba->H0/V->pba->h,2)/2./V->ptsz->Rho_crit_0,-1)*pow((1.+z),1.)*chi_star/chi/(chi_star-chi);

V->pvectsz[V->ptsz->index_lensing_Sigma_crit] =  sqrt(sigma_crit_lensmag*sigma_crit_kappa);

  }
  else {

// CMB lensing
// Eq. 6 of https://arxiv.org/pdf/1312.4525.pdf

double chi = sqrt(V->pvectsz[V->ptsz->index_chi2]);
double chi_star =  V->ptsz->chi_star;  // in Mpc/h
double sigma_crit_kappa =  pow(3.*pow(V->pba->H0/V->pba->h,2)/2./V->ptsz->Rho_crit_0,-1)*pow((1.+z),1.)*chi_star/chi/(chi_star-chi);
V->pvectsz[V->ptsz->index_lensing_Sigma_crit] = sigma_crit_kappa;


  }

  // critical density in (Msun/h)/(Mpc/h)^3
  V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                        *pow(_Mpc_over_m_,1)
                                        *pow(_c_,2)
                                        *V->pvecback[V->pba->index_bg_rho_crit]
                                        /pow(V->pba->h,2);

  double Eh = V->pvecback[V->pba->index_bg_H]/V->pba->H0;
  double omega = V->pvecback[V->pba->index_bg_Omega_m];//pow(Eh,2.);
  V->pvectsz[V->ptsz->index_Delta_c] = Delta_c_of_Omega_m(omega);





if (((V->ptsz->has_tSZ_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_gal_1h))
     || ((V->ptsz->has_tSZ_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_gal_2h))
     || ((V->ptsz->has_gal_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_1h))
     || ((V->ptsz->has_gal_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_2h))
     || ((V->ptsz->has_pk_gg_at_z_1h == _TRUE_) && (index_md == V->ptsz->index_md_pk_gg_at_z_1h))
     || ((V->ptsz->has_pk_gg_at_z_2h == _TRUE_) && (index_md == V->ptsz->index_md_pk_gg_at_z_2h))
     || ((V->ptsz->has_gal_lens_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_1h))
     || ((V->ptsz->has_gal_lens_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_2h))
     || ((V->ptsz->has_gal_cib_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_cib_1h))
     || ((V->ptsz->has_gal_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_cib_2h))
     || ((V->ptsz->has_gal_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_1h))
     || ((V->ptsz->has_gal_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_2h))
     || ((V->ptsz->has_lensmag_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_lensmag_lensmag_1h))
     || ((V->ptsz->has_lensmag_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lensmag_lensmag_2h))
     || ((V->ptsz->has_lens_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_1h))
     || ((V->ptsz->has_lens_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_2h))
     //|| ((V->ptsz->has_kSZ_kSZ_lensmag_1halo == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_lensmag_1halo))
     || ((V->ptsz->has_kSZ_kSZ_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_1h))
     || ((V->ptsz->has_kSZ_kSZ_gal_1h_fft == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_1h_fft))
     || ((V->ptsz->has_kSZ_kSZ_gal_2h_fft == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_2h_fft))
     || ((V->ptsz->has_kSZ_kSZ_gal_3h_fft == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_3h_fft))
     || ((V->ptsz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_2h))
     || ((V->ptsz->has_kSZ_kSZ_gal_3h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_3h))
     // || ((V->ptsz->has_kSZ_kSZ_gal_hf == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_hf))

    ) {

 V->pvectsz[V->ptsz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->ptsz);
}


  double result = 0.;

  // // first deal with quantities that does not require mass integration:
  // if ((V->ptsz->has_pk_at_z_2h == _TRUE_) && (index_md == V->ptsz->index_md_pk_at_z_2h)) {
  //
  //   // eq. 10 of https://arxiv.org/pdf/1505.07833.pdf
  //   result = 1.;
  // }
  // else
  if ((V->ptsz->has_gal_gal_hf == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_hf)) {
    int index_l = (int) V->pvectsz[V->ptsz->index_multipole];
    double l = V->ptsz->ell[index_l];
    V->pvectsz[V->ptsz->index_multipole_for_pk] = l; // l1,l2 or l3
    evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
    double pk1 = V->pvectsz[V->ptsz->index_pk_for_halo_bias];
    result = pk1;
    evaluate_effective_galaxy_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
    result *= V->pvectsz[V->ptsz->index_halo_bias]*V->pvectsz[V->ptsz->index_halo_bias];



  }
  else if ((V->ptsz->has_kSZ_kSZ_gal_hf == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_hf)) {
  //
  // int index_l_1 = (int) V->pvectsz[V->ptsz->index_multipole_1];
  // int index_l_2 = (int) V->pvectsz[V->ptsz->index_multipole_2];
  // int index_l_3 = (int) V->pvectsz[V->ptsz->index_multipole_3];

  int index_theta_1 = (int) V->pvectsz[V->ptsz->index_multipole_1];
  int index_l_2 = (int) V->pvectsz[V->ptsz->index_multipole_2];
  int index_l_3 = (int) V->pvectsz[V->ptsz->index_multipole_3];

  // printf("index_l_1 = %d index_l_2 = %d index_l_3 = %d\n",index_theta_1,index_l_2,index_l_3);


  double l2 = exp(V->ptsz->ell_kSZ2_gal_multipole_grid[index_l_2]);
  double l3 = V->ptsz->ell[index_l_3];
  double ell = l3;
  double ell_prime = l2;
  double theta_1 = V->ptsz->theta_kSZ2_gal_theta_grid[index_theta_1];
  double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos(theta_1));
  // double l1 = sqrt(ell*ell+ell_prime*ell_prime-2.*ell*ell_prime*theta_1);


  //
  // double l1 = V->pvectsz[V->ptsz->index_ell_1];
  // double l2 = V->pvectsz[V->ptsz->index_ell_2];
  // double l3 = V->pvectsz[V->ptsz->index_ell_3];

  // printf("l1 = %.3e l2 = %.3e l3 = %3.e\n",l1,l2,l3);


  double z = V->pvectsz[V->ptsz->index_z];
  double d_A = V->pvecback[V->pba->index_bg_ang_distance]*V->pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi

  double k1 = (l1+0.5)/d_A;
  double k2 = (l2+0.5)/d_A;
  double k3 = (l3+0.5)/d_A;

  // printf("k1 = %.3e k2 = %.3e k3 = %3.e\n",k1,k2,k3);


  V->pvectsz[V->ptsz->index_multipole_for_pk] = l1; // l1,l2 or l3
  evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
  double pk1 = V->pvectsz[V->ptsz->index_pk_for_halo_bias];

  V->pvectsz[V->ptsz->index_multipole_for_pk] = l2; // l1,l2 or l3
  evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
  double pk2 = V->pvectsz[V->ptsz->index_pk_for_halo_bias];

  V->pvectsz[V->ptsz->index_multipole_for_pk] = l3; // l1,l2 or l3
  evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
  double pk3 = V->pvectsz[V->ptsz->index_pk_for_halo_bias];

  // printf("pk1 = %.3e pk2 = %.3e pk3 = %3.e\n",pk1,pk2,pk3);


  double z_asked = log(1.+z);
  double R_asked = log(8./V->pba->h); //log(R) in Mpc
  double sigma8_at_z =  exp(pwl_interp_2d(
                             V->ptsz->n_arraySZ,
                             V->ptsz->ndimSZ,
                             V->ptsz->array_redshift,
                             V->ptsz->array_radius,
                             V->ptsz->array_sigma_at_z_and_R,
                             1,
                             &z_asked,
                             &R_asked));

  double knl = get_knl_at_z(z,V->ptsz);

  double n1 = get_nl_index_at_z_and_k(z,k1,V->ptsz,V->pnl);
  double n2 = get_nl_index_at_z_and_k(z,k2,V->ptsz,V->pnl);
  double n3 = get_nl_index_at_z_and_k(z,k3,V->ptsz,V->pnl);

  // printf("n1 = %.3e \t n2 = %.3e \t n3 = %.3e\n",n1,n2,n3);
  // printf("sigma8_at_z = %.3e \t knl = %.3e\n",sigma8_at_z,knl);

  double f2_eff_12 =  bispectrum_f2_kernel_eff(k1,k2,k3,n1,n2,sigma8_at_z,knl);
  double f2_eff_23 =  bispectrum_f2_kernel_eff(k2,k3,k1,n2,n3,sigma8_at_z,knl);
  double f2_eff_31 =  bispectrum_f2_kernel_eff(k3,k1,k2,n3,n1,sigma8_at_z,knl);
  // printf("f2_eff_12 = %.3e \t f2_eff_23 = %.3e \t f2_eff_31 = %.3e\n",
  // f2_eff_12,
  // f2_eff_23,
  // f2_eff_31);

  double b123 = 2.*pk1*pk2*f2_eff_12 + 2.*pk2*pk3*f2_eff_23 + 2.*pk3*pk1*f2_eff_31;


  double H_over_c_in_h_over_Mpc = V->pvecback[V->pba->index_bg_H]/V->pba->h;
  double galaxy_normalisation = 1.;//H_over_c_in_h_over_Mpc; // here is normally also the bias but we take it out and multiply afterward
  //galaxy_normalisation *= pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,-2.);///(1.+z);
  double sigmaT_over_mp = 8.305907197761162e-17 * pow(V->pba->h,2)/V->pba->h; // !this is sigmaT / m_prot in (Mpc/h)**2/(Msun/h)
  double a = 1. / (1. + z);
  double tau_fac = a*sigmaT_over_mp*V->pba->Omega0_b/V->ptsz->Omega_m_0/V->ptsz->mu_e*V->ptsz->f_free;
  double tau_normalisation =tau_fac
                            *pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,-2.)
                            *V->pvecback[V->pba->index_bg_Omega_m]*V->pvectsz[V->ptsz->index_Rho_crit];
                            // note: Omega_m*Rho_crit = rho_m(z)

  // simone ferraro's reio implementation:
  // double ya = pow(1./a,1.5);//**1.5
  // double zre = 10.3 ;    // # reionization  redshift
  // double are = 1. / (1 + zre);
  // double yre = pow(1./are,1.5);//**1.5
  // double deltay = 7.56;
  // double xe = (0.5*(1.-tanh((ya-yre)/deltay)));
  // tau_normalisation *= pow(1. + z,3) * xe;

  result = tau_normalisation*tau_normalisation*galaxy_normalisation*b123;
  // printf("result = %.3e tau_norm = %.3e gal_norm = %.3e b123 = %.3e\n",
  // result,
  // tau_normalisation,
  // galaxy_normalisation,
  // b123);

  if (isnan(b123) || isinf(b123) || result == 0.){
    // result = 1.e-300;
  printf("z = %.3e b123 = %.3e\n",z,b123);
  printf("k1 = %.3e k2 = %.3e k3 = %.3e\n",k1,k2,k3);
  printf("l1 = %.3e l2 = %.3e l3 = %.3e\n",l1,l2,l3);
  printf("pk1 = %.3e pk2 = %.3e pk3 = %.3e\n",pk1,pk2,pk3);
  printf("n1 = %.3e \t n2 = %.3e \t n3 = %.3e\n",n1,n2,n3);
  printf("\n\n");
}


  }
  else if ((V->ptsz->has_isw_lens == _TRUE_) && (index_md == V->ptsz->index_md_isw_lens)) {

  double delta_ell_lens =  delta_ell_lens_at_ell_and_z(V->pvecback,
                                                  V->pvectsz,
                                                  V->pba,
                                                  V->ppm,
                                                  V->pnl,
                                                  V->ptsz);

  double delta_ell_isw = delta_ell_isw_at_ell_and_z(V->pvecback,
                                                          V->pvectsz,
                                                          V->pba,
                                                          V->ppm,
                                                          V->pnl,
                                                          V->ptsz);
  result = delta_ell_lens*delta_ell_isw;

  }


  else if ((V->ptsz->has_isw_tsz == _TRUE_) && (index_md == V->ptsz->index_md_isw_tsz)){

  double delta_ell_isw = delta_ell_isw_at_ell_and_z(V->pvecback,
                                                          V->pvectsz,
                                                          V->pba,
                                                          V->ppm,
                                                          V->pnl,
                                                          V->ptsz);
  double delta_ell_y = integrate_over_m_at_z(V->pvecback,
                                              V->pvectsz,
                                              V->pba,
                                              V->pnl,
                                              V->ppm,
                                              V->ptsz);

  result = delta_ell_isw*delta_ell_y;


  }
  else if ((V->ptsz->has_isw_auto == _TRUE_) && (index_md == V->ptsz->index_md_isw_auto)){

  double delta_ell_isw = delta_ell_isw_at_ell_and_z(V->pvecback,
                                                    V->pvectsz,
                                                    V->pba,
                                                    V->ppm,
                                                    V->pnl,
                                                    V->ptsz);

  result = delta_ell_isw*delta_ell_isw;

  }

  // Halofit approach
  else if ((V->ptsz->has_gal_lens_hf == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_hf))
  {

  int index_l = (int) V->pvectsz[V->ptsz->index_multipole];
  double l = V->ptsz->ell[index_l];
  V->pvectsz[V->ptsz->index_multipole_for_pk] = l; // l1,l2 or l3
  evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
  double pk1 = V->pvectsz[V->ptsz->index_pk_for_halo_bias];
  result = pk1;
  evaluate_effective_galaxy_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);

  result *= V->pvectsz[V->ptsz->index_halo_bias];

  double W_lens =  radial_kernel_W_lensing_at_z(V->pvecback,
                                                V->pvectsz,
                                                V->pba,
                                                V->ppm,
                                                V->pnl,
                                                V->ptsz);
  // this is needed only in  the approximate calculation...
  // for the exact calculation in HOD, this comes out of Sigma_crit
  result *= W_lens;

  }
  // Halofit approach
  else if ((V->ptsz->has_gal_lensmag_hf == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_hf))
  {

//printf("ok\n");
  int index_l = (int) V->pvectsz[V->ptsz->index_multipole];
  double l = V->ptsz->ell[index_l];
  V->pvectsz[V->ptsz->index_multipole_for_pk] = l; // l1,l2 or l3
  evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
  double pk1 = V->pvectsz[V->ptsz->index_pk_for_halo_bias];
  result = pk1;
  evaluate_effective_galaxy_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);

  result *= V->pvectsz[V->ptsz->index_halo_bias];

  double W_lensmag =  radial_kernel_W_lensing_magnification_at_z(V->pvecback,
                                                                V->pvectsz,
                                                                V->pba,
                                                                V->ppm,
                                                                V->pnl,
                                                                V->ptsz);

    // this is needed only in  the approximate calculation
    // for the exact calculation in HOD, this comes out of Sigma_crit
    result *= W_lensmag;

  }
  // Halofit approach
else if ((V->ptsz->has_lensmag_lensmag_hf == _TRUE_) && (index_md == V->ptsz->index_md_lensmag_lensmag_hf))
  {

  int index_l = (int) V->pvectsz[V->ptsz->index_multipole];
  double l = V->ptsz->ell[index_l];
  V->pvectsz[V->ptsz->index_multipole_for_pk] = l; // l1,l2 or l3
  evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
  double pk1 = V->pvectsz[V->ptsz->index_pk_for_halo_bias];
  result = pk1;

  double W_lensmag =  radial_kernel_W_lensing_magnification_at_z(V->pvecback,
                                                                V->pvectsz,
                                                                V->pba,
                                                                V->ppm,
                                                                V->pnl,
                                                                V->ptsz);

    // this is needed only in  the approximate calculation
    // for the exact calculation in HOD, this comes out of Sigma_crit
    result *= W_lensmag*W_lensmag;

  }
  // Halofit approach
else if ((V->ptsz->has_lens_lensmag_hf == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_hf))
{

  int index_l = (int) V->pvectsz[V->ptsz->index_multipole];
  double l = V->ptsz->ell[index_l];
  V->pvectsz[V->ptsz->index_multipole_for_pk] = l; // l1,l2 or l3
  evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
  double pk1 = V->pvectsz[V->ptsz->index_pk_for_halo_bias];
  result = pk1;

  double W_lensmag =  radial_kernel_W_lensing_magnification_at_z(V->pvecback,
                                                                V->pvectsz,
                                                                V->pba,
                                                                V->ppm,
                                                                V->pnl,
                                                                V->ptsz);
  double W_lens =  radial_kernel_W_lensing_at_z(V->pvecback,
                                                  V->pvectsz,
                                                  V->pba,
                                                  V->ppm,
                                                  V->pnl,
                                                  V->ptsz);

    // this is needed only in  the approximate calculation
    // for the exact calculation in HOD, this comes out of Sigma_crit
    //printf("%.3e \t %.3e\n",W_lensmag,W_lens);
    result *= W_lensmag*W_lens;

  }

else if (
  ((V->ptsz->use_hod == 0) && (V->ptsz->has_lens_lens_2h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lens_2h))
  || ((V->ptsz->use_hod == 0) && (V->ptsz->has_lens_lens_1h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lens_1h))
){


if (index_md == V->ptsz->index_md_lens_lens_1h) {
  result = 0.;
}
else {
  double W_lens =  radial_kernel_W_lensing_at_z(V->pvecback,
                                                  V->pvectsz,
                                                  V->pba,
                                                  V->ppm,
                                                  V->pnl,
                                                  V->ptsz);
// this is needed only in  the approximate calculation
// for the exact calculation in halo model, this comes out of Sigma_crit
result = W_lens*W_lens;

}

}

  // then quantities that require mass integration
  else {
  // printf("integrating over mass\n");
  result = integrate_over_m_at_z(V->pvecback,
                                 V->pvectsz,
                                 V->pba,
                                 V->pnl,
                                 V->ppm,
                                 V->ptsz);
  // exit(0);
  }

if (((V->ptsz->has_sz_2halo == _TRUE_) && (index_md == V->ptsz->index_md_2halo))
 || ((V->ptsz->has_gal_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_2h)) //## BB debug
 || ((V->ptsz->has_cib_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_cib_cib_2h))
 || ((V->ptsz->has_tSZ_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_cib_2h))
 || ((V->ptsz->has_gal_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_cib_2h))
 || ((V->ptsz->has_lens_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_lens_cib_2h))
 || ((V->ptsz->has_tSZ_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_gal_2h))
 || ((V->ptsz->has_gal_lens_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_2h))
 || ((V->ptsz->has_gal_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_2h))
 || ((V->ptsz->has_lens_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_2h))
 || ((V->ptsz->has_lensmag_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lensmag_lensmag_2h))
 || ((V->ptsz->has_lens_lens_2h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lens_2h))
 || ((V->ptsz->has_tSZ_lens_2h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_lens_2h))
 || ((V->ptsz->has_tSZ_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_lensmag_2h))
 || ((V->ptsz->has_sz_m_y_y_2h == _TRUE_) && (index_md == V->ptsz->index_md_m_y_y_2h))
 || ((V->ptsz->has_pk_at_z_2h == _TRUE_) && (index_md == V->ptsz->index_md_pk_at_z_2h))
 || ((V->ptsz->has_pk_gg_at_z_2h == _TRUE_) && (index_md == V->ptsz->index_md_pk_gg_at_z_2h))
 // || ((V->ptsz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_2h))
    ){


  int index_l = (int) V->pvectsz[V->ptsz->index_multipole];
  V->pvectsz[V->ptsz->index_multipole_for_pk] = V->ptsz->ell[index_l];
  evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);

  // For all the above cases we multiply the linear matter power spectrum to the redshift integrand
  // evaluated at (ell+1/2)/Chi and redshift z
  result *= V->pvectsz[V->ptsz->index_pk_for_halo_bias];

}

if ((V->ptsz->has_bk_at_z_2h == _TRUE_) && (index_md == V->ptsz->index_md_bk_at_z_2h)){
  // printf("result 2h = %.3e\n",result);
  evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
  result *= 3.*V->pvectsz[V->ptsz->index_pk_for_halo_bias];
  // printf("result 2h = %.3e\n",result);
}

// if ((V->ptsz->has_bk_at_z_3h == _TRUE_) && (index_md == V->ptsz->index_md_bk_at_z_3h)){
//   double bh;
//   double pk;
//   double b0;
//   double f2;
//   int index_k = (int) V->pvectsz[V->ptsz->index_k_for_pk_hm];
//   double k = V->ptsz->k_for_pk_hm[index_k];
//   // printf("result 3h = %.3e\n",result);
//   evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
//   pk = V->pvectsz[V->ptsz->index_pk_for_halo_bias];
//   f2 = bispectrum_f2_kernel(k,k,k);
//   b0 = 3.*(2.*pk*pk*f2);
//   bh = (b0 + 3.*pk*pk);
//   result *= bh;
// }


// Power spectrum today : needed for ISW  stuff
if ( ((V->ptsz->has_isw_auto == _TRUE_) && (index_md == V->ptsz->index_md_isw_auto))
 ||  ((V->ptsz->has_isw_tsz == _TRUE_) && (index_md == V->ptsz->index_md_isw_tsz))
 ||  ((V->ptsz->has_isw_lens == _TRUE_) && (index_md == V->ptsz->index_md_isw_lens))
    ){

  evaluate_pk_at_ell_plus_one_half_over_chi_today(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);

  // For all the above cases we add the linear matter power spectrum to the redshift integrand
  // evaluated at (ell+1/2)/Chi and redshift z=0
  result *= V->pvectsz[V->ptsz->index_pk_for_halo_bias];


}

// galaxy radial kernel
if  (((V->ptsz->has_tSZ_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_gal_1h))
  || ((V->ptsz->has_tSZ_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_gal_2h))
  || ((V->ptsz->has_kSZ_kSZ_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_1h))
  || ((V->ptsz->has_kSZ_kSZ_gal_1h_fft == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_1h_fft))
  || ((V->ptsz->has_kSZ_kSZ_gal_2h_fft == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_2h_fft))
  || ((V->ptsz->has_kSZ_kSZ_gal_3h_fft == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_3h_fft))
  || ((V->ptsz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_2h))
  || ((V->ptsz->has_kSZ_kSZ_gal_3h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_3h))
  || ((V->ptsz->has_kSZ_kSZ_gal_hf == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_hf))
  || ((V->ptsz->has_gal_lens_hf == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_hf))
  || ((V->ptsz->has_gal_lens_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_2h))
  || ((V->ptsz->has_gal_lens_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_1h))
  || ((V->ptsz->has_gal_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_cib_2h))
  || ((V->ptsz->has_gal_cib_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_cib_1h))
  || ((V->ptsz->has_gal_lensmag_hf == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_hf))
  || ((V->ptsz->has_gal_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_2h))
  || ((V->ptsz->has_gal_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_1h))
    ){
// multiply by radial kernel for galaxies
double Wg = radial_kernel_W_galaxy_at_z(V->pvecback,V->pvectsz,V->pba,V->ptsz);
result *= Wg/V->pvectsz[V->ptsz->index_chi2];
}

// multiply by velocity dispersion
if ((V->ptsz->has_kSZ_kSZ_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_1h)
 || (V->ptsz->has_kSZ_kSZ_gal_1h_fft == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_1h_fft)
 || (V->ptsz->has_kSZ_kSZ_gal_2h_fft == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_2h_fft)
 || (V->ptsz->has_kSZ_kSZ_gal_3h_fft == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_3h_fft)
 || (V->ptsz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_2h)
 || (V->ptsz->has_kSZ_kSZ_gal_3h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_3h)
 || (V->ptsz->has_kSZ_kSZ_gal_hf == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_hf)
 || (V->ptsz->has_kSZ_kSZ_lensmag_1halo == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_lensmag_1halo)){
  // printf("evaluating vrms2\n");
  evaluate_vrms2(V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);
  result *= V->pvectsz[V->ptsz->index_vrms2]/3./pow(_c_*1e-3,2.);
}

// multiply by dsigma2_hsv
if ((V->ptsz->has_sz_cov_N_N_hsv == _TRUE_) && (index_md == V->ptsz->index_md_cov_N_N_hsv)){
  evaluate_sigma2_hsv(V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);
  result *= V->pvectsz[V->ptsz->index_sigma2_hsv];
}


// gxg needs Wg^2:
if ( ((V->ptsz->has_gal_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_1h))
   ||((V->ptsz->has_gal_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_2h))
   ||((V->ptsz->has_gal_gal_hf == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_hf))

  ){
// multiply by radial kernel for galaxies (squared for gxg quantities)
double Wg = radial_kernel_W_galaxy_at_z(V->pvecback,V->pvectsz,V->pba,V->ptsz);
result *= pow(Wg/V->pvectsz[V->ptsz->index_chi2],2.);
}


if  (((V->ptsz->has_cib_cib_1h == _TRUE_) && (index_md == V->ptsz->index_md_cib_cib_1h))
   ||((V->ptsz->has_cib_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_cib_cib_2h))
    ){
// cib redshift kernel, see McCarthy and Madhavacheril 2020
result *= 1./(1.+z)*1./(1.+z)*pow(1./V->pvectsz[V->ptsz->index_chi2],2.);
}

if  (((V->ptsz->has_tSZ_cib_1h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_cib_1h))
    ||((V->ptsz->has_tSZ_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_cib_2h))
    ||((V->ptsz->has_gal_cib_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_cib_1h))
    ||((V->ptsz->has_gal_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_cib_2h))
    ||((V->ptsz->has_lens_cib_1h == _TRUE_) && (index_md == V->ptsz->index_md_lens_cib_1h))
    ||((V->ptsz->has_lens_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_lens_cib_2h))
    ||((V->ptsz->has_cib_monopole == _TRUE_) && (index_md == V->ptsz->index_md_cib_monopole))
    ){

  result *= 1./(1.+z)*pow(1./V->pvectsz[V->ptsz->index_chi2],1.);
}

// if computing 3d matter power spectrum P(k):
if( ((V->ptsz->has_pk_at_z_1h == _TRUE_) && (index_md == V->ptsz->index_md_pk_at_z_1h))
 || ((V->ptsz->has_pk_at_z_2h == _TRUE_) && (index_md == V->ptsz->index_md_pk_at_z_2h))){
   // result *= pow((V->pba->Omega0_cdm+V->pba->Omega0_b)*V->ptsz->Rho_crit_0,-2);
   return result;
 }// end P(k)

if( ((V->ptsz->has_pk_gg_at_z_1h == _TRUE_) && (index_md == V->ptsz->index_md_pk_gg_at_z_1h))
 || ((V->ptsz->has_pk_gg_at_z_2h == _TRUE_) && (index_md == V->ptsz->index_md_pk_gg_at_z_2h))){
   return result;
 }// end P(k)
else if( ((V->ptsz->has_bk_at_z_1h == _TRUE_) && (index_md == V->ptsz->index_md_bk_at_z_1h))
      || ((V->ptsz->has_bk_at_z_2h == _TRUE_) && (index_md == V->ptsz->index_md_bk_at_z_2h))
      || ((V->ptsz->has_bk_at_z_3h == _TRUE_) && (index_md == V->ptsz->index_md_bk_at_z_3h))){
      // result *= pow((V->pba->Omega0_cdm+V->pba->Omega0_b)*V->ptsz->Rho_crit_0,-3);


   return result;
 }// end b(k,k,k)

else{

  // finally multiply by volume element Chi^2 dChi
  result *= V->pvectsz[V->ptsz->index_chi2];
  // printf("multiplying by volume\n");

  // integrate w.r.t ln(1+z); dz =  (1+z)dln(1+z)
  // volume element in units h^-3 Mpc^3
  // volume = dv/(dzdOmega)
  // Chi^2 dChi = dV/(dzdOmega) dz
  // Chi^2 dChi = dV/(dzdOmega)*(1+z) dln(1+z)
  // dChi = (c/H) *(1+z) dln(1+z) ---> this is used
  // dChi = (c/H) dz
  // d/dchi = H/c d/dz
  double H_over_c_in_h_over_Mpc = V->pvecback[V->pba->index_bg_H]/V->pba->h;
  result = (1.+V->pvectsz[V->ptsz->index_z])*result/H_over_c_in_h_over_Mpc;


  return result;
}

}


int integrate_over_redshift(struct background * pba,
                            struct nonlinear * pnl,
                            struct primordial * ppm,
                            struct tszspectrum * ptsz,
                            double * Pvecback,
                            double * Pvectsz)
{


  double z_min = ptsz->z1SZ;
  double z_max = ptsz->z2SZ;


  struct Parameters_for_integrand_redshift V;
  V.pnl = pnl;
  V.ppm = ppm;
  V.ptsz = ptsz;
  V.pba = pba;
  V.pvectsz = Pvectsz;
  V.pvecback = Pvecback;

  void * params = &V;
  double r; //result of the integral

  double epsrel= ptsz->redshift_epsrel;
  double epsabs= ptsz->redshift_epsabs;
  int show_neval = ptsz->patterson_show_neval;

  int index_md = (int) Pvectsz[ptsz->index_md];
if(_pk_at_z_1h_
|| _pk_at_z_2h_
|| _pk_gg_at_z_1h_
|| _pk_gg_at_z_2h_
|| _bk_at_z_1h_
|| _bk_at_z_2h_
|| _bk_at_z_3h_)
{
  r = integrand_redshift(log(1. + ptsz->z_for_pk_hm),params);
}
else{

  // printf("er = %.4e ea = %.4e\n",epsrel,epsabs);
  r = Integrate_using_Patterson_adaptive(log(1. + z_min), log(1. + z_max),
                                         epsrel, epsabs,
                                         integrand_redshift,
                                         params,show_neval);
    }

  Pvectsz[ptsz->index_integral] = r;

  return _SUCCESS_;
}





struct Parameters_for_integrand_patterson{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvecback;
  double * pvectsz;
};


//
double integrand_patterson_test(double logM, void *p){

  struct Parameters_for_integrand_patterson *V = ((struct Parameters_for_integrand_patterson *) p);

  double result = integrand_at_m_and_z(logM,
                                        V->pvecback,
                                        V->pvectsz,
                                        V->pba,
                                        V->ppm,
                                        V->pnl,
                                        V->ptsz);

  return result;

}



//Integration over the mass range at a given redshift
 double integrate_over_m_at_z(double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct nonlinear * pnl,
                             struct primordial * ppm,
                             struct tszspectrum * ptsz)
{

// if ( ((int) pvectsz[ptsz->index_md] == ptsz->index_md_dndlnM)
//   return integrand_at_m_and_z(log(1e16),
//                               pvecback,
//                               pvectsz,
//                               pba,
//                               ppm,
//                               pnl,
//                               ptsz);



  double epsrel=ptsz->mass_epsrel;
  double epsabs=ptsz->mass_epsabs;

  double m_min;
  double m_max;

  if ( ((int) pvectsz[ptsz->index_md] == ptsz->index_md_cov_Y_N )
    || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_cov_Y_N_next_order )
    || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_cov_N_N )
    || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_cov_N_N_hsv )){
    int index_m_1 = (int) pvectsz[ptsz->index_mass_bin_1];
    m_min = ptsz->cov_Y_N_mass_bin_edges[index_m_1];
    m_max = ptsz->cov_Y_N_mass_bin_edges[index_m_1+1];
  }

  else {
    m_min = ptsz->M1SZ;
    m_max = ptsz->M2SZ;
  }

  struct Parameters_for_integrand_patterson V;
  V.pnl = pnl;
  V.ppm = ppm;
  V.ptsz = ptsz;
  V.pba = pba;
  V.pvectsz = pvectsz;
  V.pvecback = pvecback;
  void * params = &V;


  double r; //store result of mass integral




  if ( ((int) pvectsz[ptsz->index_md] == ptsz->index_md_cov_Y_N_next_order )){

  double r_cov_Y_N_next_order_1; // for cov_Y_N_next_order: first part of redshift integrand
  double r_cov_Y_N_next_order_2; // for cov_Y_N_next_order: second part of redshift integrand

  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate within the mass bin ('N' part)
  r_cov_Y_N_next_order_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                             epsrel, epsabs,
                                                             integrand_patterson_test,
                                                             params,ptsz->patterson_show_neval);

  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  m_min = ptsz->M1SZ;
  m_max = ptsz->M2SZ;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Y' part)
  r_cov_Y_N_next_order_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                             epsrel, epsabs,
                                                             integrand_patterson_test,
                                                             params,ptsz->patterson_show_neval);
  r = r_cov_Y_N_next_order_1*r_cov_Y_N_next_order_2;
                                     }

  else if ( ((int) pvectsz[ptsz->index_md] == ptsz->index_md_cov_N_N_hsv )){

  double r_cov_N_N_hsv_1; // for cov_Y_N_next_order: first part of redshift integrand
  double r_cov_N_N_hsv_2; // for cov_Y_N_next_order: second part of redshift integrand

  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate within the mass bin ('N' part)
  r_cov_N_N_hsv_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                     epsrel, epsabs,
                                                     integrand_patterson_test,
                                                     params,ptsz->patterson_show_neval);

  pvectsz[ptsz->index_part_id_cov_hsv] = 2;

  int index_m_2 = (int) pvectsz[ptsz->index_mass_bin_2];
  m_min = ptsz->cov_Y_N_mass_bin_edges[index_m_2];
  m_max = ptsz->cov_Y_N_mass_bin_edges[index_m_2+1];


  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Y' part)
  r_cov_N_N_hsv_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                     epsrel, epsabs,
                                                     integrand_patterson_test,
                                                     params,ptsz->patterson_show_neval);
  r = r_cov_N_N_hsv_1*r_cov_N_N_hsv_2;
                                     }

  else if ( ((int) pvectsz[ptsz->index_md] == ptsz->index_md_cov_Y_Y_ssc )){

  double r_cov_Y_Y_ssc_1; // for cov_Y_Y_ssc: first part of redshift integrand
  double r_cov_Y_Y_ssc_2; // for cov_Y_Y_ssc: second part of redshift integrand

  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  pvectsz[ptsz->index_multipole] =  pvectsz[ptsz->index_multipole_1];
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_cov_Y_Y_ssc_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                     epsrel, epsabs,
                                                     integrand_patterson_test,
                                                     params,ptsz->patterson_show_neval);

  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  pvectsz[ptsz->index_multipole] =  pvectsz[ptsz->index_multipole_2];
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Y' part)
  r_cov_Y_Y_ssc_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                     epsrel, epsabs,
                                                     integrand_patterson_test,
                                                     params,ptsz->patterson_show_neval);
  r = r_cov_Y_Y_ssc_1*r_cov_Y_Y_ssc_2;
                                     }

  else if ( ((int) pvectsz[ptsz->index_md] == ptsz->index_md_tSZ_lens_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Phi' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);
  r = r_m_1*r_m_2;
                                     }

  else if ( ((int) pvectsz[ptsz->index_md] == ptsz->index_md_tSZ_cib_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);
  r = r_m_1*r_m_2;
                                     }

  else if ( ((int) pvectsz[ptsz->index_md] == ptsz->index_md_gal_cib_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
   double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);
   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
   double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
                                     }


  else if ( ((int) pvectsz[ptsz->index_md] == ptsz->index_md_lens_cib_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Phi' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);
  r = r_m_1*r_m_2;
                                     }

  else if ( ((int) pvectsz[ptsz->index_md] == ptsz->index_md_cib_cib_2h ) && (pvectsz[ptsz->index_frequency_for_cib_profile] != pvectsz[ptsz->index_frequency_prime_for_cib_profile])){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate for frequency nu
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate for frequency nu_prime
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);
  r = r_m_1*r_m_2;
                                     }



  else if ((int) pvectsz[ptsz->index_md] == ptsz->index_md_tSZ_gal_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

 if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
   double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_1\n");
}


  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('galaxy' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
   double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
    }


  else if ((int) pvectsz[ptsz->index_md] == ptsz->index_md_tSZ_lensmag_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

 if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
   double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_1\n");
}


  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('galaxy' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
   double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
    }



  else if ((int) pvectsz[ptsz->index_md] == ptsz->index_md_gal_lensmag_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('gal' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Phi' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  // corrected to match Mat M. consistency treatment in hmvec
  // in hmvec :
  // (integral+b1-consistency1)
  // with integral = r_m_2, b1 = 1, consistency1 = low mass part
  r_m_2 = r_m_2+ptsz->Omega_m_0*ptsz->Rho_crit_0*pow(pvecback[pba->index_bg_ang_distance]*pba->h,-2.)/pvectsz[ptsz->index_lensing_Sigma_crit];


  r = r_m_1*r_m_2;
  }

  else if ((int) pvectsz[ptsz->index_md] == ptsz->index_md_gal_lens_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('gal' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);


   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Phi' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
   double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  // corrected to match Mat M. consistency treatment in hmvec
  // r_m_2 = r_m_2 + ptsz->Omega_m_0*ptsz->Rho_crit_0*pow(pvecback[pba->index_bg_ang_distance]*pba->h,-2.)/pvectsz[ptsz->index_lensing_Sigma_crit];


  r = r_m_1*r_m_2;
  }

  else if ((int) pvectsz[ptsz->index_md] == ptsz->index_md_kSZ_kSZ_gal_2h){
  double r_m_11; // first part of redshift integrand
  double r_m_21; // second part of redshift integrand
  double r_m_12; // first part of redshift integrand
  double r_m_22; // second part of redshift integrand
  double r_m_13; // first part of redshift integrand
  double r_m_23; // second part of redshift integrand



  // r_m_11*r_m_21
  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_11=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);



   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_11 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_21=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_21 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  // r_m_12*r_m_22
  pvectsz[ptsz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_12=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);


   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_12 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[ptsz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_22=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_22 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  // r_m_13*r_m_23
  pvectsz[ptsz->index_part_id_cov_hsv] = 5;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_13=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_13 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[ptsz->index_part_id_cov_hsv] = 6;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_23=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);


   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_23 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  //int index_l_1 = (int) pvectsz[ptsz->index_multipole_1];
  int index_theta_1 = (int) pvectsz[ptsz->index_multipole_1];
  double theta_1 = ptsz->theta_kSZ2_gal_theta_grid[index_theta_1];
  int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
  int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
  double l2 = exp(ptsz->ell_kSZ2_gal_multipole_grid[index_l_2]);
  double l3 = ptsz->ell[index_l_3];
  double ell = l3;
  double ell_prime = l2;
  double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos(theta_1));
  pvectsz[ptsz->index_multipole_for_pk] = l1;
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk1 = pvectsz[ptsz->index_pk_for_halo_bias];


  pvectsz[ptsz->index_multipole_for_pk] = l2;
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk2 = pvectsz[ptsz->index_pk_for_halo_bias];


  pvectsz[ptsz->index_multipole_for_pk] = l3;
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk3 = pvectsz[ptsz->index_pk_for_halo_bias];


  r = pk3*r_m_11*r_m_21  +  pk2*r_m_12*r_m_22  +  pk1*r_m_13*r_m_23;
  }




  else if ((int) pvectsz[ptsz->index_md] == ptsz->index_md_kSZ_kSZ_gal_2h_fft){
  double r_m_11; // first part of redshift integrand
  double r_m_21; // second part of redshift integrand
  double r_m_12; // first part of redshift integrand
  double r_m_22; // second part of redshift integrand
  double r_m_13; // first part of redshift integrand
  double r_m_23; // second part of redshift integrand

  double z = pvectsz[ptsz->index_z];


  int index_l_3 = (int) pvectsz[ptsz->index_multipole];
  double l3 = ptsz->ell[index_l_3];
  pvectsz[ptsz->index_multipole_for_pk] = l3;
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk3 = pvectsz[ptsz->index_pk_for_halo_bias];




  // r_m_11*r_m_21
  // pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_11=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);


   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_11 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  r_m_21 = pk3*get_psi_b1g_at_l_and_z(l3,z,ptsz);

  // printf("%.5e %.5e\n",r_m_11,r_m_21);



//// r_m_12 and r_m_22
//// no mass integral
//// apply convolution theorem


  /// set-up:

double l_min = 1e-2;
double l_max = 2e4; // this is a precision parameter
// tabulate the integrand in the "l" dimension:
const int N = ptsz->N_samp_fftw;
double k[N], Pk1[N],Pk2[N], Pkr[N];
double lnk[N],lnpk[N];
int ik;
double fl;
// double taul;
double l;
// double m = exp(logM);

// printf("z = %.5e l = %.5e\n",z,l3);

for (ik=0; ik<N; ik++){
k[ik] = exp(log(l_min)+ik/(N-1.)*(log(l_max)-log(l_min)));
lnk[ik] = log(l_min)+ik/(N-1.)*(log(l_max)-log(l_min));
l = k[ik];
fl = get_ksz_filter_at_l(l,ptsz);

Pk1[ik] = fl*get_psi_b1gt_at_l1_l2_and_z(l3,l,z,ptsz);

pvectsz[ptsz->index_multipole_for_pk] = l;
evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
double pkl = pvectsz[ptsz->index_pk_for_halo_bias];
Pk2[ik] = fl*get_psi_b1t_at_l_and_z(l,z,ptsz)*pkl;
// if(l>3e3)
  // printf("k = %.5e pk = %.5e\n",l,Pk2[ik]);
}
// printf("k pk done\n");

double rp[N], xi1[N], xi2[N], xi12[N];

// go to Fourier space:
xi2pk(N,k,Pk1,rp,xi1,ptsz);
xi2pk(N,k,Pk2,rp,xi2,ptsz);
for (ik=0; ik<N; ik++){
// convolution:
xi12[ik] = xi1[ik]*xi2[ik];
}
// printf("xi pi done\n");

// move back to position space:
pk2xi(N,rp,xi12,k,Pkr,ptsz);

// evaluate at l3
double f_psi_f_psi = pwl_value_1d(N,lnk,Pkr,log(l3));

r_m_12 = f_psi_f_psi;
r_m_22 = 1.;

r_m_13 = r_m_12;
r_m_23 = 1.;

r = r_m_11*r_m_21 +  r_m_12*r_m_22  +  r_m_13*r_m_23;
// printf("xi pd done r=%.5e\n",r);
  }



  else if ((int) pvectsz[ptsz->index_md] == ptsz->index_md_kSZ_kSZ_gal_3h_fft){
  double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12;

  double z = pvectsz[ptsz->index_z];


  int index_l_3 = (int) pvectsz[ptsz->index_multipole];
  double l3 = ptsz->ell[index_l_3];
  pvectsz[ptsz->index_multipole_for_pk] = l3;
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk3 = pvectsz[ptsz->index_pk_for_halo_bias];


  double psi_bg = get_psi_b1g_at_l_and_z(l3,z,ptsz);
  double psi_b2g = get_psi_b2g_at_l_and_z(l3,z,ptsz);
  // double psi_b2t = get_psi_b2t_at_l_and_z(l3,z,ptsz);

  // printf("%.5e %.5e\n",r_m_11,r_m_21);



//// r_m_12 and r_m_22
//// no mass integral
//// apply convolution theorem


  /// set-up:

double l_min = 1e-8;
double l_max = 1e8; // this is a precision parameter
// tabulate the integrand in the "l" dimension:
const int N = ptsz->N_samp_fftw;
double k[N];

double t1_xi12[N],t1_Pkr[N];
double t2_xi12[N],t2_Pkr[N];
double t3_xi12[N],t3_Pkr[N];
double t4_xi12[N],t4_Pkr[N];
double t5_xi12[N],t5_Pkr[N];
double t6_xi12[N],t6_Pkr[N];
double t7_xi12[N],t7_Pkr[N];
double t8_xi12[N],t8_Pkr[N];
double t9_xi12[N],t9_Pkr[N];
double t10_xi12[N],t10_Pkr[N];
double t11_xi12[N],t11_Pkr[N];
double t12_xi12[N],t12_Pkr[N];


// double  xi1[N], xi2[N], xi12[N];

double lnk[N];
int ik;
double fl;
// double taul;
double l,pkl;

double pk_phi_0[N],pk_phi_m2[N],pk_phi_4[N],pk_phi_2[N];
double pk_tilde_phi_0[N],pk_tilde_phi_m2[N],pk_tilde_phi_2[N];
double pk_tilde_phi_b20[N];


double xi_phi_0[N],xi_phi_m2[N],xi_phi_4[N],xi_phi_2[N];
double xi_tilde_phi_0[N],xi_tilde_phi_m2[N],xi_tilde_phi_2[N];
double xi_tilde_phi_b20[N];





double psi_bt;
double psi_b2t;
// double m = exp(logM);

// printf("z = %.5e l = %.5e\n",z,l3);

for (ik=0; ik<N; ik++){
k[ik] = exp(log(l_min)+ik/(N-1.)*(log(l_max)-log(l_min)));
lnk[ik] = log(l_min)+ik/(N-1.)*(log(l_max)-log(l_min));
l = k[ik];
pvectsz[ptsz->index_multipole_for_pk] = l;
evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
pkl = pvectsz[ptsz->index_pk_for_halo_bias];
fl = get_ksz_filter_at_l(l,ptsz);
psi_bt = get_psi_b1t_at_l_and_z(l,z,ptsz);
psi_b2t = get_psi_b2t_at_l_and_z(l,z,ptsz);

pk_phi_0[ik] = fl*psi_bt;
pk_phi_m2[ik] = pow(l,-2)*fl*psi_bt;
pk_phi_4[ik] = pow(l,4)*fl*psi_bt;
pk_phi_2[ik] = pow(l,2)*fl*psi_bt;

pk_tilde_phi_0[ik] = fl*pkl*psi_bt;
pk_tilde_phi_m2[ik] = pow(l,-2)*fl*pkl*psi_bt;
pk_tilde_phi_2[ik] = pow(l,2)*fl*pkl*psi_bt;
pk_tilde_phi_b20[ik] =  fl*pkl*psi_b2t;





//
// t1_Pk1[ik] = tilde_phi_0;
// t1_Pk2[ik] = tilde_phi_0;
//
//


// if(l>3e3)
  // printf("k = %.5e pk = %.5e\n",l,Pk2[ik]);
}
// printf("k pk done\n");

double rp[N];

// go to Fourier space:
xi2pk(N,k,pk_phi_0,rp,xi_phi_0,ptsz);
xi2pk(N,k,pk_phi_2,rp,xi_phi_2,ptsz);
xi2pk(N,k,pk_phi_m2,rp,xi_phi_m2,ptsz);
xi2pk(N,k,pk_tilde_phi_0,rp,xi_tilde_phi_0,ptsz);
xi2pk(N,k,pk_tilde_phi_b20,rp,xi_tilde_phi_b20,ptsz);
xi2pk(N,k,pk_tilde_phi_2,rp,xi_tilde_phi_2,ptsz);
xi2pk(N,k,pk_tilde_phi_m2,rp,xi_tilde_phi_m2,ptsz);
xi2pk(N,k,pk_phi_4,rp,xi_phi_4,ptsz);

for (ik=0; ik<N; ik++){
// convolution:
t1_xi12[ik] = xi_tilde_phi_0[ik]*xi_tilde_phi_0[ik];

t2_xi12[ik] = xi_tilde_phi_2[ik]*xi_tilde_phi_m2[ik];

t3_xi12[ik] = xi_tilde_phi_0[ik]*xi_tilde_phi_m2[ik];

t4_xi12[ik] = xi_tilde_phi_m2[ik]*xi_tilde_phi_m2[ik];

t5_xi12[ik] = xi_tilde_phi_0[ik]*xi_phi_0[ik];

t6_xi12[ik] = xi_tilde_phi_2[ik]*xi_phi_0[ik];


t7_xi12[ik] = xi_tilde_phi_m2[ik]*xi_phi_0[ik];

t8_xi12[ik] = xi_tilde_phi_0[ik]*xi_phi_2[ik];

t9_xi12[ik] = xi_tilde_phi_m2[ik]*xi_phi_2[ik];

t10_xi12[ik] = xi_tilde_phi_m2[ik]*xi_phi_4[ik];

t11_xi12[ik] = xi_tilde_phi_b20[ik]* xi_tilde_phi_0[ik];

t12_xi12[ik] = xi_tilde_phi_b20[ik]* xi_phi_0[ik];

}
// printf("xi pi done\n");

// move back to position space:
pk2xi(N,rp,t1_xi12,k,t1_Pkr,ptsz);
pk2xi(N,rp,t2_xi12,k,t2_Pkr,ptsz);
pk2xi(N,rp,t3_xi12,k,t3_Pkr,ptsz);
pk2xi(N,rp,t4_xi12,k,t4_Pkr,ptsz);
pk2xi(N,rp,t5_xi12,k,t5_Pkr,ptsz);
pk2xi(N,rp,t6_xi12,k,t6_Pkr,ptsz);
pk2xi(N,rp,t7_xi12,k,t7_Pkr,ptsz);
pk2xi(N,rp,t8_xi12,k,t8_Pkr,ptsz);
pk2xi(N,rp,t9_xi12,k,t9_Pkr,ptsz);
pk2xi(N,rp,t10_xi12,k,t10_Pkr,ptsz);
pk2xi(N,rp,t11_xi12,k,t11_Pkr,ptsz);
pk2xi(N,rp,t12_xi12,k,t12_Pkr,ptsz);



r = 19./7.*psi_bg*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
+9./7.*psi_bg*pwl_value_1d(N,lnk,t2_Pkr,log(l3))
-11./7.*pow(l3,2.)*psi_bg*pwl_value_1d(N,lnk,t3_Pkr,log(l3))
+1./7.*psi_bg*pow(l3,4.)*pwl_value_1d(N,lnk,t4_Pkr,log(l3))
+24./7.*psi_bg*pk3*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
+2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
+2./7.*pow(l3,2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
-4./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
-4./7.*psi_bg*pk3*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
+2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t10_Pkr,log(l3))
// b2 terms:
+psi_b2g*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
+psi_bg*pwl_value_1d(N,lnk,t11_Pkr,log(l3))
+psi_bg*pk3*pwl_value_1d(N,lnk,t12_Pkr,log(l3));

  }

  else if ((int) pvectsz[ptsz->index_md] == ptsz->index_md_kSZ_kSZ_gal_3h){
  double r_m_b1t1;
  double r_m_b1t2;
  double r_m_b1g3;
  double r_m_b2g3;

  // r_m_11*r_m_21
  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_b1t1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }



  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_b1t2 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }



  // r_m_12*r_m_22
  pvectsz[ptsz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_b1g3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }




  pvectsz[ptsz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias_b2];
     r_m_b2g3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  // int index_l_1 = (int) pvectsz[ptsz->index_multipole_1];
  int index_theta_1 = (int) pvectsz[ptsz->index_multipole_1];
  double theta_1 = ptsz->theta_kSZ2_gal_theta_grid[index_theta_1];
  int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
  int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
  double l2 = exp(ptsz->ell_kSZ2_gal_multipole_grid[index_l_2]);
  double l3 = ptsz->ell[index_l_3];
  double ell = l3;
  double ell_prime = l2;
  double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos(theta_1));

  pvectsz[ptsz->index_multipole_for_pk] = l1;//ptsz->ell_kSZ2_gal_multipole_grid[index_l_1];
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk1 = pvectsz[ptsz->index_pk_for_halo_bias];

  // int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
  pvectsz[ptsz->index_multipole_for_pk] = l2;//ptsz->ell_kSZ2_gal_multipole_grid[index_l_2];
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk2 = pvectsz[ptsz->index_pk_for_halo_bias];

  // int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
  pvectsz[ptsz->index_multipole_for_pk] = l3;//ptsz->ell[index_l_3];
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk3 = pvectsz[ptsz->index_pk_for_halo_bias];


  double z = pvectsz[ptsz->index_z];
  double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z);
  double k1 = (l1 + 0.5)/d_A;
  double k2 = (l2 + 0.5)/d_A;
  double k3 = (l3 + 0.5)/d_A;
  double Fk1k2 = bispectrum_f2_kernel(k1,k2,k3);
  double Fk1k3 = bispectrum_f2_kernel(k3,k1,k2);
  double Fk2k3 = bispectrum_f2_kernel(k2,k3,k1);

  double comb_pks = pk1*pk2+pk1*pk3+pk2*pk3;
  double comb_pks_fks = 2.*pk1*pk2*Fk1k2+2.*pk1*pk3*Fk1k3+2.*pk2*pk3*Fk2k3;


  r = r_m_b1t1*r_m_b1t2*r_m_b1g3*comb_pks_fks+r_m_b1t1*r_m_b1t2*r_m_b2g3*comb_pks;

  }


  else if ((int) pvectsz[ptsz->index_md] == ptsz->index_md_bk_at_z_2h){
  double r_m_b1t1;
  double r_m_b1t2;


  // r_m_11*r_m_21
  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_b1t1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }



  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_b1t2 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  r = r_m_b1t1*r_m_b1t2;

  }

  else if ((int) pvectsz[ptsz->index_md] == ptsz->index_md_bk_at_z_3h){
  double r_m_b1t1;
  double r_m_b2t1;


  // r_m_11*r_m_21
  pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r_m_b1t1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias_b2];
     r_m_b2t1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[ptsz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  double r_m_mean=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double nmin_umin = nmin*I0/pvectsz[ptsz->index_hmf];
     r_m_mean += nmin_umin;
     // printf("counter terms done r_m_1\n");
  }
  int index_k = (int) pvectsz[ptsz->index_k_for_pk_hm];
  double k = ptsz->k_for_pk_hm[index_k];
  printf("k = %.8e int = %.8e rhom = %.8e\n",k,r_m_mean,pow((pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0,1));



  double bh;
  double pk;
  double b0;
  double f2;

  // printf("result 3h = %.3e\n",result);
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  pk = pvectsz[ptsz->index_pk_for_halo_bias];
  f2 = bispectrum_f2_kernel(k,k,k);
  b0 = 3.*(2.*pk*pk*f2)*r_m_b1t1*r_m_b1t1*r_m_b1t1;
  bh = 3.*pk*pk*r_m_b1t1*r_m_b1t1*r_m_b2t1;
  r  = (bh+b0);


  }


  else {
  r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                       epsrel, epsabs,
                                       integrand_patterson_test,
                                       params,ptsz->patterson_show_neval);

// halo model consistency:

if ( (int) pvectsz[ptsz->index_md] != ptsz->index_md_cov_N_N
  && (int) pvectsz[ptsz->index_md] != ptsz->index_md_cov_N_N_hsv){
  if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
 // 'pure' 2-halo cases - correlation of same fields.
 if (( (int) pvectsz[ptsz->index_md] == ptsz->index_md_2halo)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_m_y_y_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_pk_at_z_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_pk_gg_at_z_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_lens_lens_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_lensmag_lensmag_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_lens_lensmag_2h)
 || (((int) pvectsz[ptsz->index_md] == ptsz->index_md_cib_cib_2h)  && (pvectsz[ptsz->index_frequency_for_cib_profile] == pvectsz[ptsz->index_frequency_prime_for_cib_profile]) )
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_gal_gal_2h)){

     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r += bmin_umin;
    }
// 1-halo cases
else {
       double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
       double I0 = integrand_patterson_test(log(ptsz->m_min_counter_terms),params);
       double nmin_umin = nmin*I0/pvectsz[ptsz->index_hmf];
       r += nmin_umin;
       }
     }

                                   }
  }

if (( (int) pvectsz[ptsz->index_md] == ptsz->index_md_2halo)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_m_y_y_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_pk_at_z_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_pk_gg_at_z_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_lens_lens_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_lensmag_lensmag_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_lens_lensmag_2h)
 || (((int) pvectsz[ptsz->index_md] == ptsz->index_md_cib_cib_2h)
      && (pvectsz[ptsz->index_frequency_for_cib_profile] == pvectsz[ptsz->index_frequency_prime_for_cib_profile]) )
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_gal_gal_2h)){
 pvectsz[ptsz->index_integral_over_m] = r*r;
 }
else
pvectsz[ptsz->index_integral_over_m] = r;



//}

return pvectsz[ptsz->index_integral_over_m];

}




//This routine reads the tabulated
//Planck noise map

int read_Planck_noise_map(struct tszspectrum * ptsz)
{
  ///read theta file for completeness
  /////////////////////////////start read theta file
  if (ptsz->sz_verbose >= 3){
    printf("Loading theta file\n");
  }

  //read the thetas
  char line[_LINE_LENGTH_MAX_];
  FILE *process;
  int n_data_guess, n_data = 0;
  //double *thetas = NULL,
  double *tmp = NULL;
  double this_lnx;
  int status;

  n_data = 0;
  n_data_guess = 100;
  ptsz->thetas   = (double *)malloc(n_data_guess*sizeof(double));

  char Filepath[_ARGUMENT_LENGTH_MAX_];
  // //printf("%s\n",ptsz->path_to_class);
  // sprintf(Filepath,
  //         // "%s%s%s",
  //         "%s%s",
  //         "cat ",
  //         // ptsz->path_to_class,
  //         "/sz_auxiliary_files/SZ_thetas.txt");
  //
  // process = popen(Filepath, "r");
  if (ptsz->experiment == 0){
  // class_open(process,"sz_auxiliary_files/SZ_thetas.txt", "r",ptsz->error_message);
  class_open(process,ptsz->Planck_thetas_file, "r",ptsz->error_message);

}
  else if (ptsz->experiment == 1){
  // class_open(process,"sz_auxiliary_files/so_3freqs_020621_thetas.txt", "r",ptsz->error_message);
  class_open(process,ptsz->SO_thetas_file, "r",ptsz->error_message);

}

  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

  if (ptsz->sz_verbose >= 3)
    printf("%lf\n", this_lnx);

    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(ptsz->thetas,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 ptsz->error_message,
                 "Error allocating memory to read the thetas.\n");
      ptsz->thetas = tmp;
    };


    /* Store */
    ptsz->thetas[n_data]   = this_lnx;
    n_data++;
  }

  // status = pclose(process);
  status = fclose(process);
    // printf("Loading theta file 2\n");
  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");
    // printf("Loading theta file 3\n");

  ptsz->nthetas = n_data;

  ptsz->theta_bin_min = ptsz->thetas[0];
  int c;
  for (c = 1; c < ptsz->nthetas; c++)
  {
    if (ptsz->thetas[c] < ptsz->theta_bin_min)
    {
      ptsz->theta_bin_min = ptsz->thetas[c];
      //location = c+1;
    }
  }
  ptsz->theta_bin_max = ptsz->thetas[0];
  for (c = 1; c < ptsz->nthetas; c++)
  {
    if (ptsz->thetas[c] > ptsz->theta_bin_max)
    {
      ptsz->theta_bin_max = ptsz->thetas[c];
      //location = c+1;
    }
  }
  // printf("theta_bin_max:=%e\n",ptsz->theta_bin_max);
  // printf("theta_bin_min:=%e\n",ptsz->theta_bin_min);

  ptsz->Nth = ptsz->nthetas;
  class_alloc(ptsz->erfs_2d_to_1d_th_array,ptsz->Nth*sizeof(double *),ptsz->error_message);


  for (c = 0; c < ptsz->Nth; c++){
    ptsz->erfs_2d_to_1d_th_array[c] = log(ptsz->thetas[c]);

  }

  ///////////////////////////end read theta file

  //end read theta file for completeness
  // start read noise map for completeness
  //read skyfracs

  //double *skyfracs = NULL;

  if (ptsz->sz_verbose >= 3){
    printf("theta file loaded with ntheta = %d\n",ptsz->nthetas);
  }
  // printf("Loading theta file 2\n");

  n_data = 0;
  n_data_guess = 100;
  ptsz->skyfracs   = (double *)malloc(n_data_guess*sizeof(double));

  // sprintf(Filepath,
  //         // "%s%s%s",
  //         "%s%s",
  //         "cat ",
  //         // ptsz->path_to_class,
  //         "/sz_auxiliary_files/SZ_skyfracs.txt");
  //
  // process = popen(Filepath, "r");
  if (ptsz->sz_verbose >= 3){
    printf("Loading skyfrac file\n");
  }
  //class_open(process,"sz_auxiliary_files/SZ_skyfracs.txt", "r",ptsz->error_message);

  if (ptsz->experiment == 0){
  // class_open(process,"sz_auxiliary_files/SZ_skyfracs.txt", "r",ptsz->error_message);
  class_open(process,ptsz->Planck_skyfracs_file, "r",ptsz->error_message);

}
  else if (ptsz->experiment == 1){
  // class_open(process,"sz_auxiliary_files/so_3freqs_020621_skyfracs.txt", "r",ptsz->error_message);
  class_open(process,ptsz->SO_skyfracs_file, "r",ptsz->error_message);

}



  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(ptsz->skyfracs,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 ptsz->error_message,
                 "Error allocating memory to read the thetas.\n");
      ptsz->skyfracs = tmp;
    };


    /* Store */
    ptsz->skyfracs[n_data]   = this_lnx;
    n_data++;
  }

  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  ptsz->nskyfracs = n_data;
  if (ptsz->sz_verbose >= 3){
    printf("sky frac file loaded with nskyfracs = %d\n",ptsz->nskyfracs);
  }

  //end read skyfracs

  if (ptsz->sz_verbose >= 3){
    printf("Loading noise map\n");
  }
  ////////////////////////read the ylims
  int index_patches;
  //double ** ylims = NULL;

  class_alloc(ptsz->ylims,
              ptsz->nskyfracs*sizeof(double *),
              ptsz->error_message);



  for (index_patches=0;
       index_patches<ptsz->nskyfracs;
       index_patches++)
  {
    class_alloc(ptsz->ylims[index_patches],
                ptsz->nthetas*sizeof(double),
                ptsz->error_message);
  }



  if (ptsz->experiment == 0){
    // class_open(process,"sz_auxiliary_files/SZ_ylims.txt", "r",ptsz->error_message);
    class_open(process,ptsz->Planck_ylims_file, "r",ptsz->error_message);

  }
  else if (ptsz->experiment == 1){
    // class_open(process,"sz_auxiliary_files/so_3freqs_020621_ylims.txt", "r",ptsz->error_message);
    class_open(process,ptsz->SO_ylims_file, "r",ptsz->error_message);

  }

  //printf("ok\n");
  // printf("noise map loaded 0\n");
  int id_patches=0;
  int index_thetas = 0;

  for (index_thetas =0;index_thetas<ptsz->nthetas;index_thetas++){
    for (id_patches =0;id_patches<ptsz->nskyfracs;id_patches++){
    fgets(line, sizeof(line)-1, process);
    sscanf(line, "%lf", &this_lnx);
    // printf("%.3e id_p = %d id_t = %d n_p = %d n_t = %d\n",this_lnx,index_patches,index_thetas,ptsz->nskyfracs,ptsz->nthetas);
    ptsz->ylims[id_patches][index_thetas]=this_lnx;
    }
  }
  // printf("noise map loaded 1\n");

  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  ///end read the files
  //end reads noise map for completeness
  class_alloc(ptsz->sky_averaged_ylims,
              ptsz->nthetas*sizeof(double),
              ptsz->error_message);

double sum_skyfracs = 0.;
for (index_patches=0;
     index_patches<ptsz->nskyfracs;
     index_patches++)
     sum_skyfracs += ptsz->skyfracs[index_patches];

for (index_thetas = 0; index_thetas<ptsz->nthetas; index_thetas ++){
  ptsz->sky_averaged_ylims[index_thetas] = 0.;
  for (index_patches=0;
       index_patches<ptsz->nskyfracs;
       index_patches++)
  {
    ptsz->sky_averaged_ylims[index_thetas] += ptsz->skyfracs[index_patches]*ptsz->ylims[index_patches][index_thetas]/sum_skyfracs;
  }
}

if (ptsz->sz_verbose >= 3){
  printf("noise map loaded\n");
}
  return  _SUCCESS_;
}

  int read_SO_Qfit(struct tszspectrum * ptsz){
      //read the Q file
      char line[_LINE_LENGTH_MAX_];
      FILE *process;
      int n_data_guess, n_data = 0;
      //double *thetas = NULL,
      double *tmp = NULL;
      double this_lnx,this_lny;
      int status;

      n_data = 0;
      n_data_guess = 100;
      ptsz->SO_Qfit   = (double *)malloc(n_data_guess*sizeof(double));
      ptsz->SO_thetas   = (double *)malloc(n_data_guess*sizeof(double));

      // char Filepath[_ARGUMENT_LENGTH_MAX_];
      // sprintf(Filepath,
      //         "%s%s",
      //         // "%s%s%s",
      //         "cat ",
      //         // ptsz->path_to_class,
      //         "/sz_auxiliary_files/SO_files/SOSim_3freq_small_Qfit_comp_test.txt");
      //
      // process = popen(Filepath, "r");
      class_open(process,"sz_auxiliary_files/SO_files/SOSim_3freq_small_Qfit_comp_test.txt", "r",ptsz->error_message);

      while (fgets(line, sizeof(line)-1, process) != NULL) {
        sscanf(line, "%lf %lf", &this_lnx, &this_lny);

        if((n_data+1) > n_data_guess) {
          n_data_guess *= 2;
          tmp = (double *)realloc(ptsz->SO_Qfit,   n_data_guess*sizeof(double));
          class_test(tmp == NULL,
                     ptsz->error_message,
                     "Error allocating memory to read SO_Qfit.\n");
          ptsz->SO_Qfit = tmp;
          tmp = (double *)realloc(ptsz->SO_thetas,   n_data_guess*sizeof(double));
          class_test(tmp == NULL,
                     ptsz->error_message,
                     "Error allocating memory to read SO_Qfit.\n");
          ptsz->SO_thetas = tmp;
        };


        /* Store */
        ptsz->SO_thetas[n_data]   = this_lnx;
        ptsz->SO_Qfit[n_data]   = this_lny;
        n_data++;
      }

      // status = pclose(process);
      status = fclose(process);
      class_test(status != 0.,
                 ptsz->error_message,
                 "The attempt to launch the external command was unsuccessful. "
                 "Try doing it by hand to check for errors.");

      ptsz->SO_Q_size = n_data;

      ///////////////////////////end read Q file


  return  _SUCCESS_;
  }



int read_SO_noise(struct tszspectrum * ptsz){
        //read the Q file
        char line[_LINE_LENGTH_MAX_];
        FILE *process;
        int n_data_guess, n_data = 0;
        double *tmp = NULL;
        double this_lnx,this_lny;
        int status;

        n_data = 0;
        n_data_guess = 100;
        ptsz->SO_RMS   = (double *)malloc(n_data_guess*sizeof(double));
        ptsz->SO_skyfrac   = (double *)malloc(n_data_guess*sizeof(double));


        class_open(process,"sz_auxiliary_files/SO_files/SOSim_3freq_small_RMSTab_comp_test.txt", "r",ptsz->error_message);

        while (fgets(line, sizeof(line)-1, process) != NULL) {
          sscanf(line, "%lf %lf", &this_lnx, &this_lny);

          if((n_data+1) > n_data_guess) {
            n_data_guess *= 2;
            tmp = (double *)realloc(ptsz->SO_RMS,   n_data_guess*sizeof(double));
            class_test(tmp == NULL,
                       ptsz->error_message,
                       "Error allocating memory to read SO_Qfit.\n");
            ptsz->SO_RMS = tmp;
            tmp = (double *)realloc(ptsz->SO_skyfrac,   n_data_guess*sizeof(double));
            class_test(tmp == NULL,
                       ptsz->error_message,
                       "Error allocating memory to read SO_Qfit.\n");
            ptsz->SO_skyfrac = tmp;
          };


          /* Store */
          ptsz->SO_skyfrac[n_data]   = this_lnx;
          ptsz->SO_RMS[n_data]   = this_lny;
          n_data++;
        }

        // status = pclose(process);
        status = fclose(process);
        class_test(status != 0.,
                   ptsz->error_message,
                   "The attempt to launch the external command was unsuccessful. "
                   "Try doing it by hand to check for errors.");

        ptsz->SO_RMS_size = n_data;

        ///////////////////////////end read Q file


    return  _SUCCESS_;
    }


//Tabulate vrms2 as functions of redshift
 int tabulate_vrms2_from_pk(struct background * pba,
                            struct nonlinear * pnl,
                            struct primordial * ppm,
                            struct tszspectrum * ptsz){

double z_min,z_max;
if (ptsz->need_sigma==0){
  class_alloc(ptsz->array_redshift,sizeof(double *)*ptsz->n_arraySZ,ptsz->error_message);
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,ptsz->z_for_pk_hm);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,ptsz->z_for_pk_hm);
}



double * vrms2_var;
class_alloc(vrms2_var,
            sizeof(double *),
            ptsz->error_message);


class_alloc(ptsz->array_vrms2_at_z,sizeof(double *)*ptsz->n_arraySZ,ptsz->error_message);

int index_z;



    for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
    {
      if (ptsz->need_sigma== 0){
      ptsz->array_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_arraySZ-1.); // log(1+z)
                                    }

            spectra_vrms2(pba,
                          ppm,
                          pnl,
                          ptsz,
                          exp(ptsz->array_redshift[index_z])-1.,
                          vrms2_var
                          );
          ptsz->array_vrms2_at_z[index_z] = log(*vrms2_var);
          // printf("z = %.3e vrms2 = %.3e\n",ptsz->array_redshift[index_z],ptsz->array_vrms2_at_z[index_z]);

       }

free(vrms2_var);

return _SUCCESS_;
    }





struct Parameters_for_integrand_mean_galaxy_number{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
};



double integrand_mean_galaxy_number(double lnM_halo, void *p){

  struct Parameters_for_integrand_mean_galaxy_number *V = ((struct Parameters_for_integrand_mean_galaxy_number *) p);

    double M_halo = exp(lnM_halo);

    double z = V->z;



      double tau;
      int first_index_back = 0;


      class_call(background_tau_of_z(V->pba,z,&tau),
                 V->pba->error_message,
                 V->pba->error_message);

      class_call(background_at_tau(V->pba,
                                   tau,
                                   V->pba->long_info,
                                   V->pba->inter_normal,
                                   &first_index_back,
                                   V->pvecback),
                 V->pba->error_message,
                 V->pba->error_message);




      V->pvectsz[V->ptsz->index_z] = z;
      V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];

      V->pvectsz[V->ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);

      V->pvectsz[V->ptsz->index_has_galaxy] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->ptsz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];

      double M_min;
      double M0;
      double M1_prime;
      double sigma_log10M;
      double nc,ns;

      M_min = V->ptsz->M_min_HOD;
      M0 = V->ptsz->M0_HOD;
      M1_prime = V->ptsz->M1_prime_HOD;
      sigma_log10M = V->ptsz->sigma_log10M_HOD;
      // }
      nc = HOD_mean_number_of_central_galaxies(z,V->pvectsz[V->ptsz->index_mass_for_galaxies],M_min,sigma_log10M,V->ptsz,V->pba);
      ns = HOD_mean_number_of_satellite_galaxies(z,V->pvectsz[V->ptsz->index_mass_for_galaxies],nc,M0,V->ptsz->alpha_s_HOD,M1_prime,V->ptsz,V->pba);

      double result = hmf*(ns+nc);

  return result;

}



int tabulate_mean_galaxy_number_density(struct background * pba,
                                        struct nonlinear * pnl,
                                        struct primordial * ppm,
                                        struct tszspectrum * ptsz){

class_alloc(ptsz->array_mean_galaxy_number_density,sizeof(double *)*ptsz->n_arraySZ,ptsz->error_message);

int index_z;
double r;
double m_min,m_max;

// here we should always integrate over the full mass range,
// since this is a normalization term

if (ptsz->hm_consistency == 0){
  m_min = 1e10; // this has to be the same as the minimal mass at whch the counter terms are tabulated
  m_max = 1e16; // this has to be the same as the maximal mass at whch the counter terms are tabulated
}
else{
m_min = ptsz->M1SZ;
m_max = ptsz->M2SZ;
}

double * pvecback;
double * pvectsz;


 class_alloc(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);


 // printf("tabulating dndlnM quantities0\n");

for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
        {
          double z = exp(ptsz->array_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_mean_galaxy_number V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ptsz = ptsz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;

          void * params = &V;
          double epsrel=ptsz->mass_epsrel;
          double epsabs=ptsz->mass_epsabs;

          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_mean_galaxy_number,
                                               params,ptsz->patterson_show_neval);

        // here we always impose the consistency condition.
        // add counter terms:
         double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
         double I0 = integrand_mean_galaxy_number(log(m_min),params);
         double nmin_umin = nmin*I0/pvectsz[ptsz->index_hmf];
         r += nmin_umin;



          ptsz->array_mean_galaxy_number_density[index_z] = log(r);

       }
 free(pvecback);
 free(pvectsz);

return _SUCCESS_;
    }


struct Parameters_for_integrand_hmf_counter_terms_b1min{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
};

double integrand_hmf_counter_terms_b1min(double lnM_halo, void *p){

  struct Parameters_for_integrand_hmf_counter_terms_b1min *V = ((struct Parameters_for_integrand_hmf_counter_terms_b1min *) p);

    //double x=exp(ln_x);
    double z = V->z;

    double M_halo = exp(lnM_halo);

      double tau;
      int first_index_back = 0;


      class_call(background_tau_of_z(V->pba,z,&tau),
                 V->pba->error_message,
                 V->pba->error_message);

      class_call(background_at_tau(V->pba,
                                   tau,
                                   V->pba->long_info,
                                   V->pba->inter_normal,
                                   &first_index_back,
                                   V->pvecback),
                 V->pba->error_message,
                 V->pba->error_message);




      V->pvectsz[V->ptsz->index_z] = z;
      V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);

      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];

      double z_asked = z;
      double  m_asked = M_halo;

      double rho_crit_at_z = V->ptsz->Rho_crit_0;
      double Omega_cb = (V->pba->Omega0_cdm + V->pba->Omega0_b);
      double rho_cb = rho_crit_at_z*Omega_cb;


      double result = hmf*M_halo/rho_cb;

      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
      double b1 = V->pvectsz[V->ptsz->index_halo_bias];
      result *= b1;




  return result;

}


struct Parameters_for_integrand_psi_b2t{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l;
};



struct Parameters_for_integrand_psi_b2g{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l;
};



struct Parameters_for_integrand_psi_b1g{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l;
};


struct Parameters_for_integrand_psi_b1t{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l;
};


struct Parameters_for_integrand_psi_b1gt{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l1;
  double l2;
};


struct Parameters_for_integrand_dydz{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
};

struct Parameters_for_integrand_dcib0dz{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  int index_nu;
};


double integrand_dydz(double lnM_halo, void *p){

  struct Parameters_for_integrand_dydz *V = ((struct Parameters_for_integrand_dydz *) p);

    double z = V->z;



    double M_halo = exp(lnM_halo);

      double tau;
      int first_index_back = 0;


      class_call(background_tau_of_z(V->pba,z,&tau),
                 V->pba->error_message,
                 V->pba->error_message);

      class_call(background_at_tau(V->pba,
                                   tau,
                                   V->pba->long_info,
                                   V->pba->inter_normal,
                                   &first_index_back,
                                   V->pvecback),
                 V->pba->error_message,
                 V->pba->error_message);




      V->pvectsz[V->ptsz->index_z] = z;
      V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->ptsz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);


      // request appropriate mass conversion
      V->pvectsz[V->ptsz->index_has_electron_pressure] = 1 ;

      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->ptsz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];
      V->pvectsz[V->ptsz->index_md] = V->ptsz->index_md_dydz;
      V->pvectsz[V->ptsz->index_multipole_for_pressure_profile] = 1.;
      evaluate_pressure_profile(V->pvecback,V->pvectsz,V->pba,V->ptsz);


      double result = hmf*V->pvectsz[V->ptsz->index_pressure_profile];

      // multiply by volume element:
      double H_over_c_in_h_over_Mpc = V->pvecback[V->pba->index_bg_H]/V->pba->h;
      result *= V->pvectsz[V->ptsz->index_chi2]/H_over_c_in_h_over_Mpc;

  return result;

}




double integrand_dcib0dz(double lnM_halo, void *p){

  struct Parameters_for_integrand_dcib0dz *V = ((struct Parameters_for_integrand_dcib0dz *) p);

    double z = V->z;
    int index_nu = V->index_nu;


    double M_halo = exp(lnM_halo);

      double tau;
      int first_index_back = 0;


      class_call(background_tau_of_z(V->pba,z,&tau),
                 V->pba->error_message,
                 V->pba->error_message);

      class_call(background_at_tau(V->pba,
                                   tau,
                                   V->pba->long_info,
                                   V->pba->inter_normal,
                                   &first_index_back,
                                   V->pvecback),
                 V->pba->error_message,
                 V->pba->error_message);




      V->pvectsz[V->ptsz->index_z] = z;
      V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->ptsz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);


      V->pvectsz[V->ptsz->index_has_cib] = 1;

      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->ptsz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];
      V->pvectsz[V->ptsz->index_frequency_for_cib_profile] = index_nu;
      V->pvectsz[V->ptsz->index_md] = V->ptsz->index_md_dcib0dz;
      V->pvectsz[V->ptsz->index_multipole_for_cib_profile] = 1.;
      evaluate_cib_profile(V->pvectsz[V->ptsz->index_mass_for_cib],
                           V->pvectsz[V->ptsz->index_radius_for_cib],
                           V->pvectsz[V->ptsz->index_concentration_for_cib],
                           V->pvecback,V->pvectsz,V->pba,V->ptsz);

      double result = hmf*V->pvectsz[V->ptsz->index_cib_profile];

      // multiply by volume element:
      double H_over_c_in_h_over_Mpc = V->pvecback[V->pba->index_bg_H]/V->pba->h;
      result *= V->pvectsz[V->ptsz->index_chi2]/H_over_c_in_h_over_Mpc;

  return result;

}






double integrand_psi_b1g(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b1g *V = ((struct Parameters_for_integrand_psi_b1g *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell = V->l;


    double M_halo = exp(lnM_halo);

      double tau;
      int first_index_back = 0;


      class_call(background_tau_of_z(V->pba,z,&tau),
                 V->pba->error_message,
                 V->pba->error_message);

      class_call(background_at_tau(V->pba,
                                   tau,
                                   V->pba->long_info,
                                   V->pba->inter_normal,
                                   &first_index_back,
                                   V->pvecback),
                 V->pba->error_message,
                 V->pba->error_message);




      V->pvectsz[V->ptsz->index_z] = z;
      V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->ptsz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);

      V->pvectsz[V->ptsz->index_has_galaxy] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->ptsz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];

      V->pvectsz[V->ptsz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->ptsz);
      V->pvectsz[V->ptsz->index_multipole_for_galaxy_profile] = ell;//ptsz->ell[index_l_3];
      double chi = sqrt(V->pvectsz[V->ptsz->index_chi2]);
      double kl = (ell+0.5)/chi;
      evaluate_galaxy_profile_2h(kl,V->pvectsz[V->ptsz->index_mass_for_galaxies],
                                 V->pvectsz[V->ptsz->index_radius_for_galaxies],
                                 V->pvectsz[V->ptsz->index_concentration_for_galaxies],
                                 V->pvecback,V->pvectsz,V->pba,V->ptsz);
      double g = V->pvectsz[V->ptsz->index_galaxy_profile];


      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
      double b1 = V->pvectsz[V->ptsz->index_halo_bias];
      double result = hmf*b1*g;



  return result;

}




double integrand_psi_b2g(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b2g *V = ((struct Parameters_for_integrand_psi_b2g *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell = V->l;


    double M_halo = exp(lnM_halo);

      double tau;
      int first_index_back = 0;


      class_call(background_tau_of_z(V->pba,z,&tau),
                 V->pba->error_message,
                 V->pba->error_message);

      class_call(background_at_tau(V->pba,
                                   tau,
                                   V->pba->long_info,
                                   V->pba->inter_normal,
                                   &first_index_back,
                                   V->pvecback),
                 V->pba->error_message,
                 V->pba->error_message);




      V->pvectsz[V->ptsz->index_z] = z;
      V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->ptsz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);


      V->pvectsz[V->ptsz->index_has_galaxy] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->ptsz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];

      V->pvectsz[V->ptsz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->ptsz);
      V->pvectsz[V->ptsz->index_multipole_for_galaxy_profile] = ell;
      double chi = sqrt(V->pvectsz[V->ptsz->index_chi2]);
      double kl = (ell+0.5)/chi;
      evaluate_galaxy_profile_2h(kl,V->pvectsz[V->ptsz->index_mass_for_galaxies],
                                 V->pvectsz[V->ptsz->index_radius_for_galaxies],
                                 V->pvectsz[V->ptsz->index_concentration_for_galaxies],
                                 V->pvecback,V->pvectsz,V->pba,V->ptsz);
      double g = V->pvectsz[V->ptsz->index_galaxy_profile];

      evaluate_halo_bias_b2(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
      double b2 = V->pvectsz[V->ptsz->index_halo_bias_b2];


      double result = hmf*b2*g;


  return result;

}



double integrand_psi_b2t(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b2t *V = ((struct Parameters_for_integrand_psi_b2t *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell = V->l;


    double M_halo = exp(lnM_halo);

      double tau;
      int first_index_back = 0;


      class_call(background_tau_of_z(V->pba,z,&tau),
                 V->pba->error_message,
                 V->pba->error_message);

      class_call(background_at_tau(V->pba,
                                   tau,
                                   V->pba->long_info,
                                   V->pba->inter_normal,
                                   &first_index_back,
                                   V->pvecback),
                 V->pba->error_message,
                 V->pba->error_message);




      V->pvectsz[V->ptsz->index_z] = z;
      V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->ptsz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);


      V->pvectsz[V->ptsz->index_has_electron_density] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->ptsz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];

      V->pvectsz[V->ptsz->index_multipole_for_tau_profile] = ell;
      evaluate_tau_profile(V->pvecback,V->pvectsz,V->pba,V->ptsz);
      double t = V->pvectsz[V->ptsz->index_tau_profile];

      evaluate_halo_bias_b2(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
      double b2 = V->pvectsz[V->ptsz->index_halo_bias_b2];


      double result = hmf*b2*t;


  return result;

}





double integrand_psi_b1t(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b1t *V = ((struct Parameters_for_integrand_psi_b1t *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell = V->l;

    double M_halo = exp(lnM_halo);

    double tau;
    int first_index_back = 0;


      class_call(background_tau_of_z(V->pba,z,&tau),
                 V->pba->error_message,
                 V->pba->error_message);

      class_call(background_at_tau(V->pba,
                                   tau,
                                   V->pba->long_info,
                                   V->pba->inter_normal,
                                   &first_index_back,
                                   V->pvecback),
                 V->pba->error_message,
                 V->pba->error_message);




      V->pvectsz[V->ptsz->index_z] = z;
      V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->ptsz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);

      V->pvectsz[V->ptsz->index_has_electron_density] = 1;
      // V->pvectsz[V->ptsz->index_has_matter_density] = 1;

      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->ptsz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];


      V->pvectsz[V->ptsz->index_multipole_for_tau_profile] = ell;
      evaluate_tau_profile(V->pvecback,V->pvectsz,V->pba,V->ptsz);

      // double r_delta_matter = V->pvectsz[V->ptsz->index_radius_for_matter_density];
      // double c_delta_matter = V->pvectsz[V->ptsz->index_concentration_for_matter_density];
      // double k = (ell+0.5)/sqrt(V->pvectsz[V->ptsz->index_chi2]);
      // evaluate_matter_density_profile(k,r_delta_matter,c_delta_matter,V->pvecback,V->pvectsz,V->pba,V->ptsz);

      double t = V->pvectsz[V->ptsz->index_tau_profile];
      // double rhom =  V->pvectsz[V->ptsz->index_density_profile];


      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
      double b1 = V->pvectsz[V->ptsz->index_halo_bias];
      double result = hmf*b1*t;
      // double result = hmf*b1*rhom;



  return result;

}

double integrand_psi_b1gt(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b1gt *V = ((struct Parameters_for_integrand_psi_b1gt *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell1 = V->l1;
    double ell2 = V->l2;

    double M_halo = exp(lnM_halo);

      double tau;
      int first_index_back = 0;


      class_call(background_tau_of_z(V->pba,z,&tau),
                 V->pba->error_message,
                 V->pba->error_message);

      class_call(background_at_tau(V->pba,
                                   tau,
                                   V->pba->long_info,
                                   V->pba->inter_normal,
                                   &first_index_back,
                                   V->pvecback),
                 V->pba->error_message,
                 V->pba->error_message);




      V->pvectsz[V->ptsz->index_z] = z;
      V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->ptsz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);

      V->pvectsz[V->ptsz->index_has_galaxy] = 1;
      V->pvectsz[V->ptsz->index_has_electron_density] = 1;

      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->ptsz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];

      V->pvectsz[V->ptsz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->ptsz);
      V->pvectsz[V->ptsz->index_multipole_for_galaxy_profile] = ell1;
      double chi = sqrt(V->pvectsz[V->ptsz->index_chi2]);
      double kl = (ell1+0.5)/chi;
      evaluate_galaxy_profile_2h(kl,V->pvectsz[V->ptsz->index_mass_for_galaxies],
                                 V->pvectsz[V->ptsz->index_radius_for_galaxies],
                                 V->pvectsz[V->ptsz->index_concentration_for_galaxies],
                                 V->pvecback,V->pvectsz,V->pba,V->ptsz);
      double g = V->pvectsz[V->ptsz->index_galaxy_profile];



      V->pvectsz[V->ptsz->index_multipole_for_tau_profile] = ell2;

      evaluate_tau_profile(V->pvecback,V->pvectsz,V->pba,V->ptsz);
      double t = V->pvectsz[V->ptsz->index_tau_profile];




      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
      double b1 = V->pvectsz[V->ptsz->index_halo_bias];
      double result = hmf*b1*g*t;



  return result;

}




int tabulate_psi_b1g(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct tszspectrum * ptsz){

class_alloc(ptsz->array_psi_b1g_redshift,sizeof(double *)*ptsz->n_z_psi_b1g,ptsz->error_message);
class_alloc(ptsz->array_psi_b1g_multipole,sizeof(double *)*ptsz->n_l_psi_b1g,ptsz->error_message);

class_alloc(ptsz->array_psi_b1g_psi,sizeof(double *)*ptsz->n_l_psi_b1g*ptsz->n_z_psi_b1g,ptsz->error_message);


int index_z, index_l;
double r;
double m_min,m_max;


m_min = ptsz->M1SZ; // for the mass integral
m_max = ptsz->M2SZ; // for the mass integral
// m_min = ptsz->m_min_counter_terms;
// m_max = ptsz->m_max_counter_terms;
double z_min = ptsz->z1SZ;
double z_max = ptsz->z2SZ;

// printf("ptsz->n_z_psi_b1g = %d\n",ptsz->n_z_psi_b1g);

for (index_z=0; index_z<ptsz->n_z_psi_b1g; index_z++)
        {

          ptsz->array_psi_b1g_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_psi_b1g-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 3e4;


for (index_l=0; index_l<ptsz->n_l_psi_b1g; index_l++)
        {

          ptsz->array_psi_b1g_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(ptsz->n_l_psi_b1g-1.); // log(l)
        }


double * pvecback;
double * pvectsz;

double tstart, tstop;
int abort;
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,\
pba,ptsz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);


#pragma omp for schedule (dynamic)
for (index_l=0; index_l<ptsz->n_l_psi_b1g; index_l++)
{
#pragma omp flush(abort)

double l = exp(ptsz->array_psi_b1g_multipole[index_l]);

for (index_z=0; index_z<ptsz->n_z_psi_b1g; index_z++)
        {

          int index_l_z = index_l * ptsz->n_z_psi_b1g + index_z;


          double z = exp(ptsz->array_psi_b1g_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b1g V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ptsz = ptsz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;
          V.l = l;

          void * params = &V;
          double epsrel=1e-3;
          double epsabs=1e-100;


          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_psi_b1g,
                                               params,
                                               ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_psi_b1g(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

          ptsz->array_psi_b1g_psi[index_l_z] = log(r);
       }
     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (ptsz->sz_verbose > 0)
         printf("In %s: time spent in parallel region b1g (loop over l's) = %e s for thread %d\n",
                __func__,tstop-tstart,omp_get_thread_num());
     #endif
 free(pvecback);
 free(pvectsz);
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;
    }




int tabulate_psi_b2g(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct tszspectrum * ptsz){

class_alloc(ptsz->array_psi_b2g_redshift,sizeof(double *)*ptsz->n_z_psi_b2g,ptsz->error_message);
class_alloc(ptsz->array_psi_b2g_multipole,sizeof(double *)*ptsz->n_l_psi_b2g,ptsz->error_message);

class_alloc(ptsz->array_psi_b2g_psi,sizeof(double *)*ptsz->n_l_psi_b2g*ptsz->n_z_psi_b2g,ptsz->error_message);


int index_z, index_l;
double r;
double m_min,m_max;


m_min = ptsz->M1SZ; // for the mass integral
m_max = ptsz->M2SZ; // for the mass integral
// m_min = ptsz->m_min_counter_terms;
// m_max = ptsz->m_max_counter_terms;
double z_min = ptsz->z1SZ;
double z_max = ptsz->z2SZ;


for (index_z=0; index_z<ptsz->n_z_psi_b2g; index_z++)
        {

          ptsz->array_psi_b2g_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_psi_b2g-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 3e4;


for (index_l=0; index_l<ptsz->n_l_psi_b2g; index_l++)
        {

          ptsz->array_psi_b2g_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(ptsz->n_l_psi_b2g-1.); // log(l)
        }


double * pvecback;
double * pvectsz;

double tstart, tstop;
int abort;
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,\
pba,ptsz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);


#pragma omp for schedule (dynamic)
for (index_l=0; index_l<ptsz->n_l_psi_b2g; index_l++)
{
#pragma omp flush(abort)

double l = exp(ptsz->array_psi_b2g_multipole[index_l]);

for (index_z=0; index_z<ptsz->n_z_psi_b2g; index_z++)
        {

          int index_l_z = index_l * ptsz->n_z_psi_b2g + index_z;


          double z = exp(ptsz->array_psi_b2g_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b2g V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ptsz = ptsz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;
          V.l = l;

          void * params = &V;
          double epsrel=1e-3;
          double epsabs=1e-100;



          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_psi_b2g,
                                               params,
                                               ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_psi_b2g(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias_b2];
     r += bmin_umin;
  }

          ptsz->array_psi_b2g_psi[index_l_z] = log(r);
       }
     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (ptsz->sz_verbose > 0)
         printf("In %s: time spent in parallel region b2g (loop over l's) = %e s for thread %d\n",
                __func__,tstop-tstart,omp_get_thread_num());
     #endif
 free(pvecback);
 free(pvectsz);
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;
    }




int tabulate_psi_b2t(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct tszspectrum * ptsz){

class_alloc(ptsz->array_psi_b2t_redshift,sizeof(double *)*ptsz->n_z_psi_b2t,ptsz->error_message);
class_alloc(ptsz->array_psi_b2t_multipole,sizeof(double *)*ptsz->n_l_psi_b2t,ptsz->error_message);

class_alloc(ptsz->array_psi_b2t_psi,sizeof(double *)*ptsz->n_l_psi_b2t*ptsz->n_z_psi_b2t,ptsz->error_message);


int index_z, index_l;
double r;
double m_min,m_max;


m_min = ptsz->M1SZ; // for the mass integral
m_max = ptsz->M2SZ; // for the mass integral
// m_min = ptsz->m_min_counter_terms;
// m_max = ptsz->m_max_counter_terms;
double z_min = ptsz->z1SZ;
double z_max = ptsz->z2SZ;


for (index_z=0; index_z<ptsz->n_z_psi_b2t; index_z++)
        {

          ptsz->array_psi_b2t_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_psi_b2t-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 3e4;


for (index_l=0; index_l<ptsz->n_l_psi_b2t; index_l++)
        {

          ptsz->array_psi_b2t_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(ptsz->n_l_psi_b2t-1.); // log(l)
        }


double * pvecback;
double * pvectsz;

double tstart, tstop;
int abort;
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,\
pba,ptsz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);


#pragma omp for schedule (dynamic)
for (index_l=0; index_l<ptsz->n_l_psi_b2t; index_l++)
{
#pragma omp flush(abort)

double l = exp(ptsz->array_psi_b2t_multipole[index_l]);

for (index_z=0; index_z<ptsz->n_z_psi_b2t; index_z++)
        {

          int index_l_z = index_l * ptsz->n_z_psi_b2t + index_z;


          double z = exp(ptsz->array_psi_b2t_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b2t V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ptsz = ptsz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;
          V.l = l;

          void * params = &V;
          double epsrel=1e-3;
          double epsabs=1e-100;



          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_psi_b2t,
                                               params,
                                               ptsz->patterson_show_neval);
          // printf("%.8e %.8e %.8e\n",z,l,r);

double ct = 0.;
   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_psi_b2t(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias_b2];
     double ct = bmin_umin;
     r += bmin_umin;

     if (r<-1.){
       printf("%.8e %.8e int = %.8e ct = %.8e %.8e %.8e\n",z,l,r-ct,ct,nmin,r);
       printf("sort out the tabulation of b2g.\n");
     }
  }


          ptsz->array_psi_b2t_psi[index_l_z] = log(1.+r);
       }
     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (ptsz->sz_verbose > 0)
         printf("In %s: time spent in parallel region b2g (loop over l's) = %e s for thread %d\n",
                __func__,tstop-tstart,omp_get_thread_num());
     #endif
 free(pvecback);
 free(pvectsz);
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;
    }



int tabulate_psi_b1t(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct tszspectrum * ptsz){

class_alloc(ptsz->array_psi_b1t_redshift,sizeof(double *)*ptsz->n_z_psi_b1t,ptsz->error_message);
class_alloc(ptsz->array_psi_b1t_multipole,sizeof(double *)*ptsz->n_l_psi_b1t,ptsz->error_message);

class_alloc(ptsz->array_psi_b1t_psi,sizeof(double *)*ptsz->n_l_psi_b1t*ptsz->n_z_psi_b1t,ptsz->error_message);


int index_z, index_l;
double r;
double m_min,m_max;


m_min = ptsz->M1SZ; // for the mass integral
m_max = ptsz->M2SZ; // for the mass integral
// m_min = ptsz->m_min_counter_terms;
// m_max = ptsz->m_max_counter_terms;
double z_min = ptsz->z1SZ;
double z_max = ptsz->z2SZ;

// printf("ptsz->n_z_psi_b1g = %d\n",ptsz->n_z_psi_b1g);

for (index_z=0; index_z<ptsz->n_z_psi_b1t; index_z++)
        {

          ptsz->array_psi_b1t_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_psi_b1t-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 3e4;


for (index_l=0; index_l<ptsz->n_l_psi_b1t; index_l++)
        {

          ptsz->array_psi_b1t_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(ptsz->n_l_psi_b1t-1.); // log(l)
        }


double * pvecback;
double * pvectsz;

double tstart, tstop;
int abort;
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,\
pba,ptsz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);


#pragma omp for schedule (dynamic)
for (index_l=0; index_l<ptsz->n_l_psi_b1t; index_l++)
{
#pragma omp flush(abort)

double l = exp(ptsz->array_psi_b1t_multipole[index_l]);

for (index_z=0; index_z<ptsz->n_z_psi_b1t; index_z++)
        {

          int index_l_z = index_l * ptsz->n_z_psi_b1t + index_z;


          double z = exp(ptsz->array_psi_b1t_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b1t V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ptsz = ptsz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;
          V.l = l;

          void * params = &V;
          double epsrel=1e-3;
          double epsabs=1e-100;

          pvectsz[ptsz->index_has_electron_density] = 1; //

          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_psi_b1t,
                                               params,
                                               ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_psi_b1t(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r += bmin_umin;
     double ct_over_int  = bmin_umin/(r-bmin_umin);
     if (ct_over_int>0.1)
     printf("z = %.8e l = %.8e int = %.8e ct = %.8e ct/int = %.8e\n",
     z,l,r-bmin_umin,bmin_umin,ct_over_int);

  }

          ptsz->array_psi_b1t_psi[index_l_z] = log(r);

       }
     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (ptsz->sz_verbose > 0)
         printf("In %s: time spent in parallel region b1t (loop over l's) = %e s for thread %d\n",
                __func__,tstop-tstart,omp_get_thread_num());
     #endif
 free(pvecback);
 free(pvectsz);
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;
    }



int tabulate_psi_b1gt(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct tszspectrum * ptsz){

class_alloc(ptsz->array_psi_b1gt_redshift,sizeof(double *)*ptsz->n_z_psi_b1gt,ptsz->error_message);
class_alloc(ptsz->array_psi_b1gt_multipole,sizeof(double *)*ptsz->n_l_psi_b1gt,ptsz->error_message);

class_alloc(ptsz->array_psi_b1gt_psi,sizeof(double *)*ptsz->n_z_psi_b1gt,ptsz->error_message);
int index_z, index_l1,index_l2;
for (index_z=0;index_z<ptsz->n_z_psi_b1gt;index_z++){
  class_alloc(ptsz->array_psi_b1gt_psi[index_z],sizeof(double *)*ptsz->n_l_psi_b1gt*ptsz->n_l_psi_b1gt,ptsz->error_message);

}



double r;
double m_min,m_max;


m_min = ptsz->M1SZ; // for the mass integral
m_max = ptsz->M2SZ; // for the mass integral
// m_min = ptsz->m_min_counter_terms;
// m_max = ptsz->m_max_counter_terms;
double z_min = ptsz->z1SZ;
double z_max = ptsz->z2SZ;

// printf("ptsz->n_z_psi_b1g = %d\n",ptsz->n_z_psi_b1g);

for (index_z=0; index_z<ptsz->n_z_psi_b1gt; index_z++)
        {

          ptsz->array_psi_b1gt_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_psi_b1gt-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 3e4;

int index_l;
for (index_l=0; index_l<ptsz->n_l_psi_b1gt; index_l++)
        {

          ptsz->array_psi_b1gt_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(ptsz->n_l_psi_b1gt-1.); // log(l)
        }


double * pvecback;
double * pvectsz;

double tstart, tstop;
int abort;
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,\
pba,ptsz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l1,index_l2,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);


// printf("doing well\n");

#pragma omp for collapse(3)
//#pragma omp for schedule (dynamic)
for (index_z=0; index_z<ptsz->n_z_psi_b1gt; index_z++)
{
//#pragma omp flush(abort)
  for (index_l1=0; index_l1<ptsz->n_l_psi_b1gt; index_l1++)
  {
    for (index_l2=0; index_l2<ptsz->n_l_psi_b1gt; index_l2++)
    {
// #pragma omp flush(abort)


          double l1 = exp(ptsz->array_psi_b1gt_multipole[index_l1]);
          double l2 = exp(ptsz->array_psi_b1gt_multipole[index_l2]);

          int index_l1_l2 = index_l2 * ptsz->n_l_psi_b1gt + index_l1;


          double z = exp(ptsz->array_psi_b1gt_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b1gt V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ptsz = ptsz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;
          V.l1 = l1;
          V.l2 = l2;

          void * params = &V;
          double epsrel=1e-3;
          double epsabs=1e-100;


          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_psi_b1gt,
                                               params,
                                               ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double I0 = integrand_psi_b1gt(log(ptsz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

          ptsz->array_psi_b1gt_psi[index_z][index_l1_l2] = log(r);

          // printf("ptsz->array_psi_b1t_psi[%d] = %.5e\n",index_l_z,r);

       }
     }
   }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (ptsz->sz_verbose > 0)
         printf("In %s: time spent in parallel region b1gt (loop over llz's) = %e s for thread %d\n",
                __func__,tstop-tstart,omp_get_thread_num());
     #endif
 free(pvecback);
 free(pvectsz);
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;
    }




int tabulate_dydz(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct tszspectrum * ptsz){

class_alloc(ptsz->array_dydz_redshift,sizeof(double *)*ptsz->n_z_dydz,ptsz->error_message);
class_alloc(ptsz->array_dydz_at_z,sizeof(double *)*ptsz->n_z_dydz,ptsz->error_message);




double r;
double m_min,m_max;


m_min = ptsz->M1SZ; // for the mass integral
m_max = ptsz->M2SZ; // for the mass integral

double z_min = ptsz->z1SZ;
double z_max = ptsz->z2SZ;
int index_z;


for (index_z=0; index_z<ptsz->n_z_dydz; index_z++)
        {

          ptsz->array_dydz_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_dydz-1.); // log(1+z)
        }


double * pvecback;
double * pvectsz;

double tstart, tstop;
int abort;
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,\
pba,ptsz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);

#pragma omp for collapse(1)
for (index_z=0; index_z<ptsz->n_z_dydz; index_z++)
{

          double z = exp(ptsz->array_dydz_redshift[index_z])-1.;;

          // at each z, perform the mass integral
          struct Parameters_for_integrand_dydz V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ptsz = ptsz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;


          void * params = &V;
          double epsrel=1e-3;
          double epsabs=1e-100;



          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_dydz,
                                               params,
                                               ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
       double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
       double I0 = integrand_dydz(log(ptsz->m_min_counter_terms),params);
       double nmin_umin = nmin*I0/pvectsz[ptsz->index_hmf];
       r += nmin_umin;
  }

          ptsz->array_dydz_at_z[index_z] = log(r/pow(ptsz->Tcmb_gNU,1)/1.e6);

     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (ptsz->sz_verbose > 0)
         printf("In %s: time spent in parallel region (loop over z nu's) = %e s for thread %d\n",
                __func__,tstop-tstart,omp_get_thread_num());
     #endif
 free(pvecback);
 free(pvectsz);
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;
    }




int tabulate_dcib0dz(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct tszspectrum * ptsz){

class_alloc(ptsz->array_dcib0dz_redshift,sizeof(double *)*ptsz->n_z_dcib0dz,ptsz->error_message);
class_alloc(ptsz->array_dcib0dz_nu,sizeof(double *)*ptsz->n_nu_dcib0dz,ptsz->error_message);

class_alloc(ptsz->array_dcib0dz_at_z_nu,sizeof(double *)*ptsz->n_z_dcib0dz*ptsz->n_nu_dcib0dz,ptsz->error_message);




double r;
double m_min,m_max;


m_min = ptsz->M1SZ; // for the mass integral
m_max = ptsz->M2SZ; // for the mass integral

double z_min = ptsz->z1SZ;
double z_max = ptsz->z2SZ;
int index_z;


for (index_z=0; index_z<ptsz->n_z_dcib0dz; index_z++)
        {

          ptsz->array_dcib0dz_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_dcib0dz-1.); // log(1+z)
        }

double nu_min = ptsz->freq_min;
double nu_max = ptsz->freq_max;

int index_nu;
for (index_nu=0; index_nu<ptsz->n_nu_dcib0dz; index_nu++)
        {

          ptsz->array_dcib0dz_nu[index_nu] =
                                      log(nu_min)
                                      +index_nu*(log(nu_max)-log(nu_min))
                                      /(ptsz->n_nu_dcib0dz-1.); // log(nu)
        }


double * pvecback;
double * pvectsz;

double tstart, tstop;
int abort;
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif
// number_of_threads= 1;
#pragma omp parallel \
shared(abort,\
pba,ptsz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_nu,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);

#pragma omp for collapse(2)
for (index_z=0; index_z<ptsz->n_z_dcib0dz; index_z++)
{
  for (index_nu=0; index_nu<ptsz->n_nu_dcib0dz; index_nu++)
  {
          double z = exp(ptsz->array_dcib0dz_redshift[index_z])-1.;;
          double nu = exp(ptsz->array_dcib0dz_nu[index_nu]);

          int index_z_nu = index_nu * ptsz->n_z_dcib0dz + index_z;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_dcib0dz V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ptsz = ptsz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;
          V.index_nu = index_nu;

          void * params = &V;
          double epsrel=1e-3;
          double epsabs=1e-100;



          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_dcib0dz,
                                               params,
                                               ptsz->patterson_show_neval);

   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
       double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
       double I0 = integrand_dcib0dz(log(ptsz->m_min_counter_terms),params);
       double nmin_umin = nmin*I0/pvectsz[ptsz->index_hmf];
       r += nmin_umin;
     // printf("counter terms done r_m_1\n");
  }

          ptsz->array_dcib0dz_at_z_nu[index_z_nu] = log(r);

       }
     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (ptsz->sz_verbose > 0)
         printf("In %s: time spent in parallel region (loop over z nu's) = %e s for thread %d\n",
                __func__,tstop-tstart,omp_get_thread_num());
     #endif
 free(pvecback);
 free(pvectsz);
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;
    }





int tabulate_hmf_counter_terms_b1min(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct tszspectrum * ptsz){

class_alloc(ptsz->array_hmf_counter_terms_b1min,sizeof(double *)*ptsz->n_z_hmf_counter_terms,ptsz->error_message);

int index_z;
double r;
double m_min,m_max;


m_min = ptsz->M1SZ;
m_max = ptsz->M2SZ;

double z_min = ptsz->z1SZ;
double z_max = ptsz->z2SZ;

double * pvecback;
double * pvectsz;




 class_alloc(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);




for (index_z=0; index_z<ptsz->n_z_hmf_counter_terms; index_z++)
        {
          double z = exp(ptsz->array_redshift_hmf_counter_terms[index_z])-1.;

          // at each z, perform the mass integral
          struct Parameters_for_integrand_hmf_counter_terms_b1min V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ptsz = ptsz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;

          void * params = &V;
          double epsrel=ptsz->mass_epsrel;
          double epsabs=ptsz->mass_epsabs;

          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_hmf_counter_terms_b1min,
                                               params,ptsz->patterson_show_neval);

          // here is (1-r) = n_min*m_min/rho_cb at z
          double rho_crit_at_z = ptsz->Rho_crit_0;//pvectsz[ptsz->index_Rho_crit];
          double Omega_cb = (pba->Omega0_cdm + pba->Omega0_b);//*pow(1.+z,3.);
          double rho_cb = rho_crit_at_z*Omega_cb;
          double n_min =  get_hmf_counter_term_nmin_at_z(z,ptsz);
          double b1_min = (1.-r)*rho_cb/m_min/n_min;
          ptsz->array_hmf_counter_terms_b1min[index_z] = b1_min; // is actually b1_min*m_min

       }
 free(pvecback);
 free(pvectsz);

return _SUCCESS_;
    }

///// b2 min


struct Parameters_for_integrand_hmf_counter_terms_b2min{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
};

double integrand_hmf_counter_terms_b2min(double lnM_halo, void *p){

  struct Parameters_for_integrand_hmf_counter_terms_b2min *V = ((struct Parameters_for_integrand_hmf_counter_terms_b2min *) p);

    //double x=exp(ln_x);
    double z = V->z;

    double M_halo = exp(lnM_halo);

      double tau;
      int first_index_back = 0;


      class_call(background_tau_of_z(V->pba,z,&tau),
                 V->pba->error_message,
                 V->pba->error_message);

      class_call(background_at_tau(V->pba,
                                   tau,
                                   V->pba->long_info,
                                   V->pba->inter_normal,
                                   &first_index_back,
                                   V->pvecback),
                 V->pba->error_message,
                 V->pba->error_message);




      V->pvectsz[V->ptsz->index_z] = z;
      V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);

      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];

      double z_asked = z;
      double  m_asked = M_halo;
      double rho_crit_at_z = V->ptsz->Rho_crit_0;
      double Omega_cb = (V->pba->Omega0_cdm + V->pba->Omega0_b);
      double rho_cb = rho_crit_at_z*Omega_cb;


      double result = hmf*M_halo/rho_cb;

      evaluate_halo_bias_b2(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
      double b2 = V->pvectsz[V->ptsz->index_halo_bias_b2];
      result *= b2;

  return result;

}



int tabulate_hmf_counter_terms_b2min(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct tszspectrum * ptsz){

class_alloc(ptsz->array_hmf_counter_terms_b2min,sizeof(double *)*ptsz->n_z_hmf_counter_terms,ptsz->error_message);

int index_z;
double r;
double m_min,m_max;
// m_min = ptsz->m_min_counter_terms;
// m_max = ptsz->m_max_counter_terms;
m_min = ptsz->M1SZ;//r8_min(ptsz->M1SZ,ptsz->m_min_counter_terms);
m_max = ptsz->M2SZ;//r8_max(ptsz->M2SZ,ptsz->m_max_counter_terms);
// double z_min = ptsz->z1SZ;
// double z_max = ptsz->z2SZ;
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,ptsz->z_for_pk_hm);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,ptsz->z_for_pk_hm);

double * pvecback;
double * pvectsz;




 class_alloc(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);




for (index_z=0; index_z<ptsz->n_z_hmf_counter_terms; index_z++)
        {
          double z = exp(ptsz->array_redshift_hmf_counter_terms[index_z])-1.;

          // at each z, perform the mass integral
          struct Parameters_for_integrand_hmf_counter_terms_b2min V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ptsz = ptsz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;

          void * params = &V;
          double epsrel=ptsz->mass_epsrel;
          double epsabs=ptsz->mass_epsabs;

          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_hmf_counter_terms_b2min,
                                               params,ptsz->patterson_show_neval);

          // here is (1-r) = n_min*m_min/rho_cb at z
          double rho_crit_at_z = ptsz->Rho_crit_0;//pvectsz[ptsz->index_Rho_crit];
          double Omega_cb = (pba->Omega0_cdm + pba->Omega0_b);//*pow(1.+z,3.);
          double rho_cb = rho_crit_at_z*Omega_cb;
          double n_min =  get_hmf_counter_term_nmin_at_z(z,ptsz);
          double b2_min = -r*rho_cb/m_min/n_min;
          ptsz->array_hmf_counter_terms_b2min[index_z] = b2_min; // is actually b2_min*m_min

       }
 free(pvecback);
 free(pvectsz);

return _SUCCESS_;
    }




//// b2 min


struct Parameters_for_integrand_hmf_counter_terms_nmin{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct tszspectrum * ptsz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
};

double integrand_hmf_counter_terms_nmin(double lnM_halo, void *p){

  struct Parameters_for_integrand_hmf_counter_terms_nmin *V = ((struct Parameters_for_integrand_hmf_counter_terms_nmin *) p);


    double z = V->z;

    double M_halo = exp(lnM_halo);

      double tau;
      int first_index_back = 0;


      class_call(background_tau_of_z(V->pba,z,&tau),
                 V->pba->error_message,
                 V->pba->error_message);

      class_call(background_at_tau(V->pba,
                                   tau,
                                   V->pba->long_info,
                                   V->pba->inter_normal,
                                   &first_index_back,
                                   V->pvecback),
                 V->pba->error_message,
                 V->pba->error_message);




      V->pvectsz[V->ptsz->index_z] = z;
      V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);


      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];

      double z_asked = z;
      double  m_asked = M_halo;

      double rho_crit_at_z = V->ptsz->Rho_crit_0;
      double Omega_cb = (V->pba->Omega0_cdm + V->pba->Omega0_b);
      double rho_cb = rho_crit_at_z*Omega_cb;


      double result = hmf*M_halo/rho_cb;

  return result;

}


int tabulate_hmf_counter_terms_nmin(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct tszspectrum * ptsz){

class_alloc(ptsz->array_hmf_counter_terms_nmin,sizeof(double *)*ptsz->n_z_hmf_counter_terms,ptsz->error_message);
class_alloc(ptsz->array_redshift_hmf_counter_terms,sizeof(double *)*ptsz->n_z_hmf_counter_terms,ptsz->error_message);

int index_z;
double r;
double m_min,m_max;

if (ptsz->hm_consistency == 0){
  m_min = 1e10;
  m_max = 1e16;
}
else{
m_min = ptsz->M1SZ;
m_max = ptsz->M2SZ;
}


double z_min = ptsz->z1SZ;
double z_max = ptsz->z2SZ;

// printf("z_min = %.8e\n",z_min);

double * pvecback;
double * pvectsz;




 class_alloc(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);




for (index_z=0; index_z<ptsz->n_z_hmf_counter_terms; index_z++)
        {

          ptsz->array_redshift_hmf_counter_terms[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_hmf_counter_terms-1.); // log(1+z)

          double z = exp(ptsz->array_redshift_hmf_counter_terms[index_z])-1.;

          // at each z, perform the mass integral
          struct Parameters_for_integrand_hmf_counter_terms_nmin V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ptsz = ptsz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;

          void * params = &V;
          double epsrel=ptsz->mass_epsrel;
          double epsabs=ptsz->mass_epsabs;

          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_hmf_counter_terms_nmin,
                                               params,ptsz->patterson_show_neval);

          // here is (1-r) = n_min*m_min/rho_cb at z
          double rho_crit_at_z = ptsz->Rho_crit_0;
          double Omega_cb = (pba->Omega0_cdm + pba->Omega0_b);
          double rho_cb = rho_crit_at_z*Omega_cb;
          double n_min = (1.-r)*rho_cb/m_min;
          ptsz->array_hmf_counter_terms_nmin[index_z] = n_min; // is actually n_min*m_min

       }
 free(pvecback);
 free(pvectsz);

return _SUCCESS_;
    }




struct Parameters_for_integrand_patterson_L_sat{
  double nu;
  double z;
  double M_host;
  struct tszspectrum * ptsz;
};


double integrand_patterson_L_sat(double lnM_sub, void *p){
  struct Parameters_for_integrand_patterson_L_sat *V = ((struct Parameters_for_integrand_patterson_L_sat *) p);

  double M_sub = exp(lnM_sub);
  double nu = V->nu;
  double z = V->z;
  double M_host = V->M_host;

  double L_gal_at_nu = evaluate_galaxy_luminosity(z, M_sub, nu, V->ptsz);
  double dNdlnMs = subhalo_hmf_dndlnMs(M_host,M_sub);
  double result = L_gal_at_nu*dNdlnMs;


  return result;
}


int tabulate_L_sat_at_z_m_nu(struct background * pba,
                             struct tszspectrum * ptsz){
if ( ptsz->has_cib_monopole + ptsz->has_dcib0dz == _FALSE_)
return 0;

class_alloc(ptsz->array_L_sat_at_M_z_nu,sizeof(double *)*ptsz->n_nu_L_sat,ptsz->error_message);
int index_nu, index_M,index_z;
class_alloc(ptsz->array_nu_L_sat,sizeof(double *)*ptsz->n_nu_L_sat,ptsz->error_message);
double nu_min  = ptsz->freq_min;
double nu_max  = ptsz->freq_max;

for (index_nu=0;index_nu<ptsz->n_nu_L_sat;index_nu++){
  class_alloc(ptsz->array_L_sat_at_M_z_nu[index_nu],sizeof(double *)*ptsz->n_z_L_sat*ptsz->n_m_L_sat,ptsz->error_message);
      ptsz->array_nu_L_sat[index_nu] =
                                      log(nu_min)
                                      +index_nu*(log(nu_max)-log(nu_min))
                                      /(ptsz->n_nu_L_sat-1.); // log(1+z)

}

  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_L_sat);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_L_sat);
  double logM_min = r8_min(log(ptsz->M1SZ/pba->h),log(ptsz->M1SZ_L_sat)); //in Msun
  double logM_max = r8_max(log(ptsz->M2SZ/pba->h),log(ptsz->M2SZ_L_sat)); //in Msun

if (
      ptsz->has_tSZ_cib_1h
    + ptsz->has_tSZ_cib_2h
    + ptsz->has_cib_cib_1h
    // + ptsz->has_cib_monopole
    + ptsz->has_cib_cib_2h
    + ptsz->has_lens_cib_1h
    + ptsz->has_lens_cib_2h
    + ptsz->has_gal_cib_1h
    + ptsz->has_gal_cib_2h
    == _FALSE_
  ){



  class_alloc(ptsz->array_z_L_sat,sizeof(double *)*ptsz->n_z_L_sat,ptsz->error_message);
  class_alloc(ptsz->array_m_L_sat,sizeof(double *)*ptsz->n_m_L_sat,ptsz->error_message);


for (index_z=0; index_z<ptsz->n_z_L_sat; index_z++)
{
      ptsz->array_z_L_sat[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_L_sat-1.); // log(1+z)

}

for (index_M=0; index_M<ptsz->n_m_L_sat; index_M++)
{
      ptsz->array_m_L_sat[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(ptsz->n_m_L_sat-1.); //log(R)


}

}




double r;

double * pvecback;
double * pvectsz;

double tstart, tstop;
int abort;
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif
// number_of_threads = 1;
#pragma omp parallel \
shared(abort,\
pba,ptsz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_z,index_nu,index_M,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);

#pragma omp for collapse(3)
for (index_nu=0; index_nu<ptsz->n_nu_L_sat; index_nu++)
{
  for (index_M=0; index_M<ptsz->n_m_L_sat; index_M++)
  {
    for (index_z=0; index_z<ptsz->n_z_L_sat; index_z++)
    {

          int index_M_z = index_z * ptsz->n_m_L_sat + index_M;


      ptsz->array_z_L_sat[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_L_sat-1.); // log(1+z)

      ptsz->array_m_L_sat[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(ptsz->n_m_L_sat-1.); //log(R)



      double z =   exp(ptsz->array_z_L_sat[index_z])-1.;
      double logM =   ptsz->array_m_L_sat[index_M];


      double lnMs_min = log(1e10);
      double lnMs_max = logM;//log(1e11);

      double epsrel = ptsz->epsrel_L_sat;
      double epsabs = ptsz->epsabs_L_sat;

      struct Parameters_for_integrand_patterson_L_sat V;
      V.nu = exp(ptsz->array_nu_L_sat[index_nu]);
      V.z = z;
      V.ptsz = ptsz;
      V.M_host = exp(logM);


      void * params = &V;
      params = &V;


          r=Integrate_using_Patterson_adaptive(lnMs_min, lnMs_max,
                                               epsrel, epsabs,
                                               integrand_patterson_L_sat,
                                               params,ptsz->patterson_show_neval);

          if (r==0.){
            r = 1e-100;
          }
          ptsz->array_L_sat_at_M_z_nu[index_nu][index_M_z] = log(1.+r);
          if (isnan(ptsz->array_L_sat_at_M_z_nu[index_nu][index_M_z])){
            printf("nan in interp L_sat table\n");
            exit(0);
          }
          if (isinf(ptsz->array_L_sat_at_M_z_nu[index_nu][index_M_z])){
            printf("inf in interp L_sat table\n");
            printf("r = %.3e\n",r);
            exit(0);
          }


       }
     }
   }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (ptsz->sz_verbose > 0)
         printf("In %s: time spent in parallel region Lsat (loop over mnuz's) = %e s for thread %d\n",
                __func__,tstop-tstart,omp_get_thread_num());
     #endif
 free(pvecback);
 free(pvectsz);
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;





                             }

// tabulate L_nu^sat as a function of M (M_host) and z at frequency nu

int tabulate_L_sat_at_nu_and_nu_prime(struct background * pba,
                                      struct tszspectrum * ptsz){

if (
      ptsz->has_tSZ_cib_1h
    + ptsz->has_tSZ_cib_2h
    + ptsz->has_cib_cib_1h
    // + ptsz->has_cib_monopole
    + ptsz->has_cib_cib_2h
    + ptsz->has_lens_cib_1h
    + ptsz->has_lens_cib_2h
    + ptsz->has_gal_cib_1h
    + ptsz->has_gal_cib_2h
    == _FALSE_
    )
return 0;

  //Array of z
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_L_sat);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_L_sat);
  int index_z;

  double tstart, tstop;
  int index_l;

  // double * pvecback;
  // double * pvectsz;
  int abort;

  //Array of M in Msun
  double logM_min = r8_min(log(ptsz->M1SZ/pba->h),log(ptsz->M1SZ_L_sat)); //in Msun
  double logM_max = r8_max(log(ptsz->M2SZ/pba->h),log(ptsz->M2SZ_L_sat)); //in Msun
  int index_M;

  int index_z_M = 0;

  double ** array_L_sat_at_z_and_M_at_nu;
  // double ** array_L_sat_at_z_and_M_at_nu_prime;

  class_alloc(ptsz->array_z_L_sat,sizeof(double *)*ptsz->n_z_L_sat,ptsz->error_message);
  class_alloc(ptsz->array_m_L_sat,sizeof(double *)*ptsz->n_m_L_sat,ptsz->error_message);


  class_alloc(ptsz->array_L_sat_at_z_and_M_at_nu,
              ptsz->cib_frequency_list_num*sizeof(double *),
              ptsz->error_message);

int index_nu;
for (index_nu=0;index_nu<ptsz->cib_frequency_list_num;index_nu++){

class_alloc(ptsz->array_L_sat_at_z_and_M_at_nu[index_nu],
            ptsz->n_z_L_sat*ptsz->n_m_L_sat*sizeof(double),
            ptsz->error_message);


class_alloc(array_L_sat_at_z_and_M_at_nu,
            ptsz->n_z_L_sat*sizeof(double *),
            ptsz->error_message);



for (index_l=0;
     index_l<ptsz->n_z_L_sat;
     index_l++)
{
  class_alloc(array_L_sat_at_z_and_M_at_nu[index_l],
              ptsz->n_m_L_sat*sizeof(double),
              ptsz->error_message);
}

/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */


int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,index_nu,index_z_M,\
pba,ptsz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


#pragma omp for schedule (dynamic)
for (index_z=0; index_z<ptsz->n_z_L_sat; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<ptsz->n_m_L_sat; index_M++)
{
      ptsz->array_z_L_sat[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_L_sat-1.); // log(1+z)

      ptsz->array_m_L_sat[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(ptsz->n_m_L_sat-1.); //log(R)


      double z =   exp(ptsz->array_z_L_sat[index_z])-1.;
      double logM =   ptsz->array_m_L_sat[index_M];

      double lnMs_min = log(1e10);
      double lnMs_max = logM;//log(1e11);

      double epsrel = ptsz->epsrel_L_sat;
      double epsabs = ptsz->epsabs_L_sat;

      struct Parameters_for_integrand_patterson_L_sat V;
      V.nu = ptsz->cib_frequency_list[index_nu];
      V.z = z;
      V.ptsz = ptsz;
      V.M_host = exp(logM);


      void * params = &V;
      params = &V;

      double L_sat_at_nu = Integrate_using_Patterson_adaptive(lnMs_min, lnMs_max,
                                                               epsrel, epsabs,
                                                               integrand_patterson_L_sat,
                                                               params,ptsz->patterson_show_neval);

      array_L_sat_at_z_and_M_at_nu[index_z][index_M] = log(1.+L_sat_at_nu);
      index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (L_sat) at %.3e GHz = %e s for thread %d\n",
           __func__,ptsz->cib_frequency_list[index_nu],tstop-tstart,omp_get_thread_num());
#endif

    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<ptsz->n_m_L_sat; index_M++)
{
  for (index_z=0; index_z<ptsz->n_z_L_sat; index_z++)
  {
    ptsz->array_L_sat_at_z_and_M_at_nu[index_nu][index_z_M] = array_L_sat_at_z_and_M_at_nu[index_z][index_M];
    index_z_M += 1;
  }
}

for (index_z=0;index_z<ptsz->n_z_L_sat;index_z++){
  free(array_L_sat_at_z_and_M_at_nu[index_z]);
}
  free(array_L_sat_at_z_and_M_at_nu);


}

//exit(0);

return _SUCCESS_;

                                      }



//Tabulate Sigma2(R,z) and dSigma2dR
//as functions of z and logR
int tabulate_sigma_and_dsigma_from_pk(struct background * pba,
                                      struct nonlinear * pnl,
                                      struct primordial * ppm,
                                      struct tszspectrum * ptsz){

if (ptsz->need_sigma == 0)
    return 0;

// printf("tabulating sigma\n");

   // bounds array of radii for sigma computations:
   ptsz->logR1SZ = log(pow(3.*0.1*1e0/(4*_PI_*ptsz->Omega_m_0*ptsz->Rho_crit_0),1./3.));
   ptsz->logR2SZ = log(pow(3.*10.*1e17/(4*_PI_*ptsz->Omega_m_0*ptsz->Rho_crit_0),1./3.));


   // ptsz->logR1SZ = -5.684; // 0.0034Mpc/h, 1.8e4  solar mass
   // ptsz->logR2SZ = 4.; //default =4 , i.e., 54.9Mpc/h, 7.5e16 solar mass


  //Array of z
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,ptsz->z_for_pk_hm);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,ptsz->z_for_pk_hm);
  // z_min = 0.99*z_min;

  int index_z;

  double tstart, tstop;
  int index_l;
  // double * sigma_var;
  // double * dsigma_var;
  int abort;

  //Array of R in Mpc
  double logR_min = log(exp(ptsz->logR1SZ)/pba->h); //in Mpc
  double logR_max = log(exp(ptsz->logR2SZ)/pba->h); //in Mpc
  int index_R;

  int index_z_R = 0;

  double ** array_sigma_at_z_and_R;
  double ** array_dsigma2dR_at_z_and_R;

class_alloc(ptsz->array_redshift,sizeof(double *)*ptsz->n_arraySZ,ptsz->error_message);
class_alloc(ptsz->array_radius,sizeof(double *)*ptsz->ndimSZ,ptsz->error_message);


class_alloc(ptsz->array_sigma_at_z_and_R,
            sizeof(double *)*ptsz->n_arraySZ*ptsz->ndimSZ,
            ptsz->error_message);

class_alloc( ptsz->array_dsigma2dR_at_z_and_R,
            sizeof(double *)*ptsz->n_arraySZ*ptsz->ndimSZ,
            ptsz->error_message);

class_alloc(array_sigma_at_z_and_R,
            ptsz->n_arraySZ*sizeof(double *),
            ptsz->error_message);

class_alloc(array_dsigma2dR_at_z_and_R,
            ptsz->n_arraySZ*sizeof(double *),
            ptsz->error_message);

for (index_l=0;
     index_l<ptsz->n_arraySZ;
     index_l++)
{
  class_alloc(array_sigma_at_z_and_R[index_l],
              ptsz->ndimSZ*sizeof(double),
              ptsz->error_message);

  class_alloc(array_dsigma2dR_at_z_and_R[index_l],
              ptsz->ndimSZ*sizeof(double),
              ptsz->error_message);
}


//Parallelization of Sigma2(R,z) computation
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,index_z_R,\
pba,ptsz,ppm,pnl,z_min,z_max,logR_min,logR_max)\
private(tstart, tstop,index_R,index_z) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

  // class_alloc_parallel(sigma_var,
  //                      sizeof(double *),
  //                      ptsz->error_message);
  //
  // class_alloc_parallel(dsigma_var,
  //                      sizeof(double *),
  //                      ptsz->error_message);


#pragma omp for schedule (dynamic)
  for (index_R=0; index_R<ptsz->ndimSZ; index_R++)
  {
#pragma omp flush(abort)

    double sigma_var,dsigma_var;
    ptsz->array_radius[index_R] =
                                logR_min
                                +index_R*(logR_max-logR_min)
                                /(ptsz->ndimSZ-1.); //log(R)

    for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
    {
      ptsz->array_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_arraySZ-1.); // log(1+z)


  if (ptsz->need_sigma == 1){
      if (ptsz->HMF_prescription_NCDM == 2) //No-pres
        spectra_sigma_for_tSZ(pba,
                              ppm,
                              pnl,
                              ptsz,
                              exp(ptsz->array_radius[index_R]),
                              exp(ptsz->array_redshift[index_z])-1.,
                              &sigma_var//&sigma_at_z_and_R
                              );
      else
        spectra_sigma_ncdm( pba,
                           // spectra_sigma_ncdm( pba,
                           ppm,
                           pnl,
                           ptsz,
                           exp(ptsz->array_radius[index_R]),
                           exp(ptsz->array_redshift[index_z])-1.,
                           &sigma_var//&sigma_at_z_and_R
                           );



      //ptsz->array_sigma_at_z_and_R[index_z_R] = log(*sigma_var);//sigma_at_z_and_R); //log(sigma)
      array_sigma_at_z_and_R[index_z][index_R] = log(sigma_var);//sigma_at_z_and_R); //log(sigma)
      // printf("s=%.6e r=%.6e z=%.6e\n",
      // array_sigma_at_z_and_R[index_z][index_R],
      // exp(ptsz->array_radius[index_R]),
      // exp(ptsz->array_redshift[index_z])-1.);
      // exit(0);


      if (ptsz->HMF_prescription_NCDM == 2) //No-pres
        spectra_sigma_prime( pba,
                            ppm,
                            pnl,
                            ptsz,
                            exp(ptsz->array_radius[index_R]),
                            exp(ptsz->array_redshift[index_z])-1.,
                            &dsigma_var//&dsigma2dR_at_z_and_R
                            );
      else
        spectra_sigma_ncdm_prime( pba,
                                 ppm,
                                 pnl,
                                 ptsz,
                                 exp(ptsz->array_radius[index_R]),
                                 exp(ptsz->array_redshift[index_z])-1.,
                                 &dsigma_var
                                 );



      array_dsigma2dR_at_z_and_R[index_z][index_R] = dsigma_var;

      // array_sigma_at_z_and_R[index_z][index_R] = 0.;
      // array_dsigma2dR_at_z_and_R[index_z][index_R] =  0.;
                         }
      else {
        array_sigma_at_z_and_R[index_z][index_R] = 0.;
        array_dsigma2dR_at_z_and_R[index_z][index_R] =  0.;
      }
      index_z_R += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over R's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    // free(sigma_var);
    // free(dsigma_var);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_R = 0;
for (index_R=0; index_R<ptsz->ndimSZ; index_R++)
{
  for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
  {

    ptsz->array_sigma_at_z_and_R[index_z_R] = array_sigma_at_z_and_R[index_z][index_R];
    double sigma =   ptsz->array_sigma_at_z_and_R[index_z_R];

    // printf("z = %.3e sig = %.3e\n",
    // exp(ptsz->array_redshift[index_z])-1.,
    // ptsz->array_sigma_at_z_and_R[index_z_R]);
    ptsz->array_dsigma2dR_at_z_and_R[index_z_R]=array_dsigma2dR_at_z_and_R[index_z][index_R];
    double dsigma = ptsz->array_sigma_at_z_and_R[index_z_R];
    if (isnan(sigma+dsigma) || isinf(sigma+dsigma) ){
     printf("z=%.3e R=%.3e sigma=%.3e dsigma2dR=%.3e\n",exp(ptsz->array_redshift[index_z]),exp(ptsz->array_radius[index_R]),sigma,dsigma);
     exit(0);
   }
    index_z_R += 1;
  }
}



// freeing memory
for (index_l=0;
     index_l<ptsz->n_arraySZ;
     index_l++)
{
  free(array_sigma_at_z_and_R[index_l]);
  free(array_dsigma2dR_at_z_and_R[index_l]);
}


  free(array_sigma_at_z_and_R);
  free(array_dsigma2dR_at_z_and_R);

return _SUCCESS_;
}


// Tabulate redshift_int_lensmag
// as functions of z
int tabulate_redshift_int_lensmag(struct tszspectrum * ptsz,
                                  struct background * pba){

if (ptsz->has_kSZ_kSZ_lensmag_1halo
  + ptsz->has_gal_lensmag_1h
  + ptsz->has_gal_lensmag_2h
  + ptsz->has_gal_lensmag_hf
  + ptsz->has_lens_lensmag_1h
  + ptsz->has_lens_lensmag_2h
  + ptsz->has_lens_lensmag_hf
  + ptsz->has_tSZ_lensmag_1h
  + ptsz->has_tSZ_lensmag_2h
  + ptsz->has_lensmag_lensmag_1h
  + ptsz->has_lensmag_lensmag_2h
  + ptsz->has_lensmag_lensmag_hf
   == _FALSE_){
    // if (ptsz->sz_verbose>=1)
    // printf("-> Not tabulating Wz for lensing magnification\n");
    return 0;
  }

if (ptsz->sz_verbose>=1){
printf("-> Tabulating Wz for lensing magnification\n");
}
  //Array of z
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  int index_z;
  double ln1pz,z;
  class_alloc(ptsz->array_W_lensmag,sizeof(double *)*ptsz->n_z_W_lensmag,ptsz->error_message);
  class_alloc(ptsz->array_z_W_lensmag,sizeof(double *)*ptsz->n_z_W_lensmag,ptsz->error_message);

  double * pvectsz;
  class_alloc(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);

  double * pvecback;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  for (index_z=0; index_z<ptsz->n_z_W_lensmag; index_z++)
  {
    ln1pz =  log(1.+z_min)
              +index_z*(log(1.+z_max)-log(1.+z_min))
              /(ptsz->n_z_W_lensmag-1.); // log(1+z)

    z = exp(ln1pz) - 1.;

    // set redshift z
    pvectsz[ptsz->index_z] = z;

    int first_index_back = 0;
    double tau;

    // printf("-> start tabulating Wz for lensing magnification\n");

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

    // set chi at redshift z in Mpc/h
    pvectsz[ptsz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);

    // printf("-> doing tabulating Wz for lensing magnification\n");

    // printf("-> Computing integral at z=%.3e\n",z);
    double result;
    redshift_int_lensmag(ptsz,pba,pvectsz,&result);
      // printf("-> 2 doing tabulating Wz for lensing magnification\n");
    if (result <= 0.)
      result = 1e-100;
    ptsz->array_W_lensmag[index_z] = log(result);
    ptsz->array_z_W_lensmag[index_z] = ln1pz;
    // printf("-> integral z = %.3e W = =%.3e\n",z,exp(ptsz->array_W_lensmag[index_z]));
}
if (ptsz->sz_verbose>=1)
printf("-> end tabulating Wz for lensing magnification\n");
 free(pvectsz);
 free(pvecback);
}



int evaluate_redshift_int_lensmag(double * pvectsz,
                                  struct tszspectrum * ptsz)
  {

   double z = pvectsz[ptsz->index_z];
   double z_asked = log(1.+z);

   if (z<exp(ptsz->array_z_W_lensmag[0])-1.)
      z_asked = ptsz->array_z_W_lensmag[0];
   if (z>exp(ptsz->array_z_W_lensmag[ptsz->n_z_W_lensmag-1])-1.)
      z_asked =  ptsz->array_z_W_lensmag[ptsz->n_z_W_lensmag-1];


   pvectsz[ptsz->index_W_lensmag] =  exp(pwl_value_1d(ptsz->n_z_W_lensmag,
                                                        ptsz->array_z_W_lensmag,
                                                        ptsz->array_W_lensmag,
                                                        z_asked));

return _SUCCESS_;
}




// Tabulate dndlnM
// as functions of z and M
int tabulate_dndlnM(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct tszspectrum * ptsz){

  //Array of z
  // double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,ptsz->z_for_pk_hm);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,ptsz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(ptsz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(ptsz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_dndlnM_at_z_and_M;

  class_alloc(ptsz->array_z_dndlnM,sizeof(double *)*ptsz->n_z_dndlnM,ptsz->error_message);
  class_alloc(ptsz->array_m_dndlnM,sizeof(double *)*ptsz->n_m_dndlnM,ptsz->error_message);


class_alloc(ptsz->array_dndlnM_at_z_and_M,
            sizeof(double *)*ptsz->n_z_dndlnM*ptsz->n_m_dndlnM,
            ptsz->error_message);


class_alloc(array_dndlnM_at_z_and_M,
            ptsz->n_z_dndlnM*sizeof(double *),
            ptsz->error_message);


for (index_l=0;
     index_l<ptsz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_dndlnM_at_z_and_M[index_l],
              ptsz->n_m_dndlnM*sizeof(double),
              ptsz->error_message);
}

/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,\
pba,ptsz,ppm,pnl,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_z,index_M,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);

class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);


int i;
for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
{
#pragma omp flush(abort)

      double tau;
      int first_index_back = 0;
      ptsz->array_z_dndlnM[index_z] =
                                log(1.+z_min)
                                +index_z*(log(1.+z_max)-log(1.+z_min))
                                /(ptsz->n_z_dndlnM-1.); // log(1+z)
      double z =   exp(ptsz->array_z_dndlnM[index_z])-1.;

      class_call_parallel(background_tau_of_z(pba,z,&tau),
                 pba->error_message,
                 pba->error_message);

      class_call_parallel(background_at_tau(pba,
                                   tau,
                                   pba->long_info,
                                   pba->inter_normal,
                                   &first_index_back,
                                   pvecback),
                 pba->error_message,
                 pba->error_message);





for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{


      ptsz->array_m_dndlnM[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(ptsz->n_m_dndlnM-1.); //log(R)

      //background quantities @ z:

      double logM =   ptsz->array_m_dndlnM[index_M];


      pvectsz[ptsz->index_z] = z;
      pvectsz[ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);

      double omega = pvecback[pba->index_bg_Omega_m];
      pvectsz[ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      evaluate_HMF_at_logM_and_z(logM,z,pvecback,pvectsz,pba,pnl,ptsz);
      array_dndlnM_at_z_and_M[index_z][index_M] = log(pvectsz[ptsz->index_hmf]);
    }

  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    // free(dndlnM_var);
    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
  {
    ptsz->array_dndlnM_at_z_and_M[index_z_M] = array_dndlnM_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}


// freeing memory:
for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++){
  free(array_dndlnM_at_z_and_M[index_z]);
}
free(array_dndlnM_at_z_and_M);

return _SUCCESS_;
}



///Tabulate m200c_to_m500c conversion
//as functions of z and M
int tabulate_m200c_to_m500c(struct background * pba,
                            struct tszspectrum * ptsz){

  //Array of z
  // double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,ptsz->z_for_pk_hm);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,ptsz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(ptsz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(ptsz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_m200c_to_m500c_at_z_and_M;

  class_alloc(ptsz->array_ln_1pz_m200c_to_m500c,sizeof(double *)*ptsz->n_z_dndlnM,ptsz->error_message);
  class_alloc(ptsz->array_m_m200c_to_m500c,sizeof(double *)*ptsz->n_m_dndlnM,ptsz->error_message);


class_alloc(ptsz->array_m200c_to_m500c_at_z_and_M,
            sizeof(double *)*ptsz->n_z_dndlnM*ptsz->n_m_dndlnM,
            ptsz->error_message);


class_alloc(array_m200c_to_m500c_at_z_and_M,
            ptsz->n_z_dndlnM*sizeof(double *),
            ptsz->error_message);


for (index_l=0;
     index_l<ptsz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_m200c_to_m500c_at_z_and_M[index_l],
              ptsz->n_m_dndlnM*sizeof(double),
              ptsz->error_message);
}

/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,index_z_M,\
pba,ptsz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{
      ptsz->array_ln_1pz_m200c_to_m500c[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_dndlnM-1.);

      ptsz->array_m_m200c_to_m500c[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(ptsz->n_m_dndlnM-1.);

      //background quantities @ z:
      double z =   exp(ptsz->array_ln_1pz_m200c_to_m500c[index_z])-1.;
      double logM =   ptsz->array_m_m200c_to_m500c[index_M];
      pvectsz[ptsz->index_m200c] = exp(logM);
      double tau;
      int first_index_back = 0;


      class_call_parallel(background_tau_of_z(pba,z,&tau),
                 pba->error_message,
                 pba->error_message);

      class_call_parallel(background_at_tau(pba,
                                   tau,
                                   pba->long_info,
                                   pba->inter_normal,
                                   &first_index_back,
                                   pvecback),
                 pba->error_message,
                 pba->error_message);




      pvectsz[ptsz->index_z] = z;
      pvectsz[ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[ptsz->index_Rho_crit];
    double delrho = 200.*rhoc; // 200m
    double delrho_prime = 500.*rhoc; //500c
    double mdel = pvectsz[ptsz->index_m200c];

    double mdel_prime;
    class_call_parallel(mDEL_to_mDELprime(mdel,
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
    pvectsz[ptsz->index_m500c] = mdel_prime;

    array_m200c_to_m500c_at_z_and_M[index_z][index_M] = log(pvectsz[ptsz->index_m500c]);

    index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
  {
    ptsz->array_m200c_to_m500c_at_z_and_M[index_z_M] = array_m200c_to_m500c_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}

// freeing memory
for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
{
free(array_m200c_to_m500c_at_z_and_M[index_z]);
}
free(array_m200c_to_m500c_at_z_and_M);

return _SUCCESS_;
}



///Tabulate m500c_to_m200c conversion
//as functions of z and M
int tabulate_m500c_to_m200c(struct background * pba,
                            struct tszspectrum * ptsz){

  //Array of z
  // double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,ptsz->z_for_pk_hm);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,ptsz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(ptsz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(ptsz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_m200c_to_m500c_at_z_and_M;

  class_alloc(ptsz->array_ln_1pz_m500c_to_m200c,sizeof(double *)*ptsz->n_z_dndlnM,ptsz->error_message);
  class_alloc(ptsz->array_m_m500c_to_m200c,sizeof(double *)*ptsz->n_m_dndlnM,ptsz->error_message);


class_alloc(ptsz->array_m500c_to_m200c_at_z_and_M,
            sizeof(double *)*ptsz->n_z_dndlnM*ptsz->n_m_dndlnM,
            ptsz->error_message);


class_alloc(array_m200c_to_m500c_at_z_and_M,
            ptsz->n_z_dndlnM*sizeof(double *),
            ptsz->error_message);


for (index_l=0;
     index_l<ptsz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_m200c_to_m500c_at_z_and_M[index_l],
              ptsz->n_m_dndlnM*sizeof(double),
              ptsz->error_message);
}


//Parallelization of Sigma2(R,z) computation
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,index_z_M,\
pba,ptsz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{
      ptsz->array_ln_1pz_m500c_to_m200c[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_dndlnM-1.); // log(1+z)

      ptsz->array_m_m500c_to_m200c[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(ptsz->n_m_dndlnM-1.); //log(R)

      //background quantities @ z:
      double z =   exp(ptsz->array_ln_1pz_m500c_to_m200c[index_z])-1.;
      double logM =   ptsz->array_m_m500c_to_m200c[index_M];
      pvectsz[ptsz->index_m500c] = exp(logM);
      double tau;
      int first_index_back = 0;


      class_call_parallel(background_tau_of_z(pba,z,&tau),
                 pba->error_message,
                 pba->error_message);

      class_call_parallel(background_at_tau(pba,
                                   tau,
                                   pba->long_info,
                                   pba->inter_normal,
                                   &first_index_back,
                                   pvecback),
                 pba->error_message,
                 pba->error_message);




      pvectsz[ptsz->index_z] = z;
      pvectsz[ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[ptsz->index_Rho_crit];
    double delrho = 500.*rhoc; // 200m
    double delrho_prime = 200.*rhoc; //500c
    double mdel = pvectsz[ptsz->index_m500c];

    double mdel_prime;
    class_call_parallel(mDEL_to_mDELprime(mdel,
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
    pvectsz[ptsz->index_m200c] = mdel_prime;

    array_m200c_to_m500c_at_z_and_M[index_z][index_M] = log(pvectsz[ptsz->index_m200c]);
    // printf("m = %.3e\n",array_m200m_to_m500c_at_z_and_M[index_z][index_M]);

    index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
  {
    ptsz->array_m500c_to_m200c_at_z_and_M[index_z_M] = array_m200c_to_m500c_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}

// freeing memory
for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
{
free(array_m200c_to_m500c_at_z_and_M[index_z]);
}
free(array_m200c_to_m500c_at_z_and_M);

return _SUCCESS_;
}




///Tabulate m200m_to_m500c conversion
//as functions of z and M
int tabulate_m200m_to_m500c(struct background * pba,
                            struct tszspectrum * ptsz){

  //Array of z
  // double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,ptsz->z_for_pk_hm);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,ptsz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(ptsz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(ptsz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_m200m_to_m500c_at_z_and_M;

  class_alloc(ptsz->array_ln_1pz_m200m_to_m500c,sizeof(double *)*ptsz->n_z_dndlnM,ptsz->error_message);
  class_alloc(ptsz->array_m_m200m_to_m500c,sizeof(double *)*ptsz->n_m_dndlnM,ptsz->error_message);


class_alloc(ptsz->array_m200m_to_m500c_at_z_and_M,
            sizeof(double *)*ptsz->n_z_dndlnM*ptsz->n_m_dndlnM,
            ptsz->error_message);


class_alloc(array_m200m_to_m500c_at_z_and_M,
            ptsz->n_z_dndlnM*sizeof(double *),
            ptsz->error_message);


for (index_l=0;
     index_l<ptsz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_m200m_to_m500c_at_z_and_M[index_l],
              ptsz->n_m_dndlnM*sizeof(double),
              ptsz->error_message);
}


//Parallelization of Sigma2(R,z) computation
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,index_z_M,\
pba,ptsz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{
      ptsz->array_ln_1pz_m200m_to_m500c[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_dndlnM-1.); // log(1+z)

      ptsz->array_m_m200m_to_m500c[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(ptsz->n_m_dndlnM-1.); //log(R)

      //background quantities @ z:
      double z =   exp(ptsz->array_ln_1pz_m200m_to_m500c[index_z])-1.;
      double logM =   ptsz->array_m_m200m_to_m500c[index_M];
      pvectsz[ptsz->index_m200m] = exp(logM);
      double tau;
      int first_index_back = 0;


      class_call_parallel(background_tau_of_z(pba,z,&tau),
                 pba->error_message,
                 pba->error_message);

      class_call_parallel(background_at_tau(pba,
                                   tau,
                                   pba->long_info,
                                   pba->inter_normal,
                                   &first_index_back,
                                   pvecback),
                 pba->error_message,
                 pba->error_message);




      pvectsz[ptsz->index_z] = z;
      pvectsz[ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[ptsz->index_Rho_crit];
    double delrho = 200.*omega*rhoc; // 200m
    double delrho_prime = 500.*rhoc; //500c
    double mdel = pvectsz[ptsz->index_m200m];

    double mdel_prime;
    class_call_parallel(mDEL_to_mDELprime(mdel,
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
    pvectsz[ptsz->index_m500c] = mdel_prime;

    array_m200m_to_m500c_at_z_and_M[index_z][index_M] = log(pvectsz[ptsz->index_m500c]);
    // printf("m = %.3e\n",array_m200m_to_m500c_at_z_and_M[index_z][index_M]);

    index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
  {
    ptsz->array_m200m_to_m500c_at_z_and_M[index_z_M] = array_m200m_to_m500c_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}

// freeing memory
for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
{
free(array_m200m_to_m500c_at_z_and_M[index_z]);
}
free(array_m200m_to_m500c_at_z_and_M);

return _SUCCESS_;
}


///Tabulate m200m_to_m200c conversion
//as functions of z and M
int tabulate_m200m_to_m200c(struct background * pba,
                            struct tszspectrum * ptsz){

  //Array of z
  // double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,ptsz->z_for_pk_hm);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,ptsz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(ptsz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(ptsz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_m200m_to_m200c_at_z_and_M;

  class_alloc(ptsz->array_ln_1pz_m200m_to_m200c,sizeof(double *)*ptsz->n_z_dndlnM,ptsz->error_message);
  class_alloc(ptsz->array_m_m200m_to_m200c,sizeof(double *)*ptsz->n_m_dndlnM,ptsz->error_message);


class_alloc(ptsz->array_m200m_to_m200c_at_z_and_M,
            sizeof(double *)*ptsz->n_z_dndlnM*ptsz->n_m_dndlnM,
            ptsz->error_message);


class_alloc(array_m200m_to_m200c_at_z_and_M,
            ptsz->n_z_dndlnM*sizeof(double *),
            ptsz->error_message);


for (index_l=0;
     index_l<ptsz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_m200m_to_m200c_at_z_and_M[index_l],
              ptsz->n_m_dndlnM*sizeof(double),
              ptsz->error_message);
}


//Parallelization of Sigma2(R,z) computation
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,index_z_M,\
pba,ptsz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{
      ptsz->array_ln_1pz_m200m_to_m200c[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_dndlnM-1.); // log(1+z)

      ptsz->array_m_m200m_to_m200c[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(ptsz->n_m_dndlnM-1.); //log(R)

      //background quantities @ z:
      double z =   exp(ptsz->array_ln_1pz_m200m_to_m200c[index_z])-1.;
      double logM =   ptsz->array_m_m200m_to_m200c[index_M];
      pvectsz[ptsz->index_m200m] = exp(logM);
      double tau;
      int first_index_back = 0;


      class_call_parallel(background_tau_of_z(pba,z,&tau),
                 pba->error_message,
                 pba->error_message);

      class_call_parallel(background_at_tau(pba,
                                   tau,
                                   pba->long_info,
                                   pba->inter_normal,
                                   &first_index_back,
                                   pvecback),
                 pba->error_message,
                 pba->error_message);




      pvectsz[ptsz->index_z] = z;
      pvectsz[ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[ptsz->index_Rho_crit];
    double delrho = 200.*omega*rhoc; // 200m
    double delrho_prime = 200.*rhoc; //500c
    double mdel = pvectsz[ptsz->index_m200m];

    double mdel_prime;
    class_call_parallel(mDEL_to_mDELprime(mdel,
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
    pvectsz[ptsz->index_m200c] = mdel_prime;

    array_m200m_to_m200c_at_z_and_M[index_z][index_M] = log(pvectsz[ptsz->index_m200c]);
    // printf("m = %.3e\n",array_m200m_to_m200c_at_z_and_M[index_z][index_M]);

    index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
  {
    ptsz->array_m200m_to_m200c_at_z_and_M[index_z_M] = array_m200m_to_m200c_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}

// freeing memory
for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
{
free(array_m200m_to_m200c_at_z_and_M[index_z]);
}
free(array_m200m_to_m200c_at_z_and_M);

return _SUCCESS_;
}



///Tabulate m200c_to_m200m conversion
//as functions of z and M
int tabulate_m200c_to_m200m(struct background * pba,
                            struct tszspectrum * ptsz){

  //Array of z
  // double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,ptsz->z_for_pk_hm);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,ptsz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(ptsz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(ptsz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  // int index_z_M = 0;

  // double ** array_m200c_to_m200m_at_z_and_M;

  class_alloc(ptsz->array_ln_1pz_m200c_to_m200m,sizeof(double *)*ptsz->n_z_dndlnM,ptsz->error_message);
  class_alloc(ptsz->array_m_m200c_to_m200m,sizeof(double *)*ptsz->n_m_dndlnM,ptsz->error_message);


class_alloc(ptsz->array_m200c_to_m200m_at_z_and_M,
            sizeof(double *)*ptsz->n_z_dndlnM*ptsz->n_m_dndlnM,
            ptsz->error_message);


// class_alloc(array_m200c_to_m200m_at_z_and_M,
//             ptsz->n_z_dndlnM*sizeof(double *),
//             ptsz->error_message);


// for (index_l=0;
//      index_l<ptsz->n_z_dndlnM;
//      index_l++)
// {
//   class_alloc(array_m200c_to_m200m_at_z_and_M[index_l],
//               ptsz->n_m_dndlnM*sizeof(double),
//               ptsz->error_message);
// }

for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
{
      ptsz->array_ln_1pz_m200c_to_m200m[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_dndlnM-1.); // log(1+z)
}
for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{
      ptsz->array_m_m200c_to_m200m[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(ptsz->n_m_dndlnM-1.); //log(R)
}
//Parallelization of Sigma2(R,z) computation
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,\
pba,ptsz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for collapse(2)
for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
{
for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{

  int index_z_M = index_M * ptsz->n_z_dndlnM + index_z;
  // printf("%d\n",index_z_M);

      //background quantities @ z:
      double z =   exp(ptsz->array_ln_1pz_m200c_to_m200m[index_z])-1.;
      double logM =   ptsz->array_m_m200c_to_m200m[index_M];
      pvectsz[ptsz->index_m200c] = exp(logM);
      double tau;
      int first_index_back = 0;


      class_call_parallel(background_tau_of_z(pba,z,&tau),
                 pba->error_message,
                 pba->error_message);

      class_call_parallel(background_at_tau(pba,
                                   tau,
                                   pba->long_info,
                                   pba->inter_normal,
                                   &first_index_back,
                                   pvecback),
                 pba->error_message,
                 pba->error_message);




      pvectsz[ptsz->index_z] = z;
      pvectsz[ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[ptsz->index_Rho_crit];
    double delrho = 200.*rhoc; // 200c
    double delrho_prime = 200.*omega*rhoc; //200m
    double mdel = pvectsz[ptsz->index_m200c];

    double mdel_prime;
    class_call_parallel(mDEL_to_mDELprime(mdel,
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
    // printf("%.8e\n",mdel_prime);

    ptsz->array_m200c_to_m200m_at_z_and_M[index_z_M] = log(pvectsz[ptsz->index_m200m]);


    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;
}



struct Parameters_for_nl_fitting_function{
  gsl_interp_accel *acc;
  gsl_spline *spline;
};

 int tabulate_nl_index(struct background * pba,
                       struct nonlinear * pnl,
                       struct primordial * ppm,
                       struct tszspectrum * ptsz){

  //Array of z
  double z_min = ptsz->array_redshift[0];
  double z_max = ptsz->array_redshift[ptsz->n_arraySZ-1];
  int index_z;

  double tstart, tstop;
  int index_l;

  // double * pvecback;
  // double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logk_min = ptsz->ln_k_for_tSZ[0]; //in Msun/h
  double logk_max = ptsz->ln_k_for_tSZ[ptsz->ln_k_size_for_tSZ-1]; //in Msun/h
  int index_k;

  int index_z_k = 0;

  double ** array_nl_index_at_z_and_k;
  double ** array_nl_index_at_z_and_k_no_wiggles;



class_alloc(ptsz->array_nl_index_at_z_and_k,
            sizeof(double *)*ptsz->n_arraySZ*ptsz->ln_k_size_for_tSZ,
            ptsz->error_message);
class_alloc(ptsz->array_nl_index_at_z_and_k_no_wiggles,
            sizeof(double *)*ptsz->n_arraySZ*ptsz->ln_k_size_for_tSZ,
            ptsz->error_message);


class_alloc(array_nl_index_at_z_and_k,
            ptsz->n_arraySZ*sizeof(double *),
            ptsz->error_message);
class_alloc(array_nl_index_at_z_and_k_no_wiggles,
            ptsz->n_arraySZ*sizeof(double *),
            ptsz->error_message);


// class_alloc(ptsz->array_nl_index_at_z_and_k_splined,
//             sizeof(double *)*ptsz->n_arraySZ*ptsz->ln_k_size_for_tSZ,
//             ptsz->error_message);
//
//
// class_alloc(array_nl_index_at_z_and_k_splined,
//             ptsz->n_arraySZ*sizeof(double *),
//             ptsz->error_message);

for (index_l=0;
     index_l<ptsz->n_arraySZ;
     index_l++)
{
  class_alloc(array_nl_index_at_z_and_k[index_l],
              ptsz->ln_k_size_for_tSZ*sizeof(double),
              ptsz->error_message);
  class_alloc(array_nl_index_at_z_and_k_no_wiggles[index_l],
              ptsz->ln_k_size_for_tSZ*sizeof(double),
              ptsz->error_message);
}


//Parallelization of Sigma2(R,z) computation
/* initialize error management flag */
abort = _FALSE_;
/* beginning of parallel region */

int number_of_threads= 1;
#ifdef _OPENMP
#pragma omp parallel
  {
    number_of_threads = omp_get_num_threads();
  }
#endif

#pragma omp parallel \
shared(abort,index_z_k,\
pba,pnl,ppm,ptsz,z_min,z_max,logk_min,logk_max)\
private(tstart, tstop,index_k,index_z) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

  //
  // class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
  //
  // class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
{

#pragma omp flush(abort)

for (index_k=0; index_k<ptsz->ln_k_size_for_tSZ; index_k++)
{

      //background quantities @ z:
      double z =   exp(ptsz->array_redshift[index_z])-1.;
      double logk =   ptsz->ln_k_for_tSZ[index_k];


  enum pk_outputs pk_for_nl_index;
  pk_for_nl_index = pk_linear;




    ///////////////////////////////


double result;
double tol=1.e-6;

double lnk1,lnk2;
double pkl1,pkl2;


  lnk1 = logk - tol;
  lnk2 = logk + tol;

  double * pk_ic = NULL;
  double pk;
  double k;

  k = exp(lnk1);
    //Input: wavenumber in 1/Mpc
    //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
   class_call_parallel(nonlinear_pk_at_k_and_z(
                                     pba,
                                     ppm,
                                     pnl,
                                     pk_for_nl_index,
                                     k*pba->h,
                                     z,
                                     pnl->index_pk_cb,
                                     &pk, // number *out_pk_l
                                     pk_ic // array out_pk_ic_l[index_ic_ic]
                                   ),
                                   pnl->error_message,
                                   pnl->error_message);
  pkl1 = pk;

  k = exp(lnk2);
    //Input: wavenumber in 1/Mpc
    //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
   class_call_parallel(nonlinear_pk_at_k_and_z(
                                     pba,
                                     ppm,
                                     pnl,
                                     pk_for_nl_index,
                                     k*pba->h,
                                     z,
                                     pnl->index_pk_cb,
                                     &pk, // number *out_pk_l
                                     pk_ic // array out_pk_ic_l[index_ic_ic]
                                   ),
                                   pnl->error_message,
                                   pnl->error_message);
  pkl2 = pk;

  double dlnpkldlnk = (log(pkl2)-log(pkl1))/2./tol;;

  double nl_index = dlnpkldlnk;
  array_nl_index_at_z_and_k[index_z][index_k] = nl_index;
  // printf("m = %.3e\n",array_m200m_to_m500c_at_z_and_M[index_z][index_M]);
  // printf("z = %.4e \t k = %.4e n = %.3e\n",z,exp(logk),nl_index);

  index_z_k += 1;
    }

  // int i;
  // int index_num;
  // int index_x;
  // int index_y;
  // int index_ddy;
  // i=0;
  // index_x=i;
  // i++;
  // index_y=i;
  // i++;
  // index_ddy=i;
  // i++;
  // index_num=i;

// interpolate the data:
gsl_interp_accel *acc = gsl_interp_accel_alloc ();
gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, ptsz->ln_k_size_for_tSZ);
gsl_spline_init (spline, ptsz->ln_k_for_tSZ, array_nl_index_at_z_and_k[index_z], ptsz->ln_k_size_for_tSZ);
// find where the second derivative changes sign:
int it;
// for (it=0;it<10;it++){
// printf ("%g %g %g\n", ptsz->ln_k_for_tSZ[it], array_nl_index_at_z_and_k[index_z][it],gsl_spline_eval (spline, ptsz->ln_k_for_tSZ[it], acc));
// }
// printf(" ####### ");
// printf("nk = %d\n",ptsz->ln_k_size_for_tSZ);
// exit(0);


struct Parameters_for_nl_fitting_function V;
V.acc = acc;
V.spline = spline;
void * params = &V;

// for (it=0;it<10;it++){
// printf ("nl_fitting_function %g %g %g\n", ptsz->ln_k_for_tSZ[it], array_nl_index_at_z_and_k[index_z][it],nl_fitting_function (ptsz->ln_k_for_tSZ[it], params));
// }
//
// exit(0);

// set-up table of 0's of second derivatives:
int n_data_guess, n_data = 0;
int *lnx_zeros = NULL, *tmp = NULL;
n_data = 0;
n_data_guess = 1;
lnx_zeros   = (int *)malloc(n_data_guess*sizeof(double));

double tol=1.e-6;
double lnkl,lnkc,lnkr;
double nlkl,nlkc,nlkr;
double ddndlnk;
double n_sign_l = -1;
double n_sign_r = -1;

// compute second derivative at the start;
index_k = 0;
// lnkl = ptsz->ln_k_for_tSZ[index_k]+2.*tol-tol;
// lnkr = ptsz->ln_k_for_tSZ[index_k]+2.*tol+tol;
// lnkc = ptsz->ln_k_for_tSZ[index_k]+2.*tol;
tol = 1.e-6;//fabs(ptsz->ln_k_for_tSZ[index_k+1] - ptsz->ln_k_for_tSZ[index_k])/100.;
lnkl = ptsz->ln_k_for_tSZ[index_k]+2.*tol;
lnkc = ptsz->ln_k_for_tSZ[index_k]+tol;
lnkr = ptsz->ln_k_for_tSZ[index_k];

nlkl = nl_fitting_function(lnkl,params);
nlkr = nl_fitting_function(lnkr,params);
nlkc = nl_fitting_function(lnkc,params);


// ddndlnk = (nlkl - 2.*nlkc + nlkr) / tol/tol;


ddndlnk = (nlkc - nlkr)/tol;
ddndlnk =  array_nl_index_at_z_and_k[index_z][index_k+1]-array_nl_index_at_z_and_k[index_z][index_k];
// printf("index_z = %d index_k =  %d dndlnk = %.3e\n",index_z, index_k, ddndlnk );

// printf("fitting = %.3e\n",nl_fitting_function(0.,params));
// //testing
//        int i;
//        double xi, yi, x[10], y[10];
//
//        printf ("#m=0,S=2\n");
//
//        for (i = 0; i < 10; i++)
//          {
//            x[i] = i + 0.5 * sin (i);
//            y[i] = i + cos (i * i);
//            // printf ("%g %g\n", x[i], y[i]);
//          }
//
//        printf ("#m=1,S=0\n");
//
//        {
//          gsl_interp_accel *acc2
//            = gsl_interp_accel_alloc ();
//          gsl_spline *spline2
//            = gsl_spline_alloc (gsl_interp_cspline, 10);
//
//          gsl_spline_init (spline2, x, y, 10);
//          i= 0;
//          for (xi = x[0]; xi < x[9]; xi += 0.01)
//            {
//              yi = gsl_spline_eval (spline2, xi, acc2);
//              // printf ("%g %g\n", xi, yi);
//
//            }
//        for (i = 0; i < 10; i++)
//          {
//            x[i] = i + 0.5 * sin (i);
//            y[i] = i + cos (i * i);
//            printf ("%g %g %g\n", x[i], y[i],gsl_spline_eval (spline2, x[i], acc2));
//          }
//          gsl_spline_free (spline2);
//          gsl_interp_accel_free (acc2);
//        }
//
// exit(0);
if (ddndlnk < 0){
  n_sign_l = -1;
}
else{
  n_sign_l = +1;
}


for (index_k=1; index_k<ptsz->ln_k_size_for_tSZ-1; index_k++){
tol = 1.e-6;//fabs(ptsz->ln_k_for_tSZ[index_k+1] - ptsz->ln_k_for_tSZ[index_k])/100.;
lnkl = ptsz->ln_k_for_tSZ[index_k]+2.*tol;
lnkc = ptsz->ln_k_for_tSZ[index_k]+tol;
lnkr = ptsz->ln_k_for_tSZ[index_k];



nlkl = nl_fitting_function(lnkl,params);
nlkr = nl_fitting_function(lnkr,params);
nlkc = nl_fitting_function(lnkc,params);

// ddndlnk = (nlkl - 2.*nlkc + nlkr) / tol/tol;

ddndlnk = (nlkc - nlkr)/tol;
ddndlnk = array_nl_index_at_z_and_k[index_z][index_k+1]-array_nl_index_at_z_and_k[index_z][index_k];


if (ddndlnk < 0){
  n_sign_r = -1;
}
else{
  n_sign_r = +1;
}

if (n_sign_r*n_sign_l == -1){
  if((n_data+1) > n_data_guess){
    n_data_guess *= 2;
    tmp = (int *)realloc(lnx_zeros,n_data_guess*sizeof(int));
    lnx_zeros = tmp;
    // printf("reallocating memory\n");
  }
  // store the point
  // lnx_zeros[n_data] = ptsz->ln_k_for_tSZ[index_k];
  lnx_zeros[n_data] = index_k;

  n_data++;
}
// printf("index_z = %d index_k =  %d dndlnk = %.3e\n",index_z, index_k, ddndlnk );
  n_sign_l = n_sign_r;
} // continue loop over k
// exit(0);

double array_nl_no_wiggles = 0.;






if (n_data >= 3){
// before entering the oscillations:
int id_zero = 0;
int id_zero_next = 1;
double nod_1_x = (ptsz->ln_k_for_tSZ[lnx_zeros[id_zero_next]] + ptsz->ln_k_for_tSZ[lnx_zeros[id_zero]])/2.;
double nod_1_y = (array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero_next]] + array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero]])/2.;

double nod_2_x = (ptsz->ln_k_for_tSZ[lnx_zeros[id_zero_next+2]] + ptsz->ln_k_for_tSZ[lnx_zeros[id_zero+2]])/2.;
double nod_2_y = (array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero_next+2]] + array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero+2]])/2.;




  for (index_k=0; index_k<ptsz->ln_k_size_for_tSZ; index_k++){
    // if (ptsz->ln_k_for_tSZ[index_k]<=lnx_zeros[0]){
    if (index_k<=lnx_zeros[0]){
      // if (ptsz->ln_k_for_tSZ[index_k] >ptsz->ln_k_for_tSZ[lnx_zeros[id_zero]] - (nod_2_x-nod_1_x)/2. ){

      double  array_nl_no_wiggles_interp = nod_1_y + (nod_2_y-nod_1_y)/(nod_2_x-nod_1_x)*(ptsz->ln_k_for_tSZ[index_k]-nod_1_x);
      double nl_exact = array_nl_index_at_z_and_k[index_z][index_k];
      // }

      // else{
      if (nl_exact<array_nl_no_wiggles_interp & (ptsz->ln_k_for_tSZ[index_k] >ptsz->ln_k_for_tSZ[lnx_zeros[id_zero]] - (nod_2_x-nod_1_x)))
      array_nl_no_wiggles = array_nl_no_wiggles_interp;
      else
      array_nl_no_wiggles = nl_exact;
    // }
    }
    // else if (ptsz->ln_k_for_tSZ[index_k]>=lnx_zeros[n_data-1]){
    else if (index_k+1>=lnx_zeros[n_data-1]){
      array_nl_no_wiggles = array_nl_index_at_z_and_k[index_z][index_k];
    }
    else{

      // linear interpolation:
      // array_nl_no_wiggles = nl_fitting_function(lnx_zeros[id_zero],params)
      //                               + (nl_fitting_function(lnx_zeros[id_zero_next],params)
      //                                  -nl_fitting_function(lnx_zeros[id_zero],params))
      //                               /(lnx_zeros[id_zero_next]-lnx_zeros[id_zero])*(ptsz->ln_k_for_tSZ[index_k]-lnx_zeros[id_zero]);
      nod_1_x = (ptsz->ln_k_for_tSZ[lnx_zeros[id_zero_next]] + ptsz->ln_k_for_tSZ[lnx_zeros[id_zero]])/2.;
      nod_1_y = (array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero_next]] + array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero]])/2.;

      nod_2_x = (ptsz->ln_k_for_tSZ[lnx_zeros[id_zero_next+2]] + ptsz->ln_k_for_tSZ[lnx_zeros[id_zero+2]])/2.;
      nod_2_y = (array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero_next+2]] + array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero+2]])/2.;


      // array_nl_no_wiggles = array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero]]
      //                               + (array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero_next]]
      //                                  -array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero]])
      //                               /(ptsz->ln_k_for_tSZ[lnx_zeros[id_zero_next]]-ptsz->ln_k_for_tSZ[lnx_zeros[id_zero]])*(ptsz->ln_k_for_tSZ[index_k]-ptsz->ln_k_for_tSZ[lnx_zeros[id_zero]]);

      array_nl_no_wiggles = nod_1_y + (nod_2_y-nod_1_y)/(nod_2_x-nod_1_x)*(ptsz->ln_k_for_tSZ[index_k]-nod_1_x);

    // if k larger than next node -> switch range:
    // if (ptsz->ln_k_for_tSZ[index_k+1]==lnx_zeros[id_zero_next]){
    if (index_k+1==lnx_zeros[id_zero_next+1]){
      id_zero_next += 2;
      id_zero += 2;
    }

      // array_nl_no_wiggles = nl_fitting_function(lnx_zeros[0],params)
      //                               + (nl_fitting_function(lnx_zeros[n_data-2],params)-nl_fitting_function(lnx_zeros[0],params))
      //                               /(lnx_zeros[n_data-2]-lnx_zeros[0])*(ptsz->ln_k_for_tSZ[index_k]-lnx_zeros[0]);

    }
   array_nl_index_at_z_and_k_no_wiggles[index_z][index_k] = array_nl_no_wiggles;
  }

}
else{
  for (index_k=0; index_k<ptsz->ln_k_size_for_tSZ; index_k++){
    array_nl_index_at_z_and_k_no_wiggles[index_z][index_k] = array_nl_index_at_z_and_k[index_z][index_k];
  }
}


gsl_spline_free (spline);
gsl_interp_accel_free (acc);
free(lnx_zeros);





  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    // free(pvecback);
    // free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_k = 0;
for (index_k=0; index_k<ptsz->ln_k_size_for_tSZ; index_k++)
{
  for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
  {
    ptsz->array_nl_index_at_z_and_k[index_z_k] = array_nl_index_at_z_and_k[index_z][index_k];
    ptsz->array_nl_index_at_z_and_k_no_wiggles[index_z_k] = array_nl_index_at_z_and_k_no_wiggles[index_z][index_k];
    // printf("index_z = %d index_k =  %d dndlnk_nw = %.3e dndlnk = %.3e\n",index_z, index_k, ptsz->array_nl_index_at_z_and_k_no_wiggles[index_z_k], ptsz->array_nl_index_at_z_and_k[index_z_k]);

    index_z_k += 1;

  }
}

for (index_z=0; index_z<ptsz->n_arraySZ; index_z++){
free(array_nl_index_at_z_and_k[index_z]);
}
  free(array_nl_index_at_z_and_k);
// exit(0);
return _SUCCESS_;
}





double nl_fitting_function(double lnk,void *p){
  struct Parameters_for_nl_fitting_function *V = ((struct Parameters_for_nl_fitting_function *) p);
  double result = gsl_spline_eval(V->spline, lnk, V->acc);
  return result;
}



//Tabulate vrms2 as functions of redshift
int tabulate_sigma2_hsv_from_pk(struct background * pba,
                                struct nonlinear * pnl,
                                struct primordial * ppm,
                                struct tszspectrum * ptsz){


double * sigma2_hsv_var;
class_alloc(sigma2_hsv_var,
            sizeof(double *),
            ptsz->error_message);


class_alloc(ptsz->array_sigma2_hsv_at_z,sizeof(double *)*ptsz->n_arraySZ,ptsz->error_message);

int index_z;

for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
        {

            spectra_sigma2_hsv(pba,
                                ppm,
                                pnl,
                                ptsz,
                                exp(ptsz->array_redshift[index_z])-1.,
                                sigma2_hsv_var
                                );
          ptsz->array_sigma2_hsv_at_z[index_z] = log(*sigma2_hsv_var);
          //printf("%.4e \t %.4e\n",exp(ptsz->array_redshift[index_z])-1.,ptsz->array_sigma2_hsv_at_z[index_z]);
       }

free(sigma2_hsv_var);

return _SUCCESS_;
    }





//Tabulate k non linear as functions of redshift
int tabulate_knl(struct background * pba,
                 struct nonlinear * pnl,
                 struct primordial * ppm,
                 struct tszspectrum * ptsz){



double knl_var;
class_alloc(ptsz->array_knl_at_z,sizeof(double *)*ptsz->n_arraySZ,ptsz->error_message);

int index_z;
double z;
for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
        {
          z = exp(ptsz->array_redshift[index_z])-1.;
          solve_pkl_to_knl(&knl_var,
          z,
          ptsz,
          pba,
          pnl,
          ppm);

          ptsz->array_knl_at_z[index_z] = log(knl_var);
          // printf("z = %.4e \t knl = %.4e\n",z,knl_var);
       }


return _SUCCESS_;
    }





double get_knl_at_z(double z, struct tszspectrum * ptsz){
   double z_asked = log(1.+z);
 if (z<exp(ptsz->array_redshift[0])-1.)
    z_asked = ptsz->array_redshift[0];
 if (z>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
    z_asked =  ptsz->array_redshift[ptsz->n_arraySZ-1];
 return  exp(pwl_value_1d(ptsz->n_arraySZ,
                          ptsz->array_redshift,
                          ptsz->array_knl_at_z,
                          z_asked));
}

double get_nl_index_at_z_and_k(double z_asked, double k_asked, struct tszspectrum * ptsz, struct nonlinear * pnl){
  double z = log(1.+z_asked);
  double k = log(k_asked); // in h/Mpc

 if (z_asked<exp(ptsz->array_redshift[0])-1.)
    z = ptsz->array_redshift[0];
 if (z_asked>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
    z =  ptsz->array_redshift[ptsz->n_arraySZ-1];

 if (k_asked<exp(ptsz->ln_k_for_tSZ[0]))
    k =  ptsz->ln_k_for_tSZ[0];
 if (k_asked>exp(ptsz->ln_k_for_tSZ[ptsz->ln_k_size_for_tSZ-1]))
    k =  ptsz->ln_k_for_tSZ[ptsz->ln_k_size_for_tSZ-1];

 return pwl_interp_2d(ptsz->n_arraySZ,
                      ptsz->ln_k_size_for_tSZ,
                      ptsz->array_redshift,
                      ptsz->ln_k_for_tSZ,
                      ptsz->array_nl_index_at_z_and_k,
                      1,
                      &z,
                      &k);
}
//
double get_nl_index_at_z_and_k_no_wiggles(double z_asked, double k_asked, struct tszspectrum * ptsz, struct nonlinear * pnl){
  double z = log(1.+z_asked);
  double k = log(k_asked); // in h/Mpc

 if (z_asked<exp(ptsz->array_redshift[0])-1.)
    z = ptsz->array_redshift[0];
 if (z_asked>exp(ptsz->array_redshift[ptsz->n_arraySZ-1])-1.)
    z =  ptsz->array_redshift[ptsz->n_arraySZ-1];

 if (k_asked<exp(ptsz->ln_k_for_tSZ[0]))
    k =  ptsz->ln_k_for_tSZ[0];
 if (k_asked>exp(ptsz->ln_k_for_tSZ[ptsz->ln_k_size_for_tSZ-1]))
    k =  ptsz->ln_k_for_tSZ[ptsz->ln_k_size_for_tSZ-1];

 return pwl_interp_2d(ptsz->n_arraySZ,
                      ptsz->ln_k_size_for_tSZ,
                      ptsz->array_redshift,
                      ptsz->ln_k_for_tSZ,
                      ptsz->array_nl_index_at_z_and_k_no_wiggles,
                      1,
                      &z,
                      &k);
}

double get_completeness_at_z_and_M(double z_asked, double m_asked, double * completeness_2d_to_1d, struct tszspectrum * ptsz){
  if (ptsz->has_completeness == 0){
    return 1.;
  }
  else{

  double z = z_asked;

  double m = log(m_asked);
  // printf("z = %.3e m = %.3e logm = %.4e mmin = %.4e mmax = %,4e\n",z_asked,m_asked,m,ptsz->steps_m[0],ptsz->steps_m[ptsz->nsteps_m-1]);
  if (m < ptsz->steps_m[0])
    m = ptsz->steps_m[0];
  if (m > ptsz->steps_m[ptsz->nsteps_m-1])
    m = ptsz->steps_m[ptsz->nsteps_m-1];
  if (z < ptsz->steps_z[0])
    z = ptsz->steps_z[0];
  if (z > ptsz->steps_z[ptsz->nsteps_z-1])
    z = ptsz->steps_z[ptsz->nsteps_z-1];
 return exp(pwl_interp_2d(ptsz->nsteps_m,
                          ptsz->nsteps_z,
                          ptsz->steps_m,
                          ptsz->steps_z,
                          completeness_2d_to_1d,
                          1,
                          &m,
                          &z));

  }
}



double get_detection_proba_at_y_and_theta(double y_asked, double th_asked, double * erfs_2d_to_1d, struct tszspectrum * ptsz){
  double y = log(y_asked);
  double th = log(th_asked);
  double r = 0.;
  // printf("z = %.3e m = %.3e logm = %.4e mmin = %.4e mmax = %,4e\n",z_asked,m_asked,m,ptsz->steps_m[0],ptsz->steps_m[ptsz->nsteps_m-1]);
  if (y < ptsz->erfs_2d_to_1d_y_array[0]){
    r = 0.;}
  else if (y > ptsz->erfs_2d_to_1d_y_array[ptsz->Ny-1]){
    r = 0.;}

  else if (th < ptsz->erfs_2d_to_1d_th_array[0]){
    double th1 = ptsz->erfs_2d_to_1d_th_array[0];
    double th2 = ptsz->erfs_2d_to_1d_th_array[1];
    double r1 = exp(pwl_interp_2d( ptsz->Ny,
                          ptsz->Nth,
                          ptsz->erfs_2d_to_1d_y_array,
                          ptsz->erfs_2d_to_1d_th_array,
                          erfs_2d_to_1d,
                          1,
                          &y,
                          &th1));
    double r2 = exp(pwl_interp_2d( ptsz->Ny,
                          ptsz->nthetas,
                          ptsz->erfs_2d_to_1d_y_array,
                          ptsz->erfs_2d_to_1d_th_array,
                          erfs_2d_to_1d,
                          1,
                          &y,
                          &th2));
    r = r1 + (r2 - r1)/(exp(th2)-exp(th1))*(exp(th)-exp(th1));

  }
  else if (th > ptsz->erfs_2d_to_1d_th_array[ptsz->Nth-1]){
    double th1 = ptsz->erfs_2d_to_1d_th_array[ptsz->Nth-1];
    double th2 = ptsz->erfs_2d_to_1d_th_array[ptsz->Nth-2];
    double r1 = exp(pwl_interp_2d( ptsz->Ny,
                          ptsz->Nth,
                          ptsz->erfs_2d_to_1d_y_array,
                          ptsz->erfs_2d_to_1d_th_array,
                          erfs_2d_to_1d,
                          1,
                          &y,
                          &th1));
    double r2 = exp(pwl_interp_2d( ptsz->Ny,
                          ptsz->Nth,
                          ptsz->erfs_2d_to_1d_y_array,
                          ptsz->erfs_2d_to_1d_th_array,
                          erfs_2d_to_1d,
                          1,
                          &y,
                          &th2));
    r = r1 + (r2 - r1)/(exp(th2)-exp(th1))*(exp(th)-exp(th1));

}

  else{
  // r = exp(pwl_interp_2d(ptsz->nthetas,
  //                         ptsz->Ny,
  //                         ptsz->thetas,
  //                         ptsz->erfs_2d_to_1d_y_array,
  //                         erfs_2d_to_1d,
  //                         1,
  //                         &th,
  //                         &y));
  r = exp(pwl_interp_2d( ptsz->Ny,
                          ptsz->Nth,
                          ptsz->erfs_2d_to_1d_y_array,
                          ptsz->erfs_2d_to_1d_th_array,
                          erfs_2d_to_1d,
                          1,
                          &y,
                          &th));
                        }
  return r;
}

double get_psi_b2t_at_l_and_z(double l_asked, double z_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double l = log(l_asked);

 if (z<ptsz->array_psi_b2t_redshift[0])
    z = ptsz->array_psi_b2t_redshift[0];
 if (z>ptsz->array_psi_b2t_redshift[ptsz->n_z_psi_b2t-1])
    z = ptsz->array_psi_b2t_redshift[ptsz->n_z_psi_b2t-1];

 if (l<ptsz->array_psi_b2t_multipole[0])
    l = ptsz->array_psi_b2t_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l>ptsz->array_psi_b2t_multipole[ptsz->n_l_psi_b2t-1])
    l =  ptsz->array_psi_b2t_multipole[ptsz->n_l_psi_b2t-1];


// printf("l=%.3e\n",l);

 return exp(pwl_interp_2d(

                          ptsz->n_z_psi_b2t,
                          ptsz->n_l_psi_b2t,

                          ptsz->array_psi_b2t_redshift,
                          ptsz->array_psi_b2t_multipole,
                          ptsz->array_psi_b2t_psi,
                          1,
                          &z,
                          &l))-1.;
}



double get_psi_b2g_at_l_and_z(double l_asked, double z_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double l = log(l_asked);

 if (z<ptsz->array_psi_b2g_redshift[0])
    z = ptsz->array_psi_b2g_redshift[0];
 if (z>ptsz->array_psi_b2g_redshift[ptsz->n_z_psi_b2g-1])
    z = ptsz->array_psi_b2g_redshift[ptsz->n_z_psi_b2g-1];

 if (l<ptsz->array_psi_b2g_multipole[0])
    l = ptsz->array_psi_b2g_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l>ptsz->array_psi_b2g_multipole[ptsz->n_l_psi_b2g-1])
    l =  ptsz->array_psi_b2g_multipole[ptsz->n_l_psi_b2g-1];


// printf("l=%.3e\n",l);

 return exp(pwl_interp_2d(

                          ptsz->n_z_psi_b2g,
                          ptsz->n_l_psi_b2g,

                          ptsz->array_psi_b2g_redshift,
                          ptsz->array_psi_b2g_multipole,
                          ptsz->array_psi_b2g_psi,
                          1,
                          &z,
                          &l));
}


double get_psi_b1g_at_l_and_z(double l_asked, double z_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double l = log(l_asked);

 if (z<ptsz->array_psi_b1g_redshift[0])
    z = ptsz->array_psi_b1g_redshift[0];
 if (z>ptsz->array_psi_b1g_redshift[ptsz->n_z_psi_b1g-1])
    z = ptsz->array_psi_b1g_redshift[ptsz->n_z_psi_b1g-1];

 if (l<ptsz->array_psi_b1g_multipole[0])
    l = ptsz->array_psi_b1g_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l>ptsz->array_psi_b1g_multipole[ptsz->n_l_psi_b1g-1])
    l =  ptsz->array_psi_b1g_multipole[ptsz->n_l_psi_b1g-1];


// printf("l=%.3e\n",l);

 return exp(pwl_interp_2d(

                          ptsz->n_z_psi_b1g,
                          ptsz->n_l_psi_b1g,

                          ptsz->array_psi_b1g_redshift,
                          ptsz->array_psi_b1g_multipole,
                          ptsz->array_psi_b1g_psi,
                          1,
                          &z,
                          &l));
}


double get_psi_b1t_at_l_and_z(double l_asked, double z_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double l = log(l_asked);

 if (z<ptsz->array_psi_b1t_redshift[0])
    z = ptsz->array_psi_b1t_redshift[0];
 if (z>ptsz->array_psi_b1t_redshift[ptsz->n_z_psi_b1t-1])
    z = ptsz->array_psi_b1t_redshift[ptsz->n_z_psi_b1t-1];

 if (l<ptsz->array_psi_b1t_multipole[0])
    l = ptsz->array_psi_b1t_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l>ptsz->array_psi_b1t_multipole[ptsz->n_l_psi_b1t-1])
    l =  ptsz->array_psi_b1t_multipole[ptsz->n_l_psi_b1t-1];


 return exp(pwl_interp_2d(

                          ptsz->n_z_psi_b1t,
                          ptsz->n_l_psi_b1t,

                          ptsz->array_psi_b1t_redshift,
                          ptsz->array_psi_b1t_multipole,
                          ptsz->array_psi_b1t_psi,
                          1,
                          &z,
                          &l));
}




double get_psi_b1gt_at_l1_l2_and_z(double l_asked,double l_asked2, double z_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double l1 = log(l_asked);
  double l2 = log(l_asked2);
// printf("z=%.8e\n",z);
 if (z<ptsz->array_psi_b1gt_redshift[0])
    z = ptsz->array_psi_b1gt_redshift[0];
 if (z>ptsz->array_psi_b1gt_redshift[ptsz->n_z_psi_b1gt-1])
    z = ptsz->array_psi_b1gt_redshift[ptsz->n_z_psi_b1gt-1];

 if (l1<ptsz->array_psi_b1gt_multipole[0])
    l1 = ptsz->array_psi_b1gt_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l1>ptsz->array_psi_b1gt_multipole[ptsz->n_l_psi_b1gt-1])
    l1 =  ptsz->array_psi_b1gt_multipole[ptsz->n_l_psi_b1gt-1];

 if (l2<ptsz->array_psi_b1gt_multipole[0])
    l2 = ptsz->array_psi_b1gt_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l2>ptsz->array_psi_b1gt_multipole[ptsz->n_l_psi_b1gt-1])
    l2 =  ptsz->array_psi_b1gt_multipole[ptsz->n_l_psi_b1gt-1];

  // find the closest z's in the grid:
  int id_z_low;
  int id_z_up;
  r8vec_bracket(ptsz->n_z_psi_b1gt,ptsz->array_psi_b1gt_redshift,z,&id_z_low,&id_z_up);


 double ln_rho_low = pwl_interp_2d(ptsz->n_l_psi_b1gt,
                                  ptsz->n_l_psi_b1gt,
                                  ptsz->array_psi_b1gt_multipole,
                                  ptsz->array_psi_b1gt_multipole,
                                  ptsz->array_psi_b1gt_psi[id_z_low-1],
                                  1,
                                  &l1,
                                  &l2);

 double ln_rho_up = pwl_interp_2d(ptsz->n_l_psi_b1gt,
                                  ptsz->n_l_psi_b1gt,
                                  ptsz->array_psi_b1gt_multipole,
                                  ptsz->array_psi_b1gt_multipole,
                                  ptsz->array_psi_b1gt_psi[id_z_up-1],
                                  1,
                                  &l1,
                                  &l2);
 double ln_l_low = ptsz->array_psi_b1gt_redshift[id_z_low-1];
 double ln_l_up = ptsz->array_psi_b1gt_redshift[id_z_up-1];

 return exp(ln_rho_low + ((z - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low));


}


double get_dydz_at_z(double z_asked, struct tszspectrum * ptsz)
{
  double z = log(1.+z_asked);


 if (z<ptsz->array_dydz_redshift[0])
    z = ptsz->array_dydz_redshift[0];
 if (z>ptsz->array_dydz_redshift[ptsz->n_z_dydz-1])
    z = ptsz->array_dydz_redshift[ptsz->n_z_dydz-1];




 return exp(pwl_value_1d(ptsz->n_z_dydz,
                         ptsz->array_dydz_redshift,
                         ptsz->array_dydz_at_z,
                         z));

}




double get_dcib0dz_at_z_and_nu(double z_asked, double nu_asked, struct tszspectrum * ptsz)
{
  double z = log(1.+z_asked);
  double m = log(nu_asked);

 if (z<ptsz->array_dcib0dz_redshift[0])
    z = ptsz->array_dcib0dz_redshift[0];
 if (z>ptsz->array_dcib0dz_redshift[ptsz->n_z_dcib0dz-1])
    z = ptsz->array_dcib0dz_redshift[ptsz->n_z_dcib0dz-1];

 if (m<ptsz->array_dcib0dz_nu[0])
    m = ptsz->array_dcib0dz_nu[0];

 if (m>ptsz->array_dcib0dz_nu[ptsz->n_nu_dcib0dz-1])
    m =  ptsz->array_dcib0dz_nu[ptsz->n_nu_dcib0dz-1];




 return exp(pwl_interp_2d(ptsz->n_z_dcib0dz,
                          ptsz->n_nu_dcib0dz,
                          ptsz->array_dcib0dz_redshift,
                          ptsz->array_dcib0dz_nu,
                          ptsz->array_dcib0dz_at_z_nu,
                          1,
                          &z,
                          &m));

}


double get_dndlnM_at_z_and_M(double z_asked, double m_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);

 if (z<ptsz->array_z_dndlnM[0])
    z = ptsz->array_z_dndlnM[0];
 if (z>ptsz->array_z_dndlnM[ptsz->n_z_dndlnM-1])
    z = ptsz->array_z_dndlnM[ptsz->n_z_dndlnM-1];

 if (m<ptsz->array_m_dndlnM[0])
    m = ptsz->array_m_dndlnM[0];

 if (m>ptsz->array_m_dndlnM[ptsz->n_m_dndlnM-1])
    m =  ptsz->array_m_dndlnM[ptsz->n_m_dndlnM-1];


 return exp(pwl_interp_2d(ptsz->n_z_dndlnM,
                          ptsz->n_m_dndlnM,
                          ptsz->array_z_dndlnM,
                          ptsz->array_m_dndlnM,
                          ptsz->array_dndlnM_at_z_and_M,
                          1,
                          &z,
                          &m));
}

double get_m200m_to_m500c_at_z_and_M(double z_asked, double m_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
 return exp(pwl_interp_2d(ptsz->n_z_dndlnM,
                          ptsz->n_m_dndlnM,
                          ptsz->array_ln_1pz_m200m_to_m500c,
                          ptsz->array_m_m200m_to_m500c,
                          ptsz->array_m200m_to_m500c_at_z_and_M,
                          1,
                          &z,
                          &m));
}


double get_m200m_to_m200c_at_z_and_M(double z_asked, double m_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
 return exp(pwl_interp_2d(ptsz->n_z_dndlnM,
                          ptsz->n_m_dndlnM,
                          ptsz->array_ln_1pz_m200m_to_m200c,
                          ptsz->array_m_m200m_to_m200c,
                          ptsz->array_m200m_to_m200c_at_z_and_M,
                          1,
                          &z,
                          &m));
}


double get_m200c_to_m200m_at_z_and_M(double z_asked, double m_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);

 if (z<ptsz->array_ln_1pz_m200c_to_m200m[0])
    z = ptsz->array_ln_1pz_m200c_to_m200m[0];
 if (z>ptsz->array_ln_1pz_m200c_to_m200m[ptsz->n_z_dndlnM-1])
    z = ptsz->array_ln_1pz_m200c_to_m200m[ptsz->n_z_dndlnM-1];

 if (m<ptsz->array_m_m200c_to_m200m[0])
    m = ptsz->array_m_m200c_to_m200m[0];

 if (m>ptsz->array_m_m200c_to_m200m[ptsz->n_m_dndlnM-1])
    m =  ptsz->array_m_m200c_to_m200m[ptsz->n_m_dndlnM-1];



 return exp(pwl_interp_2d(ptsz->n_z_dndlnM,
                          ptsz->n_m_dndlnM,
                          ptsz->array_ln_1pz_m200c_to_m200m,
                          ptsz->array_m_m200c_to_m200m,
                          ptsz->array_m200c_to_m200m_at_z_and_M,
                          1,
                          &z,
                          &m));
}


double get_m200c_to_m500c_at_z_and_M(double z_asked, double m_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
 return exp(pwl_interp_2d(ptsz->n_z_dndlnM,
                          ptsz->n_m_dndlnM,
                          ptsz->array_ln_1pz_m200c_to_m500c,
                          ptsz->array_m_m200c_to_m500c,
                          ptsz->array_m200c_to_m500c_at_z_and_M,
                          1,
                          &z,
                          &m));
}


double get_m500c_to_m200c_at_z_and_M(double z_asked, double m_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
 return exp(pwl_interp_2d(ptsz->n_z_dndlnM,
                          ptsz->n_m_dndlnM,
                          ptsz->array_ln_1pz_m500c_to_m200c,
                          ptsz->array_m_m500c_to_m200c,
                          ptsz->array_m500c_to_m200c_at_z_and_M,
                          1,
                          &z,
                          &m));
}


double get_T10_alpha_at_z(double z_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  if (z<ptsz->T10_ln1pz[0])
   z = ptsz->T10_ln1pz[0];
  else if (z>ptsz->T10_ln1pz[ptsz->T10_lnalpha_size-1])
   z = ptsz->T10_ln1pz[ptsz->T10_lnalpha_size-1];

double result = exp(pwl_value_1d(ptsz->T10_lnalpha_size,
                          ptsz->T10_ln1pz,
                          ptsz->T10_lnalpha,
                          z));
// printf("z = %.3e  alpha = %.3e\n",z_asked,result);
return result;
}

double get_hmf_counter_term_nmin_at_z(double z_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  if (z<ptsz->array_redshift_hmf_counter_terms[0])
   z = ptsz->array_redshift_hmf_counter_terms[0];
  else if (z>ptsz->array_redshift_hmf_counter_terms[ptsz->n_z_hmf_counter_terms-1])
   z = ptsz->array_redshift_hmf_counter_terms[ptsz->n_z_hmf_counter_terms-1];

double result = pwl_value_1d( ptsz->n_z_hmf_counter_terms,
                              ptsz->array_redshift_hmf_counter_terms,
                              ptsz->array_hmf_counter_terms_nmin,
                              z);
//printf("z = %.3e  alpha = %.3e\n",z_asked,result);
return result;
}



double get_hmf_counter_term_b1min_at_z(double z_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  if (z<ptsz->array_redshift_hmf_counter_terms[0])
   z = ptsz->array_redshift_hmf_counter_terms[0];
  else if (z>ptsz->array_redshift_hmf_counter_terms[ptsz->n_z_hmf_counter_terms-1])
   z = ptsz->array_redshift_hmf_counter_terms[ptsz->n_z_hmf_counter_terms-1];

double result = pwl_value_1d( ptsz->n_z_hmf_counter_terms,
                              ptsz->array_redshift_hmf_counter_terms,
                              ptsz->array_hmf_counter_terms_b1min,
                              z);
//printf("z = %.3e  alpha = %.3e\n",z_asked,result);
return result;
}


double get_hmf_counter_term_b2min_at_z(double z_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  if (z<ptsz->array_redshift_hmf_counter_terms[0])
   z = ptsz->array_redshift_hmf_counter_terms[0];
  else if (z>ptsz->array_redshift_hmf_counter_terms[ptsz->n_z_hmf_counter_terms-1])
   z = ptsz->array_redshift_hmf_counter_terms[ptsz->n_z_hmf_counter_terms-1];

double result = pwl_value_1d( ptsz->n_z_hmf_counter_terms,
                              ptsz->array_redshift_hmf_counter_terms,
                              ptsz->array_hmf_counter_terms_b2min,
                              z);
//printf("z = %.3e  alpha = %.3e\n",z_asked,result);
return result;
}


double get_theta_at_m_and_z(double m, double z, struct tszspectrum * ptsz, struct background * pba){

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

// Eq. 9 of https://arxiv.org/pdf/1303.5080.pdf
double H0 = pba->h*100.;
double Eh = pvecback[pba->index_bg_H]/pba->H0;
double d_A = pvecback[pba->index_bg_ang_distance]*pba->h;
double mp_bias = m/ptsz->HSEbias;
double thetastar2 = ptsz->thetastar * pow(H0/70.,-2./3.);
double theta500_for_mp_at_zp =  thetastar2 * pow(mp_bias/3.e14* (100./H0),ptsz->alpha_theta);
theta500_for_mp_at_zp *=    pow(Eh,-2./3) *pow(100.*d_A/(500.0*H0),-1.);
double thp = theta500_for_mp_at_zp;
free(pvecback);

return thp;
}


double get_volume_at_z(double z, struct background * pba){

  double tau;
  int first_index_back = 0;
  double * pvecback;
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              pba->error_message);

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


double H0 = pba->h*100.;
double Eh = pvecback[pba->index_bg_H]/pba->H0;
double d_A = pvecback[pba->index_bg_ang_distance]*pba->h;
double rz = d_A*(1.+z);
double volume = 3.0e8/1.0e5*rz*rz/Eh;
free(pvecback);

return volume;
}



double get_y_at_m_and_z(double m, double z, struct tszspectrum * ptsz, struct background * pba){

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


double H0 = pba->h*100.;
double Eh = pvecback[pba->index_bg_H]/pba->H0;
double d_A = pvecback[pba->index_bg_ang_distance]*pba->h;
double mp_bias = m/ptsz->HSEbias;
double yp;

if (ptsz->y_m_relation == 1){
        double A = ptsz->A_ym;
        double B = ptsz->B_ym;
        double t = -0.00848*pow(mp_bias/(3.e14*70./(pba->h*100.))*Eh,-0.585);
        double f_rel = 1. + 3.79*t -28.2*t*t;
        yp = A*pow(Eh,2.)*pow(mp_bias/(3.e14*pba->h),1.+B)*f_rel;

        // double a = -1.29389e-01;
        // double b = 9.991387e-01;
        // double c = -3.403211e+01;
        // double d = 1.279992e-01;
        // double e = 8.441768e-01;
        // double M200m_ws = mp_bias/pba->h;
        // yp = 1e-4*exp(a*exp(-pow(log(M200m_ws/3e14/e),2.))
        // +(b*(log(M200m_ws))+c)*(1.-exp(-pow(log(M200m_ws/3e14/d),2.))));



      }
else if (ptsz->y_m_relation == 0){

        double ystar2 = pow(10.,ptsz->ystar_ym)/pow(2., ptsz->alpha_ym)*0.00472724; ////8.9138435358806980e-004;

        ystar2 *=  pow(H0/70.,-2.+ptsz->alpha_ym);
        double y500_for_mp_at_zp =  ystar2 * pow(mp_bias/3.e14* (100./H0),ptsz->alpha_ym);
        y500_for_mp_at_zp *=   pow(Eh,ptsz->beta_ym) *pow(100.*d_A/(500.0*H0),-2.);


        yp = y500_for_mp_at_zp;



}

free(pvecback);

return yp;
}




double  get_L_sat_at_z_M_nu(double z_asked, double m_asked, double nu_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
  double nu = log(nu_asked);


  // double z = log(1.+z_asked);
  // double m = log(m_asked);
   if (z<ptsz->array_z_L_sat[0])
      z = ptsz->array_z_L_sat[0];
        // printf("dealing with mass conversion in hmf\n");
   if (z>ptsz->array_z_L_sat[ptsz->n_z_L_sat-1])
      z =  ptsz->array_z_L_sat[ptsz->n_z_L_sat-1];

   if (m<ptsz->array_m_L_sat[0])
    m = ptsz->array_m_L_sat[0];
      // printf("dealing with mass conversion in hmf\n");
   if (m>ptsz->array_m_L_sat[ptsz->n_m_L_sat-1])
      m =  ptsz->array_m_L_sat[ptsz->n_m_L_sat-1];

   if (nu<ptsz->array_nu_L_sat[0])
    nu = ptsz->array_nu_L_sat[0];
      // printf("dealing with mass conversion in hmf\n");
   if (nu>ptsz->array_nu_L_sat[ptsz->n_nu_L_sat-1])
      nu =  ptsz->array_nu_L_sat[ptsz->n_nu_L_sat-1];

  // if (ptsz->tau_profile == 1){
  // find the closest l's in the grid:
  int id_l_low;
  int id_l_up;
  int n_nu = ptsz->n_nu_L_sat;
  int n_m = ptsz->n_m_L_sat;
  int n_z = ptsz->n_z_L_sat;
  r8vec_bracket(n_nu,ptsz->array_nu_L_sat,nu,&id_l_low,&id_l_up);

  // interpolate 2d at l_low:

 double ln_rho_low = pwl_interp_2d(n_m,
                                n_z,
                                ptsz->array_m_L_sat,
                                ptsz->array_z_L_sat,
                                ptsz->array_L_sat_at_M_z_nu[id_l_low-1],
                                1,
                                &m,
                                &z);

 double ln_rho_up = pwl_interp_2d(n_m,
                                n_z,
                                ptsz->array_m_L_sat,
                                ptsz->array_z_L_sat,
                                ptsz->array_L_sat_at_M_z_nu[id_l_up-1],
                                1,
                                &m,
                                &z);
 double ln_l_low = ptsz->array_nu_L_sat[id_l_low-1];
 double ln_l_up = ptsz->array_nu_L_sat[id_l_up-1];

 // printf("lnrh %.5e %.5e %d %d\n",ln_rho_low,ln_rho_up,id_l_low,id_l_up);

 return exp(ln_rho_low + ((nu - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low))-1.;
 // return ln_rho_low + ((l - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low);


}




double get_L_sat_at_z_and_M_at_nu(double z_asked,
                                  double m_asked,
                                  int index_nu,
                                  struct background * pba,
                                  struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
   if (z<ptsz->array_z_L_sat[0])
      z = ptsz->array_z_L_sat[0];
        // printf("dealing with mass conversion in hmf\n");
   if (z>ptsz->array_z_L_sat[ptsz->n_z_L_sat-1])
      z =  ptsz->array_z_L_sat[ptsz->n_z_L_sat-1];
   if (m<ptsz->array_m_L_sat[0])
    m = ptsz->array_m_L_sat[0];
      // printf("dealing with mass conversion in hmf\n");
   if (m>ptsz->array_m_L_sat[ptsz->n_m_L_sat-1])
      m =  ptsz->array_m_L_sat[ptsz->n_m_L_sat-1];
 // printf("index_nu = %d\n",index_nu);
 return exp(pwl_interp_2d(ptsz->n_z_L_sat,
                          ptsz->n_m_L_sat,
                          ptsz->array_z_L_sat,
                          ptsz->array_m_L_sat,
                          ptsz->array_L_sat_at_z_and_M_at_nu[index_nu],
                          1,
                          &z,
                          &m))-1.;
}

// double get_L_sat_at_z_and_M_at_nu_prime(double z_asked,
//                                   double m_asked,
//                                   struct background * pba,
//                                   struct tszspectrum * ptsz){
//   double z = log(1.+z_asked);
//   double m = log(m_asked);
//  return exp(pwl_interp_2d(ptsz->n_z_L_sat,
//                           ptsz->n_m_L_sat,
//                           ptsz->array_z_L_sat,
//                           ptsz->array_m_L_sat,
//                           ptsz->array_L_sat_at_z_and_M_at_nu_prime,
//                           1,
//                           &z,
//                           &m))-1.;
// }


int bispectrum_condition(double ell_1, double ell_2, double ell_3){
int ell_1_min = abs(ell_2-ell_3);
// TBC
if (ell_1_min<2)
 ell_1_min  = 2;



// full sky:
//if ( (ell_1 >= ell_1_min) && (( (int)floor(ell_1 + ell_2 + ell_3)) % 2 == 0) )
// flat sky:
if ( (ell_1 >= ell_1_min))
  return 1;
else
  return 0;
}
