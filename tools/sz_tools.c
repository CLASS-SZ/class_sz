# include "szpowerspectrum.h"
# include "sz_tools.h"
# include "Patterson.h"
# include "r8lib.h"


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

  double mvir_test;
  double cvir_test;
  double rvir_test;



  EPS2=3.e-8;
  a =x1;
  b =x2;


  mvir_test = exp(a);
  cvir_test = evaluate_cvir_of_mvir(mvir_test,z,ptsz);
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
  cvir_test = evaluate_cvir_of_mvir(mvir_test,z,ptsz);
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
  cvir_test = evaluate_cvir_of_mvir(mvir_test,z,ptsz);
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

  // char Filepath[_ARGUMENT_LENGTH_MAX_];
  // sprintf(Filepath,
  //         "%s%s",
  //         // "%s%s%s",
  //         "cat ",
  //         // ptsz->path_to_class,
  //         "/sz_auxiliary_files/C-M_Zhao09/lnconcentration_vs_z_and_lnm-redshits.txt");
  //
  // process = popen(Filepath, "r");

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
  fclose(process);
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


  // sprintf(Filepath,
  //         "%s%s",
  //         // "%s%s%s",
  //         "cat ",
  //         // ptsz->path_to_class,
  //         "/sz_auxiliary_files/C-M_Zhao09/lnconcentration_vs_z_and_lnm-masses.txt");
  //
  // process = popen(Filepath, "r");

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
  fclose(process);
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

  // sprintf(Filepath,
  //         "%s%s",
  //         // "%s%s%s",
  //         "cat ",
  //         // ptsz->path_to_class,
  //         "/sz_auxiliary_files/C-M_Zhao09/lnconcentration_vs_z_and_lnm.txt");
  //
  // process = popen(Filepath, "r");

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
  //printf("index C-M = %d\n", index);
  // printf("index C-M = %d\n", ptsz->CM_redshift_size*ptsz->CM_logM_size);

  // status = pclose(process);
  fclose(process);


  //for (index_redshift=0;
  //index_redshift<ptsz->CM_redshift_size*ptsz->CM_logM_size;
  // index_redshift++)
  // printf("C-M = %e\n",ptsz->CM_logC[index_redshift]);
  ///printf("C-M = %e\n", logC[ptsz->CM_redshift_size-1][ptsz->CM_logM_size-1]);

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
  //Tabulated values
  //given in Table 3 of 1002.3660v4
  //double z_tab[6] =
  //{
  //0.,0.5,
  //1., 2.,
  //3., 5.
  //};
  //double c0_tab[6] =
  //{
  //9.60, 7.08,
  //5.45, 3.67,
  //2.83, 2.34
  //};
  //double M0_tab[6] =
  //{
  //1.e50, 1.5e17,
  //2.5e15, 6.8e13,
  //6.3e12, 6.6e11
  //};

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
                            struct tszspectrum * ptsz){
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
               struct tszspectrum * ptsz
             ){

double mvir;
//first go from mDEL to mVIR:
class_call(mDEL_to_mVIR(mDEL,
                        delrho,
                        delc,
                        rhoc,
                        z,
                        &mvir,
                        ptsz),
                ptsz->error_message,
                ptsz->error_message);

//then go from mvir to mdel_prime
double rvir = evaluate_rvir_of_mvir(mvir,delc,rhoc,ptsz);
double cvir = evaluate_cvir_of_mvir(mvir,z,ptsz);

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
              struct tszspectrum * ptsz
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
  cvir_test = evaluate_cvir_of_mvir(mvir_test,z,ptsz);
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
      cvir_test = evaluate_cvir_of_mvir(mvir_test,z,ptsz);
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
      cvir_test = evaluate_cvir_of_mvir(mvir_test,z,ptsz);
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
                       ptsz
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

  class_call(background_tau_of_z(V->pba,V->z,&tau),
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

double k_min = 1e-5;
double k_max = 1e1;


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
              pnl->ln_k_size_for_tSZ*index_num*sizeof(double),
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
    //printf("ok z = %e\n",W);

    free(pvecback);



      for (i=0;i<pnl->ln_k_size_for_tSZ;i++) {
        k=exp(pnl->ln_k_for_tSZ[i]);
        if (i == (pnl->ln_k_size_for_tSZ-1)) k *= 0.9999999;

    // //Input: wavenumber in 1/Mpc
    // //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
  enum pk_outputs pk_for_vrms2;
  if (ptsz->pk_nonlinear_for_vrms2 == 1){
    pk_for_vrms2 = pk_nonlinear;
  }
  else {
    pk_for_vrms2 = pk_linear;
  }

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



    array_for_sigma[i*index_num+index_k]=k;
    array_for_sigma[i*index_num+index_y]=pk*W*W;
    //printf("ok k = %e I = %e\n",k,pk*W*W);
  }
//printf("ok z = %e\n",W);
  class_call(array_spline(array_for_sigma,
                          index_num,
                          pnl->ln_k_size_for_tSZ,
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
                                        pnl->ln_k_size_for_tSZ,
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
//printf("ok z = %e\n",W);
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
                       double R,
                       double z,
                       double * sigma
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
              pnl->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);




  for (i=0;i<pnl->ln_k_size_for_tSZ;i++) {
    k=exp(pnl->ln_k_for_tSZ[i]);
    if (i == (pnl->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
    W=3./x/x/x*(sin(x)-x*cos(x));



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


    array_for_sigma[i*index_num+index_k]=k;
    array_for_sigma[i*index_num+index_y]=k*k*pk*W*W;
  }

  class_call(array_spline(array_for_sigma,
                          index_num,
                          pnl->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        pnl->ln_k_size_for_tSZ,
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
                             //struct tszspectrum * ptsz,
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
              pnl->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);

  for (i=0;i<pnl->ln_k_size_for_tSZ;i++) {
    k=exp(pnl->ln_k_for_tSZ[i]);
    if (i == (pnl->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
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
                          pnl->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        pnl->ln_k_size_for_tSZ,
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
  int flag_matter_type;
};


double integrand_nfw_profile(double x, void *p){

  struct Parameters_for_integrand_nfw_profile *V = ((struct Parameters_for_integrand_nfw_profile *) p);

    double nfw_profile_at_x = 0.;

    if (V->flag_matter_type == 1 && V->ptsz->tau_profile == 1){
    rho_gnfw(&nfw_profile_at_x,x,V->pvectsz,V->pba,V->ptsz);
    }
    else{
    rho_nfw(&nfw_profile_at_x,x,V->pvectsz,V->pba,V->ptsz);
    }

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
                          double * result,
                          int flag_matter_type
                          ) {

  struct Parameters_for_integrand_nfw_profile V;
  V.ptsz = ptsz;
  V.pba = pba;
  V.pvectsz = pvectsz;
  V.flag_matter_type = flag_matter_type;

  void * params = &V;

  gsl_function F;
  F.function = &integrand_nfw_profile;
  F.params = params;

  double eps_abs = ptsz->nfw_profile_epsabs;
  double eps_rel = ptsz->nfw_profile_epsrel;

  double result_gsl, error;

  double xin = 1.e-5;
  double rvir = pvectsz[ptsz->index_rVIR]; //in Mpc/h
  double rs = pvectsz[ptsz->index_rs]; //in Mpc/h
  double c_nfw;

  // def at delta
  if (ptsz->MF == 1){

    c_nfw = pvectsz[ptsz->index_c200m];
  }
  else{
    c_nfw = pvectsz[ptsz->index_cVIR];
  }



  double c_nfw_prime;
  if (flag_matter_type == 1){
    // for tau profile, option of rescaling concentration
    c_nfw_prime = ptsz->cvir_tau_profile_factor*c_nfw;
  }
  else{
    c_nfw_prime = c_nfw;
  }

  // with xout = 2.5*rvir/rs the halo model cl^phi^phi matches class cl phi_phi
  // in the settings of KFSW20
  //if ()
  double xout = ptsz->x_out_nfw_profile*c_nfw_prime; //rvir/rs = cvir

  // Battaglia 16 case:
  if (flag_matter_type == 1 && ptsz->tau_profile == 1){
  double rvir = pvectsz[ptsz->index_rVIR]; //in Mpc/h
  // double r200c = pvectsz[ptsz->index_r200c]; //in Mpc/h
  double rs = pvectsz[ptsz->index_rs]; //in Mpc/h
  xout = 50.*rvir/rs; // as in hmvec (default 20, but set to 50 in example file)
}

  //double delta = ptsz->x_out_nfw_profile;
  //double delta_prime = delta_to_delta_prime_nfw(delta,cvir,cvir_prime,ptsz);
  //xout = delta_prime*ptsz->cvir_tau_profile_factor*rvir/rs; //delta_prime*cvir_prime

  //xout = 5.;


// QAWO



  double delta_l = xout - xin;

  gsl_integration_workspace * w;
  gsl_integration_qawo_table * wf;

  int size_w = 500;
  w = gsl_integration_workspace_alloc(size_w);

  //int index_l = (int) pvectsz[ptsz->index_multipole_for_nfw_profile];
  double w0;

  int index_md = (int) pvectsz[ptsz->index_md];
  double y_eff;
  if (_pk_at_z_1h_){
    int index_k = (int) pvectsz[ptsz->index_k_for_pk_hm];
    double k = ptsz->k_for_pk_hm[index_k];
    y_eff = k*rs;
  }
  else{
    y_eff = (pvectsz[ptsz->index_multipole_for_nfw_profile]+0.5)
             /pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile];
  }

  //w0 = (pvectsz[ptsz->index_multipole_for_nfw_profile]+0.5)/pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile];
  //y_eff = 1.e-10;
  w0 = y_eff;


  wf = gsl_integration_qawo_table_alloc(w0, delta_l,GSL_INTEG_SINE,30);


  int limit = size_w; //number of sub interval
  gsl_integration_qawo(&F,xin,eps_abs,eps_rel,limit,w,wf,&result_gsl,&error);

  *result = result_gsl;

  gsl_integration_qawo_table_free(wf);
  gsl_integration_workspace_free(w);

//
// //ROMBERG
// int n_subintervals_gsl = 100;
// gsl_integration_romberg_workspace * w = gsl_integration_romberg_alloc (n_subintervals_gsl);
//
// size_t neval;
// gsl_integration_romberg(&F,xin,xout,eps_abs,eps_rel,&result_gsl,&neval,w);
// gsl_integration_romberg_free(w);
// *result = result_gsl;


}

// Code sample from Colin Hill
// DOUBLE PRECISION FUNCTION rhoNFW(x,y,cvir) ! x=r/rs & y=(l+1/2)/ls & cvir=cvir
//   IMPLICIT none
//   double precision :: x,y,cvir
//
//   rhoNFW=x**(-1d0)*(1d0+x)**(-2d0)*x**2d0*dsin(y*x)/(y*x)
//   return
// END FUNCTION rhoNFW

int rho_nfw(double * rho_nfw_x,
            double x ,
            double * pvectsz,
            struct background * pba,
            struct tszspectrum * ptsz)
{

 int index_md = (int) pvectsz[ptsz->index_md];
 double y_eff;
 if (_pk_at_z_1h_){
   int index_k = (int) pvectsz[ptsz->index_k_for_pk_hm];
   double k = ptsz->k_for_pk_hm[index_k];
   y_eff = k*pvectsz[ptsz->index_rs];
 }
 else{
   y_eff = (pvectsz[ptsz->index_multipole_for_nfw_profile]+0.5)
            /pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile];//*ptsz->cvir_tau_profile_factor;
 }

  //y_eff = 1.e-10;
  *rho_nfw_x = 1./x*1./pow(1.+x,2)*pow(x,2)/(x*y_eff);
               //*sin(x*(pvectsz[ptsz->index_multipole_for_nfw_profile]+0.5)
               ///pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile]);
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
            /pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile];//*ptsz->cvir_tau_profile_factor;


    double A_rho0;
    double A_alpha;
    double A_beta;

    double alpha_m_rho0;
    double alpha_m_alpha;
    double alpha_m_beta;

    double alpha_z_rho0;
    double alpha_z_alpha;
    double alpha_z_beta;

  // Battaglia 16 -- https://arxiv.org/pdf/1607.02442.pdf
  // Table 2
  if (ptsz->tau_profile_mode == 0){
    // agn feedback
    A_rho0 = 4.e3;
    A_alpha = 0.88;
    A_beta = 3.83;

    alpha_m_rho0 = 0.29;
    alpha_m_alpha = -0.03;
    alpha_m_beta = 0.04;

    alpha_z_rho0 = -0.66;
    alpha_z_alpha = 0.19;
    alpha_z_beta = -0.025;
    }
  else if (ptsz->tau_profile_mode == 1){
    // shock heating
    A_rho0 = 1.9e4;
    A_alpha = 0.70;
    A_beta = 4.43;

    alpha_m_rho0 = 0.09;
    alpha_m_alpha = -0.017;
    alpha_m_beta = 0.005;

    alpha_z_rho0 = -0.95;
    alpha_z_alpha = 0.27;
    alpha_z_beta = 0.037;
  }

  // Eq. A1 and A2:
  double m200_over_msol = pvectsz[ptsz->index_m200c]/pba->h; // convert to Msun
  // double rho0  = 1.;
  double rho0 = A_rho0*pow(m200_over_msol/1e14,alpha_m_rho0)*pow(1.+z,alpha_z_rho0);
  double alpha = A_alpha*pow(m200_over_msol/1e14,alpha_m_alpha)*pow(1.+z,alpha_z_alpha);
  double beta = A_beta*pow(m200_over_msol/1e14,alpha_m_beta)*pow(1.+z,alpha_z_beta);

  double gamma = -0.2;
  double xc = 0.5;

  *rho_nfw_x = rho0*pow(x/xc,gamma)*pow(1.+ pow(x/xc,alpha),-(beta+gamma)/alpha)*pow(x,2)/(x*y_eff);

}


double get_density_profile_at_l_M_z(double l_asked, double m_asked, double z_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
  double l = log(l_asked);


  // if (ptsz->tau_profile == 1){
  // find the closest l's in the grid:
  int id_l_low;
  int id_l_up;
  int n_ell = 50;
  int n_m = 100;
  int n_z = 100;
  r8vec_bracket(n_ell,ptsz->array_profile_ln_l,l,&id_l_low,&id_l_up);

  // interpolate 2d at l_low:

 double ln_rho_low = pwl_interp_2d(n_m,
                                n_z,
                                ptsz->array_profile_ln_m,
                                ptsz->array_profile_ln_1pz,
                                ptsz->array_profile_ln_rho_at_lnl_lnM_z[id_l_low],
                                1,
                                &m,
                                &z);

 double ln_rho_up = pwl_interp_2d(n_m,
                                n_z,
                                ptsz->array_profile_ln_m,
                                ptsz->array_profile_ln_1pz,
                                ptsz->array_profile_ln_rho_at_lnl_lnM_z[id_l_up],
                                1,
                                &m,
                                &z);
 double ln_l_low = ptsz->array_profile_ln_l[id_l_low];
 double ln_l_up = ptsz->array_profile_ln_l[id_l_up];

 // return exp(ln_rho_low + ((l - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low));
 return ln_rho_low + ((l - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low);

// }
// else if (){
//   pvectsz[ptsz->index_multipole_for_galaxy_profile] = pvectsz[ptsz->index_multipole_for_nfw_profile]; // this is the multipole going into truncated nfw... TBD: needs to be renamed
//   // set 1 for matter_type = tau
//   result =  evaluate_truncated_nfw_profile(pvectsz,pba,ptsz,1);
//
// }

}




// Tabulate 2D Fourier transform of density profile on a [z - ln_M - ln_ell] grid
int tabulate_density_profile(struct background * pba,
                             struct tszspectrum * ptsz){

if (ptsz->has_kSZ_kSZ_lensmag_1halo
+ ptsz->has_kSZ_kSZ_gal_1h
+ ptsz->has_kSZ_kSZ_gal_2h
+ ptsz->has_kSZ_kSZ_gal_3h
// + ptsz->has_lensmag_lensmag_1h
// + ptsz->has_lensmag_lensmag_2h
// + ptsz->has_lens_lensmag_1h
// + ptsz->has_lens_lensmag_2h
// + ptsz->has_gal_lensmag_1h
// + ptsz->has_gal_lensmag_2h
// + ptsz->has_tSZ_lensmag_1h
// + ptsz->has_tSZ_lensmag_2h
== _FALSE_
)
  return 0;


 // array of multipoles:
 double ln_ell_min = log(0.5);
 double ln_ell_max = log(5e4);
 int n_ell = 50;
 int n_m = 100;
 int n_z = 100;

 class_alloc(ptsz->array_profile_ln_l,sizeof(double *)*n_ell,ptsz->error_message);

 // array of masses:
 double ln_m_min = log(1e8);
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


class_alloc(ptsz->array_profile_ln_rho_at_lnl_lnM_z,n_ell*sizeof(double **),ptsz->error_message);
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
class_alloc(ptsz->array_profile_ln_rho_at_lnl_lnM_z[index_l],n_m*n_z*sizeof(double *),ptsz->error_message);
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
  ptsz->array_profile_ln_rho_at_lnl_lnM_z[index_l][index_m_z] = 1e-100; // initialize with super small number
  index_m_z += 1;
}

     }
}

int has_ksz_bkp = ptsz->has_kSZ_kSZ_gal_1h;
ptsz->has_kSZ_kSZ_gal_1h = _TRUE_; //pretend we need the tau_profile

//Parallelization of profile computation
/* initialize error management flag */

double * pvectsz;
double * pvecback;
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
private(tstart, tstop,index_l,index_z,index_m,index_m_z,pvecback,pvectsz) \
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
  pvectsz[ptsz->index_md] = ptsz->index_md_pk_at_z_1h + 2; // avoid the if condition in rho_nfw

  pvectsz[ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);
  pvectsz[ptsz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
  double omega = pvecback[pba->index_bg_Omega_m];
  pvectsz[ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);


  double result;
   // only  do the tabulation of Battaglia profile
 if (ptsz->tau_profile == 1){
  pvectsz[ptsz->index_m200c] = exp(lnM);
  class_call_parallel(mDEL_to_mVIR(pvectsz[ptsz->index_m200c],
                                   200.*(pvectsz[ptsz->index_Rho_crit]),
                                   pvectsz[ptsz->index_Delta_c],
                                   pvectsz[ptsz->index_Rho_crit],
                                   z,
                                   &pvectsz[ptsz->index_mVIR],
                                   ptsz),
                  ptsz->error_message,
                  ptsz->error_message);

  // rvir needed to cut off the integral --> e.g., xout = 50.*rvir/r200c
  pvectsz[ptsz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[ptsz->index_mVIR],pvectsz[ptsz->index_Delta_c],pvectsz[ptsz->index_Rho_crit],ptsz);
 //compute concentration_parameter using mVIR
  //pvectsz[ ptsz->index_cVIR] = evaluate_cvir_of_mvir(pvectsz[ptsz->index_mVIR],z,ptsz);

  pvectsz[ptsz->index_r200c] = pow(3.*pvectsz[ptsz->index_m200c]/(4.*_PI_*200.*pvectsz[ptsz->index_Rho_crit]),1./3.);
  pvectsz[ptsz->index_l200c] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_r200c];
  // pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile] = pvectsz[ptsz->index_l200c];
  evaluate_c200c_D08(pvecback,pvectsz,pba,ptsz);
  pvectsz[ptsz->index_rs] = pvectsz[ptsz->index_r200c]/pvectsz[ptsz->index_c200c];
  pvectsz[ptsz->index_ls] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_rs];
  pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile] = pvectsz[ptsz->index_ls];


  // matter type = 1 for electron profile
  class_call_parallel(two_dim_ft_nfw_profile(ptsz,pba,pvectsz,&result,1),
                                     ptsz->error_message,
                                     ptsz->error_message);
  result *= pvectsz[ptsz->index_Rho_crit];// normalisation here

}
else if (ptsz->tau_profile == 0){

    pvectsz[ptsz->index_m200m] = exp(lnM);
    class_call_parallel(mDEL_to_mVIR(pvectsz[ptsz->index_m200m],
                                     200.*pvecback[pba->index_bg_Omega_m]*(pvectsz[ptsz->index_Rho_crit]),
                                     pvectsz[ptsz->index_Delta_c],
                                     pvectsz[ptsz->index_Rho_crit],
                                     z,
                                     &pvectsz[ptsz->index_mVIR],
                                     ptsz),
                    ptsz->error_message,
                    ptsz->error_message);

    pvectsz[ptsz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[ptsz->index_mVIR],pvectsz[ptsz->index_Delta_c],pvectsz[ptsz->index_Rho_crit],ptsz);
    //compute concentration_parameter using mVIR
    //pvectsz[ ptsz->index_cVIR] = evaluate_cvir_of_mvir(pvectsz[ptsz->index_mVIR],z,ptsz);

    pvectsz[ptsz->index_r200m] = pow(3.*pvectsz[ptsz->index_m200m]/(4.*_PI_*200.*pvecback[pba->index_bg_Omega_m]*pvectsz[ptsz->index_Rho_crit]),1./3.);
    // pvectsz[ptsz->index_l200m] = sqrt(pvectsz[ptsz->index_chi2])/(1.+z)/pvectsz[ptsz->index_r200m];
    // pvectsz[ptsz->index_characteristic_multipole_for_nfw_profile] = pvectsz[ptsz->index_l200m];
   pvectsz[ptsz->index_multipole_for_galaxy_profile] = pvectsz[ptsz->index_multipole_for_nfw_profile]; // this is the multipole going into truncated nfw... TBD: needs to be renamed
   pvectsz[ptsz->index_md] = ptsz->index_md_kSZ_kSZ_gal_1h; // make sure the mode is set up properly
   evaluate_c200m_D08(pvecback,pvectsz,pba,ptsz);
   pvectsz[ptsz->index_rs] =  pvectsz[ptsz->index_r200m]/pvectsz[ptsz->index_c200m];
   // set 1 for matter_type = tau
   result =  evaluate_truncated_nfw_profile(pvectsz,pba,ptsz,1);
   double tau_normalisation = pvectsz[ptsz->index_m200m]/(4.*_PI_*pow(pvectsz[ptsz->index_rs],3.));
   tau_normalisation *= pba->Omega0_b/ptsz->Omega_m_0/ptsz->mu_e*ptsz->f_free;
   result *= tau_normalisation;
 }

  // ptsz->array_profile_ln_rho_at_lnl_lnM_z[index_l][index_m_z] = log(result);
  ptsz->array_profile_ln_rho_at_lnl_lnM_z[index_l][index_m_z] = result;
  // printf("ell = %.3e z = %.3e m = %.3e lnrho = %.3e\n",ell,z,exp(lnM),log(result));
  index_m_z += 1;
     }


     }



}

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in tab profile parallel region (loop over ell's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif
free(pvectsz);
free(pvecback);
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

// restore initial state:
ptsz->has_kSZ_kSZ_gal_1h = has_ksz_bkp;


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


  ////////////////////////////////
  //Patterson [from Jens Chluba]
  ////////////////////////////////

  if (ptsz->integration_method_pressure_profile==0){

  struct Parameters_for_integrand_patterson_pp V;
  V.ptsz = ptsz;
  V.pba = pba;
  V.pvectsz = pvectsz;
  void * params = &V;
  //
  //
  double epsrel=ptsz->pressure_profile_epsrel;
  double epsabs=ptsz->pressure_profile_epsabs;
  int show_neval = ptsz->patterson_show_neval;
  //

  int id_max = ptsz->x_size_for_pp-1;



  double r=Integrate_using_Patterson_adaptive(ptsz->x_for_pp[0],
                                              ptsz->x_for_pp[id_max],
                                              epsrel, epsabs,
                                              integrand_patterson_test_pp,
                                              params,show_neval);

*result = r;

}
  ///////////////////////////////////
  //end Patterson [from Jens Chluba]
  //////////////////////////////////

  //GSL
  else if (ptsz->integration_method_pressure_profile==1){

  // QAWO
  //int id_max = ptsz->x_size_for_pp-1;
  //double delta_l = ptsz->x_for_pp[id_max] - ptsz->x_for_pp[0];

  double xin = ptsz->x_inSZ;
  double xout = 0.;
    if (ptsz->pressure_profile == 4) { //for Battaglia et al 2012 pressure profile
  double rvir = pvectsz[ptsz->index_rVIR]; //in Mpc/h
  double r200c = pvectsz[ptsz->index_r200c]; //in Mpc/h
  xout = 1.5*rvir/r200c;
}

  else
  xout = ptsz->x_outSZ;

  double delta_l = xout - xin;

  gsl_integration_workspace * w;
  gsl_integration_qawo_table * wf;

  int size_w = 3000;
  w = gsl_integration_workspace_alloc(size_w);

  //int index_l = (int) pvectsz[ptsz->index_multipole_for_pressure_profile];
  double w0;
  if (ptsz->pressure_profile == 4) //for Battaglia et al 2012 pressure profile
  w0 = (pvectsz[ptsz->index_multipole_for_pressure_profile]+0.5)/pvectsz[ptsz->index_l200c];
  else
  w0 =  (pvectsz[ptsz->index_multipole_for_pressure_profile]+0.5)/pvectsz[ptsz->index_l500];

  wf = gsl_integration_qawo_table_alloc(w0, delta_l,GSL_INTEG_SINE,10);

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


/////////

  // ROMBERG
  //
  // gsl_integration_romberg_workspace * w
  //   = gsl_integration_romberg_alloc (30);
  //
  //   double result_gsl, error;
  //
  //
  //   struct Parameters_for_integrand_patterson_pp V;
  //   V.ptsz = ptsz;
  //   V.pba = pba;
  //   V.pvectsz = pvectsz;
  //
  //   void * params = &V;
  //
  //   gsl_function F;
  //   F.function = &integrand_patterson_test_pp;
  //   F.params = params;
  //
  //   int id_max = ptsz->x_size_for_pp-1;
  //
  //   double eps_abs = 1e-8;
  //   double eps_rel = 1e-10;
  //
  //   size_t neval;
  //   gsl_integration_romberg (&F, ptsz->x_for_pp[0],ptsz->x_for_pp[id_max], eps_abs, eps_rel, &result_gsl,&neval,
  //                         w);
  //
  //   gsl_integration_romberg_free (w);
  //
  //   *result = result_gsl;

}
  ///////////////////////////////////
  //end GSL
  //////////////////////////////////

//spline
else if (ptsz->integration_method_pressure_profile==2){
  double * array_for_integral;
  int index_num;
  int index_x;
  int index_y;
  int index_ddy;
  int i;

  double x;


  i=0;
  index_x=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

 double plc_gnfw_at_x = 0.;

  class_alloc(array_for_integral,
              ptsz->x_size_for_pp*index_num*sizeof(double),
              ptsz->error_message);

  for (i=0;i<ptsz->x_size_for_pp;i++) {
    x=ptsz->x_for_pp[i];

    plc_gnfw(&plc_gnfw_at_x,x,pvectsz,pba,ptsz);



    double pp_at_x_and_ell_over_ell_char = plc_gnfw_at_x;
    array_for_integral[i*index_num+index_x]= x;
    array_for_integral[i*index_num+index_y]= pp_at_x_and_ell_over_ell_char;
  }




  class_call(array_spline(array_for_integral,
                          index_num,
                          ptsz->x_size_for_pp,
                          index_x,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          ptsz->error_message),
             ptsz->error_message,
             ptsz->error_message);

  class_call(array_integrate_all_spline(array_for_integral,
                                        index_num,
                                        ptsz->x_size_for_pp,
                                        index_x,
                                        index_y,
                                        index_ddy,
                                        result,
                                        ptsz->error_message),
             ptsz->error_message,
             ptsz->error_message);

  free(array_for_integral);
}

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
              pnl->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);

  for (i=0;i<pnl->ln_k_size_for_tSZ;i++) {
    k=exp(pnl->ln_k_for_tSZ[i]);
    if (i == (pnl->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
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
                          pnl->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        pnl->ln_k_size_for_tSZ,
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
              pnl->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);

  for (i=0;i<pnl->ln_k_size_for_tSZ;i++) {
    k=exp(pnl->ln_k_for_tSZ[i]);
    if (i == (pnl->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
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
                          pnl->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        pnl->ln_k_size_for_tSZ,
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
    class_open(process,ppr->P13_file, "r",ptsz->error_message);
  }
  else if (ptsz->pressure_profile==2){
    if (ptsz->sz_verbose > 0)
      printf("-> Openning the pressure profile file for A10\n");
    //class_open(process,"sz_auxiliary_files/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_A10.txt", "r",ptsz->error_message);
    class_open(process,ppr->A10_file, "r",ptsz->error_message);
    printf("-> File Name: %s\n",ppr->A10_file);
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
  if (ptsz->sz_verbose>=1)
    printf("-> noise curve for yxy covariance not requested\n");
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
    && (ptsz->has_kSZ_kSZ_gal_2h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_3h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_lensmag_1halo != _TRUE_ )
    && (ptsz->has_gal_gal_1h != _TRUE_ )
    && (ptsz->has_gal_lens_1h != _TRUE_ )
    && (ptsz->has_gal_lens_2h != _TRUE_ )
    && (ptsz->has_gal_lensmag_1h != _TRUE_ )
    && (ptsz->has_gal_lensmag_2h != _TRUE_ )
    && (ptsz->has_tSZ_lensmag_1h != _TRUE_ )
    && (ptsz->has_tSZ_lensmag_2h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_1h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_2h != _TRUE_ )
    && (ptsz->has_lens_lensmag_1h != _TRUE_ )
    && (ptsz->has_lens_lensmag_2h != _TRUE_ )
    && (ptsz->has_gal_gal_2h != _TRUE_ ))
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
  class_open(process,ptsz->UNWISE_dndz_file, "r",ptsz->error_message);
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
  sprintf(Filepath,
          "%s%s",
          "cat ",
          ptsz->full_path_to_dndz_gal);
          //"/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/run_scripts/yxg/data/dndz/unwise_red.txt");
  process = popen(Filepath, "r");
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
  if (ptsz->galaxy_sample == 2){
  // if (ptsz->galaxy_sample == 2 || ptsz->galaxy_sample == 0 ){
    status = pclose(process);
  }
  else{
    status = fclose(process);
  }
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
    && (ptsz->has_kSZ_kSZ_gal_2h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_3h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_lensmag_1halo != _TRUE_ )
    && (ptsz->has_gal_gal_1h != _TRUE_ )
    && (ptsz->has_gal_lens_1h != _TRUE_ )
    && (ptsz->has_gal_lens_2h != _TRUE_ )
    && (ptsz->has_gal_lensmag_1h != _TRUE_ )
    && (ptsz->has_gal_lensmag_2h != _TRUE_ )
    && (ptsz->has_tSZ_lensmag_1h != _TRUE_ )
    && (ptsz->has_tSZ_lensmag_2h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_1h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_2h != _TRUE_ )
    && (ptsz->has_lens_lensmag_1h != _TRUE_ )
    && (ptsz->has_lens_lensmag_2h != _TRUE_ )
    && (ptsz->has_gal_gal_2h != _TRUE_ ))
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
  char Filepath[_ARGUMENT_LENGTH_MAX_];

  //unwise

    sprintf(Filepath,
            "%s%s",
            "cat ",
            //ptsz->path_to_class,
            "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/UNWISE_galaxy_ditributions/normalised_fdndz.txt");



  process = popen(Filepath, "r");

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
  status = pclose(process);
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
    && (ptsz->has_kSZ_kSZ_gal_2h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_gal_3h != _TRUE_ )
    && (ptsz->has_kSZ_kSZ_lensmag_1halo != _TRUE_ )
    && (ptsz->has_gal_gal_1h != _TRUE_ )
    && (ptsz->has_gal_lens_1h != _TRUE_ )
    && (ptsz->has_gal_lens_2h != _TRUE_ )
    && (ptsz->has_gal_lensmag_1h != _TRUE_ )
    && (ptsz->has_gal_lensmag_2h != _TRUE_ )
    && (ptsz->has_tSZ_lensmag_1h != _TRUE_ )
    && (ptsz->has_tSZ_lensmag_2h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_1h != _TRUE_ )
    && (ptsz->has_lensmag_lensmag_2h != _TRUE_ )
    && (ptsz->has_lens_lensmag_1h != _TRUE_ )
    && (ptsz->has_lens_lensmag_2h != _TRUE_ )
    && (ptsz->has_gal_gal_2h != _TRUE_ ))
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




int load_unwise_filter(struct tszspectrum * ptsz)
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

    sprintf(Filepath,
            "%s%s",
            "cat ",
            //ptsz->path_to_class,
            "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/UNWISE_galaxy_ditributions/unwise_filter_functions_l_fl.txt");
  process = popen(Filepath, "r");

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
  status = pclose(process);
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
  //   sprintf(Filepath,
  //           "%s%s",
  //           "cat ",
  //           // ptsz->path_to_class,
  //           "/sz_auxiliary_files/Tinker_et_al_10_alpha_consistency_msyriac.txt");
  // process = popen(Filepath, "r");

  class_open(process,"sz_auxiliary_files/Tinker_et_al_10_alpha_consistency_msyriac.txt", "r",ptsz->error_message);


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
  fclose(process);
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


//This routine reads the tabulated
//nfw profiles,
//and stores the tabulated values.

int load_rho_nfw_profile(struct tszspectrum * ptsz)
{

// don't load the lensing profile if lensing/kSZ observables not required
if (ptsz->has_gal_lens_2h != _TRUE_
  && ptsz->has_gal_lens_1h != _TRUE_
  // && ptsz->has_gal_lensmag_2h != _TRUE_
  // && ptsz->has_gal_lensmag_1h != _TRUE_
  // && ptsz->has_lensmag_lensmag_2h != _TRUE_
  // && ptsz->has_lensmag_lensmag_1h != _TRUE_
  && ptsz->has_tSZ_lens_1h != _TRUE_
  && ptsz->has_tSZ_lens_2h != _TRUE_
  && ptsz->has_lens_lens_1h != _TRUE_
  && ptsz->has_lens_lens_2h != _TRUE_
  && ptsz->has_lens_cib_1h != _TRUE_
  && ptsz->has_lens_cib_2h != _TRUE_
//  && ptsz->has_kSZ_kSZ_lensmag_1halo != _TRUE_
  && ptsz->has_kSZ_kSZ_gal_1h != _TRUE_)
  return 0;


  class_alloc(ptsz->RNFW_lnx,sizeof(double *)*100,ptsz->error_message);
  class_alloc(ptsz->RNFW_lnI,sizeof(double *)*100,ptsz->error_message);
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
  // char Filepath[_ARGUMENT_LENGTH_MAX_];
  //
  //   sprintf(Filepath,
  //           "%s%s%s",
  //           "cat ",
  //           ptsz->path_to_class,
  //           "/sz_auxiliary_files/class_sz_lnInfw-vs-lnell-over-ells_xout_5.txt");
  // process = popen(Filepath, "r");

  class_open(process,"sz_auxiliary_files/class_sz_lnInfw-vs-lnell-over-ells_xout_5.txt", "r",ptsz->error_message);

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
  fclose(process);

  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  ptsz->RNFW_lnx_size = n_data;
  /** Make room */

  class_realloc(ptsz->RNFW_lnx,
                ptsz->RNFW_lnx,
                ptsz->RNFW_lnx_size*sizeof(double),
                ptsz->error_message);
  class_realloc(ptsz->RNFW_lnI,
                ptsz->RNFW_lnI,
                ptsz->RNFW_lnx_size*sizeof(double),
                ptsz->error_message);



  /** Store them */
  for (index_x=0; index_x<ptsz->RNFW_lnx_size; index_x++) {
    ptsz->RNFW_lnx[index_x] = lnx[index_x];
    ptsz->RNFW_lnI[index_x] = lnI[index_x];
  };

  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}




//HMF Tinker et al 2008
//interpolated at m500c
int MF_T08_m500(
                double * result,
                double * lognu ,
                double z ,
                struct tszspectrum * ptsz
                )
{
  //T08@m500
  if(z>3.) z=3.;
  double om0 = ptsz->Omega_m_0;
  double ol0 = 1.-ptsz->Omega_m_0;
  double delta_crit = 500.;
  double Omega_m_z = om0*pow(1.+z,3.)/(om0*pow(1.+z,3.)+ ol0);
  double  delta_mean =
  delta_crit/Omega_m_z;

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
      *plc_gnfw_x =  (1./(pow(ptsz->c500*x,ptsz->gammaGNFW)
                    *pow(1.+ pow(ptsz->c500*x,ptsz->alphaGNFW),
                         (ptsz->betaGNFW-ptsz->gammaGNFW)/ptsz->alphaGNFW)))
                    *pow(x,2)
                    //*sin(x*(ptsz->ell[index_l]+0.5)/pvectsz[ptsz->index_l500])
                    /(x*(pvectsz[ptsz->index_multipole_for_pressure_profile]+0.5)/pvectsz[ptsz->index_l500]);


    if (_mean_y_)
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
    if (_mean_y_)
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
  if(z>3.) z=3.;

  double alpha;
  if(ptsz->hm_consistency==0 || ptsz->hm_consistency==1)
  alpha = ptsz->alphaSZ;
  else if (ptsz->hm_consistency==2)
  alpha = get_T10_alpha_at_z(z,ptsz);


  *result = 0.5
            *alpha
            *(1.+pow(pow(ptsz->beta0SZ*pow(1.+z,0.2),2.)
                   *exp(*lognu),
                   -ptsz->phi0SZ
                   *pow(1.+z,-0.08)))*pow(
                    exp(*lognu),
                    ptsz->eta0SZ
                    *pow(1.+z,0.27))
            *exp(-ptsz->gamma0SZ
                 *pow(1.+z,-0.01)
                 *exp(*lognu)/2.)
            *sqrt(exp(*lognu));

  return _SUCCESS_;
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

double next_z(double z_i, double dz){
  // Compute redshift bins where bins are defined with higher resolution for z<0.2
  double dz_i;
  double highres_z_cutoff = 0.2;
  if (z_i < highres_z_cutoff)
    dz_i = 1.e-3;
  //else if ((z_i >= highres_z_cutoff) && (z_i <= 1.))
  else if ((z_i >= highres_z_cutoff) && (z_i <= 1.))
    dz_i = 1.e-2;
  else
    dz_i = dz;

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
    //printf("redshift_int_lensmag = %.3e at z = %.3e\n",redshift_int_lensmag,z);
    V->pvectsz[V->ptsz->index_lensing_Sigma_crit] = _c_/_Mpc_over_m_*_c_/_Mpc_over_m_*pow((1.+z),1.)
                                                   /(4.*_PI_*_G_*_M_sun_/pow(_Mpc_over_m_,3.)
                                                   *sqrt(V->pvectsz[V->ptsz->index_chi2])*redshift_int_lensmag);

  // printf("analytical nfw z = %.3e lm = %.3e\n",z,redshift_int_lensmag);


  }
  else if (
      ((V->ptsz->has_lens_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_1h))
    ||((V->ptsz->has_lens_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_2h))
  )
  {
    evaluate_redshift_int_lensmag(V->pvectsz,V->ptsz);
    double redshift_int_lensmag = V->pvectsz[V->ptsz->index_W_lensmag];
    V->pvectsz[V->ptsz->index_lensing_Sigma_crit] = _c_/_Mpc_over_m_*_c_/_Mpc_over_m_*V->ptsz->chi_star*pow((1.+z),1.)
                                                   /(4.*_PI_*_G_*_M_sun_/pow(_Mpc_over_m_,3.)
                                                   *sqrt(V->pvectsz[V->ptsz->index_chi2])
                                                   *sqrt((V->ptsz->chi_star-sqrt(V->pvectsz[V->ptsz->index_chi2]))*redshift_int_lensmag)
                                                 );
  }
  else {
    // CMB lensing
    // Eq. 6 of https://arxiv.org/pdf/1312.4525.pdf
    V->pvectsz[V->ptsz->index_lensing_Sigma_crit] = _c_/_Mpc_over_m_*_c_/_Mpc_over_m_*V->ptsz->chi_star*pow((1.+z),1.)
                                                   /(4.*_PI_*_G_*_M_sun_/pow(_Mpc_over_m_,3.)
                                                   *sqrt(V->pvectsz[V->ptsz->index_chi2])
                                                   *(V->ptsz->chi_star-sqrt(V->pvectsz[V->ptsz->index_chi2])));

  }


  // V->pvectsz[V->ptsz->index_lensing_Sigma_crit] = 1.6625e18*V->ptsz->chi_star*pow((1.+z),1.) //there is an issue somewhere with (1+z)...
  //                                                 /(sqrt(V->pvectsz[V->ptsz->index_chi2])*(V->ptsz->chi_star-sqrt(V->pvectsz[V->ptsz->index_chi2])));

  // critical density in (Msun/h)/(Mpc/h)^3
  V->pvectsz[V->ptsz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                        *pow(_Mpc_over_m_,1)
                                        *pow(_c_,2)
                                        *V->pvecback[V->pba->index_bg_rho_crit]
                                        /pow(V->pba->h,2);

  double Eh = V->pvecback[V->pba->index_bg_H]/V->ptsz->H0_in_class_units;
  double omega = V->pvecback[V->pba->index_bg_Omega_m];//pow(Eh,2.);
  V->pvectsz[V->ptsz->index_Delta_c] = Delta_c_of_Omega_m(omega);





if (((V->ptsz->has_tSZ_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_gal_1h))
     || ((V->ptsz->has_tSZ_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_gal_2h))
     || ((V->ptsz->has_gal_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_1h))
     || ((V->ptsz->has_gal_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_2h))
     || ((V->ptsz->has_gal_lens_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_1h))
     || ((V->ptsz->has_gal_lens_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_2h))
     || ((V->ptsz->has_gal_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_1h))
     || ((V->ptsz->has_gal_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_2h))
     || ((V->ptsz->has_lensmag_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_lensmag_lensmag_1h))
     || ((V->ptsz->has_lensmag_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lensmag_lensmag_2h))
     || ((V->ptsz->has_lens_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_1h))
     || ((V->ptsz->has_lens_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_2h))
     //|| ((V->ptsz->has_kSZ_kSZ_lensmag_1halo == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_lensmag_1halo))
     || ((V->ptsz->has_kSZ_kSZ_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_1h))
     || ((V->ptsz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_2h))
     || ((V->ptsz->has_kSZ_kSZ_gal_3h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_3h))
    ) {

 evaluate_mean_galaxy_number_density_at_z(V->pvectsz,V->ptsz);

 //printf("z = %.3e ng = %.3e\n",z,V->pvectsz[V->ptsz->index_mean_galaxy_number_density]);
}


  double result = 0.;

  // first deal with quantities that does not require mass integration:
  if ((V->ptsz->has_pk_at_z_2h == _TRUE_) && (index_md == V->ptsz->index_md_pk_at_z_2h)) {

    result = 1.;
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
  else if (
    // galaxy_sample = 1 : unwise
    ((V->ptsz->galaxy_sample==1 && V->ptsz->use_hod == 0) && (V->ptsz->has_gal_lens_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_2h))
  ||((V->ptsz->galaxy_sample==1 && V->ptsz->use_hod == 0) && (V->ptsz->has_gal_lens_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_1h))
  ||((V->ptsz->galaxy_sample==2 && V->ptsz->use_hod == 0) && (V->ptsz->has_gal_lens_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_2h))
  ||((V->ptsz->galaxy_sample==2 && V->ptsz->use_hod == 0) && (V->ptsz->has_gal_lens_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_1h))
  ){

  if (index_md == V->ptsz->index_md_gal_lens_2h){
  evaluate_effective_galaxy_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);

  result = V->pvectsz[V->ptsz->index_halo_bias];

  double W_lens =  radial_kernel_W_lensing_at_z(V->pvecback,
                                                V->pvectsz,
                                                V->pba,
                                                V->ppm,
                                                V->pnl,
                                                V->ptsz);
  // this is needed only in  the approximate calculation
  // for the exact calculation in HOD, this comes out of Sigma_crit
  result *= W_lens;
  }
  else {
    result = 0.;
  }
  }
  // Halofit approach
  else if (
    // galaxy_sample = 1 : unwise
    ((V->ptsz->galaxy_sample==1 && V->ptsz->use_hod == 0) && (V->ptsz->has_gal_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_2h))
  ||((V->ptsz->galaxy_sample==1 && V->ptsz->use_hod == 0) && (V->ptsz->has_gal_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_1h))
  ||((V->ptsz->galaxy_sample==2 && V->ptsz->use_hod == 0) && (V->ptsz->has_gal_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_2h))
  ||((V->ptsz->galaxy_sample==2 && V->ptsz->use_hod == 0) && (V->ptsz->has_gal_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_1h))
  ){

  if (index_md == V->ptsz->index_md_gal_lensmag_2h){
//printf("ok\n");

  evaluate_effective_galaxy_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);

  result = V->pvectsz[V->ptsz->index_halo_bias];

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
    else {
      result = 0.;
    }
  }
  // Halofit approach
  else if (
    // galaxy_sample = 1 : unwise
    ((V->ptsz->galaxy_sample==1 && V->ptsz->use_hod == 0) && (V->ptsz->has_lensmag_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lensmag_lensmag_2h))
  ||((V->ptsz->galaxy_sample==1 && V->ptsz->use_hod == 0) && (V->ptsz->has_lensmag_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_lensmag_lensmag_1h))
  ||((V->ptsz->galaxy_sample==2 && V->ptsz->use_hod == 0) && (V->ptsz->has_lensmag_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lensmag_lensmag_2h))
  ||((V->ptsz->galaxy_sample==2 && V->ptsz->use_hod == 0) && (V->ptsz->has_lensmag_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_lensmag_lensmag_1h))
  ){

  if (index_md == V->ptsz->index_md_lensmag_lensmag_2h){


  // evaluate_effective_galaxy_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
  //
  // result = V->pvectsz[V->ptsz->index_halo_bias];

  double W_lensmag =  radial_kernel_W_lensing_magnification_at_z(V->pvecback,
                                                                V->pvectsz,
                                                                V->pba,
                                                                V->ppm,
                                                                V->pnl,
                                                                V->ptsz);

    // this is needed only in  the approximate calculation
    // for the exact calculation in HOD, this comes out of Sigma_crit
    result = W_lensmag;
    result *= W_lensmag;
    }
    else {
      result = 0.;
    }
  }
  // Halofit approach
  else if (
    // galaxy_sample = 1 : unwise
    ((V->ptsz->galaxy_sample==1 && V->ptsz->use_hod == 0) && (V->ptsz->has_lens_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_2h))
  ||((V->ptsz->galaxy_sample==1 && V->ptsz->use_hod == 0) && (V->ptsz->has_lens_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_1h))
  ||((V->ptsz->galaxy_sample==2 && V->ptsz->use_hod == 0) && (V->ptsz->has_lens_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_2h))
  ||((V->ptsz->galaxy_sample==2 && V->ptsz->use_hod == 0) && (V->ptsz->has_lens_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_lens_lensmag_1h))
  ){

  if (index_md == V->ptsz->index_md_lens_lensmag_2h){


  // evaluate_effective_galaxy_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
  //
  // result = V->pvectsz[V->ptsz->index_halo_bias];

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
    result = W_lensmag*W_lens;

    }
    else {
      result = 0.;
    }
  }

// halofit approach
else if (
  // galaxy_sample = 1 : unwise
  ((V->ptsz->galaxy_sample==1 && V->ptsz->use_hod == 0) && (V->ptsz->has_gal_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_2h))
||((V->ptsz->galaxy_sample==1 && V->ptsz->use_hod == 0) && (V->ptsz->has_gal_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_1h))
||((V->ptsz->galaxy_sample==2 && V->ptsz->use_hod == 0) && (V->ptsz->has_gal_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_2h))
||((V->ptsz->galaxy_sample==2 && V->ptsz->use_hod == 0) && (V->ptsz->has_gal_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_1h))
){
  if (index_md == V->ptsz->index_md_gal_gal_1h) {
    result = 0.;
  }
  else {

evaluate_effective_galaxy_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);

result = V->pvectsz[V->ptsz->index_halo_bias]*V->pvectsz[V->ptsz->index_halo_bias];
}
}
// halofit approach
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
  //printf("integrating over mass\n");
  result = integrate_over_m_at_z(V->pvecback,
                                 V->pvectsz,
                                 V->pba,
                                 V->pnl,
                                 V->ppm,
                                 V->ptsz);
  }



if (((V->ptsz->has_sz_2halo == _TRUE_) && (index_md == V->ptsz->index_md_2halo))
 || ((V->ptsz->has_gal_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_2h)) //## BB debug
 || ((V->ptsz->has_cib_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_cib_cib_2h))
 || ((V->ptsz->has_tSZ_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_tSZ_cib_2h))
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
 // || ((V->ptsz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_2h))
    ){


  int index_l = (int) V->pvectsz[V->ptsz->index_multipole];
  V->pvectsz[V->ptsz->index_multipole_for_pk] = V->ptsz->ell[index_l];
  evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);

  // For all the above cases we multiply the linear matter power spectrum to the redshift integrand
  // evaluated at (ell+1/2)/Chi and redshift z
  result *= V->pvectsz[V->ptsz->index_pk_for_halo_bias];

}

// // Halofit correction facto as in KAP20:
// // [commented out for now: this was incorrect!!
// if (((V->ptsz->has_gal_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_2h) && (V->ptsz->galaxy_sample==0)) // WIxSC
//     || ((V->ptsz->has_gal_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_1h) && (V->ptsz->galaxy_sample==0)) // WIxSC
//     )
// {
//
// result *= evaluate_pk_halofit_over_pk_linear_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
//
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
  || ((V->ptsz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_2h))
  || ((V->ptsz->has_kSZ_kSZ_gal_3h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_3h))
  || ((V->ptsz->has_gal_lens_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_2h))
  || ((V->ptsz->has_gal_lens_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lens_1h))
  || ((V->ptsz->has_gal_lensmag_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_2h))
  || ((V->ptsz->has_gal_lensmag_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_lensmag_1h))
    ){
// multiply by radial kernel for galaxies
double Wg = radial_kernel_W_galaxy_at_z(V->pvecback,V->pvectsz,V->pba,V->ptsz);
result *= Wg/V->pvectsz[V->ptsz->index_chi2];

}

// multiply by velocity dispersion
if ((V->ptsz->has_kSZ_kSZ_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_1h)
 || (V->ptsz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_2h)
 || (V->ptsz->has_kSZ_kSZ_gal_3h == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_gal_3h)
 || (V->ptsz->has_kSZ_kSZ_lensmag_1halo == _TRUE_) && (index_md == V->ptsz->index_md_kSZ_kSZ_lensmag_1halo)){
  evaluate_vrms2(V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);
  result *= V->pvectsz[V->ptsz->index_vrms2]/3./pow(_c_*1e-3,2.);
}

// gxg needs Wg^2:
if ( ((V->ptsz->has_gal_gal_1h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_1h))
   ||((V->ptsz->has_gal_gal_2h == _TRUE_) && (index_md == V->ptsz->index_md_gal_gal_2h))
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
    ||((V->ptsz->has_lens_cib_1h == _TRUE_) && (index_md == V->ptsz->index_md_lens_cib_1h))
    ||((V->ptsz->has_lens_cib_2h == _TRUE_) && (index_md == V->ptsz->index_md_lens_cib_2h))
    ){

  result *= 1./(1.+z)*pow(1./V->pvectsz[V->ptsz->index_chi2],1.);
}

// if computing 3d matter power spectrum P(k):
if( ((V->ptsz->has_pk_at_z_1h == _TRUE_) && (index_md == V->ptsz->index_md_pk_at_z_1h))
 || ((V->ptsz->has_pk_at_z_2h == _TRUE_) && (index_md == V->ptsz->index_md_pk_at_z_2h))){

    if( ((V->ptsz->has_pk_at_z_1h == _TRUE_) && (index_md == V->ptsz->index_md_pk_at_z_1h))){
      //result *= pow((V->pba->Omega0_cdm+V->pba->Omega0_b)*pow(1.+z,3.)*V->pvectsz[V->ptsz->index_Rho_crit],-2);
      result *= pow((V->pba->Omega0_cdm+V->pba->Omega0_b)*V->ptsz->Rho_crit_0,-2);
      // proper to comoving volume?
      result *= pow(1.+z,-2);
      //result *= pow((V->pba->Omega0_cdm+V->pba->Omega0_b)*V->ptsz->Rho_crit_0,-2);
    }

   return result;
 }// end P(k)

else{

  // finally multiply by volume element Chi^2 dChi
  result *= V->pvectsz[V->ptsz->index_chi2];

  // integrate w.r.t ln(1+z); dz =  (1+z)dln(1+z)
  // volume element in units h^-3 Mpc^3
  // volume = dv/(dzdOmega)
  // Chi^2 dChi = dV/(dzdOmega) dz
  // Chi^2 dChi = dV/(dzdOmega)*(1+z) dln(1+z)
  // dChi = (c/H) *(1+z) dln(1+z) ---> this is used
  // dChi = (c/H) dz
  double H_over_c_in_h_over_Mpc = V->pvecback[V->pba->index_bg_H]/V->pba->h;
  result = (1.+V->pvectsz[V->ptsz->index_z])*result/H_over_c_in_h_over_Mpc;



  // test agains KA20
  //
  // result = (1.+V->pvectsz[V->ptsz->index_z])/H_over_c_in_h_over_Mpc/V->pba->h/(V->pvectsz[V->ptsz->index_chi2]/V->pba->h/V->pba->h);
  // double Wg = radial_kernel_W_galaxy_at_z(V->pvecback,V->pvectsz,V->pba,V->ptsz);
  //
  // result *= pow(Wg*V->pba->h,2.);
  //
  //
  // evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
  //
  // // For all the above cases we multiply the linear matter power spectrum to the redshift integrand
  // // evaluated at (ell+1/2)/Chi and redshift z
  // result *= V->pvectsz[V->ptsz->index_pk_for_halo_bias]/pow(V->pba->h,3.);
  //  evaluate_mean_galaxy_number_density_at_z(V->pvectsz,V->ptsz);
  // result *= 1./(V->pvectsz[V->ptsz->index_mean_galaxy_number_density]*V->pvectsz[V->ptsz->index_mean_galaxy_number_density]);
  // //result *= pow(V->pba->h,-6.);
  //
  //
  // result *= integrate_over_m_at_z(V->pvecback,
  //                                 V->pvectsz,
  //                                 V->pba,
  //                                 V->pnl,
  //                                 V->ppm,
  //                                 V->ptsz);

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

  //printf("integrating between z_min=%.3e and z_max=%.3e\n",z_min,z_max);

  struct Parameters_for_integrand_redshift V;
  V.pnl = pnl;
  V.ppm = ppm;
  V.ptsz = ptsz;
  V.pba = pba;
  V.pvectsz = Pvectsz;
  V.pvecback = Pvecback;
  void * params = &V;
  double r; //result of the integral

  double epsrel= ptsz->redshift_epsrel;//ptsz->patterson_epsrel;
  double epsabs= ptsz->redshift_epsabs;//ptsz->patterson_epsabs;
  int show_neval = ptsz->patterson_show_neval;

  int index_md = (int) Pvectsz[ptsz->index_md];
if(_pk_at_z_1h_ || _pk_at_z_2h_)
{
  r = integrand_redshift(log(1. + ptsz->z_for_pk_hm),params);
}
else{
  r=Integrate_using_Patterson_adaptive(log(1. + z_min), log(1. + z_max),
                                       epsrel, epsabs,
                                       integrand_redshift,
                                       params,show_neval);
    }

//   ///////
// gsl_function F;
// F.function = &integrand_redshift;
// F.params = params;
//
// int n_subintervals_gsl = 30;
//
//
// gsl_integration_romberg_workspace * w = gsl_integration_romberg_alloc (n_subintervals_gsl);
//
// double result_gsl;
// size_t neval;
// gsl_integration_romberg(&F,log(1. + z_min), log(1. + z_max),epsabs,epsrel,&result_gsl,&neval,w);
// gsl_integration_romberg_free(w);
//
// r = result_gsl;
//
//
//   /////////

  Pvectsz[ptsz->index_integral] = r;
  //int index_l = (int)  Pvectsz[ptsz->index_multipole];
  //printf("result=%.4e z_min = %.3e z_max=%.3e\n",r,log(1. + z_min),log(1. + z_max));
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

  if ( ((int) pvectsz[ptsz->index_md] == ptsz->index_md_cov_Y_N )|| ((int) pvectsz[ptsz->index_md] == ptsz->index_md_cov_Y_N_next_order )|| ((int) pvectsz[ptsz->index_md] == ptsz->index_md_cov_N_N ) || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_cov_N_N_hsv )){
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

//   if ((int) pvectsz[ptsz->index_md] == ptsz->index_md_isw_lens ){
//
//   double logM = 1.;
//   r = integrand_patterson_test(logM, params);
//
//     pvectsz[ptsz->index_integral_over_m] = r;
//
// }


   // else {

  // //Patterson [Jens Chluba]
  // if (ptsz->integration_method_mass==0){


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

 // // galaxy cases
 //  else if ( ((int) pvectsz[ptsz->index_md] == ptsz->index_md_tSZ_gal_1h )
 //         || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_gal_gal_1h)
 //         || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_tSZ_gal_2h)
 //         || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_gal_gal_2h)
 //         || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_gal_lens_1h)
 //         || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_gal_lens_2h)
 //         || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_kSZ_kSZ_gal_1h)){


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
   double bmin_umin = bmin*integrand_patterson_test(log(ptsz->m_min_counter_terms),params)/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
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
   double bmin_umin = bmin*integrand_patterson_test(log(ptsz->m_min_counter_terms),params)/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
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

//  if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
//    double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
//    double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
//    double bmin_umin = bmin*integrand_patterson_test(log(ptsz->m_min_counter_terms),params)/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
//    r_m_1 += bmin_umin;
//    // printf("counter terms done r_m_1\n");
// }


  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('galaxy' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

 //   if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
 //   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
 //   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
 //   double bmin_umin = bmin*integrand_patterson_test(log(ptsz->m_min_counter_terms),params)/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
 //   r_m_2 += bmin_umin;
 //   // printf("counter terms done r_m_2\n");
 // }


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

  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Phi' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  // corrected to match Mat M. consistency treatment in hmvec
  r_m_2 = r_m_2+ptsz->Omega_m_0*ptsz->Rho_crit_0*pow(pvecback[pba->index_bg_ang_distance]*pba->h,-2.)/pvectsz[ptsz->index_lensing_Sigma_crit];


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

  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_21=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  // r_m_12*r_m_22
  pvectsz[ptsz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_12=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  pvectsz[ptsz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_22=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  // r_m_13*r_m_23
  pvectsz[ptsz->index_part_id_cov_hsv] = 5;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_12=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  pvectsz[ptsz->index_part_id_cov_hsv] = 6;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_22=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);




  int index_l_1 = (int) pvectsz[ptsz->index_multipole_1];
  pvectsz[ptsz->index_multipole_for_pk] = ptsz->ell_kSZ2_gal_multipole_grid[index_l_1];
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk1 = pvectsz[ptsz->index_pk_for_halo_bias];

  int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
  pvectsz[ptsz->index_multipole_for_pk] = ptsz->ell_kSZ2_gal_multipole_grid[index_l_2];
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk2 = pvectsz[ptsz->index_pk_for_halo_bias];

  int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
  pvectsz[ptsz->index_multipole_for_pk] = ptsz->ell[index_l_3];
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk3 = pvectsz[ptsz->index_pk_for_halo_bias];


  r = pk3*r_m_11*r_m_21  +  pk2*r_m_12*r_m_22  +  pk1*r_m_13*r_m_23;
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

  pvectsz[ptsz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  // r_m_12*r_m_22
  pvectsz[ptsz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);

  pvectsz[ptsz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_patterson_test,
                                           params,ptsz->patterson_show_neval);



  int index_l_1 = (int) pvectsz[ptsz->index_multipole_1];
  pvectsz[ptsz->index_multipole_for_pk] = ptsz->ell_kSZ2_gal_multipole_grid[index_l_1];
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk1 = pvectsz[ptsz->index_pk_for_halo_bias];

  int index_l_2 = (int) pvectsz[ptsz->index_multipole_2];
  pvectsz[ptsz->index_multipole_for_pk] = ptsz->ell_kSZ2_gal_multipole_grid[index_l_2];
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk2 = pvectsz[ptsz->index_pk_for_halo_bias];

  int index_l_3 = (int) pvectsz[ptsz->index_multipole_3];
  pvectsz[ptsz->index_multipole_for_pk] = ptsz->ell[index_l_3];
  evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,ptsz);
  double pk3 = pvectsz[ptsz->index_pk_for_halo_bias];


  double z = pvectsz[ptsz->index_z];
  double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z);
  double k1 = (ptsz->ell_kSZ2_gal_multipole_grid[index_l_1] + 0.5)/d_A;
  double k2 = (ptsz->ell_kSZ2_gal_multipole_grid[index_l_2] + 0.5)/d_A;
  double k3 = (ptsz->ell[index_l_3] + 0.5)/d_A;
  double Fk1k2 = bispectrum_f2_kernel(k1,k2,k3);
  double Fk1k3 = bispectrum_f2_kernel(k3,k1,k2);
  double Fk2k3 = bispectrum_f2_kernel(k2,k3,k1);

  double comb_pks = pk1*pk2+pk1*pk3+pk2*pk3;
  double comb_pks_fks = 2.*pk1*pk2*Fk1k2+2.*pk1*pk3*Fk1k3+2.*pk2*pk3*Fk2k3;

  // printf("comb_pks = %.3e comb_pks_fks = %.3e\n",comb_pks,comb_pks_fks);
  // printf("r_m_b1t1 = %.3e r_m_b1t2 = %.3e r_m_b1g3 = %.3e r_m_b2g3 = %.3e\n",r_m_b1t1,r_m_b1t2,r_m_b1g3,r_m_b2g3);

  r = r_m_b1t1*r_m_b1t2*r_m_b1g3*comb_pks_fks+r_m_b1t1*r_m_b1t2*r_m_b2g3*comb_pks;

  }






  //pvectsz[ptsz->index_part_id_cov_hsv] = 1;
  //V.pvectsz = pvectsz;
  //params = &V;
  // else {
  //   // integrate over the mass range
  //   r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
  //                                            epsrel, epsabs,
  //                                            integrand_patterson_test,
  //                                            params,ptsz->patterson_show_neval);}
  //
  //                                }

  else {
  r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                       epsrel, epsabs,
                                       integrand_patterson_test,
                                       params,ptsz->patterson_show_neval);

// evaluate low mass part:
  if (ptsz->M1SZ == ptsz->m_min_counter_terms)  {
 double rho0 = (pba->Omega0_cdm+pba->Omega0_b)*ptsz->Rho_crit_0;
 double nmin;
 if (( (int) pvectsz[ptsz->index_md] == ptsz->index_md_2halo)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_m_y_y_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_lens_lens_2h)
 || (((int) pvectsz[ptsz->index_md] == ptsz->index_md_cib_cib_2h)  && (pvectsz[ptsz->index_frequency_for_cib_profile] == pvectsz[ptsz->index_frequency_prime_for_cib_profile]) )
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_gal_gal_2h)){

     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[ptsz->index_z],ptsz)*nmin;
     double bmin_umin = bmin*integrand_patterson_test(log(ptsz->m_min_counter_terms),params)/pvectsz[ptsz->index_hmf]/pvectsz[ptsz->index_halo_bias];
     r += bmin_umin;
    }
else {
       // printf("counter terms at low M\n");
       double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
       double nmin_umin = nmin*integrand_patterson_test(log(ptsz->m_min_counter_terms),params)/pvectsz[ptsz->index_hmf];
       r += nmin_umin;
       }

       // printf("counter terms done\n");
                                   }
  }

if (( (int) pvectsz[ptsz->index_md] == ptsz->index_md_2halo)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_m_y_y_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_lens_lens_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_lensmag_lensmag_2h)
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_lens_lensmag_2h)
 || (((int) pvectsz[ptsz->index_md] == ptsz->index_md_cib_cib_2h)
      && (pvectsz[ptsz->index_frequency_for_cib_profile] == pvectsz[ptsz->index_frequency_prime_for_cib_profile]) )
 || ((int) pvectsz[ptsz->index_md] == ptsz->index_md_gal_gal_2h)){


 // if ((int) pvectsz[ptsz->index_md] == ptsz->index_md_lens_lens_2h){
 //   //printf("ok\n");
 //   r = r+ptsz->Omega_m_0*ptsz->Rho_crit_0*pow(pvecback[pba->index_bg_ang_distance]*pba->h,-2.)/pvectsz[ptsz->index_lensing_Sigma_crit];
 // }

 pvectsz[ptsz->index_integral_over_m] = r*r; //BB commented for debug
//pvectsz[ptsz->index_integral_over_m] = r/r; // BB debug
//printf("2halo_int = %.3e at z = %.3e\n",r,pvectsz[ptsz->index_z]); // BB debug
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
  class_open(process,"sz_auxiliary_files/SZ_thetas.txt", "r",ptsz->error_message);


  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

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
  fclose(process);
  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

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
  //printf("theta_bin_max:=%e\n",theta_bin_max);
  //printf("theta_bin_min:=%e\n",theta_bin_min);

  ///////////////////////////end read theta file

  //end read theta file for completeness
  // start read noise map for completeness
  //read skyfracs

  //double *skyfracs = NULL;


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
  class_open(process,"sz_auxiliary_files/SZ_skyfracs.txt", "r",ptsz->error_message);

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
  fclose(process);
  class_test(status != 0.,
             ptsz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  ptsz->nskyfracs = n_data;


  //end read skyfracs


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

  // sprintf(Filepath,
  //         "%s%s",
  //         // "%s%s%s",
  //         "cat ",
  //         // ptsz->path_to_class,
  //         "/sz_auxiliary_files/SZ_ylims.txt");
  //
  // process = popen(Filepath, "r");
  class_open(process,"sz_auxiliary_files/SZ_ylims.txt", "r",ptsz->error_message);

  //printf("ok\n");

  index_patches=0;
  int index_thetas = 0;
  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);
    ptsz->ylims[index_patches][index_thetas]=this_lnx;
    index_patches+=1;
    if (index_patches >= ptsz->nskyfracs){
      index_patches=0;
      index_thetas+=1;
    }

  }

  // status = pclose(process);
  fclose(process);
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
      fclose(process);
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

        // char Filepath[_ARGUMENT_LENGTH_MAX_];
        // sprintf(Filepath,
        //         "%s%s",
        //         // "%s%s%s",
        //         "cat ",
        //         // ptsz->path_to_class,
        //         "/sz_auxiliary_files/SO_files/SOSim_3freq_small_RMSTab_comp_test.txt");//MFMF_SOSim_3freq_small_RMSTab.txt");
        //
        // process = popen(Filepath, "r");
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
        fclose(process);
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


double * vrms2_var;
class_alloc(vrms2_var,
            sizeof(double *),
            ptsz->error_message);


class_alloc(ptsz->array_vrms2_at_z,sizeof(double *)*ptsz->n_arraySZ,ptsz->error_message);

int index_z;

for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
        {

            spectra_vrms2(pba,
                          ppm,
                          pnl,
                          ptsz,
                          exp(ptsz->array_redshift[index_z])-1.,
                          vrms2_var
                          );
          ptsz->array_vrms2_at_z[index_z] = log(*vrms2_var);

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
    //double x=exp(ln_x);
    double z = V->z;
    // M_halo = 1e16;
    // z = .5;
    // nc = HOD_mean_number_of_central_galaxies(z,M_halo,V->ptsz->M_min_HOD,V->ptsz->sigma_lnM_HOD,V->ptsz);
    // ns =  HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,V->ptsz->M_min_HOD,V->ptsz->alpha_s_HOD,V->ptsz->M1_prime_HOD,V->ptsz);
    // printf("z=%.3e Mhalo=%.3e  nc=%.3e ns=%.3e\n",z,M_halo,nc,ns);
    // exit(0);
    //exit(0);
    // if (V->ptsz->galaxy_sample ==1){ //unwise case, at the moment just nc=1, ns=0
    // nc = 1.;
    // ns = 0.;
    // }


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

      evaluate_HMF(lnM_halo,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];//*ptsz->Rho_crit_0/pba->h/pba->h;//pvectsz[ptsz->index_hmf];

      double z_asked = z;
      double  m_asked = M_halo;
      // this  also works:
      //double hmf = get_dndlnM_at_z_and_M(z_asked,m_asked,V->ptsz);



      V->pvectsz[V->ptsz->index_md] = V->ptsz->index_md_gal_gal_1h;
      double nc = HOD_mean_number_of_central_galaxies(z,M_halo,V->ptsz->M_min_HOD,V->ptsz->sigma_lnM_HOD,V->pvectsz,V->ptsz,V->pba);
      double ns = HOD_mean_number_of_satellite_galaxies(z,M_halo,nc,V->ptsz->M_min_HOD,V->ptsz->alpha_s_HOD,V->ptsz->M1_prime_HOD,V->ptsz,V->pba);


      double result = hmf*(ns+nc);
      //double result = (ns+nc)/log(10.);

      // printf("hmf = %.3e\n",hmf);
      // printf("ns = %.3e\n",ns);
      // printf("nc = %.3e\n",nc);

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
m_min = ptsz->m_min_counter_terms;//ptsz->M1SZ;
m_max = ptsz->m_max_counter_terms;//ptsz->M2SZ;
// m_min = ptsz->M1SZ;
// m_max = ptsz->M2SZ;
double * pvecback;
double * pvectsz;


 class_alloc(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);
   int i;
   for(i = 0; i<ptsz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),ptsz->error_message);




for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
        {
          double z = exp(ptsz->array_redshift[index_z])-1.;

          // if using unwise we set a specific lower bound
          // if (ptsz->galaxy_sample==1 &&  (ptsz->use_hod ==1)){
          // m_min = evaluate_unwise_m_min_cut(pvectsz[ptsz->index_z],ptsz->unwise_galaxy_sample_id);
          // if (m_min>=m_max)
          // m_min=m_max;
          //  }
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

           double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[ptsz->index_z],ptsz);
           double nmin_umin = nmin*integrand_mean_galaxy_number(log(m_min),params)/pvectsz[ptsz->index_hmf];
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

      evaluate_HMF(lnM_halo,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];//*ptsz->Rho_crit_0/pba->h/pba->h;//pvectsz[ptsz->index_hmf];

      double z_asked = z;
      double  m_asked = M_halo;
      // this  also works:
      //double hmf = get_dndlnM_at_z_and_M(z_asked,m_asked,V->ptsz);

      // printf("hmf = %.3e\n",hmf);
      // printf("ns = %.3e\n",ns);
      // printf("nc = %.3e\n",nc);
      double rho_crit_at_z = V->ptsz->Rho_crit_0;// V->pvectsz[V->ptsz->index_Rho_crit];
      double Omega_cb = (V->pba->Omega0_cdm + V->pba->Omega0_b);//*pow(1.+z,3.);
      double rho_cb = rho_crit_at_z*Omega_cb;


      double result = hmf*M_halo/rho_cb;

      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
      double b1 = V->pvectsz[V->ptsz->index_halo_bias];
      result *= b1;
      //double result = (ns+nc)/log(10.);



  return result;

}



int tabulate_hmf_counter_terms_b1min(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct tszspectrum * ptsz){

class_alloc(ptsz->array_hmf_counter_terms_b1min,sizeof(double *)*ptsz->n_z_hmf_counter_terms,ptsz->error_message);
// class_alloc(ptsz->array_redshift_hmf_counter_terms,sizeof(double *)*ptsz->n_z_hmf_counter_terms,ptsz->error_message);

int index_z;
double r;
double m_min,m_max;
// m_min = ptsz->M1SZ;
// m_max = ptsz->M2SZ;
m_min = ptsz->m_min_counter_terms;
m_max = ptsz->m_max_counter_terms;
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

          // ptsz->array_redshift_hmf_counter_terms[index_z] =
          //                             log(1.+z_min)
          //                             +index_z*(log(1.+z_max)-log(1.+z_min))
          //                             /(ptsz->n_z_hmf_counter_terms-1.); // log(1+z)

          double z = exp(ptsz->array_redshift_hmf_counter_terms[index_z])-1.;

          // if using unwise we set a specific lower bound
          // if (ptsz->galaxy_sample==1 &&  (ptsz->use_hod ==1)){
          // m_min = evaluate_unwise_m_min_cut(pvectsz[ptsz->index_z],ptsz->unwise_galaxy_sample_id);
          // if (m_min>=m_max)
          // m_min=m_max;
          //  }
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
          ptsz->array_hmf_counter_terms_b1min[index_z] = b1_min;

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

      evaluate_HMF(lnM_halo,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];//*ptsz->Rho_crit_0/pba->h/pba->h;//pvectsz[ptsz->index_hmf];

      double z_asked = z;
      double  m_asked = M_halo;
      // this  also works:
      //double hmf = get_dndlnM_at_z_and_M(z_asked,m_asked,V->ptsz);

      // printf("hmf = %.3e\n",hmf);
      // printf("ns = %.3e\n",ns);
      // printf("nc = %.3e\n",nc);
      double rho_crit_at_z = V->ptsz->Rho_crit_0;// V->pvectsz[V->ptsz->index_Rho_crit];
      double Omega_cb = (V->pba->Omega0_cdm + V->pba->Omega0_b);//*pow(1.+z,3.);
      double rho_cb = rho_crit_at_z*Omega_cb;


      double result = hmf*M_halo/rho_cb;

      evaluate_halo_bias_b2(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ptsz);
      double b2 = V->pvectsz[V->ptsz->index_halo_bias_b2];
      result *= b2;
      //double result = (ns+nc)/log(10.);



  return result;

}



int tabulate_hmf_counter_terms_b2min(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct tszspectrum * ptsz){

class_alloc(ptsz->array_hmf_counter_terms_b2min,sizeof(double *)*ptsz->n_z_hmf_counter_terms,ptsz->error_message);
// class_alloc(ptsz->array_redshift_hmf_counter_terms,sizeof(double *)*ptsz->n_z_hmf_counter_terms,ptsz->error_message);

int index_z;
double r;
double m_min,m_max;
m_min = ptsz->m_min_counter_terms;
m_max = ptsz->m_max_counter_terms;
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

          // ptsz->array_redshift_hmf_counter_terms[index_z] =
          //                             log(1.+z_min)
          //                             +index_z*(log(1.+z_max)-log(1.+z_min))
          //                             /(ptsz->n_z_hmf_counter_terms-1.); // log(1+z)

          double z = exp(ptsz->array_redshift_hmf_counter_terms[index_z])-1.;

          // if using unwise we set a specific lower bound
          // if (ptsz->galaxy_sample==1 &&  (ptsz->use_hod ==1)){
          // m_min = evaluate_unwise_m_min_cut(pvectsz[ptsz->index_z],ptsz->unwise_galaxy_sample_id);
          // if (m_min>=m_max)
          // m_min=m_max;
          //  }
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
          ptsz->array_hmf_counter_terms_b2min[index_z] = b2_min;

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

      evaluate_HMF(lnM_halo,V->pvecback,V->pvectsz,V->pba,V->pnl,V->ptsz);

      double hmf = V->pvectsz[V->ptsz->index_hmf];//*ptsz->Rho_crit_0/pba->h/pba->h;//pvectsz[ptsz->index_hmf];

      double z_asked = z;
      double  m_asked = M_halo;
      // this  also works:
      //double hmf = get_dndlnM_at_z_and_M(z_asked,m_asked,V->ptsz);

      // printf("hmf = %.3e\n",hmf);
      // printf("ns = %.3e\n",ns);
      // printf("nc = %.3e\n",nc);
      double rho_crit_at_z = V->ptsz->Rho_crit_0;// V->pvectsz[V->ptsz->index_Rho_crit];
      double Omega_cb = (V->pba->Omega0_cdm + V->pba->Omega0_b);//*pow(1.+z,3.);
      double rho_cb = rho_crit_at_z*Omega_cb;


      double result = hmf*M_halo/rho_cb;
      //double result = (ns+nc)/log(10.);



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
m_min = ptsz->m_min_counter_terms;
m_max = ptsz->m_max_counter_terms;
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

          ptsz->array_redshift_hmf_counter_terms[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_hmf_counter_terms-1.); // log(1+z)

          double z = exp(ptsz->array_redshift_hmf_counter_terms[index_z])-1.;

          // if using unwise we set a specific lower bound
          // if (ptsz->galaxy_sample==1 &&  (ptsz->use_hod ==1)){
          // m_min = evaluate_unwise_m_min_cut(pvectsz[ptsz->index_z],ptsz->unwise_galaxy_sample_id);
          // if (m_min>=m_max)
          // m_min=m_max;
          //  }
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
          double rho_crit_at_z = ptsz->Rho_crit_0;//pvectsz[ptsz->index_Rho_crit];
          double Omega_cb = (pba->Omega0_cdm + pba->Omega0_b);//*pow(1.+z,3.);
          double rho_cb = rho_crit_at_z*Omega_cb;
          double n_min = (1.-r)*rho_cb/m_min;
          ptsz->array_hmf_counter_terms_nmin[index_z] = n_min;

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


// tabulate L_nu^sat as a function of M (M_host) and z at frequency nu

int tabulate_L_sat_at_nu_and_nu_prime(struct background * pba,
                                      struct tszspectrum * ptsz){

if (
      ptsz->has_tSZ_cib_1h
    + ptsz->has_tSZ_cib_2h
    + ptsz->has_cib_cib_1h
    + ptsz->has_cib_cib_2h
    + ptsz->has_lens_cib_1h
    + ptsz->has_lens_cib_2h
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

// class_alloc(ptsz->array_L_sat_at_z_and_M_at_nu_prime,
//             sizeof(double *)*ptsz->n_z_L_sat*ptsz->n_m_L_sat,
//             ptsz->error_message);
//
//
// class_alloc(array_L_sat_at_z_and_M_at_nu_prime,
//             ptsz->n_z_L_sat*sizeof(double *),
//             ptsz->error_message);
//

for (index_l=0;
     index_l<ptsz->n_z_L_sat;
     index_l++)
{
  class_alloc(array_L_sat_at_z_and_M_at_nu[index_l],
              ptsz->n_m_L_sat*sizeof(double),
              ptsz->error_message);
  // class_alloc(array_L_sat_at_z_and_M_at_nu_prime[index_l],
  //             ptsz->n_m_L_sat*sizeof(double),
  //             ptsz->error_message);
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

      // V.nu = ptsz->nu_prime_cib_GHz;
      // params = &V;
      // double  L_sat_at_nu_prime = Integrate_using_Patterson_adaptive(lnMs_min, lnMs_max,
      //                                                          epsrel, epsabs,
      //                                                          integrand_patterson_L_sat,
      //                                                          params,ptsz->patterson_show_neval);


      // double dlnM = (lnMs_max - lnMs_min)/50.;
      //
      // double L_sat_at_nu = 0.;
      // double L_sat_at_nu_prime = 0.;
      // double L_gal_at_nu;
      // double L_gal_at_nu_prime;
      // double lnMs = lnMs_min;
      // double M_sub;
      // double M_host = exp(logM);
      // double dNdlnMs;
      //
      // double nu  = ptsz->nu_cib_GHz;
      // double nu_prime = ptsz->nu_prime_cib_GHz;
      //
      // while (lnMs<lnMs_max){
      // M_sub = exp(lnMs);
      // L_gal_at_nu = evaluate_galaxy_luminosity(z, M_sub, nu, ptsz);
      // L_gal_at_nu_prime = evaluate_galaxy_luminosity(z, M_sub, nu_prime, ptsz);
      // dNdlnMs = subhalo_hmf_dndlnMs(M_host,M_sub);
      // L_sat_at_nu += L_gal_at_nu*dNdlnMs;
      // L_sat_at_nu_prime += L_gal_at_nu_prime*dNdlnMs;
      // lnMs += dlnM;
      // }


      array_L_sat_at_z_and_M_at_nu[index_z][index_M] = log(1.+L_sat_at_nu);
      //array_L_sat_at_z_and_M_at_nu_prime[index_z][index_M] = log(1.+L_sat_at_nu_prime);
      //printf("%.3e %.3e %.3e %.3e %.3e %.3e\n",nu, nu_prime, z,logM,log(1.+L_sat_at_nu),log(1.+L_sat_at_nu_prime));

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
    //ptsz->array_L_sat_at_z_and_M_at_nu_prime[index_nu][index_z_M] = array_L_sat_at_z_and_M_at_nu_prime[index_z][index_M];
    index_z_M += 1;
  }
}

  free(array_L_sat_at_z_and_M_at_nu);
  //free(array_L_sat_at_z_and_M_at_nu_prime);

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

   // bounds array of radii for sigma computations:
   ptsz->logR1SZ = log(pow(3.*0.1*1e6/(4*_PI_*ptsz->Omega_m_0*ptsz->Rho_crit_0),1./3.));
   ptsz->logR2SZ = log(pow(3.*10.*1e17/(4*_PI_*ptsz->Omega_m_0*ptsz->Rho_crit_0),1./3.));


   // ptsz->logR1SZ = -5.684; // 0.0034Mpc/h, 1.8e4  solar mass
   // ptsz->logR2SZ = 4.; //default =4 , i.e., 54.9Mpc/h, 7.5e16 solar mass


  //Array of z
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM)+1e-7;
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  int index_z;

  double tstart, tstop;
  int index_l;
  double * sigma_var;
  double * dsigma_var;
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
private(tstart, tstop,index_R,index_z,sigma_var,dsigma_var) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

  class_alloc_parallel(sigma_var,
                       sizeof(double *),
                       ptsz->error_message);

  class_alloc_parallel(dsigma_var,
                       sizeof(double *),
                       ptsz->error_message);


#pragma omp for schedule (dynamic)
  for (index_R=0; index_R<ptsz->ndimSZ; index_R++)
  {
#pragma omp flush(abort)

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



      if (ptsz->HMF_prescription_NCDM == 2) //No-pres
        spectra_sigma_for_tSZ( pba,
                              ppm,
                              pnl,
                              exp(ptsz->array_radius[index_R]),
                              exp(ptsz->array_redshift[index_z])-1.,
                              sigma_var//&sigma_at_z_and_R
                              );
      else
        spectra_sigma_ncdm( pba,
                           // spectra_sigma_ncdm( pba,
                           ppm,
                           pnl,
                           //ptsz,
                           exp(ptsz->array_radius[index_R]),
                           exp(ptsz->array_redshift[index_z])-1.,
                           sigma_var//&sigma_at_z_and_R
                           );


      //ptsz->array_sigma_at_z_and_R[index_z_R] = log(*sigma_var);//sigma_at_z_and_R); //log(sigma)
      array_sigma_at_z_and_R[index_z][index_R] = log(*sigma_var);//sigma_at_z_and_R); //log(sigma)

      if (ptsz->HMF_prescription_NCDM == 2) //No-pres
        spectra_sigma_prime( pba,
                            ppm,
                            pnl,
                            exp(ptsz->array_radius[index_R]),
                            exp(ptsz->array_redshift[index_z])-1.,
                            dsigma_var//&dsigma2dR_at_z_and_R
                            );
      else
        spectra_sigma_ncdm_prime( pba,
                                 ppm,
                                 pnl,
                                 exp(ptsz->array_radius[index_R]),
                                 exp(ptsz->array_redshift[index_z])-1.,
                                 dsigma_var
                                 );


      array_dsigma2dR_at_z_and_R[index_z][index_R] = *dsigma_var;

      index_z_R += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over R's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(sigma_var);
    free(dsigma_var);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_R = 0;
for (index_R=0; index_R<ptsz->ndimSZ; index_R++)
{
  for (index_z=0; index_z<ptsz->n_arraySZ; index_z++)
  {
    ptsz->array_sigma_at_z_and_R[index_z_R] = array_sigma_at_z_and_R[index_z][index_R];
    ptsz->array_dsigma2dR_at_z_and_R[index_z_R]=array_dsigma2dR_at_z_and_R[index_z][index_R];
    index_z_R += 1;
  }
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
  + ptsz->has_lens_lensmag_1h
  + ptsz->has_lens_lensmag_2h
  + ptsz->has_tSZ_lensmag_1h
  + ptsz->has_tSZ_lensmag_2h
  + ptsz->has_lensmag_lensmag_1h
  + ptsz->has_lensmag_lensmag_2h == _FALSE_){
    if (ptsz->sz_verbose>=1)
    printf("-> Not tabulating Wz for lensing magnification\n");
    return 0;
  }

if (ptsz->sz_verbose>=1)
printf("-> Tabulating Wz for lensing magnification\n");
  //Array of z
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  int index_z;
  double ln1pz,z;

printf("-> start tabulating Wz for lensing magnification\n");
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

    //printf("-> Computing integral at z=%.3e\n",z);
    double result;
    redshift_int_lensmag(ptsz,pba,pvectsz,&result);
      // printf("-> 2 doing tabulating Wz for lensing magnification\n");
    if (result <= 0.)
      result = 1e-100;
    ptsz->array_W_lensmag[index_z] = log(result);
    ptsz->array_z_W_lensmag[index_z] = ln1pz;
    // printf("-> integral z = %.3e W = =%.3e\n",z,exp(ptsz->array_W_lensmag[index_z]));
}

// printf("-> end tabulating Wz for lensing magnification\n");
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




///Tabulate dndlnM
//as functions of z and M
int tabulate_dndlnM(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct tszspectrum * ptsz){

  //Array of z
  double z_min = r8_min(ptsz->z1SZ,ptsz->z1SZ_dndlnM);
  double z_max = r8_max(ptsz->z2SZ,ptsz->z2SZ_dndlnM);
  int index_z;

  double tstart, tstop;
  int index_l;
  double * dndlnM_var;
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
pba,ptsz,ppm,pnl,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,dndlnM_var,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

  class_alloc_parallel(dndlnM_var,
                       sizeof(double *),
                       ptsz->error_message);
  class_alloc_parallel(pvectsz,ptsz->tsz_size*sizeof(double),ptsz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<ptsz->n_z_dndlnM; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<ptsz->n_m_dndlnM; index_M++)
{
      ptsz->array_z_dndlnM[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(ptsz->n_z_dndlnM-1.); // log(1+z)

      ptsz->array_m_dndlnM[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(ptsz->n_m_dndlnM-1.); //log(R)

      //background quantities @ z:
      double z =   exp(ptsz->array_z_dndlnM[index_z])-1.;
      double logM =   ptsz->array_m_dndlnM[index_M];
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

      double omega = pvecback[pba->index_bg_Omega_m];
      pvectsz[ptsz->index_Delta_c]= Delta_c_of_Omega_m(omega);

      evaluate_HMF(logM,pvecback,pvectsz,pba,pnl,ptsz);


      *dndlnM_var = pvectsz[ptsz->index_hmf];//*ptsz->Rho_crit_0/pba->h/pba->h;//pvectsz[ptsz->index_hmf];



      array_dndlnM_at_z_and_M[index_z][index_M] = log(*dndlnM_var);

      index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (ptsz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(dndlnM_var);
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

  free(array_dndlnM_at_z_and_M);

return _SUCCESS_;
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



double get_dndlnM_at_z_and_M(double z_asked, double m_asked, struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
 return exp(pwl_interp_2d(ptsz->n_z_dndlnM,
                          ptsz->n_m_dndlnM,
                          ptsz->array_z_dndlnM,
                          ptsz->array_m_dndlnM,
                          ptsz->array_dndlnM_at_z_and_M,
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
//printf("z = %.3e  alpha = %.3e\n",z_asked,result);
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




double get_L_sat_at_z_and_M_at_nu(double z_asked,
                                  double m_asked,
                                  int index_nu,
                                  struct background * pba,
                                  struct tszspectrum * ptsz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
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
