# include "class_sz.h"
# include "class_sz_tools.h"
# include "class_sz_custom_profiles.h"
# include "class_sz_custom_bias.h"
# include "Patterson.h"
# include "r8lib.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "fft.h"

/////////////////////////////////CLASS_SZ-TOOLS//////////

int zbrent_sz_delta_to_delta_prime_nfw(double x1,
                                        double x2,
                                        double tol,
                                        double cvir,
                                        double cvir_prime,
                                        double delta,
                                        double fa,
                                        double fb,
                                        double * delta_prime,
                                       struct class_sz_structure * pclass_sz){
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
             pclass_sz->error_message,
             pclass_sz->error_message);

  class_call(
               dtod_prime_nfw(
                      b,
                      delta,
                      cvir,
                      cvir_prime,
                      &fb
                      ),
             pclass_sz->error_message,
             pclass_sz->error_message);


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
               pclass_sz->error_message,
               pclass_sz->error_message);
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
  struct class_sz_structure * pclass_sz
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
              pclass_sz->error_message);



  mTEST[0] = delta;


  class_call(
             dtod_prime_nfw(
                    mTEST[0],
                    delta,
                    cvir,
                    cvir_prime,
                    &lTEST
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message
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
                 pclass_sz->error_message,
                 pclass_sz->error_message
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
                 pclass_sz->error_message,
                 pclass_sz->error_message);

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
                      pclass_sz),
             pclass_sz->error_message,
             pclass_sz->error_message);





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
              struct class_sz_structure * pclass_sz
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
                    pclass_sz
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message);

  class_call(
             mVtomD(
                    b,
                    VAR1,
                    VAR2,
                    VAR3,
                    VAR4,
                    &fb,
                    pclass_sz
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message);


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
                      pclass_sz
                      ),
               pclass_sz->error_message,
               pclass_sz->error_message);
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
              struct class_sz_structure * pclass_sz,
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
  cvir_test = evaluate_cvir_of_mvir(mvir_test,z,pclass_sz,pba);
  rvir_test = evaluate_rvir_of_mvir(mvir_test,delc,rhoc,pclass_sz);



  class_call(
             mDtomV(
                    a,
                    mDEL,
                    rvir_test,
                    cvir_test,
                    delrho,
                    &fa,
                    pclass_sz
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message);

  mvir_test = exp(b);
  cvir_test = evaluate_cvir_of_mvir(mvir_test,z,pclass_sz,pba);
  rvir_test = evaluate_rvir_of_mvir(mvir_test,delc,rhoc,pclass_sz);


  class_call(
             mDtomV(
                    b,
                    mDEL,
                    rvir_test,
                    cvir_test,
                    delrho,
                    &fb,
                    pclass_sz
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message);


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
  cvir_test = evaluate_cvir_of_mvir(mvir_test,z,pclass_sz,pba);
  rvir_test = evaluate_rvir_of_mvir(mvir_test,delc,rhoc,pclass_sz);

    class_call(
               mDtomV(
                      b,
                      mDEL,
                      rvir_test,
                      cvir_test,
                      delrho,
                      &fb,
                      pclass_sz
                      ),
               pclass_sz->error_message,
               pclass_sz->error_message);
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
              struct class_sz_structure * pclass_sz,
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
                        pclass_sz,
                        pba,
                        pnl,
                        ppm
                        ),
             pclass_sz->error_message,
             pclass_sz->error_message);

  // knl_test = exp(b);


  class_call(
             pkl_to_knl(
                    b,
                    &fb,
                    z,
                    pclass_sz,
                    pba,
                    pnl,
                    ppm
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message);


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
                      pclass_sz,
                      pba,
                      pnl,
                      ppm
                      ),
               pclass_sz->error_message,
               pclass_sz->error_message);
  }

  printf("Max. num. of ite. exceeded in ZBRENT\n");

  return _FAILURE_;
}




//Root finding algorithm
//for finding z of chi
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
             chi_to_z(
                        a,
                        &fa,
                        chi,
                        pclass_sz,
                        pba
                        ),
             pclass_sz->error_message,
             pclass_sz->error_message);

  // knl_test = exp(b);


  class_call(
             chi_to_z(
                    b,
                    &fb,
                    chi,
                    pclass_sz,
                    pba
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message);


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
               chi_to_z(
                      b,
                      &fb,
                      chi,
                      pclass_sz,
                      pba
                      ),
               pclass_sz->error_message,
               pclass_sz->error_message);
  }

  printf("Max. num. of ite. exceeded in ZBRENT\n");

  return _FAILURE_;
}




//Root finding algorithm
//for the inverting ym relation
int zbrent_y_to_m(
              double x1,
              double x2,
              double tol,
              double fa,
              double fb,
              double * knl,
              double z,
              double y,
              // double rd,
              struct class_sz_structure * pclass_sz,
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
             y_to_m(
                        a,
                        &fa,
                        z,
                        y,
                        // rd,
                        pclass_sz,
                        pba,
                        pnl,
                        ppm
                        ),
             pclass_sz->error_message,
             pclass_sz->error_message);

  // knl_test = exp(b);


  class_call(
             y_to_m(
                    b,
                    &fb,
                    z,
                    y,
                    // rd,
                    pclass_sz,
                    pba,
                    pnl,
                    ppm
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message);


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
               y_to_m(
                      b,
                      &fb,
                      z,
                      y,
                      // rd,
                      pclass_sz,
                      pba,
                      pnl,
                      ppm
                      ),
               pclass_sz->error_message,
               pclass_sz->error_message);
  }

  printf("Max. num. of ite. exceeded in ZBRENT\n");

  return _FAILURE_;
}





//Root finding algorithm
//for the outer radius of gas density profile
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
             m_to_xout(
                        a,
                        &fa,
                        z,
                        m,
                        rd,
                        pclass_sz,
                        pba,
                        pnl,
                        ppm
                        ),
             pclass_sz->error_message,
             pclass_sz->error_message);

  // knl_test = exp(b);


  class_call(
             m_to_xout(
                    b,
                    &fb,
                    z,
                    m,
                    rd,
                    pclass_sz,
                    pba,
                    pnl,
                    ppm
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message);


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
               m_to_xout(
                      b,
                      &fb,
                      z,
                      m,
                      rd,
                      pclass_sz,
                      pba,
                      pnl,
                      ppm
                      ),
               pclass_sz->error_message,
               pclass_sz->error_message);
  }

  printf("Max. num. of ite. exceeded in ZBRENT\n");

  return _FAILURE_;
}



//This routine reads the tabulated
//pk(z,k) for n5k challenge
int load_n5k_pk_zk(
                      struct class_sz_structure * pclass_sz
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
  int index_k;
  int index_z;


  class_alloc(pclass_sz->n5k_pk_z,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->n5k_pk_k,sizeof(double *)*100,pclass_sz->error_message);



  n_data = 0;
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));



class_open(process,pclass_sz->full_path_to_n5k_z, "r",pclass_sz->error_message);

  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
                 "Error allocating memory to read the n5k files.\n");
      lnx = tmp;
    };


    /* Store */
    lnx[n_data]   = this_lnx;
    n_data++;
  }

  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  pclass_sz->n5k_pk_z_size = n_data;

  class_realloc(pclass_sz->n5k_pk_z,
                pclass_sz->n5k_pk_z,
                pclass_sz->n5k_pk_z_size*sizeof(double),
                pclass_sz->error_message);


  /** Store them */
  for (index_x=0; index_x<pclass_sz->n5k_pk_z_size; index_x++) {
    pclass_sz->n5k_pk_z[index_x] = lnx[index_x];
  };


  //Masses

  n_data = 0;
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));


  class_open(process,pclass_sz->full_path_to_n5k_k, "r",pclass_sz->error_message);

  // printf("-> %s\n",Filepath);

  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
                 "Error allocating memory to read the n5k file.\n");
      lnx = tmp;
    };


    /* Store */
    lnx[n_data]   = this_lnx;
    n_data++;
  }

  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  pclass_sz->n5k_pk_k_size = n_data;

  class_realloc(pclass_sz->n5k_pk_k,
                pclass_sz->n5k_pk_k,
                pclass_sz->n5k_pk_k_size*sizeof(double),
                pclass_sz->error_message);


  /** Store them */
  for (index_x=0; index_x<pclass_sz->n5k_pk_k_size; index_x++) {
    pclass_sz->n5k_pk_k[index_x] = log(lnx[index_x]);
  };


  /** Release the memory used locally */
  free(lnx);

  //Read pk

  class_alloc(pclass_sz->n5k_pk_pk,
              sizeof(double *)*pclass_sz->n5k_pk_z_size*pclass_sz->n5k_pk_k_size,
              pclass_sz->error_message);

  class_alloc(logC,
              pclass_sz->n5k_pk_z_size*sizeof(double *),
              pclass_sz->error_message);


  for (index_z=0;
       index_z<pclass_sz->n5k_pk_z_size;
       index_z++)
  {
    class_alloc(logC[index_z],
                pclass_sz->n5k_pk_k_size*sizeof(double),
                pclass_sz->error_message);
  }


  class_open(process,pclass_sz->full_path_to_n5k_pk_nl, "r",pclass_sz->error_message);


  int z =0;
  while (fgets(line, sizeof(line)-1, process) != NULL) {
    // printf("%s", line);
    // exit(0);
    int i=0;
    char *err, *p = line;
    double val;
    while (*p) {
      val = strtod(p, &err);
      logC[z][i] = log(val); //printf("%d-%.3e ",i,val);
      p = err + 1;
      i+=1;
    }
    // printf("\n %d \n",z);
    z+=1;
  }

  // printf("storing");
  int index = 0;
  for (index_z=0;
       index_z<pclass_sz->n5k_pk_z_size;
       index_z++){
    for (index_k=0;
         index_k<pclass_sz->n5k_pk_k_size;
         index_k++){

      pclass_sz->n5k_pk_pk[index] = logC[index_z][index_k];
      // printf("pk %.5e\n", logC[index_z][index_k]);//pclass_sz->n5k_pk_pk[index]);
      index += 1;
    }
  }

  status = fclose(process);


  for (index_z=0;
       index_z<pclass_sz->n5k_pk_z_size;
       index_z++){
         free(logC[index_z]);
       }
  free(logC);

  printf("n5k pk loaded with %d z and %d k\n",pclass_sz->n5k_pk_z_size,pclass_sz->n5k_pk_k_size);

  return _SUCCESS_;
}



//This routine reads the tabulated
//C-M relation Zhao2009,
//and stores the tabulated values.
//C(redshift[i],lnmass[j])

int read_Zhao_CM_init(
                      struct class_sz_structure * pclass_sz
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


  class_alloc(pclass_sz->CM_redshift,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->CM_logM,sizeof(double *)*100,pclass_sz->error_message);



  n_data = 0;
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));

  class_open(process,"class_sz_auxiliary_files/C-M_Zhao09/lnconcentration_vs_z_and_lnm-redshits.txt", "r",pclass_sz->error_message);

  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  pclass_sz->CM_redshift_size = n_data;

  class_realloc(pclass_sz->CM_redshift,
                pclass_sz->CM_redshift,
                pclass_sz->CM_redshift_size*sizeof(double),
                pclass_sz->error_message);


  /** Store them */
  for (index_x=0; index_x<pclass_sz->CM_redshift_size; index_x++) {
    pclass_sz->CM_redshift[index_x] = lnx[index_x];
  };


  //Masses

  n_data = 0;
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));


  class_open(process,"class_sz_auxiliary_files/C-M_Zhao09/lnconcentration_vs_z_and_lnm-masses.txt", "r",pclass_sz->error_message);



  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  pclass_sz->CM_logM_size = n_data;

  class_realloc(pclass_sz->CM_logM,
                pclass_sz->CM_logM,
                pclass_sz->CM_logM_size*sizeof(double),
                pclass_sz->error_message);


  /** Store them */
  for (index_x=0; index_x<pclass_sz->CM_logM_size; index_x++) {
    pclass_sz->CM_logM[index_x] = lnx[index_x];
  };


  /** Release the memory used locally */
  free(lnx);

  //Read concentration (lnC)

  class_alloc(pclass_sz->CM_logC,
              sizeof(double *)*pclass_sz->CM_redshift_size*pclass_sz->CM_logM_size,
              pclass_sz->error_message);

  class_alloc(logC,
              pclass_sz->CM_redshift_size*sizeof(double *),
              pclass_sz->error_message);


  for (index_redshift=0;
       index_redshift<pclass_sz->CM_redshift_size;
       index_redshift++)
  {
    class_alloc(logC[index_redshift],
                pclass_sz->CM_logM_size*sizeof(double),
                pclass_sz->error_message);
  }


  class_open(process,"class_sz_auxiliary_files/C-M_Zhao09/lnconcentration_vs_z_and_lnm.txt", "r",pclass_sz->error_message);



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
       index_redshift<pclass_sz->CM_redshift_size;
       index_redshift++){
    for (index_mass=0;
         index_mass<pclass_sz->CM_logM_size;
         index_mass++){

      pclass_sz->CM_logC[index] = logC[index_redshift][index_mass];
      index += 1;
    }
  }

  status = fclose(process);


  for (index_redshift=0;
       index_redshift<pclass_sz->CM_redshift_size;
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
                  struct class_sz_structure * pclass_sz
                  )
{


  double logz_asked = logz;
  double logM_asked = logM;

  if (logz<pclass_sz->CM_redshift[0])
    logz_asked = pclass_sz->CM_redshift[0];
  if (logz>pclass_sz->CM_redshift[pclass_sz->CM_redshift_size-1])
    logz_asked =  pclass_sz->CM_redshift[pclass_sz->CM_redshift_size-1];
  if (logM<pclass_sz->CM_logM[0])
    logM_asked = pclass_sz->CM_logM[0];
  if (logM>pclass_sz->CM_logM[pclass_sz->CM_logM_size-1])
    logM_asked =  pclass_sz->CM_logM[pclass_sz->CM_logM_size-1];

  *result = exp(pwl_interp_2d(
                              pclass_sz->CM_redshift_size,
                              pclass_sz->CM_logM_size,
                              pclass_sz->CM_redshift,
                              pclass_sz->CM_logM,
                              pclass_sz->CM_logC,
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
                  struct class_sz_structure * pclass_sz
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

int  CvirMvirKLYPIN(
                   double * result,
                   double logM ,
                   double z,
                   struct class_sz_structure * pclass_sz
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
              pclass_sz->error_message);
  class_alloc(M0z,
              1*sizeof(double),
              pclass_sz->error_message);

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
                            struct class_sz_structure * pclass_sz){

return pow(3.*mvir/(4*_PI_*delc*rhoc),1./3.);
                            }


double evaluate_cvir_of_mvir(double mvir,
                            double z,
                            struct class_sz_structure * pclass_sz,
                            struct background * pba){
double cvir;
//D08 c-m relation
if (pclass_sz->concentration_parameter==0){
cvir = 7.85*pow(mvir/2.e12,-0.081)*pow(1.+z,-0.71);
// cvir = 7.; // websky uses 7
}

//S00 c-m relation
else if (pclass_sz->concentration_parameter==1){
  cvir =10.*pow(mvir/3.42e12,-0.2)/(1.+z);
}

//K10 c-m relation
else if (pclass_sz->concentration_parameter==2){
   class_call(CvirMvirKLYPIN(&cvir,log(mvir),z,pclass_sz),
                   pclass_sz->error_message,
                   pclass_sz->error_message);
}


//SC14 c-m relation
// TBD m200c, r200c
else if (pclass_sz->concentration_parameter==3){
  printf("Warning: implementation of this concentration needs check.\n");
  exit(0);
   class_call(C200M200SC14(&cvir,
                           log(mvir),
                           z,
                           pclass_sz),
                   pclass_sz->error_message,
                   pclass_sz->error_message);
}


 //Z09 interpolated c-m relation
 else if (pclass_sz->concentration_parameter==4){
    printf("Warning: implementation of this concentration needs check.\n");
    exit(0);
    class_call(CvirMvirZHAO(&cvir,log(mvir),log(z),pclass_sz),
                    pclass_sz->error_message,
                    pclass_sz->error_message);
 }

// Dutton and Maccio 2014 (https://arxiv.org/pdf/1402.7073.pdf)
else if (pclass_sz->concentration_parameter==5){
  // here for virial mass in Msun/h:
  // see ea. 7 of 1402.7073
  double a =  0.537 + (1.025-0.537)*exp(-0.718*pow(z,1.08));
  double b = -0.097 + 0.024*z;
  double log10cvir = a + b*log10(mvir/1.e12);
  cvir = pow(10.,log10cvir);
}

// else if (pclass_sz->HMF==1 && pclass_sz->tau_profile == 1){
//
// }

else if (pclass_sz->concentration_parameter==6){ // Battacharya et al 2013

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

  // double sig = get_sigma_at_z_and_m(z,mvir,pclass_sz,pba);

double nu = sqrt(get_nu_at_z_and_m(z,mvir,pclass_sz,pba));

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
            struct class_sz_structure * pclass_sz
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
            struct class_sz_structure * pclass_sz,
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



struct Parameters_for_integrand_n5k_at_k{
  // struct nonlinear * pnl;
  // struct primordial * ppm;
  // struct perturbs * ppt;
  struct class_sz_structure * pclass_sz;
  // struct background * pba;
  // double * pvecback;
  // double * pvectsz;
  //double * llprime_grid;
  // double m;
  // double z;
  // double rd;
  double l;
};


double integrand_n5k_at_k(double lk, void *p){
  struct Parameters_for_integrand_n5k_at_k *V = ((struct Parameters_for_integrand_n5k_at_k *) p);
  double k = exp(lk);
  // double k = exp(pclass_sz->array_n5k_F1_k[index_k]);
  // int l = pclass_sz->array_n5k_F1_l[index_l];


  double chi_min = V->pclass_sz->chi_min_n5k_samp_fftw;//1e0;//pclass_sz->l_min_samp_fftw; //precision parameter
  double chi_max = V->pclass_sz->chi_max_n5k_samp_fftw;//7e3;//pclass_sz->l_max_samp_fftw; //precision parameter

  const int N = V->pclass_sz->N_samp_fftw; //precision parameter
  int ichi;
  double chi[N], Pchi[N];

  // printf("k = %.5e, chi_min = %.5e, chi_max = %.5e, N = %d\n", k, chi_min, chi_max, N);


  for (ichi=0; ichi<N; ichi++){
    chi[ichi] =  exp(log(chi_min)+ichi/(N-1.)*(log(chi_max)-log(chi_min)));
    double zchi = get_n5k_z_of_chi(chi[ichi],V->pclass_sz);
    Pchi[ichi] = sqrt(get_n5k_pk_at_z_and_k(zchi,k,V->pclass_sz))*get_n5k_cl_K1_at_chi(chi[ichi],V->pclass_sz);
    // if (ichi==101)
    // printf("ichi = %d Pchi = %.3e chi = %.3e zchi = %.5e K1 = %.5e\n",ichi,Pchi[ichi],chi[ichi],zchi,get_n5k_cl_K1_at_chi(chi[ichi],V->pclass_sz));
  }

  double chit[N], Pchit[N];
  /* Compute the function
  *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
  * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
  * in this notation.  The input k-values must be logarithmically spaced.  The
  * resulting xi_l^m(r) will be evaluated at the dual r-values
  *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
  fftlog_ComputeXiLMsloz(V->l, 0, N, chi,  Pchi, chit, Pchit,V->pclass_sz);
  double F1 = 2.*_PI_*_PI_*pwl_value_1d(N,chit,Pchit,k);
  fftlog_ComputeXiLMsloz(V->l, 0, N, chi,  Pchi, chit, Pchit,V->pclass_sz);
  double F2 = 2.*_PI_*_PI_*pwl_value_1d(N,chit,Pchit,k);
  double intk = F1*F2*k*k;
  intk *= k; // integrate in logk
  return intk;
  // double dlk = (log(k_max)-log(k_min))/(pclass_sz->n_k_n5k-1.);
  // sumk +=  intk*k*dlk;
}


struct Parameters_for_integrand_m_to_xout{
  struct nonlinear * pnl;
  struct primordial * ppm;
  // struct perturbs * ppt;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  // double * pvecback;
  // double * pvectsz;
  //double * llprime_grid;
  double m;
  double z;
  double rd;
  double c;
};


double integrand_m_to_xout(double x, void *p){

struct Parameters_for_integrand_m_to_xout *V = ((struct Parameters_for_integrand_m_to_xout *) p);

double r = 0.;
double cd = 1.;
// double rd = pow(3.*m/(4.*_PI_*200.*rho_crit),1./3.);
double rs = V->rd/cd;
r = 4.*_PI_*pow(rs,3.)*get_gas_profile_at_x_M_z_b16_200c(x,
                                               V->m,
                                               V->z,
                                               V->pclass_sz->c_B16, // TBC
                                               V->pclass_sz->A_rho0,
                                               V->pclass_sz->A_alpha,
                                               V->pclass_sz->A_beta,
                                               V->pclass_sz->alpha_m_rho0,
                                               V->pclass_sz->alpha_m_alpha,
                                               V->pclass_sz->alpha_m_beta,
                                               V->pclass_sz->alpha_z_rho0,
                                               V->pclass_sz->alpha_z_alpha,
                                               V->pclass_sz->alpha_z_beta,
                                               // break model param
					                                     V->pclass_sz->mcut,
					                                     V->pclass_sz->alphap_m_rho0,
                                               V->pclass_sz->alphap_m_alpha,
                                               V->pclass_sz->alphap_m_beta,
					                                     V->pclass_sz->alpha_c_rho0,
                                               V->pclass_sz->alpha_c_alpha,
                                               V->pclass_sz->alpha_c_beta,
                                               // end break model param
                                               V->pclass_sz->gamma_B16,
                                               V->pclass_sz->xc_B16,
                                               V->pba,
                                               V->pclass_sz)*pow(x,2);


return r;
}

//Routine used for
//finding the non-linear scale
int y_to_m(
            double xout,
            double * mRES,
            double z,
            double y,
            // double rd,
            struct class_sz_structure * pclass_sz,
            struct background * pba,
            struct nonlinear * pnl,
            struct primordial * ppm
            )
{

  // struct Parameters_for_integrand_y_to_m V;
  // V.pnl = pnl;
  // V.ppm = ppm;
  // V.pclass_sz = pclass_sz;
  // V.pba = pba;
  // V.m = m;
  // V.z = z;
  // // V.rd = rd;
  // // V.c = 0.; // TBC!
  // // V.pvectsz = Pvectsz;
  // // V.pvecback = Pvecback;
  //
  // void * params = &V;
  //
  //
  // double epsrel= pclass_sz->m_to_xout_epsrel;
  // double epsabs= pclass_sz->m_to_xout_epsabs;
  // int show_neval = pclass_sz->patterson_show_neval;
  // //integral of density profile.
  // double m_profile = Integrate_using_Patterson_adaptive(1e-5, xout,
  //                                                       epsrel, epsabs,
  //                                                       integrand_m_to_xout,
  //                                                       params,show_neval);
  //

  *mRES = get_y_at_m_and_z(xout,z,pclass_sz,pba) - y;

  return _SUCCESS_;
}



//Routine used for
//finding z of chi
int chi_to_z(
            double xout,
            double * mRES,
            double chi,
            struct class_sz_structure * pclass_sz,
            struct background * pba
            )
{


  double z = xout;
  double chitest;
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

  chitest = pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h;
  *mRES = chitest - chi;

  free(pvecback);

  return _SUCCESS_;
}








//Routine used for
//finding the non-linear scale
int m_to_xout(
            double xout,
            double * mRES,
            double z,
            double m,
            double rd,
            struct class_sz_structure * pclass_sz,
            struct background * pba,
            struct nonlinear * pnl,
            struct primordial * ppm
            )
{

  struct Parameters_for_integrand_m_to_xout V;
  V.pnl = pnl;
  V.ppm = ppm;
  V.pclass_sz = pclass_sz;
  V.pba = pba;
  V.m = m;
  V.z = z;
  V.rd = rd;
  V.c = pclass_sz->c_B16; // TBC!
  // V.pvectsz = Pvectsz;
  // V.pvecback = Pvecback;

  void * params = &V;


  double epsrel= pclass_sz->m_to_xout_epsrel;
  double epsabs= pclass_sz->m_to_xout_epsabs;
  int show_neval = pclass_sz->patterson_show_neval;
  //integral of density profile.
  double m_profile = Integrate_using_Patterson_adaptive(1e-5, xout,
                                                        epsrel, epsabs,
                                                        integrand_m_to_xout,
                                                        params,show_neval);


  *mRES =m_profile-pclass_sz->f_b_gas*m;

  return _SUCCESS_;
}




int tabulate_y_to_m(struct background * pba,
                   struct nonlinear * pnl,
                   struct primordial * ppm,
                   struct class_sz_structure * pclass_sz){

if (pclass_sz->sz_verbose > 0)
 printf("->SZ_counts tabulating y to m.\n");


double r;
double y_min,y_max;
y_min = exp(pclass_sz->lnymin); // for the mass integral
y_max = exp(pclass_sz->lnymax); // for the mass integral
int index_y;
for (index_y=0; index_y<pclass_sz->n_y_y_to_m; index_y++)
        {

          pclass_sz->array_y_to_m_y[index_y] =
                                      log(y_min)
                                      +index_y*(log(y_max)-log(y_min))
                                      /(pclass_sz->n_y_y_to_m-1.); // log(nu)
        }

double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;
int index_z;
for (index_z=0; index_z<pclass_sz->n_z_y_to_m; index_z++)
        {

          pclass_sz->array_y_to_m_redshift[index_z] =
                                                  log(1.+z_min)
                                                  +index_z*(log(1.+z_max)-log(1.+z_min))
                                                  /(pclass_sz->n_z_y_to_m-1.); // log(1+z)
        }



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
pba,pclass_sz,ppm,pnl,z_min,z_max,y_min,y_max)\
private(tstart, tstop,index_z,index_y,r) \
num_threads(number_of_threads)
{

  #pragma omp for collapse(2)
  for (index_z=0; index_z<pclass_sz->n_z_y_to_m; index_z++)
  {
    for (index_y=0; index_y<pclass_sz->n_y_y_to_m; index_y++)
    {

  #ifdef _OPENMP
    tstart = omp_get_wtime();
  #endif

// double xout_var; // in multiples of 200c
double z = exp(pclass_sz->array_y_to_m_redshift[index_z])-1.;;
double y = exp(pclass_sz->array_y_to_m_y[index_y]);

int index_z_y = index_y * pclass_sz->n_z_y_to_m + index_z;

solve_y_to_m(&r,
             z,
             y,
             pclass_sz,
             pba,
             pnl,
             ppm);


if (isinf(r)){
  printf("z = %.5e y=%.5e r = %.5e\n",z,y,r);
  exit(0);
}
  pclass_sz->array_y_to_m_at_z_y[index_z_y] = r;
}
}

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z m's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif
}
if (abort == _TRUE_) return _FAILURE_;
}

double get_z_of_chi(
                    double chi,
                    struct class_sz_structure * pclass_sz,
                    struct background * pba){

double z; 
solve_chi_to_z(&z,chi,pclass_sz,pba);
return z;

  }


int tabulate_m_to_xout(struct background * pba,
                       struct nonlinear * pnl,
                       struct primordial * ppm,
                       struct class_sz_structure * pclass_sz){

class_alloc(pclass_sz->array_m_to_xout_redshift,sizeof(double *)*pclass_sz->n_z_m_to_xout,pclass_sz->error_message);
class_alloc(pclass_sz->array_m_to_xout_mass,sizeof(double *)*pclass_sz->n_mass_m_to_xout,pclass_sz->error_message);
class_alloc(pclass_sz->array_m_to_xout_at_z_m,sizeof(double *)*pclass_sz->n_z_m_to_xout*pclass_sz->n_mass_m_to_xout,pclass_sz->error_message);


double r;
double m_min,m_max;
m_min = pclass_sz->M1SZ; // for the mass integral
m_max = pclass_sz->M2SZ; // for the mass integral
int index_m;
for (index_m=0; index_m<pclass_sz->n_mass_m_to_xout; index_m++)
        {

          pclass_sz->array_m_to_xout_mass[index_m] =
                                      log(m_min)
                                      +index_m*(log(m_max)-log(m_min))
                                      /(pclass_sz->n_mass_m_to_xout-1.); // log(nu)
        }

double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;
int index_z;
for (index_z=0; index_z<pclass_sz->n_z_m_to_xout; index_z++)
        {

          pclass_sz->array_m_to_xout_redshift[index_z] =
                                                  log(1.+z_min)
                                                  +index_z*(log(1.+z_max)-log(1.+z_min))
                                                  /(pclass_sz->n_z_m_to_xout-1.); // log(1+z)
        }



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
pba,pclass_sz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_m,r) \
num_threads(number_of_threads)
{

  #pragma omp for collapse(2)
  for (index_z=0; index_z<pclass_sz->n_z_m_to_xout; index_z++)
  {
    for (index_m=0; index_m<pclass_sz->n_mass_m_to_xout; index_m++)
    {

  #ifdef _OPENMP
    tstart = omp_get_wtime();
  #endif

double xout_var; // in multiples of 200c
double z = exp(pclass_sz->array_m_to_xout_redshift[index_z])-1.;;
double m = exp(pclass_sz->array_m_to_xout_mass[index_m]);

int index_z_m = index_m * pclass_sz->n_z_m_to_xout + index_z;

solve_m_to_xout(&r,
                 z,
                 m,
                 pclass_sz,
                 pba,
                 pnl,
                 ppm);

// printf("z = %.5e m=%.5e xout = %.5e\n",z,m,r);
  pclass_sz->array_m_to_xout_at_z_m[index_z_m] = r;
}
}

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z m's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif
}
if (abort == _TRUE_) return _FAILURE_;
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
            struct class_sz_structure * pclass_sz
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
               struct class_sz_structure * pclass_sz,
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
                        pclass_sz,
                        pba),
                pclass_sz->error_message,
                pclass_sz->error_message);

//then go from mvir to mdel_prime
double rvir = evaluate_rvir_of_mvir(mvir,delc,rhoc,pclass_sz);
double cvir = evaluate_cvir_of_mvir(mvir,z,pclass_sz,pba);

class_call(mVIR_to_mDEL(mvir,
                     rvir,
                     cvir,
                     delrho_prime,
                     mDELprime,
                     pclass_sz),
                pclass_sz->error_message,
                pclass_sz->error_message);

return _SUCCESS_;
             }




//Routine used for
//the non linear scale
int solve_pkl_to_knl(
              double * result,
              double z,
              struct class_sz_structure * pclass_sz,
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
              pclass_sz->error_message);



  mTEST[0] = 1.;

 // printf("res 0 ini : %.3e\n",lTEST);
  class_call(
             pkl_to_knl(
                    mTEST[0],
                    &lTEST,
                    z,
                    pclass_sz,
                    pba,
                    pnl,
                    ppm
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message
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
                        pclass_sz,
                        pba,
                        pnl,
                        ppm
                        ),
                 pclass_sz->error_message,
                 pclass_sz->error_message
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
                        pclass_sz,
                        pba,
                        pnl,
                        ppm
                        ),
                 pclass_sz->error_message,
                 pclass_sz->error_message);

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
                               pclass_sz,
                               pba,
                               pnl,
                               ppm
                               ),
             pclass_sz->error_message,
             pclass_sz->error_message);

  mDEL = logMDEL;
  *result = mDEL;


  free(mTEST);


  return _SUCCESS_;
}



//Routine used for
//the cut-off radius of the gnfw profile
int solve_m_to_xout(
                    double * result,
                    double z,
                    double m,
                    struct class_sz_structure * pclass_sz,
                    struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm
                    )
{

  // printf("z = %.5e m=%.5e xout = %.5e\n",z,m,r);
  // printf("z = %.5e m=%.5e\n",z,m);


/// get rhoc and rd

double * pvecback;
double * pvectsz;

double tau;
int first_index_back = 0;

class_alloc(pvecback,
            pba->bg_size*sizeof(double),
            pclass_sz->error_message);

class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
 int i;
 for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

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




// pvectsz[pclass_sz->index_z] = z;
pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);

double rho_crit = pvectsz[pclass_sz->index_Rho_crit];
double delta = 200.;//*pvecback[pba->index_bg_Omega_m];
double c_delta = get_c200c_at_m_and_z(m,z,pba,pclass_sz);
double rd = pow(3.*m/(4.*_PI_*delta*rho_crit),1./3.); //in units of h^-1 Mpc


free(pvecback);
free(pvectsz);
/////



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



  // int  i;
  int iMAX = 50;

  double * mTEST;
  class_alloc(mTEST,
              iMAX*sizeof( double ),
              pclass_sz->error_message);



  mTEST[0] = 1.;

 // printf("res 0 ini : %.3e\n",lTEST);
  class_call(
             m_to_xout(
                    mTEST[0],
                    &lTEST,
                    z,
                    m,
                    rd,
                    pclass_sz,
                    pba,
                    pnl,
                    ppm
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message
             );
 // printf("res 0 : %.3e\n",lTEST);
 //exit(0);

  if (lTEST <= 0.) {
    for (i=1;i<iMAX;i++ ) {

      mTEST[i] = 2.*mTEST[i-1];

      class_call(
                 m_to_xout(
                        mTEST[i],
                        &lTEST,
                        z,
                        m,
                        rd,
                        pclass_sz,
                        pba,
                        pnl,
                        ppm
                        ),
                 pclass_sz->error_message,
                 pclass_sz->error_message
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
                 m_to_xout(
                        mTEST[i],
                        &lTEST,
                        z,
                        m,
                        rd,
                        pclass_sz,
                        pba,
                        pnl,
                        ppm
                        ),
                 pclass_sz->error_message,
                 pclass_sz->error_message);

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

  class_call(zbrent_m_to_xout(
                               mLO,
                               mUP,
                               1.e-4,
                               fa,
                               fb,
                               &logMDEL,
                               z,
                               m,
                               rd,
                               pclass_sz,
                               pba,
                               pnl,
                               ppm
                               ),
             pclass_sz->error_message,
             pclass_sz->error_message);

  mDEL = logMDEL;
  *result = mDEL;

  free(mTEST);
  return _SUCCESS_;
}




//Routine used for
//the invert ym relation
int solve_y_to_m(
                    double * result,
                    double z,
                    double y,
                    struct class_sz_structure * pclass_sz,
                    struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm
                    )
{

  // printf("z = %.5e m=%.5e xout = %.5e\n",z,m,r);
  // printf("z = %.5e m=%.5e\n",z,m);


// /// get rhoc and rd
//
// double * pvecback;
// double * pvectsz;
//
// double tau;
// int first_index_back = 0;
//
// class_alloc(pvecback,
//             pba->bg_size*sizeof(double),
//             pclass_sz->error_message);
//
// class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
 int i;
//  for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;
//
// class_call(background_tau_of_z(pba,z,&tau),
//            pba->error_message,
//            pba->error_message);
//
// class_call(background_at_tau(pba,
//                              tau,
//                              pba->long_info,
//                              pba->inter_normal,
//                              &first_index_back,
//                              pvecback),
//            pba->error_message,
//            pba->error_message);
//
//
//
//
// // pvectsz[pclass_sz->index_z] = z;
// pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
//                                 *pow(_Mpc_over_m_,1)
//                                 *pow(_c_,2)
//                                 *pvecback[pba->index_bg_rho_crit]
//                                 /pow(pba->h,2);
//
// double rho_crit = pvectsz[pclass_sz->index_Rho_crit];
// double delta = 200.;//*pvecback[pba->index_bg_Omega_m];
// double c_delta = get_c200c_at_m_and_z(m,z,pba,pclass_sz);
// double rd = pow(3.*m/(4.*_PI_*delta*rho_crit),1./3.); //in units of h^-1 Mpc
//
//
// free(pvecback);
// free(pvectsz);
/////



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



  // int  i;
  int iMAX = 50;

  double * mTEST;
  class_alloc(mTEST,
              iMAX*sizeof( double ),
              pclass_sz->error_message);



  mTEST[0] = 1.e11;

 // printf("res 0 ini : %.3e\n",y);
  class_call(
             y_to_m(
                    mTEST[0],
                    &lTEST,
                    z,
                    y,
                    // rd,
                    pclass_sz,
                    pba,
                    pnl,
                    ppm
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message
             );
 // printf("res 0 : %.3e\n",lTEST);
 // exit(0);

  if (lTEST <= 0.) {
    for (i=1;i<iMAX;i++ ) {

      mTEST[i] = 2.*mTEST[i-1];

      class_call(
                 y_to_m(
                        mTEST[i],
                        &lTEST,
                        z,
                        y,
                        // rd,
                        pclass_sz,
                        pba,
                        pnl,
                        ppm
                        ),
                 pclass_sz->error_message,
                 pclass_sz->error_message
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
                 y_to_m(
                        mTEST[i],
                        &lTEST,
                        z,
                        y,
                        // rd,
                        pclass_sz,
                        pba,
                        pnl,
                        ppm
                        ),
                 pclass_sz->error_message,
                 pclass_sz->error_message);

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

  class_call(zbrent_y_to_m(
                               mLO,
                               mUP,
                               1.e-4,
                               fa,
                               fb,
                               &logMDEL,
                               z,
                               y,
                               // rd,
                               pclass_sz,
                               pba,
                               pnl,
                               ppm
                               ),
             pclass_sz->error_message,
             pclass_sz->error_message);

  mDEL = logMDEL;
  *result = mDEL;

  free(mTEST);
  return _SUCCESS_;
}





//Routine used for
//the invert chi(z)relation
int solve_chi_to_z(
                    double * result,
                    double chi,
                    struct class_sz_structure * pclass_sz,
                    struct background * pba
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
              pclass_sz->error_message);



  mTEST[0] = 1.;
 
 // printf("res 0 ini : %.3e\n",y);
  class_call(
             chi_to_z(
                    mTEST[0],
                    &lTEST,
                    chi,
                    pclass_sz,
                    pba
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message
             );
 // printf("res 0 : %.3e\n",lTEST);
 // exit(0);

  if (lTEST <= 0.) {
    for (i=1;i<iMAX;i++ ) {

      mTEST[i] = 2.*mTEST[i-1];

      class_call(
                 chi_to_z(
                        mTEST[i],
                        &lTEST,
                        chi,
                        pclass_sz,
                        pba
                        ),
                 pclass_sz->error_message,
                 pclass_sz->error_message
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
                 chi_to_z(
                        mTEST[i],
                        &lTEST,
                        chi,
                        pclass_sz,
                        pba
                        ),
                 pclass_sz->error_message,
                 pclass_sz->error_message);

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

  class_call(zbrent_chi_to_z(
                               mLO,
                               mUP,
                               1.e-4,
                               fa,
                               fb,
                               &logMDEL,
                               chi,
                               pclass_sz,
                               pba
                               ),
             pclass_sz->error_message,
             pclass_sz->error_message);

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
              struct class_sz_structure * pclass_sz
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
              pclass_sz->error_message);



  mTEST[0] = mVIR;


  class_call(
             mVtomD(
                    log(mTEST[0]),
                    mVIR,
                    rvir,
                    c,
                    delrho,
                    &lTEST,
                    pclass_sz
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message
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
                        pclass_sz
                        ),
                 pclass_sz->error_message,
                 pclass_sz->error_message
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
                        pclass_sz
                        ),
                 pclass_sz->error_message,
                 pclass_sz->error_message);

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
                       pclass_sz
                       ),
             pclass_sz->error_message,
             pclass_sz->error_message);

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
              struct class_sz_structure * pclass_sz,
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
              pclass_sz->error_message);


  // var[0] = mDEL;
  // // var[1] = rs;
  // // var[2] = c;
  // var[1] = delrho;
  double mvir_test,cvir_test,rvir_test;

  mTEST[0] = mDEL;

  mvir_test = mTEST[0];
  cvir_test = evaluate_cvir_of_mvir(mvir_test,z,pclass_sz,pba);
  rvir_test = evaluate_rvir_of_mvir(mvir_test,delc,rhoc,pclass_sz);




  class_call(
             mDtomV(
                    log(mTEST[0]),
                    mDEL,
                    rvir_test,
                    cvir_test,
                    delrho,
                    &lTEST,
                    pclass_sz
                    ),
             pclass_sz->error_message,
             pclass_sz->error_message
             );

  if (lTEST <= 0.) {
    for (i=1;i<iMAX;i++ ) {

      mTEST[i] = 2.*mTEST[i-1];

      mvir_test = mTEST[i];
      cvir_test = evaluate_cvir_of_mvir(mvir_test,z,pclass_sz,pba);
      rvir_test = evaluate_rvir_of_mvir(mvir_test,delc,rhoc,pclass_sz);


      class_call(
                 mDtomV(
                        log(mTEST[i]),
                        mDEL,
                        rvir_test,
                        cvir_test,
                        delrho,
                        &lTEST,
                        pclass_sz
                        ),
                 pclass_sz->error_message,
                 pclass_sz->error_message
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
      cvir_test = evaluate_cvir_of_mvir(mvir_test,z,pclass_sz,pba);
      rvir_test = evaluate_rvir_of_mvir(mvir_test,delc,rhoc,pclass_sz);

      class_call(
                 mDtomV(
                        log(mTEST[i]),
                        mDEL,
                        rvir_test,
                        cvir_test,
                        delrho,
                        &lTEST,
                        pclass_sz
                        ),
                 pclass_sz->error_message,
                 pclass_sz->error_message);

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
                       pclass_sz,
                       pba
                       ),
             pclass_sz->error_message,
             pclass_sz->error_message);

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
  struct class_sz_structure * pclass_sz;
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


  double Theta_s = sqrt(V->pclass_sz->Omega_survey/_PI_); // see below Eq. 45 of Takada and Spergel 2013
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
                   struct class_sz_structure * pclass_sz,
                   double z,
                   double * sigma2_hsv
                   ) {

double k_min = exp(pclass_sz->ln_k_for_tSZ[0]);
double k_max = exp(pclass_sz->ln_k_for_tSZ[pclass_sz->ln_k_size_for_tSZ-1]);


struct Parameters_for_integrand_sigma2_hsv V;
  V.pnl = pnl;
  V.ppm = ppm;
  V.pclass_sz = pclass_sz;
  V.pba = pba;
  V.z = z;

  void * params = &V;
  double r; //result of the integral

  double epsrel = 1e-6;
  double epsabs = 1e-30;
  //int show_neval = pclass_sz->patterson_show_neval;

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
  // // double epsrel=pclass_sz->mass_epsrel;
  // // double epsabs=pclass_sz->mass_epsabs;
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




struct Parameters_for_integrand_gallens_sources{
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
};



double integrand_gallens_sources(double ln1pzs, void *p){

  struct Parameters_for_integrand_gallens_sources *V = ((struct Parameters_for_integrand_gallens_sources *) p);

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
  double Chi_at_z = sqrt(V->pvectsz[V->pclass_sz->index_chi2])/V->pba->h;  //'Chi' comoving distance in Mpc



  free(pvecback);



  W = (Chi_at_zs-Chi_at_z)/Chi_at_zs;

  double dndzs = 0.;

/////////////////////////////////
  double z_asked  = zs;
  double phig = 0.;

phig = get_source_galaxy_number_counts(z_asked,V->pclass_sz);
 dndzs = phig;
////////////////////////////////

  integrand = dndzs*W;

  integrand *= (1.+zs);

  // printf("-> integrand z = %.3e phig = =%.3e\n",z_asked,phig);

  return integrand;

}

int redshift_int_gallens_sources(
                  struct class_sz_structure * pclass_sz,
                  struct background * pba,
                  double * pvectsz,
                  double * result
                   ) {

double z =  pvectsz[pclass_sz->index_z];

double zs_min = z;
double zs_max = pclass_sz->z2SZ;


struct Parameters_for_integrand_gallens_sources V;
  V.pvectsz = pvectsz;
  V.pclass_sz = pclass_sz;
  V.pba = pba;


  void * params = &V;
  double r; //result of the integral

  double epsrel = 1e-6;
  double epsabs = 1e-30;
  //int show_neval = pclass_sz->patterson_show_neval;

  r=Integrate_using_Patterson_adaptive(log(1.+zs_min),
                                        log(1.+zs_max),
                                        epsrel, epsabs,
                                        integrand_gallens_sources,
                                        params,0);


  *result = r;
                     }




struct Parameters_for_integrand_lensmag{
  struct class_sz_structure * pclass_sz;
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
  double Chi_at_z = sqrt(V->pvectsz[V->pclass_sz->index_chi2])/V->pba->h;  //'Chi' comoving distance in Mpc



  free(pvecback);



  W = (Chi_at_zs-Chi_at_z)/Chi_at_zs;

  double dndzs = 0.;

/////////////////////////////////
  double z_asked  = zs;
  double phig = 0.;
//unwise: use Cosmos cross-match dndz
if (V->pclass_sz->galaxy_sample==1){
if(z_asked<V->pclass_sz->normalized_cosmos_dndz_z[0])
   phig = 1e-100;
else if (z_asked>V->pclass_sz->normalized_cosmos_dndz_z[V->pclass_sz->normalized_cosmos_dndz_size-1])
   phig = 1e-100;
else  phig =  pwl_value_1d(V->pclass_sz->normalized_cosmos_dndz_size,
                             V->pclass_sz->normalized_cosmos_dndz_z,
                             V->pclass_sz->normalized_cosmos_dndz_phig,
                             z_asked);
// printf("integrand ok phig\n");
}
else{
  // if(z_asked<V->pclass_sz->normalized_dndz_z[0])
  //    phig = 1e-100;
  // else if (z_asked>V->pclass_sz->normalized_dndz_z[V->pclass_sz->normalized_dndz_size-1])
  //    phig = 1e-100;
  // else  phig =  pwl_value_1d(V->pclass_sz->normalized_dndz_size,
  //                              V->pclass_sz->normalized_dndz_z,
  //                              V->pclass_sz->normalized_dndz_phig,
  //                              z_asked);

phig = get_galaxy_number_counts(z_asked,V->pclass_sz);
}
 dndzs = phig;
////////////////////////////////

  integrand = dndzs*W;

  integrand *= (1.+zs);

  // printf("-> integrand z = %.3e phig = =%.3e\n",z_asked,phig);

  return integrand;

}


struct Parameters_for_integrand_nlensmag{
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
  int index_g;
};



double integrand_nlensmag(double ln1pzs, void *p){

  struct Parameters_for_integrand_nlensmag *V = ((struct Parameters_for_integrand_nlensmag *) p);

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
  double Chi_at_z = sqrt(V->pvectsz[V->pclass_sz->index_chi2])/V->pba->h;  //'Chi' comoving distance in Mpc



  free(pvecback);



  W = (Chi_at_zs-Chi_at_z)/Chi_at_zs;

  double dndzs = 0.;

/////////////////////////////////
  double z_asked  = zs;
  double phig = 0.;

  // if(z_asked<V->pclass_sz->normalized_dndz_z[0])
  //    phig = 1e-100;
  // else if (z_asked>V->pclass_sz->normalized_dndz_z[V->pclass_sz->normalized_dndz_size-1])
  //    phig = 1e-100;
  // else  phig =  pwl_value_1d(V->pclass_sz->normalized_dndz_size,
  //                              V->pclass_sz->normalized_dndz_z,
  //                              V->pclass_sz->normalized_dndz_phig,
  //                              z_asked);
if(z_asked<V->pclass_sz->normalized_dndz_ngal_z[V->index_g][0])
   phig = 0.;//1e-100;
else if (z_asked>V->pclass_sz->normalized_dndz_ngal_z[V->index_g][V->pclass_sz->normalized_dndz_ngal_size[V->index_g]-1])
   phig = 0.;//1e-100;
else  phig =  pwl_value_1d(V->pclass_sz->normalized_dndz_ngal_size[V->index_g],
                         V->pclass_sz->normalized_dndz_ngal_z[V->index_g],
                         V->pclass_sz->normalized_dndz_ngal_phig[V->index_g],
                         z_asked);



// phig = get_galaxy_number_counts(z_asked,V->pclass_sz);
 dndzs = phig;
////////////////////////////////

  integrand = dndzs*W;

  integrand *= (1.+zs);

  // printf("-> integrand z = %.3e phig = =%.3e\n",z_asked,phig);

  return integrand;

}



int redshift_int_nlensmag(
                  int index_g,
                  struct class_sz_structure * pclass_sz,
                  struct background * pba,
                  double * pvectsz,
                  double * result
                   ) {

double z =  pvectsz[pclass_sz->index_z];

double zs_min = z;
double zs_max = pclass_sz->z2SZ;


struct Parameters_for_integrand_nlensmag V;
  V.pvectsz = pvectsz;
  V.pclass_sz = pclass_sz;
  V.pba = pba;
  V.index_g = index_g;


  void * params = &V;
  double r; //result of the integral

  double epsrel = 1e-6;
  double epsabs = 1e-30;
  //int show_neval = pclass_sz->patterson_show_neval;

  r=Integrate_using_Patterson_adaptive(log(1.+zs_min),
                                        log(1.+zs_max),
                                        epsrel, epsabs,
                                        integrand_nlensmag,
                                        params,0);


  *result = r;
                     }


int redshift_int_lensmag(
                  struct class_sz_structure * pclass_sz,
                  struct background * pba,
                  double * pvectsz,
                  double * result
                   ) {

double z =  pvectsz[pclass_sz->index_z];

double zs_min = z;
double zs_max = pclass_sz->z2SZ;


struct Parameters_for_integrand_lensmag V;
  V.pvectsz = pvectsz;
  V.pclass_sz = pclass_sz;
  V.pba = pba;


  void * params = &V;
  double r; //result of the integral

  double epsrel = 1e-6;
  double epsabs = 1e-30;
  //int show_neval = pclass_sz->patterson_show_neval;

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
                   struct class_sz_structure * pclass_sz,
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
              pclass_sz->ln_k_size_for_tSZ*index_num*sizeof(double),
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



      for (i=0;i<pclass_sz->ln_k_size_for_tSZ;i++) {
        k=exp(pclass_sz->ln_k_for_tSZ[i]);
        if (i == (pclass_sz->ln_k_size_for_tSZ-1)) k *= 0.9999999;
// printf("ok k = %e I = %e\n",k,pk*W*W);
    // //Input: wavenumber in 1/Mpc
    // //Output: total matter power spectrum P(k) in \f$ Mpc^3 \f$
  enum pk_outputs pk_for_vrms2;
  if (pclass_sz->pk_nonlinear_for_vrms2 == 1){
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
                          pclass_sz->ln_k_size_for_tSZ,
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
                                        pclass_sz->ln_k_size_for_tSZ,
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
                       struct class_sz_structure * pclass_sz,
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
              pclass_sz->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);


  // printf("entering2 \n");

  for (i=0;i<pclass_sz->ln_k_size_for_tSZ;i++) {
    k=exp(pclass_sz->ln_k_for_tSZ[i]);
    if (i == (pclass_sz->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
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
                          pclass_sz->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        pclass_sz->ln_k_size_for_tSZ,
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
                             struct class_sz_structure * pclass_sz,
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
              pclass_sz->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);

  for (i=0;i<pclass_sz->ln_k_size_for_tSZ;i++) {
    k=exp(pclass_sz->ln_k_for_tSZ[i]);
    if (i == (pclass_sz->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
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
                          pclass_sz->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_spline(array_for_sigma,
                                        index_num,
                                        pclass_sz->ln_k_size_for_tSZ,
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
// int integrate_over_m_at_z_spline(struct class_sz_structure * pclass_sz,
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
//               pclass_sz->ln_M_size*index_num*sizeof(double),
//               pclass_sz->error_message);
//
//   for (i=0;i<pclass_sz->ln_M_size;i++) {
//     x=exp(pclass_sz->ln_x_for_pp[i]);
//
//     p_gnfw(&p_gnfw_at_x,x,pvectsz,pba,pclass_sz);
//
//
//
//     double pp_at_x_and_ell_over_ell_char = x*p_gnfw_at_x;
//     array_for_integral[i*index_num+index_x]= log(x);
//     array_for_integral[i*index_num+index_y]= pp_at_x_and_ell_over_ell_char;
//   }
//
//
//
//
//   class_call(array_spline(array_for_integral,
//                           index_num,
//                           pclass_sz->ln_x_size_for_pp,
//                           index_x,
//                           index_y,
//                           index_ddy,
//                           _SPLINE_EST_DERIV_,
//                           pclass_sz->error_message),
//              pclass_sz->error_message,
//              pclass_sz->error_message);
//
//   class_call(array_integrate_all_spline(array_for_integral,
//                                         index_num,
//                                         pclass_sz->ln_x_size_for_pp,
//                                         index_x,
//                                         index_y,
//                                         index_ddy,
//                                         result,
//                                         pclass_sz->error_message),
//              pclass_sz->error_message,
//              pclass_sz->error_message);
//
//   free(array_for_integral);
//
//   return _SUCCESS_;
//
// }
//


struct Parameters_for_integrand_gas_pressure_profile{
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
  double kl;
};


struct Parameters_for_integrand_nfw_profile{
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
};


double integrand_nfw_profile(double x, void *p){

  struct Parameters_for_integrand_nfw_profile *V = ((struct Parameters_for_integrand_nfw_profile *) p);

    double nfw_profile_at_x = 0.;
    rho_gnfw(&nfw_profile_at_x,x,V->pvectsz,V->pba,V->pclass_sz);
    double result = nfw_profile_at_x;
  return result;

}

double integrand_gas_pressure_profile(double x, void *p){

  struct Parameters_for_integrand_gas_pressure_profile *V = ((struct Parameters_for_integrand_gas_pressure_profile *) p);

    //double x=exp(ln_x);

    double p_gnfw_at_x = 0.;
    p_gnfw(&p_gnfw_at_x,x,V->kl,V->pvectsz,V->pba,V->pclass_sz);

    double result = p_gnfw_at_x;

  return result;

}

int two_dim_ft_nfw_profile(struct class_sz_structure * pclass_sz,
                          struct background * pba,
                          double * pvectsz,
                          double * result
                          ) {

  struct Parameters_for_integrand_nfw_profile V;
  V.pclass_sz = pclass_sz;
  V.pba = pba;
  V.pvectsz = pvectsz;


  void * params = &V;

  gsl_function F;
  F.function = &integrand_nfw_profile;
  F.params = params;

  double eps_abs = pclass_sz->nfw_profile_epsabs;
  double eps_rel = pclass_sz->nfw_profile_epsrel;

  double result_gsl, error;

  double xin = 1.e-5;
  double c_nfw;

  //Battaglia 16 case:

  // double rvir = pvectsz[pclass_sz->index_rVIR]; //in Mpc/h
  // double r200c = pvectsz[pclass_sz->index_r200c]; //in Mpc/h
  // double rs = pvectsz[pclass_sz->index_rs]; //in Mpc/h
  // xout = 50.*rvir/rs; // as in hmvec (default 20, but set to 50 in example file)
  double xout = pclass_sz->x_out_truncated_nfw_profile_electrons; // as in hmvec (default 20, but set to 50 in example file) // is this value ok?

  if (pclass_sz->use_xout_in_density_profile_from_enclosed_mass){
    xout = get_m_to_xout_at_z_and_m(pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_m200c],pclass_sz);

    // printf("xout = %.5e\n",xout);
  }
  c_nfw = 1.;


// QAWO

  double delta_l = xout - xin;

  gsl_integration_workspace * w;
  gsl_integration_qawo_table * wf;

  int size_w = 20000; // was 3000... not sure if it matters
  w = gsl_integration_workspace_alloc(size_w);


  double w0;

  int index_md = (int) pvectsz[pclass_sz->index_md];
  double y_eff;
  // y_eff = (pvectsz[pclass_sz->index_multipole_for_nfw_profile]+0.5)
  //            /pvectsz[pclass_sz->index_characteristic_multipole_for_nfw_profile];

  //pclass_sz->index_multipole_for_nfw_profile is k
  y_eff = pvectsz[pclass_sz->index_multipole_for_nfw_profile]*pvectsz[pclass_sz->index_r200c]*(1.+pvectsz[pclass_sz->index_z]);

  w0 = y_eff;


  wf = gsl_integration_qawo_table_alloc(w0, delta_l,GSL_INTEG_SINE,300); // default 30


  int limit = size_w; //number of sub interval
  gsl_integration_qawo(&F,xin,eps_abs,eps_rel,limit,w,wf,&result_gsl,&error);

  *result = result_gsl;

  gsl_integration_qawo_table_free(wf);
  gsl_integration_workspace_free(w);

// QAWO --->> end
// try fft:



}


struct Parameters_for_integrand_bcm_profile_norm{
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double m;
  double z;
};
double integrand_bcm_profile_norm(double x, void *p)
{
  // double x = exp(lnx);
  // printf("being integrated\n");
  struct Parameters_for_integrand_bcm_profile_norm *V = ((struct Parameters_for_integrand_bcm_profile_norm *) p);
  double xout = V->pclass_sz->x_out_truncated_nfw_profile_electrons;
      if (x>xout){
        return 0.;
      }
      else{
      return  get_gas_profile_at_x_M_z_bcm_200c(x,
                                                V->m,
                                                V->z,
                                                V->pba,
                                                V->pclass_sz)*pow(x,2);
      }

}

struct Parameters_for_integrand_matter_density_profile_norm{
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double m;
  double z;
  double c_delta;
};
//
double integrand_matter_density_profile_norm(double x, void *p)
{
  // double x = exp(lnx);
  // printf("being integrated\n");
  struct Parameters_for_integrand_matter_density_profile_norm *V = ((struct Parameters_for_integrand_matter_density_profile_norm *) p);
  // double xout = V->pclass_sz->x_out_truncated_nfw_profile_electrons;
      // if (x>xout){
      //   return 0.;
      // }
      // else{
      return  get_nfw_with_power_law_profile_at_x(x,
                                                  V->pclass_sz->matter_nfw_power_law_index,
                                                  // V->m,
                                                  // V->z,
                                                  // V->pba,
                                                  V->c_delta*V->pclass_sz->x_out_matter_density_profile_normalization)*pow(x,2);
      // }

}



int rho_gnfw(double * rho_nfw_x,
            double x ,
            double * pvectsz,
            struct background * pba,
            struct class_sz_structure * pclass_sz)
{

 int index_md = (int) pvectsz[pclass_sz->index_md];

 double z = pvectsz[pclass_sz->index_z];

 // double y_eff;
 //   y_eff = (pvectsz[pclass_sz->index_multipole_for_nfw_profile]+0.5)
 //            /pvectsz[pclass_sz->index_characteristic_multipole_for_nfw_profile];

 double y_eff;
 // y_eff = (pvectsz[pclass_sz->index_multipole_for_nfw_profile]+0.5)
 //            /pvectsz[pclass_sz->index_characteristic_multipole_for_nfw_profile];
 y_eff = pvectsz[pclass_sz->index_multipole_for_nfw_profile]*pvectsz[pclass_sz->index_r200c]*(1.+pvectsz[pclass_sz->index_z]);

    double A_rho0 = pclass_sz->A_rho0;
    double A_alpha = pclass_sz->A_alpha;
    double A_beta = pclass_sz->A_beta;

    double alpha_m_rho0 = pclass_sz->alpha_m_rho0;
    double alpha_m_alpha = pclass_sz->alpha_m_alpha;
    double alpha_m_beta = pclass_sz->alpha_m_beta;

    double alphap_m_rho0 = pclass_sz->alphap_m_rho0;
    double alphap_m_alpha = pclass_sz->alphap_m_alpha;
    double alphap_m_beta = pclass_sz->alphap_m_beta;

    double alpha_z_rho0 = pclass_sz->alpha_z_rho0;
    double alpha_z_alpha = pclass_sz->alpha_z_alpha;
    double alpha_z_beta = pclass_sz->alpha_z_beta;

  // Eq. A1 and A2:
  // double m200_over_msol = pvectsz[pclass_sz->index_m200c]/pba->h; // convert to Msun
  // double rho0 = 1.;
  // double alpha = A_alpha*pow(m200_over_msol/1e14,alpha_m_alpha)*pow(1.+z,alpha_z_alpha);
  // double beta = A_beta*pow(m200_over_msol/1e14,alpha_m_beta)*pow(1.+z,alpha_z_beta);

  double gamma = pclass_sz->gamma_B16;
  double xc = pclass_sz->xc_B16;

  double c_asked = pclass_sz->c_B16;
  *rho_nfw_x = get_gas_profile_at_x_M_z_b16_200c(x,
                                                 pvectsz[pclass_sz->index_m200c],
                                                 z,
                                                 c_asked,
                                                 A_rho0,
                                                 A_alpha,
                                                 A_beta,
                                                 alpha_m_rho0,
                                                 alpha_m_alpha,
                                                 alpha_m_beta,
                                                 alpha_z_rho0,
                                                 alpha_z_alpha,
                                                 alpha_z_beta,
                                                 pclass_sz->mcut,
					                                       pclass_sz->alphap_m_rho0,
                                                 pclass_sz->alphap_m_alpha,
                                                 pclass_sz->alphap_m_beta,
					                                       pclass_sz->alpha_c_rho0,
                                                 pclass_sz->alpha_c_alpha,
                                                 pclass_sz->alpha_c_beta,
                                                 gamma,
                                                 xc,
                                                 pba,
                                                 pclass_sz)*pow(x,2)/(x*y_eff);

}






double get_rho_crit_at_z(double z_asked,
                         struct background * pba,
                         struct class_sz_structure * pclass_sz){
double result;
double rho_crit;




double * pvecback;
double * pvectsz;



double tau;
double z = z_asked;
int first_index_back = 0;

class_alloc(pvecback,
            pba->bg_size*sizeof(double),
            pclass_sz->error_message);

// class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
//  int i;
//  for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

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




// pvectsz[pclass_sz->index_z] = z;
// pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
//                                 *pow(_Mpc_over_m_,1)
//                                 *pow(_c_,2)
//                                 *pvecback[pba->index_bg_rho_crit]
//                                 /pow(pba->h,2);
//
// rho_crit = pvectsz[pclass_sz->index_Rho_crit];

rho_crit = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);


free(pvecback);
// free(pvectsz);



result = rho_crit;

return result;
}


double get_gas_profile_at_x_M_z_nfw_200c(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct class_sz_structure * pclass_sz){
double result;
double r_asked = 0.;
double rho_s;
double delta;
double c_delta;
double r_delta;
double r_s;
double rho_crit;
double f_b = pclass_sz->f_b_gas;//pba->Omega0_b/pclass_sz->Omega_m_0;
double x;
double p_x;



double * pvecback;
double * pvectsz;



double tau;
double z = z_asked;
int first_index_back = 0;

class_alloc(pvecback,
            pba->bg_size*sizeof(double),
            pclass_sz->error_message);

class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
 int i;
 for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

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




pvectsz[pclass_sz->index_z] = z;
pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);

rho_crit = pvectsz[pclass_sz->index_Rho_crit];
delta = 200.;//*pvecback[pba->index_bg_Omega_m];
c_delta = get_c200c_at_m_and_z(m_asked,z,pba,pclass_sz);
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


double get_gas_profile_at_x_M_z_nfw_200m(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct class_sz_structure * pclass_sz){
double result;
double r_asked = 0.;
double rho_s;
double delta;
double c_delta;
double r_delta;
double r_s;
double rho_crit;
double f_b = pclass_sz->f_b_gas;//pba->Omega0_b/pclass_sz->Omega_m_0;
double x;
double p_x;



double * pvecback;
double * pvectsz;



double tau;
double z = z_asked;
int first_index_back = 0;

class_alloc(pvecback,
            pba->bg_size*sizeof(double),
            pclass_sz->error_message);

class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
 int i;
 for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

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




pvectsz[pclass_sz->index_z] = z;
pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);

rho_crit = pvectsz[pclass_sz->index_Rho_crit];
delta = 200.*pvecback[pba->index_bg_Omega_m];
c_delta = get_c200m_at_m_and_z(m_asked,z,pba,pclass_sz);
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


double get_mass_profile_at_x_M_z_nfw_200m(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         struct background * pba,
                                         struct class_sz_structure * pclass_sz){
double result;
double r_asked = 0.;
double rho_s;
double delta;
double c_delta;
double r_delta;
double r_s;
double rho_crit;
// double f_b = pclass_sz->f_b_gas;//pba->Omega0_b/pclass_sz->Omega_m_0;
double x;
double p_x;



double * pvecback;
double * pvectsz;



double tau;
double z = z_asked;
int first_index_back = 0;

class_alloc(pvecback,
            pba->bg_size*sizeof(double),
            pclass_sz->error_message);

class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
 int i;
 for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

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




pvectsz[pclass_sz->index_z] = z;
pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);

rho_crit = pvectsz[pclass_sz->index_Rho_crit];
delta = 200.*pvecback[pba->index_bg_Omega_m];
c_delta = get_c200m_at_m_and_z(m_asked,z,pba,pclass_sz);
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
result = rho_s*p_x;
return result;
}

double get_rvir_of_m200c_at_z(//double x_asked, // this is just radius
                              double m_asked,
                              double z,
                              struct background * pba,
                              struct class_sz_structure * pclass_sz){

    double result;
    double rvir;

    double * pvectsz;
    double * pvecback;
    double tau;
    int first_index_back;

    class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
    class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);


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


  pvectsz[pclass_sz->index_z] = z;



  pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
  double chi = sqrt(pvectsz[pclass_sz->index_chi2]);
  // pvectsz[pclass_sz->index_multipole_for_pressure_profile] = k*chi;
  // pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in p_gnfw for the pk mode computation

  pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);

  double omega = pvecback[pba->index_bg_Omega_m];
  pvectsz[pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
  pvectsz[pclass_sz->index_m200c] = m_asked;
  pvectsz[pclass_sz->index_r200c] = pow(3.*pvectsz[pclass_sz->index_m200c]/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc

  // double r_asked = x_asked*pvectsz[pclass_sz->index_r200c];

  class_call(mDEL_to_mVIR(pvectsz[pclass_sz->index_m200c],
                          200.*(pvectsz[pclass_sz->index_Rho_crit]),
                          pvectsz[pclass_sz->index_Delta_c],
                          pvectsz[pclass_sz->index_Rho_crit],
                          z,
                          &pvectsz[pclass_sz->index_mVIR],
                          pclass_sz,
                          pba),
                  pclass_sz->error_message,
                  pclass_sz->error_message);
 //
 //  // rvir needed in bcm model
  rvir = evaluate_rvir_of_mvir(pvectsz[pclass_sz->index_mVIR],pvectsz[pclass_sz->index_Delta_c],pvectsz[pclass_sz->index_Rho_crit],pclass_sz);


free(pvectsz);
free(pvecback);

return rvir;


}






double get_gas_profile_at_x_M_z_bcm_200c(double x_asked, // this is just radius
                                         double m_asked,
                                         double z,
                                         struct background * pba,
                                         struct class_sz_structure * pclass_sz){
    double result;
    double rvir;

    double * pvectsz;
    double * pvecback;
    double tau;
    int first_index_back;

    class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
    class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);


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


  pvectsz[pclass_sz->index_z] = z;



  pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
  double chi = sqrt(pvectsz[pclass_sz->index_chi2]);
  // pvectsz[pclass_sz->index_multipole_for_pressure_profile] = k*chi;
  // pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in p_gnfw for the pk mode computation

  pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);

  double omega = pvecback[pba->index_bg_Omega_m];
  pvectsz[pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
  pvectsz[pclass_sz->index_m200c] = m_asked;
  pvectsz[pclass_sz->index_r200c] = pow(3.*pvectsz[pclass_sz->index_m200c]/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.); //in units of h^-1 Mpc

  double r_asked = x_asked*pvectsz[pclass_sz->index_r200c];

  class_call(mDEL_to_mVIR(pvectsz[pclass_sz->index_m200c],
                          200.*(pvectsz[pclass_sz->index_Rho_crit]),
                          pvectsz[pclass_sz->index_Delta_c],
                          pvectsz[pclass_sz->index_Rho_crit],
                          z,
                          &pvectsz[pclass_sz->index_mVIR],
                          pclass_sz,
                          pba),
                  pclass_sz->error_message,
                  pclass_sz->error_message);
 //
 //  // rvir needed in bcm model
  rvir = evaluate_rvir_of_mvir(pvectsz[pclass_sz->index_mVIR],pvectsz[pclass_sz->index_Delta_c],pvectsz[pclass_sz->index_Rho_crit],pclass_sz);


double omega_b_over_omega_m = pclass_sz->f_b_gas;
double fstar = get_fstar_of_m_at_z(m_asked,z,pclass_sz);
double num = omega_b_over_omega_m-fstar;
double delta = pclass_sz->delta_bcm;
double gamma = pclass_sz->gamma_bcm;
double thetaej = pclass_sz->theta_ej_bcm;
double mu = pclass_sz->mu_bcm;
double mc;// = pow(10.,pclass_sz->log10Mc_bcm);

// include redshift dependence of Mc:
double mc_z = pclass_sz->log10Mc_bcm*pow(1.+z,pclass_sz->nu_log10Mc_bcm);
mc = pow(10.,mc_z);

double betam = 3.*pow(m_asked/mc,mu)/(1.+pow(m_asked/mc,mu));
double den1 = pow(1.+10.*r_asked/rvir,betam);
double den2 = pow(1.+pow(r_asked/thetaej/rvir,gamma),(delta-betam)/gamma);
result = num/den1/den2;

free(pvectsz);
free(pvecback);

return result;


}




double get_gas_profile_at_x_M_z_b16_200c(double x_asked,
                                         double m_asked,
                                         double z_asked,
                                         double c_asked,
                                         double A_rho0,
                                         double A_alpha,
                                         double A_beta,
                                         double alpha_m_rho0,
                                         double alpha_m_alpha,
                                         double alpha_m_beta,
                                         double alpha_z_rho0,
                                         double alpha_z_alpha,
                                         double alpha_z_beta,
                            						 double mcut,
                            						 double alphap_m_rho0,
                            						 double alphap_m_alpha,
                            						 double alphap_m_beta,
                            						 double alpha_c_rho0,
                            						 double alpha_c_alpha,
                            						 double alpha_c_beta,
                                         double gamma,
                                         double xc,
                                         struct background * pba,
                                         struct class_sz_structure * pclass_sz){
double result;
double r_asked = 0.;
double rho_s;
double delta;
double c_delta;
double r_delta;
double r_s;
double rho_crit;
double f_b = pclass_sz->f_b_gas;//pba->Omega0_b/pclass_sz->Omega_m_0;
double x;
double p_x;



double * pvecback;
double * pvectsz;



double tau;
double z = z_asked;
int first_index_back = 0;

class_alloc(pvecback,
            pba->bg_size*sizeof(double),
            pclass_sz->error_message);

class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
 int i;
 for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

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




pvectsz[pclass_sz->index_z] = z;
pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);

rho_crit = pvectsz[pclass_sz->index_Rho_crit];
delta = 200.;
//pvectsz[pclass_sz->index_m200c] = m_asked;

c_delta = 1.;
r_delta = pow(3.*m_asked/(4.*_PI_*delta*rho_crit),1./3.); //in units of h^-1 Mpc

// rho_s = pow(c_delta,3.)*delta*rho_crit/3./m_nfw(c_delta);

r_s = r_delta/c_delta;
// x = r_asked/r_s;
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
  // if (pclass_sz->tau_profile_mode == 0){
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
  // else if (pclass_sz->tau_profile_mode == 1){
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
  double rho0;
  double alpha;
  double beta;
  // printf("mcut = %.5e %.5e\n",mcut,pclass_sz->mcut);
  if (m200_over_msol > mcut) {
  // rho0 = A_rho0*pow(m200_over_msol/1e14,alpha_m_rho0)*pow(1.+z,alpha_z_rho0);
  // alpha = A_alpha*pow(m200_over_msol/1e14,alpha_m_alpha)*pow(1.+z,alpha_z_alpha);
  // beta = A_beta*pow(m200_over_msol/1e14,alpha_m_beta)*pow(1.+z,alpha_z_beta);
  rho0 = A_rho0*pow(m200_over_msol/mcut,alpha_m_rho0)*pow(1.+z,alpha_z_rho0)*pow(1.+c_asked,alpha_c_rho0);
  alpha = A_alpha*pow(m200_over_msol/mcut,alpha_m_alpha)*pow(1.+z,alpha_z_alpha)*pow(1.+c_asked,alpha_c_alpha);
  beta = A_beta*pow(m200_over_msol/mcut,alpha_m_beta)*pow(1.+z,alpha_z_beta)*pow(1.+c_asked,alpha_c_beta);
  }
  else{

  rho0 = A_rho0*pow(m200_over_msol/mcut,alphap_m_rho0)*pow(1.+z,alpha_z_rho0)*pow(1.+c_asked,alpha_c_rho0);
  alpha = A_alpha*pow(m200_over_msol/mcut,alphap_m_alpha)*pow(1.+z,alpha_z_alpha)*pow(1.+c_asked,alpha_c_alpha);
  beta = A_beta*pow(m200_over_msol/mcut,alphap_m_beta)*pow(1.+z,alpha_z_beta)*pow(1.+c_asked,alpha_c_beta);
  }
  // double gamma = -0.2;
  // double xc = 0.5;

  p_x = pow(x/xc,gamma)*pow(1.+ pow(x/xc,alpha),-(beta+gamma)/alpha);
  // p_x = m200_over_msol;

  //
  // double rho0 = A_rho0*pow(m200_over_msol/1e14,alpha_m_rho0)*pow(1.+z,alpha_z_rho0);
  // tau_normalisation = rho0*4.*_PI_*pow(pvectsz[pclass_sz->index_r200c],3)
  //                     *pvectsz[pclass_sz->index_Rho_crit];


// result = rho0*rho_crit*f_b*p_x/rho_crit/f_b;
result = rho0*rho_crit*p_x*f_b;
return result;
}

double get_HI_density_profile_at_k_M_z(double k_asked, double m_asked, double z_asked, struct class_sz_structure * pclass_sz){

return 0.;
}





double get_gas_density_profile_at_k_M_z(double k_asked, double m_asked, double z_asked, int normalize, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
  double k = log(k_asked);

   if (z<pclass_sz->array_profile_ln_1pz[0])
    return 0.;//z = pclass_sz->array_profile_ln_1pz[0];
 if (z>pclass_sz->array_profile_ln_1pz[pclass_sz->n_z_density_profile-1])
    return 0.;//z = pclass_sz->array_profile_ln_1pz[pclass_sz->n_z_density_profile-1];

 if (m<pclass_sz->array_profile_ln_m[0])
    return 0.;//m = pclass_sz->array_profile_ln_m[0];
 if (m>pclass_sz->array_profile_ln_m[pclass_sz->n_m_density_profile-1])
    return 0.;//m =  pclass_sz->array_profile_ln_m[pclass_sz->n_m_density_profile-1];

if (k<pclass_sz->array_profile_ln_k[0])
    return 0.;//l = pclass_sz->array_profile_ln_k[0];
 if (k>pclass_sz->array_profile_ln_k[pclass_sz->n_k_density_profile-1])
    return 0.;//l =  pclass_sz->array_profile_ln_k[pclass_sz->n_ell_density_profile-1];



  // if (pclass_sz->tau_profile == 1){
  // find the closest l's in the grid:
  int id_k_low;
  int id_k_up;
  int n_k = pclass_sz->n_k_density_profile;
  int n_m = pclass_sz->n_m_density_profile;
  int n_z = pclass_sz->n_z_density_profile;
  r8vec_bracket(n_k,pclass_sz->array_profile_ln_k,k,&id_k_low,&id_k_up);

  // interpolate 2d at l_low:

 double ln_rho_low = pwl_interp_2d(
                                n_z,
                                n_m,

                                pclass_sz->array_profile_ln_1pz,
                                pclass_sz->array_profile_ln_m,
                                pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[id_k_low-1],
                                1,
                                &z,
                                &m);

 double ln_rho_up = pwl_interp_2d(
                                n_z,
                                n_m,
                                pclass_sz->array_profile_ln_1pz,
                                pclass_sz->array_profile_ln_m,
                                pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[id_k_up-1],
                                1,
                                &z,
                                &m);
 double ln_k_low = pclass_sz->array_profile_ln_k[id_k_low-1];
 double ln_k_up = pclass_sz->array_profile_ln_k[id_k_up-1];

 double result = exp(ln_rho_low + ((k - ln_k_low) / (ln_k_up - ln_k_low)) * (ln_rho_up - ln_rho_low));

// BCM needs to be normalized
if (normalize == 1){
  double norm = get_normalization_gas_density_profile(z_asked,m_asked,pclass_sz);
  result *= 1./norm;
  // printf("norm = %.5e\n",norm);
}

 if (isnan(result) || isinf(result)){
 printf("in get gas: z %.8e m %.8e k %.8e  ln_rho_low  %.8e ln_rho_low  %.8e id_k_low %d\n",z_asked,m_asked,k_asked,ln_rho_low,ln_rho_up,id_k_low);
 exit(0);
}
 return result;


}




double integrand_gas_pressure_profile_2h(double lnM_halo, void *p){

  struct Parameters_for_integrand_gas_pressure_profile_2h *V = ((struct Parameters_for_integrand_gas_pressure_profile_2h *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double kl = V->k;


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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);


      if (V->pclass_sz->sz_verbose>12)
        printf("doing mass conversion and hmf.\n");
      V->pvectsz[V->pclass_sz->index_has_electron_pressure] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      if (V->pclass_sz->sz_verbose>12)
        printf("mass conversion done.\n");

      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      V->pvectsz[V->pclass_sz->index_md] = 0;
      evaluate_pressure_profile(kl,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      double gas_profile_at_k_1 = V->pvectsz[V->pclass_sz->index_pressure_profile];
      if (isnan(gas_profile_at_k_1) || isinf(gas_profile_at_k_1)){
        printf("nan in pressure profile tranfform for 2halo.\n");
        exit(0);
      }

      // double gas_profile_at_k_1 = get_gas_pressure_profile_at_k_m_z(kl,
      //                                                               V->pvectsz[V->pclass_sz->index_mass_for_electron_pressure],
      //                                                               z,
      //                                                               V->pclass_sz);


      // here we need to convert back to P rather than yl:
      double conv_pe_to_y;
      // double sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb = 283./0.5176; //1./0.5176=1.932=(5Xh+3)/2(Xh+1) with Xh = 0.76 and Pth=1.932Pe
      // (Xh is the primodial hydrogen mass fraction)
      // more accurate version (see explanation below):
      // in units of Mpc^-1*micro Kelvins
      // double sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb = 283.2980000259841/0.5176*V->pba->T_cmb/2.725;

      conv_pe_to_y =  V->pclass_sz->sigmaT_over_mec2_times_50eV_per_cm3_times_Tcmb // here Tcmb is in muK
                       /50. // to cancel the factor 50 above 50eV/cm^3
                       /V->pba->T_cmb
                       // *pressure_normalisation // what we get in get_pressure
                       // *pvectsz[pclass_sz->index_pressure_profile] // what we get in get_pressure
                       // *(4*_PI_) // fourier transform factors
                       // *pow(characteristic_multipole,-2) // fourier transform factors
                       // *characteristic_radius //rs in Mpc // fourier transform factors
                       *V->pclass_sz->Tcmb_gNU;

      double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      double d_A = chi/(1.+z);

      gas_profile_at_k_1 = gas_profile_at_k_1/conv_pe_to_y*pow(d_A,2);


      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ppt,V->pclass_sz);
      double b1 = V->pvectsz[V->pclass_sz->index_halo_bias];
      double result = hmf*b1*gas_profile_at_k_1;



  return result;

}





double integrand_gas_density_profile_2h(double lnM_halo, void *p){

  struct Parameters_for_integrand_gas_density_profile_2h *V = ((struct Parameters_for_integrand_gas_density_profile_2h *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double kl = V->k;


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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);

      V->pvectsz[V->pclass_sz->index_has_electron_density] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];


      int normalize = 0;
      if (V->pclass_sz->tau_profile == 2)
        normalize = 1;
      double k_asked = kl*(1.+z)*V->pvectsz[V->pclass_sz->index_radius_for_electron_density];
      // double k_asked = kl;
      double gas_profile_at_k_1 = get_gas_density_profile_at_k_M_z(k_asked,
                                                                   V->pvectsz[V->pclass_sz->index_mass_for_electron_density],
                                                                   z,
                                                                   normalize,
                                                                   V->pclass_sz);


      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ppt,V->pclass_sz);
      double b1 = V->pvectsz[V->pclass_sz->index_halo_bias];
      double result = hmf*b1*gas_profile_at_k_1;



  return result;

}


struct Parameters_for_integrand_rho_2h_qawo{
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  struct primordial * ppm;
  struct nonlinear * pnl;
  double z;
  double r;
};


int rho_2h_qawo(double * rho_nfw_x,
                   double x ,
                   double r,
                   double z,
                   struct background * pba,
                   struct nonlinear * pnl,
                   struct primordial * ppm,
                   struct class_sz_structure * pclass_sz){

double rho2h = get_rho_2h_at_k_and_z(x,z,pclass_sz);
double pklin = get_pk_lin_at_k_and_z(x,z,pba,ppm,pnl,pclass_sz);

*rho_nfw_x = rho2h*pklin*pow(x,2)/(r*x)/2./_PI_/_PI_;
                   }

double integrand_rho_2h_qawo(double x, void *p){

  struct Parameters_for_integrand_rho_2h_qawo *V = ((struct Parameters_for_integrand_rho_2h_qawo *) p);

  double result = 0.;
  rho_2h_qawo(&result,x,V->r,
                        V->z,
                        V->pba,
                        V->pnl,
                        V->ppm,
                        V->pclass_sz);
  return result;

}



int tabulate_gas_density_profile_2h_fft_at_z_and_r(struct background * pba,
                                                   struct nonlinear * pnl,
                                                   struct primordial * ppm,
                                                   struct class_sz_structure * pclass_sz){

  if (pclass_sz->sz_verbose>10) printf("-> ftabulateing density profile kmz.\n");
  // exit(0);

double k_min = pclass_sz->k_min_samp_fftw;
double k_max = pclass_sz->k_max_samp_fftw; // this is a precision parameter
// tabulate the integrand in the "l" dimension:
const int N = pclass_sz->N_samp_fftw;


class_alloc(pclass_sz->array_profile_rho_2h_at_r_and_z,
            N*pclass_sz->n_z_density_profile*sizeof(double),
            pclass_sz->error_message);
class_alloc(pclass_sz->array_profile_ln_r,
            N*sizeof(double),
            pclass_sz->error_message);

int index_z;

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
pclass_sz,pba,ppm,pnl)\
private(tstart, tstop,index_z) \
num_threads(number_of_threads)
{
#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_density_profile; index_z++){
#pragma omp flush(abort)
  double z = exp(pclass_sz->array_profile_ln_1pz[index_z])-1.;
  double k[N], Pk1[N];
  int index_k;
  for (index_k=0; index_k<N; index_k++)
  {

    k[index_k] = exp(log(k_min)+index_k/(N-1.)*(log(k_max)-log(k_min)));
    Pk1[index_k] = get_rho_2h_at_k_and_z(k[index_k],z,pclass_sz);
    Pk1[index_k] *= get_pk_lin_at_k_and_z(k[index_k],z,pba,ppm,pnl,pclass_sz);
    // printf("z = %.3e k = %.5e pk1 = %.5e\n",z,k[index_k],Pk1[index_k]);
  }

  double rp[N], xi1[N];
  // pk2xi(N,k,Pk1,rp,xi1,pclass_sz);
  /* Compute the function
   *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
   * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
   * in this notation.  The input k-values must be logarithmically spaced.  The
   * resulting xi_l^m(r) will be evaluated at the dual r-values
   *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
  //void fftlog_ComputeXiLM(int l, int m, int N, const double k[],  const double pk[], double r[], double xi[]);
  fftlog_ComputeXiLMsloz(0, 2, N, k,  Pk1, rp, xi1,pclass_sz);
  // printf("\n##############\n");



  for (index_k=0; index_k<N; index_k++){
    int index_k_z = index_k * pclass_sz->n_z_density_profile + index_z;
    pclass_sz->array_profile_rho_2h_at_r_and_z[index_k_z] = xi1[index_k];
    pclass_sz->array_profile_ln_r[index_k] = log(rp[index_k]);

  // //// try alternative integration to double check result of fft:
  // // QAWO -- this currently gives the correct result
  // gsl_function F;
  // struct Parameters_for_integrand_rho_2h_qawo V;
  // V.pclass_sz = pclass_sz;
  // V.pba = pba;
  // V.ppm = ppm;
  // V.pnl = pnl;
  // V.z = z;
  // V.r = rp[index_k];
  // void * params = &V;
  //
  // F.function = &integrand_rho_2h_qawo;
  // F.params = params;
  // gsl_integration_workspace * w;
  // gsl_integration_qawo_table * wf;
  // int size_w = 8000;
  // double w0;
  // w = gsl_integration_workspace_alloc(size_w);
  // double xout = 100.;
  // double xin = 1e-2;
  // double delta_l = xout - xin;
  // w0 = rp[index_k];
  // wf = gsl_integration_qawo_table_alloc(w0, delta_l,GSL_INTEG_SINE,200);
  // int limit = size_w;
  // double result_gsl, error;
  // double eps_abs = 1e-10;
  // double eps_rel = 1e-3;
  // gsl_integration_qawo(&F,xin,eps_abs,eps_rel,limit,w,wf,&result_gsl,&error);
  // // *result = result_gsl;
  // gsl_integration_qawo_table_free(wf);
  // gsl_integration_workspace_free(w);
  //
  //
  // // printf("z = %.3e r = %.5e xi1 = %.5e gsl = %.5e ration gsl/xi1 = %.3e\n",z,
  // // rp[index_k],xi1[index_k],result_gsl,result_gsl/xi1[index_k]);
  //
  // pclass_sz->array_profile_rho_2h_at_r_and_z[index_k_z] = result_gsl;
   }


}
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in tab density profile 2h parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

}




int tabulate_gas_pressure_profile_2h_fft_at_z_and_r(struct background * pba,
                                                   struct nonlinear * pnl,
                                                   struct primordial * ppm,
                                                   struct class_sz_structure * pclass_sz){


double k_min = pclass_sz->k_min_samp_fftw;
double k_max = pclass_sz->k_max_samp_fftw; // this is a precision parameter
// tabulate the integrand in the "l" dimension:
const int N = pclass_sz->N_samp_fftw;


class_alloc(pclass_sz->array_pressure_profile_pressure_2h_at_r_and_z,
            N*pclass_sz->n_z_pressure_profile*sizeof(double),
            pclass_sz->error_message);
class_alloc(pclass_sz->array_pressure_profile_ln_r,
            N*sizeof(double),
            pclass_sz->error_message);

int index_z;

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
pclass_sz,pba,ppm,pnl)\
private(tstart, tstop,index_z) \
num_threads(number_of_threads)
{
#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_pressure_profile; index_z++){
#pragma omp flush(abort)
  double z = exp(pclass_sz->array_pressure_profile_ln_1pz[index_z])-1.;
  double k[N], Pk1[N];
  int index_k;
  for (index_k=0; index_k<N; index_k++)
  {

    k[index_k] = exp(log(k_min)+index_k/(N-1.)*(log(k_max)-log(k_min)));
    Pk1[index_k] = get_gas_pressure_2h_at_k_and_z(k[index_k],z,pclass_sz);
    Pk1[index_k] *= get_pk_lin_at_k_and_z(k[index_k],z,pba,ppm,pnl,pclass_sz);
    // printf("z = %.3e k = %.5e pk1 = %.5e\n",z,k[index_k],Pk1[index_k]);
  }

  double rp[N], xi1[N];
  // pk2xi(N,k,Pk1,rp,xi1,pclass_sz);
  /* Compute the function
   *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
   * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
   * in this notation.  The input k-values must be logarithmically spaced.  The
   * resulting xi_l^m(r) will be evaluated at the dual r-values
   *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
  //void fftlog_ComputeXiLM(int l, int m, int N, const double k[],  const double pk[], double r[], double xi[]);
  fftlog_ComputeXiLMsloz(0, 2, N, k,  Pk1, rp, xi1,pclass_sz);
  // printf("\n##############\n");



  for (index_k=0; index_k<N; index_k++){
    int index_k_z = index_k * pclass_sz->n_z_pressure_profile + index_z;
    pclass_sz->array_pressure_profile_pressure_2h_at_r_and_z[index_k_z] = xi1[index_k];
    pclass_sz->array_pressure_profile_ln_r[index_k] = log(rp[index_k]);

  // //// try alternative integration to double check result of fft:
  // // QAWO -- this currently gives the correct result
  // gsl_function F;
  // struct Parameters_for_integrand_rho_2h_qawo V;
  // V.pclass_sz = pclass_sz;
  // V.pba = pba;
  // V.ppm = ppm;
  // V.pnl = pnl;
  // V.z = z;
  // V.r = rp[index_k];
  // void * params = &V;
  //
  // F.function = &integrand_rho_2h_qawo;
  // F.params = params;
  // gsl_integration_workspace * w;
  // gsl_integration_qawo_table * wf;
  // int size_w = 8000;
  // double w0;
  // w = gsl_integration_workspace_alloc(size_w);
  // double xout = 100.;
  // double xin = 1e-2;
  // double delta_l = xout - xin;
  // w0 = rp[index_k];
  // wf = gsl_integration_qawo_table_alloc(w0, delta_l,GSL_INTEG_SINE,200);
  // int limit = size_w;
  // double result_gsl, error;
  // double eps_abs = 1e-10;
  // double eps_rel = 1e-3;
  // gsl_integration_qawo(&F,xin,eps_abs,eps_rel,limit,w,wf,&result_gsl,&error);
  // // *result = result_gsl;
  // gsl_integration_qawo_table_free(wf);
  // gsl_integration_workspace_free(w);
  //
  //
  // // printf("z = %.3e r = %.5e xi1 = %.5e gsl = %.5e ration gsl/xi1 = %.3e\n",z,
  // // rp[index_k],xi1[index_k],result_gsl,result_gsl/xi1[index_k]);
  //
  // pclass_sz->array_profile_rho_2h_at_r_and_z[index_k_z] = result_gsl;
   }


} // end z loop

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in tab pressure profile 2h parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

}


// Tabulate 2-halo term of density profile on a [z - k] grid
int tabulate_gas_density_profile_2h(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct perturbs * ppt,
                                    struct class_sz_structure * pclass_sz){

// if (pclass_sz->has_pk_b_at_z_2h
// == _FALSE_
// )
//   return 0;


 // array of multipoles:

 int n_z = pclass_sz->n_z_density_profile;
 int n_k = pclass_sz->n_k_density_profile; // dimension of pclass_sz->k_for_pk_hm

 // array of redshifts:
 double ln_1pz_min = log(1.+pclass_sz->z1SZ);
 double ln_1pz_max = log(1.+pclass_sz->z2SZ);


//  class_alloc(pclass_sz->array_profile_2h_ln_1pz,sizeof(double *)*n_z,pclass_sz->error_message);
  int index_k;
  int index_z;
//
// for (index_z=0;
//      index_z<n_z;
//      index_z++)
// {
//   pclass_sz->array_profile_2h_ln_1pz[index_z] = ln_1pz_min
//                                    +index_z*(ln_1pz_max-ln_1pz_min)
//                                    /(n_z-1.);
// }


class_alloc(pclass_sz->array_profile_ln_rho_2h_at_k_and_z,
            n_k*n_z*sizeof(double),
            pclass_sz->error_message);
// for (index_k=0;
//      index_k<n_k;
//      index_k++)
//     {
//     for (index_z=0;
//          index_z<n_z;
//          index_z++)
//          {
//           pclass_sz->array_profile_ln_rho_2h_at_k_and_z[index_k_z] = log(1e-100); // initialize with super small number
//           index_k_z += 1;
//          }
//     }


    double r;


    double m_min,m_max;
    m_min = pclass_sz->M1SZ; // for the mass integral
    m_max = pclass_sz->M2SZ; // for the mass integral

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
    pba,pclass_sz,ppm,ppt,pnl,m_min,m_max,n_k,n_z)\
    private(tstart, tstop,index_z,index_k,r) \
    num_threads(number_of_threads)
    {


    #ifdef _OPENMP
      tstart = omp_get_wtime();
    #endif

    double * pvecback;
    double * pvectsz;


     class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
       int i;
       for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

     class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


#pragma omp for collapse(2)
    for (index_k=0; index_k<n_k; index_k++)
    {
      for (index_z=0; index_z<n_z; index_z++)
              {
              double k = exp(pclass_sz->array_profile_ln_k[index_k]);
              int index_k_z = index_k * n_z + index_z;


              double z = exp(pclass_sz->array_profile_ln_1pz[index_z])-1.;


              // at each z, perform the mass integral
              struct Parameters_for_integrand_gas_density_profile_2h V;
              V.pnl = pnl;
              V.ppm = ppm;
              V.ppt = ppt;
              V.pclass_sz = pclass_sz;
              V.pba = pba;
              V.pvectsz = pvectsz;
              V.pvecback = pvecback;
              V.z = z;
              V.k = k;

              void * params = &V;
              double epsrel=1e-3;
              double epsabs=1e-100;

              // pvectsz[pclass_sz->index_has_electron_density] = 1; //


              r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                   epsrel, epsabs,
                                                   integrand_gas_density_profile_2h,
                                                   params,
                                                   pclass_sz->patterson_show_neval);


      // printf("k = %.5e, z = %.5e r = %.5e  zz = %.5e beofre ct\n",k,z,r,pvectsz[pclass_sz->index_z]);

       if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
         double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
         double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
         double I0 = integrand_gas_density_profile_2h(log(pclass_sz->m_min_counter_terms),params);
         double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
         r += bmin_umin;
         // printf("counter terms done r_m_1\n");
      }

      // printf("k = %.5e, z = %.5e r = %.5e after ct\n",k,z,r);

              pclass_sz->array_profile_ln_rho_2h_at_k_and_z[index_k_z] = log(r);
           }
         }

         free(pvecback);
         free(pvectsz);

         #ifdef _OPENMP
           tstop = omp_get_wtime();
           if (pclass_sz->sz_verbose > 0)
             printf("In %s: time spent in parallel region rho electons 2h (loop over z,k's) = %e s for thread %d\n",
                    __func__,tstop-tstart,omp_get_thread_num());
         #endif
    }
    if (abort == _TRUE_) return _FAILURE_;
    //end of parallel region
    return _SUCCESS_;

}




// Tabulate 2-halo term of pressure profile on a [z - k] grid
// currently for battaglia profile only (10.02.23)
int tabulate_gas_pressure_profile_2h(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct perturbs * ppt,
                                    struct class_sz_structure * pclass_sz){


 int n_z = pclass_sz->n_z_pressure_profile;
 int n_k = pclass_sz->n_k_pressure_profile_2h; // dimension of pclass_sz->k_for_pk_hm

  if (pclass_sz->sz_verbose>2)
    printf("setting up grid...\n");
    // printf("n_z = %d n_k = %d\n",n_z,n_k);


 // array of redshifts:
 double ln_1pz_min = log(1.+pclass_sz->z1SZ);
 double ln_1pz_max = log(1.+pclass_sz->z2SZ);


  int index_k;
  int index_z;


if (pclass_sz->sz_verbose>2)
  printf("assigning lnp array\n");

class_alloc(pclass_sz->array_pressure_profile_ln_pressure_2h_at_k_and_z,
            n_k*n_z*sizeof(double),
            pclass_sz->error_message);


// if (pclass_sz->sz_verbose>2)
//   printf("assigning lnp array done\n");

double ln_k_min = log(pclass_sz->k_min_gas_pressure_profile_2h);
double ln_k_max = log(pclass_sz->k_max_gas_pressure_profile_2h);

class_alloc(pclass_sz->array_pressure_profile_2h_ln_k,
            n_k*sizeof(double),
            pclass_sz->error_message);

if (pclass_sz->sz_verbose>2)
  printf("assigning lmnk array\n");

for (index_k=0;
     index_k<n_k;
     index_k++)
{
  pclass_sz->array_pressure_profile_2h_ln_k[index_k] = ln_k_min
                                              +index_k*(ln_k_max-ln_k_min)
                                              /(n_k-1.);
}



    double r;


    double m_min,m_max;
    m_min = pclass_sz->M1SZ; // for the mass integral
    m_max = pclass_sz->M2SZ; // for the mass integral

    double tstart, tstop;
    int abort;
    /* initialize error management flag */
    abort = _FALSE_;
    /* beginning of parallel region */
 // printf("-> start parallel n_z = %d n_k =%d\n",n_z,n_k);

if (pclass_sz->sz_verbose>2)
  printf("starting parallel block\n");

    int number_of_threads= 1;
    #ifdef _OPENMP
    #pragma omp parallel
      {
        number_of_threads = omp_get_num_threads();
      }
    #endif

    #pragma omp parallel \
    shared(abort,\
    pba,pclass_sz,ppm,ppt,pnl,m_min,m_max,n_k,n_z)\
    private(tstart, tstop,index_z,index_k,r) \
    num_threads(number_of_threads)
    {


    #ifdef _OPENMP
      tstart = omp_get_wtime();
    #endif

    double * pvecback;
    double * pvectsz;


     class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
       int i;
       for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

     class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


#pragma omp for collapse(2)
    for (index_k=0; index_k<n_k; index_k++)
    {
      for (index_z=0; index_z<n_z; index_z++)
              {

              double k = exp(pclass_sz->array_pressure_profile_2h_ln_k[index_k]);
              int index_k_z = index_k * n_z + index_z;


              double z = exp(pclass_sz->array_pressure_profile_ln_1pz[index_z])-1.;




              // at each z, perform the mass integral
              struct Parameters_for_integrand_gas_pressure_profile_2h V;
              V.pnl = pnl;
              V.ppm = ppm;
              V.ppt = ppt;
              V.pclass_sz = pclass_sz;
              V.pba = pba;
              V.pvectsz = pvectsz;
              V.pvecback = pvecback;
              V.z = z;
              V.k = k;

              void * params = &V;
              double epsrel=1e-3;
              double epsabs=1e-100;

              // pvectsz[pclass_sz->index_has_electron_density] = 1; //

              if (pclass_sz->sz_verbose>2)
                printf("-> starting integration of pressure profile 2h at z k %.5e %.5e\n",z,k);

              r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                   epsrel, epsabs,
                                                   integrand_gas_pressure_profile_2h,
                                                   params,
                                                   pclass_sz->patterson_show_neval);

if (pclass_sz->sz_verbose>2)
      printf("k = %.5e, z = %.5e r = %.5e  zz = %.5e beofre ct\n",k,z,r,pvectsz[pclass_sz->index_z]);

if (r<=0){
  // printf("getting r<0 after integrand_gas_pressure_profile_2h.\n");
  // exit(0);
  r = 1e-100;
}

       if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
         double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
         double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
         double I0 = integrand_gas_pressure_profile_2h(log(pclass_sz->m_min_counter_terms),params);
         double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
         r += bmin_umin;
         // printf("counter terms done r_m_1\n");
      }

      // printf("k = %.5e, z = %.5e r = %.5e after ct\n",k,z,r);

              pclass_sz->array_pressure_profile_ln_pressure_2h_at_k_and_z[index_k_z] = log(r);
           }
         }

         free(pvecback);
         free(pvectsz);

         #ifdef _OPENMP
           tstop = omp_get_wtime();
           if (pclass_sz->sz_verbose > 0)
             printf("In %s: time spent in parallel region electons pressure 2h (loop over z,k's) = %e s for thread %d\n",
                    __func__,tstop-tstart,omp_get_thread_num());
         #endif
    }
    if (abort == _TRUE_) return _FAILURE_;
    //end of parallel region
    return _SUCCESS_;

}







// double get_gas_density_profile_2h_at_r_M_z(double l_asked, double m_asked, double z_asked, struct class_sz_structure * pclass_sz){
//   double z = log(1.+z_asked);
//   double m = log(m_asked);
//   double l = log(l_asked);
//
//    if (z<pclass_sz->array_profile_ln_1pz[0])
//     return 0.;//z = pclass_sz->array_profile_ln_1pz[0];
//  if (z>pclass_sz->array_profile_ln_1pz[pclass_sz->n_z_density_profile-1])
//     return 0.;//z = pclass_sz->array_profile_ln_1pz[pclass_sz->n_z_density_profile-1];
//
//  if (m<pclass_sz->array_profile_ln_m[0])
//     return 0.;//m = pclass_sz->array_profile_ln_m[0];
//  if (m>pclass_sz->array_profile_ln_m[pclass_sz->n_m_density_profile-1])
//     return 0.;//m =  pclass_sz->array_profile_ln_m[pclass_sz->n_m_density_profile-1];
//
// if (l<pclass_sz->array_profile_ln_r[0])
//     return 0.;//l = pclass_sz->array_profile_ln_k[0];
//  if (l>pclass_sz->array_profile_ln_r[pclass_sz->n_r_density_profile-1])
//     return 0.;//l =  pclass_sz->array_profile_ln_k[pclass_sz->n_ell_density_profile-1];
//
//
//
//   // if (pclass_sz->tau_profile == 1){
//   // find the closest l's in the grid:
//   int id_l_low;
//   int id_l_up;
//   int n_ell = pclass_sz->n_k_density_profile;
//   int n_m = pclass_sz->n_m_density_profile;
//   int n_z = pclass_sz->n_z_density_profile;
//   r8vec_bracket(n_ell,pclass_sz->array_profile_ln_k,l,&id_l_low,&id_l_up);
//
//   // interpolate 2d at l_low:
//
//  double ln_rho_low = pwl_interp_2d(
//                                 n_z,
//                                 n_m,
//
//                                 pclass_sz->array_profile_ln_1pz,
//                                 pclass_sz->array_profile_ln_m,
//                                 pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[id_l_low-1],
//                                 1,
//                                 &z,
//                                 &m);
//
//  double ln_rho_up = pwl_interp_2d(
//                                 n_z,
//                                 n_m,
//                                 pclass_sz->array_profile_ln_1pz,
//                                 pclass_sz->array_profile_ln_m,
//                                 pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[id_l_up-1],
//                                 1,
//                                 &z,
//                                 &m);
//  double ln_l_low = pclass_sz->array_profile_ln_k[id_l_low-1];
//  double ln_l_up = pclass_sz->array_profile_ln_k[id_l_up-1];
//
//  double result = exp(ln_rho_low + ((l - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low));
//
// if (pclass_sz->normalize_gas_density_profile == 1){
//   double norm = get_normalization_gas_density_profile(z_asked,m_asked,pclass_sz)/pclass_sz->f_b_gas;
//   result *= 1./norm;
//   // printf("norm = %.5e\n",norm);
// }
//
//  if (isnan(result) || isinf(result)){
//  printf("in get gas: z %.8e m %.8e l %.8e  ln_rho_low  %.8e ln_rho_low  %.8e id_l_low %d\n",z_asked,m_asked,l_asked,ln_rho_low,ln_rho_up,id_l_low);
//  exit(0);
// }
//  return result;
//
//
// }



int tabulate_normalization_gas_density_profile(struct class_sz_structure *pclass_sz,struct background * pba){
  int n_m = pclass_sz->n_m_density_profile;
  int n_z = pclass_sz->n_z_density_profile;
  int index_m;
  int index_z;

  class_alloc(pclass_sz->array_ln_density_norm_at_z_and_m,
              sizeof(double *)*n_z*n_m,
              pclass_sz->error_message);

  double tstart, tstop;
  int abort;

  ///////////////////////////////////////////////
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
  pclass_sz,pba,n_z,n_m)\
  private(tstart, tstop, index_z,index_m) \
  num_threads(number_of_threads)
  {

  #ifdef _OPENMP
    tstart = omp_get_wtime();
  #endif


#pragma omp for collapse(2)
for (index_z=0; index_z<n_z; index_z++)
{
for (index_m=0; index_m<n_m; index_m++)
{
int index_z_m = index_m * n_z + index_z;
    // printf("index z = %d index_m = %d\n",index_z,index_m);

double z = exp(pclass_sz->array_profile_ln_1pz[index_z])-1.;
double m = exp(pclass_sz->array_profile_ln_m[index_m]);

pclass_sz->array_ln_density_norm_at_z_and_m[index_z_m] = 0.;

if (pclass_sz->tau_profile == 2){ // BCM

struct Parameters_for_integrand_bcm_profile_norm V;
V.pclass_sz = pclass_sz;
V.pba = pba;
V.z = z;
V.m = m;

void * params = &V;
double eps_rel = pclass_sz->density_norm_epsrel;
double eps_abs = pclass_sz->density_norm_epsabs;

// in the BCM paper this is taken very large
double x_out = pclass_sz->x_out_truncated_gas_density_profile_normalization;
double x_in = 0.;

double r = Integrate_using_Patterson_adaptive(x_in, x_out,
                                              eps_rel, eps_abs,
                                              integrand_bcm_profile_norm,
                                              params,pclass_sz->patterson_show_neval);

// need to multiply by 4pi r200c^3:
double r200c =  pow(m*3./4./_PI_/200./get_rho_crit_at_z(z,pba,pclass_sz),1./3.);
double omega_b_over_omega_m = pclass_sz->f_b_gas;
double fstar = get_fstar_of_m_at_z(m,z,pclass_sz);
double num = omega_b_over_omega_m-fstar;
double norm = 4.*_PI_*pow(r200c,3.)*r/m/num;
pclass_sz->array_ln_density_norm_at_z_and_m[index_z_m] = log(norm);

if (index_z == 71 && index_m == 58)
printf("index z = %d index_m = %d  z = %.4e m = %.4e lognorm = %.5e\n",index_z,index_m,z,m,pclass_sz->array_ln_density_norm_at_z_and_m[index_z_m]);

}
else{
  printf("normalization for this gas density profile not implemented yet.\n");
  exit(0);
}


    }
  }

  #ifdef _OPENMP
    tstop = omp_get_wtime();
    if (pclass_sz->sz_verbose > 0)
      printf("In %s: time spent in parallel region (loop over zm's) = %e s for thread %d\n",
             __func__,tstop-tstart,omp_get_thread_num());
  #endif


  // free(data);
  }


if (abort == _TRUE_) return _FAILURE_;

  return _SUCCESS_;
}



int tabulate_normalization_matter_density_profile(struct class_sz_structure *pclass_sz,struct background * pba){
  // printf("starting tabulating norm matter profile 0. \n");

  int n_m = pclass_sz->n_m_matter_density_profile;
  int n_z = pclass_sz->n_z_matter_density_profile;
  int index_m;
  int index_z;
  // printf("starting tabulating norm matter profile 1. \n");
  // printf("starting tabulating norm matter profile 1. %d %d\n",n_m,n_m);

  // exit(0);


  // printf("%.5e %d\n",pclass_sz->fixed_c200m,pclass_sz->n_m_matter_density_profile);

  // printf("starting tabulating norm matter profile 1. %d %d\n",n_m,n_z);


 class_alloc(pclass_sz->array_matter_density_profile_ln_m,sizeof(double *)*n_m,pclass_sz->error_message);
 // printf("starting tabulating norm matter profile 2a.\n");

 class_alloc(pclass_sz->array_matter_density_profile_ln_1pz,sizeof(double *)*n_z,pclass_sz->error_message);
 // printf("starting tabulating norm matter profile 2.\n");
//
 double ln_m_min = log(pclass_sz->M1SZ);
 double ln_m_max = log(pclass_sz->M2SZ);

 // array of redshifts:
 double ln_1pz_min = log(1.+pclass_sz->z1SZ);
 double ln_1pz_max = log(1.+pclass_sz->z2SZ);

// int index_m;
for (index_m=0;
     index_m<n_m;
     index_m++)
{
  pclass_sz->array_matter_density_profile_ln_m[index_m] = ln_m_min
                                      +index_m*(ln_m_max-ln_m_min)
                                      /(n_m-1.);
}

// int index_z;
for (index_z=0;
     index_z<n_z;
     index_z++)
{
  pclass_sz->array_matter_density_profile_ln_1pz[index_z] = ln_1pz_min
                                                       +index_z*(ln_1pz_max-ln_1pz_min)
                                                       /(n_z-1.);
}



  class_alloc(pclass_sz->array_ln_matter_density_norm_at_z_and_m,
              sizeof(double *)*n_z*n_m,
              pclass_sz->error_message);

  double tstart, tstop;
  int abort;

  ///////////////////////////////////////////////
  //Parallelization of computation
  /* initialize error management flag */
  abort = _FALSE_;
  /* beginning of parallel region */

  // printf("starting tabulating norm matter profile.\n");

  int number_of_threads= 1;
  #ifdef _OPENMP
  #pragma omp parallel
    {
      number_of_threads = omp_get_num_threads();
    }
  #endif

  #pragma omp parallel \
  shared(abort,\
  pclass_sz,pba,n_z,n_m)\
  private(tstart, tstop, index_z,index_m) \
  num_threads(number_of_threads)
  {

  #ifdef _OPENMP
    tstart = omp_get_wtime();
  #endif


#pragma omp for collapse(2)
for (index_z=0; index_z<n_z; index_z++)
{
for (index_m=0; index_m<n_m; index_m++)
{
int index_z_m = index_m * n_z + index_z;

double z = exp(pclass_sz->array_matter_density_profile_ln_1pz[index_z])-1.;
double m = exp(pclass_sz->array_matter_density_profile_ln_m[index_m]);

pclass_sz->array_ln_matter_density_norm_at_z_and_m[index_z_m] = 0.;

if (pclass_sz->profile_matter_density == 1){ // nfw with power law

double c_delta_matter;
  if (pclass_sz->delta_def_matter_density == 0){
    c_delta_matter = get_c200m_at_m_and_z(m,z,pba,pclass_sz);
  }
  else if (pclass_sz->delta_def_matter_density == 1){
    c_delta_matter = get_c200c_at_m_and_z(m,z,pba,pclass_sz);
  }
  else if (pclass_sz->delta_def_matter_density == 2){
    c_delta_matter = get_c500c_at_m_and_z(m,z,pba,pclass_sz);
  }
  else if (pclass_sz->delta_def_matter_density == 3){
    c_delta_matter = evaluate_cvir_of_mvir(m,z,pclass_sz,pba);
  }


struct Parameters_for_integrand_matter_density_profile_norm V;
V.pclass_sz = pclass_sz;
V.pba = pba;
V.z = z;
V.m = m;
V.c_delta = c_delta_matter;

void * params = &V;
double eps_rel = pclass_sz->matter_density_norm_epsrel;
double eps_abs = pclass_sz->matter_density_norm_epsabs;

// in for power law nfw this should be 2
// for truncated nfw, this should be c
// double x_out = pclass_sz->x_out_truncated_matter_density_profile_normalization;
double x_out = pclass_sz->x_out_matter_density_profile_normalization*c_delta_matter;//pclass_sz->x_out_truncated_matter_density_profile_normalization;

double x_in = 0.;



double r = Integrate_using_Patterson_adaptive(x_in, x_out,
                                              eps_rel, eps_abs,
                                              integrand_matter_density_profile_norm,
                                              params,
                                              pclass_sz->patterson_show_neval);

// printf("norm matter: index z = %d index_m = %d xout = %.3e r = %.5e m = %.5e z = %.5e index_z_m = %d\n",index_z,index_m,x_out,r,m,z, index_z_m);

// BCM: need to multiply by 4pi r200c^3:
// double r200c =  pow(m*3./4./_PI_/200./get_rho_crit_at_z(z,pba,pclass_sz),1./3.);
// double omega_b_over_omega_m = pclass_sz->f_b_gas;
// double fstar = get_fstar_of_m_at_z(m,z,pclass_sz);
// double num = omega_b_over_omega_m-fstar;
// double norm = 4.*_PI_*pow(r200c,3.)*r/m/num;
double norm =  r;
pclass_sz->array_ln_matter_density_norm_at_z_and_m[index_z_m] = log(norm);

if (index_z == 71 && index_m == 58)
printf("matter_density_profile_norm: index z = %d index_m = %d  z = %.4e m = %.4e lognorm = %.5e\n",
index_z,index_m,z,m,pclass_sz->array_ln_matter_density_norm_at_z_and_m[index_z_m]);

}
else{
  printf("normalization for this gas density profile not implemented yet.\n");
  exit(0);
}


    }
  }

  #ifdef _OPENMP
    tstop = omp_get_wtime();
    if (pclass_sz->sz_verbose > 0)
      printf("In %s: time spent in parallel region (loop over zm's) = %e s for thread %d\n",
             __func__,tstop-tstart,omp_get_thread_num());
  #endif


  // free(data);
  }


if (abort == _TRUE_) return _FAILURE_;

  return _SUCCESS_;
}



double get_dkappacmbdz_pklin_at_l_and_z(double l,
                                  double z,
                                  struct background * pba,
                                  struct primordial * ppm,
                                  struct nonlinear * pnl,
                                  struct class_sz_structure * pclass_sz){

  double tau;
  int first_index_back = 0;


  double * pvecback;
  double * pvectsz;
 class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);
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




      pvectsz[pclass_sz->index_z] = z;
      pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *pvecback[pba->index_bg_rho_crit]
                                            /pow(pba->h,2);

      double omega = pvecback[pba->index_bg_Omega_m];
      pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
      double chi = sqrt(pvectsz[pclass_sz->index_chi2]);
      double kl = (l+0.5)/chi;



      double pk1 = get_pk_lin_at_k_and_z(kl,z,pba,ppm,pnl,pclass_sz); // volume
      double result = pk1;
      double Wg = radial_kernel_W_cmb_lensing_at_z(z,pvectsz,pba,pclass_sz); // dimensionless
      result *= pow(Wg,2.);
      double Omega_m = pclass_sz->Omega_m_0;
      result *= pow(3.*pow(Omega_m,1.)*pow(pba->H0/pba->h,2)/2./chi*pow(1.+z,1.),2.); // volume^-2
      result *= get_volume_at_z(pvectsz[pclass_sz->index_z],pba); // volume


      free(pvecback);
      free(pvectsz);
      return result;
                                  }

double get_dkappacmbdz_at_l_and_z(double l,
                                  double z,
                                  struct background * pba,
                                  struct primordial * ppm,
                                  struct nonlinear * pnl,
                                  struct class_sz_structure * pclass_sz){

  double tau;
  int first_index_back = 0;


  double * pvecback;
  double * pvectsz;
 class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);
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




      pvectsz[pclass_sz->index_z] = z;
      pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *pvecback[pba->index_bg_rho_crit]
                                            /pow(pba->h,2);

      double omega = pvecback[pba->index_bg_Omega_m];
      pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
      double chi = sqrt(pvectsz[pclass_sz->index_chi2]);
      double kl = (l+0.5)/chi;



      double pk1 = get_pk_nonlin_at_k_and_z(kl,z,pba,ppm,pnl,pclass_sz); // volume
      double result = pk1;
      double Wg = radial_kernel_W_cmb_lensing_at_z(z,pvectsz,pba,pclass_sz); // dimensionless
      result *= pow(Wg,2.);
      double Omega_m = pclass_sz->Omega_m_0;
      result *= pow(3.*pow(Omega_m,1.)*pow(pba->H0/pba->h,2)/2./chi*pow(1.+z,1.),2.); // volume^-2
      result *= get_volume_at_z(pvectsz[pclass_sz->index_z],pba); // volume


      free(pvecback);
      free(pvectsz);
      return result;
                                  }

double get_dyldzdlnm_at_l_z_and_m(double l,
                                  double z,
                                  double m,
                                  struct background * pba,
                                  struct nonlinear * pnl,
                                  struct class_sz_structure * pclass_sz){


// double result = get_dndlnM_at_z_and_M(z_asked,m,pclass_sz)
//                 *get_volume_at_z(z,pba)
//                 *evaluate_pressure_profile(kl,pvecback,pvectsz,pba,pclass_sz);



  // double M_halo = m;

  double tau;
  int first_index_back = 0;


  double * pvecback;
  double * pvectsz;
 class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);
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




      pvectsz[pclass_sz->index_z] = z;
      pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *pvecback[pba->index_bg_rho_crit]
                                            /pow(pba->h,2);

      double omega = pvecback[pba->index_bg_Omega_m];
      pvectsz[pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);


      // request appropriate mass conversion
      pvectsz[pclass_sz->index_has_electron_pressure] = 1 ;

      do_mass_conversions(log(m),z,pvecback,pvectsz,pba,pclass_sz);
      evaluate_HMF_at_logM_and_z(log(m),z,pvecback,pvectsz,pba,pnl,pclass_sz);

      double hmf = pvectsz[pclass_sz->index_hmf];
      pvectsz[pclass_sz->index_md] = -1;//pclass_sz->index_md_dydz;


      double kl;
      if (l==0)
        kl = 0.;
      else
        kl = (l+0.5)/sqrt(pvectsz[pclass_sz->index_chi2]);

      evaluate_pressure_profile(kl,pvecback,pvectsz,pba,pclass_sz);


      double result = hmf*pvectsz[pclass_sz->index_pressure_profile];

      // multiply by volume element:
      double H_over_c_in_h_over_Mpc = pvecback[pba->index_bg_H]/pba->h;
      result *= pvectsz[pclass_sz->index_chi2]/H_over_c_in_h_over_Mpc;
      result *= 1./pow(pclass_sz->Tcmb_gNU,1)/1.e6;
      free(pvecback);
      free(pvectsz);

return result;
                                }

double get_normalization_gas_density_profile(double z_asked, double m_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);

  double result = 1e100;
  if (z<pclass_sz->array_profile_ln_1pz[0]){
    result = 1e100;
  }
  else if (z>pclass_sz->array_profile_ln_1pz[pclass_sz->n_z_density_profile-1]){
    result = 1e100;
  }
  else if (m<pclass_sz->array_profile_ln_m[0]){
    result = 1e100;
  }
  else if (m>pclass_sz->array_profile_ln_m[pclass_sz->n_m_density_profile-1]){
    result = 1e100;
  }
  else{
    result = exp(pwl_interp_2d(
                              pclass_sz->n_z_density_profile,
                               pclass_sz->n_m_density_profile,
                               pclass_sz->array_profile_ln_1pz,
                               pclass_sz->array_profile_ln_m,
                               pclass_sz->array_ln_density_norm_at_z_and_m,
                               1,
                               &z,
                               &m));

     // result *= 1/exp(m);

  // if (pclass_sz->tau_profile == 2){
  //   double omega_b_over_omega_m = pclass_sz->f_b_gas;
  //   double fstar = get_fstar_of_m_at_z(m_asked,z,pclass_sz);
  //   double num = omega_b_over_omega_m-fstar;
  //   result *= 1./num;
  // }

  // nfw case already normalized.
  // if (pclass_sz->tau_profile == 2){
  // do nothing
  // }

  }
  return result;
}


double get_normalization_matter_density_profile(double z_asked, double m_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);

  double result = 1e100;
  if (z<pclass_sz->array_matter_density_profile_ln_1pz[0]){
    result = 1e100;
  }
  else if (z>pclass_sz->array_matter_density_profile_ln_1pz[pclass_sz->n_z_matter_density_profile-1]){
    result = 1e100;
  }
  else if (m<pclass_sz->array_matter_density_profile_ln_m[0]){
    result = 1e100;
  }
  else if (m>pclass_sz->array_matter_density_profile_ln_m[pclass_sz->n_m_matter_density_profile-1]){
    result = 1e100;
  }
  else{
    result = exp(pwl_interp_2d(
                              pclass_sz->n_z_matter_density_profile,
                               pclass_sz->n_m_matter_density_profile,
                               pclass_sz->array_matter_density_profile_ln_1pz,
                               pclass_sz->array_matter_density_profile_ln_m,
                               pclass_sz->array_ln_matter_density_norm_at_z_and_m,
                               1,
                               &z,
                               &m));

     // result *= 1/exp(m);

  // if (pclass_sz->tau_profile == 2){
  //   double omega_b_over_omega_m = pclass_sz->f_b_gas;
  //   double fstar = get_fstar_of_m_at_z(m_asked,z,pclass_sz);
  //   double num = omega_b_over_omega_m-fstar;
  //   result *= 1./num;
  // }

  // nfw case already normalized.
  // if (pclass_sz->tau_profile == 2){
  // do nothing
  // }

  }
  return result;
}


double get_n5k_pk_at_z_and_k(double z_asked, double k_asked, struct class_sz_structure * pclass_sz){
  double z = z_asked;
  double k = log(k_asked);

  double result = 0.;
  if (z<pclass_sz->n5k_pk_z[0]){
    printf("z too small\n");
    result = 0.;
  }
  else if (z>pclass_sz->n5k_pk_z[pclass_sz->n5k_pk_z_size-1]){
    printf("z too big\n");
    result = 0.;
  }
  else if (k<pclass_sz->n5k_pk_k[0]){
    printf("k too small\n");
    result = 0.;
  }
  else if (k>pclass_sz->n5k_pk_k[pclass_sz->n5k_pk_k_size-1]){
    printf("k too big\n");
    result = 0.;
  }
  else{
    result = exp(pwl_interp_2d(pclass_sz->n5k_pk_k_size,
                               pclass_sz->n5k_pk_z_size,
                               pclass_sz->n5k_pk_k,
                               pclass_sz->n5k_pk_z,
                               pclass_sz->n5k_pk_pk,
                               1,
                               &k,
                               &z));
  }
  return result;
}




double get_cib_Snu_z_and_nu(double z_asked, double nu_asked, struct class_sz_structure * pclass_sz){
  double z = z_asked;
  double k = nu_asked;

  double result = 0.;
  if (z<pclass_sz->cib_Snu_z[0]){
    printf("z too small in Snu z = %.5e snuzmax = %.5e \n",z,pclass_sz->cib_Snu_z[0]);
    result = 0.;
  }
  else if (z>pclass_sz->cib_Snu_z[pclass_sz->cib_Snu_z_size-1]){
    printf("z too big in Snu z = %.5e snuzmax = %.5e \n",z,pclass_sz->cib_Snu_z[pclass_sz->cib_Snu_z_size-1]);
    result = 0.;
  }
  else if (k<pclass_sz->cib_Snu_nu[0]){
    printf("nu too small\n");
    result = 0.;
  }
  else if (k>pclass_sz->cib_Snu_nu[pclass_sz->cib_Snu_nu_size-1]){
    printf("nu too big nu = %.5e numax = %.3e\n",k,pclass_sz->cib_Snu_nu[pclass_sz->cib_Snu_nu_size-1]);
    result = 0.;
  }
  else{
    result = exp(pwl_interp_2d(pclass_sz->cib_Snu_nu_size,
                               pclass_sz->cib_Snu_z_size,
                               pclass_sz->cib_Snu_nu,
                               pclass_sz->cib_Snu_z,
                               pclass_sz->cib_Snu_snu,
                               1,
                               &k,
                               &z));

    if (isnan(result) || isinf(result)) {
      printf("Error in get_cib_Snu_z_and_nu: Non-finite result detected!\n");
      printf("Parameters: z = %.5e, nu = %.5e, result = %.5e\n", z, k, result);
      exit(1);
    }
  }
  return result;
}






// Tabulate 2D Fourier transform of density profile on a [z - ln_M - ln_ell] grid
// this is the tau profile for kSZ
int tabulate_gas_density_profile(struct background * pba,
                                struct class_sz_structure * pclass_sz){

if (pclass_sz->has_kSZ_kSZ_lensmag_1halo

+ pclass_sz->has_kSZ_kSZ_gal_1h_fft
+ pclass_sz->has_kSZ_kSZ_gal_2h_fft
+ pclass_sz->has_kSZ_kSZ_gal_3h_fft
+ pclass_sz->has_kSZ_kSZ_gal_1h
+ pclass_sz->has_kSZ_kSZ_gal_2h
+ pclass_sz->has_kSZ_kSZ_gal_3h
+ pclass_sz->has_kSZ_kSZ_tSZ_1h
+ pclass_sz->has_kSZ_kSZ_tSZ_2h
+ pclass_sz->has_tau_gal_1h
+ pclass_sz->has_tau_gal_2h
+ pclass_sz->has_tau_tau_1h
+ pclass_sz->has_tau_tau_2h
+ pclass_sz->has_kSZ_kSZ_1h
+ pclass_sz->has_kSZ_kSZ_2h
+ pclass_sz->has_pk_bb_at_z_1h
+ pclass_sz->has_pk_bb_at_z_2h
+ pclass_sz->has_pk_b_at_z_2h
+ pclass_sz->has_gas_density_profile_2h
+ pclass_sz->has_pk_em_at_z_1h
+ pclass_sz->has_pk_em_at_z_2h
+ pclass_sz->has_kSZ_kSZ_tSZ_3h
+ pclass_sz->has_bk_ttg_at_z_1h
+ pclass_sz->has_bk_ttg_at_z_2h
+ pclass_sz->has_bk_ttg_at_z_3h
+ pclass_sz->has_kSZ_kSZ_gallens_1h_fft
+ pclass_sz->has_kSZ_kSZ_gallens_2h_fft
+ pclass_sz->has_kSZ_kSZ_gallens_3h_fft
+ pclass_sz->has_kSZ_kSZ_lens_1h_fft
+ pclass_sz->has_kSZ_kSZ_lens_2h_fft
+ pclass_sz->has_kSZ_kSZ_lens_3h_fft
== _FALSE_
)
  return 0;


 // array of multipoles:

 double ln_ell_min = log(pclass_sz->k_min_gas_density_profile);
 double ln_ell_max = log(pclass_sz->k_max_gas_density_profile);
 int n_ell = pclass_sz->n_k_density_profile;
 int n_m = pclass_sz->n_m_density_profile;
 int n_z = pclass_sz->n_z_density_profile;

 class_alloc(pclass_sz->array_profile_ln_k,sizeof(double *)*n_ell,pclass_sz->error_message);

 // array of masses:
 double ln_m_min = log(5e8);
 double ln_m_max = log(1e16);


 class_alloc(pclass_sz->array_profile_ln_m,sizeof(double *)*n_m,pclass_sz->error_message);


 // array of redshifts:
 double ln_1pz_min = log(1.+pclass_sz->z1SZ);
 double ln_1pz_max = log(1.+pclass_sz->z2SZ);


 class_alloc(pclass_sz->array_profile_ln_1pz,sizeof(double *)*n_z,pclass_sz->error_message);
int index_m_z;

int index_l;
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
  // this is k
  pclass_sz->array_profile_ln_k[index_l] = ln_ell_min
                                      +index_l*(ln_ell_max-ln_ell_min)
                                      /(n_ell-1.);
}

int index_m;
for (index_m=0;
     index_m<n_m;
     index_m++)
{
  pclass_sz->array_profile_ln_m[index_m] = ln_m_min
                                      +index_m*(ln_m_max-ln_m_min)
                                      /(n_m-1.);
}

int index_z;
for (index_z=0;
     index_z<n_z;
     index_z++)
{
  pclass_sz->array_profile_ln_1pz[index_z] = ln_1pz_min
                                   +index_z*(ln_1pz_max-ln_1pz_min)
                                   /(n_z-1.);
}


class_alloc(pclass_sz->array_profile_ln_rho_at_lnk_lnM_z,n_ell*sizeof(double *),pclass_sz->error_message);
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
class_alloc(pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l],n_m*n_z*sizeof(double *),pclass_sz->error_message);
index_m_z = 0;
for (index_m=0;
     index_m<n_m;
     index_m++){

for (index_z=0;
     index_z<n_z;
     index_z++)
{
  pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z] = log(1e-100); // initialize with super small number
  index_m_z += 1;
}

     }
}

int has_ksz_bkp = pclass_sz->has_kSZ_kSZ_gal_1h;
pclass_sz->has_kSZ_kSZ_gal_1h = _TRUE_; //pretend we need the tau_profile

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
pclass_sz,pba)\
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
class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
int index_pvectsz;
for (index_pvectsz=0;
     index_pvectsz<pclass_sz->tsz_size;
     index_pvectsz++){
       pvectsz[index_pvectsz] = 0.; // set everything to 0.
     }
index_m_z = 0;
for (index_m=0;
     index_m<n_m;
     index_m++){
for (index_z=0;
     index_z<n_z;
     index_z++){




  double z = exp(pclass_sz->array_profile_ln_1pz[index_z])-1.;
  double lnM = pclass_sz->array_profile_ln_m[index_m];
  double ell = exp(pclass_sz->array_profile_ln_k[index_l]);


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
  pvectsz[pclass_sz->index_z] = z;



  pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in rho_nfw for the pk mode computation

  pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);
  pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
  double omega = pvecback[pba->index_bg_Omega_m];
  pvectsz[pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);



   double r_delta;// = pvectsz[pclass_sz->index_radius_for_electron_density];
   double c_delta;// = pvectsz[pclass_sz->index_concentration_for_electron_density];
   double m_delta;// = pvectsz[pclass_sz->index_mass_for_electron_density];
   // printf("de = %d\n",pclass_sz->delta_def_electron_density);
   // exit(0);

  if (pclass_sz->delta_def_electron_density == 0){
    m_delta = exp(lnM);
    r_delta = pow(3.*m_delta/(4.*_PI_*200.*pvecback[pba->index_bg_Omega_m]*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
    c_delta = get_c200m_at_m_and_z(m_delta,z,pba,pclass_sz);
  }
  else if (pclass_sz->delta_def_electron_density == 1){
    m_delta = exp(lnM);
    r_delta = pow(3.*m_delta/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
    c_delta = get_c200c_at_m_and_z(m_delta,z,pba,pclass_sz);
  }
  else if (pclass_sz->delta_def_electron_density == 2){
    m_delta = exp(lnM);
    r_delta = pow(3.*m_delta/(4.*_PI_*500.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
    c_delta = get_c500c_at_m_and_z(m_delta,z,pba,pclass_sz);
  }

  ell = ell/(1.+z)/r_delta;

  pvectsz[pclass_sz->index_multipole_for_nfw_profile] = ell;

  double result;


  // pvectsz[pclass_sz->index_has_electron_density] = 1;
  // do_mass_conversions(lnM,z,pvecback,pvectsz,pba,pclass_sz);

 // only  do the integration of Battaglia and BCM models
 // nfw has an analytical formula
 if (pclass_sz->tau_profile == 1 || pclass_sz->tau_profile == 2){
  pvectsz[pclass_sz->index_m200c] = exp(lnM);
  // class_call_parallel(mDEL_to_mVIR(pvectsz[pclass_sz->index_m200c],
  //                                  200.*(pvectsz[pclass_sz->index_Rho_crit]),
  //                                  pvectsz[pclass_sz->index_Delta_c],
  //                                  pvectsz[pclass_sz->index_Rho_crit],
  //                                  z,
  //                                  &pvectsz[pclass_sz->index_mVIR],
  //                                  pclass_sz,
  //                                  pba),
  //                 pclass_sz->error_message,
  //                 pclass_sz->error_message);
 //
 //  // rvir needed to cut off the integral --> e.g., xout = 50.*rvir/r200c
  // pvectsz[pclass_sz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[pclass_sz->index_mVIR],pvectsz[pclass_sz->index_Delta_c],pvectsz[pclass_sz->index_Rho_crit],pclass_sz);

  pvectsz[pclass_sz->index_r200c] = pow(3.*pvectsz[pclass_sz->index_m200c]/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
  pvectsz[pclass_sz->index_l200c] = sqrt(pvectsz[pclass_sz->index_chi2])/(1.+z)/pvectsz[pclass_sz->index_r200c];
  pvectsz[pclass_sz->index_characteristic_multipole_for_nfw_profile] = pvectsz[pclass_sz->index_l200c];
  // pvectsz[pclass_sz->index_rs] = pvectsz[pclass_sz->index_r200c];

 double result_int;


two_dim_ft_nfw_profile(pclass_sz,pba,pvectsz,&result_int);


 result = result_int;
 double tau_normalisation = 1.;
 tau_normalisation = 4.*_PI_*pow(pvectsz[pclass_sz->index_r200c],3);
 // printf("In tab gas: k %.4e z %.8e rt %.8e mt %.8e res = %.4e\n",ell,pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_r200c],pvectsz[pclass_sz->index_m200c],result);


 if (result<=0 || isnan(result) || isinf(result)){
 // printf("ERROR: In tab gas: k %.4e z %.8e rt %.8e mt %.8e res = %.4e\n",ell,pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_r200c],pvectsz[pclass_sz->index_m200c],result);
 // printf("check precision and input parameters?\n");
 // exit(0);
 result = 1e-200;
}

 result *= tau_normalisation;
}
else if (pclass_sz->tau_profile == 0){ // truncated nfw profile
   //
   // pvectsz[pclass_sz->index_m200c] = exp(lnM);
   // pvectsz[pclass_sz->index_r200c] = pow(3.*pvectsz[pclass_sz->index_m200c]/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
   // pvectsz[pclass_sz->index_c200c] = get_c200c_at_m_and_z(pvectsz[pclass_sz->index_m200c],z,pba,pclass_sz);
   //
   // double r_delta = pvectsz[pclass_sz->index_r200c];
   // double c_delta = pvectsz[pclass_sz->index_c200c];
   // double m_delta = pvectsz[pclass_sz->index_m200c];


   double xout = pclass_sz->x_out_truncated_nfw_profile_electrons;



   // pvectsz[pclass_sz->index_rs] = r_delta/c_delta;

  // pvectsz[pclass_sz->index_multipole_for_truncated_nfw_profile] = pvectsz[pclass_sz->index_multipole_for_nfw_profile];
  // double l = pvectsz[pclass_sz->index_multipole_for_truncated_nfw_profile];
  double chi = sqrt(pvectsz[pclass_sz->index_chi2]);
  double k = ell;//*(1.+z)*r_delta;
   result =  evaluate_truncated_nfw_profile(pvectsz[pclass_sz->index_z],k,r_delta,c_delta,xout);
   //result *= 1.;//m_delta;///(4.*_PI_*pow(pvectsz[pclass_sz->index_rs],3));
   double f_b = pclass_sz->f_b_gas;//pba->Omega0_b/pclass_sz->Omega_m_0;
   result *= f_b*m_delta;//*pow((pba->Omega0_cdm+pba->Omega0_b)*pclass_sz->Rho_crit_0,-1);

    if (isnan(result) || isinf(result)){
    printf("z %.8e rt %.8e ct %.8e mt %.8e\n",pvectsz[pclass_sz->index_z],r_delta,c_delta,m_delta);
    exit(0);
  }
   // double tau_normalisation = 4.*_PI_*pow(pvectsz[pclass_sz->index_rs],3);
   // result *= tau_normalisation;
 }

  pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z] = log(result);
  // printf("l = %.8e m = %.8e z = %.8e lnrho = %.8e\n",ell,exp(lnM),z,log(result));

  index_m_z += 1;
     }


     }

     free(pvectsz);
     free(pvecback);
}

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in tab profile parallel region (loop over ell's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

// restore initial state:
pclass_sz->has_kSZ_kSZ_gal_1h = has_ksz_bkp;


                                      }




//
// // Tabulate 2D Fourier transform of density profile on a [z - ln_M - ln_ell] grid
// // this is the tau profile for kSZ
// int tabulate_gas_density_profile_fft(struct background * pba,
//                                 struct class_sz_structure * pclass_sz){
//
// if (pclass_sz->has_kSZ_kSZ_lensmag_1halo
//
// + pclass_sz->has_kSZ_kSZ_gal_1h_fft
// + pclass_sz->has_kSZ_kSZ_gal_2h_fft
// + pclass_sz->has_kSZ_kSZ_gal_3h_fft
// + pclass_sz->has_kSZ_kSZ_gal_1h
// + pclass_sz->has_kSZ_kSZ_gal_2h
// + pclass_sz->has_kSZ_kSZ_gal_3h
// + pclass_sz->has_kSZ_kSZ_tSZ_1h
// + pclass_sz->has_kSZ_kSZ_tSZ_2h
// + pclass_sz->has_kSZ_kSZ_1h
// + pclass_sz->has_kSZ_kSZ_2h
// + pclass_sz->has_pk_bb_at_z_1h
// + pclass_sz->has_pk_bb_at_z_2h
// + pclass_sz->has_pk_b_at_z_2h
// + pclass_sz->has_pk_em_at_z_1h
// + pclass_sz->has_pk_em_at_z_2h
// + pclass_sz->has_kSZ_kSZ_tSZ_3h
// + pclass_sz->has_bk_ttg_at_z_1h
// + pclass_sz->has_bk_ttg_at_z_2h
// + pclass_sz->has_bk_ttg_at_z_3h
// + pclass_sz->has_kSZ_kSZ_gallens_1h_fft
// + pclass_sz->has_kSZ_kSZ_gallens_2h_fft
// + pclass_sz->has_kSZ_kSZ_gallens_3h_fft
// + pclass_sz->has_kSZ_kSZ_lens_1h_fft
// + pclass_sz->has_kSZ_kSZ_lens_2h_fft
// + pclass_sz->has_kSZ_kSZ_lens_3h_fft
// == _FALSE_
// )
//   return 0;
//
//
//  // array of multipoles:
//
//  double ln_ell_min = log(pclass_sz->k_min_gas_density_profile);
//  double ln_ell_max = log(pclass_sz->k_max_gas_density_profile);
//  int n_ell = pclass_sz->n_k_density_profile;
//  int n_m = pclass_sz->n_m_density_profile;
//  int n_z = pclass_sz->n_z_density_profile;
//
//  class_alloc(pclass_sz->array_profile_ln_k,sizeof(double *)*n_ell,pclass_sz->error_message);
//
//  // array of masses:
//  double ln_m_min = log(5e8);
//  double ln_m_max = log(1e16);
//
//
//  class_alloc(pclass_sz->array_profile_ln_m,sizeof(double *)*n_m,pclass_sz->error_message);
//
//
//  // array of redshifts:
//  double ln_1pz_min = log(1.+pclass_sz->z1SZ);
//  double ln_1pz_max = log(1.+pclass_sz->z2SZ);
//
//
//  class_alloc(pclass_sz->array_profile_ln_1pz,sizeof(double *)*n_z,pclass_sz->error_message);
// int index_m_z;
// int index_m_z_tab[n_m][n_z];
// int index_l;
// for (index_l=0;
//      index_l<n_ell;
//      index_l++)
// {
//   // this is k
//   pclass_sz->array_profile_ln_k[index_l] = ln_ell_min
//                                       +index_l*(ln_ell_max-ln_ell_min)
//                                       /(n_ell-1.);
// }
//
// int index_m;
// for (index_m=0;
//      index_m<n_m;
//      index_m++)
// {
//   pclass_sz->array_profile_ln_m[index_m] = ln_m_min
//                                       +index_m*(ln_m_max-ln_m_min)
//                                       /(n_m-1.);
//
// }
//
// int index_z;
// for (index_z=0;
//      index_z<n_z;
//      index_z++)
// {
//   pclass_sz->array_profile_ln_1pz[index_z] = ln_1pz_min
//                                    +index_z*(ln_1pz_max-ln_1pz_min)
//                                    /(n_z-1.);
// }
//
//
// class_alloc(pclass_sz->array_profile_ln_rho_at_lnk_lnM_z,n_ell*sizeof(double *),pclass_sz->error_message);
// for (index_l=0;
//      index_l<n_ell;
//      index_l++)
// {
// class_alloc(pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l],n_m*n_z*sizeof(double *),pclass_sz->error_message);
// index_m_z = 0;
// for (index_m=0;
//      index_m<n_m;
//      index_m++){
//
// for (index_z=0;
//      index_z<n_z;
//      index_z++)
// {
//   pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z] = log(1e-100); // initialize with super small number
//
//   index_m_z_tab[index_m][index_z] = index_m_z;
//   index_m_z += 1;
// }
//
//      }
// }
//
// int has_ksz_bkp = pclass_sz->has_kSZ_kSZ_gal_1h;
// pclass_sz->has_kSZ_kSZ_gal_1h = _TRUE_; //pretend we need the tau_profile
//
// //Parallelization of profile computation
// /* initialize error management flag */
//
//
// int abort;
// double tstart, tstop;
// abort = _FALSE_;
// /* beginning of parallel region */
//
// int number_of_threads= 1;
// #ifdef _OPENMP
// #pragma omp parallel
//   {
//     number_of_threads = omp_get_num_threads();
//   }
// #endif
//
// #pragma omp parallel \
// shared(abort,\
// pclass_sz,pba,index_m_z_tab)\
// private(tstart, tstop,index_l,index_z,index_m,index_m_z) \
// num_threads(number_of_threads)
// {
//
// #ifdef _OPENMP
//   tstart = omp_get_wtime();
// #endif
//
// #pragma omp for schedule (dynamic)
// for (index_z=0;
//      index_z<n_z;
//      index_z++){
// #pragma omp flush(abort)
// double * pvectsz;
// double * pvecback;
//
// class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);
// class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
// int index_pvectsz;
// for (index_pvectsz=0;
//      index_pvectsz<pclass_sz->tsz_size;
//      index_pvectsz++){
//        pvectsz[index_pvectsz] = 0.; // set everything to 0.
//      }
// // index_m_z = 0;
// for (index_m=0;
//      index_m<n_m;
//      index_m++){
//
//
// // index_m_z = index_z + index_m*n_m;
// index_m_z  = index_m_z_tab[index_m][index_z];
//        for (index_l=0;
//             index_l<n_ell;
//             index_l++)
//        {
//   double z = exp(pclass_sz->array_profile_ln_1pz[index_z])-1.;
//   double lnM = pclass_sz->array_profile_ln_m[index_m];
//   double ell = exp(pclass_sz->array_profile_ln_k[index_l]);
//
//
//   double tau;
//   int first_index_back = 0;
//
//
//   class_call_parallel(background_tau_of_z(pba,z,&tau),
//              pba->error_message,
//              pba->error_message);
//
//   class_call_parallel(background_at_tau(pba,
//                                tau,
//                                pba->long_info,
//                                pba->inter_normal,
//                                &first_index_back,
//                                pvecback),
//              pba->error_message,
//              pba->error_message);
//
//
//   // fill relevant entries
//   pvectsz[pclass_sz->index_z] = z;
//
//   pvectsz[pclass_sz->index_multipole_for_nfw_profile] = ell;
//   pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in rho_nfw for the pk mode computation
//
//   pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
//                                 *pow(_Mpc_over_m_,1)
//                                 *pow(_c_,2)
//                                 *pvecback[pba->index_bg_rho_crit]
//                                 /pow(pba->h,2);
//   pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
//   double omega = pvecback[pba->index_bg_Omega_m];
//   pvectsz[pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
//
//
//   double result;
//
//
//   // pvectsz[pclass_sz->index_has_electron_density] = 1;
//   // do_mass_conversions(lnM,z,pvecback,pvectsz,pba,pclass_sz);
//
//  // only  do the integration of Battaglia profile
//  // nfw has an analytical formula
//  if (pclass_sz->tau_profile == 1){
//   pvectsz[pclass_sz->index_m200c] = exp(lnM);
//   // class_call_parallel(mDEL_to_mVIR(pvectsz[pclass_sz->index_m200c],
//   //                                  200.*(pvectsz[pclass_sz->index_Rho_crit]),
//   //                                  pvectsz[pclass_sz->index_Delta_c],
//   //                                  pvectsz[pclass_sz->index_Rho_crit],
//   //                                  z,
//   //                                  &pvectsz[pclass_sz->index_mVIR],
//   //                                  pclass_sz,
//   //                                  pba),
//   //                 pclass_sz->error_message,
//   //                 pclass_sz->error_message);
//  //
//  //  // rvir needed to cut off the integral --> e.g., xout = 50.*rvir/r200c
//   // pvectsz[pclass_sz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[pclass_sz->index_mVIR],pvectsz[pclass_sz->index_Delta_c],pvectsz[pclass_sz->index_Rho_crit],pclass_sz);
//
//   pvectsz[pclass_sz->index_r200c] = pow(3.*pvectsz[pclass_sz->index_m200c]/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
//   pvectsz[pclass_sz->index_l200c] = sqrt(pvectsz[pclass_sz->index_chi2])/(1.+z)/pvectsz[pclass_sz->index_r200c];
//   pvectsz[pclass_sz->index_characteristic_multipole_for_nfw_profile] = pvectsz[pclass_sz->index_l200c];
//   // pvectsz[pclass_sz->index_rs] = pvectsz[pclass_sz->index_r200c];
//
//  double result_int;
//
//
// two_dim_ft_nfw_profile(pclass_sz,pba,pvectsz,&result_int);
//
//
//  result = result_int;
//  double tau_normalisation = 1.;
//  tau_normalisation = 4.*_PI_*pow(pvectsz[pclass_sz->index_r200c],3);
//  // printf("In tab gas: k %.4e z %.8e rt %.8e mt %.8e res = %.4e\n",ell,pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_r200c],pvectsz[pclass_sz->index_m200c],result);
//
//
//  if (result<=0 || isnan(result) || isinf(result)){
//  // printf("ERROR: In tab gas: k %.4e z %.8e rt %.8e mt %.8e res = %.4e\n",ell,pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_r200c],pvectsz[pclass_sz->index_m200c],result);
//  // printf("check precision and input parameters?\n");
//  // exit(0);
//  result = 1e-200;
// }
//
//  result *= tau_normalisation;
// }
// else if (pclass_sz->tau_profile == 0){ // truncated nfw profile
//    //
//    // pvectsz[pclass_sz->index_m200c] = exp(lnM);
//    // pvectsz[pclass_sz->index_r200c] = pow(3.*pvectsz[pclass_sz->index_m200c]/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
//    // pvectsz[pclass_sz->index_c200c] = get_c200c_at_m_and_z(pvectsz[pclass_sz->index_m200c],z,pba,pclass_sz);
//    //
//    // double r_delta = pvectsz[pclass_sz->index_r200c];
//    // double c_delta = pvectsz[pclass_sz->index_c200c];
//    // double m_delta = pvectsz[pclass_sz->index_m200c];
//
//    double r_delta;// = pvectsz[pclass_sz->index_radius_for_electron_density];
//    double c_delta;// = pvectsz[pclass_sz->index_concentration_for_electron_density];
//    double m_delta;// = pvectsz[pclass_sz->index_mass_for_electron_density];
//    // printf("de = %d\n",pclass_sz->delta_def_electron_density);
//    // exit(0);
//
//   if (pclass_sz->delta_def_electron_density == 0){
//     m_delta = exp(lnM);
//     r_delta = pow(3.*m_delta/(4.*_PI_*200.*pvecback[pba->index_bg_Omega_m]*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
//     c_delta = get_c200m_at_m_and_z(m_delta,z,pba,pclass_sz);
//   }
//   else if (pclass_sz->delta_def_electron_density == 1){
//     m_delta = exp(lnM);
//     r_delta = pow(3.*m_delta/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
//     c_delta = get_c200c_at_m_and_z(m_delta,z,pba,pclass_sz);
//   }
//   else if (pclass_sz->delta_def_electron_density == 2){
//     m_delta = exp(lnM);
//     r_delta = pow(3.*m_delta/(4.*_PI_*500.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
//     c_delta = get_c500c_at_m_and_z(m_delta,z,pba,pclass_sz);
//   }
//    double xout = pclass_sz->x_out_truncated_nfw_profile_electrons;
//
//
//
//    // pvectsz[pclass_sz->index_rs] = r_delta/c_delta;
//
//   // pvectsz[pclass_sz->index_multipole_for_truncated_nfw_profile] = pvectsz[pclass_sz->index_multipole_for_nfw_profile];
//   // double l = pvectsz[pclass_sz->index_multipole_for_truncated_nfw_profile];
//   double chi = sqrt(pvectsz[pclass_sz->index_chi2]);
//   double k = ell;
//    result =  evaluate_truncated_nfw_profile(pvectsz[pclass_sz->index_z],k,r_delta,c_delta,xout);
//    //result *= 1.;//m_delta;///(4.*_PI_*pow(pvectsz[pclass_sz->index_rs],3));
//    double f_b = pclass_sz->f_b_gas;//pba->Omega0_b/pclass_sz->Omega_m_0;
//    result *= f_b*m_delta;//*pow((pba->Omega0_cdm+pba->Omega0_b)*pclass_sz->Rho_crit_0,-1);
//
//     if (isnan(result) || isinf(result)){
//     printf("z %.8e rt %.8e ct %.8e mt %.8e\n",pvectsz[pclass_sz->index_z],r_delta,c_delta,m_delta);
//     exit(0);
//   }
//    // double tau_normalisation = 4.*_PI_*pow(pvectsz[pclass_sz->index_rs],3);
//    // result *= tau_normalisation;
//  }
//
//   pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z] = log(result);
//   // printf("l = %.8e m = %.8e z = %.8e lnrho = %.8e\n",ell,exp(lnM),z,log(result));
//
//   // index_m_z += 1;
//      }
//
//
//      }
//
//      free(pvectsz);
//      free(pvecback);
// }
//
// #ifdef _OPENMP
//   tstop = omp_get_wtime();
//   if (pclass_sz->sz_verbose > 0)
//     printf("In %s: time spent in tab profile parallel region (loop over ell's) = %e s for thread %d\n",
//            __func__,tstop-tstart,omp_get_thread_num());
// #endif
// }
// if (abort == _TRUE_) return _FAILURE_;
// //end of parallel region
//
// // restore initial state:
// pclass_sz->has_kSZ_kSZ_gal_1h = has_ksz_bkp;
//
//
//                                       }
//
//
//

int tabulate_matter_nfw_with_power_law_profile_fft(struct background * pba,
                                                   struct class_sz_structure * pclass_sz){

 int n_k = pclass_sz->N_samp_fftw;
 int n_m = pclass_sz->n_m_matter_density_profile;
 int n_z = pclass_sz->n_z_matter_density_profile;
 int index_m_z_tab[n_m][n_z];
 // array of masses:
 // double ln_m_min = log(pclass_sz->M1SZ);
 // double ln_m_max = log(pclass_sz->M2SZ);
 //
 // // array of redshifts:
 // double ln_1pz_min = log(1.+pclass_sz->z1SZ);
 // double ln_1pz_max = log(1.+pclass_sz->z2SZ);

 class_alloc(pclass_sz->array_matter_density_profile_ln_k,sizeof(double *)*n_k,pclass_sz->error_message);
 // class_alloc(pclass_sz->array_matter_density_profile_ln_m,sizeof(double *)*n_m,pclass_sz->error_message);
 // class_alloc(pclass_sz->array_matter_density_profile_ln_1pz,sizeof(double *)*n_z,pclass_sz->error_message);


 class_alloc(pclass_sz->array_profile_ln_rho_matter_at_lnk,n_k*sizeof(double *),pclass_sz->error_message);

int index_m;
// for (index_m=0;
//      index_m<n_m;
//      index_m++)
// {
//   pclass_sz->array_matter_density_profile_ln_m[index_m] = ln_m_min
//                                       +index_m*(ln_m_max-ln_m_min)
//                                       /(n_m-1.);
// }

int index_z;
// for (index_z=0;
//      index_z<n_z;
//      index_z++)
// {
//   pclass_sz->array_matter_density_profile_ln_1pz[index_z] = ln_1pz_min
//                                                        +index_z*(ln_1pz_max-ln_1pz_min)
//                                                        /(n_z-1.);
// }


int index_m_z;
int index_k;

for (index_k=0;
     index_k<n_k;
     index_k++)
{
class_alloc(pclass_sz->array_profile_ln_rho_matter_at_lnk[index_k],n_m*n_z*sizeof(double *),pclass_sz->error_message);
index_m_z = 0;
for (index_z=0;
     index_z<n_z;
     index_z++)
{

for (index_m=0;
     index_m<n_m;
     index_m++){
//class_alloc(pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z],n_z*sizeof(double ),pclass_sz->error_message);


  // pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z] = -100.; // initialize with super small number
  pclass_sz->array_profile_ln_rho_matter_at_lnk[index_k][index_m_z] = 1e-100; // initialize with super small number
  index_m_z_tab[index_m][index_z] = index_m_z;
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
pclass_sz,pba,index_m_z_tab)\
private(tstart, tstop,index_k,index_z,index_m,index_m_z) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
for (index_z=0;
     index_z<n_z;
     index_z++){
#pragma omp flush(abort)
double * pvectsz;
double * pvecback;
class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);
class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);


int index_pvectsz;
for (index_pvectsz=0;
     index_pvectsz<pclass_sz->tsz_size;
     index_pvectsz++){
       pvectsz[index_pvectsz] = 0.; // set everything to 0.
     }

double z = exp(pclass_sz->array_matter_density_profile_ln_1pz[index_z])-1.;
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
  pvectsz[pclass_sz->index_z] = z;


  pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
  double chi = sqrt(pvectsz[pclass_sz->index_chi2]);
  // pvectsz[pclass_sz->index_multipole_for_pressure_profile] = k*chi;
  pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in p_gnfw for the pk mode computation

  pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);

  double omega = pvecback[pba->index_bg_Omega_m];
  pvectsz[pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);



for (index_m=0;
     index_m<n_m;
     index_m++){
  double lnM = pclass_sz->array_matter_density_profile_ln_m[index_m];
  // index_m_z = index_m+ index_z*n_z;
  index_m_z  = index_m_z_tab[index_m][index_z];



// here we FFT the profile ===== commented
const int N = pclass_sz->N_samp_fftw; //precision parameter
int ix;
double x[N], Px[N];
double x_min = pclass_sz->x_min_matter_density_fftw;
double x_max = pclass_sz->x_max_matter_density_fftw;

double m_delta = exp(lnM);
double c_delta_matter;
double r_delta_matter;
  if (pclass_sz->delta_def_matter_density == 0){
    c_delta_matter = get_c200m_at_m_and_z(m_delta,z,pba,pclass_sz);
    r_delta_matter = pow(3.*m_delta/(4.*_PI_*200.*pvecback[pba->index_bg_Omega_m]*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
  }
  else if (pclass_sz->delta_def_matter_density == 1){
    c_delta_matter = get_c200c_at_m_and_z(m_delta,z,pba,pclass_sz);
    r_delta_matter = pow(3.*m_delta/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);

  }
  else if (pclass_sz->delta_def_matter_density == 2){
    c_delta_matter = get_c500c_at_m_and_z(m_delta,z,pba,pclass_sz);
    r_delta_matter = pow(3.*m_delta/(4.*_PI_*500.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
  }
  else if (pclass_sz->delta_def_matter_density == 3){
    c_delta_matter = evaluate_cvir_of_mvir(m_delta,z,pclass_sz,pba);
    r_delta_matter = evaluate_rvir_of_mvir(m_delta,pvectsz[pclass_sz->index_Delta_c],pvectsz[pclass_sz->index_Rho_crit],pclass_sz);
  }

for (ix=0; ix<N; ix++){
   x[ix] = exp(log(x_min)+ix/(N-1.)*(log(x_max)-log(x_min)));
  // if (x[ix]>x_out){
  //     Px[ix] = 0.;
  //   }
  // else{
    Px[ix] =  get_nfw_with_power_law_profile_at_x(x[ix],
                                                  pclass_sz->matter_nfw_power_law_index,
                                                  pclass_sz->x_out_matter_density_profile*c_delta_matter);
                                              // }
  // printf("x = %.3e Px = %.3e\n",x[ix],Px[ix]);
  }

  double kp[N], Pkp[N];
    // pk2xi(N,k,Pk1,rp,xi1,pclass_sz);
  /* Compute the function
   *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
   * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
   * in this notation.  The input k-values must be logarithmically spaced.  The
   * resulting xi_l^m(r) will be evaluated at the dual r-values
   *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
  //void fftlog_ComputeXiLM(int l, int m, int N, const double k[],  const double pk[], double r[], double xi[]);
  fftlog_ComputeXiLMsloz(0, 2, N, x, Px, kp, Pkp,pclass_sz);

/// ---- > FFT end commented

  double result;
  int index_k;
//
// if (index_z == 10 && index_m == 23){
//
// char Filepath[_ARGUMENT_LENGTH_MAX_];
//
// FILE *fp;
// sprintf(Filepath,"%s%s%s",pclass_sz->root,"","test_nfw.txt");
// fp=fopen(Filepath, "w");
//
//   for (index_k=0;
//        index_k<n_k;
//        index_k++)
//   {
//
//     double k;
//     double  result_fft;
//
//     pclass_sz->array_matter_density_profile_ln_k[index_k] = log(kp[index_k]);//(pvectsz[pclass_sz->index_r200c]*(1.+pvectsz[pclass_sz->index_z])));
//     result_fft = 2.*_PI_*_PI_*Pkp[index_k];
//
//
//     double result = result_fft;
//
//
//    if (result<=0 || isnan(result) || isinf(result)){
//          // printf("ERROR: In tab gas: k %.4e z %.8e rt %.8e mt %.8e res = %.4e\n",ell,pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_r200c],pvectsz[pclass_sz->index_m200c],result);
//          // printf("check precision and input parameters?\n");
//          // exit(0);
//          result = 1e-200;
//         }
//     double result_trunc = evaluate_truncated_nfw_profile(z,
//                                                          kp[index_k]/(r_delta_matter/c_delta_matter*(1.+z)),
//                                                          r_delta_matter,
//                                                          c_delta_matter,
//                                                          pclass_sz->x_out_truncated_nfw_profile);
//     printf("matter nfw: k = %.5e r = %.5e ratio = %.5e\n",kp[index_k],result,result/result_trunc);
//
//     double norm = get_normalization_matter_density_profile(z,m_delta,pclass_sz);
// fprintf(fp,"%.5e \t %.5e \t %.5e\n",kp[index_k],result/norm,result_trunc);
//
//    pclass_sz->array_profile_ln_rho_matter_at_lnk[index_k][index_m_z] = log(result);
//  } // k loop
// fclose(fp);
//
// exit(0);
// }
// else{


  for (index_k=0;
       index_k<n_k;
       index_k++)
  {

    double k;
    double  result_fft;

    pclass_sz->array_matter_density_profile_ln_k[index_k] = log(kp[index_k]);//(pvectsz[pclass_sz->index_r200c]*(1.+pvectsz[pclass_sz->index_z])));
    result_fft = 2.*_PI_*_PI_*Pkp[index_k];


    double result = result_fft;


   if (result<=0 || isnan(result) || isinf(result)){
         // printf("ERROR: In tab gas: k %.4e z %.8e rt %.8e mt %.8e res = %.4e\n",ell,pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_r200c],pvectsz[pclass_sz->index_m200c],result);
         // printf("check precision and input parameters?\n");
         // exit(0);
         result = 1e-200;
        }
    // double result_trunc = evaluate_truncated_nfw_profile(0.,
    //                                                      kp[index_k],
    //                                                      1.,
    //                                                      1.,
    //                                                      1.)*m_nfw(1.);


   pclass_sz->array_profile_ln_rho_matter_at_lnk[index_k][index_m_z] = log(result);
 } // k loop

// } // else condition


} // m loop

free(pvectsz);
free(pvecback);

} // z loop

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in tab pressure profile parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

}
if (abort == _TRUE_) return _FAILURE_;

return _SUCCESS_;
                                     }

// Tabulate 2D Fourier transform of density profile on a [z - ln_M - ln_ell] grid
// this is the tau profile for kSZ
int tabulate_gas_density_profile_fft(struct background * pba,
                                     struct class_sz_structure * pclass_sz){

if (pclass_sz->has_kSZ_kSZ_lensmag_1halo

+ pclass_sz->has_kSZ_kSZ_gal_1h_fft
+ pclass_sz->has_kSZ_kSZ_gal_2h_fft
+ pclass_sz->has_kSZ_kSZ_gal_3h_fft
+ pclass_sz->has_kSZ_kSZ_gal_1h
+ pclass_sz->has_kSZ_kSZ_gal_2h
+ pclass_sz->has_kSZ_kSZ_gal_3h
+ pclass_sz->has_kSZ_kSZ_tSZ_1h
+ pclass_sz->has_kSZ_kSZ_tSZ_2h
+ pclass_sz->has_tau_gal_1h
+ pclass_sz->has_tau_gal_2h
+ pclass_sz->has_tau_tau_1h
+ pclass_sz->has_tau_tau_2h
+ pclass_sz->has_kSZ_kSZ_1h
+ pclass_sz->has_kSZ_kSZ_2h
+ pclass_sz->has_pk_bb_at_z_1h
+ pclass_sz->has_pk_bb_at_z_2h
+ pclass_sz->has_pk_b_at_z_2h
+ pclass_sz->has_gas_density_profile_2h
+ pclass_sz->has_pk_em_at_z_1h
+ pclass_sz->has_pk_em_at_z_2h
+ pclass_sz->has_kSZ_kSZ_tSZ_3h
+ pclass_sz->has_bk_ttg_at_z_1h
+ pclass_sz->has_bk_ttg_at_z_2h
+ pclass_sz->has_bk_ttg_at_z_3h
+ pclass_sz->has_kSZ_kSZ_gallens_1h_fft
+ pclass_sz->has_kSZ_kSZ_gallens_2h_fft
+ pclass_sz->has_kSZ_kSZ_gallens_3h_fft
+ pclass_sz->has_kSZ_kSZ_lens_1h_fft
+ pclass_sz->has_kSZ_kSZ_lens_2h_fft
+ pclass_sz->has_kSZ_kSZ_lens_3h_fft
== _FALSE_
)
  return 0;


 // array of multipoles:

 // double ln_k_min = log(pclass_sz->k_min_gas_density_profile);
 // double ln_k_max = log(pclass_sz->k_max_gas_density_profile);
 int n_k = pclass_sz->n_k_density_profile;
 int n_m = pclass_sz->n_m_density_profile;
 int n_z = pclass_sz->n_z_density_profile;

 class_alloc(pclass_sz->array_profile_ln_k,sizeof(double *)*n_k,pclass_sz->error_message);

 // array of masses:
 double ln_m_min = log(pclass_sz->M1SZ);
 double ln_m_max = log(pclass_sz->M2SZ);


 class_alloc(pclass_sz->array_profile_ln_m,sizeof(double *)*n_m,pclass_sz->error_message);


 // array of redshifts:
 double ln_1pz_min = log(1.+pclass_sz->z1SZ);
 double ln_1pz_max = log(1.+pclass_sz->z2SZ);


 class_alloc(pclass_sz->array_profile_ln_1pz,sizeof(double *)*n_z,pclass_sz->error_message);
int index_m_z;
int index_m_z_tab[n_m][n_z];
int index_k;
// for (index_k=0;
//      index_k<n_k;
//      index_k++)
// {
//   // this is k
//   pclass_sz->array_profile_ln_k[index_k] = ln_k_min
//                                       +index_k*(ln_k_max-ln_k_min)
//                                       /(n_k-1.);
// }

int index_m;
for (index_m=0;
     index_m<n_m;
     index_m++)
{
  pclass_sz->array_profile_ln_m[index_m] = ln_m_min
                                      +index_m*(ln_m_max-ln_m_min)
                                      /(n_m-1.);
}

int index_z;
for (index_z=0;
     index_z<n_z;
     index_z++)
{
  pclass_sz->array_profile_ln_1pz[index_z] = ln_1pz_min
                                   +index_z*(ln_1pz_max-ln_1pz_min)
                                   /(n_z-1.);
}


class_alloc(pclass_sz->array_profile_ln_rho_at_lnk_lnM_z,n_k*sizeof(double *),pclass_sz->error_message);
for (index_k=0;
     index_k<n_k;
     index_k++)
{
class_alloc(pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_k],n_m*n_z*sizeof(double *),pclass_sz->error_message);
index_m_z = 0;
for (index_m=0;
     index_m<n_m;
     index_m++){

for (index_z=0;
     index_z<n_z;
     index_z++)
{
  pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_k][index_m_z] = log(1e-100); // initialize with super small number
  index_m_z_tab[index_m][index_z] = index_m_z;
  index_m_z += 1;
}

     }
}

int has_ksz_bkp = pclass_sz->has_kSZ_kSZ_gal_1h;
pclass_sz->has_kSZ_kSZ_gal_1h = _TRUE_; //pretend we need the tau_profile

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


// index_m_z = 0;
#pragma omp parallel \
shared(abort,\
pclass_sz,pba,index_m_z_tab)\
private(tstart, tstop,index_k,index_z,index_m,index_m_z) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
for (index_z=0;
     index_z<n_z;
     index_z++){
#pragma omp flush(abort)

double * pvectsz;
double * pvecback;

class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);
class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);

int index_pvectsz;
for (index_pvectsz=0;
     index_pvectsz<pclass_sz->tsz_size;
     index_pvectsz++){
       pvectsz[index_pvectsz] = 0.; // set everything to 0.
     }


  double z = exp(pclass_sz->array_profile_ln_1pz[index_z])-1.;

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
  pvectsz[pclass_sz->index_z] = z;


  pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in rho_nfw for the pk mode computation

  pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);
  pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
  double omega = pvecback[pba->index_bg_Omega_m];
  pvectsz[pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);




for (index_m=0;
     index_m<n_m;
     index_m++){

  double lnM = pclass_sz->array_profile_ln_m[index_m];
  // index_m_z = index_z + index_m*n_m;
  index_m_z  = index_m_z_tab[index_m][index_z];
  // index_m_z = index_m+ index_z*n_z; // no

    pvectsz[pclass_sz->index_m200c] = exp(lnM);
    pvectsz[pclass_sz->index_r200c] = pow(3.*pvectsz[pclass_sz->index_m200c]/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
    pvectsz[pclass_sz->index_l200c] = sqrt(pvectsz[pclass_sz->index_chi2])/(1.+z)/pvectsz[pclass_sz->index_r200c];
    pvectsz[pclass_sz->index_characteristic_multipole_for_nfw_profile] = pvectsz[pclass_sz->index_l200c];


// only  do the integration of Battaglia profile
// nfw has an analytical formula
if (pclass_sz->tau_profile == 1){
// here we FFT the profile ===== commented
const int N = pclass_sz->N_samp_fftw; //precision parameter
int ix;
double x[N], Px[N];
double x_min = pclass_sz->x_min_gas_density_fftw;
double x_max = pclass_sz->x_max_gas_density_fftw;
// double x_max = pclass_sz->x_out_truncated_nfw_profile_electrons;
double x_out = pclass_sz->x_out_truncated_nfw_profile_electrons;
if (pclass_sz->use_xout_in_density_profile_from_enclosed_mass){
    x_out = get_m_to_xout_at_z_and_m(pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_m200c],pclass_sz);
    }

double A_rho0 = pclass_sz->A_rho0;
double A_alpha = pclass_sz->A_alpha;
double A_beta = pclass_sz->A_beta;

double alpha_m_rho0 = pclass_sz->alpha_m_rho0;
double alpha_m_alpha = pclass_sz->alpha_m_alpha;
double alpha_m_beta = pclass_sz->alpha_m_beta;

double alpha_z_rho0 = pclass_sz->alpha_z_rho0;
double alpha_z_alpha = pclass_sz->alpha_z_alpha;
double alpha_z_beta = pclass_sz->alpha_z_beta;
double gamma = pclass_sz->gamma_B16;
double xc = pclass_sz->xc_B16;


for (ix=0; ix<N; ix++){
  x[ix] = exp(log(x_min)+ix/(N-1.)*(log(x_max)-log(x_min)));
  if (x[ix]>x_out){
      Px[ix] = 0.;
    }
  else{
    double c_asked = pclass_sz->c_B16;
    Px[ix] =  get_gas_profile_at_x_M_z_b16_200c(x[ix],
                                                pvectsz[pclass_sz->index_m200c],
                                                z,
                                                c_asked,
                                                A_rho0,
                                                A_alpha,
                                                A_beta,
                                                alpha_m_rho0,
                                                alpha_m_alpha,
                                                alpha_m_beta,
                                                alpha_z_rho0,
                                                alpha_z_alpha,
                                                alpha_z_beta,
					                                      pclass_sz->mcut,
					                                      pclass_sz->alphap_m_rho0,
                                                pclass_sz->alphap_m_alpha,
                                                pclass_sz->alphap_m_beta,
					                                      pclass_sz->alpha_c_rho0,
                                                pclass_sz->alpha_c_alpha,
                                                pclass_sz->alpha_c_beta,
                                                gamma,
                                                xc,
                                                pba,
                                                pclass_sz);
                                              }
  // printf("x = %.3e Px = %.3e\n",x[ix],Px[ix]);
  }

  double kp[N], Pkp[N];
    // pk2xi(N,k,Pk1,rp,xi1,pclass_sz);
  /* Compute the function
   *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
   * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
   * in this notation.  The input k-values must be logarithmically spaced.  The
   * resulting xi_l^m(r) will be evaluated at the dual r-values
   *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
  //void fftlog_ComputeXiLM(int l, int m, int N, const double k[],  const double pk[], double r[], double xi[]);
  fftlog_ComputeXiLMsloz(0, 2, N, x, Px, kp, Pkp,pclass_sz);

/// ---- > FFT end commented

  double result;

  for (index_k=0;
       index_k<n_k;
       index_k++)
  {

    double k;
    double  result_fft;

  //   k = exp(pclass_sz->array_profile_ln_k[index_k]);
  //   pvectsz[pclass_sz->index_multipole_for_nfw_profile] = k;
  //
  // // pvectsz[pclass_sz->index_has_electron_density] = 1;
  // // do_mass_conversions(lnM,z,pvecback,pvectsz,pba,pclass_sz);
  //
  //   result_fft = 2.*_PI_*_PI_*pwl_value_1d(N,
  //                                          kp,
  //                                          Pkp,
  //                                          k*pvectsz[pclass_sz->index_r200c]*(1.+pvectsz[pclass_sz->index_z]));
  //

    pclass_sz->array_profile_ln_k[index_k] = log(kp[index_k]);//(pvectsz[pclass_sz->index_r200c]*(1.+pvectsz[pclass_sz->index_z])));
    result_fft = 2.*_PI_*_PI_*Pkp[index_k];


   // double result_int;
   // two_dim_ft_nfw_profile(pclass_sz,pba,pvectsz,&result_int);

   // printf("result_fft = %.5e result_qawo = %.5e ratio = %.5e \n",result_fft,result_int,result_int/result_fft);


   double result = result_fft;
   // double result = result_int;
   double tau_normalisation = 4.*_PI_*pow(pvectsz[pclass_sz->index_r200c],3);
   // printf("In tab gas: k %.4e z %.8e rt %.8e mt %.8e res = %.4e\n",ell,pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_r200c],pvectsz[pclass_sz->index_m200c],result);




   result *= tau_normalisation;

   if (result<=0 || isnan(result) || isinf(result)){
         // printf("ERROR: In tab gas: k %.4e z %.8e rt %.8e mt %.8e res = %.4e\n",ell,pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_r200c],pvectsz[pclass_sz->index_m200c],result);
         // printf("check precision and input parameters?\n");
         // exit(0);
         result = 1e-200;
        }

   pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_k][index_m_z] = log(result);
 } // k loop
} // density mode
//BCM model
if (pclass_sz->tau_profile == 2){
// here we FFT the profile ===== commented
const int N = pclass_sz->N_samp_fftw; //precision parameter
int ix;
double x[N], Px[N];
double x_min = pclass_sz->x_min_gas_density_fftw;
double x_max = pclass_sz->x_max_gas_density_fftw;
// double x_max = pclass_sz->x_out_truncated_nfw_profile_electrons;
double x_out = pclass_sz->x_out_truncated_nfw_profile_electrons;
// if (pclass_sz->use_xout_in_density_profile_from_enclosed_mass){
//     x_out = get_m_to_xout_at_z_and_m(pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_m200c],pclass_sz);
//     }

// double A_rho0 = pclass_sz->A_rho0;
// double A_alpha = pclass_sz->A_alpha;
// double A_beta = pclass_sz->A_beta;
//
// double alpha_m_rho0 = pclass_sz->alpha_m_rho0;
// double alpha_m_alpha = pclass_sz->alpha_m_alpha;
// double alpha_m_beta = pclass_sz->alpha_m_beta;
//
// double alpha_z_rho0 = pclass_sz->alpha_z_rho0;
// double alpha_z_alpha = pclass_sz->alpha_z_alpha;
// double alpha_z_beta = pclass_sz->alpha_z_beta;
// double gamma = pclass_sz->gamma_B16;
// double xc = pclass_sz->xc_B16;


for (ix=0; ix<N; ix++){
  x[ix] = exp(log(x_min)+ix/(N-1.)*(log(x_max)-log(x_min)));
  if (x[ix]>x_out){
      Px[ix] = 0.;
    }
  else{
    // double c_asked = pclass_sz->c_B16;
    Px[ix] =  get_gas_profile_at_x_M_z_bcm_200c(x[ix],
                                                pvectsz[pclass_sz->index_m200c],
                                                z,
                                                pba,
                                                pclass_sz);
                                              }
  // printf("x = %.3e Px = %.3e\n",x[ix],Px[ix]);
  }

  double kp[N], Pkp[N];
    // pk2xi(N,k,Pk1,rp,xi1,pclass_sz);
  /* Compute the function
   *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
   * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
   * in this notation.  The input k-values must be logarithmically spaced.  The
   * resulting xi_l^m(r) will be evaluated at the dual r-values
   *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
  //void fftlog_ComputeXiLM(int l, int m, int N, const double k[],  const double pk[], double r[], double xi[]);
  fftlog_ComputeXiLMsloz(0, 2, N, x, Px, kp, Pkp,pclass_sz);

/// ---- > FFT end commented

  double result;

  for (index_k=0;
       index_k<n_k;
       index_k++)
  {

  //   double k = exp(pclass_sz->array_profile_ln_k[index_k]);
  //   pvectsz[pclass_sz->index_multipole_for_nfw_profile] = k;
  //
  // // pvectsz[pclass_sz->index_has_electron_density] = 1;
  // // do_mass_conversions(lnM,z,pvecback,pvectsz,pba,pclass_sz);
  //
  // double  result_fft = 2.*_PI_*_PI_*pwl_value_1d(N,
  //                                               kp,
  //                                               Pkp,
  //                                               k*pvectsz[pclass_sz->index_r200c]*(1.+pvectsz[pclass_sz->index_z]));
  //
  //
  //

    double k;
    double  result_fft;

  //   k = exp(pclass_sz->array_profile_ln_k[index_k]);
  //   pvectsz[pclass_sz->index_multipole_for_nfw_profile] = k;
  //
  // // pvectsz[pclass_sz->index_has_electron_density] = 1;
  // // do_mass_conversions(lnM,z,pvecback,pvectsz,pba,pclass_sz);
  //
  //   result_fft = 2.*_PI_*_PI_*pwl_value_1d(N,
  //                                          kp,
  //                                          Pkp,
  //                                          k*pvectsz[pclass_sz->index_r200c]*(1.+pvectsz[pclass_sz->index_z]));
  //

    pclass_sz->array_profile_ln_k[index_k] = log(kp[index_k]);//(pvectsz[pclass_sz->index_r200c]*(1.+pvectsz[pclass_sz->index_z])));
    result_fft = 2.*_PI_*_PI_*Pkp[index_k];


   // double result_int;
   // two_dim_ft_nfw_profile(pclass_sz,pba,pvectsz,&result_int);

   // printf("result_fft = %.5e result_qawo = %.5e ratio = %.5e \n",result_fft,result_int,result_int/result_fft);


   double result = result_fft;
   // double result = result_int;
   double tau_normalisation = 4.*_PI_*pow(pvectsz[pclass_sz->index_r200c],3);
   // printf("In tab gas: k %.4e z %.8e rt %.8e mt %.8e res = %.4e\n",ell,pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_r200c],pvectsz[pclass_sz->index_m200c],result);




   result *= tau_normalisation;

   if (result<=0 || isnan(result) || isinf(result)){
         // printf("ERROR: In tab gas: k %.4e z %.8e rt %.8e mt %.8e res = %.4e\n",ell,pvectsz[pclass_sz->index_z],pvectsz[pclass_sz->index_r200c],pvectsz[pclass_sz->index_m200c],result);
         // printf("check precision and input parameters?\n");
         // exit(0);
         result = 1e-200;
        }


   pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_k][index_m_z] = log(result);
 } // k loop
} // density mode

// else if (pclass_sz->tau_profile == 0){ // truncated nfw profile
//
//
//     for (index_k=0;
//          index_k<n_k;
//          index_k++)
//     {
//
//       double k = exp(pclass_sz->array_profile_ln_k[index_k]);
//       // pvectsz[pclass_sz->index_multipole_for_nfw_profile] = k;
//    //
//    // pvectsz[pclass_sz->index_m200c] = exp(lnM);
//    // pvectsz[pclass_sz->index_r200c] = pow(3.*pvectsz[pclass_sz->index_m200c]/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
//    // pvectsz[pclass_sz->index_c200c] = get_c200c_at_m_and_z(pvectsz[pclass_sz->index_m200c],z,pba,pclass_sz);
//    //
//    // double r_delta = pvectsz[pclass_sz->index_r200c];
//    // double c_delta = pvectsz[pclass_sz->index_c200c];
//    // double m_delta = pvectsz[pclass_sz->index_m200c];
//
//    double r_delta;// = pvectsz[pclass_sz->index_radius_for_electron_density];
//    double c_delta;// = pvectsz[pclass_sz->index_concentration_for_electron_density];
//    double m_delta;// = pvectsz[pclass_sz->index_mass_for_electron_density];
//    // printf("de = %d\n",pclass_sz->delta_def_electron_density);
//    // exit(0);
//
//   if (pclass_sz->delta_def_electron_density == 0){
//     m_delta = exp(lnM);
//     r_delta = pow(3.*m_delta/(4.*_PI_*200.*pvecback[pba->index_bg_Omega_m]*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
//     c_delta = get_c200m_at_m_and_z(m_delta,z,pba,pclass_sz);
//   }
//   else if (pclass_sz->delta_def_electron_density == 1){
//     m_delta = exp(lnM);
//     r_delta = pow(3.*m_delta/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
//     c_delta = get_c200c_at_m_and_z(m_delta,z,pba,pclass_sz);
//   }
//   else if (pclass_sz->delta_def_electron_density == 2){
//     m_delta = exp(lnM);
//     r_delta = pow(3.*m_delta/(4.*_PI_*500.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
//     c_delta = get_c500c_at_m_and_z(m_delta,z,pba,pclass_sz);
//   }
//    double xout = pclass_sz->x_out_truncated_nfw_profile_electrons;
//
//
//
//    // pvectsz[pclass_sz->index_rs] = r_delta/c_delta;
//
//   // pvectsz[pclass_sz->index_multipole_for_truncated_nfw_profile] = pvectsz[pclass_sz->index_multipole_for_nfw_profile];
//   // double l = pvectsz[pclass_sz->index_multipole_for_truncated_nfw_profile];
//   double chi = sqrt(pvectsz[pclass_sz->index_chi2]);
//   // double k = ell;
//   double result =  evaluate_truncated_nfw_profile(pvectsz[pclass_sz->index_z],k,r_delta,c_delta,xout);
//    //result *= 1.;//m_delta;///(4.*_PI_*pow(pvectsz[pclass_sz->index_rs],3));
//    double f_b = pclass_sz->f_b_gas;//pba->Omega0_b/pclass_sz->Omega_m_0;
//    result *= f_b*m_delta;//*pow((pba->Omega0_cdm+pba->Omega0_b)*pclass_sz->Rho_crit_0,-1);
//
//     if (isnan(result) || isinf(result)){
//     printf("z %.8e rt %.8e ct %.8e mt %.8e\n",pvectsz[pclass_sz->index_z],r_delta,c_delta,m_delta);
//     exit(0);
//   }
//    // double tau_normalisation = 4.*_PI_*pow(pvectsz[pclass_sz->index_rs],3);
//    // result *= tau_normalisation;
//
//
//   pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_k][index_m_z] = log(result);
//   // printf("l = %.8e m = %.8e z = %.8e lnrho = %.8e\n",ell,exp(lnM),z,log(result));
//
//   // index_m_z += 1;
// } // k loop
//
//
// } // tau mode

} // m loop


     free(pvectsz);
     free(pvecback);
} // z loop

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in tab profile parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

// restore initial state:
pclass_sz->has_kSZ_kSZ_gal_1h = has_ksz_bkp;


                                      }




double get_gas_pressure_profile_at_k(double k_asked,
                                     struct class_sz_structure * pclass_sz){

  if (log(k_asked)<pclass_sz->array_pressure_profile_ln_k[0])
    return 0.;
  if (log(k_asked)>pclass_sz->array_pressure_profile_ln_k[pclass_sz->n_k_pressure_profile-1])
    return 0.;


    return  pwl_value_1d(pclass_sz->n_k_pressure_profile,
                         pclass_sz->array_pressure_profile_ln_k,
                         pclass_sz->array_pressure_profile_ln_p_at_lnk,
                         log(k_asked));
                                     }



// double get_gas_pressure_profile_at_k_m_z(double k_asked,
//                                          double m_asked,
//                                          double z_asked,
//                                          struct class_sz_structure * pclass_sz){
//   double z = log(1.+z_asked);
//   double m = log(m_asked);
//   double k = log(k_asked);


//   // if (pclass_sz->tau_profile == 1){
//   // find the closest l's in the grid:
//   int id_k_low;
//   int id_k_up;
//   int n_k = pclass_sz->n_k_pressure_profile;
//   int n_m = pclass_sz->n_m_pressure_profile;
//   int n_z = pclass_sz->n_z_pressure_profile;
//   r8vec_bracket(n_k,pclass_sz->array_pressure_profile_ln_k,k,&id_k_low,&id_k_up);

//   if (id_k_low == id_k_up){
//     printf("class_sz_tools.c bug in get_gas_pressure_profile_at_k_m_z");
//     exit(0);
//   }

//   if (m<pclass_sz->array_pressure_profile_ln_m[0])
//     return 0.;
//   if (m>pclass_sz->array_pressure_profile_ln_m[n_m-1])
//     return 0.;

//   // interpolate 2d at l_low:

//  double ln_rho_low = pwl_interp_2d(n_m,
//                                 n_z,
//                                 pclass_sz->array_pressure_profile_ln_m,
//                                 pclass_sz->array_pressure_profile_ln_1pz,
//                                 pclass_sz->array_pressure_profile_ln_p_at_lnk_lnm_z[id_k_low-1],
//                                 1,
//                                 &m,
//                                 &z);

//  double ln_rho_up = pwl_interp_2d(n_m,
//                                 n_z,
//                                 pclass_sz->array_pressure_profile_ln_m,
//                                 pclass_sz->array_pressure_profile_ln_1pz,
//                                 pclass_sz->array_pressure_profile_ln_p_at_lnk_lnm_z[id_k_up-1],
//                                 1,
//                                 &m,
//                                 &z);
//  double ln_k_low = pclass_sz->array_pressure_profile_ln_k[id_k_low-1];
//  double ln_k_up = pclass_sz->array_pressure_profile_ln_k[id_k_up-1];


//  double res = ln_rho_low + ((k - ln_k_low) / (ln_k_up - ln_k_low)) * (ln_rho_up - ln_rho_low);
// //  printf("interpolating between lmin  = %.5e and lmax = %.5e got %.5e\n",exp(ln_k_low),exp(ln_k_up), res );

//  return res;



// }




double get_gas_pressure_profile_at_l_m_z(double l_asked,
                                         double m_asked,
                                         double z_asked,
                                         struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
  double l = log(l_asked);


  // if (pclass_sz->tau_profile == 1){
  // find the closest l's in the grid:
  int id_l_low;
  int id_l_up;
  int n_l = pclass_sz->n_l_pressure_profile;
  int n_m = pclass_sz->n_m_pressure_profile;
  int n_z = pclass_sz->n_z_pressure_profile;
  r8vec_bracket(n_l,pclass_sz->array_pressure_profile_ln_l,l,&id_l_low,&id_l_up);

  if (id_l_low == id_l_up){
    printf("bug in get_gas_pressure_profile_at_l_m_z");
    exit(0);
  }

  if (m<pclass_sz->array_pressure_profile_ln_m[0])
    return 0.;
  if (m>pclass_sz->array_pressure_profile_ln_m[n_m-1])
    return 0.;

  // interpolate 2d at l_low:

 double ln_rho_low = pwl_interp_2d(n_m,
                                n_z,
                                pclass_sz->array_pressure_profile_ln_m,
                                pclass_sz->array_pressure_profile_ln_1pz,
                                pclass_sz->array_pressure_profile_ln_p_at_lnl_lnm_z[id_l_low-1],
                                1,
                                &m,
                                &z);

 double ln_rho_up = pwl_interp_2d(n_m,
                                n_z,
                                pclass_sz->array_pressure_profile_ln_m,
                                pclass_sz->array_pressure_profile_ln_1pz,
                                pclass_sz->array_pressure_profile_ln_p_at_lnl_lnm_z[id_l_up-1],
                                1,
                                &m,
                                &z);
 double ln_l_low = pclass_sz->array_pressure_profile_ln_l[id_l_low-1];
 double ln_l_up = pclass_sz->array_pressure_profile_ln_l[id_l_up-1];

 return ln_rho_low + ((l - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low);



}



// Tabulate 2D Fourier transform of custom1 on a [ln1pz - ln_M - ln_k] grid
int tabulate_custom1_profile_fft(struct background * pba,
                                 struct class_sz_structure * pclass_sz){


 int n_k = pclass_sz->n_k_custom1_profile;
 int n_m = pclass_sz->n_m_custom1_profile;
 int n_z = pclass_sz->n_z_custom1_profile;
 int index_m_z_tab[n_m][n_z];


int index_k;
int index_m;
int index_z;

int index_m_z = 0;


for (index_z=0;
     index_z<n_z;
     index_z++)
{


for (index_m=0;
     index_m<n_m;
     index_m++){

  index_m_z_tab[index_m][index_z] = index_m_z;
  index_m_z += 1;
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
pclass_sz,pba,index_m_z_tab)\
private(tstart, tstop,index_k,index_z,index_m,index_m_z) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
for (index_z=0;
     index_z<n_z;
     index_z++){
#pragma omp flush(abort)




double z = exp(pclass_sz->array_custom1_profile_ln_1pz[index_z])-1.;

for (index_m=0;
     index_m<n_m;
     index_m++){

  double m = exp(pclass_sz->array_custom1_profile_ln_m[index_m]);

  index_m_z  = index_m_z_tab[index_m][index_z];

  int ix;
  double x[n_k], Px[n_k];



  double x_out = pclass_sz->x_out_custom1;


for (ix=0; ix<n_k; ix++){
  x[ix] = exp(pclass_sz->array_custom1_profile_ln_x[ix]);
  if (x[ix]>x_out){
      Px[ix] = 0.;
    }
  else{
    Px[ix] = get_custom1_profile_at_x_m_z(x[ix],m,z,pclass_sz);
    }
  }

  double kp[n_k], Pkp[n_k];
    // pk2xi(N,k,Pk1,rp,xi1,pclass_sz);
  /* Compute the function
   *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
   * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
   * in this notation.  The input k-values must be logarithmically spaced.  The
   * resulting xi_l^m(r) will be evaluated at the dual r-values
   *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
  //void fftlog_ComputeXiLM(int l, int m, int N, const double k[],  const double pk[], double r[], double xi[]);
  fftlog_ComputeXiLMsloz(0, 2, n_k, x, Px, kp, Pkp,pclass_sz);



  for (index_k=0;
       index_k<n_k;
       index_k++)
  {


  pclass_sz->array_custom1_profile_ln_k[index_k] = log(kp[index_k]);

  double  result_fft = 2.*_PI_*_PI_*Pkp[index_k];

  pclass_sz->array_custom1_profile_u_at_lnk_lnm_ln1pz[index_k][index_m_z] = log(result_fft);

  } // k loop

} // m loop


} // z loop

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in tab fft custom1 profile parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;

                                      }




// Tabulate 2D Fourier transform of pressure profile on a [z - ln_M - ln_k] grid
int tabulate_gas_pressure_profile_B12_fft(struct background * pba,
                                          struct class_sz_structure * pclass_sz){

 // // array of multipoles:
 // double ln_k_min = log(pclass_sz->k_min_gas_pressure_profile);
 // double ln_k_max = log(pclass_sz->k_max_gas_pressure_profile);


double lnl_min = log(pclass_sz->l_min_gas_pressure_profile);
double lnl_max = log(pclass_sz->l_max_gas_pressure_profile);


 int n_l = pclass_sz->n_l_pressure_profile;
 int n_m = pclass_sz->n_m_pressure_profile;
 int n_z = pclass_sz->n_z_pressure_profile;

class_alloc(pclass_sz->array_pressure_profile_ln_l,sizeof(double *)*n_l,pclass_sz->error_message);


//  printf("allocate arrays.\n");

//  class_alloc(pclass_sz->array_pressure_profile_ln_k,sizeof(double *)*n_k,pclass_sz->error_message);

 // array of masses:
 double ln_m_min = log(pclass_sz->M1SZ);
 double ln_m_max = log(pclass_sz->M2SZ);


 class_alloc(pclass_sz->array_pressure_profile_ln_m,sizeof(double *)*n_m,pclass_sz->error_message);


 // array of redshifts:
 double ln_1pz_min = log(1.+pclass_sz->z1SZ);
 double ln_1pz_max = log(1.+pclass_sz->z2SZ);


 class_alloc(pclass_sz->array_pressure_profile_ln_1pz,sizeof(double *)*n_z,pclass_sz->error_message);
int index_m_z;
int index_m_z_tab[n_m][n_z];
int index_l;

//  printf("allocate arrays2.\n");
// for (index_k=0;
//      index_k<n_k;
//      index_k++)
// {
//   pclass_sz->array_pressure_profile_ln_k[index_k] = ln_k_min
//                                               +index_k*(ln_k_max-ln_k_min)
//                                               /(n_k-1.);
// }

    for (index_l=0;
         index_l<n_l;
         index_l++)
      {
        pclass_sz->array_pressure_profile_ln_l[index_l] = lnl_min
                                            +index_l*(lnl_max-lnl_min)
                                            /(n_l-1.);
      }



int index_m;
for (index_m=0;
     index_m<n_m;
     index_m++)
{
  pclass_sz->array_pressure_profile_ln_m[index_m] = ln_m_min
                                      +index_m*(ln_m_max-ln_m_min)
                                      /(n_m-1.);
}

int index_z;
for (index_z=0;
     index_z<n_z;
     index_z++)
{
  pclass_sz->array_pressure_profile_ln_1pz[index_z] = ln_1pz_min
                                   +index_z*(ln_1pz_max-ln_1pz_min)
                                   /(n_z-1.);
}

//  printf("allocate arrays 3.\n");

class_alloc(pclass_sz->array_pressure_profile_ln_p_at_lnl_lnm_z,n_l*sizeof(double *),pclass_sz->error_message);
for (index_l=0;
     index_l<n_l;
     index_l++)
      {
      class_alloc(pclass_sz->array_pressure_profile_ln_p_at_lnl_lnm_z[index_l],n_m*n_z*sizeof(double *),pclass_sz->error_message);
      index_m_z = 0;


                  //class_alloc(pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z],n_z*sizeof(double ),pclass_sz->error_message);


            for (index_z=0;
                index_z<n_z;
                index_z++)
            {



      for (index_m=0;
          index_m<n_m;
          index_m++){





                    // pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z] = -100.; // initialize with super small number
                    pclass_sz->array_pressure_profile_ln_p_at_lnl_lnm_z[index_l][index_m_z] = 1e-100; // initialize with super small number
                    index_m_z_tab[index_m][index_z] = index_m_z;
                    index_m_z += 1;
                  }

                }
      }


// printf("starting parallel block.\n");

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
pclass_sz,pba,index_m_z_tab)\
private(tstart, index_z,index_m,tstop,index_m_z) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
for (index_z=0;
     index_z<n_z;
     index_z++){
#pragma omp flush(abort)

// int index_z,index_m;
double * pvectsz;
double * pvecback;
class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);
class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);


int index_pvectsz;

// int index_l;
// int n_l = pclass_sz->n_l_pressure_profile;

for (index_pvectsz=0;
     index_pvectsz<pclass_sz->tsz_size;
     index_pvectsz++){
       pvectsz[index_pvectsz] = 0.; // set everything to 0.
     }

  double z = exp(pclass_sz->array_pressure_profile_ln_1pz[index_z])-1.;
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
  pvectsz[pclass_sz->index_z] = z;


  pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
  double chi = sqrt(pvectsz[pclass_sz->index_chi2]);
  // pvectsz[pclass_sz->index_multipole_for_pressure_profile] = k*chi;
  pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in p_gnfw for the pk mode computation

  pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);

  double omega = pvecback[pba->index_bg_Omega_m];
  pvectsz[pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);



for (index_m=0;
     index_m<n_m;
     index_m++){
  double lnM = pclass_sz->array_pressure_profile_ln_m[index_m];
  // index_m_z = index_m+ index_z*n_z;
  index_m_z  = index_m_z_tab[index_m][index_z];






  pvectsz[pclass_sz->index_m200c] = exp(lnM);

  pvectsz[pclass_sz->index_r200c] = pow(3.*pvectsz[pclass_sz->index_m200c]/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
  pvectsz[pclass_sz->index_l200c] = sqrt(pvectsz[pclass_sz->index_chi2])/(1.+z)/pvectsz[pclass_sz->index_r200c];
  // pvectsz[pclass_sz->index_characteristic_multipole_for_profile] = pvectsz[pclass_sz->index_l200c];
  pvectsz[pclass_sz->index_rs] = pvectsz[pclass_sz->index_r200c];///pvectsz[pclass_sz->index_c200c];




  const int N = pclass_sz->N_samp_fftw; //precision parameter
  int ix;
  double x[N], Px[N];
  double x_min = pclass_sz->x_min_gas_pressure_fftw;
  double x_max = pclass_sz->x_max_gas_pressure_fftw;


  double r200c = pvectsz[pclass_sz->index_r200c]; //in Mpc/h
  double x_out;
  if (pclass_sz->truncate_gas_pressure_wrt_rvir == 1){

          // class_call_parallel(mDEL_to_mVIR(pvectsz[pclass_sz->index_m200c],
          //                                   200.*(pvectsz[pclass_sz->index_Rho_crit]),
          //                                   pvectsz[pclass_sz->index_Delta_c],
          //                                   pvectsz[pclass_sz->index_Rho_crit],
          //                                   z,
          //                                   &pvectsz[pclass_sz->index_mVIR],
          //                                   pclass_sz,
          //                                   pba),
          //                 pclass_sz->error_message,
          //                 pclass_sz->error_message);
          // //
          // //  // rvir needed to cut off the integral --> e.g., xout = 50.*rvir/r200c
          //  pvectsz[pclass_sz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[pclass_sz->index_mVIR],pvectsz[pclass_sz->index_Delta_c],pvectsz[pclass_sz->index_Rho_crit],pclass_sz);

          double mvir = pvectsz[pclass_sz->index_mVIR];

          pvectsz[pclass_sz->index_mVIR] = get_m200c_to_mvir_at_z_and_M(z,pvectsz[pclass_sz->index_m200c],pclass_sz);

          // printf("mvir/mvir = %.5e\n",mvir/pvectsz[pclass_sz->index_mVIR]);

          pvectsz[pclass_sz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[pclass_sz->index_mVIR],
                                                            pvectsz[pclass_sz->index_Delta_c],
                                                            pvectsz[pclass_sz->index_Rho_crit],
                                                            pclass_sz);

          double rvir = pvectsz[pclass_sz->index_rVIR]; //in Mpc/h

          x_out = pclass_sz->x_outSZ*rvir/r200c; // the truncation radius is in multiples of rvir

          }
  else{
          x_out = pclass_sz->x_outSZ;
      }

for (ix=0; ix<N; ix++){
    x[ix] = exp(log(x_min)+ix/(N-1.)*(log(x_max)-log(x_min)));
    if (x[ix]>x_out){
        Px[ix] = 0.;
      }
    else{

          // double xc;
          // double beta;
          // double P0;
          //
          double m200_over_msol = pvectsz[pclass_sz->index_m200c]/pba->h; // convert to Msun
          double z = pvectsz[pclass_sz->index_z];
          //
          //
          double P0 = pclass_sz->P0_B12*pow(m200_over_msol/1e14,pclass_sz->alpha_m_P0_B12)*pow(1+z,pclass_sz->alpha_z_P0_B12);
          double xc = pclass_sz->xc_B12*pow(m200_over_msol/1e14,pclass_sz->alpha_m_xc_B12)*pow(1+z,pclass_sz->alpha_z_xc_B12);
          double beta = pclass_sz->beta_B12*pow(m200_over_msol/1e14,pclass_sz->alpha_m_beta_B12)*pow(1+z,pclass_sz->alpha_z_beta_B12);

          double gamma = pclass_sz->gamma_B12;
          double alpha = pclass_sz->alpha_B12;

          double p_gnfw_x = P0*pow(x[ix]/xc,gamma)*pow(1.+ pow(x[ix]/xc,alpha),-beta);
          Px[ix] = p_gnfw_x;

          // double c_asked = pclass_sz->c_B12;//what we pass there?
          // Px[ix] = get_pressure_P_over_P_delta_at_x_M_z_b12_200c(x[ix],m200_over_msol,z,
          //                                                       c_asked,pclass_sz->P0_B12,
          //                                                       pclass_sz->xc_B12,pclass_sz->beta_B12,
          //                                                       pclass_sz->alpha_m_P0_B12,pclass_sz->alpha_m_xc_B12,
          //                                                       pclass_sz->alpha_m_beta_B12,pclass_sz->alpha_z_P0_B12,
          //                                                       pclass_sz->alpha_z_xc_B12,pclass_sz->alpha_z_beta_B12,
          //                                                       // break model
          //                                                       pclass_sz->mcut_B12,pclass_sz->alphap_m_P0_B12,
          //                                                       pclass_sz->alphap_m_xc_B12,pclass_sz->alphap_m_beta_B12,
          //                                                       pclass_sz->alpha_c_P0_B12,
          //                                                       pclass_sz->alpha_c_xc_B12,
          //                                                       pclass_sz->alpha_c_beta_B12,
          //                                                             // end break model
          //                                                       pclass_sz->alpha_B12,
          //                                                       pclass_sz->gamma_B12,
          //                                                       pba,pclass_sz);




      }
  // printf("x = %.3e Px = %.3e\n",x[ix],Px[ix]);

}

  double kp[N], Pkp[N];
    // pk2xi(N,k,Pk1,rp,xi1,pclass_sz);
  /* Compute the function
   *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
   * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
   * in this notation.  The input k-values must be logarithmically spaced.  The
   * resulting xi_l^m(r) will be evaluated at the dual r-values
   *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
  //void fftlog_ComputeXiLM(int l, int m, int N, const double k[],  const double pk[], double r[], double xi[]);
  fftlog_ComputeXiLMsloz(0, 2, N, x, Px, kp, Pkp,pclass_sz);


  // printf("FFT done.\n");
  int index_l;
  int n_l = pclass_sz->n_l_pressure_profile;
  for (index_l=0;
       index_l<n_l;
       index_l++)
  {


  // double k = exp(pclass_sz->array_pressure_profile_ln_k[index_k]); // l/ls
  // printf("calling ft\n");

  double l = exp(pclass_sz->array_pressure_profile_ln_l[index_l]);
  // printf("l = %.5e\n",l);
  double kl = (l+0.5)/pvectsz[pclass_sz->index_l200c];

  double result_fft;
  // printf("interpolating....\n");

  if ((kl<kp[0])||(kl>kp[N-1]))
    result_fft = 0.;
  else
    result_fft = 2.*_PI_*_PI_*pwl_value_1d(N,kp,Pkp,kl);

  if (result_fft<0.)  result_fft = 0.;


  // printf("z = %.5e m = %.5e l = %.5e result_fft = %.5e\n",z,pvectsz[pclass_sz->index_m200c],l,result_fft);


  // double lp = kp[index_k]*pvectsz[pclass_sz->index_l200c]-0.5;
  // if (lp<0.){
  //   Pkp[index_k] = 0.;
  //   // lp = 1.e-100;
  // }

  //
  // pclass_sz->array_pressure_profile_ln_k[index_k] = log(kp[index_k]);

  // pclass_sz->array_pressure_profile_ln_k[index_k] = log(lp);


  // double  result_fft = 2.*_PI_*_PI_*Pkp[index_k];


   // double result_int;
   //
   // // printf("calling ft\n");
   // pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in p_gnfw for the pk mode computation
   // class_call_parallel(two_dim_ft_pressure_profile(k,
   //                                                 pclass_sz,pba,pvectsz,&result_int),
   //                                                 pclass_sz->error_message,
   //                                                 pclass_sz->error_message);
   // double result = result_int;
   //
   // printf("result_fft = %.5e result_qawo = %.5e ratio = %.5e \n",result_fft,result_int,result_int/result_fft);


 // double tau_normalisation = 1.;//pvectsz[pclass_sz->index_m200m];///(4.*_PI_*pow(pvectsz[pclass_sz->index_rs],3.));
 // tau_normalisation = 4.*_PI_*pow(pvectsz[pclass_sz->index_r200c],3);//*pvectsz[pclass_sz->index_Rho_crit];
 // result *= tau_normalisation;


  // pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z] = log(result);
  pclass_sz->array_pressure_profile_ln_p_at_lnl_lnm_z[index_l][index_m_z] = result_fft;
  // printf("ell = %.3e z = %.3e m = %.3e lnrho = %.3e\n",l,z,exp(lnM),pclass_sz->array_pressure_profile_ln_p_at_lnk_lnm_z[index_k][index_m_z]);
  // index_m_z += 1;
} // k loop


} // m loop

free(pvectsz);
free(pvecback);

} // z loop

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in tab pressure profile parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;

// exit(0);

                                      }


// // Tabulate 2D Fourier transform of pressure profile on a [z - ln_M - ln_ell] grid
int tabulate_gas_pressure_profile_B12_l(struct background * pba,
                                  struct class_sz_structure * pclass_sz){

 // array of multipoles:
 double lnl_min = log(pclass_sz->l_min_gas_pressure_profile);
 double lnl_max = log(pclass_sz->l_max_gas_pressure_profile);
 int n_l = pclass_sz->n_l_pressure_profile;
 int n_m = pclass_sz->n_m_pressure_profile;
 int n_z = pclass_sz->n_z_pressure_profile;

 class_alloc(pclass_sz->array_pressure_profile_ln_l,sizeof(double *)*n_l,pclass_sz->error_message);

 // array of masses:
 double ln_m_min = log(pclass_sz->M1SZ);
 double ln_m_max = log(pclass_sz->M2SZ);


 class_alloc(pclass_sz->array_pressure_profile_ln_m,sizeof(double *)*n_m,pclass_sz->error_message);


 // array of redshifts:
 double ln_1pz_min = log(1.+pclass_sz->z1SZ);
 double ln_1pz_max = log(1.+pclass_sz->z2SZ);


 class_alloc(pclass_sz->array_pressure_profile_ln_1pz,sizeof(double *)*n_z,pclass_sz->error_message);
int index_m_z;

int index_l;
for (index_l=0;
     index_l<n_l;
     index_l++)
  {
    pclass_sz->array_pressure_profile_ln_l[index_l] = lnl_min
                                        +index_l*(lnl_max-lnl_min)
                                        /(n_l-1.);
  }

int index_m;
for (index_m=0;
     index_m<n_m;
     index_m++)
  {
    pclass_sz->array_pressure_profile_ln_m[index_m] = ln_m_min
                                        +index_m*(ln_m_max-ln_m_min)
                                        /(n_m-1.);
  }

int index_z;
for (index_z=0;
     index_z<n_z;
     index_z++)
  {
    pclass_sz->array_pressure_profile_ln_1pz[index_z] = ln_1pz_min
                                    +index_z*(ln_1pz_max-ln_1pz_min)
                                    /(n_z-1.);
  }


class_alloc(pclass_sz->array_pressure_profile_ln_p_at_lnl_lnm_z,n_l*sizeof(double *),pclass_sz->error_message);
for (index_l=0;
     index_l<n_l;
     index_l++)
      {
      class_alloc(pclass_sz->array_pressure_profile_ln_p_at_lnl_lnm_z[index_l],n_m*n_z*sizeof(double *),pclass_sz->error_message);
      index_m_z = 0;


      // int index_z;
      for (index_z=0;
          index_z<n_z;
          index_z++)
          {



          for (index_m=0;
              index_m<n_m;
              index_m++){



      //class_alloc(pclass_sz->array_profile_ln_rho_at_lnl_lnM_z[index_l][index_m_z],n_z*sizeof(double ),pclass_sz->error_message);


                  // pclass_sz->array_profile_ln_rho_at_lnl_lnM_z[index_l][index_m_z] = -100.; // initialize with super small number
                  pclass_sz->array_pressure_profile_ln_p_at_lnl_lnm_z[index_l][index_m_z] = 1e-100; // initialize with super small number
                  index_m_z += 1;
                }

          }
      }

// printf("starting para...\n");
// exit(0);

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
pclass_sz,pba)\
private(tstart, tstop,index_l,index_z,index_m,index_m_z) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)
for (index_l=0;
     index_l<n_l;
     index_l++)
{
#pragma omp flush(abort)
double * pvectsz;
double * pvecback;
class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);
class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
int index_pvectsz;
for (index_pvectsz=0;
     index_pvectsz<pclass_sz->tsz_size;
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



  double z = exp(pclass_sz->array_pressure_profile_ln_1pz[index_z])-1.;
  double lnM = pclass_sz->array_pressure_profile_ln_m[index_m];
  double ell = exp(pclass_sz->array_pressure_profile_ln_l[index_l]);


// if (pclass_sz->sz_verbose>1)
//   printf("--->In parallel region B12 (class_sz_tool.c) l = %.3e m = %.3e z = %.3e\n",
//         ell,
//         exp(lnM),
//         z);

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
  pvectsz[pclass_sz->index_z] = z;

  pvectsz[pclass_sz->index_multipole_for_pressure_profile] = ell;
  pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in rho_nfw for the pk mode computation

  pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                *pow(_Mpc_over_m_,1)
                                *pow(_c_,2)
                                *pvecback[pba->index_bg_rho_crit]
                                /pow(pba->h,2);
  pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
  double omega = pvecback[pba->index_bg_Omega_m];
  pvectsz[pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);


  double result;
  pvectsz[pclass_sz->index_m200c] = exp(lnM);
  // class_call_parallel(mDEL_to_mVIR(pvectsz[pclass_sz->index_m200c],
  //                                  200.*(pvectsz[pclass_sz->index_Rho_crit]),
  //                                  pvectsz[pclass_sz->index_Delta_c],
  //                                  pvectsz[pclass_sz->index_Rho_crit],
  //                                  z,
  //                                  &pvectsz[pclass_sz->index_mVIR],
  //                                  pclass_sz,
  //                                  pba),
  //                 pclass_sz->error_message,
  //                 pclass_sz->error_message);
  if (pclass_sz->sz_verbose>2)
    printf("----> tab B12 getting mvir\n");

  if (pclass_sz->truncate_gas_pressure_wrt_rvir)
    pvectsz[pclass_sz->index_mVIR] = get_m200c_to_mvir_at_z_and_M(z,pvectsz[pclass_sz->index_m200c],pclass_sz);
  else
    pvectsz[pclass_sz->index_mVIR] = pvectsz[pclass_sz->index_m200c];

  if (pclass_sz->sz_verbose>2)
    printf("----> tab B12 got mvir = %.3e\n",pvectsz[pclass_sz->index_mVIR]);
 //
 //  // rvir needed to cut off the integral --> e.g., xout = 50.*rvir/r200c
  pvectsz[pclass_sz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[pclass_sz->index_mVIR],pvectsz[pclass_sz->index_Delta_c],pvectsz[pclass_sz->index_Rho_crit],pclass_sz);
  pvectsz[pclass_sz->index_r200c] = pow(3.*pvectsz[pclass_sz->index_m200c]/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
  pvectsz[pclass_sz->index_l200c] = sqrt(pvectsz[pclass_sz->index_chi2])/(1.+z)/pvectsz[pclass_sz->index_r200c];
  // pvectsz[pclass_sz->index_characteristic_multipole_for_profile] = pvectsz[pclass_sz->index_l200c];
  pvectsz[pclass_sz->index_rs] = pvectsz[pclass_sz->index_r200c];///pvectsz[pclass_sz->index_c200c];
 double result_int;
 double kl = (ell+1./2.)/pvectsz[pclass_sz->index_l200c];
 class_call_parallel(two_dim_ft_pressure_profile(kl,pclass_sz,pba,pvectsz,&result_int),
                                    pclass_sz->error_message,
                                    pclass_sz->error_message);
 result = result_int;



 // double tau_normalisation = 1.;//pvectsz[pclass_sz->index_m200m];///(4.*_PI_*pow(pvectsz[pclass_sz->index_rs],3.));
 // tau_normalisation = 4.*_PI_*pow(pvectsz[pclass_sz->index_r200c],3);//*pvectsz[pclass_sz->index_Rho_crit];
 // result *= tau_normalisation;


  // pclass_sz->array_profile_ln_rho_at_lnl_lnM_z[index_l][index_m_z] = log(result);
  pclass_sz->array_pressure_profile_ln_p_at_lnl_lnm_z[index_l][index_m_z] = result;
  // printf("ell = %.3e z = %.3e m = %.3e rho = %.3e\n",ell,z,exp(lnM),log(result));

  // printf("--------->In parallel region B12 l (class_sz_tool.c)  l = %.5e m = %.5e z = %.5e res = %.5e\n",
  //       ell,
  //       exp(lnM),
  //       z,
  //       result);


  index_m_z += 1;
     }


     }

free(pvectsz);
free(pvecback);

}
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in tab profile parallel region (loop over ell's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;

                                      }



// // Tabulate 2D Fourier transform of pressure profile on a [z - ln_M - ln_ell] grid
// int tabulate_gas_pressure_profile_B12(struct background * pba,
//                                       struct class_sz_structure * pclass_sz){


// if (pclass_sz->sz_verbose>1)
//   printf("---> Start tabulating B12 (class_sz_tool.c)\n");

//  // array of multipoles:
//  double ln_k_min = log(pclass_sz->k_min_gas_pressure_profile);
//  double ln_k_max = log(pclass_sz->k_max_gas_pressure_profile);
//  int n_k = pclass_sz->n_k_pressure_profile;
//  int n_m = pclass_sz->n_m_pressure_profile;
//  int n_z = pclass_sz->n_z_pressure_profile;

//  class_alloc(pclass_sz->array_pressure_profile_ln_k,sizeof(double *)*n_k,pclass_sz->error_message);

//  // array of masses:
//  double ln_m_min = log(pclass_sz->M1SZ);
//  double ln_m_max = log(pclass_sz->M2SZ);


//  class_alloc(pclass_sz->array_pressure_profile_ln_m,sizeof(double *)*n_m,pclass_sz->error_message);


//  // array of redshifts:
//  double ln_1pz_min = log(1.+pclass_sz->z1SZ);
//  double ln_1pz_max = log(1.+pclass_sz->z2SZ);


//  class_alloc(pclass_sz->array_pressure_profile_ln_1pz,sizeof(double *)*n_z,pclass_sz->error_message);
// int index_m_z;

// int index_k;
// for (index_k=0;
//      index_k<n_k;
//      index_k++)
// {
//   pclass_sz->array_pressure_profile_ln_k[index_k] = ln_k_min
//                                               +index_k*(ln_k_max-ln_k_min)
//                                               /(n_k-1.);
// }

// int index_m;
// for (index_m=0;
//      index_m<n_m;
//      index_m++)
// {
//   pclass_sz->array_pressure_profile_ln_m[index_m] = ln_m_min
//                                       +index_m*(ln_m_max-ln_m_min)
//                                       /(n_m-1.);
// }

// int index_z;
// for (index_z=0;
//      index_z<n_z;
//      index_z++)
// {
//   pclass_sz->array_pressure_profile_ln_1pz[index_z] = ln_1pz_min
//                                    +index_z*(ln_1pz_max-ln_1pz_min)
//                                    /(n_z-1.);
// }


// class_alloc(pclass_sz->array_pressure_profile_ln_p_at_lnk_lnm_z,n_k*sizeof(double *),pclass_sz->error_message);
// for (index_k=0;
//      index_k<n_k;
//      index_k++)
// {
// class_alloc(pclass_sz->array_pressure_profile_ln_p_at_lnk_lnm_z[index_k],n_m*n_z*sizeof(double *),pclass_sz->error_message);
// index_m_z = 0;
// for (index_m=0;
//      index_m<n_m;
//      index_m++){
//     //class_alloc(pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z],n_z*sizeof(double ),pclass_sz->error_message);

//     for (index_z=0;
//         index_z<n_z;
//         index_z++)
//     {
//       // pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z] = -100.; // initialize with super small number
//       pclass_sz->array_pressure_profile_ln_p_at_lnk_lnm_z[index_k][index_m_z] = 1e-100; // initialize with super small number
//       index_m_z += 1;
//     }

//         }
// }


// if (pclass_sz->sz_verbose>1)
//   printf("---> Start parallel region B12 (class_sz_tool.c)\n");
// //Parallelization of profile computation
// /* initialize error management flag */


// int abort;
// double tstart, tstop;
// abort = _FALSE_;
// /* beginning of parallel region */

// int number_of_threads= 1;
// #ifdef _OPENMP
// #pragma omp parallel
//   {
//     number_of_threads = omp_get_num_threads();
//   }
// #endif

// #pragma omp parallel \
// shared(abort,\
// pclass_sz,pba)\
// private(tstart, tstop,index_k,index_z,index_m,index_m_z) \
// num_threads(number_of_threads)
// {

// #ifdef _OPENMP
//   tstart = omp_get_wtime();
// #endif

// #pragma omp for schedule (dynamic)
// for (index_k=0;
//      index_k<n_k;
//      index_k++)
// {
// #pragma omp flush(abort)
// double * pvectsz;
// double * pvecback;
// class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);
// class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
// int index_pvectsz;
// for (index_pvectsz=0;
//      index_pvectsz<pclass_sz->tsz_size;
//      index_pvectsz++){
//        pvectsz[index_pvectsz] = 0.; // set everything to 0.
//      }
// index_m_z = 0;
// for (index_z=0;
//      index_z<n_z;
//      index_z++){
// for (index_m=0;
//      index_m<n_m;
//      index_m++){



//   double z = exp(pclass_sz->array_pressure_profile_ln_1pz[index_z])-1.;
//   double lnM = pclass_sz->array_pressure_profile_ln_m[index_m];
//   double k = exp(pclass_sz->array_pressure_profile_ln_k[index_k]); // l/ls
// // printf("calling ft\n");

// // if (pclass_sz->sz_verbose>1)
// //   printf("--->In parallel region B12 (class_sz_tool.c) k = %.3e m = %.3e z = %.3e\n",
// //         k,
// //         exp(lnM),
// //         z);

//   // double l = sqrt(pvectsz[pclass_sz->index_chi2])*k-0.5;
//   // printf("l = %.5e\n",l);

//   double tau;
//   int first_index_back = 0;


//   class_call_parallel(background_tau_of_z(pba,z,&tau),
//              pba->error_message,
//              pba->error_message);

//   class_call_parallel(background_at_tau(pba,
//                                tau,
//                                pba->long_info,
//                                pba->inter_normal,
//                                &first_index_back,
//                                pvecback),
//              pba->error_message,
//              pba->error_message);


//   // fill relevant entries
//   pvectsz[pclass_sz->index_z] = z;


//   pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);
//   double chi = sqrt(pvectsz[pclass_sz->index_chi2]);
//   // pvectsz[pclass_sz->index_multipole_for_pressure_profile] = k*chi;
//   pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in p_gnfw for the pk mode computation

//   pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
//                                 *pow(_Mpc_over_m_,1)
//                                 *pow(_c_,2)
//                                 *pvecback[pba->index_bg_rho_crit]
//                                 /pow(pba->h,2);

//   double omega = pvecback[pba->index_bg_Omega_m];
//   pvectsz[pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);


//   double result;
//   pvectsz[pclass_sz->index_m200c] = exp(lnM);
//   class_call_parallel(mDEL_to_mVIR(pvectsz[pclass_sz->index_m200c],
//                                    200.*(pvectsz[pclass_sz->index_Rho_crit]),
//                                    pvectsz[pclass_sz->index_Delta_c],
//                                    pvectsz[pclass_sz->index_Rho_crit],
//                                    z,
//                                    &pvectsz[pclass_sz->index_mVIR],
//                                    pclass_sz,
//                                    pba),
//                   pclass_sz->error_message,
//                   pclass_sz->error_message);
//   // if (pclass_sz->sz_verbose>1)
//   //   printf("----> tab B12 getting mvir = %.5e\n",pvectsz[pclass_sz->index_mVIR]);
//   double mvir  = pvectsz[pclass_sz->index_mVIR];

//   // if (pclass_sz->truncate_gas_pressure_wrt_rvir){

//   //   pvectsz[pclass_sz->index_mVIR] = get_m200c_to_mvir_at_z_and_M(z,pvectsz[pclass_sz->index_m200c],pclass_sz);
//   //   if (pclass_sz->sz_verbose>1){
//   //       printf("truncating wrt rvir\n");
//   //     }
//   // }
//   // else
//   //   pvectsz[pclass_sz->index_mVIR] = pvectsz[pclass_sz->index_m200c];

//   // if (pclass_sz->sz_verbose>1)
//   //   printf("----> tab B12 got mvir = %.3e ratio = %.18e\n",pvectsz[pclass_sz->index_mVIR], pvectsz[pclass_sz->index_mVIR]/mvir);
//  //
//  //  // rvir needed to cut off the integral --> e.g., xout = 50.*rvir/r200c
//   pvectsz[pclass_sz->index_rVIR] = evaluate_rvir_of_mvir(pvectsz[pclass_sz->index_mVIR],pvectsz[pclass_sz->index_Delta_c],pvectsz[pclass_sz->index_Rho_crit],pclass_sz);
//   pvectsz[pclass_sz->index_r200c] = pow(3.*pvectsz[pclass_sz->index_m200c]/(4.*_PI_*200.*pvectsz[pclass_sz->index_Rho_crit]),1./3.);
//   pvectsz[pclass_sz->index_l200c] = sqrt(pvectsz[pclass_sz->index_chi2])/(1.+z)/pvectsz[pclass_sz->index_r200c];
//   // pvectsz[pclass_sz->index_characteristic_multipole_for_profile] = pvectsz[pclass_sz->index_l200c];
//   pvectsz[pclass_sz->index_rs] = pvectsz[pclass_sz->index_r200c];///pvectsz[pclass_sz->index_c200c];
//    double result_int;
//    double kl =  k; // (l+0.5)/l200c
//   //  printf("calling ft\n");
//    pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in p_gnfw for the pk mode computation
//    class_call_parallel(two_dim_ft_pressure_profile(kl,
//                                                    pclass_sz,pba,pvectsz,&result_int),
//                                                    pclass_sz->error_message,
//                                                    pclass_sz->error_message);
//    result = result_int;



//  // double tau_normalisation = 1.;//pvectsz[pclass_sz->index_m200m];///(4.*_PI_*pow(pvectsz[pclass_sz->index_rs],3.));
//  // tau_normalisation = 4.*_PI_*pow(pvectsz[pclass_sz->index_r200c],3);//*pvectsz[pclass_sz->index_Rho_crit];
//  // result *= tau_normalisation;


//   // pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z] = log(result);
//   pclass_sz->array_pressure_profile_ln_p_at_lnk_lnm_z[index_k][index_m_z] = result;
//   // printf("ell = %.3e z = %.3e m = %.3e lnrho = %.3e\n",ell,z,exp(lnM),log(result));

//   // printf("--------->In parallel region B12 (class_sz_tool.c) k = %.5e l = %.5e m = %.5e z = %.5e res = %.5e\n",
//   //       k,
//   //       l,
//   //       exp(lnM),
//   //       z,
//   //       result);


//   index_m_z += 1;
//      }


//      }

// free(pvectsz);
// free(pvecback);

// }

// #ifdef _OPENMP
//   tstop = omp_get_wtime();
//   if (pclass_sz->sz_verbose > 0)
//     printf("In %s: time spent in tab pressure profile parallel region (loop over k's) = %e s for thread %d\n",
//            __func__,tstop-tstart,omp_get_thread_num());
// #endif

// }
// if (abort == _TRUE_) return _FAILURE_;
// //end of parallel region
// return _SUCCESS_;

//                                       }




// Tabulate 2D Fourier transform of density profile on a [ln_ell_over_ell_char] grid
int tabulate_gas_pressure_profile_gNFW_fft(struct background * pba,
                                           struct class_sz_structure * pclass_sz){


pclass_sz->n_k_pressure_profile = pclass_sz->N_samp_fftw;
// if (pclass_sz->has_kSZ_kSZ_gal_1h + pclass_sz->has_kSZ_kSZ_gal_2h + pclass_sz->has_kSZ_kSZ_gal_3h == _FALSE_)
//   return 0;


//  // array of multipoles:
//  // this is (l+0.5)/ls
//  double ln_ell_min = log(1e-2);
//  double ln_ell_max = log(50.);
//  int n_ell = pclass_sz->array_profile_ln_PgNFW_at_lnl_over_ls_size;


//  class_alloc(pclass_sz->array_profile_ln_l_over_ls,sizeof(double *)*n_ell,pclass_sz->error_message);

//  // // array of masses:
//  // double ln_m_min = log(1e8);
//  // double ln_m_max = log(1e18);

//  //
//  // class_alloc(pclass_sz->array_profile_ln_m,sizeof(double *)*n_m,pclass_sz->error_message);
//  //

//  // // array of redshifts:
//  // double ln_1pz_min = log(1.+pclass_sz->z1SZ);
//  // double ln_1pz_max = log(1.+pclass_sz->z2SZ);


// //  class_alloc(pclass_sz->array_profile_ln_1pz,sizeof(double *)*n_z,pclass_sz->error_message);
// // int index_m_z;

// int index_l;
// for (index_l=0;
//      index_l<n_ell;
//      index_l++)
// {
//   pclass_sz->array_profile_ln_l_over_ls[index_l] = ln_ell_min
//                                               +index_l*(ln_ell_max-ln_ell_min)
//                                               /(n_ell-1.);
// }


// class_alloc(pclass_sz->array_profile_ln_PgNFW_at_lnl_over_ls,n_ell*sizeof(double *),pclass_sz->error_message);
// for (index_l=0;
//      index_l<n_ell;
//      index_l++)
// {

//   pclass_sz->array_profile_ln_PgNFW_at_lnl_over_ls[index_l] = 1e-100; // initialize with super small number

// }

// int has_ksz_bkp = pclass_sz->has_kSZ_kSZ_gal_1h;
// pclass_sz->has_kSZ_kSZ_gal_1h = _TRUE_; //pretend we need the tau_profile

//Parallelization of profile computation
/* initialize error management flag */



////FFT part
  int n_k = pclass_sz->n_k_pressure_profile; //
  class_alloc(pclass_sz->array_pressure_profile_ln_k,sizeof(double *)*n_k,pclass_sz->error_message);
  class_alloc(pclass_sz->array_pressure_profile_ln_p_at_lnk,n_k*sizeof(double *),pclass_sz->error_message);

  const int N = pclass_sz->N_samp_fftw; //precision parameter
  int ix;
  int index_k;
  double x[N], Px[N];
  double x_min = pclass_sz->x_min_gas_pressure_fftw;
  double x_max = pclass_sz->x_max_gas_pressure_fftw;


  double x_out = pclass_sz->x_outSZ;

  for (ix=0; ix<N; ix++){
        x[ix] = exp(log(x_min)+ix/(N-1.)*(log(x_max)-log(x_min)));
        if (x[ix]>x_out){
            Px[ix] = 0.;
          }
        else{

            double P0GNFW = 1.; // doesnt matter
            double alphaGNFW = pclass_sz->alphaGNFW;
            double betaGNFW = pclass_sz->betaGNFW;
            double gammaGNFW = pclass_sz->gammaGNFW;
            double c500 = pclass_sz->c500;
            Px[ix] = get_pressure_P_over_P_delta_at_x_gnfw_500c(x[ix],
                                                                P0GNFW,
                                                                alphaGNFW,
                                                                betaGNFW,
                                                                gammaGNFW,
                                                                c500,
                                                                pba,
                                                                pclass_sz
                                                                );

        }
  }

  double kp[N], Pkp[N];
  /* Compute the function
   *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
   * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
   * in this notation.  The input k-values must be logarithmically spaced.  The
   * resulting xi_l^m(r) will be evaluated at the dual r-values
   *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
  //void fftlog_ComputeXiLM(int l, int m, int N, const double k[],  const double pk[], double r[], double xi[]);
  fftlog_ComputeXiLMsloz(0, 2, N, x, Px, kp, Pkp,pclass_sz);
  // printf("fft computed.\n");

  for (index_k=0;
       index_k<n_k;
       index_k++)
  {


  // double k = exp(pclass_sz->array_pressure_profile_ln_k[index_k]); // l/ls
  // printf("calling ft\n");

  // double  result_fft = 2.*_PI_*_PI_*pwl_value_1d(N,
  //                                               kp,
  //                                               Pkp,
  //                                               k);
  //
  pclass_sz->array_pressure_profile_ln_k[index_k] = log(kp[index_k]);

  double  result_fft = 2.*_PI_*_PI_*Pkp[index_k];

  pclass_sz->array_pressure_profile_ln_p_at_lnk[index_k] = result_fft;
    // printf("ell = %.3e z = %.3e m = %.3e lnrho = %.3e\n",ell,z,exp(lnM),log(result));
    // index_m_z += 1;
  } // k loop

///FFT part done



// int abort;
// double tstart, tstop;
// abort = _FALSE_;
// /* beginning of parallel region */

// int number_of_threads= 1;
// #ifdef _OPENMP
// #pragma omp parallel
//   {
//     number_of_threads = 1; //omp_get_num_threads();
//   }
// #endif

// #pragma omp parallel \
// shared(abort,\
// pclass_sz,pba)\
// private(tstart, tstop,index_l) \
// num_threads(number_of_threads)
// {

// #ifdef _OPENMP
//   tstart = omp_get_wtime();
// #endif

// #pragma omp for schedule (dynamic)
// for (index_l=0;
//      index_l<n_ell;
//      index_l++)
// {
// #pragma omp flush(abort)
// double * pvectsz;
// double * pvecback;
// class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);
// class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
// int index_pvectsz;
// for (index_pvectsz=0;
//      index_pvectsz<pclass_sz->tsz_size;
//      index_pvectsz++){
//        pvectsz[index_pvectsz] = 0.; // set everything to 0.
//      }




//   double result;
//   double kl = exp(pclass_sz->array_profile_ln_l_over_ls[index_l]);

//   pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in p_gnfw for the pk mode computation
//   class_call_parallel(two_dim_ft_pressure_profile(kl,pclass_sz,pba,pvectsz,&result),
//                                                  pclass_sz->error_message,
//                                                  pclass_sz->error_message);


//   double fft_result = get_gas_pressure_profile_at_k(kl,pclass_sz);


//   printf("starting fft here to check ----- \n");



//   printf("wanting kl = %.3e\n",kl);
//   printf("original result = %.3e\n",result);
//   printf("fft result = %.3e\n",fft_result);
//   printf("ratio = %.3e\n",fft_result/result);

//   printf("\n\n\n\n");
//   // exit(0);



//   // pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z] = log(result);
//   if (result < 0.){
//     result = 1e-10;
//     // printf("negative result.!!\n");
//     // exit(0);
//   }
//   pclass_sz->array_profile_ln_PgNFW_at_lnl_over_ls[index_l] = log(result);
//   // printf("ell/ells = %.3e ln_pgnfw = %.3e result = %.3e\n",exp(pclass_sz->array_profile_ln_l_over_ls[index_l]),log(result),result);
//   // printf("ell/ells = %.3e ln_pgnfw = %.3e\n",exp(pclass_sz->array_profile_ln_k_over_ls[index_l]),log(result));



//   // printf("freeing pp pvectsz, pvecback\n");
//   // free(pvectsz);
//   // free(pvecback);
//   // printf("freed\n");


// }

// #ifdef _OPENMP
//   tstop = omp_get_wtime();
//   if (pclass_sz->sz_verbose > 0)
//     printf("In %s: time spent in tab profile parallel region (loop over ell's) = %e s for thread %d\n",
//            __func__,tstop-tstart,omp_get_thread_num());
// #endif

// }
// if (abort == _TRUE_) return _FAILURE_;
// //end of parallel region

return _SUCCESS_;

                                      }



// Tabulate 2D Fourier transform of density profile on a [ln_ell_over_ell_char] grid
int tabulate_gas_pressure_profile_gNFW(struct background * pba,
                                   struct class_sz_structure * pclass_sz){

// if (pclass_sz->has_kSZ_kSZ_gal_1h + pclass_sz->has_kSZ_kSZ_gal_2h + pclass_sz->has_kSZ_kSZ_gal_3h == _FALSE_)
//   return 0;


 // array of multipoles:
 // this is (l+0.5)/ls
 double ln_ell_min = log(1e-2);
 double ln_ell_max = log(50.);
 int n_ell = pclass_sz->array_profile_ln_PgNFW_at_lnl_over_ls_size;


 class_alloc(pclass_sz->array_profile_ln_l_over_ls,sizeof(double *)*n_ell,pclass_sz->error_message);

 // // array of masses:
 // double ln_m_min = log(1e8);
 // double ln_m_max = log(1e18);

 //
 // class_alloc(pclass_sz->array_profile_ln_m,sizeof(double *)*n_m,pclass_sz->error_message);
 //

 // // array of redshifts:
 // double ln_1pz_min = log(1.+pclass_sz->z1SZ);
 // double ln_1pz_max = log(1.+pclass_sz->z2SZ);


//  class_alloc(pclass_sz->array_profile_ln_1pz,sizeof(double *)*n_z,pclass_sz->error_message);
// int index_m_z;

int index_l;
for (index_l=0;
     index_l<n_ell;
     index_l++)
{
  pclass_sz->array_profile_ln_l_over_ls[index_l] = ln_ell_min
                                              +index_l*(ln_ell_max-ln_ell_min)
                                              /(n_ell-1.);
}


class_alloc(pclass_sz->array_profile_ln_PgNFW_at_lnl_over_ls,n_ell*sizeof(double *),pclass_sz->error_message);
for (index_l=0;
     index_l<n_ell;
     index_l++)
{

  pclass_sz->array_profile_ln_PgNFW_at_lnl_over_ls[index_l] = 1e-100; // initialize with super small number

}

// int has_ksz_bkp = pclass_sz->has_kSZ_kSZ_gal_1h;
// pclass_sz->has_kSZ_kSZ_gal_1h = _TRUE_; //pretend we need the tau_profile

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
pclass_sz,pba)\
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
class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
int index_pvectsz;
for (index_pvectsz=0;
     index_pvectsz<pclass_sz->tsz_size;
     index_pvectsz++){
       pvectsz[index_pvectsz] = 0.; // set everything to 0.
     }




  double result;
  double kl = exp(pclass_sz->array_profile_ln_l_over_ls[index_l]);

  pvectsz[pclass_sz->index_md] = 0; // avoid the if condition in p_gnfw for the pk mode computation
  class_call_parallel(two_dim_ft_pressure_profile(kl,pclass_sz,pba,pvectsz,&result),
                                                 pclass_sz->error_message,
                                                 pclass_sz->error_message);



  // pclass_sz->array_profile_ln_rho_at_lnk_lnM_z[index_l][index_m_z] = log(result);
  if (result < 0.){
    result = 1e-10;
    // printf("negative result.!!\n");
    // exit(0);
  }
  pclass_sz->array_profile_ln_PgNFW_at_lnl_over_ls[index_l] = log(result);
  // printf("ell/ells = %.3e ln_pgnfw = %.3e result = %.3e\n",exp(pclass_sz->array_profile_ln_l_over_ls[index_l]),log(result),result);
  // printf("ell/ells = %.3e ln_pgnfw = %.3e\n",exp(pclass_sz->array_profile_ln_k_over_ls[index_l]),log(result));



  // printf("freeing pp pvectsz, pvecback\n");
  // free(pvectsz);
  // free(pvecback);
  // printf("freed\n");


}

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
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

int two_dim_ft_pressure_profile(double kl,
                                struct class_sz_structure * pclass_sz,
                                struct background * pba,
                                double * pvectsz,
                                double * result
                          ) {




  double xin = pclass_sz->x_inSZ;
  double xout = 0.;
if (pclass_sz->pressure_profile == 4) { //for Battaglia et al 2012 pressure profile
      double rvir = pvectsz[pclass_sz->index_rVIR]; //in Mpc/h
      double r200c = pvectsz[pclass_sz->index_r200c]; //in Mpc/h
      if (pclass_sz->truncate_gas_pressure_wrt_rvir == 1){
          xout = pclass_sz->x_outSZ*rvir/r200c; // the truncation radius is in multiples of rvir
          }
      else{
        xout = pclass_sz->x_outSZ;
        }
    }// end Battaglia et al 2012 pressure profile
else{
    xout = pclass_sz->x_outSZ; // in all other cases the truncation radius is in multiples of rs=r_delta/c_delta
    }


  //GSL
  // QAWO --> START
  double delta_l = xout - xin;

  gsl_integration_workspace * w;
  gsl_integration_qawo_table * wf;

  int size_w = 3000;
  w = gsl_integration_workspace_alloc(size_w);


  double w0;
  // if (pclass_sz->pressure_profile == 4) //for Battaglia et al 2012 pressure profile
  // w0 = (pvectsz[pclass_sz->index_multipole_for_pressure_profile])/pvectsz[pclass_sz->index_l200c];
  // else
  // w0 =  (pvectsz[pclass_sz->index_multipole_for_pressure_profile]);


  // correct:
  w0 = kl; // this is (l+1/2)/ls (see eq. 2 in komatsu & seljak)

  // debugging
  // w0 = (kl+0.5)/pvectsz[pclass_sz->index_l200c];

  // printf("w0 = %.5e\n",w0);
  // exit(0);

  wf = gsl_integration_qawo_table_alloc(w0, delta_l,GSL_INTEG_SINE,50);

  struct Parameters_for_integrand_gas_pressure_profile V;
  V.pclass_sz = pclass_sz;
  V.pba = pba;
  V.pvectsz = pvectsz;
  V.kl =kl;

  void * params = &V;

  gsl_function F;
  F.function = &integrand_gas_pressure_profile;
  F.params = params;

  double eps_abs = pclass_sz->pressure_profile_epsabs;
  double eps_rel = pclass_sz->pressure_profile_epsrel;

  double result_gsl, error;
  int limit = size_w; //number of sub interval
  gsl_integration_qawo(&F,pclass_sz->x_for_pp[0],eps_abs,eps_rel,limit,w,wf,&result_gsl,&error);

  *result = result_gsl;

  gsl_integration_qawo_table_free(wf);
  gsl_integration_workspace_free(w);

  // GSL
  // QAWO --> END

  // FFTLog

  // printf("class_sz.c res = %.3e\n",result_gsl);


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
                          struct class_sz_structure * pclass_sz,
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
  double t;



  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              pclass_sz->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);

  for (i=0;i<pclass_sz->ln_k_size_for_tSZ;i++) {
    k=exp(pclass_sz->ln_k_for_tSZ[i]);
    t = 1./(1.+k);
    if (i == (pclass_sz->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;

    if (x<0.01)
      W = 1.-x*x/10.;
    else
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


    array_for_sigma[i*index_num+index_k]=t;
    array_for_sigma[i*index_num+index_y]=k*k*k*pk*W*W/(t*(1.-t));
  }




  class_call(array_spline(array_for_sigma,
                          index_num,
                          pclass_sz->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_trapzd_or_spline(array_for_sigma,
                                        index_num,
                                        pclass_sz->ln_k_size_for_tSZ,
                                        0,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);


  double sigmat = *sigma;

  //
  // for (i=0;i<pclass_sz->ln_k_size_for_tSZ;i++) {
  //   k=exp(pclass_sz->ln_k_for_tSZ[i]);
  //   if (i == (pclass_sz->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
  //   x=k*R;
  //   W=3./x/x/x*(sin(x)-x*cos(x));
  //
  //   /*
  //   class_call(spectra_pk_at_k_and_z(pba,ppm,psp,k,z,&pk,pk_ic),
  //              psp->error_message,
  //              psp->error_message);*/
  //
  //    class_call(nonlinear_pk_at_k_and_z(
  //                                      pba,
  //                                      ppm,
  //                                      pnl,
  //                                      pk_linear,
  //                                      k,
  //                                      z,
  //                                      pnl->index_pk_m,
  //                                      &pk, // number *out_pk_l
  //                                      pk_ic // array out_pk_ic_l[index_ic_ic]
  //                                    ),
  //                                    pnl->error_message,
  //                                    pnl->error_message);
  //
  //
  //   array_for_sigma[i*index_num+index_k]=k;
  //   array_for_sigma[i*index_num+index_y]=k*k*pk*W*W;
  // }
  //
  //
  //
  //
  // class_call(array_spline(array_for_sigma,
  //                         index_num,
  //                         pclass_sz->ln_k_size_for_tSZ,
  //                         index_k,
  //                         index_y,
  //                         index_ddy,
  //                         _SPLINE_EST_DERIV_,
  //                         pnl->error_message),
  //            pnl->error_message,
  //            pnl->error_message);
  //
  // class_call(array_integrate_all_trapzd_or_spline(array_for_sigma,
  //                                       index_num,
  //                                       pclass_sz->ln_k_size_for_tSZ,
  //                                       0,
  //                                       index_k,
  //                                       index_y,
  //                                       index_ddy,
  //                                       sigma,
  //                                       pnl->error_message),
  //            pnl->error_message,
  //            pnl->error_message);
  //
  // // printf("sigma t = %.5e n = %.5e\n",sigmat,*sigma);

  free(array_for_sigma);


  // *sigma = sqrt(*sigma/(2.*_PI_*_PI_));
  *sigma = sqrt(-sigmat/(2.*_PI_*_PI_));

  return _SUCCESS_;

}


//This routine computes dSigma2/dR
//at R and z

int spectra_sigma_prime(
                        struct background * pba,
                        struct primordial * ppm,
                        struct nonlinear *pnl,
                        struct class_sz_structure * pclass_sz,
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

  double k,W,x,W_prime,t;



  i=0;
  index_k=i;
  i++;
  index_y=i;
  i++;
  index_ddy=i;
  i++;
  index_num=i;

  class_alloc(array_for_sigma,
              pclass_sz->ln_k_size_for_tSZ*index_num*sizeof(double),
              pnl->error_message);

  for (i=0;i<pclass_sz->ln_k_size_for_tSZ;i++) {
    k=exp(pclass_sz->ln_k_for_tSZ[i]);
    t = 1./(1.+k);
    if (i == (pclass_sz->ln_k_size_for_tSZ-1)) k *= 0.9999999; // to prevent rounding error leading to k being bigger than maximum value
    x=k*R;
    if (x<0.01) {
      W = 1.-x*x/10.;
      W_prime = -0.2*x;
    }
    else {
    W=3./x/x/x*(sin(x)-x*cos(x));
    W_prime=3./x/x*sin(x)-9./x/x/x/x*(sin(x)-x*cos(x));
    }

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


    array_for_sigma[i*index_num+index_k]=t;//k;
    array_for_sigma[i*index_num+index_y]=k*k*k*pk*2.*k*W*W_prime/(t*(1.-t));//k*k*pk*k*2.*W*W_prime;
  }

  class_call(array_spline(array_for_sigma,
                          index_num,
                          pclass_sz->ln_k_size_for_tSZ,
                          index_k,
                          index_y,
                          index_ddy,
                          _SPLINE_EST_DERIV_,
                          pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  class_call(array_integrate_all_trapzd_or_spline(array_for_sigma,
                                        index_num,
                                        pclass_sz->ln_k_size_for_tSZ,
                                        0,
                                        index_k,
                                        index_y,
                                        index_ddy,
                                        sigma_prime,
                                        pnl->error_message),
             pnl->error_message,
             pnl->error_message);

  free(array_for_sigma);



  *sigma_prime = -*sigma_prime/(2.*_PI_*_PI_);

  return _SUCCESS_;

}



//This routine reads the tabulated
//pressure profiles,
//and stores the tabulated values.

int external_pressure_profile_init(struct precision * ppr, struct class_sz_structure * pclass_sz)
{

if (pclass_sz->pressure_profile != 0 && pclass_sz->pressure_profile != 2 )
  return 0;
  if (pclass_sz->has_tSZ_gal_1h
     +pclass_sz->has_tSZ_gal_2h
     +pclass_sz->has_sz_te_y_y
     +pclass_sz->has_sz_trispec
     +pclass_sz->has_sz_m_y_y_1h
     +pclass_sz->has_sz_m_y_y_2h
     +pclass_sz->has_sz_cov_Y_N
     +pclass_sz->has_sz_cov_Y_Y_ssc
     +pclass_sz->has_sz_cov_Y_N_next_order
     +pclass_sz->has_tSZ_lensmag_2h
     +pclass_sz->has_tSZ_lensmag_1h
     +pclass_sz->has_tSZ_gallens_1h
     +pclass_sz->has_tSZ_gallens_2h
     +pclass_sz->has_tSZ_gal_1h
     +pclass_sz->has_tSZ_gal_2h
     +pclass_sz->has_ngal_tsz_1h
     +pclass_sz->has_ngal_tsz_2h
     +pclass_sz->has_nlensmag_tsz_1h //ola
     +pclass_sz->has_nlensmag_tsz_2h
     +pclass_sz->has_tSZ_cib_1h
     +pclass_sz->has_tSZ_cib_2h
     +pclass_sz->has_tSZ_lens_1h
     +pclass_sz->has_tSZ_lens_2h
     +pclass_sz->has_tSZ_tSZ_tSZ_1halo
     +pclass_sz->has_tSZ_tSZ_tSZ_2h
     +pclass_sz->has_tSZ_tSZ_tSZ_3h
     +pclass_sz->has_kSZ_kSZ_tSZ_1h
     +pclass_sz->has_kSZ_kSZ_tSZ_2h
     +pclass_sz->has_kSZ_kSZ_tSZ_3h
     +pclass_sz->has_sz_2halo
     +pclass_sz->has_sz_ps
     +pclass_sz->has_mean_y
     +pclass_sz->has_dydz
     +pclass_sz->has_isw_tsz
     == 0)
     return 0;

  if (pclass_sz->sz_verbose > 0)
    printf("-> Using tabulated pressure profile transform\n");

  class_alloc(pclass_sz->PP_lnx,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->PP_lnI,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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
  // if(strncmp("cat ", pclass_sz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];
  if (pclass_sz->pressure_profile==0){
    class_open(process,pclass_sz->P13_file, "r",pclass_sz->error_message);
  }
  else if (pclass_sz->pressure_profile==2){
    if (pclass_sz->sz_verbose > 0)
      printf("-> Openning the pressure profile file for A10\n");
    //class_open(process,"class_sz_auxiliary_files/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_A10.txt", "r",pclass_sz->error_message);
    class_open(process,pclass_sz->A10_file, "r",pclass_sz->error_message);
    if (pclass_sz->sz_verbose > 0)
      printf("-> File Name: %s\n",pclass_sz->A10_file);
    // printf("-> File Name: %s\n",ppr->sBBN_file);

  }

    // sprintf(Filepath,
    //         "%s%s",
    //         // "%s%s%s",
    //         "cat ",
    //         // pclass_sz->path_to_class,
    //         "/class_sz_auxiliary_files/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_P13.txt");
  //   sprintf(Filepath,
  //           "%s%s",
  //           // "%s%s%s",
  //           "cat ",
  //           pclass_sz->path_to_class,
  //           "/class_sz_auxiliary_files/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_A10.txt");
  // process = popen(Filepath, "r");
  if (pclass_sz->sz_verbose > 0)
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
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnI = tmp;
      tmp = (double *)realloc(d2lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ell500's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  // printf("pclass_sz->PP_lnI[index_x] = %e\n",lnI[0]);

  status = fclose(process);
  // printf("pclass_sz->PP_lnI[index_x] = %e\n",lnI[0]);

  // fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->PP_lnx_size = n_data;
  /** Make room */
  // printf("pclass_sz->PP_lnI[index_x] = %d\n",n_data);

  class_realloc(pclass_sz->PP_lnx,
                pclass_sz->PP_lnx,
                pclass_sz->PP_lnx_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->PP_lnI,
                pclass_sz->PP_lnI,
                pclass_sz->PP_lnx_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->PP_d2lnI,
                pclass_sz->PP_d2lnI,
                pclass_sz->PP_lnx_size*sizeof(double),
                pclass_sz->error_message);


  /** Store them */
  for (index_x=0; index_x<pclass_sz->PP_lnx_size; index_x++) {
    pclass_sz->PP_lnx[index_x] = lnx[index_x];
    pclass_sz->PP_lnI[index_x] = lnI[index_x];
    pclass_sz->PP_d2lnI[index_x] = d2lnI[index_x];
    // printf("pclass_sz->PP_lnI[index_x] = %e\n",pclass_sz->PP_lnI[index_x]);

  };

  /** Release the memory used locally */
  free(lnx);
  free(lnI);
  free(d2lnI);

   if (pclass_sz->sz_verbose>1)
   printf("-> pressure profile loaded.\n");


  return _SUCCESS_;
}



//This routine reads the tabulated
//noise curve for yxy covariance,
//and stores the tabulated values.

int load_unbinned_nl_yy(struct class_sz_structure * pclass_sz)
{


// don't load  if none of the following are required:
if ( (pclass_sz->include_noise_cov_y_y != 1 )){
  // if (pclass_sz->sz_verbose>=1)
  //   printf("-> noise curve for yxy covariance not requested\n");
  return 0;
}

if (pclass_sz->sz_verbose >= 1)
  printf("-> loading the noise curve for yxy covariance\n");


  class_alloc(pclass_sz->unbinned_nl_yy_ell,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->unbinned_nl_yy_n_ell,sizeof(double *)*100,pclass_sz->error_message);
  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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
  // if(strncmp("cat ", pclass_sz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];

    sprintf(Filepath,
            "%s%s",
            "cat ",
            pclass_sz->full_path_to_noise_curve_for_y_y);
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
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  status = pclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->unbinned_nl_yy_size = n_data;
  /** Make room */

  class_realloc(pclass_sz->unbinned_nl_yy_ell,
                pclass_sz->unbinned_nl_yy_ell,
                pclass_sz->unbinned_nl_yy_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->unbinned_nl_yy_n_ell,
                pclass_sz->unbinned_nl_yy_n_ell,
                pclass_sz->unbinned_nl_yy_size*sizeof(double),
                pclass_sz->error_message);



  /** Store them */
  for (index_x=0; index_x<pclass_sz->unbinned_nl_yy_size; index_x++) {
    pclass_sz->unbinned_nl_yy_ell[index_x] = lnx[index_x];
    pclass_sz->unbinned_nl_yy_n_ell[index_x] = lnI[index_x];
    //printf("z=%.3e phig=%.3e\n",pclass_sz->unbinned_nl_yy_ell[index_x],pclass_sz->unbinned_nl_yy_n_ell[index_x]);
  };

  //exit(0);
  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}


double get_lensing_noise_at_ell(double l,
                                struct class_sz_structure * pclass_sz){
double nl_kcmb_kcmb;
if (l<pclass_sz->l_lensing_noise[0])
  nl_kcmb_kcmb = 0.;
else if (l>pclass_sz->l_lensing_noise[pclass_sz->lensing_noise_size-1])
  nl_kcmb_kcmb = 1e100;

else  nl_kcmb_kcmb = pwl_value_1d(pclass_sz->lensing_noise_size,
                              pclass_sz->l_lensing_noise,
                              pclass_sz->nl_lensing_noise,
                              l);
return nl_kcmb_kcmb;
}


double get_n5k_cl_K1_at_chi(double chi,
                                struct class_sz_structure * pclass_sz){
double r;
if (chi<pclass_sz->n5k_cl_K1_chi[0])
  r = 0.;
else if (chi>pclass_sz->n5k_cl_K1_chi[pclass_sz->n5k_cl_K1_size-1])
  r = 0;

else  r = pwl_value_1d(pclass_sz->n5k_cl_K1_size,
                              pclass_sz->n5k_cl_K1_chi,
                              pclass_sz->n5k_cl_K1_K1,
                              chi);
return r;
}

double get_n5k_z_of_chi(double chi,
                                struct class_sz_structure * pclass_sz){
double r;
if (chi<pclass_sz->n5k_z_of_chi_chi[0])
  r = 0.;
else if (chi>pclass_sz->n5k_z_of_chi_chi[pclass_sz->n5k_z_of_chi_size-1])
  r = 0;

else  r = pwl_value_1d(pclass_sz->n5k_z_of_chi_size,
                              pclass_sz->n5k_z_of_chi_chi,
                              pclass_sz->n5k_z_of_chi_z,
                              chi);
return r;
}





int load_n5k_cl_K1(struct class_sz_structure * pclass_sz)
{



if (pclass_sz->sz_verbose >= 1)
  printf("-> loading n5k Kernel K1 cl\n");


  class_alloc(pclass_sz->n5k_cl_K1_chi,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->n5k_cl_K1_K1,sizeof(double *)*100,pclass_sz->error_message);
  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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
  // if(strncmp("cat ", pclass_sz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  // char Filepath[_ARGUMENT_LENGTH_MAX_];
  // if (pclass_sz->sz_verbose >= 1)
  //   printf("-> File Name: %s\n",pclass_sz->cmb_lensing_noise_file);
  class_open(process,pclass_sz->full_path_to_n5k_gg_chi_K0, "r",pclass_sz->error_message);
  // if (pclass_sz->sz_verbose >= 1)
  //   printf("-> File Name: %s\n",pclass_sz->cmb_lensing_noise_file);

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
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->n5k_cl_K1_size = n_data;
  /** Make room */

  class_realloc(pclass_sz->n5k_cl_K1_chi,
                pclass_sz->n5k_cl_K1_chi,
                pclass_sz->n5k_cl_K1_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->n5k_cl_K1_K1,
                pclass_sz->n5k_cl_K1_K1,
                pclass_sz->n5k_cl_K1_size*sizeof(double),
                pclass_sz->error_message);



  /** Store them */
  for (index_x=0; index_x<pclass_sz->n5k_cl_K1_size; index_x++) {
    pclass_sz->n5k_cl_K1_chi[index_x] = lnx[index_x];
    pclass_sz->n5k_cl_K1_K1[index_x] = lnI[index_x];

    // printf("%.5e %.5e\n",pclass_sz->l_lensing_noise[index_x],pclass_sz->nl_lensing_noise[index_x]);

    //printf("z=%.3e phig=%.3e\n",pclass_sz->unbinned_nl_yy_ell[index_x],pclass_sz->unbinned_nl_yy_n_ell[index_x]);
  };

  // exit(0);
  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}




int load_normalized_dndz_ngal(struct class_sz_structure * pclass_sz)
{



if (pclass_sz->sz_verbose > 0)
  printf("-> [ngal] loading ngal normalized redshift kernels.\n");



  // class_alloc(pclass_sz->n5k_cl_K1_chi,sizeof(double *)*100,pclass_sz->error_message);
  // class_alloc(pclass_sz->n5k_cl_K1_K1,sizeof(double *)*100,pclass_sz->error_message);
  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

  class_alloc(pclass_sz->normalized_dndz_ngal_z,
              sizeof(double **)*pclass_sz->galaxy_samples_list_num,
              pclass_sz->error_message);
  class_alloc(pclass_sz->normalized_dndz_ngal_phig,
              sizeof(double **)*pclass_sz->galaxy_samples_list_num,
              pclass_sz->error_message);
  class_alloc(pclass_sz->normalized_dndz_ngal_size,
              sizeof(int *)*pclass_sz->galaxy_samples_list_num,
              pclass_sz->error_message);

  int index_g;
  for (index_g=0;index_g<pclass_sz->galaxy_samples_list_num;index_g++){

  class_alloc(pclass_sz->normalized_dndz_ngal_z[index_g],sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->normalized_dndz_ngal_phig[index_g],sizeof(double *)*100,pclass_sz->error_message);
}

if (pclass_sz->sz_verbose>0)
  printf("-> [ngal] galaxy redshift kernel arrays allocated. Starting filling them.\n");



for (index_g=0;index_g<pclass_sz->galaxy_samples_list_num;index_g++){
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
  // if(strncmp("cat ", pclass_sz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  // char Filepath[_ARGUMENT_LENGTH_MAX_];
  // if (pclass_sz->sz_verbose >= 1)
  //   printf("-> File Name: %s\n",pclass_sz->cmb_lensing_noise_file);
  char Filepath[_ARGUMENT_LENGTH_MAX_];
  sprintf(Filepath,"%s%d%s",
          pclass_sz->full_path_and_prefix_to_dndz_ngal,
          pclass_sz->galaxy_samples_list[index_g],
          ".txt");
  if (pclass_sz->sz_verbose >= 1)
  printf("-> [ngal] kernel file name for sample %d: %s\n",index_g,Filepath);
  class_open(process,Filepath, "r",pclass_sz->error_message);
  if (pclass_sz->sz_verbose >= 1)
    printf("-> [ngal] kernel file for sample %d opened successfully\n");

  // exit(0);

  //int il = 0;
  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {

    //il++;
    sscanf(line, "%lf %lf ", &this_lnx, &this_lnI);

    // printf("%lf %lf\n",this_lnx,this_lnI);

    /* Standard technique in C:
     /*if too many data, double the size of the vectors */
    /* (it is faster and safer that reallocating every new line) */
    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
                 "Error allocating memory to read the file.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
                 "Error allocating memory to read the file.\n");
      lnI = tmp;
    };
    /* Store */
    lnx[n_data]   = this_lnx;
    lnI[n_data]   = this_lnI;

    n_data++;
    /* Check ascending order of the k's */
    if(n_data>1) {
      class_test(lnx[n_data-1] <= lnx[n_data-2],
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->normalized_dndz_ngal_size[index_g] = n_data;

// }
//
// printf("reallocating1\n");
// printf("reallocating\n");
// printf("reallocating\n");
// printf("reallocating3\n");

// for (index_g=0;index_g<pclass_sz->galaxy_samples_list_num;index_g++){
  /** Make room */
  if (pclass_sz->sz_verbose >= 1)
    printf("-> [ngal]reallocating %d with size %d\n",index_g,pclass_sz->normalized_dndz_ngal_size[index_g]);
  class_realloc(pclass_sz->normalized_dndz_ngal_z[index_g],
                pclass_sz->normalized_dndz_ngal_z[index_g],
                pclass_sz->normalized_dndz_ngal_size[index_g]*sizeof(double),
                pclass_sz->error_message);
  // printf("reallocating ngal z done\n");

  class_realloc(pclass_sz->normalized_dndz_ngal_phig[index_g],
                pclass_sz->normalized_dndz_ngal_phig[index_g],
                pclass_sz->normalized_dndz_ngal_size[index_g]*sizeof(double),
                pclass_sz->error_message);
  // printf("reallocating done\n");
// }
// exit(0);


  // /** Store them */
  for (index_x=0; index_x<pclass_sz->normalized_dndz_ngal_size[index_g]; index_x++) {
    pclass_sz->normalized_dndz_ngal_z[index_g][index_x] = lnx[index_x];
    pclass_sz->normalized_dndz_ngal_phig[index_g][index_x] = lnI[index_x];


  if (pclass_sz->sz_verbose >= 3)
    printf("-> [ngal] kernel sample %d z = %.3e phig = %.3e\n",
           index_g,
           pclass_sz->normalized_dndz_ngal_z[index_g][index_x],
           pclass_sz->normalized_dndz_ngal_phig[index_g][index_x]);

    // printf("%.5e %.5e\n",pclass_sz->l_lensing_noise[index_x],pclass_sz->nl_lensing_noise[index_x]);

    //printf("z=%.3e phig=%.3e\n",pclass_sz->unbinned_nl_yy_ell[index_x],pclass_sz->unbinned_nl_yy_n_ell[index_x]);
  }

  // exit(0);
  /** Release the memory used locally */
  free(lnx);
  free(lnI);
}

  return _SUCCESS_;
}






int load_n5k_z_of_chi(struct class_sz_structure * pclass_sz)
{



if (pclass_sz->sz_verbose >= 1)
  printf("-> loading n5k z_of_chi\n");


  class_alloc(pclass_sz->n5k_z_of_chi_chi,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->n5k_z_of_chi_z,sizeof(double *)*100,pclass_sz->error_message);
  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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




  class_open(process,pclass_sz->full_path_to_n5k_z_chi, "r",pclass_sz->error_message);


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
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->n5k_z_of_chi_size = n_data;
  /** Make room */

  class_realloc(pclass_sz->n5k_z_of_chi_chi,
                pclass_sz->n5k_z_of_chi_chi,
                pclass_sz->n5k_z_of_chi_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->n5k_z_of_chi_z,
                pclass_sz->n5k_z_of_chi_z,
                pclass_sz->n5k_z_of_chi_size*sizeof(double),
                pclass_sz->error_message);



  /** Store them */
  for (index_x=0; index_x<pclass_sz->n5k_cl_K1_size; index_x++) {
    pclass_sz->n5k_z_of_chi_z[index_x] = lnx[index_x];
    pclass_sz->n5k_z_of_chi_chi[index_x] = lnI[index_x];

    // printf("%.5e %.5e\n",pclass_sz->l_lensing_noise[index_x],pclass_sz->nl_lensing_noise[index_x]);

    //printf("z=%.3e phig=%.3e\n",pclass_sz->unbinned_nl_yy_ell[index_x],pclass_sz->unbinned_nl_yy_n_ell[index_x]);
  };

  // exit(0);
  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}




int load_nl_lensing_noise(struct class_sz_structure * pclass_sz)
{



if (pclass_sz->sz_verbose >= 1)
  printf("-> loading the noise curve for CMB lensing\n");


  class_alloc(pclass_sz->nl_lensing_noise,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->l_lensing_noise,sizeof(double *)*100,pclass_sz->error_message);
  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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
  // if(strncmp("cat ", pclass_sz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];
  if (pclass_sz->sz_verbose >= 1)
    printf("-> File Name: %s\n",pclass_sz->cmb_lensing_noise_file);
  class_open(process,pclass_sz->cmb_lensing_noise_file, "r",pclass_sz->error_message);
  if (pclass_sz->sz_verbose >= 1)
    printf("-> File Name: %s\n",pclass_sz->cmb_lensing_noise_file);

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
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->lensing_noise_size = n_data;
  /** Make room */

  class_realloc(pclass_sz->nl_lensing_noise,
                pclass_sz->nl_lensing_noise,
                pclass_sz->lensing_noise_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->l_lensing_noise,
                pclass_sz->l_lensing_noise,
                pclass_sz->lensing_noise_size*sizeof(double),
                pclass_sz->error_message);



  /** Store them */
  for (index_x=0; index_x<pclass_sz->lensing_noise_size; index_x++) {
    pclass_sz->l_lensing_noise[index_x] = lnx[index_x];
    pclass_sz->nl_lensing_noise[index_x] = lnI[index_x];

    // printf("%.5e %.5e\n",pclass_sz->l_lensing_noise[index_x],pclass_sz->nl_lensing_noise[index_x]);

    //printf("z=%.3e phig=%.3e\n",pclass_sz->unbinned_nl_yy_ell[index_x],pclass_sz->unbinned_nl_yy_n_ell[index_x]);
  };

  // exit(0);
  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}



//This routine reads the tabulated
//noise curve for t-t covariance,
//and stores the tabulated values.

int load_unbinned_nl_tt(struct class_sz_structure * pclass_sz)
{

if (pclass_sz->sz_verbose >= 1)
  printf("-> loading the noise curve for TT covariance\n");


  class_alloc(pclass_sz->unbinned_nl_tt_ell,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->unbinned_nl_tt_n_ell,sizeof(double *)*100,pclass_sz->error_message);
  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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
  // if(strncmp("cat ", pclass_sz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];

    sprintf(Filepath,
            "%s%s",
            "cat ",
            pclass_sz->full_path_to_noise_curve_for_t_t);
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
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  status = pclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->unbinned_nl_tt_size = n_data;
  /** Make room */

  class_realloc(pclass_sz->unbinned_nl_tt_ell,
                pclass_sz->unbinned_nl_tt_ell,
                pclass_sz->unbinned_nl_tt_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->unbinned_nl_tt_n_ell,
                pclass_sz->unbinned_nl_tt_n_ell,
                pclass_sz->unbinned_nl_tt_size*sizeof(double),
                pclass_sz->error_message);



  /** Store them */
  for (index_x=0; index_x<pclass_sz->unbinned_nl_tt_size; index_x++) {
    pclass_sz->unbinned_nl_tt_ell[index_x] = lnx[index_x];
    pclass_sz->unbinned_nl_tt_n_ell[index_x] = lnI[index_x];
    //printf("z=%.3e phig=%.3e\n",pclass_sz->unbinned_nl_yy_ell[index_x],pclass_sz->unbinned_nl_yy_n_ell[index_x]);
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

int load_normalized_source_dndz(struct class_sz_structure * pclass_sz)
{

// don't load the unwise  dndz  if none of the following are required:
// all quantities requiring galaxy or lensmag need that:
if (
     (pclass_sz->has_gal_gallens_1h != _TRUE_ )
    && (pclass_sz->has_gal_gallens_2h != _TRUE_ )
    && (pclass_sz->has_IA_gal_2h != _TRUE_ )
    && (pclass_sz->has_ngal_gallens_1h != _TRUE_ )
    && (pclass_sz->has_ngal_gallens_2h != _TRUE_ )
    && (pclass_sz->has_nlensmag_gallens_1h != _TRUE_ )
    && (pclass_sz->has_nlensmag_gallens_2h != _TRUE_ )
    && (pclass_sz->has_ngal_IA_2h != _TRUE_ )
    && (pclass_sz->has_gallens_gallens_2h != _TRUE_ )
    && (pclass_sz->has_gallens_gallens_2h != _TRUE_ )
    && (pclass_sz->has_gallens_cib_1h != _TRUE_ )
    && (pclass_sz->has_gallens_cib_2h != _TRUE_ )
    && (pclass_sz->has_gallens_lens_1h != _TRUE_ )
    && (pclass_sz->has_gallens_lens_2h != _TRUE_ )
    && (pclass_sz->has_gallens_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_gallens_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_tSZ_gallens_1h != _TRUE_ )
    && (pclass_sz->has_tSZ_gallens_2h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gallens_1h_fft != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gallens_2h_fft != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gallens_3h_fft != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gallens_hf != _TRUE_ )
)
  return 0;

if (pclass_sz->sz_verbose>=1){

    printf("-> Loading source dndz file\n");
}

  class_alloc(pclass_sz->normalized_source_dndz_z,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->normalized_source_dndz_phig,sizeof(double *)*100,pclass_sz->error_message);

  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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
  // if(strncmp("cat ", pclass_sz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];

    if (pclass_sz->sz_verbose > 0){
      printf("-> Openning the dndz file for source galaxies\n");
      printf("-> File Name: %s\n",pclass_sz->full_path_to_source_dndz_gal);
      // printf("-> File Name: %s\n",pclass_sz->UNWISE_fdndz_file);
      // printf("-> File Name: %s\n",pclass_sz->A10_file);
    }
  class_open(process,pclass_sz->full_path_to_source_dndz_gal, "r",pclass_sz->error_message);
  if (pclass_sz->sz_verbose > 0)
    printf("-> File opened successfully\n");


  // process = popen(Filepath, "r");

  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {

    sscanf(line, "%lf %lf ", &this_lnx, &this_lnI);



    /* Standard technique in C:
     /*if too many data, double the size of the vectors */
    /* (it is faster and safer that reallocating every new line) */
    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  // if (pclass_sz->galaxy_sample == 2){
  // // if (pclass_sz->galaxy_sample == 2 || pclass_sz->galaxy_sample == 0 ){
  //   status = pclose(process);
  // }
  // else{
    status = fclose(process);
  //}
  // status = pclose(process);

  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->normalized_source_dndz_size = n_data;
  /** Make room */

  class_realloc(pclass_sz->normalized_source_dndz_z,
                pclass_sz->normalized_source_dndz_z,
                pclass_sz->normalized_source_dndz_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->normalized_source_dndz_phig,
                pclass_sz->normalized_source_dndz_phig,
                pclass_sz->normalized_source_dndz_size*sizeof(double),
                pclass_sz->error_message);



  /** Store them */
  // printf("normalized_source_dndz_size = %d\n",pclass_sz->normalized_source_dndz_size);
  for (index_x=0; index_x<pclass_sz->normalized_source_dndz_size; index_x++) {
    pclass_sz->normalized_source_dndz_z[index_x] = lnx[index_x];
    pclass_sz->normalized_source_dndz_phig[index_x] = lnI[index_x];
    // printf("z=%.3e phig=%.3e\n",pclass_sz->normalized_source_dndz_z[index_x],pclass_sz->normalized_source_dndz_z[index_x]);
  };

  // exit(0);

  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}



//This routine reads the tabulated
//Snu(z,nu) for cib computations
int load_cib_Snu(
                      struct class_sz_structure * pclass_sz
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
  int index_k;
  int index_z;


  class_alloc(pclass_sz->cib_Snu_z,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->cib_Snu_nu,sizeof(double *)*100,pclass_sz->error_message);



  n_data = 0;
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));



class_open(process,pclass_sz->cib_Snu_file_z, "r",pclass_sz->error_message);

  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
                 "Error allocating memory to read the snu files.\n");
      lnx = tmp;
    };


    /* Store */
    lnx[n_data]   = this_lnx;
    n_data++;
  }

  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  pclass_sz->cib_Snu_z_size = n_data;

  class_realloc(pclass_sz->cib_Snu_z,
                pclass_sz->cib_Snu_z,
                pclass_sz->cib_Snu_z_size*sizeof(double),
                pclass_sz->error_message);


  /** Store them */
  for (index_x=0; index_x<pclass_sz->cib_Snu_z_size; index_x++) {
    pclass_sz->cib_Snu_z[index_x] = lnx[index_x];
  };


  //Masses

  n_data = 0;
  n_data_guess = 100;
  lnx   = (double *)malloc(n_data_guess*sizeof(double));

class_open(process,pclass_sz->cib_Snu_file_nu, "r",pclass_sz->error_message);
  // class_open(process,"class_sz_auxiliary_files/filtered_snu_planck_nu.txt", "r",pclass_sz->error_message);

  // printf("-> %s\n",Filepath);

  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
                 "Error allocating memory to read the Snu files.\n");
      lnx = tmp;
    };


    /* Store */
    lnx[n_data]   = this_lnx;
    n_data++;
  }

  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  pclass_sz->cib_Snu_nu_size = n_data;

  class_realloc(pclass_sz->cib_Snu_nu,
                pclass_sz->cib_Snu_nu,
                pclass_sz->cib_Snu_nu_size*sizeof(double),
                pclass_sz->error_message);


  /** Store them */
  for (index_x=0; index_x<pclass_sz->cib_Snu_nu_size; index_x++) {
    pclass_sz->cib_Snu_nu[index_x] = lnx[index_x];
  };


  /** Release the memory used locally */
  free(lnx);

  //Read pk

  class_alloc(pclass_sz->cib_Snu_snu,
              sizeof(double *)*pclass_sz->cib_Snu_z_size*pclass_sz->cib_Snu_nu_size,
              pclass_sz->error_message);

  class_alloc(logC,
              pclass_sz->cib_Snu_z_size*sizeof(double *),
              pclass_sz->error_message);


  for (index_z=0;
       index_z<pclass_sz->cib_Snu_z_size;
       index_z++)
  {
    class_alloc(logC[index_z],
                pclass_sz->cib_Snu_nu_size*sizeof(double),
                pclass_sz->error_message);
  }


class_open(process,pclass_sz->cib_Snu_file_snu, "r",pclass_sz->error_message);
  // class_open(process,"class_sz_auxiliary_files/filtered_snu_planck_90_100_143_217_353_545_857.txt", "r",pclass_sz->error_message);


  int z =0;
  while (fgets(line, sizeof(line)-1, process) != NULL) {
    // printf("%s", line);
    // exit(0);
    int i=0;
    char *err, *p = line;
    double val;
    while (*p) {
      val = strtod(p, &err);
      logC[z][i] = log(val); 
      
      if (val<0) printf("z %d i %d %.3e %.3e\n",z,i,val,logC[z][i]);
      p = err + 1;
      i+=1;
    }
    // printf("\n %d \n",z);
    z+=1;
  }

  // printf("storing");
  int index = 0;
  for (index_z=0;
       index_z<pclass_sz->cib_Snu_z_size;
       index_z++){
    for (index_k=0;
         index_k<pclass_sz->cib_Snu_nu_size;
         index_k++){

      pclass_sz->cib_Snu_snu[index] = logC[index_z][index_k];
      // printf("pk %.5e\n", logC[index_z][index_k]);//pclass_sz->n5k_pk_pk[index]);
      index += 1;
    }
  }

  status = fclose(process);


  for (index_z=0;
       index_z<pclass_sz->cib_Snu_z_size;
       index_z++){
         free(logC[index_z]);
       }
  free(logC);
if (pclass_sz->sz_verbose>=1){
  printf("cib Snu loaded with %d z and %d nu\n",pclass_sz->cib_Snu_z_size,pclass_sz->cib_Snu_nu_size);
}
  return _SUCCESS_;
}


//This routine reads the tabulated
//dndz galaxy counts,
//and stores the tabulated values.

int load_normalized_dndz(struct class_sz_structure * pclass_sz)
{

// don't load the unwise  dndz  if none of the following are required:
// all quantities requiring galaxy or lensmag need that:
if (   (pclass_sz->has_tSZ_gal_1h != _TRUE_ )
    && (pclass_sz->has_tSZ_gal_2h != _TRUE_ )
    && (pclass_sz->has_IA_gal_2h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_1h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_1h_fft != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_2h_fft != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_3h_fft != _TRUE_ )
    && (pclass_sz->has_gal_gal_lens_1h_fft != _TRUE_ )
    && (pclass_sz->has_gal_gal_lens_2h_fft != _TRUE_ )
    && (pclass_sz->has_gal_gal_lens_3h_fft != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_2h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_3h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_hf != _TRUE_ )
    && (pclass_sz->has_bk_ttg_at_z_1h != _TRUE_ )
    && (pclass_sz->has_bk_ttg_at_z_2h != _TRUE_ )
    && (pclass_sz->has_bk_ttg_at_z_3h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_lensmag_1halo != _TRUE_ )
    && (pclass_sz->has_gal_gal_1h != _TRUE_ )
    && (pclass_sz->has_tau_gal_1h != _TRUE_ )
    && (pclass_sz->has_tau_gal_2h != _TRUE_ )
    && (pclass_sz->has_gal_lens_1h != _TRUE_ )
    && (pclass_sz->has_gal_lens_2h != _TRUE_ )
    && (pclass_sz->has_gal_cib_1h != _TRUE_ )
    && (pclass_sz->has_gal_cib_2h != _TRUE_ )
    && (pclass_sz->has_gal_lens_hf != _TRUE_ )
    && (pclass_sz->has_gal_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_gal_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_gal_gallens_1h != _TRUE_ )
    && (pclass_sz->has_gal_gallens_2h != _TRUE_ )
    && (pclass_sz->has_gal_lensmag_hf != _TRUE_ )
    && (pclass_sz->has_tSZ_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_tSZ_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_lensmag_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_lensmag_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_lensmag_lensmag_hf != _TRUE_ )
    && (pclass_sz->has_gallens_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_gallens_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_lens_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_lens_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_lens_lensmag_hf != _TRUE_ )
    && (pclass_sz->has_gal_gal_2h != _TRUE_ )
    && (pclass_sz->has_gal_gal_hf != _TRUE_ ))
  return 0;

if (pclass_sz->sz_verbose>=1){
    if (pclass_sz->galaxy_sample == 0)
    printf("-> Loading dndz WIxSC\n");
    if (pclass_sz->galaxy_sample == 1)
    printf("-> Loading dndz unwise\n");
    if (pclass_sz->galaxy_sample == 2)
    printf("-> Loading dndz file\n");
    }

  class_alloc(pclass_sz->normalized_dndz_z,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->normalized_dndz_phig,sizeof(double *)*100,pclass_sz->error_message);

  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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
  // if(strncmp("cat ", pclass_sz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];

  //unwise
  if (pclass_sz->galaxy_sample == 1){
    if (pclass_sz->sz_verbose > 0){
      printf("-> Openning the dndz file for unWISE galaxies\n");
      printf("-> File Name: %s\n",pclass_sz->UNWISE_dndz_file);
      // printf("-> File Name: %s\n",pclass_sz->UNWISE_fdndz_file);
      // printf("-> File Name: %s\n",pclass_sz->A10_file);
    }
  class_open(process,pclass_sz->UNWISE_dndz_file, "r",pclass_sz->error_message);
  if (pclass_sz->sz_verbose > 0)
    printf("-> File opened successfully\n");
    // sprintf(Filepath,
    //         "%s%s",
    //         "cat ",
    //         //pclass_sz->path_to_class,
    //         "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/UNWISE_galaxy_ditributions/normalised_dndz.txt");
            }

  else if (pclass_sz->galaxy_sample == 0){
    if (pclass_sz->sz_verbose > 0){
      printf("-> Openning the dndz file for WISC3 galaxies\n");
      printf("-> File Name: %s\n",pclass_sz->WISC3_dndz_file);
      // printf("-> File Name: %s\n",pclass_sz->UNWISE_fdndz_file);
      // printf("-> File Name: %s\n",pclass_sz->A10_file);
    }
  class_open(process,pclass_sz->WISC3_dndz_file, "r",pclass_sz->error_message);
    if (pclass_sz->sz_verbose > 0)
      printf("-> File opened successfully\n");
 //  sprintf(Filepath,
 //          "%s%s",
 //          "cat ",
 //          //pclass_sz->path_to_class,
 //          "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/run_scripts/yxg/data/dndz/WISC_bin3.txt");
 // process = popen(Filepath, "r");
        }

  else if (pclass_sz->galaxy_sample == 2){
    if (pclass_sz->sz_verbose > 0){
      printf("-> Openning the dndz file for galaxies\n");
      printf("-> File Name: %s\n",pclass_sz->full_path_to_dndz_gal);
      // printf("-> File Name: %s\n",pclass_sz->UNWISE_fdndz_file);
      // printf("-> File Name: %s\n",pclass_sz->A10_file);
    }

  // char Filepath[_ARGUMENT_LENGTH_MAX_];
  sprintf(Filepath,"%s",pclass_sz->full_path_to_dndz_gal);
  class_open(process,pclass_sz->full_path_to_dndz_gal, "r",pclass_sz->error_message);
  // class_open(process,Filepath, "r",pclass_sz->error_message);

  if (pclass_sz->sz_verbose > 0)
    printf("-> File opened successfully\n");
  // sprintf(Filepath,
  //         "%s%s",
  //         "cat ",
  //         pclass_sz->full_path_to_dndz_gal);
  //         //"/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/run_scripts/yxg/data/dndz/unwise_red.txt");
  // process = popen(Filepath, "r");
        }

// exit(0);

  // process = popen(Filepath, "r");

  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {

    // unWISE load and read column depending on the requested color
      if (pclass_sz->galaxy_sample == 1){
    sscanf(line, "%lf %lf %lf %lf", &this_lnx, &this_lnI, &this_lnJ, &this_lnK);
    //sscanf(line, "%lf %lf", &this_lnx, &this_lnI);

    // red
    if (pclass_sz->unwise_galaxy_sample_id == 0)
    this_lnI = this_lnK;

    // green
    if (pclass_sz->unwise_galaxy_sample_id == 1 || pclass_sz->unwise_galaxy_sample_id == 2)
    this_lnI = this_lnJ;

    // blue
    //if (pclass_sz->unwise_galaxy_sample_id == 3)
    //this_lnI = this_lnI;

    // printf("lnx = %e\n",this_lnI);
                                    }

  // WIxSC and "other": just two columns files
  else if ((pclass_sz->galaxy_sample == 0) || (pclass_sz->galaxy_sample == 2)){

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
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  // if (pclass_sz->galaxy_sample == 2){
  // // if (pclass_sz->galaxy_sample == 2 || pclass_sz->galaxy_sample == 0 ){
  //   status = pclose(process);
  // }
  // else{
    status = fclose(process);
  //}
  // status = pclose(process);

  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->normalized_dndz_size = n_data;
  /** Make room */

  class_realloc(pclass_sz->normalized_dndz_z,
                pclass_sz->normalized_dndz_z,
                pclass_sz->normalized_dndz_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->normalized_dndz_phig,
                pclass_sz->normalized_dndz_phig,
                pclass_sz->normalized_dndz_size*sizeof(double),
                pclass_sz->error_message);



  /** Store them */
  for (index_x=0; index_x<pclass_sz->normalized_dndz_size; index_x++) {
    pclass_sz->normalized_dndz_z[index_x] = lnx[index_x];
    pclass_sz->normalized_dndz_phig[index_x] = lnI[index_x];
    
    if (pclass_sz->sz_verbose > 10) {
    printf("load_normalized_dndz z=%.3e phig=%.3e\n",pclass_sz->normalized_dndz_z[index_x],pclass_sz->normalized_dndz_z[index_x]);
    }
  };

  // exit(0);


  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}



int load_normalized_fdndz(struct class_sz_structure * pclass_sz)
{

// don't load the unwise  dndz  if none of the following are required:
if (   (pclass_sz->has_tSZ_gal_1h != _TRUE_ )
    && (pclass_sz->has_tSZ_gal_2h != _TRUE_ )
    && (pclass_sz->has_IA_gal_2h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_1h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_1h_fft != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_2h_fft != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_3h_fft != _TRUE_ )
    && (pclass_sz->has_gal_gal_lens_1h_fft != _TRUE_ )
    && (pclass_sz->has_gal_gal_lens_2h_fft != _TRUE_ )
    && (pclass_sz->has_gal_gal_lens_3h_fft != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_2h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_3h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_hf != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_lensmag_1halo != _TRUE_ )
    && (pclass_sz->has_gal_gal_1h != _TRUE_ )
    && (pclass_sz->has_tau_gal_1h != _TRUE_ )
    && (pclass_sz->has_tau_gal_2h != _TRUE_ )
    && (pclass_sz->has_gal_lens_1h != _TRUE_ )
    && (pclass_sz->has_gal_lens_2h != _TRUE_ )
    && (pclass_sz->has_gal_cib_1h != _TRUE_ )
    && (pclass_sz->has_gal_cib_2h != _TRUE_ )
    && (pclass_sz->has_gal_lens_hf != _TRUE_ )
    && (pclass_sz->has_gal_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_gal_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_gal_lensmag_hf != _TRUE_ )
    && (pclass_sz->has_tSZ_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_tSZ_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_lensmag_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_lensmag_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_lensmag_lensmag_hf != _TRUE_ )
    && (pclass_sz->has_gallens_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_gallens_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_lens_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_lens_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_lens_lensmag_hf != _TRUE_ )
    && (pclass_sz->has_gal_gal_2h != _TRUE_ )
    && (pclass_sz->has_gal_gal_hf != _TRUE_ ))
  return 0;

if ((pclass_sz->galaxy_sample == 0) || (pclass_sz->galaxy_sample == 2))
  return 0;

if (pclass_sz->sz_verbose >= 1)
printf("-> Loading fdndz unwise\n");

  class_alloc(pclass_sz->normalized_fdndz_z,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->normalized_fdndz_phig,sizeof(double *)*100,pclass_sz->error_message);

  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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
  // if(strncmp("cat ", pclass_sz->command, 4) == 0)
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
  //           //pclass_sz->path_to_class,
  //           "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/UNWISE_galaxy_ditributions/normalised_fdndz.txt");

  class_open(process,pclass_sz->UNWISE_fdndz_file, "r",pclass_sz->error_message);

  // process = popen(Filepath, "r");

  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {

    // unWISE load and read column depending on the requested color

    sscanf(line, "%lf %lf %lf %lf", &this_lnx, &this_lnI, &this_lnJ, &this_lnK);
    //sscanf(line, "%lf %lf", &this_lnx, &this_lnI);

    // red
    if (pclass_sz->unwise_galaxy_sample_id == 0)
    this_lnI = this_lnK;

    // green
    if (pclass_sz->unwise_galaxy_sample_id == 1 || pclass_sz->unwise_galaxy_sample_id == 2)
    this_lnI = this_lnJ;

    // blue
    //if (pclass_sz->unwise_galaxy_sample_id == 3)
    //this_lnI = this_lnI;

    //printf("lnx = %e\n",this_lnx);





    /* Standard technique in C:
     /*if too many data, double the size of the vectors */
    /* (it is faster and safer that reallocating every new line) */
    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->normalized_fdndz_size = n_data;
  /** Make room */

  class_realloc(pclass_sz->normalized_fdndz_z,
                pclass_sz->normalized_fdndz_z,
                pclass_sz->normalized_fdndz_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->normalized_fdndz_phig,
                pclass_sz->normalized_fdndz_phig,
                pclass_sz->normalized_fdndz_size*sizeof(double),
                pclass_sz->error_message);



  /** Store them */
  for (index_x=0; index_x<pclass_sz->normalized_fdndz_size; index_x++) {
    pclass_sz->normalized_fdndz_z[index_x] = lnx[index_x];
    pclass_sz->normalized_fdndz_phig[index_x] = lnI[index_x];
    //print("z=%.3e phig=%.3e\n",pclass_sz->normalized_dndz_z[index_x])
  };

  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}

// unwise dndz deduced from cross-match with spectroscopic surveys
int load_normalized_cosmos_dndz(struct class_sz_structure * pclass_sz)
{
  // printf("gs = %d\n",pclass_sz->galaxy_sample);
// don't load the unwise  dndz  if none of the following are required:
if (   (pclass_sz->has_tSZ_gal_1h != _TRUE_ )
    && (pclass_sz->has_tSZ_gal_2h != _TRUE_ )
    && (pclass_sz->has_IA_gal_2h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_1h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_1h_fft != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_2h_fft != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_3h_fft != _TRUE_ )
    && (pclass_sz->has_gal_gal_lens_1h_fft != _TRUE_ )
    && (pclass_sz->has_gal_gal_lens_2h_fft != _TRUE_ )
    && (pclass_sz->has_gal_gal_lens_3h_fft != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_2h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_3h != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_gal_hf != _TRUE_ )
    && (pclass_sz->has_kSZ_kSZ_lensmag_1halo != _TRUE_ )
    && (pclass_sz->has_gal_gal_1h != _TRUE_ )
    && (pclass_sz->has_tau_gal_1h != _TRUE_ )
    && (pclass_sz->has_tau_gal_2h != _TRUE_ )
    && (pclass_sz->has_gal_lens_1h != _TRUE_ )
    && (pclass_sz->has_gal_lens_2h != _TRUE_ )
    && (pclass_sz->has_gal_cib_1h != _TRUE_ )
    && (pclass_sz->has_gal_cib_2h != _TRUE_ )
    && (pclass_sz->has_gal_lens_hf != _TRUE_ )
    && (pclass_sz->has_gal_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_gal_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_gal_lensmag_hf != _TRUE_ )
    && (pclass_sz->has_tSZ_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_tSZ_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_lensmag_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_lensmag_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_lensmag_lensmag_hf != _TRUE_ )
    && (pclass_sz->has_gallens_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_gallens_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_lens_lensmag_1h != _TRUE_ )
    && (pclass_sz->has_lens_lensmag_2h != _TRUE_ )
    && (pclass_sz->has_lens_lensmag_hf != _TRUE_ )
    && (pclass_sz->has_gal_gal_2h != _TRUE_ )
    && (pclass_sz->has_gal_gal_hf != _TRUE_ ))
  return 0;



if ((pclass_sz->galaxy_sample == 0) || (pclass_sz->galaxy_sample == 2))
  return 0;

if (pclass_sz->sz_verbose >= 1)
printf("-> Loading cosmos dndz unwise\n");

  class_alloc(pclass_sz->normalized_cosmos_dndz_z,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->normalized_cosmos_dndz_phig,sizeof(double *)*100,pclass_sz->error_message);

  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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
  // if(strncmp("cat ", pclass_sz->command, 4) == 0)
  // {
  // sprintf(arguments, " ");
  // }

  /** 2. Launch the command and retrieve the output */
  /* Launch the process */
  char Filepath[_ARGUMENT_LENGTH_MAX_];

  //unwise
    class_open(process,pclass_sz->UNWISE_cosmos_dndz_file, "r",pclass_sz->error_message);
    if (pclass_sz->sz_verbose > 0){
      printf("-> Openning the cosmos dndz file for unWISE galaxies\n");
      printf("-> File Name: %s\n",pclass_sz->UNWISE_cosmos_dndz_file);
      // printf("-> File Name: %s\n",pclass_sz->UNWISE_fdndz_file);
      // printf("-> File Name: %s\n",pclass_sz->A10_file);
    }

    // sprintf(Filepath,
    //         "%s%s",
    //         "cat ",
    //         //pclass_sz->path_to_class,
    //         "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/UNWISE_galaxy_ditributions/normalised_dndz_cosmos.txt");
    //


  // process = popen(Filepath, "r");

  /* Read output and store it */
  while (fgets(line, sizeof(line)-1, process) != NULL) {

    // unWISE load and read column depending on the requested color

    sscanf(line, "%lf %lf %lf %lf", &this_lnx, &this_lnI, &this_lnJ, &this_lnK);
    // sscanf(line, "%lf %lf", &this_lnx, &this_lnI);

    // red
    if (pclass_sz->unwise_galaxy_sample_id == 0)
    this_lnI = this_lnK;

    // green
    if (pclass_sz->unwise_galaxy_sample_id == 1 || pclass_sz->unwise_galaxy_sample_id == 2)
    this_lnI = this_lnJ;

    // blue
    //if (pclass_sz->unwise_galaxy_sample_id == 3)
    //this_lnI = this_lnI;

    // printf("lnx = %e\n",this_lnI);





    /* Standard technique in C:
     /*if too many data, double the size of the vectors */
    /* (it is faster and safer that reallocating every new line) */
    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(lnx,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  // printf("closing the process\n");

  /* Close the process */
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->normalized_cosmos_dndz_size = n_data;
  /** Make room */

  class_realloc(pclass_sz->normalized_cosmos_dndz_z,
                pclass_sz->normalized_cosmos_dndz_z,
                pclass_sz->normalized_cosmos_dndz_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->normalized_cosmos_dndz_phig,
                pclass_sz->normalized_cosmos_dndz_phig,
                pclass_sz->normalized_cosmos_dndz_size*sizeof(double),
                pclass_sz->error_message);


  // printf("fillling the arrays of dim = %d\n",pclass_sz->normalized_cosmos_dndz_size);
  /** Store them */
  for (index_x=0; index_x<pclass_sz->normalized_cosmos_dndz_size; index_x++) {
    pclass_sz->normalized_cosmos_dndz_z[index_x] = lnx[index_x];
    pclass_sz->normalized_cosmos_dndz_phig[index_x] = lnI[index_x];
    // printf("z=%.3e phig=%.3e\n",pclass_sz->normalized_dndz_z[index_x]);
  };
// exit(0);
  /** Release the memory used locally */
  free(lnx);
  free(lnI);

    if (pclass_sz->sz_verbose > 0){
      printf("-> Cosmos dndz file for unWISE galaxies loaded.\n");
    }

  return _SUCCESS_;
}


int load_M_min_of_z(struct class_sz_structure * pclass_sz)
{

  if (pclass_sz->sz_verbose >= 1)
    printf("-> loading the minimal mass vs redshift\n");


  class_alloc(pclass_sz->M_min_of_z_z,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->M_min_of_z_M_min,sizeof(double *)*100,pclass_sz->error_message);
  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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

  if (pclass_sz->sz_verbose >= 1)
    printf("-> File Name: %s\n",pclass_sz->full_path_to_redshift_dependent_M_min);
  class_open(process,pclass_sz->full_path_to_redshift_dependent_M_min, "r",pclass_sz->error_message);


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
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->M_min_of_z_size = n_data;
  /** Make room */

  class_realloc(pclass_sz->M_min_of_z_z,
                pclass_sz->M_min_of_z_z,
                pclass_sz->M_min_of_z_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->M_min_of_z_M_min,
                pclass_sz->M_min_of_z_M_min,
                pclass_sz->M_min_of_z_size*sizeof(double),
                pclass_sz->error_message);



  /** Store them */
  for (index_x=0; index_x<pclass_sz->M_min_of_z_size; index_x++) {
    pclass_sz->M_min_of_z_z[index_x] = lnx[index_x];
    pclass_sz->M_min_of_z_M_min[index_x] = lnI[index_x];
  };

  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}



int load_ksz_filter(struct class_sz_structure * pclass_sz)
{

  if (pclass_sz->sz_verbose >= 1)
    printf("-> loading the filter f(l) for cl^kSZ2_gal\n");


  class_alloc(pclass_sz->l_unwise_filter,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->f_unwise_filter,sizeof(double *)*100,pclass_sz->error_message);
  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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

  class_open(process,pclass_sz->ksz_filter_file, "r",pclass_sz->error_message);
  if (pclass_sz->sz_verbose >= 1)
    printf("-> File Name: %s\n",pclass_sz->ksz_filter_file);


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
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->unwise_filter_size = n_data;
  /** Make room */

  class_realloc(pclass_sz->l_unwise_filter,
                pclass_sz->l_unwise_filter,
                pclass_sz->unwise_filter_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->f_unwise_filter,
                pclass_sz->f_unwise_filter,
                pclass_sz->unwise_filter_size*sizeof(double),
                pclass_sz->error_message);



  /** Store them */
  for (index_x=0; index_x<pclass_sz->unwise_filter_size; index_x++) {
    pclass_sz->l_unwise_filter[index_x] = lnx[index_x];
    pclass_sz->f_unwise_filter[index_x] = lnI[index_x];
  };

  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  if (pclass_sz->sz_verbose >= 1)
    printf("-> filter f(l) for cl^kSZ2_x successfully loaded.\n");


  return _SUCCESS_;
}




//This routine reads the tabulated
//alpha(z) normalisation for Tinker et al 2010 HMF
//and stores the tabulated values.

int load_T10_alpha_norm(struct class_sz_structure * pclass_sz)
{


  class_alloc(pclass_sz->T10_ln1pz,sizeof(double *)*100,pclass_sz->error_message);
  class_alloc(pclass_sz->T10_lnalpha,sizeof(double *)*100,pclass_sz->error_message);
  //class_alloc(pclass_sz->PP_d2lnI,sizeof(double *)*100,pclass_sz->error_message);

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

  // class_open(process,"class_sz_auxiliary_files/Tinker_et_al_10_alpha_consistency_msyriac.txt", "r",pclass_sz->error_message);
  class_open(process,pclass_sz->Tinker_et_al_10_alpha_consistency_msyriac_file, "r",pclass_sz->error_message);


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
                 pclass_sz->error_message,
                 "Error allocating memory to read the pressure profile.\n");
      lnx = tmp;
      tmp = (double *)realloc(lnI, n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
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
                 pclass_sz->error_message,
                 "The ell/ells's are not strictly sorted in ascending order, "
                 "as it is required for the calculation of the splines.\n");
    }
  }

  /* Close the process */
  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  /** 3. Store the read results into CLASS structures */
  pclass_sz->T10_lnalpha_size = n_data;
  /** Make room */

  class_realloc(pclass_sz->T10_ln1pz,
                pclass_sz->T10_ln1pz,
                pclass_sz->T10_lnalpha_size*sizeof(double),
                pclass_sz->error_message);
  class_realloc(pclass_sz->T10_lnalpha,
                pclass_sz->T10_lnalpha,
                pclass_sz->T10_lnalpha_size*sizeof(double),
                pclass_sz->error_message);



  /** Store them */
  for (index_x=0; index_x<pclass_sz->T10_lnalpha_size; index_x++) {
    pclass_sz->T10_ln1pz[index_x] = log(1.+lnx[index_x]);
    pclass_sz->T10_lnalpha[index_x] = log(lnI[index_x]);
  };

  /** Release the memory used locally */
  free(lnx);
  free(lnI);

  return _SUCCESS_;
}



/**
 * Calculates the mean density contrast at a given redshift based on a specified critical density contrast.
 *
 * @param delta_crit The critical density contrast threshold.
 * @param z The redshift at which the calculation is performed.
 * @param ptsz Pointer to a struct class_sz_structure containing necessary cosmological parameters.
 * @return The mean density contrast adjusted for the matter density at redshift z.
 *
 * This function computes the mean density contrast, delta_mean, using the formula:
 * delta_mean = delta_crit / Omega_m_z
 * where Omega_m_z is the matter density parameter at redshift z, calculated from the cosmological parameters
 * contained in ptsz. This is useful for cosmological studies where the evolution of density contrast relative
 * to the background cosmology at different redshifts is of interest.
 */
double get_delta_mean_from_delta_crit_at_z(double delta_crit,
                                           double z,
                                           struct class_sz_structure *  pclass_sz){
    double om0 = pclass_sz->Omega_m_0;
    double om0_nonu = pclass_sz->Omega0_cdm + pclass_sz->Omega0_b;
    double or0 = pclass_sz->Omega_r_0;
    double ol0 = 1. - om0 - or0;
    double Omega_m_z = om0_nonu * pow(1. + z, 3.) / (om0 * pow(1. + z, 3.) + ol0 + or0 * pow(1. + z, 4.)); // omega_matter without neutrinos
    double delta_mean = delta_crit / Omega_m_z;

    return delta_mean;
}


//HMF Tinker et al 2008
//interpolated at mdeltac
int MF_T08_m500(
                double * result,
                double * lognu ,
                double z ,
                double delta_crit,
                struct class_sz_structure * pclass_sz
                )
{
  //T08@m500
  if (pclass_sz->hmf_apply_zthreshold_to_hmf_and_bias){
  if(z>3.) z=3.; // ccl doesnt have this.. commenting for now.
}
  // double om0 = pclass_sz->Omega_m_0;
  // double or0 = pclass_sz->Omega_r_0;
  // double ol0 = 1.-om0-or0;
  // double Omega_m_z = om0*pow(1.+z,3.)/(om0*pow(1.+z,3.)+ ol0 + or0*pow(1.+z,4.));
  // double  delta_mean = delta_crit/Omega_m_z;
  double  delta_mean = get_delta_mean_from_delta_crit_at_z(delta_crit,z,pclass_sz);
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
              pclass_sz->error_message);
  class_alloc(a0,
              1*sizeof(double),
              pclass_sz->error_message);
  class_alloc(b0,
              1*sizeof(double),
              pclass_sz->error_message);
  class_alloc(c0,
              1*sizeof(double),
              pclass_sz->error_message);

  if (pclass_sz->no_spline_in_tinker == 1){
    // printf("interpolating without splines.\n");

    *Ap0 = pwl_value_1d(9,
                       delta_mean_tab,
                       A_tab,
                       delta_mean);

    *a0 = pwl_value_1d(9,
                       delta_mean_tab,
                       aa_tab,
                       delta_mean);
    *b0 = pwl_value_1d(9,
                       delta_mean_tab,
                       b_tab,
                       delta_mean);
    *c0 = pwl_value_1d(9,
                       delta_mean_tab,
                       c_tab,
                       delta_mean);


  }
  else{

    // printf("interpolating with splines.\n");
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

  // printf("interpolation done %.5e.\n",c0);
  }

  double alphaT08 =
  pow(10.,-pow(0.75/log10(pow(10.,delta_mean)/75.),1.2)); // pow(10.,delta_mean) is delta
  // note in ccl this is written slighlty diff:
  // def _pd(self, ld):
  //     return 10.**(-(0.75/(ld - 1.8750612633))**1.2)
  // where 1.8750612633=np.log10(75)

  double   Ap=*Ap0*pow(1.+z,-0.14);
  double   a=*a0*pow(1.+z,-0.06);
  double   b=*b0*pow(1.+z,-alphaT08);
  double   c=*c0;
  double   nu= exp(*lognu);
  double sigma= pclass_sz->delta_cSZ/sqrt(nu);

  free(Ap0);
  free(a0);
  free(b0);
  free(c0);

  // sigma = .15;

  *result = 0.5*(Ap*(pow(sigma/b,-a)+1.)*exp(-c/pow(sigma,2.)));
  // *result = 0.5*(Ap);
  // *result = 0.5*(Ap*(pow(1.0/b,-a)+1.)*exp(-c/pow(1.0,2.)));
  // *result = 0.5;
  // *result = 0.5*exp(-1./pow(sigma,2));
  // *result = 0.5*sigma;
  return _SUCCESS_;
}



int p_gnfw(double * p_gnfw_x,
             double x ,
             double kl,
             double * pvectsz,
             struct background * pba,
             struct class_sz_structure * pclass_sz)
{
  //Custom. GNFW pressure profile
  //int index_l = (int) pvectsz[pclass_sz->index_multipole_for_pressure_profile];
  int index_md = (int) pvectsz[pclass_sz->index_md]; // important for mean y and dydz computation


  //}



  //Battaglia et al 2012 pressure profile
  //Eq. 10
  if(pclass_sz->pressure_profile == 4){
        //pclass_sz->P0_B12 = 18.1;
        //pclass_sz->xc_B12 = 0.497;
        //pclass_sz->beta_B12 = 4.35;

        //pclass_sz->alpha_m_P0_B12 = 0.154;
        //pclass_sz->alpha_m_xc_B12 = -0.00865;
        //pclass_sz->alpha_m_beta_B12 = 0.0393;

        //pclass_sz->alpha_z_P0_B12 = -0.758;
        //pclass_sz->alpha_z_xc_B12 = 0.731;
        //pclass_sz->alpha_z_beta_B12 = 0.415;



        // double xc;
        // double beta;
        // double P0;
        //
        double m200_over_msol = pvectsz[pclass_sz->index_m200c]/pba->h; // convert to Msun
        // double z = pvectsz[pclass_sz->index_z];
        //
        //
        // P0 = pclass_sz->P0_B12*pow(m200_over_msol/1e14,pclass_sz->alpha_m_P0_B12)*pow(1+z,pclass_sz->alpha_z_P0_B12);
        // xc = pclass_sz->xc_B12*pow(m200_over_msol/1e14,pclass_sz->alpha_m_xc_B12)*pow(1+z,pclass_sz->alpha_z_xc_B12);
        // beta = pclass_sz->beta_B12*pow(m200_over_msol/1e14,pclass_sz->alpha_m_beta_B12)*pow(1+z,pclass_sz->alpha_z_beta_B12);
        //
        // double gamma = pclass_sz->gamma_B12;
        // double alpha = pclass_sz->alpha_B12;

double c_asked = pclass_sz->c_B12;
double Px = get_pressure_P_over_P_delta_at_x_M_z_b12_200c(x,m200_over_msol,pvectsz[pclass_sz->index_z],
                                              c_asked,pclass_sz->P0_B12,
                                              pclass_sz->xc_B12,pclass_sz->beta_B12,
                                              pclass_sz->alpha_m_P0_B12,pclass_sz->alpha_m_xc_B12,
                                              pclass_sz->alpha_m_beta_B12,pclass_sz->alpha_z_P0_B12,
                                              pclass_sz->alpha_z_xc_B12,pclass_sz->alpha_z_beta_B12,
                                              // break model
                                  						pclass_sz->mcut_B12,pclass_sz->alphap_m_P0_B12,
                                  						pclass_sz->alphap_m_xc_B12,pclass_sz->alphap_m_beta_B12,
                                  						pclass_sz->alpha_c_P0_B12,
                                  						pclass_sz->alpha_c_xc_B12,
                                  						pclass_sz->alpha_c_beta_B12,
                                                     // end break model
                                              pclass_sz->alpha_B12,
                                              pclass_sz->gamma_B12,
                                              pba,pclass_sz);

          // *p_gnfw_x = P0*pow(x/xc,gamma)*pow(1.+ pow(x/xc,alpha),-beta)
          //               *pow(x,2)
          //               /(x*kl);

        *p_gnfw_x = Px*pow(x,2)/(x*kl);

    double z = pvectsz[pclass_sz->index_z];
    double P0 = pclass_sz->P0_B12*pow(m200_over_msol/1e14,pclass_sz->alpha_m_P0_B12)*pow(1+z,pclass_sz->alpha_z_P0_B12);
    double xc = pclass_sz->xc_B12*pow(m200_over_msol/1e14,pclass_sz->alpha_m_xc_B12)*pow(1+z,pclass_sz->alpha_z_xc_B12);
    double beta = pclass_sz->beta_B12*pow(m200_over_msol/1e14,pclass_sz->alpha_m_beta_B12)*pow(1+z,pclass_sz->alpha_z_beta_B12);

    double gamma = pclass_sz->gamma_B12;
    double alpha = pclass_sz->alpha_B12;

      *p_gnfw_x = P0*pow(x/xc,gamma)*pow(1.+ pow(x/xc,alpha),-beta)
                    *pow(x,2)
                    // /(x*(kl+0.5)/pvectsz[pclass_sz->index_l200c]);
                    /(x*kl);




  if (_mean_y_ || _dydz_)
        // *p_gnfw_x = P0*pow(x/xc,gamma)*pow(1.+ pow(x/xc,alpha),-beta)*pow(x,2);
        *p_gnfw_x = Px*pow(x,2);
    }
  else{
      *p_gnfw_x = 0.;

      // Example Arnaud 2010
      // pclass_sz->P0GNFW = 8.130;
      // pclass_sz->c500 = 1.156;
      // pclass_sz->gammaGNFW = 0.3292;
      // pclass_sz->alphaGNFW = 1.0620;
      // pclass_sz->betaGNFW = 5.4807;

      //Custom. GNFW
      //if(pclass_sz->pressure_profile == 3){
        *p_gnfw_x = (1./(pow(pclass_sz->c500*x,pclass_sz->gammaGNFW)
                      *pow(1.+ pow(pclass_sz->c500*x,pclass_sz->alphaGNFW),
                          (pclass_sz->betaGNFW-pclass_sz->gammaGNFW)/pclass_sz->alphaGNFW)))
                      *pow(x,2)
                      //*sin(x*(pclass_sz->ell[index_l]+0.5)/pvectsz[pclass_sz->index_l500])
                      /(x*kl);


      if (_mean_y_ || _dydz_)
        *p_gnfw_x = (1./(pow(pclass_sz->c500*x,pclass_sz->gammaGNFW)
                         *pow(1.+ pow(pclass_sz->c500*x,pclass_sz->alphaGNFW),
                             (pclass_sz->betaGNFW-pclass_sz->gammaGNFW)/pclass_sz->alphaGNFW)))
                        *pow(x,2);
      }


  return _SUCCESS_;
}



// HMF Tinker 2010
// https://wwwmpa.mpa-garching.mpg.de/~komatsu/CRL/clusters/szpowerspectrumdks/szfastdks/mf_tinker10.f90
int MF_T10 (
            double * result,
            double * lognu ,
            double z ,
            struct class_sz_structure * pclass_sz
            )
{
  *result = get_f_tinker10_at_nu_and_z(exp(*lognu),z,pclass_sz);
  return _SUCCESS_;
}


double get_f_tinker10_at_nu_and_z(double nu, double z,struct class_sz_structure * pclass_sz){
if (pclass_sz->hmf_apply_zthreshold_to_hmf_and_bias){
  if(z>3.) z=3.;
}
  double alpha;



  // fix alpha or not
  if (pclass_sz->T10_alpha_fixed == 1){
  alpha = pclass_sz->alphaSZ;
  }
  else{
  alpha = get_T10_alpha_at_z(z,pclass_sz);
  }

  double lognu = log(nu);


  double result = 0.5
            *alpha
            *(1.+pow(pow(pclass_sz->beta0SZ*pow(1.+z,0.2),2.)
                   *exp(lognu),
                   -pclass_sz->phi0SZ
                   *pow(1.+z,-0.08)))*pow(
                    exp(lognu),
                    pclass_sz->eta0SZ
                    *pow(1.+z,0.27))
            *exp(-pclass_sz->gamma0SZ
                 *pow(1.+z,-0.01)
                 *exp(lognu)/2.)
            *sqrt(exp(lognu));
  return result;


}



// //HMF Tinker et al 2008
// //@ M200m
// ---> deprecated: all done in the one that we interpolate
// int MF_T08(
//            double * result,
//            double * lognu ,
//            double z ,
//            struct class_sz_structure * pclass_sz
//            )
// {
//   // double alphaT08 = pow(10.,-pow(0.75/log10(200./75.),1.2));
//   //
//   // double   Ap=0.186*pow(1.+z,-0.14);
//   // double   a=1.47*pow(1.+z,-0.06);
//   // double   b=2.57*pow(1.+z,-alphaT08);
//   // double   c=1.19;
//   // double   nu= exp(*lognu);
//   // double sigma= pclass_sz->delta_cSZ/sqrt(nu);
//   //
//   // *result = 0.5*(Ap*(pow(sigma/b,-a)+1.)*exp(-c/pow(sigma,2.)));
//
//   *result = get_f_tinker08_at_nu_and_z(exp(*lognu),z,pclass_sz);
//
//   return _SUCCESS_;
// }



double get_f_tinker08_at_nu_and_z(double nu, double z,  struct class_sz_structure * pclass_sz){

if (pclass_sz->hmf_apply_zthreshold_to_hmf_and_bias){
  if(z>2.5) z=2.5; // see sec 4 of https://arxiv.org/pdf/0803.2706.pdf
  // this is not in CCL
}

  double alphaT08 = pow(10.,-pow(0.75/log10(200./75.),1.2));

  // double   Ap=0.186*pow(1.+z,-0.14);
  // double   a=1.47*pow(1.+z,-0.06);
  // double   b=2.57*pow(1.+z,-alphaT08);
  // double   c=1.19;
  //
  // A_hmfcalc = 1.858659e-01
  // a_hmfcalc = 1.466904
  // b_hmfcalc = 2.571104
  // c_hmfcalc = 1.193958

  double   Ap=1.858659e-01*pow(1.+z,-0.14);
  double   a=1.466904*pow(1.+z,-0.06);
  double   b=2.571104*pow(1.+z,-alphaT08);
  double   c=1.193958;
  // double   nu= exp(*lognu);
  double   sigma= pclass_sz->delta_cSZ/sqrt(nu);


  return 0.5*(Ap*(pow(sigma/b,-a)+1.)*exp(-c/pow(sigma,2.)));
}

//HMF Tinker et al 2008
//@ M1600m
int MF_T08_M1600m(
                  double * result,
                  double * lognu ,
                  double z ,
                  struct class_sz_structure * pclass_sz
                  )
{
  double alphaT08 = pow(10.,-pow(0.75/log10(1600./75.),1.2));
  double   Ap=0.260*pow(1.+z,-0.14);
  double   a=2.30*pow(1.+z,-0.06);
  double   b=1.46*pow(1.+z,-alphaT08);
  double   c=1.97;
  double   nu= exp(*lognu);
  double sigma= pclass_sz->delta_cSZ/sqrt(nu);

  *result = 0.5*(Ap*(pow(sigma/b,-a)+1.)*exp(-c/pow(sigma,2.)));
  return _SUCCESS_;
}


//HMF Boquet et al 2015
int MF_B15(
           double * result,
           double * lognu ,
           double z ,
           struct class_sz_structure * pclass_sz
           )
{
  //B15
  double   Ap=pclass_sz->Ap0*pow(1.+z,0.285);
  double   a=pclass_sz->a0*pow(1.+z,-0.058);
  double   b=pclass_sz->b0*pow(1.+z,-0.366);
  double   c=pclass_sz->c0*pow(1.+z,-0.045);
  double   nu= exp(*lognu);
  double sigma= pclass_sz->delta_cSZ/sqrt(nu);

  *result = 0.5*(Ap*(pow(sigma/b,-a)+1.)*exp(-c/pow(sigma,2.)));

  return _SUCCESS_;
}


//HMF Boquet et al 2015 @ M500c
//TBC (11.04.19)
int MF_B15_M500c(double * result,
                 double * lognu ,
                 double z ,
                 struct class_sz_structure * pclass_sz)
{
  //B15

  double   Ap=0.180*pow(1.+z,1.088);
  double   a=2.29*pow(1.+z,0.150);
  double   b=2.44*pow(1.+z,-1.008);
  double   c=1.97*pow(1.+z,-0.322);
  double   nu= exp(*lognu);
  double sigma= pclass_sz->delta_cSZ/sqrt(nu);

  *result = 0.5*(Ap*(pow(sigma/b,-a)+1.)*exp(-c/pow(sigma,2.)));

  return _SUCCESS_;
}



int MF_J01(double * result,
           double * lognu ,
           struct class_sz_structure * pclass_sz)
{
  double   nu= exp(*lognu);
  double sigma= pclass_sz->delta_cSZ/sqrt(nu);

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

double next_z(double z_i, double dz, struct class_sz_structure * pclass_sz){
  // Compute redshift bins where bins are defined with higher resolution for z<0.2
  double dz_i;
  double highres_z_cutoff = pclass_sz->cluster_count_completeness_grid_z_cutoff_low;
  if (z_i < highres_z_cutoff)
    dz_i = pclass_sz->dz_cluster_count_completeness_grid_low_z;
  //else if ((z_i >= highres_z_cutoff) && (z_i <= 1.))
  else if ((z_i >= highres_z_cutoff) && (z_i <= pclass_sz->cluster_count_completeness_grid_z_cutoff_mid))
    dz_i = pclass_sz->dz_cluster_count_completeness_grid_mid_z;
  else
    dz_i = pclass_sz->dz_cluster_count_completeness_grid_high_z;;

  double next_z = z_i + dz_i;
  return next_z;
}

double erf_compl_nicola(double y,
                        double sn,
                        double q,
                        double ymin,
                        double ymax,
                        double szcc_dof){
  //Completeness with error function
  // double arg = (y - q * sn)/(sqrt(2.) * sn);
  double arg1;
  double ylim;
  ylim = ymax;
  arg1 = (sqrt(y/sn*y/sn+szcc_dof)-ylim)/(sqrt(2.));
  double arg2;
  if (ymin>q) ylim = ymin;
  else ylim = q;
  arg2 = (sqrt(y/sn*y/sn+szcc_dof)-ylim)/(sqrt(2.));
  double erf_compl = (erf(arg2) - erf(arg1))/2.;

  if (ymax<q) erf_compl = 1e-100;
  return erf_compl;
}



double erf_compl(double y,
                 double sn,
                 double q,
                 double dof){
  //Completeness with error function
  // double arg = (y - q * sn)/(sqrt(2.) * sn);
  // with optimization bias
  double arg = (sqrt(pow(y/sn,2.)+pow(dof,1.)) - q )/(sqrt(2.));
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
  struct perturbs * ppt;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvecback;
  double * pvectsz;
  //double * llprime_grid;
};

double integrand_redshift(double ln1pz, void *p){

  struct Parameters_for_integrand_redshift *V = ((struct Parameters_for_integrand_redshift *) p);

   V->pvectsz[V->pclass_sz->index_z] = exp(ln1pz)-1.;

  double z =  V->pvectsz[V->pclass_sz->index_z];




  //Evaluation of background quantities @ z:
  double tau;
  int first_index_back = 0;

  class_call(background_tau_of_z(V->pba,z,&tau),
             V->pclass_sz->error_message,
             V->pclass_sz->error_message);

  class_call(background_at_tau(V->pba,
                               tau,
                               V->pba->long_info,
                               V->pba->inter_normal,
                               &first_index_back,
                               V->pvecback),
             V->pclass_sz->error_message,
             V->pclass_sz->error_message);


  //volume element in units h^-3 Mpc^3
  //volume = dv/(dzdOmega)*(c/H)
  // Chi^2 dChi = dV/(dzdOmega)*(c/H) dz

  V->pvectsz[V->pclass_sz->index_volume] = pow(1.+z,2)
                                      *pow(V->pvecback[V->pba->index_bg_ang_distance]*V->pba->h,2)
                                      *_c_*1.e-5
                                      /(V->pvecback[V->pba->index_bg_H]/V->pba->H0);


  V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2); // conformal distance squared in [Mpc/h]^2
  double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
  V->pvectsz[V->pclass_sz->index_dgdz] = V->pvecback[V->pba->index_bg_D]*(1.-V->pvecback[V->pba->index_bg_f]); // d/dz(D/a)



  int index_md = (int) V->pvectsz[V->pclass_sz->index_md];


  double kl;


  // critical density at z in (Msun/h)/(Mpc/h)^3
  V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                        *pow(_Mpc_over_m_,1)
                                        *pow(_c_,2)
                                        *V->pvecback[V->pba->index_bg_rho_crit]
                                        /pow(V->pba->h,2);

  double Eh = V->pvecback[V->pba->index_bg_H]/V->pba->H0;
  double omega = V->pvecback[V->pba->index_bg_Omega_m];//pow(Eh,2.);
  V->pvectsz[V->pclass_sz->index_Delta_c] = Delta_c_of_Omega_m(omega);





if     (((V->pclass_sz->has_tSZ_gal_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_gal_1h))
     || ((V->pclass_sz->has_tSZ_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_gal_2h))
     || ((V->pclass_sz->has_IA_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_IA_gal_2h))
     || ((V->pclass_sz->has_gal_gal_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_1h))
     || ((V->pclass_sz->has_gal_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_2h))
     || ((V->pclass_sz->has_pk_gg_at_z_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_gg_at_z_1h))
     || ((V->pclass_sz->has_pk_gg_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_gg_at_z_2h))
     || ((V->pclass_sz->has_gal_lens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lens_1h))
     || ((V->pclass_sz->has_gal_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lens_2h))
     || ((V->pclass_sz->has_tau_gal_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_tau_gal_1h))
     || ((V->pclass_sz->has_tau_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tau_gal_2h))
     || ((V->pclass_sz->has_gal_cib_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_cib_1h))
     || ((V->pclass_sz->has_gal_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_cib_2h))
     || ((V->pclass_sz->has_gal_lensmag_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lensmag_1h))
     || ((V->pclass_sz->has_gal_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lensmag_2h))
     || ((V->pclass_sz->has_gal_gallens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gallens_1h))
     || ((V->pclass_sz->has_gal_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gallens_2h))
     // || ((V->pclass_sz->has_lensmag_lensmag_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_lensmag_lensmag_1h))
     // || ((V->pclass_sz->has_lensmag_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lensmag_lensmag_2h))
     // || ((V->pclass_sz->has_lens_lensmag_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lensmag_1h))
     // || ((V->pclass_sz->has_lens_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lensmag_2h))
     //|| ((V->pclass_sz->has_kSZ_kSZ_lensmag_1halo == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lensmag_1halo))
     || ((V->pclass_sz->has_kSZ_kSZ_gal_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_1h))
     || ((V->pclass_sz->has_kSZ_kSZ_gal_1h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_1h_fft))
     || ((V->pclass_sz->has_kSZ_kSZ_gal_2h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_2h_fft))
     || ((V->pclass_sz->has_kSZ_kSZ_gal_3h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_3h_fft))
     || ((V->pclass_sz->has_gal_gal_lens_1h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_lens_1h_fft))
     || ((V->pclass_sz->has_gal_gal_lens_2h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_lens_2h_fft))
     || ((V->pclass_sz->has_gal_gal_lens_3h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_lens_3h_fft))
     || ((V->pclass_sz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_2h))
     || ((V->pclass_sz->has_kSZ_kSZ_gal_3h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_3h))
     // || ((V->pclass_sz->has_kSZ_kSZ_gal_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_hf))

    ) {

 V->pvectsz[V->pclass_sz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->pclass_sz);
 // printf("z = %.5e ng = %.5e\n",z,V->pvectsz[V->pclass_sz->index_mean_galaxy_number_density]);

}


  double result = 0.;

  // // first deal with quantities that does not require mass integration:
  // if ((V->pclass_sz->has_pk_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_at_z_2h)) {
  //
  //   // eq. 10 of https://arxiv.org/pdf/1505.07833.pdf
  //   result = 1.;
  // }
  // else
  if ((V->pclass_sz->has_gal_gal_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_hf)) {
    int index_l = (int) V->pvectsz[V->pclass_sz->index_multipole];
    double l = V->pclass_sz->ell[index_l];
    double pk1;
    if (V->pclass_sz->use_pkl_in_linbias_calc){
        pk1 =  get_pk_lin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
       }
    else{
        pk1 =  get_pk_nonlin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
        }
    result = pk1;
    evaluate_effective_galaxy_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);
    result *= V->pvectsz[V->pclass_sz->index_halo_bias]*V->pvectsz[V->pclass_sz->index_halo_bias];


    // printf("b=%.5e pk=%.5e\n",V->pvectsz[V->pclass_sz->index_halo_bias],pk1);

  }
  else if ((V->pclass_sz->has_ngal_ngal_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_ngal_hf)) {

    
    int index_l = (int) V->pvectsz[V->pclass_sz->index_multipole];
    double l = V->pclass_sz->ell[index_l];
    double pk1;
    if (V->pclass_sz->use_pkl_in_linbias_calc){
        pk1 =  get_pk_lin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
       }
    else{
        pk1 =  get_pk_nonlin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
        }
    result = pk1;
    int index_g = (int) V->pvectsz[V->pclass_sz->index_ngal_for_galaxy_profile];
    int index_g_prime = (int) V->pvectsz[V->pclass_sz->index_ngal_prime_for_galaxy_profile];
    evaluate_effective_galaxy_bias_ngal(index_g,V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);
    double bg = V->pvectsz[V->pclass_sz->index_halo_bias];
    evaluate_effective_galaxy_bias_ngal(index_g_prime,V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);
    double bg_prime =  V->pvectsz[V->pclass_sz->index_halo_bias];
    result *= bg*bg_prime;

    if (V->pclass_sz->use_nl_bias){

      double pknl = get_pk_nonlin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
      double pkl = get_pk_lin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
      // double bnl = V->pclass_sz->bnl;

      double bnl =  V->pclass_sz->effective_galaxy_bias_nl_ngal[index_g];
      double bnl_prime =  V->pclass_sz->effective_galaxy_bias_nl_ngal[index_g_prime];
      result = bg*bg_prime*pkl + bnl*bnl_prime*(pknl-pkl);

    }


    // printf("result = %.5e b=%.5e pk=%.5e\n",result,V->pvectsz[V->pclass_sz->index_halo_bias],pk1);

  }
  else if ((V->pclass_sz->has_kSZ_kSZ_gal_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_hf)) {

  int index_theta_1 = (int) V->pvectsz[V->pclass_sz->index_multipole_1];
  int index_l_2 = (int) V->pvectsz[V->pclass_sz->index_multipole_2];
  int index_l_3 = (int) V->pvectsz[V->pclass_sz->index_multipole_3];


  double l2 = exp(V->pclass_sz->ell_kSZ2_gal_multipole_grid[index_l_2]);
  double l3 = V->pclass_sz->ell[index_l_3];
  double ell = l3;
  double ell_prime = l2;
  double theta_1 = V->pclass_sz->theta_kSZ2_gal_theta_grid[index_theta_1];
  double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos(theta_1));


  double z = V->pvectsz[V->pclass_sz->index_z];
  double d_A = V->pvecback[V->pba->index_bg_ang_distance]*V->pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi

  double k1 = (l1+0.5)/d_A;
  double k2 = (l2+0.5)/d_A;
  double k3 = (l3+0.5)/d_A;

  result = get_ttg_bispectrum_at_z_effective_approach(k1,k2,k3,z,V->pclass_sz,V->pba,V->pnl,V->ppm);

  double bg = 1.;
  if (V->pclass_sz->use_bg_at_z_in_ksz2g_eff==1){
    bg = get_mean_galaxy_bias_at_z(z,V->pclass_sz);
  }
  else if (V->pclass_sz->use_bg_eff_in_ksz2g_eff==1){
    bg = V->pclass_sz->effective_galaxy_bias;
  }
  result *= bg;


  }

  else if(
    ((V->pclass_sz->has_kSZ_kSZ_gallens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gallens_hf))
  ||((V->pclass_sz->has_kSZ_kSZ_lens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lens_hf))
  ) {

  int index_theta_1 = (int) V->pvectsz[V->pclass_sz->index_multipole_1];
  int index_l_2 = (int) V->pvectsz[V->pclass_sz->index_multipole_2];
  int index_l_3 = (int) V->pvectsz[V->pclass_sz->index_multipole_3];


  double l2 = exp(V->pclass_sz->ell_kSZ2_gal_multipole_grid[index_l_2]);
  double l3 = V->pclass_sz->ell[index_l_3];
  double ell = l3;
  double ell_prime = l2;
  double theta_1 = V->pclass_sz->theta_kSZ2_gal_theta_grid[index_theta_1];
  double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos(theta_1));


  double z = V->pvectsz[V->pclass_sz->index_z];
  double d_A = V->pvecback[V->pba->index_bg_ang_distance]*V->pba->h*(1.+z); //multiply by h to get in Mpc/h => conformal distance Chi

  double k1 = (l1+0.5)/d_A;
  double k2 = (l2+0.5)/d_A;
  double k3 = (l3+0.5)/d_A;

  result = get_ttg_bispectrum_at_z_effective_approach(k1,k2,k3,z,V->pclass_sz,V->pba,V->pnl,V->ppm);

  }

  else if ((V->pclass_sz->has_isw_lens == _TRUE_) && (index_md == V->pclass_sz->index_md_isw_lens)) {

  double delta_ell_lens =  delta_ell_lens_at_ell_and_z(V->pvecback,
                                                  V->pvectsz,
                                                  V->pba,
                                                  V->ppm,
                                                  V->pnl,
                                                  V->pclass_sz);

  double delta_ell_isw = delta_ell_isw_at_ell_and_z(V->pvecback,
                                                          V->pvectsz,
                                                          V->pba,
                                                          V->ppm,
                                                          V->pnl,
                                                          V->pclass_sz);
  result = delta_ell_lens*delta_ell_isw;

  }


  else if ((V->pclass_sz->has_isw_tsz == _TRUE_) && (index_md == V->pclass_sz->index_md_isw_tsz)){

  double delta_ell_isw = delta_ell_isw_at_ell_and_z(V->pvecback,
                                                          V->pvectsz,
                                                          V->pba,
                                                          V->ppm,
                                                          V->pnl,
                                                          V->pclass_sz);
  double delta_ell_y = integrate_over_m_at_z(V->pvecback,
                                              V->pvectsz,
                                              V->pba,
                                              V->pnl,
                                              V->ppm,
                                              V->ppt,
                                              V->pclass_sz);

  result = delta_ell_isw*delta_ell_y;


  }
  else if ((V->pclass_sz->has_isw_auto == _TRUE_) && (index_md == V->pclass_sz->index_md_isw_auto)){

  double delta_ell_isw = delta_ell_isw_at_ell_and_z(V->pvecback,
                                                    V->pvectsz,
                                                    V->pba,
                                                    V->ppm,
                                                    V->pnl,
                                                    V->pclass_sz);

  result = delta_ell_isw*delta_ell_isw;

  }

  // Halofit approach
  else if ((V->pclass_sz->has_gal_lens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lens_hf))
  {

  int index_l = (int) V->pvectsz[V->pclass_sz->index_multipole];
  double l = V->pclass_sz->ell[index_l];

    double pk1;
    if (V->pclass_sz->use_pkl_in_linbias_calc){
        pk1 =  get_pk_lin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
       }
    else{
        pk1 =  get_pk_nonlin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
        }

  evaluate_effective_galaxy_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);

  result *= V->pvectsz[V->pclass_sz->index_halo_bias];

  double W_lens =  radial_kernel_W_lensing_at_z(V->pvecback,
                                                V->pvectsz,
                                                V->pba,
                                                V->ppm,
                                                V->pnl,
                                                V->pclass_sz);
  // this is needed only in  the approximate calculation...
  // for the exact calculation in HOD, this comes out of Sigma_crit
  result *= W_lens;
    // printf("result gal lens hf = %.3e  b = %.3e z = %.3e\n",result,V->pvectsz[V->pclass_sz->index_halo_bias],z);
    // exit(0);
  }

  else if ((V->pclass_sz->has_ngal_lens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_lens_hf)) {
    int index_l = (int) V->pvectsz[V->pclass_sz->index_multipole];
    double l = V->pclass_sz->ell[index_l];


    double pk1;
    if (V->pclass_sz->use_pkl_in_linbias_calc){
        pk1 =  get_pk_lin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
       }
    else{
        pk1 =  get_pk_nonlin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
        }

    result = pk1;
    int index_g = (int) V->pvectsz[V->pclass_sz->index_ngal_for_galaxy_profile];
    evaluate_effective_galaxy_bias_ngal(index_g,V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);
    double bg = V->pvectsz[V->pclass_sz->index_halo_bias];
    result *= bg;


    if (V->pclass_sz->use_nl_bias){

      double pknl = get_pk_nonlin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
      double pkl = get_pk_lin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
      // double bnl = V->pclass_sz->bnl;

      double bnl = V->pclass_sz->effective_galaxy_bias_nl_ngal[index_g];
      result = bg*pkl + bnl*(pknl-pkl);

    }


    double W_lens =  radial_kernel_W_lensing_at_z(V->pvecback,
                                                  V->pvectsz,
                                                  V->pba,
                                                  V->ppm,
                                                  V->pnl,
                                                  V->pclass_sz);
    result *= W_lens;
    // printf("result ngal lens hf = %.3e  b = %.3e z = %.3e\n",result,bg,z);
    // exit(0);
  }


  // Halofit approach
  else if ((V->pclass_sz->has_gal_lensmag_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lensmag_hf))
  {

//printf("ok\n");
  int index_l = (int) V->pvectsz[V->pclass_sz->index_multipole];
  double l = V->pclass_sz->ell[index_l];

    double pk1;
    if (V->pclass_sz->use_pkl_in_linbias_calc){
        pk1 =  get_pk_lin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
       }
    else{
        pk1 =  get_pk_nonlin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
        }


  result = pk1;
  evaluate_effective_galaxy_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);

  result *= V->pvectsz[V->pclass_sz->index_halo_bias];

  double W_lensmag =  radial_kernel_W_lensing_magnification_at_z(V->pvecback,
                                                                V->pvectsz,
                                                                V->pba,
                                                                V->ppm,
                                                                V->pnl,
                                                                V->pclass_sz);

    // this is needed only in  the approximate calculation
    // for the exact calculation in HOD, this comes out of Sigma_crit
    result *= W_lensmag;

  }
  // Halofit approach
else if ((V->pclass_sz->has_lensmag_lensmag_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_lensmag_lensmag_hf))
  {

  int index_l = (int) V->pvectsz[V->pclass_sz->index_multipole];
  double l = V->pclass_sz->ell[index_l];
    double pk1;
    if (V->pclass_sz->use_pkl_in_linbias_calc){
        pk1 =  get_pk_lin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
       }
    else{
        pk1 =  get_pk_nonlin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
        }
  result = pk1;

  double W_lensmag =  radial_kernel_W_lensing_magnification_at_z(V->pvecback,
                                                                V->pvectsz,
                                                                V->pba,
                                                                V->ppm,
                                                                V->pnl,
                                                                V->pclass_sz);

    // this is needed only in  the approximate calculation
    // for the exact calculation in HOD, this comes out of Sigma_crit
    result *= W_lensmag*W_lensmag;

  }
  // Halofit approach
else if ((V->pclass_sz->has_lens_lensmag_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lensmag_hf))
{

  int index_l = (int) V->pvectsz[V->pclass_sz->index_multipole];
  double l = V->pclass_sz->ell[index_l];

    double pk1;
    if (V->pclass_sz->use_pkl_in_linbias_calc){
        pk1 =  get_pk_lin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
       }
    else{
        pk1 =  get_pk_nonlin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
        }

  result = pk1;

  double W_lensmag =  radial_kernel_W_lensing_magnification_at_z(V->pvecback,
                                                                V->pvectsz,
                                                                V->pba,
                                                                V->ppm,
                                                                V->pnl,
                                                                V->pclass_sz);
  double W_lens =  radial_kernel_W_lensing_at_z(V->pvecback,
                                                  V->pvectsz,
                                                  V->pba,
                                                  V->ppm,
                                                  V->pnl,
                                                  V->pclass_sz);

    // this is needed only in  the approximate calculation
    // for the exact calculation in HOD, this comes out of Sigma_crit
    //printf("%.3e \t %.3e\n",W_lensmag,W_lens);
    result *= W_lensmag*W_lens;

  }

  // Halofit approach
else if ((V->pclass_sz->has_lens_lens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lens_hf))
  {
  int index_l = (int) V->pvectsz[V->pclass_sz->index_multipole];
  double l = V->pclass_sz->ell[index_l];

    double pk1;
    if (V->pclass_sz->use_pkl_in_linbias_calc){
        pk1 =  get_pk_lin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
       }
    else{
        pk1 =  get_pk_nonlin_at_k_and_z((l+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
        }

  result = pk1;
  }

// else if (
//   ((V->pclass_sz->use_hod == 0) && (V->pclass_sz->has_lens_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lens_2h))
//   || ((V->pclass_sz->use_hod == 0) && (V->pclass_sz->has_lens_lens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lens_1h))
// ){
//
//
// if (index_md == V->pclass_sz->index_md_lens_lens_1h) {
//   result = 0.;
// }
// else {
//   double W_lens =  radial_kernel_W_lensing_at_z(V->pvecback,
//                                                   V->pvectsz,
//                                                   V->pba,
//                                                   V->ppm,
//                                                   V->pnl,
//                                                   V->pclass_sz);
// // this is needed only in  the approximate calculation
// // for the exact calculation in halo model, this comes out of Sigma_crit
// result = W_lens*W_lens;
//
// }
//
// }

  // then quantities that require mass integration
  else {

    if (V->pclass_sz->sz_verbose>10)
      printf("integrating over mass at z = %.3e\n",z);

  result = integrate_over_m_at_z(V->pvecback,
                                 V->pvectsz,
                                 V->pba,
                                 V->pnl,
                                 V->ppm,
                                 V->ppt,
                                 V->pclass_sz);

  if (V->pclass_sz->sz_verbose>10)
    printf("Result of mass integration at z = %.3e: %.3e\n", z, result);

  if ((V->pclass_sz->has_hmf == _TRUE_) && (index_md == V->pclass_sz->index_md_hmf)){
      // printf("returning integrated obver mass, intm = %.3e\n",result);
      result *= (1.+V->pvectsz[V->pclass_sz->index_z])*get_volume_at_z(V->pvectsz[V->pclass_sz->index_z],V->pba);
      return result;
      }


// if computing 3d matter power spectrum P(k) of bispectrum:
// this are not integrated over volume
if( ((V->pclass_sz->has_pk_at_z_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_at_z_1h))
    || ((V->pclass_sz->has_pk_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_at_z_2h))
    || ((V->pclass_sz->has_pk_gg_at_z_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_gg_at_z_1h))
    || ((V->pclass_sz->has_pk_gg_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_gg_at_z_2h))
    || ((V->pclass_sz->has_pk_bb_at_z_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_bb_at_z_1h))
    || ((V->pclass_sz->has_pk_bb_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_bb_at_z_2h))
    || ((V->pclass_sz->has_pk_b_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_b_at_z_2h))
    || ((V->pclass_sz->has_pk_em_at_z_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_em_at_z_1h))
    || ((V->pclass_sz->has_pk_em_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_em_at_z_2h))
    || ((V->pclass_sz->has_pk_HI_at_z_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_HI_at_z_1h))
    || ((V->pclass_sz->has_pk_HI_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_HI_at_z_2h))
    || ((V->pclass_sz->has_bk_at_z_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_bk_at_z_1h))
    || ((V->pclass_sz->has_bk_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_bk_at_z_2h))
    || ((V->pclass_sz->has_bk_at_z_3h == _TRUE_) && (index_md == V->pclass_sz->index_md_bk_at_z_3h))
    || ((V->pclass_sz->has_bk_ttg_at_z_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_bk_ttg_at_z_1h))
    || ((V->pclass_sz->has_bk_ttg_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_bk_ttg_at_z_2h))
    || ((V->pclass_sz->has_bk_ttg_at_z_3h == _TRUE_) && (index_md == V->pclass_sz->index_md_bk_ttg_at_z_3h))
    ){
      int index_k = (int) V->pvectsz[V->pclass_sz->index_k_for_pk_hm];
      kl = V->pclass_sz->k_for_pk_hm[index_k];

        if (((V->pclass_sz->has_pk_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_at_z_2h))
          || ((V->pclass_sz->has_pk_bb_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_bb_at_z_2h))
          || ((V->pclass_sz->has_pk_b_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_b_at_z_2h))
          || ((V->pclass_sz->has_pk_em_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_em_at_z_2h))
          || ((V->pclass_sz->has_pk_gg_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_gg_at_z_2h))
          || ((V->pclass_sz->has_pk_HI_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_HI_at_z_2h))
          ){
            result *= get_pk_lin_at_k_and_z(kl,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
            // evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);
            // result *= V->pvectsz[V->pclass_sz->index_pk_for_halo_bias];
          }



        // if ((V->pclass_sz->has_bk_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_bk_at_z_2h)){
        //   evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);
        //   result *= 3.*V->pvectsz[V->pclass_sz->index_pk_for_halo_bias];
        // }

      if (((V->pclass_sz->has_bk_ttg_at_z_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_bk_ttg_at_z_1h))
      || ((V->pclass_sz->has_bk_ttg_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_bk_ttg_at_z_2h))
      || ((V->pclass_sz->has_bk_ttg_at_z_3h == _TRUE_) && (index_md == V->pclass_sz->index_md_bk_ttg_at_z_3h))){

        evaluate_vrms2(V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);
        result *= V->pvectsz[V->pclass_sz->index_vrms2]/3./pow(_c_*1e-3,2.);


      }



        return result;
    }  

 if ((V->pclass_sz->has_sz_rates == _TRUE_) && (index_md == V->pclass_sz->index_md_szrates)){
   // if (V->pclass_sz->sz_verbose>0) printf("finnished mass integration for szrates.\n");
   return result*get_volume_at_z(V->pvectsz[V->pclass_sz->index_z],V->pba);
 }



  // exit(0);
  }// END MASS INTEGRATION
  // NOW MULTIPLY BY REDSHIFT DEPENDENT KERNELS

if (((V->pclass_sz->has_sz_2halo == _TRUE_) && (index_md == V->pclass_sz->index_md_2halo))
 || ((V->pclass_sz->has_gal_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_2h)) //## BB debug
 || ((V->pclass_sz->has_cib_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_cib_cib_2h))
 || ((V->pclass_sz->has_kSZ_kSZ_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_2h))
 || ((V->pclass_sz->has_tSZ_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_cib_2h))
 || ((V->pclass_sz->has_ngal_ngal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_ngal_2h))
 || ((V->pclass_sz->has_ngal_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_lens_2h))
 || ((V->pclass_sz->has_ngal_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_gallens_2h))
 // || ((V->pclass_sz->has_ngal_IA_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_IA_2h))
 || ((V->pclass_sz->has_ngal_tsz_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_tsz_2h))
 || ((V->pclass_sz->has_nlensmag_tsz_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_nlensmag_tsz_2h))
 || ((V->pclass_sz->has_nlensmag_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_nlensmag_gallens_2h))
 || ((V->pclass_sz->has_gallens_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_cib_2h))
 || ((V->pclass_sz->has_gal_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_cib_2h))
 || ((V->pclass_sz->has_lens_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_cib_2h))
 || ((V->pclass_sz->has_tSZ_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_gal_2h))
 // || ((V->pclass_sz->has_IA_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_IA_gal_2h))
 || ((V->pclass_sz->has_tau_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tau_gal_2h))
 || ((V->pclass_sz->has_tau_tau_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tau_tau_2h))
 || ((V->pclass_sz->has_gal_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lens_2h))
 || ((V->pclass_sz->has_gal_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lensmag_2h))
 || ((V->pclass_sz->has_tSZ_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_gallens_2h))
 || ((V->pclass_sz->has_gal_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gallens_2h))
 || ((V->pclass_sz->has_gallens_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_gallens_2h))
 || ((V->pclass_sz->has_gallens_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_lens_2h))
 || ((V->pclass_sz->has_lens_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lensmag_2h))
 || ((V->pclass_sz->has_gallens_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_lensmag_2h))
 || ((V->pclass_sz->has_lensmag_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lensmag_lensmag_2h))
 || ((V->pclass_sz->has_lens_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lens_2h))
 || ((V->pclass_sz->has_tSZ_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_lens_2h))
 || ((V->pclass_sz->has_tSZ_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_lensmag_2h))
 || ((V->pclass_sz->has_sz_m_y_y_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_m_y_y_2h))
 || ((V->pclass_sz->has_custom1_custom1_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_custom1_2h))
 || ((V->pclass_sz->has_custom1_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_lens_2h))
 || ((V->pclass_sz->has_custom1_tSZ_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_tSZ_2h))
 || ((V->pclass_sz->has_custom1_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_gal_2h))
 || ((V->pclass_sz->has_custom1_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_gallens_2h))
 || ((V->pclass_sz->has_custom1_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_cib_2h))
 // || ((V->pclass_sz->has_pk_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_at_z_2h))
 // || ((V->pclass_sz->has_pk_gg_at_z_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_pk_gg_at_z_2h))
 // || ((V->pclass_sz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_2h))
    ){


  int index_l = (int) V->pvectsz[V->pclass_sz->index_multipole];
  // V->pvectsz[V->pclass_sz->index_multipole_for_pk] = V->pclass_sz->ell[index_l];
  // evaluate_pk_at_ell_plus_one_half_over_chi(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);
  // double z = V->pvectsz[V->pclass_sz->index_z];
  //  double d_A = V->pvecback[V->pba->index_bg_ang_distance]*V->pba->h*(1.+z);
  //  double pk;

// double pkr;
// double fr = get2_pk_lin_at_k_and_z(//V->pvecback,//V->pvectsz,
//   &pkr,(V->pclass_sz->ell[index_l]+0.5)/d_A,z,V->pba,V->ppm,V->pnl,V->pclass_sz);
//   printf("k=%.3e z=%.3e pke=%.3e pklin2=%.3e pklin=%.3e fr=%.3e\n",
//          (V->pclass_sz->ell[index_l]+0.5)/chi,z,
//          V->pvectsz[V->pclass_sz->index_pk_for_halo_bias],
//          pkr,
//          pkp,
//          fr);
//
// if (fr == 1.)
//   exit(0);
// if (fr == 0.)
//   exit(0);
  // For all the above cases we multiply the linear matter power spectrum to the redshift integrand
  // evaluated at (ell+1/2)/Chi and redshift z
  if (V->pclass_sz->use_pknl_in_2hterms){
  result *= get_pk_nonlin_at_k_and_z((V->pclass_sz->ell[index_l]+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);//V->pvectsz[V->pclass_sz->index_pk_for_halo_bias];
  }
  else{
  result *= get_pk_lin_at_k_and_z((V->pclass_sz->ell[index_l]+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);//V->pvectsz[V->pclass_sz->index_pk_for_halo_bias];
  }

}
if (((V->pclass_sz->has_ngal_IA_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_IA_2h))
    || ((V->pclass_sz->has_IA_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_IA_gal_2h))
  ){
    int index_l = (int) V->pvectsz[V->pclass_sz->index_multipole];

    if (V->pclass_sz->use_pknl_in_2hterms_IA_only){
    result *= get_pk_nonlin_at_k_and_z((V->pclass_sz->ell[index_l]+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);//V->pvectsz[V->pclass_sz->index_pk_for_halo_bias];
    }
    else{
    result *= get_pk_lin_at_k_and_z((V->pclass_sz->ell[index_l]+0.5)/chi,z,V->pba,V->ppm,V->pnl,V->pclass_sz);//V->pvectsz[V->pclass_sz->index_pk_for_halo_bias];
    }

  }





// Power spectrum today : needed for ISW  stuff
if ( ((V->pclass_sz->has_isw_auto == _TRUE_) && (index_md == V->pclass_sz->index_md_isw_auto))
 ||  ((V->pclass_sz->has_isw_tsz == _TRUE_) && (index_md == V->pclass_sz->index_md_isw_tsz))
 ||  ((V->pclass_sz->has_isw_lens == _TRUE_) && (index_md == V->pclass_sz->index_md_isw_lens))
    ){

  // evaluate_pk_at_ell_plus_one_half_over_chi_today(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);
  int index_l = (int) V->pvectsz[V->pclass_sz->index_multipole];
  double pk1 =  get_pk_lin_at_k_and_z((V->pclass_sz->ell[index_l]+0.5)/chi,0.,V->pba,V->ppm,V->pnl,V->pclass_sz);

  // For all the above cases we add the linear matter power spectrum to the redshift integrand
  // evaluated at (ell+1/2)/Chi and redshift z=0
  result *= pk1; //V->pvectsz[V->pclass_sz->index_pk_for_halo_bias];


}

// galaxy radial kernel
if  (((V->pclass_sz->has_tSZ_gal_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_gal_1h))
  || ((V->pclass_sz->has_tSZ_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_gal_2h))
  || ((V->pclass_sz->has_custom1_gal_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_gal_1h))
  || ((V->pclass_sz->has_custom1_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_gal_2h))
  || ((V->pclass_sz->has_IA_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_IA_gal_2h))
  || ((V->pclass_sz->has_kSZ_kSZ_gal_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_1h))
  || ((V->pclass_sz->has_kSZ_kSZ_gal_1h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_1h_fft))
  || ((V->pclass_sz->has_kSZ_kSZ_gal_2h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_2h_fft))
  || ((V->pclass_sz->has_kSZ_kSZ_gal_3h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_3h_fft))
  || ((V->pclass_sz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_2h))
  || ((V->pclass_sz->has_kSZ_kSZ_gal_3h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_3h))
  || ((V->pclass_sz->has_kSZ_kSZ_gal_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_hf))
  || ((V->pclass_sz->has_gal_lens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lens_hf))
  || ((V->pclass_sz->has_tau_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tau_gal_2h))
  || ((V->pclass_sz->has_tau_gal_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_tau_gal_1h))
  || ((V->pclass_sz->has_gal_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lens_2h))
  || ((V->pclass_sz->has_gal_lens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lens_1h))
  || ((V->pclass_sz->has_gal_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_cib_2h))
  || ((V->pclass_sz->has_gal_cib_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_cib_1h))
  || ((V->pclass_sz->has_gal_lensmag_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lensmag_hf))
  || ((V->pclass_sz->has_gal_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lensmag_2h))
  || ((V->pclass_sz->has_gal_lensmag_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lensmag_1h))
  || ((V->pclass_sz->has_gal_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gallens_2h))
  || ((V->pclass_sz->has_gal_gallens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gallens_1h))
    ){
// multiply by radial kernel for galaxies
double Wg = radial_kernel_W_galaxy_at_z(V->pvecback,V->pvectsz,V->pba,V->pclass_sz);

result *= Wg/V->pvectsz[V->pclass_sz->index_chi2];
}

// gxg needs Wg^2:
if ( ((V->pclass_sz->has_gal_gal_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_1h))
   ||((V->pclass_sz->has_gal_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_2h))
   ||((V->pclass_sz->has_gal_gal_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_hf))
   || ((V->pclass_sz->has_gal_gal_lens_1h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_lens_1h_fft))
   || ((V->pclass_sz->has_gal_gal_lens_2h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_lens_2h_fft))
   || ((V->pclass_sz->has_gal_gal_lens_3h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_lens_3h_fft))

  ){
// multiply by radial kernel for galaxies (squared for gxg quantities)
double Wg = radial_kernel_W_galaxy_at_z(V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
 if (V->pclass_sz->sz_verbose>10)
    printf("multiply by radial kernel for galaxies z = %.5e Wg = %.5e\n",z,Wg);  

result *= pow(Wg/V->pvectsz[V->pclass_sz->index_chi2],2.);
}


if ( ((V->pclass_sz->has_ngal_ngal_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_ngal_1h))
   ||((V->pclass_sz->has_ngal_ngal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_ngal_2h))
   ||((V->pclass_sz->has_ngal_ngal_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_ngal_hf))
  ){
// multiply by radial kernel for galaxies (squared for gxg quantities)

int index_g = (int) V->pvectsz[V->pclass_sz->index_ngal_for_galaxy_profile];
int index_g_prime = (int) V->pvectsz[V->pclass_sz->index_ngal_prime_for_galaxy_profile];

// printf("getting kernels\n");
double Wg = radial_kernel_W_galaxy_ngal_at_z(index_g,
                                             V->pvecback,
                                             z,
                                             V->pba,
                                             V->pclass_sz);

double Wg_galprime = radial_kernel_W_galaxy_ngal_at_z(index_g_prime,
                                                      V->pvecback,
                                                      z,
                                                      V->pba,
                                                      V->pclass_sz);

// printf("getting kernels: wg = %.5e wgp = %.5e\n",Wg,Wg_galprime);

result *= Wg*Wg_galprime*pow(1/V->pvectsz[V->pclass_sz->index_chi2],2.);
}

if ( ((V->pclass_sz->has_ngal_lens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_lens_1h))
   ||((V->pclass_sz->has_ngal_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_lens_2h))
   ||((V->pclass_sz->has_ngal_lens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_lens_hf))
   ||((V->pclass_sz->has_ngal_nlensmag_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_nlensmag_hf))
   ||((V->pclass_sz->has_ngal_tsz_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_tsz_1h))
   ||((V->pclass_sz->has_ngal_tsz_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_tsz_2h))
   ||((V->pclass_sz->has_ngal_gallens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_gallens_1h))
   ||((V->pclass_sz->has_ngal_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_gallens_2h))
   ||((V->pclass_sz->has_ngal_IA_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_IA_2h))
  ){
// multiply by radial kernel for galaxies (squared for gxg quantities)

int index_g = (int) V->pvectsz[V->pclass_sz->index_ngal_for_galaxy_profile];
double Wg = radial_kernel_W_galaxy_ngal_at_z(index_g,
                                             V->pvecback,
                                             z,
                                             V->pba,
                                             V->pclass_sz);
// if (index_g<10)
// printf("index_g = %d Wg = %.5e z = %.5e\n",index_g,Wg,z);

result *= Wg/V->pvectsz[V->pclass_sz->index_chi2];
}

//  n lensing magification needs lensing kernel:

if ( ((V->pclass_sz->has_nlensmag_gallens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_nlensmag_gallens_1h))
   ||((V->pclass_sz->has_nlensmag_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_nlensmag_gallens_2h))
   ||((V->pclass_sz->has_nlensmag_tsz_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_nlensmag_tsz_1h))
   ||((V->pclass_sz->has_nlensmag_tsz_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_nlensmag_tsz_2h))
  ){

int index_g = (int) V->pvectsz[V->pclass_sz->index_ngal_for_galaxy_profile];
double Wg = radial_kernel_W_galaxy_lensing_magnification_nlensmag_at_z(index_g,
                                                                       V->pvectsz,
                                                                       z,
                                                                       V->pba,
                                                                       V->pclass_sz);

result *= Wg;
}

// lensing magification needs lensing kernel:

if (((V->pclass_sz->has_kSZ_kSZ_lensmag_1halo == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lensmag_1halo))
    ||((V->pclass_sz->has_tSZ_lensmag_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_lensmag_1h))
    ||((V->pclass_sz->has_tSZ_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_lensmag_2h))
    ||((V->pclass_sz->has_gal_lensmag_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lensmag_1h))
    ||((V->pclass_sz->has_gal_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lensmag_2h))
    ||((V->pclass_sz->has_lens_lensmag_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lensmag_1h))
    ||((V->pclass_sz->has_lens_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lensmag_2h))
    ||((V->pclass_sz->has_gallens_lensmag_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_lensmag_1h))
    ||((V->pclass_sz->has_gallens_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_lensmag_2h))
){

double Wg = radial_kernel_W_galaxy_lensing_magnification_at_z(z,V->pvectsz,V->pba,V->pclass_sz);
// printf("Wg  = %.5e\n",Wg);
result *= Wg;
}

if(
  ((V->pclass_sz->has_lensmag_lensmag_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_lensmag_lensmag_1h))
||((V->pclass_sz->has_lensmag_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lensmag_lensmag_2h))
){
double Wg = radial_kernel_W_galaxy_lensing_magnification_at_z(z,V->pvectsz,V->pba,V->pclass_sz);
result *= pow(Wg,2.);
}

// cmb lensing needs lensing kernel:

if (
    ((V->pclass_sz->has_gal_lens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lens_1h))
    ||((V->pclass_sz->has_gal_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_lens_2h))
    ||((V->pclass_sz->has_ngal_lens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_lens_1h))
    ||((V->pclass_sz->has_ngal_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_lens_2h))
    ||((V->pclass_sz->has_lens_lensmag_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lensmag_1h))
    ||((V->pclass_sz->has_lens_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lensmag_2h))
    ||((V->pclass_sz->has_lens_cib_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_cib_1h))
    ||((V->pclass_sz->has_lens_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_cib_2h))
    ||((V->pclass_sz->has_custom1_lens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_lens_1h))
    ||((V->pclass_sz->has_custom1_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_lens_2h))
    ||((V->pclass_sz->has_tSZ_lens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_lens_1h))
    ||((V->pclass_sz->has_tSZ_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_lens_2h))
    ||((V->pclass_sz->has_gallens_lens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_lens_1h))
    ||((V->pclass_sz->has_gallens_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_lens_2h))
    ||((V->pclass_sz->has_kSZ_kSZ_lens_1h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lens_1h_fft))
    ||((V->pclass_sz->has_kSZ_kSZ_lens_2h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lens_2h_fft))
    ||((V->pclass_sz->has_kSZ_kSZ_lens_3h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lens_3h_fft))
    ||((V->pclass_sz->has_gal_gal_lens_1h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_lens_1h_fft))
    ||((V->pclass_sz->has_gal_gal_lens_2h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_lens_2h_fft))
    ||((V->pclass_sz->has_gal_gal_lens_3h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gal_lens_3h_fft))
    ||((V->pclass_sz->has_kSZ_kSZ_lens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lens_hf)) // this is not correct??
){

if (((V->pclass_sz->has_kSZ_kSZ_lens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lens_hf))){
  printf("this seems not correct, check lensing kernel for kSZ_kSZ_lens_hf!\n");
  printf(" i think in the lensing hf case the kernel should be radial_kernel_W_lensing_at_z.\n");
  exit(0);
}
double Wg = radial_kernel_W_cmb_lensing_at_z(z,V->pvectsz,V->pba,V->pclass_sz);
result *= Wg;


}

if(
  ((V->pclass_sz->has_lens_lens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lens_1h))
||((V->pclass_sz->has_lens_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lens_2h))
||((V->pclass_sz->has_lens_lens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lens_hf))
){
double Wg = radial_kernel_W_cmb_lensing_at_z(z,V->pvectsz,V->pba,V->pclass_sz);
result *= pow(Wg,2.);
}


if(
  ((V->pclass_sz->has_custom1_custom1_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_custom1_1h))
||((V->pclass_sz->has_custom1_custom1_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_custom1_2h))
){
double Wg = get_radial_kernel_W_custom1_at_z(z,V->pclass_sz);
result *= pow(Wg,2.);
}

if(
  ((V->pclass_sz->has_custom1_lens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_lens_1h))
||((V->pclass_sz->has_custom1_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_lens_2h))
||((V->pclass_sz->has_custom1_tSZ_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_tSZ_1h))
||((V->pclass_sz->has_custom1_tSZ_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_tSZ_2h))
||((V->pclass_sz->has_custom1_cib_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_cib_1h))
||((V->pclass_sz->has_custom1_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_cib_2h))
||((V->pclass_sz->has_custom1_gal_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_gal_1h))
||((V->pclass_sz->has_custom1_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_gal_2h))
||((V->pclass_sz->has_custom1_gallens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_gallens_1h))
||((V->pclass_sz->has_custom1_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_gallens_2h))
){
double Wg = get_radial_kernel_W_custom1_at_z(z,V->pclass_sz);
result *= pow(Wg,1.);
}

// galaxy lensing lensing needs lensing kernel:

if (
    ((V->pclass_sz->has_gal_gallens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gallens_1h))
    ||((V->pclass_sz->has_gal_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_gallens_2h))
    ||((V->pclass_sz->has_gallens_cib_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_cib_1h))
    ||((V->pclass_sz->has_gallens_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_cib_2h))
    ||((V->pclass_sz->has_custom1_gallens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_gallens_1h))
    ||((V->pclass_sz->has_custom1_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_gallens_2h))
    ||((V->pclass_sz->has_tSZ_gallens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_gallens_1h))
    ||((V->pclass_sz->has_tSZ_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_gallens_2h))
    ||((V->pclass_sz->has_gallens_lens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_lens_1h))
    ||((V->pclass_sz->has_gallens_lens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_lens_2h))
    ||((V->pclass_sz->has_gallens_lensmag_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_lensmag_1h))
    ||((V->pclass_sz->has_gallens_lensmag_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_lensmag_2h))
    ||((V->pclass_sz->has_kSZ_kSZ_gallens_1h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gallens_1h_fft))
    ||((V->pclass_sz->has_kSZ_kSZ_gallens_2h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gallens_2h_fft))
    ||((V->pclass_sz->has_kSZ_kSZ_gallens_3h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gallens_3h_fft))
    ||((V->pclass_sz->has_kSZ_kSZ_gallens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gallens_hf))
    ||((V->pclass_sz->has_ngal_gallens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_gallens_1h))
    ||((V->pclass_sz->has_ngal_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_ngal_gallens_2h))
    ||((V->pclass_sz->has_nlensmag_gallens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_nlensmag_gallens_1h))
    ||((V->pclass_sz->has_nlensmag_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_nlensmag_gallens_2h))
){

double Wg = radial_kernel_W_galaxy_lensing_at_z(z,//V->pvectsz,V->pba,
                                                V->pclass_sz);
result *= Wg;
}

if(
  ((V->pclass_sz->has_gallens_gallens_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_gallens_1h))
||((V->pclass_sz->has_gallens_gallens_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_gallens_2h))
){
double Wg = radial_kernel_W_galaxy_lensing_at_z(z,//V->pvectsz,V->pba,
                                                V->pclass_sz);
result *= pow(Wg,2.);
}

if (
  ((V->pclass_sz->has_kSZ_kSZ_lens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lens_hf))
||((V->pclass_sz->has_kSZ_kSZ_gallens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gallens_hf))
){
  // result *= 1./(pow(3.*pow(V->pba->H0/V->pba->h,2)/2./V->pclass_sz->Rho_crit_0,-1)*pow((1.+z),1.)/chi);
  // result *= 1./(pow(3.*pow(V->pba->H0/V->pba->h,2)/2./V->pclass_sz->Rho_crit_0,-1)*pow((1.+z),1.)/chi/chi/chi);
  double Omega_m = V->pclass_sz->Omega_m_0;
  result *= 3.*pow(Omega_m,1.)*pow(V->pba->H0/V->pba->h,2)/2./chi*pow(1.+z,1.);
}

if (((V->pclass_sz->has_lens_lens_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_lens_hf))){
  double Omega_m = V->pclass_sz->Omega_m_0;
  result *= pow(3.*pow(Omega_m,1.)*pow(V->pba->H0/V->pba->h,2)/2./chi*pow(1.+z,1.),2.);
}

// multiply by velocity dispersion
if (
  ((V->pclass_sz->has_kSZ_kSZ_gal_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_1h))
 || ((V->pclass_sz->has_kSZ_kSZ_gal_1h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_1h_fft))
 || ((V->pclass_sz->has_kSZ_kSZ_gal_2h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_2h_fft))
 || ((V->pclass_sz->has_kSZ_kSZ_gal_3h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_3h_fft))
 || ((V->pclass_sz->has_kSZ_kSZ_gallens_1h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gallens_1h_fft))
 || ((V->pclass_sz->has_kSZ_kSZ_gallens_2h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gallens_2h_fft))
 || ((V->pclass_sz->has_kSZ_kSZ_gallens_3h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gallens_3h_fft))
 || ((V->pclass_sz->has_kSZ_kSZ_lens_1h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lens_1h_fft))
 || ((V->pclass_sz->has_kSZ_kSZ_lens_2h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lens_2h_fft))
 || ((V->pclass_sz->has_kSZ_kSZ_lens_3h_fft == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lens_3h_fft))
 || ((V->pclass_sz->has_kSZ_kSZ_gal_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_2h))
 || ((V->pclass_sz->has_kSZ_kSZ_gal_3h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_3h))
 || ((V->pclass_sz->has_kSZ_kSZ_tSZ_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_tSZ_1h))
 || ((V->pclass_sz->has_kSZ_kSZ_tSZ_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_tSZ_2h))
 || ((V->pclass_sz->has_kSZ_kSZ_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_1h))
 || ((V->pclass_sz->has_kSZ_kSZ_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_2h))
 || ((V->pclass_sz->has_kSZ_kSZ_tSZ_3h == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_tSZ_3h))
 // || (V->pclass_sz->has_kSZ_kSZ_gal_hf == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_gal_hf)
 || ((V->pclass_sz->has_kSZ_kSZ_lensmag_1halo == _TRUE_) && (index_md == V->pclass_sz->index_md_kSZ_kSZ_lensmag_1halo))
){
  // printf("evaluating vrms2\n");
  evaluate_vrms2(V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);
  result *= V->pvectsz[V->pclass_sz->index_vrms2]/3./pow(_c_*1e-3,2.);
}





// multiply by dsigma2_hsv
if ((V->pclass_sz->has_sz_cov_N_N_hsv == _TRUE_) && (index_md == V->pclass_sz->index_md_cov_N_N_hsv)){
  evaluate_sigma2_hsv(V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);
  result *= V->pvectsz[V->pclass_sz->index_sigma2_hsv];
}





if  (((V->pclass_sz->has_cib_cib_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_cib_cib_1h))
   ||((V->pclass_sz->has_cib_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_cib_cib_2h))
   ||((V->pclass_sz->has_cib_shotnoise == _TRUE_) && (index_md == V->pclass_sz->index_md_cib_shotnoise))
    ){

if (V->pclass_sz->use_maniyar_cib_model == 0){
// cib redshift kernel, see McCarthy and Madhavacheril 2020
 result *= 1./(1.+z)*1./(1.+z)*pow(1./V->pvectsz[V->pclass_sz->index_chi2],2.);
}
}


if   (((V->pclass_sz->has_tSZ_cib_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_cib_1h))
    ||((V->pclass_sz->has_tSZ_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_tSZ_cib_2h))
    ||((V->pclass_sz->has_gal_cib_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_cib_1h))
    ||((V->pclass_sz->has_gal_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gal_cib_2h))
    ||((V->pclass_sz->has_gallens_cib_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_cib_1h))
    ||((V->pclass_sz->has_gallens_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_gallens_cib_2h))
    ||((V->pclass_sz->has_custom1_cib_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_cib_1h))
    ||((V->pclass_sz->has_custom1_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_custom1_cib_2h))
    ||((V->pclass_sz->has_lens_cib_1h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_cib_1h))
    ||((V->pclass_sz->has_lens_cib_2h == _TRUE_) && (index_md == V->pclass_sz->index_md_lens_cib_2h))
    ||((V->pclass_sz->has_cib_monopole == _TRUE_) && (index_md == V->pclass_sz->index_md_cib_monopole))
    ){
if (V->pclass_sz->use_maniyar_cib_model == 0){
  result *= 1./(1.+z)*pow(1./V->pvectsz[V->pclass_sz->index_chi2],1.);
}
}


// else{

  // finally multiply by volume element Chi^2 dChi
  // result *= V->pvectsz[V->pclass_sz->index_chi2];


  // integrate w.r.t ln(1+z); dz =  (1+z)dln(1+z)
  // volume element in units h^-3 Mpc^3
  // volume = dv/(dzdOmega)
  // Chi^2 dChi = dV/(dzdOmega) dz
  // Chi^2 dChi = dV/(dzdOmega)*(1+z) dln(1+z)
  // dChi = (c/H) *(1+z) dln(1+z) ---> this is used
  // dChi = (c/H) dz
  // d/dchi = H/c d/dz
  // double H_over_c_in_h_over_Mpc = V->pvecback[V->pba->index_bg_H]/V->pba->h;
  //
  // printf("multiplying by volume %.3e %.3e\n",V->pvectsz[V->pclass_sz->index_chi2]/H_over_c_in_h_over_Mpc, get_volume_at_z(V->pvectsz[V->pclass_sz->index_z],V->pba));
  // result = (1.+V->pvectsz[V->pclass_sz->index_z])*result/H_over_c_in_h_over_Mpc;

  double volume = (1.+V->pvectsz[V->pclass_sz->index_z])*get_volume_at_z(V->pvectsz[V->pclass_sz->index_z],V->pba);
  

  
  result *= volume;
  // note : get_vol is c/H*chi2...dchi/dz*ch2




  if (isnan(result)||isinf(result)){
  printf("nan or inf in integrand redshift 1h, volume = %.3e\n",volume);
  exit(0);
  }
  return result;
// }

}


int integrate_over_redshift(struct background * pba,
                            struct nonlinear * pnl,
                            struct primordial * ppm,
                            struct perturbs * ppt,
                            struct class_sz_structure * pclass_sz,
                            double * Pvecback,
                            double * Pvectsz)
{


  double z_min = pclass_sz->z1SZ;
  double z_max = pclass_sz->z2SZ;


  struct Parameters_for_integrand_redshift V;
  V.pnl = pnl;
  V.ppm = ppm;
  V.ppt = ppt;
  V.pclass_sz = pclass_sz;
  V.pba = pba;
  V.pvectsz = Pvectsz;
  V.pvecback = Pvecback;

  void * params = &V;
  double r; //result of the integral

  double epsrel= pclass_sz->redshift_epsrel;
  double epsabs= pclass_sz->redshift_epsabs;
  int show_neval = pclass_sz->patterson_show_neval;


// hhere put the things that do not need integration over z:
  int index_md = (int) Pvectsz[pclass_sz->index_md];
if(_pk_at_z_1h_
|| _pk_at_z_2h_
|| _pk_gg_at_z_1h_
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
)
{
  r = integrand_redshift(log(1. + pclass_sz->z_for_pk_hm),params);
}
else if ( _szrates_ )
{
  if (pclass_sz->sz_verbose > 10) printf("evaluating rate at z = %.3e.\n",pclass_sz->szcat_z[(int)Pvectsz[pclass_sz->index_szrate]]);
    if (pclass_sz->szcat_z[(int)Pvectsz[pclass_sz->index_szrate]] <= 0){
          r = 1.; // contrinbutes to nothing to the lkl. since we take sum of ln(rates).
      }
    else {
        r = integrand_redshift(log(1. + pclass_sz->szcat_z[(int)Pvectsz[pclass_sz->index_szrate]]),params);

       if (r == 0.) r = 1e-300; // when snr cat is too high the integral becomes too small.
      }
}
else{
//   if (_gal_gal_hf_
//    || _ngal_ngal_hf_){
//
//   printf("z_min = %.4e z_max  = %.4e\n",z_min,z_max);
//   printf("currently in mode %d\n",_gal_gal_hf_);
//   exit(0);
//   z_min = pclass_sz->normalized_dndz_ngal_z[0];
//   z_max = .......
// }
if (pclass_sz->sz_verbose>10)
  printf("integrating over redshift\n");

  r = Integrate_using_Patterson_adaptive(log(1. + z_min), log(1. + z_max),
                                         epsrel, epsabs,
                                         integrand_redshift,
                                         params,show_neval);
if (pclass_sz->sz_verbose>10)
  printf("integrating over redshift got r = %.5e\n",r);
    }

  Pvectsz[pclass_sz->index_integral] = r;

  return _SUCCESS_;
}





struct Parameters_for_integrand_mass{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  struct perturbs * ppt;
  double * pvecback;
  double * pvectsz;
};


//
double integrand_mass(double logM, void *p){

  struct Parameters_for_integrand_mass *V = ((struct Parameters_for_integrand_mass *) p);

  double result = integrand_at_m_and_z(logM,
                                        V->pvecback,
                                        V->pvectsz,
                                        V->pba,
                                        V->ppm,
                                        V->pnl,
                                        V->ppt,
                                        V->pclass_sz);

  return result;

}



//Integration over the mass range at a given redshift
 double integrate_over_m_at_z(double * pvecback,
                             double * pvectsz,
                             struct background * pba,
                             struct nonlinear * pnl,
                             struct primordial * ppm,
                             struct perturbs * ppt,
                             struct class_sz_structure * pclass_sz)
{

if (pclass_sz->sz_verbose>10)
  printf("starting mass integral at z, preliminary calculations\n");
// if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_dndlnM)
//   return integrand_at_m_and_z(log(1e16),
//                               pvecback,
//                               pvectsz,
//                               pba,
//                               ppm,
//                               pnl,
//                               pclass_sz);


  double chi = sqrt(pvectsz[pclass_sz->index_chi2]);
  double epsrel=pclass_sz->mass_epsrel;
  double epsabs=pclass_sz->mass_epsabs;

  double m_min;
  double m_max;

  if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_cov_Y_N )
    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_cov_Y_N_next_order )
    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_cov_N_N )
    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_cov_N_N_hsv )){
    int index_m_1 = (int) pvectsz[pclass_sz->index_mass_bin_1];
    m_min = pclass_sz->cov_Y_N_mass_bin_edges[index_m_1];
    m_max = pclass_sz->cov_Y_N_mass_bin_edges[index_m_1+1];
  }

  else {
    m_min = pclass_sz->M1SZ;
    m_max = pclass_sz->M2SZ;

    if (pclass_sz->use_redshift_dependent_M_min){
      m_min = get_M_min_of_z(pvectsz[pclass_sz->index_z],pclass_sz);
      // printf("z = %.3e m_min = %.3e\n",pvectsz[pclass_sz->index_z],m_min);
    }
  }

  struct Parameters_for_integrand_mass V;
  V.pnl = pnl;
  V.ppm = ppm;
  V.ppt = ppt;
  V.pclass_sz = pclass_sz;
  V.pba = pba;
  V.pvectsz = pvectsz;
  V.pvecback = pvecback;
  void * params = &V;


  double r; //store result of mass integral

if (pclass_sz->sz_verbose>10){
  printf("starting mass integral at z, preliminary calculations done\n");
  }



  if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_cov_Y_N_next_order )){

  double r_cov_Y_N_next_order_1; // for cov_Y_N_next_order: first part of redshift integrand
  double r_cov_Y_N_next_order_2; // for cov_Y_N_next_order: second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate within the mass bin ('N' part)
  r_cov_Y_N_next_order_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                             epsrel, epsabs,
                                                             integrand_mass,
                                                             params,pclass_sz->patterson_show_neval);

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  m_min = pclass_sz->M1SZ;
  m_max = pclass_sz->M2SZ;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Y' part)
  r_cov_Y_N_next_order_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                             epsrel, epsabs,
                                                             integrand_mass,
                                                             params,pclass_sz->patterson_show_neval);
  r = r_cov_Y_N_next_order_1*r_cov_Y_N_next_order_2;
                                     }

  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_cov_N_N_hsv )){

  double r_cov_N_N_hsv_1; // for cov_Y_N_next_order: first part of redshift integrand
  double r_cov_N_N_hsv_2; // for cov_Y_N_next_order: second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate within the mass bin ('N' part)
  r_cov_N_N_hsv_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                     epsrel, epsabs,
                                                     integrand_mass,
                                                     params,pclass_sz->patterson_show_neval);

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;

  int index_m_2 = (int) pvectsz[pclass_sz->index_mass_bin_2];
  m_min = pclass_sz->cov_Y_N_mass_bin_edges[index_m_2];
  m_max = pclass_sz->cov_Y_N_mass_bin_edges[index_m_2+1];


  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Y' part)
  r_cov_N_N_hsv_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                     epsrel, epsabs,
                                                     integrand_mass,
                                                     params,pclass_sz->patterson_show_neval);
  r = r_cov_N_N_hsv_1*r_cov_N_N_hsv_2;
                                     }

  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_cov_Y_Y_ssc )){

  double r_cov_Y_Y_ssc_1; // for cov_Y_Y_ssc: first part of redshift integrand
  double r_cov_Y_Y_ssc_2; // for cov_Y_Y_ssc: second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  pvectsz[pclass_sz->index_multipole] =  pvectsz[pclass_sz->index_multipole_1];
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_cov_Y_Y_ssc_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                     epsrel, epsabs,
                                                     integrand_mass,
                                                     params,pclass_sz->patterson_show_neval);

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  pvectsz[pclass_sz->index_multipole] =  pvectsz[pclass_sz->index_multipole_2];
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Y' part)
  r_cov_Y_Y_ssc_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                                     epsrel, epsabs,
                                                     integrand_mass,
                                                     params,pclass_sz->patterson_show_neval);
  r = r_cov_Y_Y_ssc_1*r_cov_Y_Y_ssc_2;
                                     }

  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_tSZ_lens_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

if (pclass_sz->include_y_counterterms_in_yk){
 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }
}

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Phi' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }

  r = r_m_1*r_m_2;
                                     }

  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_tSZ_cib_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
                                     }

  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_ngal_ngal_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }

 int index_g = (int) pvectsz[pclass_sz->index_ngal_for_galaxy_profile];
 int index_g_prime = (int) pvectsz[pclass_sz->index_ngal_prime_for_galaxy_profile];
if (index_g_prime == index_g){
  r_m_2 = r_m_1;
}
else{

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }
}



  r = r_m_1*r_m_2;
                                     }



  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gallens_cib_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);
   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
                                     }

  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gal_cib_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);
   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
                                     }


  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_lens_cib_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Phi' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
                                     }



  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_custom1_gallens_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Phi' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
                                     }

  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_custom1_gal_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Phi' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
                                     }

  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_custom1_tSZ_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Phi' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
                                     }

  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_custom1_lens_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Phi' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
                                     }

  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_custom1_cib_2h )){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Phi' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('cib' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
                                     }

  else if ( ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_cib_cib_2h ) && (pvectsz[pclass_sz->index_frequency_for_cib_profile] != pvectsz[pclass_sz->index_frequency_prime_for_cib_profile])){

  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate for frequency nu
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate for frequency nu_prime
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }

  r = r_m_1*r_m_2;
                                     }


  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_IA_gal_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  // r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
  //                                          epsrel, epsabs,
  //                                          integrand_mass,
  //                                          params,pclass_sz->patterson_show_neval);
  r_m_1 = 1.;
  r_m_1 *= get_IA_of_z(pvectsz[pclass_sz->index_z],pba,pclass_sz);

//  if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
//    double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
//    double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
//    double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
//    double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
//    r_m_1 += bmin_umin;
//    // printf("counter terms done r_m_1\n");
// }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('galaxy' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
    }


  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_tSZ_gal_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_1\n");
}


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('galaxy' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
    }

  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_tSZ_gallens_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_1\n");
}


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('galaxy' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
    }



  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gal_gallens_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_1\n");
}


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('galaxy' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
    }




    else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gallens_lensmag_2h){
    double r_m_1; // first part of redshift integrand
    double r_m_2; // second part of redshift integrand

    pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
    V.pvectsz = pvectsz;
    params = &V;

    // integrate over the whole mass range ('Y' part)
    r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                             epsrel, epsabs,
                                             integrand_mass,
                                             params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


    pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
    V.pvectsz = pvectsz;
    params = &V;


    // integrate over the whole mass range ('galaxy' part)
    r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                             epsrel, epsabs,
                                             integrand_mass,
                                             params,pclass_sz->patterson_show_neval);

     if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_2 += bmin_umin;
     // printf("counter terms done r_m_2\n");
   }


    r = r_m_1*r_m_2;
      }

      else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_nlensmag_gallens_2h){
      double r_m_1; // first part of redshift integrand
      double r_m_2; // second part of redshift integrand

      pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
      V.pvectsz = pvectsz;
      params = &V;

      // integrate over the whole mass range ('Y' part)
      r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_mass,
                                               params,pclass_sz->patterson_show_neval);

     if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
       double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
       double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
       double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
       double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
       r_m_1 += bmin_umin;
       // printf("counter terms done r_m_1\n");
    }


      pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
      V.pvectsz = pvectsz;
      params = &V;


      // integrate over the whole mass range ('galaxy' part)
      r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_mass,
                                               params,pclass_sz->patterson_show_neval);

       if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
       double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
       double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
       double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
       double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
       r_m_2 += bmin_umin;
       // printf("counter terms done r_m_2\n");
     }


      r = r_m_1*r_m_2;
        }




  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_tSZ_lensmag_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_1\n");
}


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('galaxy' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
    }



  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gal_lensmag_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('gal' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

 if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_1 += bmin_umin;
   // printf("counter terms done r_m_1\n");
}


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Phi' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

  // corrected to match Mat M. consistency treatment in hmvec
  // in hmvec :
  // (integral+b1-consistency1)
  // with integral = r_m_2, b1 = 1, consistency1 = low mass part
  // r_m_2 = r_m_2+pclass_sz->Omega_m_0*pclass_sz->Rho_crit_0*pow(pvecback[pba->index_bg_ang_distance]*pba->h,-2.)/pvectsz[pclass_sz->index_lensing_Sigma_crit];

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }

  r = r_m_1*r_m_2;
  }

  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gal_lens_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('gal' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


  if (pclass_sz->include_g_counterterms_in_gk){
      if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
        double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
        double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
        double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
        double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
        r_m_1 += bmin_umin;
        // printf("counter terms done r_m_1\n");
      }
      }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Phi' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

  if (pclass_sz->include_k_counterterms_in_gk){
      if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
      double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
      double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
      double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
      double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
      r_m_2 += bmin_umin;
      // printf("counter terms done r_m_2\n");
        }
      }


  // corrected to match Mat M. consistency treatment in hmvec
  // r_m_2 = r_m_2 + pclass_sz->Omega_m_0*pclass_sz->Rho_crit_0*pow(pvecback[pba->index_bg_ang_distance]*pba->h,-2.)/pvectsz[pclass_sz->index_lensing_Sigma_crit];


  r = r_m_1*r_m_2;
  }
  // else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_ngal_lens_1h){
  //
  // }

  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_tau_gal_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('gal' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Phi' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  // corrected to match Mat M. consistency treatment in hmvec
  // r_m_2 = r_m_2 + pclass_sz->Omega_m_0*pclass_sz->Rho_crit_0*pow(pvecback[pba->index_bg_ang_distance]*pba->h,-2.)/pvectsz[pclass_sz->index_lensing_Sigma_crit];


  r = r_m_1*r_m_2;
  }
  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_ngal_lens_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('gal' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Phi' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;

 }
  r = r_m_1*r_m_2;
  }

  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_ngal_gallens_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('gal' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Phi' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;

 }
  r = r_m_1*r_m_2;
  }

  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_ngal_IA_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('Y' part)
  // r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
  //                                          epsrel, epsabs,
  //                                          integrand_mass,
  //                                          params,pclass_sz->patterson_show_neval);
  r_m_1 = 1.;
  r_m_1 *= get_IA_of_z(pvectsz[pclass_sz->index_z],pba,pclass_sz);

//  if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
//    double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
//    double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
//    double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
//    double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
//    r_m_1 += bmin_umin;
//    // printf("counter terms done r_m_1\n");
// }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('galaxy' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;
   // printf("counter terms done r_m_2\n");
 }


  r = r_m_1*r_m_2;
    }


  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_ngal_tsz_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('gal' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Phi' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;

  }
  r = r_m_1*r_m_2;
  }

  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_nlensmag_tsz_2h){
  double r_m_1; // first part of redshift integrand
  double r_m_2; // second part of redshift integrand

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;

  // integrate over the whole mass range ('gal' part)
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;


  // integrate over the whole mass range ('Phi' part)
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
   double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
   double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
   double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
   double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
   r_m_2 += bmin_umin;

  }
  r = r_m_1*r_m_2;
  }


  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_tSZ_2h){
  double r_m_11; // first part of redshift integrand
  double r_m_21; // second part of redshift integrand
  double r_m_12; // first part of redshift integrand
  double r_m_22; // second part of redshift integrand
  double r_m_13; // first part of redshift integrand
  double r_m_23; // second part of redshift integrand



  // r_m_11*r_m_21
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_11=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);



   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_11 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_21=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_21 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  // r_m_12*r_m_22
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_12=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_12 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_22=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_22 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  // r_m_13*r_m_23
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 5;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_13=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_13 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 6;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_23=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_23 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  int index_l = (int) pvectsz[pclass_sz->index_multipole];
  double l1 = pclass_sz->ell[index_l];
  double l2 = pclass_sz->bispectrum_lambda_k2*pclass_sz->ell[index_l];
  double l3 = pclass_sz->bispectrum_lambda_k3*pclass_sz->ell[index_l];
  double pk1 = get_pk_lin_at_k_and_z((l1+0.5)/chi,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  double pk2 = get_pk_lin_at_k_and_z((l2+0.5)/chi,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  double pk3 = get_pk_lin_at_k_and_z((l3+0.5)/chi,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);


  r = pk3*r_m_11*r_m_21  +  pk2*r_m_12*r_m_22  +  pk1*r_m_13*r_m_23;

  }


else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_tSZ_3h){

  double r_m_b1t1;
  double r_m_b1t2;
  double r_m_b1t3;

  double r_m_b2t1;
  double r_m_b2t2;
  double r_m_b2t3;

  double r_m_b1y1;
  double r_m_b1y2;
  double r_m_b1y3;

  double r_m_b2y1;
  double r_m_b2y2;
  double r_m_b2y3;

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1y1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1y1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1y2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1y2 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1y3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1y3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 5;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t2 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 6;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b1t3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 7;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2y1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b2y1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 8;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2y2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b2y2 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 9;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2y3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b2y3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 10;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b2t1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 11;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b2t2 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 12;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2t3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b2t3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  int index_l = (int) pvectsz[pclass_sz->index_multipole];
  double l1 = pclass_sz->ell[index_l];
  double l2 = pclass_sz->bispectrum_lambda_k2*pclass_sz->ell[index_l];
  double l3 = pclass_sz->bispectrum_lambda_k3*pclass_sz->ell[index_l];
  double k1 = (l1 + 0.5)/chi;
  double k2 = (l2 + 0.5)/chi;
  double k3 = (l3 + 0.5)/chi;
  double pk1 = 0.;
  double pk2 = 0.;
  double pk3 = 0.;

  pk1 = get_pk_lin_at_k_and_z(k1,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  pk2 = get_pk_lin_at_k_and_z(k2,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  pk3 = get_pk_lin_at_k_and_z(k3,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);

  double f2_123 = bispectrum_f2_kernel(k1,k2,k3);
  double f2_312 = bispectrum_f2_kernel(k3,k1,k2);
  double f2_231 = bispectrum_f2_kernel(k2,k3,k1);

  r = 2.*r_m_b1t1*r_m_b1t2*r_m_b1y3*f2_123*pk1*pk2
    +2.*r_m_b1t1*r_m_b1t2*r_m_b1y3*f2_312*pk3*pk1
    +2.*r_m_b1t1*r_m_b1t2*r_m_b1y3*f2_231*pk2*pk3

    +2.*r_m_b1t1*r_m_b1y2*r_m_b1t3*f2_123*pk1*pk2
    +2.*r_m_b1t1*r_m_b1y2*r_m_b1t3*f2_312*pk3*pk1
    +2.*r_m_b1t1*r_m_b1y2*r_m_b1t3*f2_231*pk2*pk3

    +2.*r_m_b1y1*r_m_b1t2*r_m_b1t3*f2_123*pk1*pk2
    +2.*r_m_b1y1*r_m_b1t2*r_m_b1t3*f2_312*pk3*pk1
    +2.*r_m_b1y1*r_m_b1t2*r_m_b1t3*f2_231*pk2*pk3

    +r_m_b1t1*r_m_b1t2*r_m_b2y3*pk1*pk2
    +r_m_b1t1*r_m_b2t2*r_m_b1y3*pk3*pk1
    +r_m_b2t1*r_m_b1t2*r_m_b1y3*pk2*pk3

    +r_m_b2t1*r_m_b1y2*r_m_b1t3*pk1*pk2
    +r_m_b1t1*r_m_b2y2*r_m_b1t3*pk3*pk1
    +r_m_b1t1*r_m_b1y2*r_m_b2t3*pk2*pk3

    +r_m_b2y1*r_m_b1t2*r_m_b1t3*pk1*pk2
    +r_m_b1y1*r_m_b2t2*r_m_b1t3*pk3*pk1
    +r_m_b1y1*r_m_b1t2*r_m_b2t3*pk2*pk3;


}




  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_tSZ_tSZ_tSZ_2h){
  double r_m_11; // first part of redshift integrand
  double r_m_21; // second part of redshift integrand
  double r_m_12; // first part of redshift integrand
  double r_m_22; // second part of redshift integrand
  double r_m_13; // first part of redshift integrand
  double r_m_23; // second part of redshift integrand



  // r_m_11*r_m_21
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_11=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);



   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_11 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_21=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_21 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  // r_m_12*r_m_22
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_12=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_12 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_22=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_22 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  // r_m_13*r_m_23
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 5;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_13=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_13 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 6;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_23=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_23 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  int index_l = (int) pvectsz[pclass_sz->index_multipole];
  double l1 = pclass_sz->ell[index_l];
  double l2 = pclass_sz->bispectrum_lambda_k2*pclass_sz->ell[index_l];
  double l3 = pclass_sz->bispectrum_lambda_k3*pclass_sz->ell[index_l];
  double pk1 = get_pk_lin_at_k_and_z((l1+0.5)/chi,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  double pk2 = get_pk_lin_at_k_and_z((l2+0.5)/chi,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  double pk3 = get_pk_lin_at_k_and_z((l3+0.5)/chi,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);


  r = pk3*r_m_11*r_m_21  +  pk2*r_m_12*r_m_22  +  pk1*r_m_13*r_m_23;

  }


else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_tSZ_tSZ_tSZ_3h){

  double r_m_b1y1;
  double r_m_b1y2;
  double r_m_b1y3;

  double r_m_b2y1;
  double r_m_b2y2;
  double r_m_b2y3;

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1y1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1y1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1y2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1y2 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1y3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1y3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2y1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b2y1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 5;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2y2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b2y2 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 6;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2y3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b2y3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }



  int index_l = (int) pvectsz[pclass_sz->index_multipole];
  double l1 = pclass_sz->ell[index_l];
  double l2 = pclass_sz->bispectrum_lambda_k2*pclass_sz->ell[index_l];
  double l3 = pclass_sz->bispectrum_lambda_k3*pclass_sz->ell[index_l];
  double k1 = (l1 + 0.5)/chi;
  double k2 = (l2 + 0.5)/chi;
  double k3 = (l3 + 0.5)/chi;
  double pk1 = 0.;
  double pk2 = 0.;
  double pk3 = 0.;

  pk1 = get_pk_lin_at_k_and_z(k1,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  pk2 = get_pk_lin_at_k_and_z(k2,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  pk3 = get_pk_lin_at_k_and_z(k3,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);

  double f2_123 = bispectrum_f2_kernel(k1,k2,k3);
  double f2_312 = bispectrum_f2_kernel(k3,k1,k2);
  double f2_231 = bispectrum_f2_kernel(k2,k3,k1);

  r = 2.*r_m_b1y1*r_m_b1y2*r_m_b1y3*f2_123*pk1*pk2
    +2.*r_m_b1y1*r_m_b1y2*r_m_b1y3*f2_312*pk3*pk1
    +2.*r_m_b1y1*r_m_b1y2*r_m_b1y3*f2_231*pk2*pk3

    +r_m_b1y1*r_m_b1y2*r_m_b2y3*pk1*pk2
    +r_m_b1y1*r_m_b2y2*r_m_b1y3*pk3*pk1
    +r_m_b2y1*r_m_b1y2*r_m_b1y3*pk2*pk3;


}



  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_gal_2h){
  double r_m_11; // first part of redshift integrand
  double r_m_21; // second part of redshift integrand
  double r_m_12; // first part of redshift integrand
  double r_m_22; // second part of redshift integrand
  double r_m_13; // first part of redshift integrand
  double r_m_23; // second part of redshift integrand



  // r_m_11*r_m_21
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_11=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);



   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_11 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_21=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_21 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  // r_m_12*r_m_22
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_12=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_12 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_22=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_22 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  // r_m_13*r_m_23
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 5;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_13=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_13 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 6;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_23=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_23 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  //int index_l_1 = (int) pvectsz[pclass_sz->index_multipole_1];
  int index_theta_1 = (int) pvectsz[pclass_sz->index_multipole_1];
  double theta_1 = pclass_sz->theta_kSZ2_gal_theta_grid[index_theta_1];
  int index_l_2 = (int) pvectsz[pclass_sz->index_multipole_2];
  int index_l_3 = (int) pvectsz[pclass_sz->index_multipole_3];
  double l2 = exp(pclass_sz->ell_kSZ2_gal_multipole_grid[index_l_2]);
  double l3 = pclass_sz->ell[index_l_3];
  double ell = l3;
  double ell_prime = l2;
  double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos(theta_1));
  // pvectsz[pclass_sz->index_multipole_for_pk] = l1;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk1 = pvectsz[pclass_sz->index_pk_for_halo_bias];
  double pk1 = get_pk_lin_at_k_and_z((l1+0.5)/chi,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);


  // pvectsz[pclass_sz->index_multipole_for_pk] = l2;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk2 = pvectsz[pclass_sz->index_pk_for_halo_bias];
  double pk2 = get_pk_lin_at_k_and_z((l2+0.5)/chi,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);


  // pvectsz[pclass_sz->index_multipole_for_pk] = l3;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk3 = pvectsz[pclass_sz->index_pk_for_halo_bias];
  double pk3 = get_pk_lin_at_k_and_z((l3+0.5)/chi,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);


  r = pk3*r_m_11*r_m_21  +  pk2*r_m_12*r_m_22  +  pk1*r_m_13*r_m_23;

  }




  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_gal_2h_fft){
  double r_m_11; // first part of redshift integrand
  double r_m_21; // second part of redshift integrand
  double r_m_12; // first part of redshift integrand
  double r_m_22; // second part of redshift integrand
  double r_m_13; // first part of redshift integrand
  double r_m_23; // second part of redshift integrand

  double z = pvectsz[pclass_sz->index_z];


  int index_l_3 = (int) pvectsz[pclass_sz->index_multipole];
  double l3 = pclass_sz->ell[index_l_3];
  // pvectsz[pclass_sz->index_multipole_for_pk] = l3;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk3 = pvectsz[pclass_sz->index_pk_for_halo_bias];
  double pk3 = get_pk_lin_at_k_and_z((l3+0.5)/chi,z,pba,ppm,pnl,pclass_sz);




  // r_m_11*r_m_21
  // pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_11=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_11 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  r_m_21 = pk3*get_psi_b1g_at_k_and_z((l3+0.5)/chi,z,pclass_sz);


  // printf("%.5e %.5e\n",r_m_11,r_m_21);



//// r_m_12 and r_m_22
//// no mass integral
//// apply convolution theorem


  /// set-up:

// double l_min = 1e-2;
// double l_max = 2e5; // this is a precision parameter
double l_min = pclass_sz->l_min_samp_fftw;
double l_max = pclass_sz->l_max_samp_fftw; // this is a precision parameter
// tabulate the integrand in the "l" dimension:
const int N = pclass_sz->N_samp_fftw;
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
fl = get_ksz_filter_at_l(l,pclass_sz);

Pk1[ik] = fl*get_psi_b1gt_at_k1_k2_and_z((l3+0.5)/chi,(l+0.5)/chi,z,pclass_sz);
if (isnan(Pk1[ik])||isinf(Pk1[ik])){
  printf("fft 2h : z %.3e k3 %.4e k' %.4e\n",z,(l3+0.5)/chi,(l+0.5)/chi);
  exit(0);
}

// pvectsz[pclass_sz->index_multipole_for_pk] = l;
// evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
double pkl = get_pk_lin_at_k_and_z((l+0.5)/chi,z,pba,ppm,pnl,pclass_sz);//pvectsz[pclass_sz->index_pk_for_halo_bias];
Pk2[ik] = fl*pkl*get_psi_b1t_at_k_and_z((l+0.5)/chi,z,pclass_sz);
// if(l>3e3)
  // printf("l = %.5e pk2 = %.5e\n",l,Pk2[ik]);
}
// printf("k pk done\n");
// exit(0);

double rp[N], xi1[N], xi2[N], xi12[N];

// go to Fourier space:
xi2pk(N,k,Pk1,rp,xi1,pclass_sz);
xi2pk(N,k,Pk2,rp,xi2,pclass_sz);
for (ik=0; ik<N; ik++){
  // if (isnan(xi1[ik]) || isinf(r)){
// convolution:
xi12[ik] = xi1[ik]*xi2[ik];
}
// printf("xi pi done\n");

// move back to position space:
pk2xi(N,rp,xi12,k,Pkr,pclass_sz);


// evaluate at l3


double f_psi_f_psi = pwl_value_1d(N,lnk,Pkr,log(l3));

r_m_12 = f_psi_f_psi;
r_m_22 = 1.;

r_m_13 = r_m_12;
r_m_23 = 1.;

r = r_m_11*r_m_21 +  r_m_12*r_m_22  +  r_m_13*r_m_23;
// printf("xi pd done r=%.5e\n",r);
if (isnan(r) || isinf(r)){
  // check transform
  for (ik=0; ik<10; ik++){
  printf("rp = %.3e xi1 = %.3e k = %.3e Pk1 = %.5e\n",rp[ik],xi1[ik],k[ik],Pk1[ik]);
  printf("rp = %.3e xi2 = %.3e k = %.3e Pk2 = %.5e\n",rp[ik],xi2[ik],k[ik],Pk2[ik]);
  printf("rp = %.3e xi12 = %.3e k = %.3e Pkr = %.5e\n",rp[ik],xi12[ik],k[ik],Pkr[ik]);
  }
  printf("in kSZ_kSZ_gal_2h_fft k %.3e z %.3e r_m_11 %.5e r_m_12 %.5e r_m_21 %.5e\n",(l3+0.5)/chi,z,r_m_11,r_m_12,r_m_21);
  exit(0);
}
  }



  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_gal_3h_fft){
  double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12;

  double z = pvectsz[pclass_sz->index_z];


  int index_l_3 = (int) pvectsz[pclass_sz->index_multipole];
  double l3 = pclass_sz->ell[index_l_3];
  // pvectsz[pclass_sz->index_multipole_for_pk] = l3;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk3 = pvectsz[pclass_sz->index_pk_for_halo_bias];
  double pk3 = get_pk_lin_at_k_and_z((l3+0.5)/chi,z,pba,ppm,pnl,pclass_sz);
  double psi_bg = get_psi_b1g_at_k_and_z((l3+0.5)/chi,z,pclass_sz);
  double psi_b2g =get_psi_b2g_at_k_and_z((l3+0.5)/chi,z,pclass_sz);
  // double psi_b2t = get_psi_b2t_at_k_and_z(l3,z,pclass_sz);

  // printf("%.5e %.5e\n",r_m_11,r_m_21);



//// r_m_12 and r_m_22
//// no mass integral
//// apply convolution theorem


  /// set-up:

  double l_min = pclass_sz->l_min_samp_fftw;
  double l_max = pclass_sz->l_max_samp_fftw; // this is a precision parameter
// tabulate the integrand in the "l" dimension:
const int N = pclass_sz->N_samp_fftw;
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
double l;
double pkl=0.;

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
// pvectsz[pclass_sz->index_multipole_for_pk] = l;
// pvectsz[pclass_sz->index_pk_for_halo_bias] = 0.;
// evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
// pkl = pvectsz[pclass_sz->index_pk_for_halo_bias];
pkl = get_pk_lin_at_k_and_z((l+0.5)/chi,z,pba,ppm,pnl,pclass_sz);
fl = get_ksz_filter_at_l(l,pclass_sz);
// if ((l+0.5)/chi>1e-2) fl = 0.;
psi_bt = get_psi_b1t_at_k_and_z((l+0.5)/chi,z,pclass_sz);

psi_b2t = get_psi_b2t_at_k_and_z((l+0.5)/chi,z,pclass_sz);
// l = 1.;

pk_phi_0[ik] = fl*psi_bt;
pk_phi_m2[ik] = pow((l+0.5)/chi,-2)*fl*psi_bt;
pk_phi_4[ik] = pow((l+0.5)/chi,4)*fl*psi_bt;
pk_phi_2[ik] = pow((l+0.5)/chi,2)*fl*psi_bt;

pk_tilde_phi_0[ik] = fl*pkl*psi_bt;
pk_tilde_phi_m2[ik] = pow((l+0.5)/chi,-2)*fl*pkl*psi_bt;
pk_tilde_phi_2[ik] = pow((l+0.5)/chi,2)*fl*pkl*psi_bt;
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
xi2pk(N,k,pk_phi_0,rp,xi_phi_0,pclass_sz);
xi2pk(N,k,pk_phi_2,rp,xi_phi_2,pclass_sz);
xi2pk(N,k,pk_phi_m2,rp,xi_phi_m2,pclass_sz);
xi2pk(N,k,pk_tilde_phi_0,rp,xi_tilde_phi_0,pclass_sz);
xi2pk(N,k,pk_tilde_phi_b20,rp,xi_tilde_phi_b20,pclass_sz);
xi2pk(N,k,pk_tilde_phi_2,rp,xi_tilde_phi_2,pclass_sz);
xi2pk(N,k,pk_tilde_phi_m2,rp,xi_tilde_phi_m2,pclass_sz);
xi2pk(N,k,pk_phi_4,rp,xi_phi_4,pclass_sz);

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
pk2xi(N,rp,t1_xi12,k,t1_Pkr,pclass_sz);
pk2xi(N,rp,t2_xi12,k,t2_Pkr,pclass_sz);
pk2xi(N,rp,t3_xi12,k,t3_Pkr,pclass_sz);
pk2xi(N,rp,t4_xi12,k,t4_Pkr,pclass_sz);
pk2xi(N,rp,t5_xi12,k,t5_Pkr,pclass_sz);
pk2xi(N,rp,t6_xi12,k,t6_Pkr,pclass_sz);
pk2xi(N,rp,t7_xi12,k,t7_Pkr,pclass_sz);
pk2xi(N,rp,t8_xi12,k,t8_Pkr,pclass_sz);
pk2xi(N,rp,t9_xi12,k,t9_Pkr,pclass_sz);
pk2xi(N,rp,t10_xi12,k,t10_Pkr,pclass_sz);
pk2xi(N,rp,t11_xi12,k,t11_Pkr,pclass_sz);
pk2xi(N,rp,t12_xi12,k,t12_Pkr,pclass_sz);


r = 10./14.*psi_bg*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
   -5./7.*psi_bg*pwl_value_1d(N,lnk,t2_Pkr,log(l3))
   +3./7.*psi_bg*pow((l3+0.5)/chi,2.)*pwl_value_1d(N,lnk,t3_Pkr,log(l3))
   +1./7.*psi_bg*pow((l3+0.5)/chi,4.)*pwl_value_1d(N,lnk,t4_Pkr,log(l3))
   // b2 terms:
   +psi_b2g*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
   +2.*psi_bg*pk3*pwl_value_1d(N,lnk,t12_Pkr,log(l3))

   +10./14.*pk3*psi_bg*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
   +3./14.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
   +3./14.*pk3*psi_bg*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
   -5./14.*pk3*psi_bg*pow((l3+0.5)/chi,2.)*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
   -5./14.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
   +1./7.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t10_Pkr,log(l3))

   +10./14.*pk3*psi_bg*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
   +3./14.*pk3*psi_bg*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
   +3./14.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
   -5./14.*pk3*psi_bg*pow((l3+0.5)/chi,2.)*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
   -5./14.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
   +1./7.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t10_Pkr,log(l3));



   //+2./7.*pwl_value_1d(N,lnk,t2_Pkr,log(l3));
// r = 19./7.*psi_bg*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
//     +9./7.*psi_bg*pwl_value_1d(N,lnk,t2_Pkr,log(l3))
//     -11./7.*pow((l3+0.5)/chi,2.)*psi_bg*pwl_value_1d(N,lnk,t3_Pkr,log(l3))
//     +1./7.*psi_bg*pow((l3+0.5)/chi,4.)*pwl_value_1d(N,lnk,t4_Pkr,log(l3));
// +24./7.*psi_bg*pk3*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
// +2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
// +2./7.*pow(l3,2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
// -4./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
// -4./7.*psi_bg*pk3*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
// +2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t10_Pkr,log(l3))
// // b2 terms:
 // psi_b2g*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
// +psi_bg*pwl_value_1d(N,lnk,t11_Pkr,log(l3))
// +psi_bg*pk3*pwl_value_1d(N,lnk,t12_Pkr,log(l3));

// r = (psi_b2g+19./7.*psi_bg)*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
// +9./7.*psi_bg*pwl_value_1d(N,lnk,t2_Pkr,log(l3))
// -11./7.*pow(l3,2.)*psi_bg*pwl_value_1d(N,lnk,t3_Pkr,log(l3))
// +1./7.*psi_bg*pow(l3,4.)*pwl_value_1d(N,lnk,t4_Pkr,log(l3))
// +(2.*psi_b2g+24./7.*psi_bg)*pk3*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
// +2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
// +2./7.*pow(l3,2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
// -4./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
// -4./7.*psi_bg*pk3*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
// +2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t10_Pkr,log(l3));




if (isnan(r) || isinf(r)){
  printf("nan in bispectrum TTG ffts\n");
  printf("pk3 = %.3e\n",pk3);
  printf("psi_bg = %.3e\n",psi_bg);
  printf("psi_b2g = %.3e\n",psi_b2g);
  printf("t1_Pkr = %.3e\n",pwl_value_1d(N,lnk,t1_Pkr,log(l3)));
  printf("t2_Pkr = %.3e\n",pwl_value_1d(N,lnk,t2_Pkr,log(l3)));
  printf("t3_Pkr = %.3e\n",pwl_value_1d(N,lnk,t3_Pkr,log(l3)));
  printf("t4_Pkr = %.3e\n",pwl_value_1d(N,lnk,t4_Pkr,log(l3)));
  printf("t5_Pkr = %.3e\n",pwl_value_1d(N,lnk,t5_Pkr,log(l3)));
  printf("t6_Pkr = %.3e\n",pwl_value_1d(N,lnk,t6_Pkr,log(l3)));
  printf("t7_Pkr = %.3e\n",pwl_value_1d(N,lnk,t7_Pkr,log(l3)));
  printf("t8_Pkr = %.3e\n",pwl_value_1d(N,lnk,t8_Pkr,log(l3)));
  printf("t9_Pkr = %.3e\n",pwl_value_1d(N,lnk,t9_Pkr,log(l3)));
  printf("t10_Pkr = %.3e\n",pwl_value_1d(N,lnk,t10_Pkr,log(l3)));
  printf("t11_Pkr = %.3e\n",pwl_value_1d(N,lnk,t11_Pkr,log(l3)));
  printf("t12_Pkr = %.3e\n",pwl_value_1d(N,lnk,t12_Pkr,log(l3)));
  exit(0);
}

  }

  else if (
          ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_gallens_2h_fft)
       || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_lens_2h_fft)
){
  double r_m_11; // first part of redshift integrand
  double r_m_21; // second part of redshift integrand
  double r_m_12; // first part of redshift integrand
  double r_m_22; // second part of redshift integrand
  double r_m_13; // first part of redshift integrand
  double r_m_23; // second part of redshift integrand

  double z = pvectsz[pclass_sz->index_z];


  int index_l_3 = (int) pvectsz[pclass_sz->index_multipole];
  double l3 = pclass_sz->ell[index_l_3];
  // pvectsz[pclass_sz->index_multipole_for_pk] = l3;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk3 = pvectsz[pclass_sz->index_pk_for_halo_bias];
  double pk3 = get_pk_lin_at_k_and_z((l3+0.5)/chi,z,pba,ppm,pnl,pclass_sz);




  // r_m_11*r_m_21
  // pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_11=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);


   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_11 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  r_m_21 = pk3*get_psi_b1kg_at_k_and_z((l3+0.5)/chi,z,pclass_sz);


  // printf("%.5e %.5e\n",r_m_11,r_m_21);



//// r_m_12 and r_m_22
//// no mass integral
//// apply convolution theorem


  /// set-up:

// double l_min = 1e-2;
// double l_max = 2e5; // this is a precision parameter
double l_min = pclass_sz->l_min_samp_fftw;
double l_max = pclass_sz->l_max_samp_fftw; // this is a precision parameter
// tabulate the integrand in the "l" dimension:
const int N = pclass_sz->N_samp_fftw;
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
fl = get_ksz_filter_at_l(l,pclass_sz);

Pk1[ik] = fl*get_psi_b1kgt_at_k1_k2_and_z((l3+0.5)/chi,(l+0.5)/chi,z,pclass_sz);
if (isnan(Pk1[ik])||isinf(Pk1[ik])){
  printf("fft 2h : z %.3e k3 %.4e k' %.4e\n",z,(l3+0.5)/chi,(l+0.5)/chi);
  exit(0);
}

// pvectsz[pclass_sz->index_multipole_for_pk] = l;
// evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
double pkl = get_pk_lin_at_k_and_z((l+0.5)/chi,z,pba,ppm,pnl,pclass_sz);//pvectsz[pclass_sz->index_pk_for_halo_bias];
Pk2[ik] = fl*pkl*get_psi_b1t_at_k_and_z((l+0.5)/chi,z,pclass_sz);
// if(l>3e3)
  // printf("k = %.5e pk = %.5e\n",l,Pk2[ik]);
}
// printf("k pk done\n");

double rp[N], xi1[N], xi2[N], xi12[N];

// go to Fourier space:
xi2pk(N,k,Pk1,rp,xi1,pclass_sz);
xi2pk(N,k,Pk2,rp,xi2,pclass_sz);
for (ik=0; ik<N; ik++){
// convolution:
xi12[ik] = xi1[ik]*xi2[ik];
}
// printf("xi pi done\n");

// move back to position space:
pk2xi(N,rp,xi12,k,Pkr,pclass_sz);

// evaluate at l3
double f_psi_f_psi = pwl_value_1d(N,lnk,Pkr,log(l3));

r_m_12 = f_psi_f_psi;
r_m_22 = 1.;

r_m_13 = r_m_12;
r_m_23 = 1.;

r = r_m_11*r_m_21 +  r_m_12*r_m_22  +  r_m_13*r_m_23;
// printf("xi pd done r=%.5e\n",r);
if (isnan(r) || isinf(r)){
  printf("in kSZ_kSZ_lens_2h_fft k %.3e z %.3e r_m_11 %.5e r_m_12 %.5e r_m_21 %.5e\n",(l3+0.5)/chi,z,r_m_11,r_m_12,r_m_21);
  exit(0);
}
  }

  else if (
          ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gal_gal_lens_2h_fft)
       //|| ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_lens_2h_fft)
){
  double r_m_11; // first part of redshift integrand
  double r_m_21; // second part of redshift integrand
  double r_m_12; // first part of redshift integrand
  double r_m_22; // second part of redshift integrand
  double r_m_13; // first part of redshift integrand
  double r_m_23; // second part of redshift integrand

  double z = pvectsz[pclass_sz->index_z];


  int index_l_3 = (int) pvectsz[pclass_sz->index_multipole];
  double l3 = pclass_sz->ell[index_l_3];
  // pvectsz[pclass_sz->index_multipole_for_pk] = l3;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk3 = pvectsz[pclass_sz->index_pk_for_halo_bias];
  double pk3 = get_pk_lin_at_k_and_z((l3+0.5)/chi,z,pba,ppm,pnl,pclass_sz);




  // r_m_11*r_m_21
  // pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_11=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

  if (isnan(r_m_11)){
    printf("rm11 nan in gal_gal_lens_2h\n");
  }

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_11 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  r_m_21 = pk3*get_psi_b1kg_at_k_and_z((l3+0.5)/chi,z,pclass_sz);


  // printf("%.5e %.5e\n",r_m_11,r_m_21);



//// r_m_12 and r_m_22
//// no mass integral
//// apply convolution theorem


  /// set-up:

// double l_min = 1e-2;
// double l_max = 2e5; // this is a precision parameter
double l_min = pclass_sz->l_min_samp_fftw;
double l_max = pclass_sz->l_max_samp_fftw; // this is a precision parameter
// tabulate the integrand in the "l" dimension:
const int N = pclass_sz->N_samp_fftw;
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
fl = get_ksz_filter_at_l(l,pclass_sz);

Pk1[ik] = fl*get_psi_b1kgg_at_k1_k2_and_z((l3+0.5)/chi,(l+0.5)/chi,z,pclass_sz);
if (isnan(Pk1[ik])||isinf(Pk1[ik])){
  printf("fft 2h : z %.3e k3 %.4e k' %.4e\n",z,(l3+0.5)/chi,(l+0.5)/chi);
  exit(0);
}

// pvectsz[pclass_sz->index_multipole_for_pk] = l;
// evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
double pkl = get_pk_lin_at_k_and_z((l+0.5)/chi,z,pba,ppm,pnl,pclass_sz);//pvectsz[pclass_sz->index_pk_for_halo_bias];
Pk2[ik] = fl*pkl*get_psi_b1g_at_k_and_z((l+0.5)/chi,z,pclass_sz);
// if(l>3e3)
  // printf("k = %.5e pk = %.5e\n",l,Pk2[ik]);
}
// printf("k pk done\n");

double rp[N], xi1[N], xi2[N], xi12[N];

// go to Fourier space:
xi2pk(N,k,Pk1,rp,xi1,pclass_sz);
xi2pk(N,k,Pk2,rp,xi2,pclass_sz);
for (ik=0; ik<N; ik++){
// convolution:
xi12[ik] = xi1[ik]*xi2[ik];
}
// printf("xi pi done\n");

// move back to position space:
pk2xi(N,rp,xi12,k,Pkr,pclass_sz);

// evaluate at l3
double f_psi_f_psi = pwl_value_1d(N,lnk,Pkr,log(l3));

r_m_12 = f_psi_f_psi;
r_m_22 = 1.;

r_m_13 = r_m_12;
r_m_23 = 1.;

r = r_m_11*r_m_21 +  r_m_12*r_m_22  +  r_m_13*r_m_23;
// printf("xi pd done r=%.5e\n",r);
if (isnan(r) || isinf(r)){
  printf("in  gal_gal_lens_2h k %.3e z %.3e r_m_11 %.5e r_m_12 %.5e r_m_21 %.5e\n",(l3+0.5)/chi,z,r_m_11,r_m_12,r_m_21);
  exit(0);
}
  }



  else if (
    ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_gallens_3h_fft)
  ||((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_lens_3h_fft)

  ){
  double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12;

  double z = pvectsz[pclass_sz->index_z];


  int index_l_3 = (int) pvectsz[pclass_sz->index_multipole];
  double l3 = pclass_sz->ell[index_l_3];
  // pvectsz[pclass_sz->index_multipole_for_pk] = l3;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk3 = pvectsz[pclass_sz->index_pk_for_halo_bias];
  double pk3 = get_pk_lin_at_k_and_z((l3+0.5)/chi,z,pba,ppm,pnl,pclass_sz);
  double psi_bg = get_psi_b1kg_at_k_and_z((l3+0.5)/chi,z,pclass_sz);
  double psi_b2g =get_psi_b2kg_at_k_and_z((l3+0.5)/chi,z,pclass_sz);
  // double psi_b2t = get_psi_b2t_at_k_and_z(l3,z,pclass_sz);

  // printf("%.5e %.5e\n",r_m_11,r_m_21);



//// r_m_12 and r_m_22
//// no mass integral
//// apply convolution theorem


  /// set-up:

  double l_min = pclass_sz->l_min_samp_fftw;
  double l_max = pclass_sz->l_max_samp_fftw; // this is a precision parameter
// tabulate the integrand in the "l" dimension:
const int N = pclass_sz->N_samp_fftw;
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
double l;
double pkl=0.;

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
// pvectsz[pclass_sz->index_multipole_for_pk] = l;
// pvectsz[pclass_sz->index_pk_for_halo_bias] = 0.;
// evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
// pkl = pvectsz[pclass_sz->index_pk_for_halo_bias];
pkl = get_pk_lin_at_k_and_z((l+0.5)/chi,z,pba,ppm,pnl,pclass_sz);
fl = get_ksz_filter_at_l(l,pclass_sz);
// if ((l+0.5)/chi>1e-2) fl = 0.;
psi_bt = get_psi_b1t_at_k_and_z((l+0.5)/chi,z,pclass_sz);

psi_b2t = get_psi_b2t_at_k_and_z((l+0.5)/chi,z,pclass_sz);
// l = 1.;

pk_phi_0[ik] = fl*psi_bt;
pk_phi_m2[ik] = pow((l+0.5)/chi,-2)*fl*psi_bt;
pk_phi_4[ik] = pow((l+0.5)/chi,4)*fl*psi_bt;
pk_phi_2[ik] = pow((l+0.5)/chi,2)*fl*psi_bt;

pk_tilde_phi_0[ik] = fl*pkl*psi_bt;
pk_tilde_phi_m2[ik] = pow((l+0.5)/chi,-2)*fl*pkl*psi_bt;
pk_tilde_phi_2[ik] = pow((l+0.5)/chi,2)*fl*pkl*psi_bt;
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
xi2pk(N,k,pk_phi_0,rp,xi_phi_0,pclass_sz);
xi2pk(N,k,pk_phi_2,rp,xi_phi_2,pclass_sz);
xi2pk(N,k,pk_phi_m2,rp,xi_phi_m2,pclass_sz);
xi2pk(N,k,pk_tilde_phi_0,rp,xi_tilde_phi_0,pclass_sz);
xi2pk(N,k,pk_tilde_phi_b20,rp,xi_tilde_phi_b20,pclass_sz);
xi2pk(N,k,pk_tilde_phi_2,rp,xi_tilde_phi_2,pclass_sz);
xi2pk(N,k,pk_tilde_phi_m2,rp,xi_tilde_phi_m2,pclass_sz);
xi2pk(N,k,pk_phi_4,rp,xi_phi_4,pclass_sz);

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
pk2xi(N,rp,t1_xi12,k,t1_Pkr,pclass_sz);
pk2xi(N,rp,t2_xi12,k,t2_Pkr,pclass_sz);
pk2xi(N,rp,t3_xi12,k,t3_Pkr,pclass_sz);
pk2xi(N,rp,t4_xi12,k,t4_Pkr,pclass_sz);
pk2xi(N,rp,t5_xi12,k,t5_Pkr,pclass_sz);
pk2xi(N,rp,t6_xi12,k,t6_Pkr,pclass_sz);
pk2xi(N,rp,t7_xi12,k,t7_Pkr,pclass_sz);
pk2xi(N,rp,t8_xi12,k,t8_Pkr,pclass_sz);
pk2xi(N,rp,t9_xi12,k,t9_Pkr,pclass_sz);
pk2xi(N,rp,t10_xi12,k,t10_Pkr,pclass_sz);
pk2xi(N,rp,t11_xi12,k,t11_Pkr,pclass_sz);
pk2xi(N,rp,t12_xi12,k,t12_Pkr,pclass_sz);


r = 10./14.*psi_bg*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
   -5./7.*psi_bg*pwl_value_1d(N,lnk,t2_Pkr,log(l3))
   +3./7.*psi_bg*pow((l3+0.5)/chi,2.)*pwl_value_1d(N,lnk,t3_Pkr,log(l3))
   +1./7.*psi_bg*pow((l3+0.5)/chi,4.)*pwl_value_1d(N,lnk,t4_Pkr,log(l3))
   // b2 terms:
   +psi_b2g*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
   +2.*psi_bg*pk3*pwl_value_1d(N,lnk,t12_Pkr,log(l3))

   +10./14.*pk3*psi_bg*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
   +3./14.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
   +3./14.*pk3*psi_bg*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
   -5./14.*pk3*psi_bg*pow((l3+0.5)/chi,2.)*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
   -5./14.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
   +1./7.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t10_Pkr,log(l3))

   +10./14.*pk3*psi_bg*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
   +3./14.*pk3*psi_bg*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
   +3./14.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
   -5./14.*pk3*psi_bg*pow((l3+0.5)/chi,2.)*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
   -5./14.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
   +1./7.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t10_Pkr,log(l3));



   //+2./7.*pwl_value_1d(N,lnk,t2_Pkr,log(l3));
// r = 19./7.*psi_bg*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
//     +9./7.*psi_bg*pwl_value_1d(N,lnk,t2_Pkr,log(l3))
//     -11./7.*pow((l3+0.5)/chi,2.)*psi_bg*pwl_value_1d(N,lnk,t3_Pkr,log(l3))
//     +1./7.*psi_bg*pow((l3+0.5)/chi,4.)*pwl_value_1d(N,lnk,t4_Pkr,log(l3));
// +24./7.*psi_bg*pk3*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
// +2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
// +2./7.*pow(l3,2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
// -4./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
// -4./7.*psi_bg*pk3*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
// +2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t10_Pkr,log(l3))
// // b2 terms:
 // psi_b2g*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
// +psi_bg*pwl_value_1d(N,lnk,t11_Pkr,log(l3))
// +psi_bg*pk3*pwl_value_1d(N,lnk,t12_Pkr,log(l3));

// r = (psi_b2g+19./7.*psi_bg)*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
// +9./7.*psi_bg*pwl_value_1d(N,lnk,t2_Pkr,log(l3))
// -11./7.*pow(l3,2.)*psi_bg*pwl_value_1d(N,lnk,t3_Pkr,log(l3))
// +1./7.*psi_bg*pow(l3,4.)*pwl_value_1d(N,lnk,t4_Pkr,log(l3))
// +(2.*psi_b2g+24./7.*psi_bg)*pk3*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
// +2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
// +2./7.*pow(l3,2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
// -4./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
// -4./7.*psi_bg*pk3*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
// +2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t10_Pkr,log(l3));




if (isnan(r) || isinf(r)){
  printf("nan in bispectrum TTG ffts\n");
  printf("pk3 = %.3e\n",pk3);
  printf("psi_bg = %.3e\n",psi_bg);
  printf("psi_b2g = %.3e\n",psi_b2g);
  printf("t1_Pkr = %.3e\n",pwl_value_1d(N,lnk,t1_Pkr,log(l3)));
  printf("t2_Pkr = %.3e\n",pwl_value_1d(N,lnk,t2_Pkr,log(l3)));
  printf("t3_Pkr = %.3e\n",pwl_value_1d(N,lnk,t3_Pkr,log(l3)));
  printf("t4_Pkr = %.3e\n",pwl_value_1d(N,lnk,t4_Pkr,log(l3)));
  printf("t5_Pkr = %.3e\n",pwl_value_1d(N,lnk,t5_Pkr,log(l3)));
  printf("t6_Pkr = %.3e\n",pwl_value_1d(N,lnk,t6_Pkr,log(l3)));
  printf("t7_Pkr = %.3e\n",pwl_value_1d(N,lnk,t7_Pkr,log(l3)));
  printf("t8_Pkr = %.3e\n",pwl_value_1d(N,lnk,t8_Pkr,log(l3)));
  printf("t9_Pkr = %.3e\n",pwl_value_1d(N,lnk,t9_Pkr,log(l3)));
  printf("t10_Pkr = %.3e\n",pwl_value_1d(N,lnk,t10_Pkr,log(l3)));
  printf("t11_Pkr = %.3e\n",pwl_value_1d(N,lnk,t11_Pkr,log(l3)));
  printf("t12_Pkr = %.3e\n",pwl_value_1d(N,lnk,t12_Pkr,log(l3)));
  exit(0);
}

  }

  else if (
    ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gal_gal_lens_3h_fft)
  // ||((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_lens_3h_fft)

  ){
  double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12;

  double z = pvectsz[pclass_sz->index_z];


  int index_l_3 = (int) pvectsz[pclass_sz->index_multipole];
  double l3 = pclass_sz->ell[index_l_3];
  // pvectsz[pclass_sz->index_multipole_for_pk] = l3;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk3 = pvectsz[pclass_sz->index_pk_for_halo_bias];
  double pk3 = get_pk_lin_at_k_and_z((l3+0.5)/chi,z,pba,ppm,pnl,pclass_sz);
  double psi_bg = get_psi_b1kg_at_k_and_z((l3+0.5)/chi,z,pclass_sz);
  double psi_b2g =get_psi_b2kg_at_k_and_z((l3+0.5)/chi,z,pclass_sz);
  // double psi_b2t = get_psi_b2t_at_k_and_z(l3,z,pclass_sz);

  // printf("%.5e %.5e\n",r_m_11,r_m_21);



//// r_m_12 and r_m_22
//// no mass integral
//// apply convolution theorem


  /// set-up:

  double l_min = pclass_sz->l_min_samp_fftw;
  double l_max = pclass_sz->l_max_samp_fftw; // this is a precision parameter
// tabulate the integrand in the "l" dimension:
const int N = pclass_sz->N_samp_fftw;
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
double l;
double pkl=0.;

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
// pvectsz[pclass_sz->index_multipole_for_pk] = l;
// pvectsz[pclass_sz->index_pk_for_halo_bias] = 0.;
// evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
// pkl = pvectsz[pclass_sz->index_pk_for_halo_bias];
pkl = get_pk_lin_at_k_and_z((l+0.5)/chi,z,pba,ppm,pnl,pclass_sz);
fl = get_ksz_filter_at_l(l,pclass_sz);
// if ((l+0.5)/chi>1e-2) fl = 0.;
psi_bt = get_psi_b1g_at_k_and_z((l+0.5)/chi,z,pclass_sz);

psi_b2t = get_psi_b2g_at_k_and_z((l+0.5)/chi,z,pclass_sz);
// l = 1.;

pk_phi_0[ik] = fl*psi_bt;
pk_phi_m2[ik] = pow((l+0.5)/chi,-2)*fl*psi_bt;
pk_phi_4[ik] = pow((l+0.5)/chi,4)*fl*psi_bt;
pk_phi_2[ik] = pow((l+0.5)/chi,2)*fl*psi_bt;

pk_tilde_phi_0[ik] = fl*pkl*psi_bt;
pk_tilde_phi_m2[ik] = pow((l+0.5)/chi,-2)*fl*pkl*psi_bt;
pk_tilde_phi_2[ik] = pow((l+0.5)/chi,2)*fl*pkl*psi_bt;
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
xi2pk(N,k,pk_phi_0,rp,xi_phi_0,pclass_sz);
xi2pk(N,k,pk_phi_2,rp,xi_phi_2,pclass_sz);
xi2pk(N,k,pk_phi_m2,rp,xi_phi_m2,pclass_sz);
xi2pk(N,k,pk_tilde_phi_0,rp,xi_tilde_phi_0,pclass_sz);
xi2pk(N,k,pk_tilde_phi_b20,rp,xi_tilde_phi_b20,pclass_sz);
xi2pk(N,k,pk_tilde_phi_2,rp,xi_tilde_phi_2,pclass_sz);
xi2pk(N,k,pk_tilde_phi_m2,rp,xi_tilde_phi_m2,pclass_sz);
xi2pk(N,k,pk_phi_4,rp,xi_phi_4,pclass_sz);

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
pk2xi(N,rp,t1_xi12,k,t1_Pkr,pclass_sz);
pk2xi(N,rp,t2_xi12,k,t2_Pkr,pclass_sz);
pk2xi(N,rp,t3_xi12,k,t3_Pkr,pclass_sz);
pk2xi(N,rp,t4_xi12,k,t4_Pkr,pclass_sz);
pk2xi(N,rp,t5_xi12,k,t5_Pkr,pclass_sz);
pk2xi(N,rp,t6_xi12,k,t6_Pkr,pclass_sz);
pk2xi(N,rp,t7_xi12,k,t7_Pkr,pclass_sz);
pk2xi(N,rp,t8_xi12,k,t8_Pkr,pclass_sz);
pk2xi(N,rp,t9_xi12,k,t9_Pkr,pclass_sz);
pk2xi(N,rp,t10_xi12,k,t10_Pkr,pclass_sz);
pk2xi(N,rp,t11_xi12,k,t11_Pkr,pclass_sz);
pk2xi(N,rp,t12_xi12,k,t12_Pkr,pclass_sz);


r = 10./14.*psi_bg*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
   -5./7.*psi_bg*pwl_value_1d(N,lnk,t2_Pkr,log(l3))
   +3./7.*psi_bg*pow((l3+0.5)/chi,2.)*pwl_value_1d(N,lnk,t3_Pkr,log(l3))
   +1./7.*psi_bg*pow((l3+0.5)/chi,4.)*pwl_value_1d(N,lnk,t4_Pkr,log(l3))
   // b2 terms:
   +psi_b2g*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
   +2.*psi_bg*pk3*pwl_value_1d(N,lnk,t12_Pkr,log(l3))

   +10./14.*pk3*psi_bg*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
   +3./14.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
   +3./14.*pk3*psi_bg*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
   -5./14.*pk3*psi_bg*pow((l3+0.5)/chi,2.)*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
   -5./14.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
   +1./7.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t10_Pkr,log(l3))

   +10./14.*pk3*psi_bg*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
   +3./14.*pk3*psi_bg*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
   +3./14.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
   -5./14.*pk3*psi_bg*pow((l3+0.5)/chi,2.)*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
   -5./14.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
   +1./7.*pk3*psi_bg*pow((l3+0.5)/chi,-2.)*pwl_value_1d(N,lnk,t10_Pkr,log(l3));



   //+2./7.*pwl_value_1d(N,lnk,t2_Pkr,log(l3));
// r = 19./7.*psi_bg*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
//     +9./7.*psi_bg*pwl_value_1d(N,lnk,t2_Pkr,log(l3))
//     -11./7.*pow((l3+0.5)/chi,2.)*psi_bg*pwl_value_1d(N,lnk,t3_Pkr,log(l3))
//     +1./7.*psi_bg*pow((l3+0.5)/chi,4.)*pwl_value_1d(N,lnk,t4_Pkr,log(l3));
// +24./7.*psi_bg*pk3*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
// +2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
// +2./7.*pow(l3,2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
// -4./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
// -4./7.*psi_bg*pk3*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
// +2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t10_Pkr,log(l3))
// // b2 terms:
 // psi_b2g*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
// +psi_bg*pwl_value_1d(N,lnk,t11_Pkr,log(l3))
// +psi_bg*pk3*pwl_value_1d(N,lnk,t12_Pkr,log(l3));

// r = (psi_b2g+19./7.*psi_bg)*pwl_value_1d(N,lnk,t1_Pkr,log(l3))
// +9./7.*psi_bg*pwl_value_1d(N,lnk,t2_Pkr,log(l3))
// -11./7.*pow(l3,2.)*psi_bg*pwl_value_1d(N,lnk,t3_Pkr,log(l3))
// +1./7.*psi_bg*pow(l3,4.)*pwl_value_1d(N,lnk,t4_Pkr,log(l3))
// +(2.*psi_b2g+24./7.*psi_bg)*pk3*pwl_value_1d(N,lnk,t5_Pkr,log(l3))
// +2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t6_Pkr,log(l3))
// +2./7.*pow(l3,2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t7_Pkr,log(l3))
// -4./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t8_Pkr,log(l3))
// -4./7.*psi_bg*pk3*pwl_value_1d(N,lnk,t9_Pkr,log(l3))
// +2./7.*pow(l3,-2.)*psi_bg*pk3*pwl_value_1d(N,lnk,t10_Pkr,log(l3));




if (isnan(r) || isinf(r)){
  printf("nan in bispectrum GGK ffts\n");
  printf("pk3 = %.3e\n",pk3);
  printf("psi_bg = %.3e\n",psi_bg);
  printf("psi_b2g = %.3e\n",psi_b2g);
  printf("t1_Pkr = %.3e\n",pwl_value_1d(N,lnk,t1_Pkr,log(l3)));
  printf("t2_Pkr = %.3e\n",pwl_value_1d(N,lnk,t2_Pkr,log(l3)));
  printf("t3_Pkr = %.3e\n",pwl_value_1d(N,lnk,t3_Pkr,log(l3)));
  printf("t4_Pkr = %.3e\n",pwl_value_1d(N,lnk,t4_Pkr,log(l3)));
  printf("t5_Pkr = %.3e\n",pwl_value_1d(N,lnk,t5_Pkr,log(l3)));
  printf("t6_Pkr = %.3e\n",pwl_value_1d(N,lnk,t6_Pkr,log(l3)));
  printf("t7_Pkr = %.3e\n",pwl_value_1d(N,lnk,t7_Pkr,log(l3)));
  printf("t8_Pkr = %.3e\n",pwl_value_1d(N,lnk,t8_Pkr,log(l3)));
  printf("t9_Pkr = %.3e\n",pwl_value_1d(N,lnk,t9_Pkr,log(l3)));
  printf("t10_Pkr = %.3e\n",pwl_value_1d(N,lnk,t10_Pkr,log(l3)));
  printf("t11_Pkr = %.3e\n",pwl_value_1d(N,lnk,t11_Pkr,log(l3)));
  printf("t12_Pkr = %.3e\n",pwl_value_1d(N,lnk,t12_Pkr,log(l3)));
  exit(0);
}

  }
  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_gal_3h){
  double r_m_b1t1;
  double r_m_b1t2;
  double r_m_b2t1;
  double r_m_b2t2;
  double r_m_b1g3;
  double r_m_b2g3;
  double r_tab;

  int index_theta_1 = (int) pvectsz[pclass_sz->index_multipole_1];
  double theta_1 = pclass_sz->theta_kSZ2_gal_theta_grid[index_theta_1];
  // double cos_theta_1 = pclass_sz->theta_kSZ2_gal_theta_grid[index_theta_1];
  int index_l_2 = (int) pvectsz[pclass_sz->index_multipole_2];
  int index_l_3 = (int) pvectsz[pclass_sz->index_multipole_3];
  double l2 = exp(pclass_sz->ell_kSZ2_gal_multipole_grid[index_l_2]);
  double l3 = pclass_sz->ell[index_l_3];
  double ell = l3;
  double ell_prime = l2;
  double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos(theta_1));
  // double l1 = sqrt(ell*ell+ell_prime*ell_prime+2.*ell*ell_prime*cos_theta_1);
  double k1 = (l1 + 0.5)/chi;
  double k2 = (l2 + 0.5)/chi;
  double k3 = (l3 + 0.5)/chi;
  double pk1 = 0.;
  double pk2 = 0.;
  double pk3 = 0.;

  double z = pvectsz[pclass_sz->index_z];

  // // r_m_11*r_m_21
  // pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  // V.pvectsz = pvectsz;
  // params = &V;
  // r_m_b1t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
  //                                          epsrel, epsabs,
  //                                          integrand_mass,
  //                                          params,pclass_sz->patterson_show_neval);
  //
  //  if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
  //    double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
  //    double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
  //    double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
  //    double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
  //    r_m_b1t1 += bmin_umin;
  //    // printf("counter terms done r_m_1\n");
  // }
  r_tab = get_psi_b1t_at_k_and_z(k1,z,pclass_sz);
  // printf("r_m_b1t1 %.8e %.8e\n",r_m_b1t1,r_tab);
  r_m_b1t1 = r_tab;



  // pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  // V.pvectsz = pvectsz;
  // params = &V;
  // r_m_b1t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
  //                                          epsrel, epsabs,
  //                                          integrand_mass,
  //                                          params,pclass_sz->patterson_show_neval);
  //
  //  if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
  //    double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
  //    double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
  //    double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
  //    double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
  //    r_m_b1t2 += bmin_umin;
  //    // printf("counter terms done r_m_1\n");
  // }
  r_tab = get_psi_b1t_at_k_and_z(k2,z,pclass_sz);
  // printf("r_m_b1t2 %.8e %.8e\n",r_m_b1t2,r_tab);
  r_m_b1t2  = r_tab;




  // // r_m_12*r_m_22
  // pvectsz[pclass_sz->index_part_id_cov_hsv] = 3;
  // V.pvectsz = pvectsz;
  // params = &V;
  // r_m_b1g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
  //                                          epsrel, epsabs,
  //                                          integrand_mass,
  //                                          params,pclass_sz->patterson_show_neval);
  //
  //  if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
  //    double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
  //    double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
  //    double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
  //    double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
  //    r_m_b1g3 += bmin_umin;
  //    // printf("counter terms done r_m_1\n");
  // }
  r_tab = get_psi_b1g_at_k_and_z(k3,z,pclass_sz);
  // printf("r_m_b1g3 %.8e %.8e\n",r_m_b1g3,r_tab);
  r_m_b1g3  = r_tab;




  // pvectsz[pclass_sz->index_part_id_cov_hsv] = 4;
  // V.pvectsz = pvectsz;
  // params = &V;
  // r_m_b2g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
  //                                          epsrel, epsabs,
  //                                          integrand_mass,
  //                                          params,pclass_sz->patterson_show_neval);
  //
  //  if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
  //    double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
  //    double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
  //    double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
  //    double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
  //    r_m_b2g3 += bmin_umin;
  //    // printf("counter terms done r_m_1\n");
  // }
  r_tab = get_psi_b2g_at_k_and_z(k3,z,pclass_sz);
  r_m_b2g3 = r_tab;

  r_tab = get_psi_b2t_at_k_and_z(k1,z,pclass_sz);
  r_m_b2t1 = r_tab;

  r_tab = get_psi_b2t_at_k_and_z(k2,z,pclass_sz);
  r_m_b2t2 = r_tab;
  //
  // // int index_l_1 = (int) pvectsz[pclass_sz->index_multipole_1];
  //
  // pvectsz[pclass_sz->index_multipole_for_pk] = l1;//pclass_sz->ell_kSZ2_gal_multipole_grid[index_l_1];
  // pvectsz[pclass_sz->index_pk_for_halo_bias] = 0.;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk1 = pvectsz[pclass_sz->index_pk_for_halo_bias];
  //
  // // int index_l_2 = (int) pvectsz[pclass_sz->index_multipole_2];
  // pvectsz[pclass_sz->index_multipole_for_pk] = l2;//pclass_sz->ell_kSZ2_gal_multipole_grid[index_l_2];
  // pvectsz[pclass_sz->index_pk_for_halo_bias] = 0.;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk2 = pvectsz[pclass_sz->index_pk_for_halo_bias];
  //
  // // int index_l_3 = (int) pvectsz[pclass_sz->index_multipole_3];
  // pvectsz[pclass_sz->index_multipole_for_pk] = l3;//pclass_sz->ell[index_l_3];
  // pvectsz[pclass_sz->index_pk_for_halo_bias] = 0.;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk3 = pvectsz[pclass_sz->index_pk_for_halo_bias];

  pk1 = get_pk_lin_at_k_and_z(k1,z,pba,ppm,pnl,pclass_sz);
  pk2 = get_pk_lin_at_k_and_z(k2,z,pba,ppm,pnl,pclass_sz);
  pk3 = get_pk_lin_at_k_and_z(k3,z,pba,ppm,pnl,pclass_sz);



  // double d_A = pvecback[pba->index_bg_ang_distance]*pba->h*(1.+z);

  // double Fk1k2 = bispectrum_f2_kernel(k1,k2,k3);
  // double Fk1k3 = bispectrum_f2_kernel(k3,k1,k2);
  // double Fk2k3 = bispectrum_f2_kernel(k2,k3,k1);

  double f2_123 = bispectrum_f2_kernel(k1,k2,k3);
  double f2_312 = bispectrum_f2_kernel(k3,k1,k2);
  double f2_231 = bispectrum_f2_kernel(k2,k3,k1);
  // printf("f2_123 = %.8e\n",f2_123);

  // double comb_pks = pk1*pk2+pk1*pk3+pk2*pk3;
  // double comb_pks_fks = 2.*pk1*pk2*Fk1k2+2.*pk1*pk3*Fk1k3+2.*pk2*pk3*Fk2k3;


  // r = r_m_b1t1*r_m_b1t2*r_m_b1g3*comb_pks_fks+r_m_b1t1*r_m_b1t2*r_m_b2g3*comb_pks;
  // r_m_b1t1 = 1.;
  // r_m_b1t2 = 1.;
  // r_m_b1g3 = 1.;
  // r_m_b2g3 = 1.;
  // r_m_b2t1 = 1.;
  // r_m_b2t2 = 1.;
  //
  //
  // pk1 = 1.;
  // pk2 = 1.;
  // pk3 = 1.;
  // f2_123 = 1.;

  r =2.*r_m_b1t1*r_m_b1t2*r_m_b1g3*f2_123*pk1*pk2
    +2.*r_m_b1t1*r_m_b1t2*r_m_b1g3*f2_312*pk3*pk1
    +2.*r_m_b1t1*r_m_b1t2*r_m_b1g3*f2_231*pk2*pk3
    +r_m_b1t1*r_m_b1t2*r_m_b2g3*pk1*pk2
    +r_m_b1t1*r_m_b2t2*r_m_b1g3*pk3*pk1
    +r_m_b2t1*r_m_b1t2*r_m_b1g3*pk2*pk3;

// printf("r = %.8e\n",r);


  }

  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_bk_at_z_2h){
  double r_m_b1t1;
  double r_m_b1t2;
  double r_m_b1t1g3;
  double r_m_b1t2g3;
  double r_m_b1g3;
  double r_m_b1t1t2;

  // r_m_11*r_m_21
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }



  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t2g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t2g3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1g3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t1t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t1t2+= bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 5;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t2+= bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 6;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t1g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t1g3+= bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  // r = r_m_b1t1*r_m_b1t2;
  int index_k = (int) pvectsz[pclass_sz->index_k_for_pk_hm];
  double k = pclass_sz->k_for_pk_hm[index_k];
  double pk1, pk2, pk3;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk = pvectsz[pclass_sz->index_pk_for_halo_bias];
  // double pk = get_pk_lin_at_k_and_z(k,z,pba,ppm,pnl,pclass_sz);

  pk1 = get_pk_lin_at_k_and_z(k,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  pk2 = get_pk_lin_at_k_and_z(pclass_sz->bispectrum_lambda_k2*k,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  pk3 = get_pk_lin_at_k_and_z(pclass_sz->bispectrum_lambda_k3*k,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);

  r = pk3*r_m_b1g3*r_m_b1t1t2
     +pk2*r_m_b1t1g3*r_m_b1t2
     +pk1*r_m_b1t1*r_m_b1t2g3;

  }

  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_bk_at_z_3h){
  double r_m_b1t1;
  double r_m_b2t1;
  double r_m_b1t2;
  double r_m_b2t2;
  double r_m_b1g3;
  double r_m_b2g3;

  // r_m_11*r_m_21
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b2t1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }



  pvectsz[pclass_sz->index_part_id_cov_hsv] = 6;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t2 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 7;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b2t2 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 8;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1g3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 9;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2g3 = Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b2g3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }


  int index_k = (int) pvectsz[pclass_sz->index_k_for_pk_hm];
  double k = pclass_sz->k_for_pk_hm[index_k];
if (pclass_sz->check_consistency_conditions == 1){
  // check consistency conditions:

  double r_m_mean;
  double r_mass;
  double r_b1;
  double r_b2;
  // mass consistency
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_mean=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double nmin_umin = nmin*I0/pvectsz[pclass_sz->index_hmf];
     r_m_mean += nmin_umin;
  }
  r_mass = r_m_mean;


  // b1 consistency
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_mean=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_mean += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }
  r_b1 = r_m_mean;


  // b2 consistency
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 5;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_mean=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_mean += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }
  r_b2 = r_m_mean;
  printf("hm consistency z = %.3e k = %.8e m = %.8e b1 = %.8e b2 %.8e\n",pvectsz[pclass_sz->index_z],k,r_mass,r_b1,1.-r_b2);
}

  // double bh;
  // double pk;
  // double b0;
  // double f2;
  //
  // // printf("result 3h = %.3e\n",result);
  // // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // // pk = pvectsz[pclass_sz->index_pk_for_halo_bias];
  // pk = get_pk_lin_at_k_and_z(k,z,pba,ppm,pnl,pclass_sz);
  // f2 = bispectrum_f2_kernel(k,k,k);
  // b0 = (2.*pk*pk*f2)*r_m_b1t1*r_m_b1t1*r_m_b1t1;
  // bh = 0.;//3.*pk*pk*r_m_b1t1*r_m_b1t1*r_m_b2t1;
  // r  = (bh+b0);

// r_m_b2t1 = 0.;
  // r = 3.*(2.*r_m_b1t1*r_m_b1t1*r_m_b1t1*f2+r_m_b1t1*r_m_b1t1*r_m_b2t1)*pk*pk;
  // double b_tree = get_matter_bispectrum_at_z_tree_level_PT(k,
  //                                                          pclass_sz->bispectrum_lambda_k2,
  //                                                          pclass_sz->bispectrum_lambda_k3,
  //                                                          pvectsz[pclass_sz->index_z],
  //                                                          pclass_sz,pba,pnl,ppm);
  // printf("bispectrum fields z = %.3e k = %.8e <bu> = %.8e <b2u> = %.8e b_hm = %.8e b_tree = %.8e\n",pvectsz[pclass_sz->index_z],k,r_m_b1t1,r_m_b2t1,r,b_tree);
  double k1,k2,k3;
  k1 = k;
  k2 = pclass_sz->bispectrum_lambda_k2*k;
  k3 = pclass_sz->bispectrum_lambda_k3*k;
  double pk1, pk2, pk3;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk = pvectsz[pclass_sz->index_pk_for_halo_bias];
  // double pk = get_pk_lin_at_k_and_z(k,z,pba,ppm,pnl,pclass_sz);
  pk1 = get_pk_lin_at_k_and_z(k1,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  pk2 = get_pk_lin_at_k_and_z(k2,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  pk3 = get_pk_lin_at_k_and_z(k3,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);

  double f2_123 = bispectrum_f2_kernel(k1,k2,k3);
  double f2_231 = bispectrum_f2_kernel(k2,k3,k1);
  double f2_312 = bispectrum_f2_kernel(k3,k1,k2);

  // r_m_b1g3 = 1.;
  //
  // r_m_b2g3 = 0.;
  // r_m_b2t2 = 0.;
  // r_m_b2t1 = 0.;



  r = 2.*r_m_b1t1*r_m_b1t2*r_m_b1g3*f2_123*pk1*pk2
     +2.*r_m_b1t1*r_m_b1t2*r_m_b1g3*f2_312*pk3*pk1
     +2.*r_m_b1t1*r_m_b1t2*r_m_b1g3*f2_231*pk2*pk3
     +r_m_b1t1*r_m_b1t2*r_m_b2g3*pk1*pk2
     +r_m_b1t1*r_m_b2t2*r_m_b1g3*pk3*pk1
     +r_m_b2t1*r_m_b1t2*r_m_b1g3*pk2*pk3;


  }

  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_bk_ttg_at_z_2h){
  double r_m_b1t1;
  double r_m_b1t2;
  double r_m_b1t1g3;
  double r_m_b1t2g3;
  double r_m_b1g3;
  double r_m_b1t1t2;
  // r_m_11*r_m_21
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t1 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }



  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t2 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t1g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t1g3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t2g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t2g3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 5;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1g3 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 6;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t1t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t1t2 += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

  int index_k = (int) pvectsz[pclass_sz->index_k_for_pk_hm];
  double k = pclass_sz->k_for_pk_hm[index_k];
  double pk1, pk2, pk3;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk = pvectsz[pclass_sz->index_pk_for_halo_bias];
  // double pk = get_pk_lin_at_k_and_z(k,z,pba,ppm,pnl,pclass_sz);

  pk1 = get_pk_lin_at_k_and_z(k,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  pk2 = get_pk_lin_at_k_and_z(pclass_sz->bispectrum_lambda_k2*k,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);
  pk3 = get_pk_lin_at_k_and_z(pclass_sz->bispectrum_lambda_k3*k,pvectsz[pclass_sz->index_z],pba,ppm,pnl,pclass_sz);

  r = pk3*r_m_b1g3*r_m_b1t1t2
     +pk2*r_m_b1t1g3*r_m_b1t2
     +pk1*r_m_b1t1*r_m_b1t2g3;

  }

  else if ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_bk_ttg_at_z_3h){
  double r_m_b1t1;
  double r_m_b2t1;
  double r_m_b1t2;
  double r_m_b2t2;
  double r_m_b1g3;
  double r_m_b2g3;

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t1 += bmin_umin;

  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1t2 += bmin_umin;

  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 3;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b1g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_b1g3 += bmin_umin;

  }

  pvectsz[pclass_sz->index_part_id_cov_hsv] = 4;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2t1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b2t1 += bmin_umin;

  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 5;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2t2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b2t2 += bmin_umin;

  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 6;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_b2g3=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r_m_b2g3 += bmin_umin;

  }



  int index_k = (int) pvectsz[pclass_sz->index_k_for_pk_hm];
  double z = pvectsz[pclass_sz->index_z];
  double k = pclass_sz->k_for_pk_hm[index_k];
  double k1,k2,k3;
  k1 = k;
  k2 = pclass_sz->bispectrum_lambda_k2*k;
  k3 = pclass_sz->bispectrum_lambda_k3*k;
  double pk1, pk2, pk3;
  // evaluate_pk_at_ell_plus_one_half_over_chi(pvecback,pvectsz,pba,ppm,pnl,pclass_sz);
  // double pk = pvectsz[pclass_sz->index_pk_for_halo_bias];
  // double pk = get_pk_lin_at_k_and_z(k,z,pba,ppm,pnl,pclass_sz);
  pk1 = get_pk_lin_at_k_and_z(k1,z,pba,ppm,pnl,pclass_sz);
  pk2 = get_pk_lin_at_k_and_z(k2,z,pba,ppm,pnl,pclass_sz);
  pk3 = get_pk_lin_at_k_and_z(k3,z,pba,ppm,pnl,pclass_sz);

  double f2_123 = bispectrum_f2_kernel(k1,k2,k3);
  double f2_231 = bispectrum_f2_kernel(k2,k3,k1);
  double f2_312 = bispectrum_f2_kernel(k3,k1,k2);

  // r_m_b1g3 = 1.;
  //
  if (pclass_sz->no_b2){
  r_m_b2g3 = 0.;
  r_m_b2t2 = 0.;
  r_m_b2t1 = 0.;
}



  r = 2.*r_m_b1t1*r_m_b1t2*r_m_b1g3*f2_123*pk1*pk2
     +2.*r_m_b1t1*r_m_b1t2*r_m_b1g3*f2_312*pk3*pk1
     +2.*r_m_b1t1*r_m_b1t2*r_m_b1g3*f2_231*pk2*pk3
     +r_m_b1t1*r_m_b1t2*r_m_b2g3*pk1*pk2
     +r_m_b1t1*r_m_b2t2*r_m_b1g3*pk3*pk1
     +r_m_b2t1*r_m_b1t2*r_m_b1g3*pk2*pk3;


// double z = pvectsz[pclass_sz->index_z];
// double r_effective = get_ttg_bispectrum_at_z_tree_level_PT(k,k,k,z,pclass_sz,pba,pnl,ppm);
// printf("bispectrum z = %.3e k = %.8e r_m_b1g3 %.8e b_hm = %.8e b_tree = %.8e\n",z,k,r_m_b1g3,r,r_effective);

  }

  else if((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_pk_em_at_z_2h){

  double r_m_1;
  pvectsz[pclass_sz->index_part_id_cov_hsv] = 1;
  V.pvectsz = pvectsz;
  params = &V;
  r_m_1=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_1 += bmin_umin;

  }


  pvectsz[pclass_sz->index_part_id_cov_hsv] = 2;
  V.pvectsz = pvectsz;
  params = &V;
  double r_m_2;
  r_m_2=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                           epsrel, epsabs,
                                           integrand_mass,
                                           params,pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r_m_2 += bmin_umin;

  }

r = r_m_1*r_m_2;
  }

// dont apply halo model consistency
else if(((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_szrates)
|| ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_hmf)
){
  // if (pclass_sz->sz_verbose>0) printf("starting mass integral for szrate id = %d.\n", (int)pvectsz[pclass_sz->index_szrate]);
  // printf("integrating over mass m_min = %.3e m_max = %.3e\n",m_min,m_max);
  r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                       epsrel, epsabs,
                                       integrand_mass,
                                       params,pclass_sz->patterson_show_neval);


// printf("found r = %.5e\n",r);
}

else {
// here we treat all the 1-halo terms and also the 2halo terms that are auto
// printf("integrating over mass\n");
  r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                       epsrel, epsabs,
                                       integrand_mass,
                                       params,pclass_sz->patterson_show_neval);
// printf("got r = %.5e\n",r);
        // halo model consistency:

        if ( (int) pvectsz[pclass_sz->index_md] != pclass_sz->index_md_cov_N_N
          && (int) pvectsz[pclass_sz->index_md] != pclass_sz->index_md_cov_N_N_hsv){

          if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
              // autocorrelation 2-halo cases (correlation of same fields).
              if (( (int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_2halo)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_m_y_y_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_pk_at_z_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_pk_gg_at_z_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_pk_bb_at_z_2h)
                    | ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_pk_b_at_z_2h)
                    // || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_pk_em_at_z_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_pk_HI_at_z_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_lens_lens_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_tau_tau_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_custom1_custom1_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_lensmag_lensmag_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_lens_lensmag_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gallens_lensmag_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_nlensmag_gallens_2h)
                    || (((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_cib_cib_2h)  && (pvectsz[pclass_sz->index_frequency_for_cib_profile] == pvectsz[pclass_sz->index_frequency_prime_for_cib_profile]) )
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gal_gal_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gallens_gallens_2h)
                    || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gallens_lens_2h)
                    ){

                  double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
                  double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
                  double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
                  double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
                  r += bmin_umin;
                  if (isnan(bmin_umin)){
                    printf("nan in concistency condition bmin = %.5e n = %.5e I0 = %.5e b = %.5e.\n",
                  bmin,pvectsz[pclass_sz->index_hmf],I0,pvectsz[pclass_sz->index_halo_bias]);
                    exit(0);
                  }
                  }
              // all of the 1-halo cases
              else {

                if (pclass_sz->sz_verbose>10)
                  printf("adding counter terms 1h\n");


                  // treat some cases seprately (e.g., useful when comparing with sims)
                  if (((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gal_lens_1h)){
                    if (pclass_sz->include_gk_counterterms_in_gk){
                        double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
                        double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
                        double nmin_umin = nmin*I0/pvectsz[pclass_sz->index_hmf];
                        r += nmin_umin;
                        }
                    }

                  else{
                    double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
                    double I0 = integrand_mass(log(pclass_sz->m_min_counter_terms),params);
                    double nmin_umin = nmin*I0/pvectsz[pclass_sz->index_hmf];
                    r += nmin_umin;
                    }
              }
                  }

                                           }
  }


// for autocorelations 2-halo we square the results
if (( (int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_2halo)
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_m_y_y_2h)
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_kSZ_kSZ_2h)
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_pk_at_z_2h)
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_pk_gg_at_z_2h)
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_pk_bb_at_z_2h)
 // || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_pk_em_at_z_2h)
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_pk_HI_at_z_2h)
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_lens_lens_2h)
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_tau_tau_2h)
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_custom1_custom1_2h)
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_lensmag_lensmag_2h)
 // || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_lens_lensmag_2h)
// || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gallens_lensmag_2h)
 || (((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_cib_cib_2h)
      && (pvectsz[pclass_sz->index_frequency_for_cib_profile] == pvectsz[pclass_sz->index_frequency_prime_for_cib_profile]) )
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gal_gal_2h)
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gallens_gallens_2h)
 || ((int) pvectsz[pclass_sz->index_md] == pclass_sz->index_md_gallens_lens_2h)
){
  // printf("squaring result\n");
 pvectsz[pclass_sz->index_integral_over_m] = r*r;
 }
else
pvectsz[pclass_sz->index_integral_over_m] = r;

// printf("pvectsz[pclass_sz->index_integral_over_m] = %.3e\n",pvectsz[pclass_sz->index_integral_over_m]);

//}

return pvectsz[pclass_sz->index_integral_over_m];

}




//This routine reads the tabulated
//Planck noise map

int read_Planck_noise_map(struct class_sz_structure * pclass_sz)
{
  ///read theta file for completeness
  /////////////////////////////start read theta file
  if (pclass_sz->sz_verbose >= 3){
    printf("Loading theta file\n");
    printf("Planck thetas file: %s\n", pclass_sz->Planck_thetas_file);

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
  pclass_sz->thetas   = (double *)malloc(n_data_guess*sizeof(double));

  char Filepath[_ARGUMENT_LENGTH_MAX_];
  // //printf("%s\n",pclass_sz->path_to_class);
  // sprintf(Filepath,
  //         // "%s%s%s",
  //         "%s%s",
  //         "cat ",
  //         // pclass_sz->path_to_class,
  //         "/class_sz_auxiliary_files/SZ_thetas.txt");
  //
  // process = popen(Filepath, "r");
  if (pclass_sz->experiment == 0){
  // class_open(process,"class_sz_auxiliary_files/SZ_thetas.txt", "r",pclass_sz->error_message);
  class_open(process,pclass_sz->Planck_thetas_file, "r",pclass_sz->error_message);

}
  else if (pclass_sz->experiment == 1){
  // class_open(process,"class_sz_auxiliary_files/so_3freqs_020621_thetas.txt", "r",pclass_sz->error_message);
  class_open(process,pclass_sz->SO_thetas_file, "r",pclass_sz->error_message);

}

  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

  if (pclass_sz->sz_verbose >= 3)
    printf("%lf\n", this_lnx);

    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(pclass_sz->thetas,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
                 "Error allocating memory to read the thetas.\n");
      pclass_sz->thetas = tmp;
    };


    /* Store */
    pclass_sz->thetas[n_data]   = this_lnx;
    n_data++;
  }

  // status = pclose(process);
  status = fclose(process);
    // printf("Loading theta file 2\n");
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");
    // printf("Loading theta file 3\n");

  pclass_sz->nthetas = n_data;

  pclass_sz->theta_bin_min = pclass_sz->thetas[0];
  int c;
  for (c = 1; c < pclass_sz->nthetas; c++)
  {
    if (pclass_sz->thetas[c] < pclass_sz->theta_bin_min)
    {
      pclass_sz->theta_bin_min = pclass_sz->thetas[c];
      //location = c+1;
    }
  }
  pclass_sz->theta_bin_max = pclass_sz->thetas[0];
  for (c = 1; c < pclass_sz->nthetas; c++)
  {
    if (pclass_sz->thetas[c] > pclass_sz->theta_bin_max)
    {
      pclass_sz->theta_bin_max = pclass_sz->thetas[c];
      //location = c+1;
    }
  }
  // printf("theta_bin_max:=%e\n",pclass_sz->theta_bin_max);
  // printf("theta_bin_min:=%e\n",pclass_sz->theta_bin_min);

  pclass_sz->Nth = pclass_sz->nthetas;
  class_alloc(pclass_sz->erfs_2d_to_1d_th_array,pclass_sz->Nth*sizeof(double *),pclass_sz->error_message);


  for (c = 0; c < pclass_sz->Nth; c++){
    pclass_sz->erfs_2d_to_1d_th_array[c] = log(pclass_sz->thetas[c]);

  }

  ///////////////////////////end read theta file

  //end read theta file for completeness
  // start read noise map for completeness
  //read skyfracs

  //double *skyfracs = NULL;

  if (pclass_sz->sz_verbose >= 3){
    printf("theta file loaded with ntheta = %d\n",pclass_sz->nthetas);
  }
  // printf("Loading theta file 2\n");

  n_data = 0;
  n_data_guess = 100;
  pclass_sz->skyfracs   = (double *)malloc(n_data_guess*sizeof(double));

  // sprintf(Filepath,
  //         // "%s%s%s",
  //         "%s%s",
  //         "cat ",
  //         // pclass_sz->path_to_class,
  //         "/class_sz_auxiliary_files/SZ_skyfracs.txt");
  //
  // process = popen(Filepath, "r");
  if (pclass_sz->sz_verbose >= 3){
    printf("Loading skyfrac file\n");
  }
  //class_open(process,"class_sz_auxiliary_files/SZ_skyfracs.txt", "r",pclass_sz->error_message);

  if (pclass_sz->experiment == 0){
  // class_open(process,"class_sz_auxiliary_files/SZ_skyfracs.txt", "r",pclass_sz->error_message);
  class_open(process,pclass_sz->Planck_skyfracs_file, "r",pclass_sz->error_message);

}
  else if (pclass_sz->experiment == 1){
  // class_open(process,"class_sz_auxiliary_files/so_3freqs_020621_skyfracs.txt", "r",pclass_sz->error_message);
  class_open(process,pclass_sz->SO_skyfracs_file, "r",pclass_sz->error_message);

}



  while (fgets(line, sizeof(line)-1, process) != NULL) {
    sscanf(line, "%lf", &this_lnx);

    if((n_data+1) > n_data_guess) {
      n_data_guess *= 2;
      tmp = (double *)realloc(pclass_sz->skyfracs,   n_data_guess*sizeof(double));
      class_test(tmp == NULL,
                 pclass_sz->error_message,
                 "Error allocating memory to read the thetas.\n");
      pclass_sz->skyfracs = tmp;
    };


    /* Store */
    pclass_sz->skyfracs[n_data]   = this_lnx;
    n_data++;
  }

  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  pclass_sz->nskyfracs = n_data;
  if (pclass_sz->sz_verbose >= 3){
    printf("sky frac file loaded with nskyfracs = %d\n",pclass_sz->nskyfracs);
  }

  //end read skyfracs

  if (pclass_sz->sz_verbose >= 3){
    printf("Loading noise map\n");
  }
  ////////////////////////read the ylims
  int index_patches;
  //double ** ylims = NULL;

  class_alloc(pclass_sz->ylims,
              pclass_sz->nskyfracs*sizeof(double *),
              pclass_sz->error_message);



  for (index_patches=0;
       index_patches<pclass_sz->nskyfracs;
       index_patches++)
  {
    class_alloc(pclass_sz->ylims[index_patches],
                pclass_sz->nthetas*sizeof(double),
                pclass_sz->error_message);
  }



  if (pclass_sz->experiment == 0){
    // class_open(process,"class_sz_auxiliary_files/SZ_ylims.txt", "r",pclass_sz->error_message);
    class_open(process,pclass_sz->Planck_ylims_file, "r",pclass_sz->error_message);

  }
  else if (pclass_sz->experiment == 1){
    // class_open(process,"class_sz_auxiliary_files/so_3freqs_020621_ylims.txt", "r",pclass_sz->error_message);
    class_open(process,pclass_sz->SO_ylims_file, "r",pclass_sz->error_message);

  }

  //printf("ok\n");
  // printf("noise map loaded 0\n");
  int id_patches=0;
  int index_thetas = 0;

  for (index_thetas =0;index_thetas<pclass_sz->nthetas;index_thetas++){
    for (id_patches =0;id_patches<pclass_sz->nskyfracs;id_patches++){
    fgets(line, sizeof(line)-1, process);
    sscanf(line, "%lf", &this_lnx);
    // printf("%.3e id_p = %d id_t = %d n_p = %d n_t = %d\n",this_lnx,index_patches,index_thetas,pclass_sz->nskyfracs,pclass_sz->nthetas);
    pclass_sz->ylims[id_patches][index_thetas]=this_lnx;
    }
  }
  // printf("noise map loaded 1\n");

  // status = pclose(process);
  status = fclose(process);
  class_test(status != 0.,
             pclass_sz->error_message,
             "The attempt to launch the external command was unsuccessful. "
             "Try doing it by hand to check for errors.");

  ///end read the files
  //end reads noise map for completeness
  class_alloc(pclass_sz->sky_averaged_ylims,
              pclass_sz->nthetas*sizeof(double),
              pclass_sz->error_message);

double sum_skyfracs = 0.;
for (index_patches=0;
     index_patches<pclass_sz->nskyfracs;
     index_patches++)
     sum_skyfracs += pclass_sz->skyfracs[index_patches];
     if (pclass_sz->sz_verbose >= 1){
       printf("sum_skyfracs =  %.3e\n",sum_skyfracs);}
      pclass_sz->fsky_from_skyfracs = sum_skyfracs;

for (index_thetas = 0; index_thetas<pclass_sz->nthetas; index_thetas ++){
  pclass_sz->sky_averaged_ylims[index_thetas] = 0.;
  for (index_patches=0;
       index_patches<pclass_sz->nskyfracs;
       index_patches++)
  {
    pclass_sz->sky_averaged_ylims[index_thetas] += pclass_sz->skyfracs[index_patches]*pclass_sz->ylims[index_patches][index_thetas]/sum_skyfracs;
  }
  if (pclass_sz->sz_verbose >= 2){
printf("sky_ave idtheta = %d sigmac = %.5e\n",index_thetas,pclass_sz->sky_averaged_ylims[index_thetas]);
  }
}

if (pclass_sz->sz_verbose >= 3){
  printf("noise map loaded\n");
}
  return  _SUCCESS_;
}

  int read_SO_Qfit(struct class_sz_structure * pclass_sz){
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
      pclass_sz->SO_Qfit   = (double *)malloc(n_data_guess*sizeof(double));
      pclass_sz->SO_thetas   = (double *)malloc(n_data_guess*sizeof(double));

      // char Filepath[_ARGUMENT_LENGTH_MAX_];
      // sprintf(Filepath,
      //         "%s%s",
      //         // "%s%s%s",
      //         "cat ",
      //         // pclass_sz->path_to_class,
      //         "/class_sz_auxiliary_files/SO_files/SOSim_3freq_small_Qfit_comp_test.txt");
      //
      // process = popen(Filepath, "r");
      class_open(process,"class_sz_auxiliary_files/SO_files/SOSim_3freq_small_Qfit_comp_test.txt", "r",pclass_sz->error_message);

      while (fgets(line, sizeof(line)-1, process) != NULL) {
        sscanf(line, "%lf %lf", &this_lnx, &this_lny);

        if((n_data+1) > n_data_guess) {
          n_data_guess *= 2;
          tmp = (double *)realloc(pclass_sz->SO_Qfit,   n_data_guess*sizeof(double));
          class_test(tmp == NULL,
                     pclass_sz->error_message,
                     "Error allocating memory to read SO_Qfit.\n");
          pclass_sz->SO_Qfit = tmp;
          tmp = (double *)realloc(pclass_sz->SO_thetas,   n_data_guess*sizeof(double));
          class_test(tmp == NULL,
                     pclass_sz->error_message,
                     "Error allocating memory to read SO_Qfit.\n");
          pclass_sz->SO_thetas = tmp;
        };


        /* Store */
        pclass_sz->SO_thetas[n_data]   = this_lnx;
        pclass_sz->SO_Qfit[n_data]   = this_lny;
        n_data++;
      }

      // status = pclose(process);
      status = fclose(process);
      class_test(status != 0.,
                 pclass_sz->error_message,
                 "The attempt to launch the external command was unsuccessful. "
                 "Try doing it by hand to check for errors.");

      pclass_sz->SO_Q_size = n_data;

      ///////////////////////////end read Q file


  return  _SUCCESS_;
  }



int read_SO_noise(struct class_sz_structure * pclass_sz){
        //read the Q file
        char line[_LINE_LENGTH_MAX_];
        FILE *process;
        int n_data_guess, n_data = 0;
        double *tmp = NULL;
        double this_lnx,this_lny;
        int status;

        n_data = 0;
        n_data_guess = 100;
        pclass_sz->SO_RMS   = (double *)malloc(n_data_guess*sizeof(double));
        pclass_sz->SO_skyfrac   = (double *)malloc(n_data_guess*sizeof(double));


        class_open(process,"class_sz_auxiliary_files/SO_files/SOSim_3freq_small_RMSTab_comp_test.txt", "r",pclass_sz->error_message);

        while (fgets(line, sizeof(line)-1, process) != NULL) {
          sscanf(line, "%lf %lf", &this_lnx, &this_lny);

          if((n_data+1) > n_data_guess) {
            n_data_guess *= 2;
            tmp = (double *)realloc(pclass_sz->SO_RMS,   n_data_guess*sizeof(double));
            class_test(tmp == NULL,
                       pclass_sz->error_message,
                       "Error allocating memory to read SO_Qfit.\n");
            pclass_sz->SO_RMS = tmp;
            tmp = (double *)realloc(pclass_sz->SO_skyfrac,   n_data_guess*sizeof(double));
            class_test(tmp == NULL,
                       pclass_sz->error_message,
                       "Error allocating memory to read SO_Qfit.\n");
            pclass_sz->SO_skyfrac = tmp;
          };


          /* Store */
          pclass_sz->SO_skyfrac[n_data]   = this_lnx;
          pclass_sz->SO_RMS[n_data]   = this_lny;
          n_data++;
        }

        // status = pclose(process);
        status = fclose(process);
        class_test(status != 0.,
                   pclass_sz->error_message,
                   "The attempt to launch the external command was unsuccessful. "
                   "Try doing it by hand to check for errors.");

        pclass_sz->SO_RMS_size = n_data;

        ///////////////////////////end read Q file


    return  _SUCCESS_;
    }


int read_sz_catalog(struct class_sz_structure * pclass_sz){

      char line[_LINE_LENGTH_MAX_];
      FILE *process;
      int n_data_guess, n_data = 0;
      //double *thetas = NULL,
      double *tmp = NULL;
      double this_lnx,this_lny,this_lnz;
      int status;

      n_data = 0;
      n_data_guess = 100;
      pclass_sz->szcat_z   = (double *)malloc(n_data_guess*sizeof(double));
      pclass_sz->szcat_snr   = (double *)malloc(n_data_guess*sizeof(double));


      class_open(process,pclass_sz->SZ_cat_file, "r",pclass_sz->error_message);

      while (fgets(line, sizeof(line)-1, process) != NULL) {
        sscanf(line, "%lf %lf %lf", &this_lnx, &this_lny, &this_lnz);

        if((n_data+1) > n_data_guess) {
          n_data_guess *= 2;
          tmp = (double *)realloc(pclass_sz->szcat_z,   n_data_guess*sizeof(double));
          class_test(tmp == NULL,
                     pclass_sz->error_message,
                     "Error allocating memory to read szcat_z.\n");
          pclass_sz->szcat_z = tmp;
          tmp = (double *)realloc(pclass_sz->szcat_snr,   n_data_guess*sizeof(double));
          class_test(tmp == NULL,
                     pclass_sz->error_message,
                     "Error allocating memory to read szcat_snr.\n");
          pclass_sz->szcat_snr = tmp;
        };


        /* Store */
        pclass_sz->szcat_z[n_data]   = this_lnx;
        pclass_sz->szcat_snr[n_data]   = this_lnz;
        n_data++;
      }

      // status = pclose(process);
      status = fclose(process);
      class_test(status != 0.,
                 pclass_sz->error_message,
                 "The attempt to launch the external command was unsuccessful. "
                 "Try doing it by hand to check for errors.");

      pclass_sz->szcat_size = n_data;

      ///////////////////////////end read Q file


  return  _SUCCESS_;
  }




int tabulate_ng_bias_contribution_at_z_and_k(struct background * pba,
                                             struct perturbs * ppt,
                                             struct class_sz_structure * pclass_sz){
double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;
int index_z;
pclass_sz->nz_ng_bias = 200; // set in parser



int index_md=ppt->index_md_scalars;
int index_k;

// double k_min = ppt->k[index_md][0]/pba->h;
// double k_max = ppt->k[ppt->k_size[index_md]-1]/pba->h;
pclass_sz->nk_ng_bias = ppt->k_size[index_md];


class_alloc(pclass_sz->array_ln_1pz_ng_bias,sizeof(double *)*pclass_sz->nz_ng_bias,pclass_sz->error_message);
class_alloc(pclass_sz->array_ln_k_ng_bias,sizeof(double *)*pclass_sz->nk_ng_bias,pclass_sz->error_message);

class_alloc(pclass_sz->array_ln_ng_bias_at_z_and_k,
            sizeof(double *)*pclass_sz->nk_ng_bias*pclass_sz->nz_ng_bias,
            pclass_sz->error_message);


for (index_z=0; index_z<pclass_sz->nz_ng_bias; index_z++)
{
      pclass_sz->array_ln_1pz_ng_bias[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->nz_ng_bias-1.); // log(1+z)
}
for (index_k=0; index_k<pclass_sz->nk_ng_bias; index_k++)
{
      pclass_sz->array_ln_k_ng_bias[index_k] = log(ppt->k[index_md][index_k]/pba->h); // in h/Mpc
}


// int index_z_k = 0;
double fNL = pclass_sz->fNL;
// double bh = get_first_order_bias_at_z_and_nu(z,nu,pclass_sz);
double beta_f = 2.*pclass_sz->delta_cSZ; // multiply by (bh-1.) in the "get functoion"
double alpha_k = 1.;


// start collecting transfer functions
double * data;
int size_data;
int number_of_titles = pclass_sz->number_of_titles;

int index_d_tot = pclass_sz->index_d_tot;
int index_phi = pclass_sz->index_phi;
int index_psi = pclass_sz->index_psi;

size_data = number_of_titles*ppt->k_size[index_md];

double tstart, tstop;
int abort;

///////////////////////////////////////////////
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
pba,pclass_sz,ppt,z_min,z_max,beta_f,fNL,size_data,number_of_titles,index_d_tot,index_phi)\
private(tstart, tstop,index_k,index_z,data) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


if (ppt->ic_size[index_md] != 1){
  printf("Please run with only one type of initial conditions to avoid confusion in class_sz.\n");
  exit(0);
}

class_alloc_parallel(data, sizeof(double)*ppt->ic_size[index_md]*size_data, pclass_sz->error_message);

#pragma omp for collapse(2)
for (index_z=0; index_z<pclass_sz->nz_ng_bias; index_z++)
{
for (index_k=0; index_k<pclass_sz->nk_ng_bias; index_k++)
  {



  int index_z_k = index_k * pclass_sz->nz_ng_bias + index_z;



      double z =   exp(pclass_sz->array_ln_1pz_ng_bias[index_z])-1.;
      double kp =  exp(pclass_sz->array_ln_k_ng_bias[index_k]);

      perturb_output_data(pba,
                          ppt,
                          class_format,
                          0., // z_pk....
                          number_of_titles,
                          data);

      // eq. 3 of this: https://arxiv.org/pdf/1810.13424.pdf
      // double alpha_kp = data[index_k*number_of_titles+index_d_tot]/data[index_k*number_of_titles+index_phi];
      //
      // double alpha_kp0 = data[0*number_of_titles+index_d_tot]/data[0*number_of_titles+index_phi];

      double om0 = pclass_sz->Omega_m_0;

      double tk_phi_plus_psi = (data[index_k*number_of_titles+index_phi]+data[index_k*number_of_titles+index_psi])
                                /(data[0*number_of_titles+index_phi]+data[0*number_of_titles+index_psi]);
      // _c_ in m/s
      double c_in_km_per_s = _c_/1000.;
      double k_in_invMpc = kp*pba->h;
      // double k0_in_invMpc = exp(pclass_sz->array_ln_k_ng_bias[0])*pba->h;


      double * pvecback;
      double tau;
      int first_index_back = 0;
      class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);


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

    double D = pvecback[pba->index_bg_D];
    free(pvecback);

    double D_normalized = D*5.*om0/2.;

    // double tk = alpha_kp*3.*om0*pow(100.*pba->h/c_in_km_per_s/k_in_invMpc,2.)/2./D_normalized;
    // double tk0 = alpha_kp0*3.*om0*pow(100.*pba->h/c_in_km_per_s/k0_in_invMpc,2.)/2./D_normalized;

  if (isnan(tk_phi_plus_psi)||isinf(tk_phi_plus_psi) || (tk_phi_plus_psi==0)){
      printf("alpha_kp = %.5e phi = %.5e psi = %.5e k = %.5e z = %.5e\n",
             tk_phi_plus_psi,
             data[index_k*number_of_titles+index_phi],
             data[index_k*number_of_titles+index_psi],
             kp,
             z
           );
      exit(0);
      }
  // else{
  //     if (alpha_kp>0){
  //       printf("alpha>0\n");
  //       exit(0);
  //     }

      // double res = fNL*3.*om0*pow(100.*pba->h/c_in_km_per_s/k_in_invMpc,2.)/tk_phi_plus_psi/D_normalized*pclass_sz->delta_cSZ;
      double res = 3.*om0*pow(100.*pba->h/c_in_km_per_s/k_in_invMpc,2.)/tk_phi_plus_psi/D_normalized*pclass_sz->delta_cSZ;


      // pclass_sz->array_ln_ng_bias_at_z_and_k[index_z_k] = log(fNL*beta_f/alpha_kp);
      // pclass_sz->array_ln_ng_bias_at_z_and_k[index_z_k] = log(res*tk*D_normalized);
      pclass_sz->array_ln_ng_bias_at_z_and_k[index_z_k] = log(res);
      // }


  }
}

#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over zk's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif


free(data);
}


if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

return _SUCCESS_;
                                             }


//Tabulate vrms2 as functions of redshift
 int tabulate_vrms2_from_pk(struct background * pba,
                            struct nonlinear * pnl,
                            struct primordial * ppm,
                            struct class_sz_structure * pclass_sz){

// double z_min,z_max;
// if (pclass_sz->need_sigma==0){
//   class_alloc(pclass_sz->array_redshift,sizeof(double *)*pclass_sz->ndim_redshifts,pclass_sz->error_message);
//   double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
//   // z_min = r8_min(z_min,pclass_sz->z_for_pk_hm);
//   double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
//   // z_max = r8_min(z_max,pclass_sz->z_for_pk_hm);
// }



double * vrms2_var;
class_alloc(vrms2_var,
            sizeof(double *),
            pclass_sz->error_message);


class_alloc(pclass_sz->array_vrms2_at_z,sizeof(double *)*pclass_sz->ndim_redshifts,pclass_sz->error_message);

int index_z;



    for (index_z=0; index_z<pclass_sz->ndim_redshifts; index_z++)
    {
      // if (pclass_sz->need_sigma== 0){
      // pclass_sz->array_redshift[index_z] =
      //                                 log(1.+z_min)
      //                                 +index_z*(log(1.+z_max)-log(1.+z_min))
      //                                 /(pclass_sz->ndim_redshifts-1.); // log(1+z)
      //                               }

            spectra_vrms2(pba,
                          ppm,
                          pnl,
                          pclass_sz,
                          exp(pclass_sz->array_redshift[index_z])-1.,
                          vrms2_var
                          );
          pclass_sz->array_vrms2_at_z[index_z] = log(*vrms2_var);
          // printf("z = %.3e vrms2 = %.3e\n",pclass_sz->array_redshift[index_z],pclass_sz->array_vrms2_at_z[index_z]);

       }

free(vrms2_var);

return _SUCCESS_;
    }



struct Parameters_for_integrand_mean_galaxy_bias{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  struct perturbs * ppt;
  double * pvectsz;
  double * pvecback;
  double z;
};



double integrand_mean_galaxy_bias(double lnM_halo, void *p){

  struct Parameters_for_integrand_mean_galaxy_bias *V = ((struct Parameters_for_integrand_mean_galaxy_bias *) p);

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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];

      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);

      V->pvectsz[V->pclass_sz->index_has_galaxy] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      double M_min;
      double M0;
      double M1_prime;
      double sigma_log10M;
      double nc,ns;

      M_min = V->pclass_sz->M_min_HOD;
      M0 = V->pclass_sz->M0_HOD;
      M1_prime = V->pclass_sz->M1_prime_HOD;
      sigma_log10M = V->pclass_sz->sigma_log10M_HOD;
      // }
      nc = HOD_mean_number_of_central_galaxies(z,V->pvectsz[V->pclass_sz->index_mass_for_galaxies],M_min,sigma_log10M,V->pclass_sz->f_cen_HOD,V->pclass_sz,V->pba);
      ns = HOD_mean_number_of_satellite_galaxies(z,V->pvectsz[V->pclass_sz->index_mass_for_galaxies],nc,M0,V->pclass_sz->alpha_s_HOD,M1_prime,V->pclass_sz,V->pba);
      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ppt,V->pclass_sz);
      double result = hmf*V->pvectsz[V->pclass_sz->index_halo_bias]*(ns+nc);

  return result;

}



struct Parameters_for_integrand_mean_galaxy_number{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  int index_g;
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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];

      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);

      V->pvectsz[V->pclass_sz->index_has_galaxy] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);

      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      double M_min;
      double M0;
      double M1_prime;
      double sigma_log10M;
      double nc,ns;

      M_min = V->pclass_sz->M_min_HOD;
      M0 = V->pclass_sz->M0_HOD;
      M1_prime = V->pclass_sz->M1_prime_HOD;
      sigma_log10M = V->pclass_sz->sigma_log10M_HOD;
      // }
      nc = HOD_mean_number_of_central_galaxies(z,V->pvectsz[V->pclass_sz->index_mass_for_galaxies],M_min,sigma_log10M,V->pclass_sz->f_cen_HOD,V->pclass_sz,V->pba);
      ns = HOD_mean_number_of_satellite_galaxies(z,V->pvectsz[V->pclass_sz->index_mass_for_galaxies],nc,M0,V->pclass_sz->alpha_s_HOD,M1_prime,V->pclass_sz,V->pba);

      if (V->pclass_sz->sz_verbose>3){
        printf("got nc ns hmf %.3e and %.3e %.3e at z = %.3e and m = %.3e\n",nc,ns,hmf,z,exp(lnM_halo));
      }

      double result = hmf*(ns+nc);




  return result;

}



double integrand_mean_galaxy_number_ngal(double lnM_halo, void *p){

  struct Parameters_for_integrand_mean_galaxy_number *V = ((struct Parameters_for_integrand_mean_galaxy_number *) p);

    double M_halo = exp(lnM_halo);

    double z = V->z;
    int index_g = V->index_g;



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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];

      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);

      V->pvectsz[V->pclass_sz->index_has_galaxy] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);

      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      double M_min;
      double M0;
      double M1_prime;
      double sigma_log10M;
      double nc,ns;

      M_min = V->pclass_sz->M_min_HOD_ngal[index_g];
      M0 = V->pclass_sz->M0_HOD_ngal[index_g];
      M1_prime = V->pclass_sz->M1_prime_HOD_ngal[index_g];
      sigma_log10M = V->pclass_sz->sigma_log10M_HOD_ngal[index_g];
      // }
      nc = HOD_mean_number_of_central_galaxies(z,V->pvectsz[V->pclass_sz->index_mass_for_galaxies],M_min,sigma_log10M,V->pclass_sz->f_cen_HOD_ngal[index_g],V->pclass_sz,V->pba);
      ns = HOD_mean_number_of_satellite_galaxies(z,V->pvectsz[V->pclass_sz->index_mass_for_galaxies],nc,M0,V->pclass_sz->alpha_s_HOD_ngal[index_g],M1_prime,V->pclass_sz,V->pba);


    //   printf("%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
    // hmf,nc,ns,M_min, M0, M1_prime, sigma_log10M,V->pclass_sz->f_cen_HOD_ngal[index_g],V->pclass_sz->alpha_s_HOD_ngal[index_g]);
    // exit(0);
      double result = hmf*(ns+nc);




  return result;

}



int tabulate_mean_galaxy_number_density(struct background * pba,
                                        struct nonlinear * pnl,
                                        struct primordial * ppm,
                                        struct class_sz_structure * pclass_sz){

class_alloc(pclass_sz->array_mean_galaxy_number_density,sizeof(double *)*pclass_sz->ndim_redshifts,pclass_sz->error_message);

int index_z;
double r;
double m_min,m_max;

// here we should always integrate over the full mass range,
// since this is a normalization term
//
// if (pclass_sz->hm_consistency == 0){
//   m_min = 1e10; // this has to be the same as the minimal mass at which the counter terms are tabulated
//   m_max = 1e16; // this has to be the same as the maximal mass at which the counter terms are tabulated
// }
// else{
m_min = pclass_sz->M_min_ng_bar;
m_max = pclass_sz->M_max_ng_bar;
// }


if (pclass_sz->sz_verbose> 3){
  printf("-> tabulating mean ngal from hod between %.3e  and %.3e Msun/h\n",m_min,m_max);
}

double * pvecback;
double * pvectsz;


 class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);



 pvectsz[pclass_sz->index_has_galaxy] = 1;
 if (pclass_sz->delta_def_galaxies == 0)
   pvectsz[pclass_sz->index_has_200m] = 1;
 else if (pclass_sz->delta_def_galaxies == 1)
   pvectsz[pclass_sz->index_has_200c] = 1;
 else if (pclass_sz->delta_def_galaxies == 2)
   pvectsz[pclass_sz->index_has_500c] = 1;

for (index_z=0; index_z<pclass_sz->ndim_redshifts; index_z++)
        {
          double z = exp(pclass_sz->array_redshift[index_z])-1.;
          if (pclass_sz->sz_verbose>3){
            printf("tabulating mean ngal from hod between %.3e  and %.3e Msun/h at z = %.3e\n",m_min,m_max,z);
          }

          // at each z, perform the mass integral
          struct Parameters_for_integrand_mean_galaxy_number V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.pclass_sz = pclass_sz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;

          void * params = &V;
          double epsrel=pclass_sz->mass_epsrel_ngbar;
          double epsabs=pclass_sz->mass_epsabs_ngbar;

          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_mean_galaxy_number,
                                               params,pclass_sz->patterson_show_neval);
           if (pclass_sz->sz_verbose>3){
             printf("got %.3e at z = %.3e\n",r,z);
           }

          // if (z< 1e-3)
          //   printf("-> [1gal] ngbar for sample at z = %.3e is %.3e.\n",z,r);
        // here we always impose the consistency condition.
        // add counter terms:
        if (pclass_sz->hm_consistency_ngbar){
         double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
         double I0 = integrand_mean_galaxy_number(log(m_min),params);
         double nmin_umin = nmin*I0/pvectsz[pclass_sz->index_hmf];

         if (pclass_sz->sz_verbose>3){
           printf("counter terms %.3e at z = %.3e\n",nmin,z);
         }

         r += nmin_umin;

         if (pclass_sz->sz_verbose>3){
           printf("after counter terms added got %.3e at z = %.3e\n",r,z);
         }
       }

          pclass_sz->array_mean_galaxy_number_density[index_z] = log(r);
          // if (z< 1e-3)
          //   printf("-> [1gal] ngbar for sample at z = %.3e is %.3e.\n",z,r);
        if (pclass_sz->sz_verbose>3){
          if (z< 1e-3)
            printf("-> [1gal] ngbar for sample at z = %.3e is %.3e.\n",z,r);
        }
       }
 free(pvecback);
 free(pvectsz);
 // exit(0);

return _SUCCESS_;
    }




int tabulate_mean_galaxy_number_density_ngal(struct background * pba,
                                             struct nonlinear * pnl,
                                             struct primordial * ppm,
                                             struct class_sz_structure * pclass_sz){

if (pclass_sz->sz_verbose>0){
  printf("-> [ngal] starting tabulation of mean galaxy number density.\n");
}
class_alloc(pclass_sz->array_mean_galaxy_number_density_ngal,sizeof(double **)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
int index_g;
for (index_g=0;index_g<pclass_sz->galaxy_samples_list_num;index_g++){
 class_alloc(pclass_sz->array_mean_galaxy_number_density_ngal[index_g],sizeof(double *)*pclass_sz->ndim_redshifts,pclass_sz->error_message);
}

if (pclass_sz->sz_verbose>0)
  printf("-> [ngal] mean galaxy number density array allocated.\n");
int index_z;
double r;
double m_min,m_max;

// here we should always integrate over the full mass range,
// since this is a normalization term
//
// if (pclass_sz->hm_consistency == 0){
//   m_min = 1e10; // this has to be the same as the minimal mass at which the counter terms are tabulated
//   m_max = 1e16; // this has to be the same as the maximal mass at which the counter terms are tabulated
// }
// else{
for (index_g=0;index_g<pclass_sz->galaxy_samples_list_num;index_g++){

  if (pclass_sz->sz_verbose>0)
    printf("-> [ngal] starting computation for sample %d.\n",index_g);

m_min = pclass_sz->M_min_ng_bar;
m_max = pclass_sz->M_max_ng_bar;
// }

double * pvecback;
double * pvectsz;


 class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);



 pvectsz[pclass_sz->index_has_galaxy] = 1;
 if (pclass_sz->delta_def_galaxies == 0)
   pvectsz[pclass_sz->index_has_200m] = 1;
 else if (pclass_sz->delta_def_galaxies == 1)
   pvectsz[pclass_sz->index_has_200c] = 1;
 else if (pclass_sz->delta_def_galaxies == 2)
   pvectsz[pclass_sz->index_has_500c] = 1;

for (index_z=0; index_z<pclass_sz->ndim_redshifts; index_z++)
        {


          double z = exp(pclass_sz->array_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_mean_galaxy_number V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.pclass_sz = pclass_sz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;
          V.index_g = index_g;

          void * params = &V;
          double epsrel=pclass_sz->mass_epsrel_ngbar;
          double epsabs=pclass_sz->mass_epsabs_ngbar;

          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_mean_galaxy_number_ngal,
                                               params,pclass_sz->patterson_show_neval);

          // if (z< 1e-3)
          //   printf("-> [ngal] ngbar for sample %d at z = %.3e is %.3e.\n",index_g,z,r);

        // here we always impose the consistency condition.
        // add counter terms:
         double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
         double I0 = integrand_mean_galaxy_number_ngal(log(m_min),params);
         double nmin_umin = nmin*I0/pvectsz[pclass_sz->index_hmf];
         r += nmin_umin;


          pclass_sz->array_mean_galaxy_number_density_ngal[index_g][index_z] = log(r);
          // if (z< 1e-3)
          //   printf("-> [ngal] ngbar for sample %d at z = %.3e is %.3e.\n",index_g,z,r);


        if (pclass_sz->sz_verbose>3){
          if (z< 1e-3)
            printf("-> [ngal] ngbar for sample %d at z = %.3e is %.3e.\n",index_g,z,r);
        }

       }
 free(pvecback);
 free(pvectsz);
 // exit(0);
}

return _SUCCESS_;
    }






int tabulate_mean_galaxy_bias(struct background * pba,
                              struct nonlinear * pnl,
                              struct primordial * ppm,
                              struct perturbs * ppt,
                              struct class_sz_structure * pclass_sz){

class_alloc(pclass_sz->array_mean_galaxy_bias,sizeof(double *)*pclass_sz->ndim_redshifts,pclass_sz->error_message);

int index_z;
double r;
double m_min,m_max;

// here we should always integrate over the full mass range,
// since this is a normalization term

if (pclass_sz->hm_consistency == 0){
  m_min = 1e10; // this has to be the same as the minimal mass at whch the counter terms are tabulated
  m_max = 1e16; // this has to be the same as the maximal mass at whch the counter terms are tabulated
}
else{
m_min = pclass_sz->M1SZ;
m_max = pclass_sz->M2SZ;
}

double * pvecback;
double * pvectsz;


 class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


 // printf("tabulating dndlnM quantities0\n");

for (index_z=0; index_z<pclass_sz->ndim_redshifts; index_z++)
        {
          double z = exp(pclass_sz->array_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_mean_galaxy_bias V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.pclass_sz = pclass_sz;
          V.pba = pba;
          V.ppt = ppt;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;

          void * params = &V;
          double epsrel=pclass_sz->mass_epsrel;
          double epsabs=pclass_sz->mass_epsabs;

          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_mean_galaxy_bias,
                                               params,pclass_sz->patterson_show_neval);

        // here we always impose the consistency condition.
        // add counter terms:
         // double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
         // double I0 = integrand_mean_galaxy_number(log(m_min),params);
         // double nmin_umin = nmin*I0/pvectsz[pclass_sz->index_hmf];
         // r += nmin_umin;
         double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
         double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
         double I0 = integrand_mean_galaxy_bias(log(pclass_sz->m_min_counter_terms),params);
         double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
         r += bmin_umin;



          pclass_sz->array_mean_galaxy_bias[index_z] = log(r/evaluate_mean_galaxy_number_density_at_z(z,pclass_sz));
          // printf("ng = %.8e\n",r);

       }
 free(pvecback);
 free(pvectsz);
 // exit(0);

return _SUCCESS_;
    }




struct Parameters_for_integrand_hmf_counter_terms_b1min{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  struct perturbs * ppt;
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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);

      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      double z_asked = z;
      double  m_asked = M_halo;

      double rho_crit_at_z = V->pclass_sz->Rho_crit_0;
      double Omega_cb = (V->pba->Omega0_cdm + V->pba->Omega0_b);
      double rho_cb = rho_crit_at_z*Omega_cb;

      // here ensure the mass is the halo mass:
      double xout = V->pclass_sz->x_out_truncated_nfw_profile;
      double c_delta_matter;
        if (V->pclass_sz->delta_def_matter_density == 0){
          c_delta_matter = get_c200m_at_m_and_z(M_halo,z,V->pba,V->pclass_sz);
        }
        else if (V->pclass_sz->delta_def_matter_density == 1){
          c_delta_matter = get_c200c_at_m_and_z(M_halo,z,V->pba,V->pclass_sz);
        }
        else if (V->pclass_sz->delta_def_matter_density == 2){
          c_delta_matter = get_c500c_at_m_and_z(M_halo,z,V->pba,V->pclass_sz);
        }
        else if (V->pclass_sz->delta_def_matter_density == 3){
          c_delta_matter = evaluate_cvir_of_mvir(M_halo,z,V->pclass_sz,V->pba);
        }
      M_halo *= 1.;//m_min_fac m_nfw(xout*c_delta_matter)/ m_nfw(c_delta_matter);
      ///done with mass consistency.

      double result = hmf*M_halo/rho_cb;

      // switch off non gaussian bias in all these cases:
      int store_ng_in_bh =  V->pclass_sz->has_ng_in_bh;
      V->pclass_sz->has_ng_in_bh = 0;
      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ppt,V->pclass_sz);
      // restore:
      V->pclass_sz->has_ng_in_bh = store_ng_in_bh;
      // done with that !


      double b1 = V->pvectsz[V->pclass_sz->index_halo_bias];
      result *= b1;




  return result;

}


struct Parameters_for_integrand_psi_b2t{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l;
};



struct Parameters_for_integrand_psi_b2g{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l;
};



struct Parameters_for_integrand_psi_b1g{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct perturbs * ppt;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l;
};


struct Parameters_for_integrand_psi_b1t{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct perturbs * ppt;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l;
};


struct Parameters_for_integrand_psi_b1gt{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct perturbs * ppt;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l1;
  double l2;
};



struct Parameters_for_integrand_psi_b2kg{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l;
};



struct Parameters_for_integrand_psi_b1kg{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct perturbs * ppt;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l;
};

struct Parameters_for_integrand_psi_b1kgg{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct perturbs * ppt;
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
  double l1;
  double l2;
};


struct Parameters_for_integrand_psi_b1kgt{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct perturbs * ppt;
  struct class_sz_structure * pclass_sz;
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
  struct class_sz_structure * pclass_sz;
  struct background * pba;
  double * pvectsz;
  double * pvecback;
  double z;
};

struct Parameters_for_integrand_dcib0dz{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct class_sz_structure * pclass_sz;
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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);


      // request appropriate mass conversion
      V->pvectsz[V->pclass_sz->index_has_electron_pressure] = 1 ;

      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];
      V->pvectsz[V->pclass_sz->index_md] = V->pclass_sz->index_md_dydz;

      evaluate_pressure_profile(0.,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);


      double result = hmf*V->pvectsz[V->pclass_sz->index_pressure_profile];

      // multiply by volume element:
      double H_over_c_in_h_over_Mpc = V->pvecback[V->pba->index_bg_H]/V->pba->h;
      result *= V->pvectsz[V->pclass_sz->index_chi2]/H_over_c_in_h_over_Mpc;

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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);


      V->pvectsz[V->pclass_sz->index_has_cib] = 1;

      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];
      V->pvectsz[V->pclass_sz->index_frequency_for_cib_profile] = index_nu;
      V->pvectsz[V->pclass_sz->index_md] = V->pclass_sz->index_md_dcib0dz;
      V->pvectsz[V->pclass_sz->index_multipole_for_cib_profile] = 1.;
      evaluate_cib_profile(V->pvectsz[V->pclass_sz->index_mass_for_cib],
                           V->pvectsz[V->pclass_sz->index_radius_for_cib],
                           V->pvectsz[V->pclass_sz->index_concentration_for_cib],
                           V->pvecback,V->pvectsz,V->pba,V->pclass_sz);

      double result = hmf*V->pvectsz[V->pclass_sz->index_cib_profile];

      // multiply by volume element:
      double H_over_c_in_h_over_Mpc = V->pvecback[V->pba->index_bg_H]/V->pba->h;
      // result *= V->pvectsz[V->pclass_sz->index_chi2]/H_over_c_in_h_over_Mpc;
      // result *= (1.+z)/H_over_c_in_h_over_Mpc;
      result *= 1./(1.+z)/H_over_c_in_h_over_Mpc;

  return result;

}






double integrand_psi_b1g(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b1g *V = ((struct Parameters_for_integrand_psi_b1g *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell = V->l; // this is actually k


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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);

      V->pvectsz[V->pclass_sz->index_has_galaxy] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      V->pvectsz[V->pclass_sz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->pclass_sz);
      V->pvectsz[V->pclass_sz->index_multipole_for_galaxy_profile] = ell;//pclass_sz->ell[index_l_3];
      // double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      // double kl = (ell+0.5)/chi;
      evaluate_galaxy_profile_2h(ell,V->pvectsz[V->pclass_sz->index_mass_for_galaxies],
                                 V->pvectsz[V->pclass_sz->index_radius_for_galaxies],
                                 V->pvectsz[V->pclass_sz->index_concentration_for_galaxies],
                                 V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      double g = V->pvectsz[V->pclass_sz->index_galaxy_profile];


      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ppt,V->pclass_sz);
      double b1 = V->pvectsz[V->pclass_sz->index_halo_bias];
      double result = hmf*b1*g;



  return result;

}


double integrand_psi_b1kg(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b1kg *V = ((struct Parameters_for_integrand_psi_b1kg *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell = V->l; // this is actually k


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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);
      double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);

      V->pvectsz[V->pclass_sz->index_has_lensing] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      // evaluate_redshift_int_gallens_sources(V->pvectsz,V->pclass_sz);
      // double redshift_int_sources = V->pvectsz[V->pclass_sz->index_W_gallens_sources];
      // // double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      // V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit] = pow(3.*pow(V->pba->H0/V->pba->h,2)/2./V->pclass_sz->Rho_crit_0,-1)*pow((1.+z),1.)/(chi*redshift_int_sources);
      // if (V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit]<0. || isnan(V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit])||isinf(V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit])){
      //   printf("%.3e\n",redshift_int_sources);
      //   printf("0, nan or inf in sigmacrit\n");
      //   exit(0);
      // }




      // V->pvectsz[V->pclass_sz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->pclass_sz);
      // V->pvectsz[V->pclass_sz->index_multipole_for_galaxy_profile] = ell;//pclass_sz->ell[index_l_3];
      // double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      // double kl = (ell+0.5)/chi;
      // V->pvectsz[pclass_sz->index_multipole_for_lensing_profile] = ell;
      evaluate_lensing_profile(ell,V->pvectsz[V->pclass_sz->index_mass_for_matter_density],
                                 V->pvectsz[V->pclass_sz->index_radius_for_matter_density],
                                 V->pvectsz[V->pclass_sz->index_concentration_for_matter_density],
                                 V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      double g = V->pvectsz[V->pclass_sz->index_lensing_profile];

      if (g<0. || isnan(g) || isinf(g)){
        printf("integrand b1kg: %.3e\n",g);
        printf("integrand b1kg: 0, nan or inf in lensing prof\n");
        exit(0);
      }


      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ppt,V->pclass_sz);
      double b1 = V->pvectsz[V->pclass_sz->index_halo_bias];
      double result = hmf*b1*g;



  return result;

}




double integrand_psi_b2g(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b2g *V = ((struct Parameters_for_integrand_psi_b2g *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell = V->l; // this is actually k


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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);


      V->pvectsz[V->pclass_sz->index_has_galaxy] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      V->pvectsz[V->pclass_sz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->pclass_sz);
      // V->pvectsz[V->pclass_sz->index_multipole_for_galaxy_profile] = ell;
      // double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      // double kl = (ell+0.5)/chi;
      evaluate_galaxy_profile_2h(ell,V->pvectsz[V->pclass_sz->index_mass_for_galaxies],
                                 V->pvectsz[V->pclass_sz->index_radius_for_galaxies],
                                 V->pvectsz[V->pclass_sz->index_concentration_for_galaxies],
                                 V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      double g = V->pvectsz[V->pclass_sz->index_galaxy_profile];

      evaluate_halo_bias_b2(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);
      double b2 = V->pvectsz[V->pclass_sz->index_halo_bias_b2];


      double result = hmf*b2*g;


  return result;

}


double integrand_psi_b2kg(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b2kg *V = ((struct Parameters_for_integrand_psi_b2kg *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell = V->l; // this is actually k


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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);
      double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);

      V->pvectsz[V->pclass_sz->index_has_lensing] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];
      //
      // evaluate_redshift_int_gallens_sources(V->pvectsz,V->pclass_sz);
      // double redshift_int_sources = V->pvectsz[V->pclass_sz->index_W_gallens_sources];
      // // double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      // V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit] = pow(3.*pow(V->pba->H0/V->pba->h,2)/2./V->pclass_sz->Rho_crit_0,-1)*pow((1.+z),1.)/(chi*redshift_int_sources);
      // if (isnan(V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit])||isinf(V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit])){
      //   printf("%.3e\n",redshift_int_sources);
      //   printf("nan or inf in sigmacrit\n");
      //   exit(0);
      // }
      // V->pvectsz[V->pclass_sz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->pclass_sz);
      // V->pvectsz[V->pclass_sz->index_multipole_for_galaxy_profile] = ell;
      // double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      // double kl = (ell+0.5)/chi;
      evaluate_lensing_profile(ell,V->pvectsz[V->pclass_sz->index_mass_for_matter_density],
                                 V->pvectsz[V->pclass_sz->index_radius_for_matter_density],
                                 V->pvectsz[V->pclass_sz->index_concentration_for_matter_density],
                                 V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      double g = V->pvectsz[V->pclass_sz->index_lensing_profile];

      evaluate_halo_bias_b2(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);
      double b2 = V->pvectsz[V->pclass_sz->index_halo_bias_b2];


      double result = hmf*b2*g;

      if(isnan(result) || isinf(result)){
        printf("nan or inf in integrand b2k %.3e %.3e %.3e\n",hmf,b2,g);
        exit(0);
      }


  return result;

}



double integrand_psi_b2t(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b2t *V = ((struct Parameters_for_integrand_psi_b2t *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell = V->l; // this is actually k


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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);


      V->pvectsz[V->pclass_sz->index_has_electron_density] = 1;
      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      // V->pvectsz[V->pclass_sz->index_multipole_for_tau_profile] = ell;

      evaluate_tau_profile(ell,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      double t = V->pvectsz[V->pclass_sz->index_tau_profile];

      evaluate_halo_bias_b2(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);
      double b2 = V->pvectsz[V->pclass_sz->index_halo_bias_b2];


      double result = hmf*b2*t;


  return result;

}





double integrand_psi_b1t(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b1t *V = ((struct Parameters_for_integrand_psi_b1t *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell = V->l; // this is actually k

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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);

      V->pvectsz[V->pclass_sz->index_has_electron_density] = 1;
      // V->pvectsz[V->pclass_sz->index_has_matter_density] = 1;

      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];


      // V->pvectsz[V->pclass_sz->index_multipole_for_tau_profile] = ell;
      evaluate_tau_profile(ell,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);

      // double r_delta_matter = V->pvectsz[V->pclass_sz->index_radius_for_matter_density];
      // double c_delta_matter = V->pvectsz[V->pclass_sz->index_concentration_for_matter_density];
      // double k = (ell+0.5)/sqrt(V->pvectsz[V->pclass_sz->index_chi2]);

      double t = V->pvectsz[V->pclass_sz->index_tau_profile];
      // double rhom =  V->pvectsz[V->pclass_sz->index_density_profile];


      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ppt,V->pclass_sz);
      double b1 = V->pvectsz[V->pclass_sz->index_halo_bias];
      double result = hmf*b1*t;
      // double result = hmf*b1*rhom;



  return result;

}

double integrand_psi_b1gt(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b1gt *V = ((struct Parameters_for_integrand_psi_b1gt *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell1 = V->l1; // this is actually k
    double ell2 = V->l2; // this is actually k

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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);

      V->pvectsz[V->pclass_sz->index_has_galaxy] = 1;
      V->pvectsz[V->pclass_sz->index_has_electron_density] = 1;

      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      V->pvectsz[V->pclass_sz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->pclass_sz);
      // V->pvectsz[V->pclass_sz->index_multipole_for_galaxy_profile] = ell1;
      // double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      // double kl = (ell1+0.5)/chi;
      evaluate_galaxy_profile_2h(ell1,V->pvectsz[V->pclass_sz->index_mass_for_galaxies],
                                 V->pvectsz[V->pclass_sz->index_radius_for_galaxies],
                                 V->pvectsz[V->pclass_sz->index_concentration_for_galaxies],
                                 V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      double g = V->pvectsz[V->pclass_sz->index_galaxy_profile];



      // V->pvectsz[V->pclass_sz->index_multipole_for_tau_profile] = ell2;

      evaluate_tau_profile(ell2,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      double t = V->pvectsz[V->pclass_sz->index_tau_profile];




      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ppt,V->pclass_sz);
      double b1 = V->pvectsz[V->pclass_sz->index_halo_bias];
      double result = hmf*b1*g*t;
      if (isnan(result)||isinf(result)){
        printf("tab b1gt : z %.3e k3 %.4e k' %.4e\n",z,ell1,ell2);
        exit(0);
      }


  return result;

}


double integrand_psi_b1kgt(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b1kgt *V = ((struct Parameters_for_integrand_psi_b1kgt *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell1 = V->l1; // this is actually k
    double ell2 = V->l2; // this is actually k

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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);
      double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);

      V->pvectsz[V->pclass_sz->index_has_lensing] = 1;
      V->pvectsz[V->pclass_sz->index_has_electron_density] = 1;

      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      // evaluate_redshift_int_gallens_sources(V->pvectsz,V->pclass_sz);
      // double redshift_int_sources = V->pvectsz[V->pclass_sz->index_W_gallens_sources];
      // // double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      // V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit] = pow(3.*pow(V->pba->H0/V->pba->h,2)/2./V->pclass_sz->Rho_crit_0,-1)*pow((1.+z),1.)/(chi*redshift_int_sources);
      // if (isnan(V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit])||isinf(V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit])){
      //   printf("%.3e\n",redshift_int_sources);
      //   printf("nan or inf in sigmacrit\n");
      //   exit(0);
      // }


      // V->pvectsz[V->pclass_sz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->pclass_sz);
      // V->pvectsz[V->pclass_sz->index_multipole_for_galaxy_profile] = ell1;
      // double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      // double kl = (ell1+0.5)/chi;
      evaluate_lensing_profile(ell1,V->pvectsz[V->pclass_sz->index_mass_for_matter_density],
                                 V->pvectsz[V->pclass_sz->index_radius_for_matter_density],
                                 V->pvectsz[V->pclass_sz->index_concentration_for_matter_density],
                                 V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      double g = V->pvectsz[V->pclass_sz->index_lensing_profile];



      // V->pvectsz[V->pclass_sz->index_multipole_for_tau_profile] = ell2;

      evaluate_tau_profile(ell2,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      double t = V->pvectsz[V->pclass_sz->index_tau_profile];




      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ppt,V->pclass_sz);
      double b1 = V->pvectsz[V->pclass_sz->index_halo_bias];
      double result = hmf*b1*g*t;
      if (isnan(result)||isinf(result)){
        printf("tab b1kgt : z %.3e k3 %.4e k %.4e hmf %.4e b1 %.4e g %.4e t %.4e\n",z,ell1,ell2,hmf,b1,g,t);
        exit(0);
      }


  return result;

}



double integrand_psi_b1kgg(double lnM_halo, void *p){

  struct Parameters_for_integrand_psi_b1kgg *V = ((struct Parameters_for_integrand_psi_b1kgg *) p);

    //double x=exp(ln_x);
    double z = V->z;
    double ell1 = V->l1; // this is actually k
    double ell2 = V->l2; // this is actually k

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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      V->pvectsz[V->pclass_sz->index_chi2] = pow(V->pvecback[V->pba->index_bg_ang_distance]*(1.+z)*V->pba->h,2);
      double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);

      V->pvectsz[V->pclass_sz->index_has_lensing] = 1;
      V->pvectsz[V->pclass_sz->index_has_galaxy] = 1;

      do_mass_conversions(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      // evaluate_redshift_int_gallens_sources(V->pvectsz,V->pclass_sz);
      // double redshift_int_sources = V->pvectsz[V->pclass_sz->index_W_gallens_sources];
      // // double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      // V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit] = pow(3.*pow(V->pba->H0/V->pba->h,2)/2./V->pclass_sz->Rho_crit_0,-1)*pow((1.+z),1.)/(chi*redshift_int_sources);
      // if (isnan(V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit])||isinf(V->pvectsz[V->pclass_sz->index_lensing_Sigma_crit])){
      //   printf("%.3e\n",redshift_int_sources);
      //   printf("nan or inf in sigmacrit\n");
      //   exit(0);
      // }


      // V->pvectsz[V->pclass_sz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->pclass_sz);
      // V->pvectsz[V->pclass_sz->index_multipole_for_galaxy_profile] = ell1;
      // double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      // double kl = (ell1+0.5)/chi;
      evaluate_lensing_profile(ell1,V->pvectsz[V->pclass_sz->index_mass_for_matter_density],
                                 V->pvectsz[V->pclass_sz->index_radius_for_matter_density],
                                 V->pvectsz[V->pclass_sz->index_concentration_for_matter_density],
                                 V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      double g = V->pvectsz[V->pclass_sz->index_lensing_profile];



      // V->pvectsz[V->pclass_sz->index_multipole_for_tau_profile] = ell2;

      // evaluate_tau_profile(ell2,V->pvecback,V->pvectsz,V->pba,V->pclass_sz);
      V->pvectsz[V->pclass_sz->index_mean_galaxy_number_density] = evaluate_mean_galaxy_number_density_at_z(z,V->pclass_sz);

      // printf("nbar = %.5e ell1 = %.5e ell2 = %.5e\n",
      //         V->pvectsz[V->pclass_sz->index_mean_galaxy_number_density],
      //         ell1,ell2);
      // V->pvectsz[V->pclass_sz->index_multipole_for_galaxy_profile] = ell1;
      // double chi = sqrt(V->pvectsz[V->pclass_sz->index_chi2]);
      // double kl = (ell1+0.5)/chi;
      evaluate_galaxy_profile_2h(ell2,V->pvectsz[V->pclass_sz->index_mass_for_galaxies],
                                 V->pvectsz[V->pclass_sz->index_radius_for_galaxies],
                                 V->pvectsz[V->pclass_sz->index_concentration_for_galaxies],
                                 V->pvecback,V->pvectsz,V->pba,V->pclass_sz);

      double t = V->pvectsz[V->pclass_sz->index_galaxy_profile];




      evaluate_halo_bias(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->ppt,V->pclass_sz);
      double b1 = V->pvectsz[V->pclass_sz->index_halo_bias];
      double result = hmf*b1*g*t;
      if (isnan(result)||isinf(result)){
        printf("tab b1kgg : z %.3e k3 %.4e k %.4e hmf %.4e b1 %.4e g %.4e t %.4e\n",z,ell1,ell2,hmf,b1,g,t);
        exit(0);
      }


  return result;

}




int tabulate_psi_b1g(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct perturbs * ppt,
                    struct class_sz_structure * pclass_sz){

class_alloc(pclass_sz->array_psi_b1g_redshift,sizeof(double *)*pclass_sz->n_z_psi_b1g,pclass_sz->error_message);
class_alloc(pclass_sz->array_psi_b1g_multipole,sizeof(double *)*pclass_sz->n_l_psi_b1g,pclass_sz->error_message);

class_alloc(pclass_sz->array_psi_b1g_psi,sizeof(double *)*pclass_sz->n_l_psi_b1g*pclass_sz->n_z_psi_b1g,pclass_sz->error_message);


int index_z, index_l;
double r;
double m_min,m_max;


m_min = pclass_sz->M1SZ; // for the mass integral
m_max = pclass_sz->M2SZ; // for the mass integral
// m_min = pclass_sz->m_min_counter_terms;
// m_max = pclass_sz->m_max_counter_terms;
double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;

// printf("pclass_sz->n_z_psi_b1g = %d\n",pclass_sz->n_z_psi_b1g);

for (index_z=0; index_z<pclass_sz->n_z_psi_b1g; index_z++)
        {

          pclass_sz->array_psi_b1g_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_psi_b1g-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 3e5;


for (index_l=0; index_l<pclass_sz->n_l_psi_b1g; index_l++)
        {

          pclass_sz->array_psi_b1g_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(pclass_sz->n_l_psi_b1g-1.); // log(l)
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
pba,pclass_sz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


#pragma omp for schedule (dynamic)
for (index_l=0; index_l<pclass_sz->n_l_psi_b1g; index_l++)
{
#pragma omp flush(abort)

double l = exp(pclass_sz->array_psi_b1g_multipole[index_l]);

for (index_z=0; index_z<pclass_sz->n_z_psi_b1g; index_z++)
        {

          int index_l_z = index_l * pclass_sz->n_z_psi_b1g + index_z;


          double z = exp(pclass_sz->array_psi_b1g_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b1g V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ppt = ppt;
          V.pclass_sz = pclass_sz;
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
                                               pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_psi_b1g(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

          pclass_sz->array_psi_b1g_psi[index_l_z] = log(r);
       }
     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
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


int tabulate_psi_b1kg(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct perturbs * ppt,
                    struct class_sz_structure * pclass_sz){
if(pclass_sz->sz_verbose>0){
  printf("->tabulating psi b1kg\n");
}


class_alloc(pclass_sz->array_psi_b1kg_redshift,sizeof(double *)*pclass_sz->n_z_psi_b1kg,pclass_sz->error_message);
class_alloc(pclass_sz->array_psi_b1kg_multipole,sizeof(double *)*pclass_sz->n_l_psi_b1kg,pclass_sz->error_message);

class_alloc(pclass_sz->array_psi_b1kg_psi,sizeof(double *)*pclass_sz->n_l_psi_b1kg*pclass_sz->n_z_psi_b1kg,pclass_sz->error_message);


int index_z, index_l;
double r;
double m_min,m_max;


m_min = pclass_sz->M1SZ; // for the mass integral
m_max = pclass_sz->M2SZ; // for the mass integral
// m_min = pclass_sz->m_min_counter_terms;
// m_max = pclass_sz->m_max_counter_terms;
double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;

// printf("pclass_sz->n_z_psi_b1kg = %d\n",pclass_sz->n_z_psi_b1kg);

for (index_z=0; index_z<pclass_sz->n_z_psi_b1kg; index_z++)
        {

          pclass_sz->array_psi_b1kg_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_psi_b1kg-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 3e5;


for (index_l=0; index_l<pclass_sz->n_l_psi_b1kg; index_l++)
        {

          pclass_sz->array_psi_b1kg_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(pclass_sz->n_l_psi_b1kg-1.); // log(l)
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
pba,pclass_sz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


#pragma omp for schedule (dynamic)
for (index_l=0; index_l<pclass_sz->n_l_psi_b1kg; index_l++)
{
#pragma omp flush(abort)

double l = exp(pclass_sz->array_psi_b1kg_multipole[index_l]);

for (index_z=0; index_z<pclass_sz->n_z_psi_b1kg; index_z++)
        {

          int index_l_z = index_l * pclass_sz->n_z_psi_b1kg + index_z;


          double z = exp(pclass_sz->array_psi_b1kg_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b1kg V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ppt = ppt;
          V.pclass_sz = pclass_sz;
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
                                               integrand_psi_b1kg,
                                               params,
                                               pclass_sz->patterson_show_neval);
if (r < 0. || isnan(r)||isinf(r)){
printf("tab b1kg after int0 : z %.3e r %.3e k1 %.3e\n",z,r,l);
exit(0);
}

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_psi_b1kg(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }

          pclass_sz->array_psi_b1kg_psi[index_l_z] = log(r);
       }
     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
         printf("In %s: time spent in parallel region b1kg (loop over l's) = %e s for thread %d\n",
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
                    struct class_sz_structure * pclass_sz){

class_alloc(pclass_sz->array_psi_b2g_redshift,sizeof(double *)*pclass_sz->n_z_psi_b2g,pclass_sz->error_message);
class_alloc(pclass_sz->array_psi_b2g_multipole,sizeof(double *)*pclass_sz->n_l_psi_b2g,pclass_sz->error_message);

class_alloc(pclass_sz->array_psi_b2g_psi,sizeof(double *)*pclass_sz->n_l_psi_b2g*pclass_sz->n_z_psi_b2g,pclass_sz->error_message);


int index_z, index_l;
double r;
double m_min,m_max;


m_min = pclass_sz->M1SZ; // for the mass integral
m_max = pclass_sz->M2SZ; // for the mass integral
// m_min = pclass_sz->m_min_counter_terms;
// m_max = pclass_sz->m_max_counter_terms;
double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;


for (index_z=0; index_z<pclass_sz->n_z_psi_b2g; index_z++)
        {

          pclass_sz->array_psi_b2g_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_psi_b2g-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 3e5;


for (index_l=0; index_l<pclass_sz->n_l_psi_b2g; index_l++)
        {

          pclass_sz->array_psi_b2g_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(pclass_sz->n_l_psi_b2g-1.); // log(l)
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
pba,pclass_sz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


#pragma omp for schedule (dynamic)
for (index_l=0; index_l<pclass_sz->n_l_psi_b2g; index_l++)
{
#pragma omp flush(abort)

double l = exp(pclass_sz->array_psi_b2g_multipole[index_l]);

for (index_z=0; index_z<pclass_sz->n_z_psi_b2g; index_z++)
        {

          int index_l_z = index_l * pclass_sz->n_z_psi_b2g + index_z;


          double z = exp(pclass_sz->array_psi_b2g_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b2g V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.pclass_sz = pclass_sz;
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
                                               pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_psi_b2g(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r += bmin_umin;
  }

          pclass_sz->array_psi_b2g_psi[index_l_z] = log(1.+r);
       }
     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
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




int tabulate_psi_b2kg(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct class_sz_structure * pclass_sz){

class_alloc(pclass_sz->array_psi_b2kg_redshift,sizeof(double *)*pclass_sz->n_z_psi_b2kg,pclass_sz->error_message);
class_alloc(pclass_sz->array_psi_b2kg_multipole,sizeof(double *)*pclass_sz->n_l_psi_b2kg,pclass_sz->error_message);

class_alloc(pclass_sz->array_psi_b2kg_psi,sizeof(double *)*pclass_sz->n_l_psi_b2kg*pclass_sz->n_z_psi_b2kg,pclass_sz->error_message);


int index_z, index_l;
double r;
double m_min,m_max;


m_min = pclass_sz->M1SZ; // for the mass integral
m_max = pclass_sz->M2SZ; // for the mass integral
// m_min = pclass_sz->m_min_counter_terms;
// m_max = pclass_sz->m_max_counter_terms;
double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;


for (index_z=0; index_z<pclass_sz->n_z_psi_b2kg; index_z++)
        {

          pclass_sz->array_psi_b2kg_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_psi_b2kg-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 3e5;


for (index_l=0; index_l<pclass_sz->n_l_psi_b2kg; index_l++)
        {

          pclass_sz->array_psi_b2kg_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(pclass_sz->n_l_psi_b2kg-1.); // log(l)
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
pba,pclass_sz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


#pragma omp for schedule (dynamic)
for (index_l=0; index_l<pclass_sz->n_l_psi_b2kg; index_l++)
{
#pragma omp flush(abort)

double l = exp(pclass_sz->array_psi_b2kg_multipole[index_l]);

for (index_z=0; index_z<pclass_sz->n_z_psi_b2kg; index_z++)
        {

          int index_l_z = index_l * pclass_sz->n_z_psi_b2kg + index_z;


          double z = exp(pclass_sz->array_psi_b2kg_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b2kg V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.pclass_sz = pclass_sz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;
          V.l = l;

          void * params = &V;
          double epsrel=1e-4;
          double epsabs=1e-100;



          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_psi_b2kg,
                                               params,
                                               pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_psi_b2kg(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     r += bmin_umin;
  }

          pclass_sz->array_psi_b2kg_psi[index_l_z] = log(1.+r);
      if(isnan(pclass_sz->array_psi_b2kg_psi[index_l_z]) || isinf(pclass_sz->array_psi_b2kg_psi[index_l_z])){
        printf("nan or inf in tabulate b2k %.3e %.3e\n",r,pvectsz[pclass_sz->index_z]);
        exit(0);
      }


       }
     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
         printf("In %s: time spent in parallel region b2kg (loop over l's) = %e s for thread %d\n",
                __func__,tstop-tstart,omp_get_thread_num());
     #endif
 free(pvecback);
 free(pvectsz);
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;
    }





int tabulate_n5k_F1(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct class_sz_structure * pclass_sz){

// printf("tabulating n5k_F\n");
class_alloc(pclass_sz->array_n5k_F1_l,sizeof(int *)*pclass_sz->n_l_n5k,pclass_sz->error_message);
class_alloc(pclass_sz->array_n5k_F1_k,sizeof(double *)*pclass_sz->n_k_n5k,pclass_sz->error_message);
class_alloc(pclass_sz->array_n5k_F1_F,sizeof(double *)*pclass_sz->n_l_n5k,pclass_sz->error_message);


int index_k, index_l;
double r;
// const int Nl_n5k = pclass_sz->n_l_n5k;
if (pclass_sz->n_l_n5k != 103)
  {
    printf("wrong l dimension for n5k\n");
    exit(0);
  }
// printf("%d ",Nl_n5k);
// exit(0);
double l_n5k[103] = {   2.,    3.,    4.,    5.,    6.,    7.,    8.,    9.,   10.,
         11.,   12.,   13.,   14.,   15.,   16.,   17.,   18.,   19.,
         20.,   21.,   23.,   24.,   25.,   27.,   28.,   30.,   32.,
         33.,   35.,   37.,   39.,   42.,   44.,   46.,   49.,   52.,
         55.,   58.,   61.,   64.,   68.,   72.,   76.,   80.,   85.,
         90.,   95.,  100.,  106.,  111.,  118.,  124.,  131.,  139.,
        146.,  155.,  163.,  172.,  182.,  192.,  203.,  215.,  227.,
        239.,  253.,  267.,  282.,  298.,  314.,  332.,  350.,  370.,
        391.,  413.,  436.,  460.,  486.,  513.,  542.,  572.,  604.,
        638.,  673.,  711.,  751.,  793.,  837.,  884.,  933.,  986.,
       1041., 1099., 1160., 1225., 1294., 1366., 1443., 1523., 1608.,
       1698., 1793., 1894., 2000.};


for (index_l=0; index_l<pclass_sz->n_l_n5k; index_l++)
        {

          pclass_sz->array_n5k_F1_l[index_l] = l_n5k[index_l];

        }

double k_min = pclass_sz->k_min_n5k;//1.e-4;
double k_max = pclass_sz->k_max_n5k;//1e2;
for (index_k=0; index_k<pclass_sz->n_k_n5k; index_k++)
        {

          pclass_sz->array_n5k_F1_k[index_k] =
                                          log(k_min)
                                          +index_k*(log(k_max)-log(k_min))
                                          /(pclass_sz->n_k_n5k-1.); // log(l)
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
pba,pclass_sz,ppm,pnl,k_max,k_min)\
private(tstart, tstop,index_k,index_l,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


 // #pragma omp for collapse(2)
 // for (index_l=0; index_l<pclass_sz->n_l_n5k; index_l++)
 // {
 // for (index_k=0; index_k<pclass_sz->n_k_n5k; index_k++)
 //   {



#pragma omp for schedule (dynamic)
for (index_l=0; index_l<pclass_sz->n_l_n5k; index_l++)
{
#pragma omp flush(abort)

double sumk = 0.;
int l = pclass_sz->array_n5k_F1_l[index_l];
// for (index_k=0; index_k<pclass_sz->n_k_n5k; index_k++)
//   {
//
//           double k = exp(pclass_sz->array_n5k_F1_k[index_k]);
//           // int l = pclass_sz->array_n5k_F1_l[index_l];
//
//
//           double chi_min = pclass_sz->chi_min_n5k_samp_fftw;//1e0;//pclass_sz->l_min_samp_fftw; //precision parameter
//           double chi_max = pclass_sz->chi_max_n5k_samp_fftw;//7e3;//pclass_sz->l_max_samp_fftw; //precision parameter
//
//           const int N = pclass_sz->N_samp_fftw; //precision parameter
//           int ichi;
//           double chi[N], Pchi[N];
//           for (ichi=0; ichi<N; ichi++){
//             chi[ichi] =  exp(log(chi_min)+ichi/(N-1.)*(log(chi_max)-log(chi_min)));
//             double zchi = get_n5k_z_of_chi(chi[ichi],pclass_sz);
//             Pchi[ichi] = sqrt(get_n5k_pk_at_z_and_k(zchi,k,pclass_sz))*get_n5k_cl_K1_at_chi(chi[ichi],pclass_sz);
//             // printf("Pchi = %.3e\n",Pchi[ichi]);
//           }
//
//           double chit[N], Pchit[N];
//         /* Compute the function
//          *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
//          * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
//          * in this notation.  The input k-values must be logarithmically spaced.  The
//          * resulting xi_l^m(r) will be evaluated at the dual r-values
//          *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
//           fftlog_ComputeXiLMsloz(l, 0, N, chi,  Pchi, chit, Pchit,pclass_sz);
//           double F1 = 2.*_PI_*_PI_*pwl_value_1d(N,chit,Pchit,k);
//           fftlog_ComputeXiLMsloz(l, 0, N, chi,  Pchi, chit, Pchit,pclass_sz);
//           double F2 = 2.*_PI_*_PI_*pwl_value_1d(N,chit,Pchit,k);
//           double intk = F1*F2*k*k;
//           double dlk = (log(k_max)-log(k_min))/(pclass_sz->n_k_n5k-1.);
//           sumk +=  intk*k*dlk;
//
//        }
    // if (pclass_sz->sz_verbose>5){
    //   printf("ell = %d sumk = %.3e\n",l,sumk);
    // }

  struct Parameters_for_integrand_n5k_at_k V;
  // V.pnl = pnl;
  // V.ppm = ppm;
  V.pclass_sz = pclass_sz;
  // V.pba = pba;
  // V.m = m;
  // V.z = z;
  // V.rd = rd;
  V.l = l; // TBC!
  // V.pvectsz = Pvectsz;
  // V.pvecback = Pvecback;

  void * params = &V;


  double epsrel= pclass_sz->integrand_n5k_epsrel;
  double epsabs= pclass_sz->integrand_n5k_epsabs;
  int show_neval = pclass_sz->patterson_show_neval;

  // printf("evaluating integrand_n5k_at_k at one k for l = %d\n",l);
  // double result = integrand_n5k_at_k(log(1e-3), params);
  // printf("Result of integrand_n5k_at_k: %.5e\n", result);
  // exit(0);


  //integral of density profile.
  double sumk_patterson= Integrate_using_Patterson_adaptive(log(pclass_sz->k_min_n5k),
                                                        log(pclass_sz->k_max_n5k),
                                                        epsrel, epsabs,
                                                        integrand_n5k_at_k,
                                                        params,show_neval);

  // printf("sumk %.8e sumkpatt %.8e\n",sumk,sumk_patterson);
sumk = sumk_patterson;

pclass_sz->array_n5k_F1_F[index_l] = sumk*2./_PI_;

     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
         printf("In %s: time spent in parallel region n5k (loop over l's) = %e s for thread %d\n",
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
                    struct class_sz_structure * pclass_sz){

class_alloc(pclass_sz->array_psi_b2t_redshift,sizeof(double *)*pclass_sz->n_z_psi_b2t,pclass_sz->error_message);
class_alloc(pclass_sz->array_psi_b2t_multipole,sizeof(double *)*pclass_sz->n_l_psi_b2t,pclass_sz->error_message);

class_alloc(pclass_sz->array_psi_b2t_psi,sizeof(double *)*pclass_sz->n_l_psi_b2t*pclass_sz->n_z_psi_b2t,pclass_sz->error_message);


int index_z, index_l;
double r;
double m_min,m_max;


m_min = pclass_sz->M1SZ; // for the mass integral
m_max = pclass_sz->M2SZ; // for the mass integral
// m_min = pclass_sz->m_min_counter_terms;
// m_max = pclass_sz->m_max_counter_terms;
double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;


for (index_z=0; index_z<pclass_sz->n_z_psi_b2t; index_z++)
        {

          pclass_sz->array_psi_b2t_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_psi_b2t-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 3e5;


for (index_l=0; index_l<pclass_sz->n_l_psi_b2t; index_l++)
        {

          pclass_sz->array_psi_b2t_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(pclass_sz->n_l_psi_b2t-1.); // log(l)
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
pba,pclass_sz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


#pragma omp for schedule (dynamic)
for (index_l=0; index_l<pclass_sz->n_l_psi_b2t; index_l++)
{
#pragma omp flush(abort)

double l = exp(pclass_sz->array_psi_b2t_multipole[index_l]);

for (index_z=0; index_z<pclass_sz->n_z_psi_b2t; index_z++)
        {

          int index_l_z = index_l * pclass_sz->n_z_psi_b2t + index_z;


          double z = exp(pclass_sz->array_psi_b2t_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b2t V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.pclass_sz = pclass_sz;
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
                                               pclass_sz->patterson_show_neval);
          // printf("%.8e %.8e %.8e\n",z,l,r);

double ct = 0.;
   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b2min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_psi_b2t(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias_b2];
     double ct = bmin_umin;
     r += bmin_umin;

     if (r<-1.){
       printf("%.8e %.8e int = %.8e ct = %.8e %.8e %.8e\n",z,l,r-ct,ct,nmin,r);
       printf("sort out the tabulation of b2g.\n");
     }
  }


          pclass_sz->array_psi_b2t_psi[index_l_z] = log(1.+r);
       }
     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
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
                    struct perturbs * ppt,
                    struct class_sz_structure * pclass_sz){

class_alloc(pclass_sz->array_psi_b1t_redshift,sizeof(double *)*pclass_sz->n_z_psi_b1t,pclass_sz->error_message);
class_alloc(pclass_sz->array_psi_b1t_multipole,sizeof(double *)*pclass_sz->n_l_psi_b1t,pclass_sz->error_message);

class_alloc(pclass_sz->array_psi_b1t_psi,sizeof(double *)*pclass_sz->n_l_psi_b1t*pclass_sz->n_z_psi_b1t,pclass_sz->error_message);


int index_z, index_l;
double r;
double m_min,m_max;


m_min = pclass_sz->M1SZ; // for the mass integral
m_max = pclass_sz->M2SZ; // for the mass integral
// m_min = pclass_sz->m_min_counter_terms;
// m_max = pclass_sz->m_max_counter_terms;
double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;

// printf("pclass_sz->n_z_psi_b1g = %d\n",pclass_sz->n_z_psi_b1g);

for (index_z=0; index_z<pclass_sz->n_z_psi_b1t; index_z++)
        {

          pclass_sz->array_psi_b1t_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_psi_b1t-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 3e5;


for (index_l=0; index_l<pclass_sz->n_l_psi_b1t; index_l++)
        {

          pclass_sz->array_psi_b1t_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(pclass_sz->n_l_psi_b1t-1.); // log(l)
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
pba,pclass_sz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


#pragma omp for schedule (dynamic)
for (index_l=0; index_l<pclass_sz->n_l_psi_b1t; index_l++)
{
#pragma omp flush(abort)

double l = exp(pclass_sz->array_psi_b1t_multipole[index_l]);

for (index_z=0; index_z<pclass_sz->n_z_psi_b1t; index_z++)
        {

          int index_l_z = index_l * pclass_sz->n_z_psi_b1t + index_z;


          double z = exp(pclass_sz->array_psi_b1t_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b1t V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ppt = ppt;
          V.pclass_sz = pclass_sz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;
          V.l = l;

          void * params = &V;
          double epsrel=1e-3;
          double epsabs=1e-100;

          pvectsz[pclass_sz->index_has_electron_density] = 1; //

          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_psi_b1t,
                                               params,
                                               pclass_sz->patterson_show_neval);

           if (isnan(r)||isinf(r)){
             printf("nan or inf in psib1t tab at z %.4e l %.3e got %.3e\n",
             z,l,r);
             exit(0);
           }

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_psi_b1t(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r += bmin_umin;
     // double ct_over_int  = bmin_umin/(r-bmin_umin);
     // if (ct_over_int>0.1)
     // printf("z = %.8e l = %.8e int = %.8e ct = %.8e ct/int = %.8e\n",
     // z,l,r-bmin_umin,bmin_umin,ct_over_int);
     if (isnan(bmin_umin)||isinf(bmin_umin)){
       printf("nan or inf in psib1t tab at z %.4e l %.3e got bmin_umin %.3e\n",
       z,l,bmin_umin);
       exit(0);
     }
  }

  if (r==0) r = 1e-200;

  if (r<0){
    printf("negative r in psib1t tab at z %.4e l %.3e got r %.3e\n",
    z,l,r);
    exit(0);
  }
          pclass_sz->array_psi_b1t_psi[index_l_z] = log(r);

       }
     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
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
                    struct perturbs * ppt,
                    struct class_sz_structure * pclass_sz){

class_alloc(pclass_sz->array_psi_b1gt_redshift,sizeof(double *)*pclass_sz->n_z_psi_b1gt,pclass_sz->error_message);
class_alloc(pclass_sz->array_psi_b1gt_multipole,sizeof(double *)*pclass_sz->n_l_psi_b1gt,pclass_sz->error_message);

class_alloc(pclass_sz->array_psi_b1gt_psi,sizeof(double *)*pclass_sz->n_z_psi_b1gt,pclass_sz->error_message);
int index_z, index_l1,index_l2;
for (index_z=0;index_z<pclass_sz->n_z_psi_b1gt;index_z++){
  class_alloc(pclass_sz->array_psi_b1gt_psi[index_z],sizeof(double *)*pclass_sz->n_l_psi_b1gt*pclass_sz->n_l_psi_b1gt,pclass_sz->error_message);

}



double r;
double m_min,m_max;


m_min = pclass_sz->M1SZ; // for the mass integral
m_max = pclass_sz->M2SZ; // for the mass integral
// m_min = pclass_sz->m_min_counter_terms;
// m_max = pclass_sz->m_max_counter_terms;
double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;

// printf("pclass_sz->n_z_psi_b1g = %d\n",pclass_sz->n_z_psi_b1g);

for (index_z=0; index_z<pclass_sz->n_z_psi_b1gt; index_z++)
        {

          pclass_sz->array_psi_b1gt_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_psi_b1gt-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 1e3;

int index_l;
for (index_l=0; index_l<pclass_sz->n_l_psi_b1gt; index_l++)
        {

          pclass_sz->array_psi_b1gt_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(pclass_sz->n_l_psi_b1gt-1.); // log(l)
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
pba,pclass_sz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l1,index_l2,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


// printf("doing well\n");

#pragma omp for collapse(3)
//#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_psi_b1gt; index_z++)
{
//#pragma omp flush(abort)
  for (index_l1=0; index_l1<pclass_sz->n_l_psi_b1gt; index_l1++)
  {
    for (index_l2=0; index_l2<pclass_sz->n_l_psi_b1gt; index_l2++)
    {
// #pragma omp flush(abort)


          double l1 = exp(pclass_sz->array_psi_b1gt_multipole[index_l1]);
          double l2 = exp(pclass_sz->array_psi_b1gt_multipole[index_l2]);

          int index_l1_l2 = index_l2 * pclass_sz->n_l_psi_b1gt + index_l1;


          double z = exp(pclass_sz->array_psi_b1gt_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b1gt V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ppt = ppt;
          V.pclass_sz = pclass_sz;
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
                                               pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_psi_b1gt(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];
     r += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }
  if (r==0){
    r = 1e-100;
  }
  if (r < 0. || isnan(r)||isinf(r)){
    printf("tab b1gt after int : z %.3e r %.3e k1 %.3e k2 %.3e\n",z,r,l1,l2);
    // exit(0);
  }
          pclass_sz->array_psi_b1gt_psi[index_z][index_l1_l2] = log(r);

          // printf("pclass_sz->array_psi_b1t_psi[%d] = %.5e\n",index_l_z,r);

       }
     }
   }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
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




int tabulate_psi_b1kgt(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct perturbs * ppt,
                    struct class_sz_structure * pclass_sz){

class_alloc(pclass_sz->array_psi_b1kgt_redshift,sizeof(double *)*pclass_sz->n_z_psi_b1kgt,pclass_sz->error_message);
class_alloc(pclass_sz->array_psi_b1kgt_multipole,sizeof(double *)*pclass_sz->n_l_psi_b1kgt,pclass_sz->error_message);

class_alloc(pclass_sz->array_psi_b1kgt_psi,sizeof(double *)*pclass_sz->n_z_psi_b1kgt,pclass_sz->error_message);
int index_z, index_l1,index_l2;
for (index_z=0;index_z<pclass_sz->n_z_psi_b1kgt;index_z++){
  class_alloc(pclass_sz->array_psi_b1kgt_psi[index_z],sizeof(double *)*pclass_sz->n_l_psi_b1kgt*pclass_sz->n_l_psi_b1kgt,pclass_sz->error_message);

}



double r;
double m_min,m_max;


m_min = pclass_sz->M1SZ; // for the mass integral
m_max = pclass_sz->M2SZ; // for the mass integral
// m_min = pclass_sz->m_min_counter_terms;
// m_max = pclass_sz->m_max_counter_terms;
double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;

// printf("pclass_sz->n_z_psi_b1g = %d\n",pclass_sz->n_z_psi_b1g);

for (index_z=0; index_z<pclass_sz->n_z_psi_b1kgt; index_z++)
        {

          pclass_sz->array_psi_b1kgt_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_psi_b1kgt-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 1e3;

int index_l;
for (index_l=0; index_l<pclass_sz->n_l_psi_b1kgt; index_l++)
        {

          pclass_sz->array_psi_b1kgt_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(pclass_sz->n_l_psi_b1kgt-1.); // log(l)
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
pba,pclass_sz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l1,index_l2,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


// printf("doing well\n");

#pragma omp for collapse(3)
//#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_psi_b1kgt; index_z++)
{
//#pragma omp flush(abort)
  for (index_l1=0; index_l1<pclass_sz->n_l_psi_b1kgt; index_l1++)
  {
    for (index_l2=0; index_l2<pclass_sz->n_l_psi_b1kgt; index_l2++)
    {
// #pragma omp flush(abort)


          double l1 = exp(pclass_sz->array_psi_b1kgt_multipole[index_l1]);
          double l2 = exp(pclass_sz->array_psi_b1kgt_multipole[index_l2]);

          int index_l1_l2 = index_l2 * pclass_sz->n_l_psi_b1kgt + index_l1;


          double z = exp(pclass_sz->array_psi_b1kgt_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b1kgt V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ppt = ppt;
          V.pclass_sz = pclass_sz;
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
                                               integrand_psi_b1kgt,
                                               params,
                                               pclass_sz->patterson_show_neval);
if (r < 0. || isnan(r)||isinf(r)){
printf("tab b1kgt after int0 : z %.3e r %.3e k1 %.3e k2 %.3e\n",z,r,l1,l2);
}

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_psi_b1kgt(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];


     r += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }
  if (r==0){
    r = 1e-100;
  }
  if (r < 0. || isnan(r)||isinf(r)){
    printf("tab b1kgt after int1 : z %.3e r %.3e k1 %.3e k2 %.3e\n",z,r,l1,l2);
    // exit(0);
  }
          pclass_sz->array_psi_b1kgt_psi[index_z][index_l1_l2] = log(r);

          // printf("pclass_sz->array_psi_b1t_psi[%d] = %.5e\n",index_l_z,r);

       }
     }
   }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
         printf("In %s: time spent in parallel region b1kgt (loop over llz's) = %e s for thread %d\n",
                __func__,tstop-tstart,omp_get_thread_num());
     #endif
 free(pvecback);
 free(pvectsz);
}
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region
return _SUCCESS_;
    }



int tabulate_psi_b1kgg(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct perturbs * ppt,
                    struct class_sz_structure * pclass_sz){
       if (pclass_sz->sz_verbose > 0)
        printf("in tabulate_psi_b1kgg\n");
class_alloc(pclass_sz->array_psi_b1kgg_redshift,sizeof(double *)*pclass_sz->n_z_psi_b1kgg,pclass_sz->error_message);
class_alloc(pclass_sz->array_psi_b1kgg_multipole,sizeof(double *)*pclass_sz->n_l_psi_b1kgg,pclass_sz->error_message);

class_alloc(pclass_sz->array_psi_b1kgg_psi,sizeof(double *)*pclass_sz->n_z_psi_b1kgg,pclass_sz->error_message);
int index_z, index_l1,index_l2;
for (index_z=0;index_z<pclass_sz->n_z_psi_b1kgg;index_z++){
  class_alloc(pclass_sz->array_psi_b1kgg_psi[index_z],sizeof(double *)*pclass_sz->n_l_psi_b1kgg*pclass_sz->n_l_psi_b1kgg,pclass_sz->error_message);

}



double r;
double m_min,m_max;


m_min = pclass_sz->M1SZ; // for the mass integral
m_max = pclass_sz->M2SZ; // for the mass integral
// m_min = pclass_sz->m_min_counter_terms;
// m_max = pclass_sz->m_max_counter_terms;
double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;

// printf("pclass_sz->n_z_psi_b1g = %d\n",pclass_sz->n_z_psi_b1g);

for (index_z=0; index_z<pclass_sz->n_z_psi_b1kgg; index_z++)
        {

          pclass_sz->array_psi_b1kgg_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_psi_b1kgg-1.); // log(1+z)
        }

// parallelize ver l
double l_min = 1.e-3;
double l_max = 1e3;

int index_l;
for (index_l=0; index_l<pclass_sz->n_l_psi_b1kgg; index_l++)
        {

          pclass_sz->array_psi_b1kgg_multipole[index_l] =
                                      log(l_min)
                                      +index_l*(log(l_max)-log(l_min))
                                      /(pclass_sz->n_l_psi_b1kgg-1.); // log(l)
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
pba,pclass_sz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_l1,index_l2,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);


// printf("doing well\n");

#pragma omp for collapse(3)
//#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_psi_b1kgg; index_z++)
{
//#pragma omp flush(abort)
  for (index_l1=0; index_l1<pclass_sz->n_l_psi_b1kgg; index_l1++)
  {
    for (index_l2=0; index_l2<pclass_sz->n_l_psi_b1kgg; index_l2++)
    {
// #pragma omp flush(abort)


          double l1 = exp(pclass_sz->array_psi_b1kgg_multipole[index_l1]);
          double l2 = exp(pclass_sz->array_psi_b1kgg_multipole[index_l2]);

          int index_l1_l2 = index_l2 * pclass_sz->n_l_psi_b1kgg + index_l1;


          double z = exp(pclass_sz->array_psi_b1kgg_redshift[index_z])-1.;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_psi_b1kgt V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ppt = ppt;
          V.pclass_sz = pclass_sz;
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
                                               integrand_psi_b1kgg,
                                               params,
                                               pclass_sz->patterson_show_neval);
if (r < 0. || isnan(r)||isinf(r)){
printf("tab b1kgg after int0 : z %.3e r %.3e k1 %.3e k2 %.3e\n",z,r,l1,l2);
}

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
     double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
     double bmin = get_hmf_counter_term_b1min_at_z(pvectsz[pclass_sz->index_z],pclass_sz)*nmin;
     double I0 = integrand_psi_b1kgg(log(pclass_sz->m_min_counter_terms),params);
     double bmin_umin = bmin*I0/pvectsz[pclass_sz->index_hmf]/pvectsz[pclass_sz->index_halo_bias];


     r += bmin_umin;
     // printf("counter terms done r_m_1\n");
  }
  if (r==0){
    r = 1e-100;
  }
  if (r < 0. || isnan(r)||isinf(r)){
    printf("tab b1kgg after int1 : z %.3e r %.3e k1 %.3e k2 %.3e\n",z,r,l1,l2);
    // exit(0);
  }
          pclass_sz->array_psi_b1kgg_psi[index_z][index_l1_l2] = log(r);

          // printf("pclass_sz->array_psi_b1t_psi[%d] = %.5e\n",index_l_z,r);

       }
     }
   }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
         printf("In %s: time spent in parallel region b1kgt (loop over llz's) = %e s for thread %d\n",
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
                    struct class_sz_structure * pclass_sz){

class_alloc(pclass_sz->array_dydz_redshift,sizeof(double *)*pclass_sz->n_z_dydz,pclass_sz->error_message);
class_alloc(pclass_sz->array_dydz_at_z,sizeof(double *)*pclass_sz->n_z_dydz,pclass_sz->error_message);




double r;
double m_min,m_max;


m_min = pclass_sz->M1SZ; // for the mass integral
m_max = pclass_sz->M2SZ; // for the mass integral

double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;
int index_z;


for (index_z=0; index_z<pclass_sz->n_z_dydz; index_z++)
        {

          pclass_sz->array_dydz_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_dydz-1.); // log(1+z)
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
pba,pclass_sz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);

#pragma omp for collapse(1)
for (index_z=0; index_z<pclass_sz->n_z_dydz; index_z++)
{

          double z = exp(pclass_sz->array_dydz_redshift[index_z])-1.;;

          // at each z, perform the mass integral
          struct Parameters_for_integrand_dydz V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.pclass_sz = pclass_sz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;


          void * params = &V;
          double epsrel=1e-6;
          double epsabs=1e-100;



          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_dydz,
                                               params,
                                               pclass_sz->patterson_show_neval);
          r *= 1./pow(pclass_sz->Tcmb_gNU,1)/1.e6;

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
       double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
       double I0 = integrand_dydz(log(pclass_sz->m_min_counter_terms),params)/pow(pclass_sz->Tcmb_gNU,1)/1.e6;
       double nmin_umin = nmin*I0/pvectsz[pclass_sz->index_hmf];
       if (nmin_umin<0.){
         printf("%.5e %.5e %.5e\n",nmin,I0,nmin_umin);
         exit(0);
       }
       r += nmin_umin;
  }

          pclass_sz->array_dydz_at_z[index_z] = log(r);

     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
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
                    struct class_sz_structure * pclass_sz){

class_alloc(pclass_sz->array_dcib0dz_redshift,sizeof(double *)*pclass_sz->n_z_dcib0dz,pclass_sz->error_message);
class_alloc(pclass_sz->array_dcib0dz_nu,sizeof(double *)*pclass_sz->n_nu_dcib0dz,pclass_sz->error_message);

class_alloc(pclass_sz->array_dcib0dz_at_z_nu,sizeof(double *)*pclass_sz->n_z_dcib0dz*pclass_sz->n_nu_dcib0dz,pclass_sz->error_message);




double r;
double m_min,m_max;


m_min = pclass_sz->M1SZ; // for the mass integral
m_max = pclass_sz->M2SZ; // for the mass integral

double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;
int index_z;


for (index_z=0; index_z<pclass_sz->n_z_dcib0dz; index_z++)
        {

          pclass_sz->array_dcib0dz_redshift[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_dcib0dz-1.); // log(1+z)
        }

double nu_min = pclass_sz->freq_min;
double nu_max = pclass_sz->freq_max;

int index_nu;
for (index_nu=0; index_nu<pclass_sz->n_nu_dcib0dz; index_nu++)
        {

          pclass_sz->array_dcib0dz_nu[index_nu] =
                                      log(nu_min)
                                      +index_nu*(log(nu_max)-log(nu_min))
                                      /(pclass_sz->n_nu_dcib0dz-1.); // log(nu)
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
pba,pclass_sz,ppm,pnl,z_min,z_max,m_min,m_max)\
private(tstart, tstop,index_z,index_nu,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);

#pragma omp for collapse(2)
for (index_z=0; index_z<pclass_sz->n_z_dcib0dz; index_z++)
{
  for (index_nu=0; index_nu<pclass_sz->n_nu_dcib0dz; index_nu++)
  {
          double z = exp(pclass_sz->array_dcib0dz_redshift[index_z])-1.;;
          double nu = exp(pclass_sz->array_dcib0dz_nu[index_nu]);

          int index_z_nu = index_nu * pclass_sz->n_z_dcib0dz + index_z;


          // at each z, perform the mass integral
          struct Parameters_for_integrand_dcib0dz V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.pclass_sz = pclass_sz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;
          V.index_nu = index_nu;

          void * params = &V;
          double epsrel=1e-6;
          double epsabs=1e-100;



          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_dcib0dz,
                                               params,
                                               pclass_sz->patterson_show_neval);

   if (pclass_sz->M1SZ == pclass_sz->m_min_counter_terms)  {
       double nmin = get_hmf_counter_term_nmin_at_z(pvectsz[pclass_sz->index_z],pclass_sz);
       double I0 = integrand_dcib0dz(log(pclass_sz->m_min_counter_terms),params);
       double nmin_umin = nmin*I0/pvectsz[pclass_sz->index_hmf];
       r += nmin_umin;
     // printf("counter terms done r_m_1\n");
  }

          pclass_sz->array_dcib0dz_at_z_nu[index_z_nu] = log(r);

       }
     }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
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
                                    struct perturbs * ppt,
                                    struct class_sz_structure * pclass_sz){
// this will only be executed if hm_consistency==1
class_alloc(pclass_sz->array_hmf_counter_terms_b1min,sizeof(double *)*pclass_sz->n_z_hmf_counter_terms,pclass_sz->error_message);

int index_z;
double r;
double m_min,m_max;


m_min = pclass_sz->m_min_counter_terms;
m_max = pclass_sz->m_max_counter_terms;

double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;

double * pvecback;
double * pvectsz;




 class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);




for (index_z=0; index_z<pclass_sz->n_z_hmf_counter_terms; index_z++)
        {
          double z = exp(pclass_sz->array_redshift_hmf_counter_terms[index_z])-1.;

          // at each z, perform the mass integral
          struct Parameters_for_integrand_hmf_counter_terms_b1min V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.ppt = ppt;
          V.pclass_sz = pclass_sz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;

          void * params = &V;
          double epsrel=pclass_sz->mass_epsrel;
          double epsabs=pclass_sz->mass_epsabs;

          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_hmf_counter_terms_b1min,
                                               params,pclass_sz->patterson_show_neval);

          // here is (1-r) = n_min*m_min/rho_cb at z
          double rho_crit_at_z = pclass_sz->Rho_crit_0;//pvectsz[pclass_sz->index_Rho_crit];
          double Omega_cb = (pba->Omega0_cdm + pba->Omega0_b);//*pow(1.+z,3.);
          double rho_cb = rho_crit_at_z*Omega_cb;
          double n_min =  get_hmf_counter_term_nmin_at_z(z,pclass_sz);

      // here ensure the mass is the halo mass:
      double xout = pclass_sz->x_out_truncated_nfw_profile;
      double c_delta_matter;
        if (pclass_sz->delta_def_matter_density == 0){
          c_delta_matter = get_c200m_at_m_and_z(m_min,z,pba,pclass_sz);
        }
        else if (pclass_sz->delta_def_matter_density == 1){
          c_delta_matter = get_c200c_at_m_and_z(m_min,z,pba,pclass_sz);
        }
        else if (pclass_sz->delta_def_matter_density == 2){
          c_delta_matter = get_c500c_at_m_and_z(m_min,z,pba,pclass_sz);
        }
        else if (pclass_sz->delta_def_matter_density == 3){
          c_delta_matter = evaluate_cvir_of_mvir(m_min,z,pclass_sz,pba);
        }
      double m_min_fac = 1.;//m_nfw(xout*c_delta_matter)/ m_nfw(c_delta_matter);
      ///done with mass consistency.


          double b1_min = (1.-r)*rho_cb/m_min/m_min_fac/n_min;
          pclass_sz->array_hmf_counter_terms_b1min[index_z] = b1_min;

       }
 free(pvecback);
 free(pvectsz);

return _SUCCESS_;
    }

///// b2 min


struct Parameters_for_integrand_hmf_counter_terms_b2min{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct class_sz_structure * pclass_sz;
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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);

      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      double z_asked = z;
      double  m_asked = M_halo;
      double rho_crit_at_z = V->pclass_sz->Rho_crit_0;
      double Omega_cb = (V->pba->Omega0_cdm + V->pba->Omega0_b);
      double rho_cb = rho_crit_at_z*Omega_cb;

      // here ensure the mass is the halo mass:
      double xout = V->pclass_sz->x_out_truncated_nfw_profile;
      double c_delta_matter;
        if (V->pclass_sz->delta_def_matter_density == 0){
          c_delta_matter = get_c200m_at_m_and_z(M_halo,z,V->pba,V->pclass_sz);
        }
        else if (V->pclass_sz->delta_def_matter_density == 1){
          c_delta_matter = get_c200c_at_m_and_z(M_halo,z,V->pba,V->pclass_sz);
        }
        else if (V->pclass_sz->delta_def_matter_density == 2){
          c_delta_matter = get_c500c_at_m_and_z(M_halo,z,V->pba,V->pclass_sz);
        }
        else if (V->pclass_sz->delta_def_matter_density == 3){
          c_delta_matter = evaluate_cvir_of_mvir(M_halo,z,V->pclass_sz,V->pba);
        }
      M_halo *= 1.;//m_min_fac m_nfw(xout*c_delta_matter)/ m_nfw(c_delta_matter);
      ///done with mass consistency.

      double result = hmf*M_halo/rho_cb;

      evaluate_halo_bias_b2(V->pvecback,V->pvectsz,V->pba,V->ppm,V->pnl,V->pclass_sz);
      double b2 = V->pvectsz[V->pclass_sz->index_halo_bias_b2];
      result *= b2;

  return result;

}



int tabulate_hmf_counter_terms_b2min(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct class_sz_structure * pclass_sz){
// this will only be executed if hm_consistency==1
class_alloc(pclass_sz->array_hmf_counter_terms_b2min,sizeof(double *)*pclass_sz->n_z_hmf_counter_terms,pclass_sz->error_message);

int index_z;
double r;
double m_min,m_max;

m_min = pclass_sz->m_min_counter_terms;
m_max = pclass_sz->m_max_counter_terms;


double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;


double * pvecback;
double * pvectsz;




 class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);




for (index_z=0; index_z<pclass_sz->n_z_hmf_counter_terms; index_z++)
        {
          double z = exp(pclass_sz->array_redshift_hmf_counter_terms[index_z])-1.;

          // at each z, perform the mass integral
          struct Parameters_for_integrand_hmf_counter_terms_b2min V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.pclass_sz = pclass_sz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;

          void * params = &V;
          double epsrel=pclass_sz->mass_epsrel;
          double epsabs=pclass_sz->mass_epsabs;

          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_hmf_counter_terms_b2min,
                                               params,pclass_sz->patterson_show_neval);

          // here is (1-r) = n_min*m_min/rho_cb at z
          double rho_crit_at_z = pclass_sz->Rho_crit_0;//pvectsz[pclass_sz->index_Rho_crit];
          double Omega_cb = (pba->Omega0_cdm + pba->Omega0_b);//*pow(1.+z,3.);
          double rho_cb = rho_crit_at_z*Omega_cb;
          double n_min =  get_hmf_counter_term_nmin_at_z(z,pclass_sz);

      // here ensure the mass is the halo mass:
      double xout = pclass_sz->x_out_truncated_nfw_profile;
      double c_delta_matter;
        if (pclass_sz->delta_def_matter_density == 0){
          c_delta_matter = get_c200m_at_m_and_z(m_min,z,pba,pclass_sz);
        }
        else if (pclass_sz->delta_def_matter_density == 1){
          c_delta_matter = get_c200c_at_m_and_z(m_min,z,pba,pclass_sz);
        }
        else if (pclass_sz->delta_def_matter_density == 2){
          c_delta_matter = get_c500c_at_m_and_z(m_min,z,pba,pclass_sz);
        }
        else if (pclass_sz->delta_def_matter_density == 3){
          c_delta_matter = evaluate_cvir_of_mvir(m_min,z,pclass_sz,pba);
        }
      double m_min_fac = 1.;//m_nfw(xout*c_delta_matter)/ m_nfw(c_delta_matter);
      ///done with mass consistency.

          double b2_min = -r*rho_cb/m_min/m_min_fac/n_min;
          pclass_sz->array_hmf_counter_terms_b2min[index_z] = b2_min;

       }
 free(pvecback);
 free(pvectsz);

return _SUCCESS_;
    }




//// b2 min


struct Parameters_for_integrand_hmf_counter_terms_nmin{
  struct nonlinear * pnl;
  struct primordial * ppm;
  struct class_sz_structure * pclass_sz;
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




      V->pvectsz[V->pclass_sz->index_z] = z;
      V->pvectsz[V->pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                            *pow(_Mpc_over_m_,1)
                                            *pow(_c_,2)
                                            *V->pvecback[V->pba->index_bg_rho_crit]
                                            /pow(V->pba->h,2);

      double omega = V->pvecback[V->pba->index_bg_Omega_m];
      V->pvectsz[V->pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);


      evaluate_HMF_at_logM_and_z(lnM_halo,z,V->pvecback,V->pvectsz,V->pba,V->pnl,V->pclass_sz);

      double hmf = V->pvectsz[V->pclass_sz->index_hmf];

      double z_asked = z;
      // double  m_asked = M_halo;

      // here ensure the mass is the halo mass:
      double xout = V->pclass_sz->x_out_truncated_nfw_profile;
      double c_delta_matter;
        if (V->pclass_sz->delta_def_matter_density == 0){
          c_delta_matter = get_c200m_at_m_and_z(M_halo,z,V->pba,V->pclass_sz);
        }
        else if (V->pclass_sz->delta_def_matter_density == 1){
          c_delta_matter = get_c200c_at_m_and_z(M_halo,z,V->pba,V->pclass_sz);
        }
        else if (V->pclass_sz->delta_def_matter_density == 2){
          c_delta_matter = get_c500c_at_m_and_z(M_halo,z,V->pba,V->pclass_sz);
        }
        else if (V->pclass_sz->delta_def_matter_density == 3){
          c_delta_matter = evaluate_cvir_of_mvir(M_halo,z,V->pclass_sz,V->pba);
        }
      M_halo *= 1.;// m_min_fac m_nfw(xout*c_delta_matter)/ m_nfw(c_delta_matter);
      ///done with mass consistency.

      double rho_crit_at_z = V->pclass_sz->Rho_crit_0;
      double Omega_cb = (V->pba->Omega0_cdm + V->pba->Omega0_b);
      double rho_cb = rho_crit_at_z*Omega_cb;


      double result = hmf*M_halo/rho_cb;

  return result;

}


int tabulate_hmf_counter_terms_nmin(struct background * pba,
                                    struct nonlinear * pnl,
                                    struct primordial * ppm,
                                    struct class_sz_structure * pclass_sz){
// this will only be executed if hm_consistency==1


class_alloc(pclass_sz->array_hmf_counter_terms_nmin,sizeof(double *)*pclass_sz->n_z_hmf_counter_terms,pclass_sz->error_message);
class_alloc(pclass_sz->array_redshift_hmf_counter_terms,sizeof(double *)*pclass_sz->n_z_hmf_counter_terms,pclass_sz->error_message);

int index_z;
double r;
double m_min,m_max;

m_min = pclass_sz->m_min_counter_terms;
m_max = pclass_sz->m_max_counter_terms;

double z_min = pclass_sz->z1SZ;
double z_max = pclass_sz->z2SZ;

// printf("z_min = %.8e\n",z_min);

double * pvecback;
double * pvectsz;




 class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);




for (index_z=0; index_z<pclass_sz->n_z_hmf_counter_terms; index_z++)
        {

          pclass_sz->array_redshift_hmf_counter_terms[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_hmf_counter_terms-1.); // log(1+z)

          double z = exp(pclass_sz->array_redshift_hmf_counter_terms[index_z])-1.;

          // at each z, perform the mass integral
          struct Parameters_for_integrand_hmf_counter_terms_nmin V;
          V.pnl = pnl;
          V.ppm = ppm;
          V.pclass_sz = pclass_sz;
          V.pba = pba;
          V.pvectsz = pvectsz;
          V.pvecback = pvecback;
          V.z = z;

          void * params = &V;
          double epsrel=pclass_sz->mass_epsrel;
          double epsabs=pclass_sz->mass_epsabs;

          r=Integrate_using_Patterson_adaptive(log(m_min), log(m_max),
                                               epsrel, epsabs,
                                               integrand_hmf_counter_terms_nmin,
                                               params,pclass_sz->patterson_show_neval);

          // here is (1-r) = n_min*m_min/rho_cb at z
          double rho_crit_at_z = pclass_sz->Rho_crit_0;
          double Omega_cb = (pba->Omega0_cdm + pba->Omega0_b);
          double rho_cb = rho_crit_at_z*Omega_cb;

      // here ensure the mass is the halo mass:
      double xout = pclass_sz->x_out_truncated_nfw_profile;
      double c_delta_matter;
        if (pclass_sz->delta_def_matter_density == 0){
          c_delta_matter = get_c200m_at_m_and_z(m_min,z,pba,pclass_sz);
        }
        else if (pclass_sz->delta_def_matter_density == 1){
          c_delta_matter = get_c200c_at_m_and_z(m_min,z,pba,pclass_sz);
        }
        else if (pclass_sz->delta_def_matter_density == 2){
          c_delta_matter = get_c500c_at_m_and_z(m_min,z,pba,pclass_sz);
        }
        else if (pclass_sz->delta_def_matter_density == 3){
          c_delta_matter = evaluate_cvir_of_mvir(m_min,z,pclass_sz,pba);
        }
      double m_min_fac = 1.;//m_nfw(xout*c_delta_matter)/ m_nfw(c_delta_matter);
      ///done with mass consistency.

          double n_min = (1.-r)*rho_cb/m_min/m_min_fac;
          pclass_sz->array_hmf_counter_terms_nmin[index_z] = n_min;

       }
 free(pvecback);
 free(pvectsz);

return _SUCCESS_;
    }




struct Parameters_for_integrand_mass_L_sat{
  double nu;
  double z;
  double M_host;
  struct class_sz_structure * pclass_sz;
};


double integrand_patterson_L_sat(double lnM_sub, void *p){
  struct Parameters_for_integrand_mass_L_sat *V = ((struct Parameters_for_integrand_mass_L_sat *) p);

  double M_sub = exp(lnM_sub);
  double nu = V->nu;
  double z = V->z;
  double M_host = V->M_host;


  double L_gal_at_nu;
  if (V->pclass_sz->use_maniyar_cib_model){
      double sfrI;
      double sfrII;
      sfrI = evaluate_galaxy_luminosity(z, M_sub, nu, V->pclass_sz);
      sfrII =  evaluate_galaxy_luminosity(z, M_host, nu, V->pclass_sz)*M_sub/M_host;

      // L_gal_at_nu = sfrI;//r8_min(sfrI,sfrII);
      L_gal_at_nu = r8_min(sfrI,sfrII);

      }
  else{
      L_gal_at_nu = evaluate_galaxy_luminosity(z, M_sub, nu, V->pclass_sz);
      }
double dNdlnMs;
if (V->pclass_sz->use_maniyar_cib_model){
  dNdlnMs = subhalo_hmf_dndlnMs(M_host/(1.-V->pclass_sz->maniyar_cib_fsub),M_sub,V->pclass_sz);
}
else{
  dNdlnMs = subhalo_hmf_dndlnMs(M_host,M_sub,V->pclass_sz);
}
  double result = L_gal_at_nu*dNdlnMs;

  if (isnan(result) ||  isinf(result)){
  printf("result z = %.5e m = %.5e nu = %.5e Lsat integrand = %.5e L_gal_at_nu = %.5e\n",
                                                                                        z,
                                                                                        exp(lnM_sub),
                                                                                        nu,
                                                                                        result,
                                                                                        L_gal_at_nu);
  exit(0);
  }

  return result;
}


int tabulate_L_sat_at_z_m_nu(struct background * pba,
                             struct class_sz_structure * pclass_sz){

if (
      pclass_sz->has_tSZ_cib_1h
    + pclass_sz->has_tSZ_cib_2h
    + pclass_sz->has_cib_cib_1h
    + pclass_sz->has_cib_monopole
    + pclass_sz->has_cib_shotnoise
    + pclass_sz->has_dcib0dz
    + pclass_sz->has_cib_cib_2h
    + pclass_sz->has_lens_cib_1h
    + pclass_sz->has_lens_cib_2h
    + pclass_sz->has_gallens_cib_1h
    + pclass_sz->has_gallens_cib_2h
    + pclass_sz->has_gal_cib_1h
    + pclass_sz->has_gal_cib_2h
    == _FALSE_
    )
return 0;

if (pclass_sz->sz_verbose>1){
  printf("Tabulating Lsat.\n");
  printf("n_nu_L_sat = %d\n",pclass_sz->n_nu_L_sat);
  printf("nu_min = %.4e\n",pclass_sz->freq_min);
  printf("nu_max = %.4e\n",pclass_sz->freq_max);
}

// printf("pclass_sz->n_nu_L_sat = %d %d\n",pclass_sz->n_nu_L_sat,pclass_sz->n_z_psi_b1gt);

class_alloc(pclass_sz->array_L_sat_at_M_z_nu,sizeof(double *)*pclass_sz->n_nu_L_sat,pclass_sz->error_message);
int index_nu, index_M,index_z;
class_alloc(pclass_sz->array_nu_L_sat,sizeof(double *)*pclass_sz->n_nu_L_sat,pclass_sz->error_message);
double nu_min  = pclass_sz->freq_min;
double nu_max  = pclass_sz->freq_max;

for (index_nu=0;index_nu<pclass_sz->n_nu_L_sat;index_nu++){
  class_alloc(pclass_sz->array_L_sat_at_M_z_nu[index_nu],sizeof(double *)*pclass_sz->n_z_L_sat*pclass_sz->n_m_L_sat,pclass_sz->error_message);
      pclass_sz->array_nu_L_sat[index_nu] =
                                      log(nu_min)
                                      +index_nu*(log(nu_max)-log(nu_min))
                                      /(pclass_sz->n_nu_L_sat-1.); // log(1+z)

}

  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_L_sat);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_L_sat);
  double logM_min = r8_min(log(pclass_sz->M1SZ/pba->h),log(pclass_sz->M1SZ_L_sat)); //in Msun
  double logM_max = r8_max(log(pclass_sz->M2SZ/pba->h),log(pclass_sz->M2SZ_L_sat)); //in Msun


if (pclass_sz->sz_verbose>1){
  printf("Tabulating Lsat.\n");
  printf("n_m_L_sat = %d\n",pclass_sz->n_m_L_sat);
  printf("M_min = %.4e\n",exp(logM_min));
  printf("M_max = %.4e\n", exp(logM_max));

  printf("Tabulating Lsat for z.\n");
  printf("n_z_L_sat = %d\n", pclass_sz->n_z_L_sat);
  printf("z_min = %.4e\n", z_min);
  printf("z_max = %.4e\n", z_max);
}

// if (
//       pclass_sz->has_tSZ_cib_1h
//     + pclass_sz->has_tSZ_cib_2h
//     + pclass_sz->has_cib_cib_1h
//     // + pclass_sz->has_cib_monopole
//     + pclass_sz->has_cib_cib_2h
//     + pclass_sz->has_lens_cib_1h
//     + pclass_sz->has_lens_cib_2h
//     + pclass_sz->has_gal_cib_1h
//     + pclass_sz->has_gal_cib_2h
//     == _FALSE_
//   ){



  class_alloc(pclass_sz->array_z_L_sat,sizeof(double *)*pclass_sz->n_z_L_sat,pclass_sz->error_message);
  class_alloc(pclass_sz->array_m_L_sat,sizeof(double *)*pclass_sz->n_m_L_sat,pclass_sz->error_message);


for (index_z=0; index_z<pclass_sz->n_z_L_sat; index_z++)
{
      pclass_sz->array_z_L_sat[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_L_sat-1.); // log(1+z)

}

for (index_M=0; index_M<pclass_sz->n_m_L_sat; index_M++)
{
      pclass_sz->array_m_L_sat[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(pclass_sz->n_m_L_sat-1.); //log(R)


}

// }




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
pba,pclass_sz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_z,index_nu,index_M,pvecback,pvectsz,r) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


 class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
   int i;
   for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

 class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pclass_sz->error_message);

#pragma omp for collapse(3)
for (index_nu=0; index_nu<pclass_sz->n_nu_L_sat; index_nu++)
{
  for (index_M=0; index_M<pclass_sz->n_m_L_sat; index_M++)
  {
    for (index_z=0; index_z<pclass_sz->n_z_L_sat; index_z++)
    {

          int index_M_z = index_z * pclass_sz->n_m_L_sat + index_M;

      //
      // pclass_sz->array_z_L_sat[index_z] =
      //                                 log(1.+z_min)
      //                                 +index_z*(log(1.+z_max)-log(1.+z_min))
      //                                 /(pclass_sz->n_z_L_sat-1.); // log(1+z)
      //
      // pclass_sz->array_m_L_sat[index_M] =
      //                               logM_min
      //                               +index_M*(logM_max-logM_min)
      //                               /(pclass_sz->n_m_L_sat-1.); //log(R)
      //


      double z =   exp(pclass_sz->array_z_L_sat[index_z])-1.;
      double logM =   pclass_sz->array_m_L_sat[index_M];
      double lnMs_min;
      if (pclass_sz->M_min_subhalo_in_Msun>0.){
      lnMs_min = log(pclass_sz->M_min_subhalo_in_Msun);
      }
      else{
      lnMs_min = log(pclass_sz->M_min_HOD_cib);
      }
      double lnMs_max;
      if (pclass_sz->use_maniyar_cib_model == 1)
        lnMs_max = log(exp(logM)*(1.-pclass_sz->maniyar_cib_fsub));
      else
        lnMs_max = logM;//log(1e11);

      if (lnMs_max<=lnMs_min){
      r = 0.;
      }

      else{
      double epsrel = pclass_sz->epsrel_L_sat;
      double epsabs = pclass_sz->epsabs_L_sat;

      struct Parameters_for_integrand_mass_L_sat V;
      V.nu = exp(pclass_sz->array_nu_L_sat[index_nu]);
      V.z = z;
      V.pclass_sz = pclass_sz;
      if (pclass_sz->use_maniyar_cib_model == 1)
        V.M_host = exp(logM)*(1.-pclass_sz->maniyar_cib_fsub);
      else
        V.M_host = exp(logM);


      void * params = &V;
      params = &V;


          r=Integrate_using_Patterson_adaptive(lnMs_min, lnMs_max,
                                               epsrel, epsabs,
                                               integrand_patterson_L_sat,
                                               params,pclass_sz->patterson_show_neval);
      }

          if (r==0. || r<1e-100){
            r = 1e-100;
          }
          pclass_sz->array_L_sat_at_M_z_nu[index_nu][index_M_z] = log(1.+r);

          if (isnan(pclass_sz->array_L_sat_at_M_z_nu[index_nu][index_M_z])){
            printf("nan in interp L_sat table\n");
          printf("-->z = %.8e M = %.8e expL = %.8e expnu = %.8e nu = %.5e r = %.5e\n",z,
                                                                                  exp(logM),pclass_sz->array_L_sat_at_M_z_nu[index_nu][index_M_z],
                                                                                  exp(pclass_sz->array_nu_L_sat[index_nu]),
                                                                                  pclass_sz->array_nu_L_sat[index_nu],
                                                                                  r);
            exit(0);
          }
          if (isinf(pclass_sz->array_L_sat_at_M_z_nu[index_nu][index_M_z])){
            printf("inf in interp L_sat table\n");
            printf("r = %.3e\n",r);
          printf("-->z = %.8e M = %.8e expL = %.8e expnu = %.8e nu = %.5e r = %.5e\n",z,
                                                                                  exp(logM),pclass_sz->array_L_sat_at_M_z_nu[index_nu][index_M_z],
                                                                                  exp(pclass_sz->array_nu_L_sat[index_nu]),
                                                                                  pclass_sz->array_nu_L_sat[index_nu],
                                                                                  r);
            exit(0);
          }


       }
     }
   }
     #ifdef _OPENMP
       tstop = omp_get_wtime();
       if (pclass_sz->sz_verbose > 0)
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

// int tabulate_L_sat_at_nu_and_nu_prime(struct background * pba,
//                                       struct class_sz_structure * pclass_sz){

// // if (
// //       pclass_sz->has_tSZ_cib_1h
// //     + pclass_sz->has_tSZ_cib_2h
// //     + pclass_sz->has_cib_cib_1h
// //     // + pclass_sz->has_cib_monopole
// //     + pclass_sz->has_cib_cib_2h
// //     + pclass_sz->has_lens_cib_1h
// //     + pclass_sz->has_lens_cib_2h
// //     + pclass_sz->has_gal_cib_1h
// //     + pclass_sz->has_gal_cib_2h
// //     == _FALSE_
// //     )
// // return 0;

//   //Array of z
//   double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_L_sat);
//   double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_L_sat);
//   int index_z;

//   double tstart, tstop;
//   int index_l;

//   // double * pvecback;
//   // double * pvectsz;
//   int abort;

//   //Array of M in Msun
//   double logM_min = r8_min(log(pclass_sz->M1SZ/pba->h),log(pclass_sz->M1SZ_L_sat)); //in Msun
//   double logM_max = r8_max(log(pclass_sz->M2SZ/pba->h),log(pclass_sz->M2SZ_L_sat)); //in Msun
//   int index_M;

//   int index_z_M = 0;

//   double ** array_L_sat_at_z_and_M_at_nu;
//   // double ** array_L_sat_at_z_and_M_at_nu_prime;

//   // class_alloc(pclass_sz->array_z_L_sat,sizeof(double *)*pclass_sz->n_z_L_sat,pclass_sz->error_message);
//   // class_alloc(pclass_sz->array_m_L_sat,sizeof(double *)*pclass_sz->n_m_L_sat,pclass_sz->error_message);


//   class_alloc(pclass_sz->array_L_sat_at_z_and_M_at_nu,
//               pclass_sz->cib_frequency_list_num*sizeof(double *),
//               pclass_sz->error_message);

// int index_nu;
// for (index_nu=0;index_nu<pclass_sz->cib_frequency_list_num;index_nu++){

// class_alloc(pclass_sz->array_L_sat_at_z_and_M_at_nu[index_nu],
//             pclass_sz->n_z_L_sat*pclass_sz->n_m_L_sat*sizeof(double),
//             pclass_sz->error_message);


// class_alloc(array_L_sat_at_z_and_M_at_nu,
//             pclass_sz->n_z_L_sat*sizeof(double *),
//             pclass_sz->error_message);



// for (index_l=0;
//      index_l<pclass_sz->n_z_L_sat;
//      index_l++)
// {
//   class_alloc(array_L_sat_at_z_and_M_at_nu[index_l],
//               pclass_sz->n_m_L_sat*sizeof(double),
//               pclass_sz->error_message);
// }

// /* initialize error management flag */
// abort = _FALSE_;
// /* beginning of parallel region */


// int number_of_threads= 1;
// #ifdef _OPENMP
// #pragma omp parallel
//   {
//     number_of_threads = omp_get_num_threads();
//   }
// #endif

// #pragma omp parallel \
// shared(abort,index_nu,index_z_M,\
// pba,pclass_sz,z_min,z_max,logM_min,logM_max)\
// private(tstart, tstop,index_M,index_z) \
// num_threads(number_of_threads)
// {

// #ifdef _OPENMP
//   tstart = omp_get_wtime();
// #endif


// #pragma omp for schedule (dynamic)
// for (index_z=0; index_z<pclass_sz->n_z_L_sat; index_z++)
// {

// #pragma omp flush(abort)

// for (index_M=0; index_M<pclass_sz->n_m_L_sat; index_M++)
// {
//       // pclass_sz->array_z_L_sat[index_z] =
//       //                                 log(1.+z_min)
//       //                                 +index_z*(log(1.+z_max)-log(1.+z_min))
//       //                                 /(pclass_sz->n_z_L_sat-1.); // log(1+z)
//       //
//       // pclass_sz->array_m_L_sat[index_M] =
//       //                               logM_min
//       //                               +index_M*(logM_max-logM_min)
//       //                               /(pclass_sz->n_m_L_sat-1.); //log(R)


//       double z =   exp(pclass_sz->array_z_L_sat[index_z])-1.;
//       double logM =   pclass_sz->array_m_L_sat[index_M];

//       double lnMs_min = log(pclass_sz->M_min_HOD_cib);
//       double lnMs_max = logM;//log(1e11);

//       double epsrel = pclass_sz->epsrel_L_sat;
//       double epsabs = pclass_sz->epsabs_L_sat;

//       struct Parameters_for_integrand_mass_L_sat V;
//       V.nu = pclass_sz->cib_frequency_list[index_nu];
//       V.z = z;
//       V.pclass_sz = pclass_sz;
//       V.M_host = exp(logM);


//       void * params = &V;
//       params = &V;

//       double L_sat_at_nu = Integrate_using_Patterson_adaptive(lnMs_min, lnMs_max,
//                                                                epsrel, epsabs,
//                                                                integrand_patterson_L_sat,
//                                                                params,pclass_sz->patterson_show_neval);

//       array_L_sat_at_z_and_M_at_nu[index_z][index_M] = log(1.+L_sat_at_nu);
//       index_z_M += 1;
//     }
//   }
// #ifdef _OPENMP
//   tstop = omp_get_wtime();
//   if (pclass_sz->sz_verbose > 0)
//     printf("In %s: time spent in parallel region (L_sat) at %.3e GHz = %e s for thread %d\n",
//            __func__,pclass_sz->cib_frequency_list[index_nu],tstop-tstart,omp_get_thread_num());
// #endif

//     }
// if (abort == _TRUE_) return _FAILURE_;
// //end of parallel region

// index_z_M = 0;
// for (index_M=0; index_M<pclass_sz->n_m_L_sat; index_M++)
// {
//   for (index_z=0; index_z<pclass_sz->n_z_L_sat; index_z++)
//   {
//     pclass_sz->array_L_sat_at_z_and_M_at_nu[index_nu][index_z_M] = array_L_sat_at_z_and_M_at_nu[index_z][index_M];
//     index_z_M += 1;
//   }
// }

// for (index_z=0;index_z<pclass_sz->n_z_L_sat;index_z++){
//   free(array_L_sat_at_z_and_M_at_nu[index_z]);
// }
//   free(array_L_sat_at_z_and_M_at_nu);


// }

// //exit(0);

// return _SUCCESS_;

//                                       }



//Tabulate Sigma2(R,z) and dSigma2dR
//as functions of z and logR
int tabulate_sigma_and_dsigma_from_pk(struct background * pba,
                                      struct nonlinear * pnl,
                                      struct primordial * ppm,
                                      struct class_sz_structure * pclass_sz){


if (pclass_sz->use_class_sz_fast_mode){

    if (pclass_sz->cosmo_model == 6)

        pclass_sz->ndim_masses = 1000;
        
    
    else
    
        pclass_sz->ndim_masses = 500;

    pclass_sz->n_m_dndlnM = pclass_sz->ndim_masses;
    
    }

if (pclass_sz->need_sigma == 0)
    return 0;



// class_alloc(pclass_sz->array_redshift,sizeof(double *)*pclass_sz->ndim_redshifts,pclass_sz->error_message);


if (pclass_sz->use_class_sz_fast_mode){
   // fixed by cosmopower emulator k-sampling
  class_alloc(pclass_sz->array_radius,sizeof(double *)*pclass_sz->ndim_masses,pclass_sz->error_message);


  class_alloc(pclass_sz->array_sigma_at_z_and_R,
              sizeof(double *)*pclass_sz->ndim_redshifts*pclass_sz->ndim_masses,
              pclass_sz->error_message);

  class_alloc( pclass_sz->array_dsigma2dR_at_z_and_R,
              sizeof(double *)*pclass_sz->ndim_redshifts*pclass_sz->ndim_masses,
              pclass_sz->error_message);


  return _SUCCESS_;
}
else{

class_alloc(pclass_sz->array_radius,sizeof(double *)*pclass_sz->ndim_masses,pclass_sz->error_message);


class_alloc(pclass_sz->array_sigma_at_z_and_R,
            sizeof(double *)*pclass_sz->ndim_redshifts*pclass_sz->ndim_masses,
            pclass_sz->error_message);

class_alloc( pclass_sz->array_dsigma2dR_at_z_and_R,
            sizeof(double *)*pclass_sz->ndim_redshifts*pclass_sz->ndim_masses,
            pclass_sz->error_message);

// printf("tabulating sigma\n");

   // bounds array of radii for sigma computations:
   pclass_sz->logR1SZ = log(pow(3.*0.1*1e4/(4*_PI_*pclass_sz->Omega_m_0*pclass_sz->Rho_crit_0),1./3.));
   pclass_sz->logR2SZ = log(pow(3.*10.*1e17/(4*_PI_*pclass_sz->Omega_m_0*pclass_sz->Rho_crit_0),1./3.));


   // pclass_sz->logR1SZ = -5.684; // 0.0034Mpc/h, 1.8e4  solar mass
   // pclass_sz->logR2SZ = 4.; //default =4 , i.e., 54.9Mpc/h, 7.5e16 solar mass


  //Array of z
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,pclass_sz->z_for_pk_hm);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,pclass_sz->z_for_pk_hm);
  // z_min = 0.99*z_min;

  int index_z;

  double tstart, tstop;
  int index_l;
  // double * sigma_var;
  // double * dsigma_var;
  int abort;

  //Array of R in Mpc
  double logR_min = log(exp(pclass_sz->logR1SZ)/pba->h); //in Mpc
  double logR_max = log(exp(pclass_sz->logR2SZ)/pba->h); //in Mpc

  // Print R_min and R_max
  double R_min = exp(logR_min);
  double R_max = exp(logR_max);
  if (pclass_sz->sz_verbose>2){
  printf("R_min = %.6e Mpc\n", R_min);
  printf("R_max = %.6e Mpc\n", R_max);
  }
// exit(0);
  
  int index_R;

  int index_z_R = 0;

  double ** array_sigma_at_z_and_R;
  double ** array_dsigma2dR_at_z_and_R;


class_alloc(array_sigma_at_z_and_R,
            pclass_sz->ndim_redshifts*sizeof(double *),
            pclass_sz->error_message);

class_alloc(array_dsigma2dR_at_z_and_R,
            pclass_sz->ndim_redshifts*sizeof(double *),
            pclass_sz->error_message);

for (index_l=0;
     index_l<pclass_sz->ndim_redshifts;
     index_l++)
{
  class_alloc(array_sigma_at_z_and_R[index_l],
              pclass_sz->ndim_masses*sizeof(double),
              pclass_sz->error_message);

  class_alloc(array_dsigma2dR_at_z_and_R[index_l],
              pclass_sz->ndim_masses*sizeof(double),
              pclass_sz->error_message);
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
pba,pclass_sz,ppm,pnl,z_min,z_max,logR_min,logR_max)\
private(tstart, tstop,index_R,index_z) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

  // class_alloc_parallel(sigma_var,
  //                      sizeof(double *),
  //                      pclass_sz->error_message);
  //
  // class_alloc_parallel(dsigma_var,
  //                      sizeof(double *),
  //                      pclass_sz->error_message);


// #pragma omp for schedule (dynamic)
#pragma omp for collapse(2)
  for (index_R=0; index_R<pclass_sz->ndim_masses; index_R++) // ndim_masses
  {
  for (index_z=0; index_z<pclass_sz->ndim_redshifts; index_z++) // ndim_redshift
  {
// #pragma omp flush(abort)

    double sigma_var,dsigma_var;
    pclass_sz->array_radius[index_R] =
                                logR_min
                                +index_R*(logR_max-logR_min)
                                /(pclass_sz->ndim_masses-1.); //log(R)

    // for (index_z=0; index_z<pclass_sz->ndim_redshifts; index_z++)
    // {
      // pclass_sz->array_redshift[index_z] =
      //                                 log(1.+z_min)
      //                                 +index_z*(log(1.+z_max)-log(1.+z_min))
      //                                 /(pclass_sz->ndim_redshifts-1.); // log(1+z)


  if (pclass_sz->need_sigma == 1){
      if (pclass_sz->HMF_prescription_NCDM == 2) //No-pres
        spectra_sigma_for_tSZ(pba,
                              ppm,
                              pnl,
                              pclass_sz,
                              exp(pclass_sz->array_radius[index_R]),
                              exp(pclass_sz->array_redshift[index_z])-1.,
                              &sigma_var//&sigma_at_z_and_R
                              );
      else
        spectra_sigma_ncdm( pba,
                           // spectra_sigma_ncdm( pba,
                           ppm,
                           pnl,
                           pclass_sz,
                           exp(pclass_sz->array_radius[index_R]),
                           exp(pclass_sz->array_redshift[index_z])-1.,
                           &sigma_var//&sigma_at_z_and_R
                           );



      //pclass_sz->array_sigma_at_z_and_R[index_z_R] = log(*sigma_var);//sigma_at_z_and_R); //log(sigma)
      array_sigma_at_z_and_R[index_z][index_R] = log(sigma_var);//sigma_at_z_and_R); //log(sigma)
      // printf("s=%.6e r=%.6e z=%.6e\n",
      // array_sigma_at_z_and_R[index_z][index_R],
      // exp(pclass_sz->array_radius[index_R]),
      // exp(pclass_sz->array_redshift[index_z])-1.);
      // exit(0);


      if (pclass_sz->HMF_prescription_NCDM == 2) //No-pres
        spectra_sigma_prime(pba,
                            ppm,
                            pnl,
                            pclass_sz,
                            exp(pclass_sz->array_radius[index_R]),
                            exp(pclass_sz->array_redshift[index_z])-1.,
                            &dsigma_var//&dsigma2dR_at_z_and_R
                            );
      else
        spectra_sigma_ncdm_prime(pba,
                                 ppm,
                                 pnl,
                                 pclass_sz,
                                 exp(pclass_sz->array_radius[index_R]),
                                 exp(pclass_sz->array_redshift[index_z])-1.,
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
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over R's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    // free(sigma_var);
    // free(dsigma_var);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_R = 0;
for (index_R=0; index_R<pclass_sz->ndim_masses; index_R++)
{
  for (index_z=0; index_z<pclass_sz->ndim_redshifts; index_z++)
  {

    pclass_sz->array_sigma_at_z_and_R[index_z_R] = array_sigma_at_z_and_R[index_z][index_R];
    double sigma =   pclass_sz->array_sigma_at_z_and_R[index_z_R];

    // printf("z = %.3e sig = %.3e\n",
    // exp(pclass_sz->array_redshift[index_z])-1.,
    // pclass_sz->array_sigma_at_z_and_R[index_z_R]);
    pclass_sz->array_dsigma2dR_at_z_and_R[index_z_R]=array_dsigma2dR_at_z_and_R[index_z][index_R];
    double dsigma = pclass_sz->array_dsigma2dR_at_z_and_R[index_z_R];
    if (isnan(sigma+dsigma) || isinf(sigma+dsigma) ){
     printf("z=%.3e R=%.3e sigma=%.3e dsigma2dR=%.3e\n",exp(pclass_sz->array_redshift[index_z]),exp(pclass_sz->array_radius[index_R]),sigma,dsigma);
     exit(0);
   }
    index_z_R += 1;
  }
}



// freeing memory
for (index_l=0;
     index_l<pclass_sz->ndim_redshifts;
     index_l++)
{
  free(array_sigma_at_z_and_R[index_l]);
  free(array_dsigma2dR_at_z_and_R[index_l]);
}


  free(array_sigma_at_z_and_R);
  free(array_dsigma2dR_at_z_and_R);

return _SUCCESS_;

}
}


// Tabulate redshift_int_lensmag
// as functions of z
int tabulate_redshift_int_lensmag(struct class_sz_structure * pclass_sz,
                                  struct background * pba){

if (pclass_sz->has_kSZ_kSZ_lensmag_1halo
  + pclass_sz->has_gal_lensmag_1h
  + pclass_sz->has_gal_lensmag_2h
  + pclass_sz->has_gal_lensmag_hf
  + pclass_sz->has_gallens_lensmag_1h
  + pclass_sz->has_gallens_lensmag_2h
  + pclass_sz->has_lens_lensmag_1h
  + pclass_sz->has_lens_lensmag_2h
  + pclass_sz->has_lens_lensmag_hf
  + pclass_sz->has_tSZ_lensmag_1h
  + pclass_sz->has_tSZ_lensmag_2h
  + pclass_sz->has_lensmag_lensmag_1h
  + pclass_sz->has_lensmag_lensmag_2h
  + pclass_sz->has_lensmag_lensmag_hf
   == _FALSE_){
    // if (pclass_sz->sz_verbose>=1)
    // printf("-> Not tabulating Wz for lensing magnification\n");
    return 0;
  }

if (pclass_sz->sz_verbose>=1){
printf("-> Tabulating Wz for lensing magnification\n");
}
  //Array of z
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  int index_z;
  double ln1pz,z;
  class_alloc(pclass_sz->array_W_lensmag,sizeof(double *)*pclass_sz->n_z_W_lensmag,pclass_sz->error_message);
  class_alloc(pclass_sz->array_z_W_lensmag,sizeof(double *)*pclass_sz->n_z_W_lensmag,pclass_sz->error_message);

  double * pvectsz;
  class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);

  double * pvecback;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  for (index_z=0; index_z<pclass_sz->n_z_W_lensmag; index_z++)
  {
    ln1pz =  log(1.+z_min)
              +index_z*(log(1.+z_max)-log(1.+z_min))
              /(pclass_sz->n_z_W_lensmag-1.); // log(1+z)

    z = exp(ln1pz) - 1.;

    // set redshift z
    pvectsz[pclass_sz->index_z] = z;

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
    pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);

    // printf("-> doing tabulating Wz for lensing magnification\n");

    // printf("-> Computing integral at z=%.3e\n",z);
    double result;
    redshift_int_lensmag(pclass_sz,pba,pvectsz,&result);

    if (result <= 0.)
      result = 1e-100;
    pclass_sz->array_W_lensmag[index_z] = log(result);
    pclass_sz->array_z_W_lensmag[index_z] = ln1pz;
    // printf("-> integral z = %.3e W = =%.3e\n",z,exp(pclass_sz->array_W_lensmag[index_z]));
}
if (pclass_sz->sz_verbose>=1)
printf("-> end tabulating Wz for lensing magnification\n");
 free(pvectsz);
 free(pvecback);
}

int evaluate_redshift_int_lensmag(double * pvectsz,
                                  struct class_sz_structure * pclass_sz)
  {

   double z = pvectsz[pclass_sz->index_z];
   double z_asked = log(1.+z);

   if (z<exp(pclass_sz->array_z_W_lensmag[0])-1.)
      z_asked = pclass_sz->array_z_W_lensmag[0];
   if (z>exp(pclass_sz->array_z_W_lensmag[pclass_sz->n_z_W_lensmag-1])-1.)
      z_asked =  pclass_sz->array_z_W_lensmag[pclass_sz->n_z_W_lensmag-1];


   pvectsz[pclass_sz->index_W_lensmag] =  exp(pwl_value_1d(pclass_sz->n_z_W_lensmag,
                                                        pclass_sz->array_z_W_lensmag,
                                                        pclass_sz->array_W_lensmag,
                                                        z_asked));

return _SUCCESS_;
}

int tabulate_redshift_int_nlensmag(struct class_sz_structure * pclass_sz,
                                  struct background * pba){

if (pclass_sz->has_nlensmag_gallens_1h
  + pclass_sz->has_nlensmag_gallens_2h
  + pclass_sz->has_nlensmag_tsz_1h
  + pclass_sz->has_nlensmag_tsz_2h
   == _FALSE_){
    // if (pclass_sz->sz_verbose>=1)
    // printf("-> Not tabulating Wz for lensing magnification\n");
    return 0;
  }

if (pclass_sz->sz_verbose>=1){
printf("-> [nlensmag] Tabulating Wz for n lensing magnification\n");
}
  //Array of z
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  int index_z;
  double ln1pz,z;

  class_alloc(pclass_sz->array_W_nlensmag,sizeof(double **)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);
  class_alloc(pclass_sz->array_z_W_nlensmag,sizeof(double **)*pclass_sz->galaxy_samples_list_num,pclass_sz->error_message);

  int index_g;
  for (index_g=0;index_g<pclass_sz->galaxy_samples_list_num;index_g++){
    class_alloc(pclass_sz->array_W_nlensmag[index_g],sizeof(double *)*pclass_sz->n_z_W_lensmag,pclass_sz->error_message);
    class_alloc(pclass_sz->array_z_W_nlensmag[index_g],sizeof(double *)*pclass_sz->n_z_W_lensmag,pclass_sz->error_message);
    }


for (index_g=0;index_g<pclass_sz->galaxy_samples_list_num;index_g++){

  if (pclass_sz->sz_verbose>0)
  printf("-> [nlensmag] starting computation for sample %d.\n",index_g);

  double * pvectsz;
  double * pvecback;
  class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  for (index_z=0; index_z<pclass_sz->n_z_W_lensmag; index_z++) //ola2 n_z_W_nlensmag
  {

    ln1pz =  log(1.+z_min)
              +index_z*(log(1.+z_max)-log(1.+z_min))
              /(pclass_sz->n_z_W_lensmag-1.); // log(1+z)

    z = exp(ln1pz) - 1.;

    // set redshift z
    pvectsz[pclass_sz->index_z] = z;

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
    pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);


    //printf("-> Computing integral at z=%.3e\n",z);
    double result;
    ///
    redshift_int_nlensmag(index_g,pclass_sz,pba,pvectsz,&result);

    if (result <= 0.)
      result = 1e-100;

    pclass_sz->array_W_nlensmag[index_g][index_z] = log(result);
    pclass_sz->array_z_W_nlensmag[index_g][index_z] = ln1pz;
//printf("-> integral z = %.3e W = %.3e\n",z,exp(pclass_sz->array_W_nlensmag[index_g][index_z]));
  }
if (pclass_sz->sz_verbose>=1)
printf("-> [nlensmag] end tabulating Wz for lensing magnification\n");
 free(pvectsz);
 free(pvecback);
}
}

// int evaluate_redshift_int_nlensmag(double * pvectsz,
//                                   struct class_sz_structure * pclass_sz)
//   {
//
//    double z = pvectsz[pclass_sz->index_z];
//    double z_asked = log(1.+z);
//
//    if (z<exp(pclass_sz->array_z_W_nlensmag[0])-1.)
//       z_asked = pclass_sz->array_z_W_nlensmag[0];
//    if (z>exp(pclass_sz->array_z_W_nlensmag[pclass_sz->n_z_W_lensmag-1])-1.)
//       z_asked =  pclass_sz->array_z_W_nlensmag[pclass_sz->n_z_W_lensmag-1];
//
//
//    pvectsz[pclass_sz->index_W_nlensmag] =  exp(pwl_value_1d(pclass_sz->n_z_W_lensmag,
//                                                         pclass_sz->array_z_W_nlensmag,
//                                                         pclass_sz->array_W_nlensmag,
//                                                         z_asked));
//
// return _SUCCESS_;
// }


// Tabulate redshift_int_lensmag
// as functions of z
int tabulate_redshift_int_gallens_sources(struct class_sz_structure * pclass_sz,
                                          struct background * pba){

if (
    pclass_sz->has_gal_gallens_2h
  + pclass_sz->has_gal_gallens_1h
  + pclass_sz->has_ngal_gallens_2h
  + pclass_sz->has_ngal_gallens_1h
  + pclass_sz->has_nlensmag_gallens_2h
  + pclass_sz->has_nlensmag_gallens_1h
  + pclass_sz->has_tSZ_gallens_2h
  + pclass_sz->has_tSZ_gallens_1h
  + pclass_sz->has_gallens_gallens_2h
  + pclass_sz->has_gallens_gallens_1h
  + pclass_sz->has_gallens_cib_2h
  + pclass_sz->has_gallens_cib_1h
  + pclass_sz->has_gallens_lens_2h
  + pclass_sz->has_gallens_lens_1h
  + pclass_sz->has_gallens_lensmag_2h
  + pclass_sz->has_gallens_lensmag_1h
  + pclass_sz->has_kSZ_kSZ_gallens_1h_fft
  + pclass_sz->has_kSZ_kSZ_gallens_2h_fft
  + pclass_sz->has_kSZ_kSZ_gallens_3h_fft
  + pclass_sz->has_kSZ_kSZ_gallens_hf
   == _FALSE_){
    // if (pclass_sz->sz_verbose>=1)
    // printf("-> Not tabulating Wz for lensing magnification\n");
    return 0;
  }

if (pclass_sz->sz_verbose>=1){
printf("-> Tabulating Wz for source galaxies\n");
}
  //Array of z
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  int index_z;
  double ln1pz,z;
  class_alloc(pclass_sz->array_W_gallens_sources,sizeof(double *)*pclass_sz->n_z_W_gallens_sources,pclass_sz->error_message);
  class_alloc(pclass_sz->array_z_W_gallens_sources,sizeof(double *)*pclass_sz->n_z_W_gallens_sources,pclass_sz->error_message);

  double * pvectsz;
  class_alloc(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);

  double * pvecback;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
// printf("-> nz=%d\n",pclass_sz->n_z_W_gallens_sources);
  for (index_z=0; index_z<pclass_sz->n_z_W_gallens_sources; index_z++)
  {
    ln1pz =  log(1.+z_min)
              +index_z*(log(1.+z_max)-log(1.+z_min))
              /(pclass_sz->n_z_W_gallens_sources-1.); // log(1+z)

    z = exp(ln1pz) - 1.;

    // set redshift z
    pvectsz[pclass_sz->index_z] = z;

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
    pvectsz[pclass_sz->index_chi2] = pow(pvecback[pba->index_bg_ang_distance]*(1.+z)*pba->h,2);

    // printf("-> doing tabulating Wz for lensing magnification\n");

    // printf("-> Computing integral at z=%.3e\n",z);
    double result;
    redshift_int_gallens_sources(pclass_sz,pba,pvectsz,&result);
      // printf("-> 2 doing tabulating Wz for lensing magnification\n");
    if (result <= 0.)
      result = 1e-100;
    pclass_sz->array_W_gallens_sources[index_z] = log(result);
    pclass_sz->array_z_W_gallens_sources[index_z] = ln1pz;
    // printf("-> integral z = %.3e W = %.3e\n",z,exp(pclass_sz->array_W_gallens_sources[index_z]));
}
if (pclass_sz->sz_verbose>=1)
printf("-> end tabulating Wz for source galaxies\n");
 free(pvectsz);
 free(pvecback);
}



double  evaluate_redshift_int_gallens_sources(double z,
                                              struct class_sz_structure * pclass_sz)
  {

   // double z = pvectsz[pclass_sz->index_z];
   double z_asked = log(1.+z);

   if (z<exp(pclass_sz->array_z_W_gallens_sources[0])-1.)
      z_asked = pclass_sz->array_z_W_gallens_sources[0];
   if (z>exp(pclass_sz->array_z_W_gallens_sources[pclass_sz->n_z_W_gallens_sources-1])-1.)
      z_asked =  pclass_sz->array_z_W_gallens_sources[pclass_sz->n_z_W_gallens_sources-1];


   // pvectsz[pclass_sz->index_W_gallens_sources]
   double result  =  exp(pwl_value_1d(pclass_sz->n_z_W_gallens_sources,
                                      pclass_sz->array_z_W_gallens_sources,
                                      pclass_sz->array_W_gallens_sources,
                                      z_asked));
// if ( pvectsz[pclass_sz->index_W_gallens_sources] == 0){
if ( result == 0){
  printf("null W gallens source %.3e\n",z);
  exit(0);
}
return result;
}




// Tabulate dndlnM
// as functions of z and M
int tabulate_dndlnM(struct background * pba,
                    struct nonlinear * pnl,
                    struct primordial * ppm,
                    struct class_sz_structure * pclass_sz){

  //Array of z
  // double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,pclass_sz->z_for_pk_hm);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,pclass_sz->z_for_pk_hm);

  if (pclass_sz->sz_verbose>10){
  printf("tabulating dndlnM between z_min=%.5e and z_max=%.5e\n",z_min,z_max);
}
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(pclass_sz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(pclass_sz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_dndlnM_at_z_and_M;

  class_alloc(pclass_sz->array_z_dndlnM,sizeof(double *)*pclass_sz->n_z_dndlnM,pclass_sz->error_message);
  class_alloc(pclass_sz->array_m_dndlnM,sizeof(double *)*pclass_sz->n_m_dndlnM,pclass_sz->error_message);


class_alloc(pclass_sz->array_dndlnM_at_z_and_M,
            sizeof(double *)*pclass_sz->n_z_dndlnM*pclass_sz->n_m_dndlnM,
            pclass_sz->error_message);


class_alloc(array_dndlnM_at_z_and_M,
            pclass_sz->n_z_dndlnM*sizeof(double *),
            pclass_sz->error_message);


for (index_l=0;
     index_l<pclass_sz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_dndlnM_at_z_and_M[index_l],
              pclass_sz->n_m_dndlnM*sizeof(double),
              pclass_sz->error_message);
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
pba,pclass_sz,ppm,pnl,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_z,index_M,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);

class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);


int i;
for(i = 0; i<pclass_sz->tsz_size;i++) pvectsz[i] = 0.;

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{
#pragma omp flush(abort)

      double tau;
      int first_index_back = 0;
      pclass_sz->array_z_dndlnM[index_z] =
                                log(1.+z_min)
                                +index_z*(log(1.+z_max)-log(1.+z_min))
                                /(pclass_sz->n_z_dndlnM-1.); // log(1+z)
      double z =   exp(pclass_sz->array_z_dndlnM[index_z])-1.;

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





for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{


      pclass_sz->array_m_dndlnM[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(pclass_sz->n_m_dndlnM-1.); //log(R)

      //background quantities @ z:

      double logM =   pclass_sz->array_m_dndlnM[index_M];


      pvectsz[pclass_sz->index_z] = z;
      pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);

      double omega = pvecback[pba->index_bg_Omega_m];
      pvectsz[pclass_sz->index_Delta_c]= Delta_c_of_Omega_m(omega);
      evaluate_HMF_at_logM_and_z(logM,z,pvecback,pvectsz,pba,pnl,pclass_sz);
      array_dndlnM_at_z_and_M[index_z][index_M] = log(pvectsz[pclass_sz->index_hmf]);
    }

  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
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
for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
  {
    pclass_sz->array_dndlnM_at_z_and_M[index_z_M] = array_dndlnM_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}


// freeing memory:
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++){
  free(array_dndlnM_at_z_and_M[index_z]);
}
free(array_dndlnM_at_z_and_M);

return _SUCCESS_;
}



///Tabulate m200c_to_m500c conversion
//as functions of z and M
int tabulate_m200c_to_m500c(struct background * pba,
                            struct class_sz_structure * pclass_sz){

  //Array of z
  // double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,pclass_sz->z_for_pk_hm);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,pclass_sz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(pclass_sz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(pclass_sz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_m200c_to_m500c_at_z_and_M;

  class_alloc(pclass_sz->array_ln_1pz_m200c_to_m500c,sizeof(double *)*pclass_sz->n_z_dndlnM,pclass_sz->error_message);
  class_alloc(pclass_sz->array_m_m200c_to_m500c,sizeof(double *)*pclass_sz->n_m_dndlnM,pclass_sz->error_message);


class_alloc(pclass_sz->array_m200c_to_m500c_at_z_and_M,
            sizeof(double *)*pclass_sz->n_z_dndlnM*pclass_sz->n_m_dndlnM,
            pclass_sz->error_message);


class_alloc(array_m200c_to_m500c_at_z_and_M,
            pclass_sz->n_z_dndlnM*sizeof(double *),
            pclass_sz->error_message);


for (index_l=0;
     index_l<pclass_sz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_m200c_to_m500c_at_z_and_M[index_l],
              pclass_sz->n_m_dndlnM*sizeof(double),
              pclass_sz->error_message);
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
pba,pclass_sz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
      pclass_sz->array_ln_1pz_m200c_to_m500c[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_dndlnM-1.);

      pclass_sz->array_m_m200c_to_m500c[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(pclass_sz->n_m_dndlnM-1.);

      //background quantities @ z:
      double z =   exp(pclass_sz->array_ln_1pz_m200c_to_m500c[index_z])-1.;
      double logM =   pclass_sz->array_m_m200c_to_m500c[index_M];
      pvectsz[pclass_sz->index_m200c] = exp(logM);
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




      pvectsz[pclass_sz->index_z] = z;
      pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[pclass_sz->index_Rho_crit];
    double delrho = 200.*rhoc; // 200m
    double delrho_prime = 500.*rhoc; //500c
    double mdel = pvectsz[pclass_sz->index_m200c];

    double mdel_prime;
    class_call_parallel(mDEL_to_mDELprime(mdel,
                           delrho,
                           delrho_prime,
                           delc,
                           rhoc,
                           z,
                           &mdel_prime,
                           pclass_sz,
                           pba),
                    pclass_sz->error_message,
                    pclass_sz->error_message);
    pvectsz[pclass_sz->index_m500c] = mdel_prime;

    array_m200c_to_m500c_at_z_and_M[index_z][index_M] = log(pvectsz[pclass_sz->index_m500c]);

    index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
  {
    pclass_sz->array_m200c_to_m500c_at_z_and_M[index_z_M] = array_m200c_to_m500c_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}

// freeing memory
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{
free(array_m200c_to_m500c_at_z_and_M[index_z]);
}
free(array_m200c_to_m500c_at_z_and_M);

return _SUCCESS_;
}



///Tabulate m500c_to_m200c conversion
//as functions of z and M
int tabulate_m500c_to_m200c(struct background * pba,
                            struct class_sz_structure * pclass_sz){

  //Array of z
  // double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,pclass_sz->z_for_pk_hm);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,pclass_sz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(pclass_sz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(pclass_sz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_m200c_to_m500c_at_z_and_M;

  class_alloc(pclass_sz->array_ln_1pz_m500c_to_m200c,sizeof(double *)*pclass_sz->n_z_dndlnM,pclass_sz->error_message);
  class_alloc(pclass_sz->array_m_m500c_to_m200c,sizeof(double *)*pclass_sz->n_m_dndlnM,pclass_sz->error_message);


class_alloc(pclass_sz->array_m500c_to_m200c_at_z_and_M,
            sizeof(double *)*pclass_sz->n_z_dndlnM*pclass_sz->n_m_dndlnM,
            pclass_sz->error_message);


class_alloc(array_m200c_to_m500c_at_z_and_M,
            pclass_sz->n_z_dndlnM*sizeof(double *),
            pclass_sz->error_message);


for (index_l=0;
     index_l<pclass_sz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_m200c_to_m500c_at_z_and_M[index_l],
              pclass_sz->n_m_dndlnM*sizeof(double),
              pclass_sz->error_message);
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
pba,pclass_sz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
      pclass_sz->array_ln_1pz_m500c_to_m200c[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_dndlnM-1.); // log(1+z)

      pclass_sz->array_m_m500c_to_m200c[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(pclass_sz->n_m_dndlnM-1.); //log(R)

      //background quantities @ z:
      double z =   exp(pclass_sz->array_ln_1pz_m500c_to_m200c[index_z])-1.;
      double logM =   pclass_sz->array_m_m500c_to_m200c[index_M];
      pvectsz[pclass_sz->index_m500c] = exp(logM);
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




      pvectsz[pclass_sz->index_z] = z;
      pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[pclass_sz->index_Rho_crit];
    double delrho = 500.*rhoc; // 200m
    double delrho_prime = 200.*rhoc; //500c
    double mdel = pvectsz[pclass_sz->index_m500c];

    double mdel_prime;
    class_call_parallel(mDEL_to_mDELprime(mdel,
                           delrho,
                           delrho_prime,
                           delc,
                           rhoc,
                           z,
                           &mdel_prime,
                           pclass_sz,
                           pba),
                    pclass_sz->error_message,
                    pclass_sz->error_message);
    pvectsz[pclass_sz->index_m200c] = mdel_prime;

    array_m200c_to_m500c_at_z_and_M[index_z][index_M] = log(pvectsz[pclass_sz->index_m200c]);
    // printf("m = %.3e\n",array_m200m_to_m500c_at_z_and_M[index_z][index_M]);

    index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
  {
    pclass_sz->array_m500c_to_m200c_at_z_and_M[index_z_M] = array_m200c_to_m500c_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}

// freeing memory
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{
free(array_m200c_to_m500c_at_z_and_M[index_z]);
}
free(array_m200c_to_m500c_at_z_and_M);

return _SUCCESS_;
}

int tabulate_m500c_to_m200m(struct background * pba,
                            struct class_sz_structure * pclass_sz){

  //Array of z
  // double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,pclass_sz->z_for_pk_hm);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,pclass_sz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(pclass_sz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(pclass_sz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_m200m_to_m500c_at_z_and_M;

  class_alloc(pclass_sz->array_ln_1pz_m500c_to_m200m,sizeof(double *)*pclass_sz->n_z_dndlnM,pclass_sz->error_message);
  class_alloc(pclass_sz->array_m_m500c_to_m200m,sizeof(double *)*pclass_sz->n_m_dndlnM,pclass_sz->error_message);


class_alloc(pclass_sz->array_m500c_to_m200m_at_z_and_M,
            sizeof(double *)*pclass_sz->n_z_dndlnM*pclass_sz->n_m_dndlnM,
            pclass_sz->error_message);


class_alloc(array_m200m_to_m500c_at_z_and_M,
            pclass_sz->n_z_dndlnM*sizeof(double *),
            pclass_sz->error_message);


for (index_l=0;
     index_l<pclass_sz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_m200m_to_m500c_at_z_and_M[index_l],
              pclass_sz->n_m_dndlnM*sizeof(double),
              pclass_sz->error_message);
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
pba,pclass_sz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
      pclass_sz->array_ln_1pz_m500c_to_m200m[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_dndlnM-1.); // log(1+z)

      pclass_sz->array_m_m500c_to_m200m[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(pclass_sz->n_m_dndlnM-1.); //log(R)

      //background quantities @ z:
      double z =   exp(pclass_sz->array_ln_1pz_m500c_to_m200m[index_z])-1.;
      double logM =   pclass_sz->array_m_m500c_to_m200m[index_M];
      pvectsz[pclass_sz->index_m500c] = exp(logM);
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




      pvectsz[pclass_sz->index_z] = z;
      pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[pclass_sz->index_Rho_crit];
    double delrho = 500.*rhoc; // 200m
    double delrho_prime = 200.*rhoc; //500c
    double mdel = pvectsz[pclass_sz->index_m500c];

    double mdel_prime;
    class_call_parallel(mDEL_to_mDELprime(mdel,
                           delrho,
                           delrho_prime,
                           delc,
                           rhoc,
                           z,
                           &mdel_prime,
                           pclass_sz,
                           pba),
                    pclass_sz->error_message,
                    pclass_sz->error_message);
    pvectsz[pclass_sz->index_m200m] = mdel_prime;

    array_m200m_to_m500c_at_z_and_M[index_z][index_M] = log(pvectsz[pclass_sz->index_m200m]);
    // printf("m = %.3e\n",array_m200m_to_m500c_at_z_and_M[index_z][index_M]);

    index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
  {
    pclass_sz->array_m500c_to_m200m_at_z_and_M[index_z_M] = array_m200m_to_m500c_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}

// freeing memory
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{
free(array_m200m_to_m500c_at_z_and_M[index_z]);
}
free(array_m200m_to_m500c_at_z_and_M);

return _SUCCESS_;
}




///Tabulate m200m_to_m500c conversion
//as functions of z and M
int tabulate_m200m_to_m500c(struct background * pba,
                            struct class_sz_structure * pclass_sz){

  //Array of z
  // double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,pclass_sz->z_for_pk_hm);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,pclass_sz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(pclass_sz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(pclass_sz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_m200m_to_m500c_at_z_and_M;

  class_alloc(pclass_sz->array_ln_1pz_m200m_to_m500c,sizeof(double *)*pclass_sz->n_z_dndlnM,pclass_sz->error_message);
  class_alloc(pclass_sz->array_m_m200m_to_m500c,sizeof(double *)*pclass_sz->n_m_dndlnM,pclass_sz->error_message);


class_alloc(pclass_sz->array_m200m_to_m500c_at_z_and_M,
            sizeof(double *)*pclass_sz->n_z_dndlnM*pclass_sz->n_m_dndlnM,
            pclass_sz->error_message);


class_alloc(array_m200m_to_m500c_at_z_and_M,
            pclass_sz->n_z_dndlnM*sizeof(double *),
            pclass_sz->error_message);


for (index_l=0;
     index_l<pclass_sz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_m200m_to_m500c_at_z_and_M[index_l],
              pclass_sz->n_m_dndlnM*sizeof(double),
              pclass_sz->error_message);
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
pba,pclass_sz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
      pclass_sz->array_ln_1pz_m200m_to_m500c[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_dndlnM-1.); // log(1+z)

      pclass_sz->array_m_m200m_to_m500c[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(pclass_sz->n_m_dndlnM-1.); //log(R)

      //background quantities @ z:
      double z =   exp(pclass_sz->array_ln_1pz_m200m_to_m500c[index_z])-1.;
      double logM =   pclass_sz->array_m_m200m_to_m500c[index_M];
      pvectsz[pclass_sz->index_m200m] = exp(logM);
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




      pvectsz[pclass_sz->index_z] = z;
      pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[pclass_sz->index_Rho_crit];
    double delrho = 200.*omega*rhoc; // 200m
    double delrho_prime = 500.*rhoc; //500c
    double mdel = pvectsz[pclass_sz->index_m200m];

    double mdel_prime;
    class_call_parallel(mDEL_to_mDELprime(mdel,
                           delrho,
                           delrho_prime,
                           delc,
                           rhoc,
                           z,
                           &mdel_prime,
                           pclass_sz,
                           pba),
                    pclass_sz->error_message,
                    pclass_sz->error_message);
    pvectsz[pclass_sz->index_m500c] = mdel_prime;

    array_m200m_to_m500c_at_z_and_M[index_z][index_M] = log(pvectsz[pclass_sz->index_m500c]);
    // printf("m = %.3e\n",array_m200m_to_m500c_at_z_and_M[index_z][index_M]);

    index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
  {
    pclass_sz->array_m200m_to_m500c_at_z_and_M[index_z_M] = array_m200m_to_m500c_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}

// freeing memory
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{
free(array_m200m_to_m500c_at_z_and_M[index_z]);
}
free(array_m200m_to_m500c_at_z_and_M);

return _SUCCESS_;
}


///Tabulate m200m_to_m200c conversion
//as functions of z and M
int tabulate_m200m_to_m200c(struct background * pba,
                            struct class_sz_structure * pclass_sz){

  //Array of z
  // double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,pclass_sz->z_for_pk_hm);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,pclass_sz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(pclass_sz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(pclass_sz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_m200m_to_m200c_at_z_and_M;

  class_alloc(pclass_sz->array_ln_1pz_m200m_to_m200c,sizeof(double *)*pclass_sz->n_z_dndlnM,pclass_sz->error_message);
  class_alloc(pclass_sz->array_m_m200m_to_m200c,sizeof(double *)*pclass_sz->n_m_dndlnM,pclass_sz->error_message);


class_alloc(pclass_sz->array_m200m_to_m200c_at_z_and_M,
            sizeof(double *)*pclass_sz->n_z_dndlnM*pclass_sz->n_m_dndlnM,
            pclass_sz->error_message);


class_alloc(array_m200m_to_m200c_at_z_and_M,
            pclass_sz->n_z_dndlnM*sizeof(double *),
            pclass_sz->error_message);


for (index_l=0;
     index_l<pclass_sz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_m200m_to_m200c_at_z_and_M[index_l],
              pclass_sz->n_m_dndlnM*sizeof(double),
              pclass_sz->error_message);
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
pba,pclass_sz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
      pclass_sz->array_ln_1pz_m200m_to_m200c[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_dndlnM-1.); // log(1+z)

      pclass_sz->array_m_m200m_to_m200c[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(pclass_sz->n_m_dndlnM-1.); //log(R)

      //background quantities @ z:
      double z =   exp(pclass_sz->array_ln_1pz_m200m_to_m200c[index_z])-1.;
      double logM =   pclass_sz->array_m_m200m_to_m200c[index_M];
      pvectsz[pclass_sz->index_m200m] = exp(logM);
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




      pvectsz[pclass_sz->index_z] = z;
      pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[pclass_sz->index_Rho_crit];
    double delrho = 200.*omega*rhoc; // 200m
    double delrho_prime = 200.*rhoc; //500c
    double mdel = pvectsz[pclass_sz->index_m200m];

    double mdel_prime;


    if (pclass_sz->use_websky_m200m_to_m200c_conversion == 1){
      // omegamz = co.omegam*(1+z)**3/(co.omegam*(1+z)**3+1-co.omegam)
      // m200c   = omegamz**0.35 * m200m # m200m to m200c conversion used for websky
      // return m200c

      pvectsz[pclass_sz->index_m200c] = pow(omega,0.35)*pvectsz[pclass_sz->index_m200m];
    }
    else{
      class_call_parallel(mDEL_to_mDELprime(mdel,
                             delrho,
                             delrho_prime,
                             delc,
                             rhoc,
                             z,
                             &mdel_prime,
                             pclass_sz,
                             pba),
                      pclass_sz->error_message,
                      pclass_sz->error_message);
      pvectsz[pclass_sz->index_m200c] = mdel_prime;
    }

    array_m200m_to_m200c_at_z_and_M[index_z][index_M] = log(pvectsz[pclass_sz->index_m200c]);
    // printf("m = %.3e\n",array_m200m_to_m200c_at_z_and_M[index_z][index_M]);

    index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
  {
    pclass_sz->array_m200m_to_m200c_at_z_and_M[index_z_M] = array_m200m_to_m200c_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}

// freeing memory
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{
free(array_m200m_to_m200c_at_z_and_M[index_z]);
}
free(array_m200m_to_m200c_at_z_and_M);

return _SUCCESS_;
}






///Tabulate m200m_to_mvir conversion
//as functions of z and M
int tabulate_m200m_to_mvir(struct background * pba,
                            struct class_sz_structure * pclass_sz){

  //Array of z
  // double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,pclass_sz->z_for_pk_hm);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,pclass_sz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(pclass_sz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(pclass_sz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_m200m_to_mvir_at_z_and_M;

  class_alloc(pclass_sz->array_ln_1pz_m200m_to_mvir,sizeof(double *)*pclass_sz->n_z_dndlnM,pclass_sz->error_message);
  class_alloc(pclass_sz->array_m_m200m_to_mvir,sizeof(double *)*pclass_sz->n_m_dndlnM,pclass_sz->error_message);


class_alloc(pclass_sz->array_m200m_to_mvir_at_z_and_M,
            sizeof(double *)*pclass_sz->n_z_dndlnM*pclass_sz->n_m_dndlnM,
            pclass_sz->error_message);


class_alloc(array_m200m_to_mvir_at_z_and_M,
            pclass_sz->n_z_dndlnM*sizeof(double *),
            pclass_sz->error_message);


for (index_l=0;
     index_l<pclass_sz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_m200m_to_mvir_at_z_and_M[index_l],
              pclass_sz->n_m_dndlnM*sizeof(double),
              pclass_sz->error_message);
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
pba,pclass_sz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
      pclass_sz->array_ln_1pz_m200m_to_mvir[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_dndlnM-1.); // log(1+z)

      pclass_sz->array_m_m200m_to_mvir[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(pclass_sz->n_m_dndlnM-1.); //log(R)

      //background quantities @ z:
      double z =   exp(pclass_sz->array_ln_1pz_m200m_to_mvir[index_z])-1.;
      double logM =   pclass_sz->array_m_m200m_to_mvir[index_M];
      pvectsz[pclass_sz->index_m200m] = exp(logM);
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




      pvectsz[pclass_sz->index_z] = z;
      pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[pclass_sz->index_Rho_crit];
    double delrho = 200.*omega*rhoc; // 200m
    double delrho_prime = delc; //vir
    double mdel = pvectsz[pclass_sz->index_m200m];

    double mdel_prime;


      class_call_parallel(mDEL_to_mVIR(mdel,
                               delrho,
                               delc,
                               pvectsz[pclass_sz->index_Rho_crit],
                               z,
                               &mdel_prime,
                               pclass_sz,
                               pba),
                      pclass_sz->error_message,
                      pclass_sz->error_message);
      pvectsz[pclass_sz->index_mVIR] = mdel_prime;

    array_m200m_to_mvir_at_z_and_M[index_z][index_M] = log(pvectsz[pclass_sz->index_mVIR]);
    // printf("m = %.3e\n",array_m200m_to_m200c_at_z_and_M[index_z][index_M]);

    index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
  {
    pclass_sz->array_m200m_to_mvir_at_z_and_M[index_z_M] = array_m200m_to_mvir_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}

// freeing memory
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{
free(array_m200m_to_mvir_at_z_and_M[index_z]);
}
free(array_m200m_to_mvir_at_z_and_M);

return _SUCCESS_;
}




///Tabulate m200c_to_mvir conversion
//as functions of z and M
int tabulate_m200c_to_mvir(struct background * pba,
                            struct class_sz_structure * pclass_sz){

  //Array of z
  // double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,pclass_sz->z_for_pk_hm);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,pclass_sz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(pclass_sz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(pclass_sz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  int index_z_M = 0;

  double ** array_m200c_to_mvir_at_z_and_M;

  class_alloc(pclass_sz->array_ln_1pz_m200c_to_mvir,sizeof(double *)*pclass_sz->n_z_dndlnM,pclass_sz->error_message);
  class_alloc(pclass_sz->array_m_m200c_to_mvir,sizeof(double *)*pclass_sz->n_m_dndlnM,pclass_sz->error_message);


class_alloc(pclass_sz->array_m200c_to_mvir_at_z_and_M,
            sizeof(double *)*pclass_sz->n_z_dndlnM*pclass_sz->n_m_dndlnM,
            pclass_sz->error_message);


class_alloc(array_m200c_to_mvir_at_z_and_M,
            pclass_sz->n_z_dndlnM*sizeof(double *),
            pclass_sz->error_message);


for (index_l=0;
     index_l<pclass_sz->n_z_dndlnM;
     index_l++)
{
  class_alloc(array_m200c_to_mvir_at_z_and_M[index_l],
              pclass_sz->n_m_dndlnM*sizeof(double),
              pclass_sz->error_message);
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
pba,pclass_sz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{

#pragma omp flush(abort)

for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
      pclass_sz->array_ln_1pz_m200c_to_mvir[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_dndlnM-1.); // log(1+z)

      pclass_sz->array_m_m200c_to_mvir[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(pclass_sz->n_m_dndlnM-1.); //log(R)

      //background quantities @ z:
      double z =   exp(pclass_sz->array_ln_1pz_m200c_to_mvir[index_z])-1.;
      double logM =   pclass_sz->array_m_m200c_to_mvir[index_M];
      pvectsz[pclass_sz->index_m200c] = exp(logM);
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




      pvectsz[pclass_sz->index_z] = z;
      pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[pclass_sz->index_Rho_crit];
    double delrho = 200.*rhoc; // 200c
    double delrho_prime = delc; //vir
    double mdel = pvectsz[pclass_sz->index_m200c];

    double mdel_prime;


      class_call_parallel(mDEL_to_mVIR(mdel,
                               delrho,
                               delc,
                               pvectsz[pclass_sz->index_Rho_crit],
                               z,
                               &mdel_prime,
                               pclass_sz,
                               pba),
                      pclass_sz->error_message,
                      pclass_sz->error_message);
      pvectsz[pclass_sz->index_mVIR] = mdel_prime;

    array_m200c_to_mvir_at_z_and_M[index_z][index_M] = log(pvectsz[pclass_sz->index_mVIR]);
    // printf("m = %.3e\n",array_m200m_to_m200c_at_z_and_M[index_z][index_M]);

    index_z_M += 1;
    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    free(pvecback);
    free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_M = 0;
for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
  for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
  {
    pclass_sz->array_m200c_to_mvir_at_z_and_M[index_z_M] = array_m200c_to_mvir_at_z_and_M[index_z][index_M];
    index_z_M += 1;
  }
}

// freeing memory
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{
free(array_m200c_to_mvir_at_z_and_M[index_z]);
}
free(array_m200c_to_mvir_at_z_and_M);

return _SUCCESS_;
}



///Tabulate m200c_to_m200m conversion
//as functions of z and M
int tabulate_m200c_to_m200m(struct background * pba,
                            struct class_sz_structure * pclass_sz){

  //Array of z
  // double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  double z_min = r8_min(pclass_sz->z1SZ,pclass_sz->z1SZ_dndlnM);
  // z_min = r8_min(z_min,pclass_sz->z_for_pk_hm);
  double z_max = r8_max(pclass_sz->z2SZ,pclass_sz->z2SZ_dndlnM);
  // z_max = r8_min(z_max,pclass_sz->z_for_pk_hm);
  int index_z;

  double tstart, tstop;
  int index_l;

  double * pvecback;
  double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logM_min = log(pclass_sz->M1SZ_dndlnM); //in Msun/h
  double logM_max = log(pclass_sz->M2SZ_dndlnM); //in Msun/h
  int index_M;

  // int index_z_M = 0;

  // double ** array_m200c_to_m200m_at_z_and_M;

  class_alloc(pclass_sz->array_ln_1pz_m200c_to_m200m,sizeof(double *)*pclass_sz->n_z_dndlnM,pclass_sz->error_message);
  class_alloc(pclass_sz->array_m_m200c_to_m200m,sizeof(double *)*pclass_sz->n_m_dndlnM,pclass_sz->error_message);


class_alloc(pclass_sz->array_m200c_to_m200m_at_z_and_M,
            sizeof(double *)*pclass_sz->n_z_dndlnM*pclass_sz->n_m_dndlnM,
            pclass_sz->error_message);


// class_alloc(array_m200c_to_m200m_at_z_and_M,
//             pclass_sz->n_z_dndlnM*sizeof(double *),
//             pclass_sz->error_message);


// for (index_l=0;
//      index_l<pclass_sz->n_z_dndlnM;
//      index_l++)
// {
//   class_alloc(array_m200c_to_m200m_at_z_and_M[index_l],
//               pclass_sz->n_m_dndlnM*sizeof(double),
//               pclass_sz->error_message);
// }

for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{
      pclass_sz->array_ln_1pz_m200c_to_m200m[index_z] =
                                      log(1.+z_min)
                                      +index_z*(log(1.+z_max)-log(1.+z_min))
                                      /(pclass_sz->n_z_dndlnM-1.); // log(1+z)
}
for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{
      pclass_sz->array_m_m200c_to_m200m[index_M] =
                                    logM_min
                                    +index_M*(logM_max-logM_min)
                                    /(pclass_sz->n_m_dndlnM-1.); //log(R)
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
pba,pclass_sz,z_min,z_max,logM_min,logM_max)\
private(tstart, tstop,index_M,index_z,pvecback,pvectsz) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif


  class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);

  class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for collapse(2)
for (index_z=0; index_z<pclass_sz->n_z_dndlnM; index_z++)
{
for (index_M=0; index_M<pclass_sz->n_m_dndlnM; index_M++)
{

  int index_z_M = index_M * pclass_sz->n_z_dndlnM + index_z;
  // printf("%d\n",index_z_M);

      //background quantities @ z:
      double z =   exp(pclass_sz->array_ln_1pz_m200c_to_m200m[index_z])-1.;
      double logM =   pclass_sz->array_m_m200c_to_m200m[index_M];
      pvectsz[pclass_sz->index_m200c] = exp(logM);
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




      pvectsz[pclass_sz->index_z] = z;
      pvectsz[pclass_sz->index_Rho_crit] = (3./(8.*_PI_*_G_*_M_sun_))
                                      *pow(_Mpc_over_m_,1)
                                      *pow(_c_,2)
                                      *pvecback[pba->index_bg_rho_crit]
                                      /pow(pba->h,2);



    double omega = pvecback[pba->index_bg_Omega_m];///pow(Eh,2.);
    double delc = Delta_c_of_Omega_m(omega);
    double rhoc = pvectsz[pclass_sz->index_Rho_crit];
    double delrho = 200.*rhoc; // 200c
    double delrho_prime = 200.*omega*rhoc; //200m
    double mdel = pvectsz[pclass_sz->index_m200c];

    double mdel_prime;
    class_call_parallel(mDEL_to_mDELprime(mdel,
                           delrho,
                           delrho_prime,
                           delc,
                           rhoc,
                           z,
                           &mdel_prime,
                           pclass_sz,
                           pba),
                    pclass_sz->error_message,
                    pclass_sz->error_message);
    pvectsz[pclass_sz->index_m200m] = mdel_prime;
    // printf("%.8e\n",mdel_prime);

    pclass_sz->array_m200c_to_m200m_at_z_and_M[index_z_M] = log(pvectsz[pclass_sz->index_m200m]);


    }
  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
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
                       struct class_sz_structure * pclass_sz){

  //Array of z
  double z_min = pclass_sz->array_redshift[0];
  double z_max = pclass_sz->array_redshift[pclass_sz->ndim_redshifts-1];
  int index_z;

  double tstart, tstop;
  int index_l;

  // double * pvecback;
  // double * pvectsz;
  int abort;

  //Array of M in Msun/h
  double logk_min = pclass_sz->ln_k_for_tSZ[0]; //in Mpc/h
  double logk_max = pclass_sz->ln_k_for_tSZ[pclass_sz->ln_k_size_for_tSZ-1]; //in Mpc/h
  int index_k;

  int index_z_k = 0;

  double ** array_nl_index_at_z_and_k;
  double ** array_nl_index_at_z_and_k_no_wiggles;



class_alloc(pclass_sz->array_nl_index_at_z_and_k,
            sizeof(double *)*pclass_sz->ndim_redshifts*pclass_sz->ln_k_size_for_tSZ,
            pclass_sz->error_message);
class_alloc(pclass_sz->array_nl_index_at_z_and_k_no_wiggles,
            sizeof(double *)*pclass_sz->ndim_redshifts*pclass_sz->ln_k_size_for_tSZ,
            pclass_sz->error_message);


class_alloc(array_nl_index_at_z_and_k,
            pclass_sz->ndim_redshifts*sizeof(double *),
            pclass_sz->error_message);
class_alloc(array_nl_index_at_z_and_k_no_wiggles,
            pclass_sz->ndim_redshifts*sizeof(double *),
            pclass_sz->error_message);


// class_alloc(pclass_sz->array_nl_index_at_z_and_k_splined,
//             sizeof(double *)*pclass_sz->ndim_redshifts*pclass_sz->ln_k_size_for_tSZ,
//             pclass_sz->error_message);
//
//
// class_alloc(array_nl_index_at_z_and_k_splined,
//             pclass_sz->ndim_redshifts*sizeof(double *),
//             pclass_sz->error_message);

for (index_l=0;
     index_l<pclass_sz->ndim_redshifts;
     index_l++)
{
  class_alloc(array_nl_index_at_z_and_k[index_l],
              pclass_sz->ln_k_size_for_tSZ*sizeof(double),
              pclass_sz->error_message);
  class_alloc(array_nl_index_at_z_and_k_no_wiggles[index_l],
              pclass_sz->ln_k_size_for_tSZ*sizeof(double),
              pclass_sz->error_message);
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
pba,pnl,ppm,pclass_sz,z_min,z_max,logk_min,logk_max)\
private(tstart, tstop,index_k,index_z) \
num_threads(number_of_threads)
{

#ifdef _OPENMP
  tstart = omp_get_wtime();
#endif

  //
  // class_alloc_parallel(pvectsz,pclass_sz->tsz_size*sizeof(double),pclass_sz->error_message);
  //
  // class_alloc_parallel(pvecback,pba->bg_size*sizeof(double),pba->error_message);

#pragma omp for schedule (dynamic)
for (index_z=0; index_z<pclass_sz->ndim_redshifts; index_z++)
{

#pragma omp flush(abort)

for (index_k=0; index_k<pclass_sz->ln_k_size_for_tSZ; index_k++)
{

      //background quantities @ z:
      double z =   exp(pclass_sz->array_redshift[index_z])-1.;
      double logk =   pclass_sz->ln_k_for_tSZ[index_k];


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

    // printf("zk loop done.\n");

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
gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, pclass_sz->ln_k_size_for_tSZ);
gsl_spline_init (spline, pclass_sz->ln_k_for_tSZ, array_nl_index_at_z_and_k[index_z], pclass_sz->ln_k_size_for_tSZ);
// find where the second derivative changes sign:
int it;
// for (it=0;it<10;it++){
// printf ("%g %g %g\n", pclass_sz->ln_k_for_tSZ[it], array_nl_index_at_z_and_k[index_z][it],gsl_spline_eval (spline, pclass_sz->ln_k_for_tSZ[it], acc));
// }
// printf(" ####### ");
// printf("nk = %d\n",pclass_sz->ln_k_size_for_tSZ);
// exit(0);


struct Parameters_for_nl_fitting_function V;
V.acc = acc;
V.spline = spline;
void * params = &V;

// for (it=0;it<10;it++){
// printf ("nl_fitting_function %g %g %g\n", pclass_sz->ln_k_for_tSZ[it], array_nl_index_at_z_and_k[index_z][it],nl_fitting_function (pclass_sz->ln_k_for_tSZ[it], params));
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
// lnkl = pclass_sz->ln_k_for_tSZ[index_k]+2.*tol-tol;
// lnkr = pclass_sz->ln_k_for_tSZ[index_k]+2.*tol+tol;
// lnkc = pclass_sz->ln_k_for_tSZ[index_k]+2.*tol;
tol = 1.e-6;//fabs(pclass_sz->ln_k_for_tSZ[index_k+1] - pclass_sz->ln_k_for_tSZ[index_k])/100.;
lnkl = pclass_sz->ln_k_for_tSZ[index_k]+2.*tol;
lnkc = pclass_sz->ln_k_for_tSZ[index_k]+tol;
lnkr = pclass_sz->ln_k_for_tSZ[index_k];

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


for (index_k=1; index_k<pclass_sz->ln_k_size_for_tSZ-1; index_k++){
tol = 1.e-6;//fabs(pclass_sz->ln_k_for_tSZ[index_k+1] - pclass_sz->ln_k_for_tSZ[index_k])/100.;
lnkl = pclass_sz->ln_k_for_tSZ[index_k]+2.*tol;
lnkc = pclass_sz->ln_k_for_tSZ[index_k]+tol;
lnkr = pclass_sz->ln_k_for_tSZ[index_k];



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
  // lnx_zeros[n_data] = pclass_sz->ln_k_for_tSZ[index_k];
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
double nod_1_x = (pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero_next]] + pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero]])/2.;
double nod_1_y = (array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero_next]] + array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero]])/2.;

double nod_2_x = (pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero_next+2]] + pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero+2]])/2.;
double nod_2_y = (array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero_next+2]] + array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero+2]])/2.;




  for (index_k=0; index_k<pclass_sz->ln_k_size_for_tSZ; index_k++){
    // if (pclass_sz->ln_k_for_tSZ[index_k]<=lnx_zeros[0]){
    if (index_k<=lnx_zeros[0]){
      // if (pclass_sz->ln_k_for_tSZ[index_k] >pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero]] - (nod_2_x-nod_1_x)/2. ){

      double  array_nl_no_wiggles_interp = nod_1_y + (nod_2_y-nod_1_y)/(nod_2_x-nod_1_x)*(pclass_sz->ln_k_for_tSZ[index_k]-nod_1_x);
      double nl_exact = array_nl_index_at_z_and_k[index_z][index_k];
      // }

      // else{
      if ((nl_exact<array_nl_no_wiggles_interp) & (pclass_sz->ln_k_for_tSZ[index_k] >pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero]] - (nod_2_x-nod_1_x)))
      array_nl_no_wiggles = array_nl_no_wiggles_interp;
      else
      array_nl_no_wiggles = nl_exact;
    // }
    }
    // else if (pclass_sz->ln_k_for_tSZ[index_k]>=lnx_zeros[n_data-1]){
    else if (index_k+1>=lnx_zeros[n_data-1]){
      array_nl_no_wiggles = array_nl_index_at_z_and_k[index_z][index_k];
    }
    else{

      // linear interpolation:
      // array_nl_no_wiggles = nl_fitting_function(lnx_zeros[id_zero],params)
      //                               + (nl_fitting_function(lnx_zeros[id_zero_next],params)
      //                                  -nl_fitting_function(lnx_zeros[id_zero],params))
      //                               /(lnx_zeros[id_zero_next]-lnx_zeros[id_zero])*(pclass_sz->ln_k_for_tSZ[index_k]-lnx_zeros[id_zero]);
      nod_1_x = (pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero_next]] + pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero]])/2.;
      nod_1_y = (array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero_next]] + array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero]])/2.;

      nod_2_x = (pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero_next+2]] + pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero+2]])/2.;
      nod_2_y = (array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero_next+2]] + array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero+2]])/2.;


      // array_nl_no_wiggles = array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero]]
      //                               + (array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero_next]]
      //                                  -array_nl_index_at_z_and_k[index_z][lnx_zeros[id_zero]])
      //                               /(pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero_next]]-pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero]])*(pclass_sz->ln_k_for_tSZ[index_k]-pclass_sz->ln_k_for_tSZ[lnx_zeros[id_zero]]);

      array_nl_no_wiggles = nod_1_y + (nod_2_y-nod_1_y)/(nod_2_x-nod_1_x)*(pclass_sz->ln_k_for_tSZ[index_k]-nod_1_x);

    // if k larger than next node -> switch range:
    // if (pclass_sz->ln_k_for_tSZ[index_k+1]==lnx_zeros[id_zero_next]){
    if (index_k+1==lnx_zeros[id_zero_next+1]){
      id_zero_next += 2;
      id_zero += 2;
    }

      // array_nl_no_wiggles = nl_fitting_function(lnx_zeros[0],params)
      //                               + (nl_fitting_function(lnx_zeros[n_data-2],params)-nl_fitting_function(lnx_zeros[0],params))
      //                               /(lnx_zeros[n_data-2]-lnx_zeros[0])*(pclass_sz->ln_k_for_tSZ[index_k]-lnx_zeros[0]);

    }
   array_nl_index_at_z_and_k_no_wiggles[index_z][index_k] = array_nl_no_wiggles;
  }

}
else{
  for (index_k=0; index_k<pclass_sz->ln_k_size_for_tSZ; index_k++){
    array_nl_index_at_z_and_k_no_wiggles[index_z][index_k] = array_nl_index_at_z_and_k[index_z][index_k];
  }
}


gsl_spline_free (spline);
gsl_interp_accel_free (acc);
free(lnx_zeros);





  }
#ifdef _OPENMP
  tstop = omp_get_wtime();
  if (pclass_sz->sz_verbose > 0)
    printf("In %s: time spent in parallel region (loop over z's) = %e s for thread %d\n",
           __func__,tstop-tstart,omp_get_thread_num());
#endif

    // free(pvecback);
    // free(pvectsz);
    }
if (abort == _TRUE_) return _FAILURE_;
//end of parallel region

index_z_k = 0;
for (index_k=0; index_k<pclass_sz->ln_k_size_for_tSZ; index_k++)
{
  for (index_z=0; index_z<pclass_sz->ndim_redshifts; index_z++)
  {
    pclass_sz->array_nl_index_at_z_and_k[index_z_k] = array_nl_index_at_z_and_k[index_z][index_k];
    pclass_sz->array_nl_index_at_z_and_k_no_wiggles[index_z_k] = array_nl_index_at_z_and_k_no_wiggles[index_z][index_k];
    // printf("index_z = %d index_k =  %d dndlnk_nw = %.3e dndlnk = %.3e\n",index_z, index_k, pclass_sz->array_nl_index_at_z_and_k_no_wiggles[index_z_k], pclass_sz->array_nl_index_at_z_and_k[index_z_k]);

    index_z_k += 1;

  }
}

for (index_z=0; index_z<pclass_sz->ndim_redshifts; index_z++){
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
                                struct class_sz_structure * pclass_sz){


double * sigma2_hsv_var;
class_alloc(sigma2_hsv_var,
            sizeof(double *),
            pclass_sz->error_message);


class_alloc(pclass_sz->array_sigma2_hsv_at_z,sizeof(double *)*pclass_sz->ndim_redshifts,pclass_sz->error_message);

int index_z;

for (index_z=0; index_z<pclass_sz->ndim_redshifts; index_z++)
        {

            spectra_sigma2_hsv(pba,
                                ppm,
                                pnl,
                                pclass_sz,
                                exp(pclass_sz->array_redshift[index_z])-1.,
                                sigma2_hsv_var
                                );
          pclass_sz->array_sigma2_hsv_at_z[index_z] = log(*sigma2_hsv_var);
          //printf("%.4e \t %.4e\n",exp(pclass_sz->array_redshift[index_z])-1.,pclass_sz->array_sigma2_hsv_at_z[index_z]);
       }

free(sigma2_hsv_var);

return _SUCCESS_;
    }





//Tabulate k non linear as functions of redshift
int tabulate_knl(struct background * pba,
                 struct nonlinear * pnl,
                 struct primordial * ppm,
                 struct class_sz_structure * pclass_sz){



double knl_var;
class_alloc(pclass_sz->array_knl_at_z,sizeof(double *)*pclass_sz->ndim_redshifts,pclass_sz->error_message);

int index_z;
double z;
for (index_z=0; index_z<pclass_sz->ndim_redshifts; index_z++)
        {
          z = exp(pclass_sz->array_redshift[index_z])-1.;
          solve_pkl_to_knl(&knl_var,
          z,
          pclass_sz,
          pba,
          pnl,
          ppm);

          pclass_sz->array_knl_at_z[index_z] = log(knl_var);
          // printf("z = %.4e \t knl = %.4e\n",z,knl_var);
       }


return _SUCCESS_;
    }



double get_planck_sigma_at_theta500(double theta500, struct class_sz_structure * pclass_sz){
  double y;
  int l1,l2;
  double th1,th2;

  if ((theta500<pclass_sz->thetas[0])){
       l1 = 0;
       l2 = 1;
       th1 = pclass_sz->thetas[l1];
       th2 = pclass_sz->thetas[l2];
    double y1 = pclass_sz->sky_averaged_ylims[l1];
    double y2 = pclass_sz->sky_averaged_ylims[l2];
    double y = y1 + (y2-y1)/(th2-th1)*(theta500-th1);
    return y;
  }
  else if ((theta500>pclass_sz->thetas[pclass_sz->nthetas-1])){
      l1 = pclass_sz->nthetas - 1;
      l2 = pclass_sz->nthetas - 2;
      th1 = pclass_sz->thetas[l1];
      th2 = pclass_sz->thetas[l2];
    double y1 = pclass_sz->sky_averaged_ylims[l1];
    double y2 = pclass_sz->sky_averaged_ylims[l2];
    double y = y1 + (y2-y1)/(th2-th1)*(theta500-th1);
    return y;
  }
  else{
  return pwl_value_1d(pclass_sz->nthetas,
                      pclass_sz->thetas,
                      pclass_sz->sky_averaged_ylims,
                      theta500);
                    }
}

double get_knl_at_z(double z, struct class_sz_structure * pclass_sz){
   double z_asked = log(1.+z);
 if (z<exp(pclass_sz->array_redshift[0])-1.)
    z_asked = pclass_sz->array_redshift[0];
 if (z>exp(pclass_sz->array_redshift[pclass_sz->ndim_redshifts-1])-1.)
    z_asked =  pclass_sz->array_redshift[pclass_sz->ndim_redshifts-1];
 return  exp(pwl_value_1d(pclass_sz->ndim_redshifts,
                          pclass_sz->array_redshift,
                          pclass_sz->array_knl_at_z,
                          z_asked));
}

double get_nl_index_at_z_and_k(double z_asked, double k_asked, struct class_sz_structure * pclass_sz, struct nonlinear * pnl){
  double z = log(1.+z_asked);
  double k = log(k_asked); // in h/Mpc

 if (z_asked<exp(pclass_sz->array_redshift[0])-1.)
    z = pclass_sz->array_redshift[0];
 if (z_asked>exp(pclass_sz->array_redshift[pclass_sz->ndim_redshifts-1])-1.)
    z =  pclass_sz->array_redshift[pclass_sz->ndim_redshifts-1];

 if (k_asked<exp(pclass_sz->ln_k_for_tSZ[0]))
    k =  pclass_sz->ln_k_for_tSZ[0];
 if (k_asked>exp(pclass_sz->ln_k_for_tSZ[pclass_sz->ln_k_size_for_tSZ-1]))
    k =  pclass_sz->ln_k_for_tSZ[pclass_sz->ln_k_size_for_tSZ-1];

 return pwl_interp_2d(pclass_sz->ndim_redshifts,
                      pclass_sz->ln_k_size_for_tSZ,
                      pclass_sz->array_redshift,
                      pclass_sz->ln_k_for_tSZ,
                      pclass_sz->array_nl_index_at_z_and_k,
                      1,
                      &z,
                      &k);
}
//
double get_nl_index_at_z_and_k_no_wiggles(double z_asked, double k_asked, struct class_sz_structure * pclass_sz, struct nonlinear * pnl){
  double z = log(1.+z_asked);
  double k = log(k_asked); // in h/Mpc

 if (z_asked<exp(pclass_sz->array_redshift[0])-1.)
    z = pclass_sz->array_redshift[0];
 if (z_asked>exp(pclass_sz->array_redshift[pclass_sz->ndim_redshifts-1])-1.)
    z =  pclass_sz->array_redshift[pclass_sz->ndim_redshifts-1];

 if (k_asked<exp(pclass_sz->ln_k_for_tSZ[0]))
    k =  pclass_sz->ln_k_for_tSZ[0];
 if (k_asked>exp(pclass_sz->ln_k_for_tSZ[pclass_sz->ln_k_size_for_tSZ-1]))
    k =  pclass_sz->ln_k_for_tSZ[pclass_sz->ln_k_size_for_tSZ-1];

 return pwl_interp_2d(pclass_sz->ndim_redshifts,
                      pclass_sz->ln_k_size_for_tSZ,
                      pclass_sz->array_redshift,
                      pclass_sz->ln_k_for_tSZ,
                      pclass_sz->array_nl_index_at_z_and_k_no_wiggles,
                      1,
                      &z,
                      &k);
}

double get_completeness_at_z_and_M(double z_asked, double m_asked, double * completeness_2d_to_1d, struct class_sz_structure * pclass_sz){
  if (pclass_sz->has_completeness == 0){
    return 1.;
  }
  else{

  double z = z_asked;

  double m = log(m_asked);
  // printf("z = %.3e m = %.3e logm = %.4e mmin = %.4e mmax = %,4e\n",z_asked,m_asked,m,pclass_sz->steps_m[0],pclass_sz->steps_m[pclass_sz->nsteps_m-1]);
  if (m < pclass_sz->steps_m[0])
    return 1e-100;//m = pclass_sz->steps_m[0];
  if (m > pclass_sz->steps_m[pclass_sz->nsteps_m-1])
    return 1e-100;//m = pclass_sz->steps_m[pclass_sz->nsteps_m-1];
  if (z < pclass_sz->steps_z[0])
    return 1e-100;//z = pclass_sz->steps_z[0];
  if (z > pclass_sz->steps_z[pclass_sz->nsteps_z-1])
    return 1e-100;//z = pclass_sz->steps_z[pclass_sz->nsteps_z-1];
 return exp(pwl_interp_2d(pclass_sz->nsteps_m,
                          pclass_sz->nsteps_z,
                          pclass_sz->steps_m,
                          pclass_sz->steps_z,
                          completeness_2d_to_1d,
                          1,
                          &m,
                          &z));

  }
}



double get_detection_proba_at_y_and_theta(double y_asked, double th_asked, double * erfs_2d_to_1d, struct class_sz_structure * pclass_sz){
  double y = log(y_asked);
  double th = log(th_asked);
  double r = 0.;
  // printf("z = %.3e m = %.3e logm = %.4e mmin = %.4e mmax = %,4e\n",z_asked,m_asked,m,pclass_sz->steps_m[0],pclass_sz->steps_m[pclass_sz->nsteps_m-1]);
  if (y < pclass_sz->erfs_2d_to_1d_y_array[0]){
    r = 0.;}
  else if (y > pclass_sz->erfs_2d_to_1d_y_array[pclass_sz->Ny-1]){
    r = 0.;}

  else if (th < pclass_sz->erfs_2d_to_1d_th_array[0]){
    double th1 = pclass_sz->erfs_2d_to_1d_th_array[0];
    double th2 = pclass_sz->erfs_2d_to_1d_th_array[1];
    double r1 = exp(pwl_interp_2d( pclass_sz->Ny,
                          pclass_sz->Nth,
                          pclass_sz->erfs_2d_to_1d_y_array,
                          pclass_sz->erfs_2d_to_1d_th_array,
                          erfs_2d_to_1d,
                          1,
                          &y,
                          &th1));
    double r2 = exp(pwl_interp_2d( pclass_sz->Ny,
                          pclass_sz->nthetas,
                          pclass_sz->erfs_2d_to_1d_y_array,
                          pclass_sz->erfs_2d_to_1d_th_array,
                          erfs_2d_to_1d,
                          1,
                          &y,
                          &th2));
    r = r1 + (r2 - r1)/(exp(th2)-exp(th1))*(exp(th)-exp(th1));

  }
  else if (th > pclass_sz->erfs_2d_to_1d_th_array[pclass_sz->Nth-1]){
    double th1 = pclass_sz->erfs_2d_to_1d_th_array[pclass_sz->Nth-1];
    double th2 = pclass_sz->erfs_2d_to_1d_th_array[pclass_sz->Nth-2];
    double r1 = exp(pwl_interp_2d( pclass_sz->Ny,
                          pclass_sz->Nth,
                          pclass_sz->erfs_2d_to_1d_y_array,
                          pclass_sz->erfs_2d_to_1d_th_array,
                          erfs_2d_to_1d,
                          1,
                          &y,
                          &th1));
    double r2 = exp(pwl_interp_2d( pclass_sz->Ny,
                          pclass_sz->Nth,
                          pclass_sz->erfs_2d_to_1d_y_array,
                          pclass_sz->erfs_2d_to_1d_th_array,
                          erfs_2d_to_1d,
                          1,
                          &y,
                          &th2));
    r = r1 + (r2 - r1)/(exp(th2)-exp(th1))*(exp(th)-exp(th1));

}

  else{
  // r = exp(pwl_interp_2d(pclass_sz->nthetas,
  //                         pclass_sz->Ny,
  //                         pclass_sz->thetas,
  //                         pclass_sz->erfs_2d_to_1d_y_array,
  //                         erfs_2d_to_1d,
  //                         1,
  //                         &th,
  //                         &y));
  r = exp(pwl_interp_2d( pclass_sz->Ny,
                          pclass_sz->Nth,
                          pclass_sz->erfs_2d_to_1d_y_array,
                          pclass_sz->erfs_2d_to_1d_th_array,
                          erfs_2d_to_1d,
                          1,
                          &y,
                          &th));
                        }
  return r;
}


double get_rho_2h_at_r_and_m_and_z(double r_asked,
                                   double m_asked,
                                   double z_asked,
                                   struct class_sz_structure * pclass_sz,
                                   struct background * pba){
  double z = log(1.+z_asked);
  double r = log(r_asked);

 if (z<pclass_sz->array_profile_ln_1pz[0]){
    // z = pclass_sz->array_profile_ln_1pz[0];
    return 1e-100;
  }
 else if (z>pclass_sz->array_profile_ln_1pz[pclass_sz->n_z_density_profile-1]){
    // z = pclass_sz->array_profile_ln_1pz[pclass_sz->n_z_density_profile-1];
    return 1e-100;
  }

 else if (r<pclass_sz->array_profile_ln_r[0]){
    // r = pclass_sz->array_profile_ln_r[0];
    return 1e-100;
  }
      // printf("dealing with mass conversion in hmf3\n");
 else if (r>pclass_sz->array_profile_ln_r[pclass_sz->N_samp_fftw-1]){
    // k =  pclass_sz->array_profile_ln_k[pclass_sz->n_k_density_profile-1]
    return 1e-100;
  }


// printf("l=%.3e\n",l);
else {

 double result = pwl_interp_2d(

                          pclass_sz->n_z_density_profile,
                          pclass_sz->N_samp_fftw,


                          pclass_sz->array_profile_ln_1pz,
                          pclass_sz->array_profile_ln_r,

                          pclass_sz->array_profile_rho_2h_at_r_and_z,

                          1,
                          &z,
                          &r);

 double nu = get_nu_at_z_and_m(z_asked,m_asked,pclass_sz,pba);
 double b_at_m = get_first_order_bias_at_z_and_nu(z_asked,nu,pclass_sz);
 result *= b_at_m;
 return result;
  }
}




double get_gas_pressure_2h_at_r_and_m_and_z(double r_asked,
                                            double m_asked,
                                            double z_asked,
                                            struct class_sz_structure * pclass_sz,
                                            struct background * pba){
  double z = log(1.+z_asked);
  double r = log(r_asked);

 if (z<pclass_sz->array_pressure_profile_ln_1pz[0]){
    // z = pclass_sz->array_profile_ln_1pz[0];
    return 1e-100;
  }
 else if (z>pclass_sz->array_pressure_profile_ln_1pz[pclass_sz->n_z_pressure_profile-1]){
    // z = pclass_sz->array_profile_ln_1pz[pclass_sz->n_z_density_profile-1];
    return 1e-100;
  }

 else if (r<pclass_sz->array_pressure_profile_ln_r[0]){
    // r = pclass_sz->array_profile_ln_r[0];
    return 1e-100;
  }
      // printf("dealing with mass conversion in hmf3\n");
 else if (r>pclass_sz->array_pressure_profile_ln_r[pclass_sz->N_samp_fftw-1]){
    // k =  pclass_sz->array_profile_ln_k[pclass_sz->n_k_density_profile-1]
    return 1e-100;
  }


// printf("l=%.3e\n",l);
else {

 double result = pwl_interp_2d(

                          pclass_sz->n_z_pressure_profile,
                          pclass_sz->N_samp_fftw,


                          pclass_sz->array_pressure_profile_ln_1pz,
                          pclass_sz->array_pressure_profile_ln_r,

                          pclass_sz->array_pressure_profile_pressure_2h_at_r_and_z,

                          1,
                          &z,
                          &r);

  if (isnan(result ) || isinf(result )){
    printf("found nan in pressure_profile_pressure_2h_at_r_and_z.\n");
    exit(0);
  }

 double nu = get_nu_at_z_and_m(z_asked,m_asked,pclass_sz,pba);
 double b_at_m = get_first_order_bias_at_z_and_nu(z_asked,nu,pclass_sz);
 result *= b_at_m;
 return result;
  }
}





double get_gas_pressure_2h_at_k_and_z(double k_asked, double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double k = log(k_asked);

 if (z<pclass_sz->array_pressure_profile_ln_1pz[0]){
    return 1e-100;
  }
 else if (z>pclass_sz->array_pressure_profile_ln_1pz[pclass_sz->n_z_pressure_profile-1]){
    return 1e-100;
  }

 else if (k<pclass_sz->array_pressure_profile_2h_ln_k[0]){
    return 1e-100;
  }
 else if (k>pclass_sz->array_pressure_profile_2h_ln_k[pclass_sz->n_k_pressure_profile_2h-1]){
    return 1e-100;
  }


// printf("l=%.3e\n",l);
else {

 double result = pwl_interp_2d(
                          pclass_sz->n_z_pressure_profile,
                          pclass_sz->n_k_pressure_profile_2h,



                          pclass_sz->array_pressure_profile_ln_1pz,
                          pclass_sz->array_pressure_profile_2h_ln_k,

                          pclass_sz->array_pressure_profile_ln_pressure_2h_at_k_and_z,
                          1,
                          &z,
                          &k);

   if (isnan(result ) || isinf(result )){
    printf("found nan in pressure_profile_ln_pressure_2h_at_k_and_z.\n");
    exit(0);
    }

 return exp(result);
  }
}



double get_rho_2h_at_k_and_z(double k_asked, double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double k = log(k_asked);

 if (z<pclass_sz->array_profile_ln_1pz[0]){
    // z = pclass_sz->array_profile_ln_1pz[0];
    return 1e-100;
  }
 else if (z>pclass_sz->array_profile_ln_1pz[pclass_sz->n_z_density_profile-1]){
    // z = pclass_sz->array_profile_ln_1pz[pclass_sz->n_z_density_profile-1];
    return 1e-100;
  }

 else if (k<pclass_sz->array_profile_ln_k[0]){
    // k = pclass_sz->array_profile_ln_k[0];
    return 1e-100;
  }
      // printf("dealing with mass conversion in hmf3\n");
 else if (k>pclass_sz->array_profile_ln_k[pclass_sz->n_k_density_profile-1]){
    // k =  pclass_sz->array_profile_ln_k[pclass_sz->n_k_density_profile-1]
    return 1e-100;
  }


// printf("l=%.3e\n",l);
else {
 return exp(pwl_interp_2d(
                          pclass_sz->n_z_density_profile,
                          pclass_sz->n_k_density_profile,



                          pclass_sz->array_profile_ln_1pz,
                          pclass_sz->array_profile_ln_k,

                          pclass_sz->array_profile_ln_rho_2h_at_k_and_z,
                          1,
                          &z,
                          &k));
  }
}



double get_psi_b2t_at_k_and_z(double l_asked, double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double l = log(l_asked);

 if (z<pclass_sz->array_psi_b2t_redshift[0])
    return 0.;//z = pclass_sz->array_psi_b2t_redshift[0];
 if (z>pclass_sz->array_psi_b2t_redshift[pclass_sz->n_z_psi_b2t-1])
    return 0.;//z = pclass_sz->array_psi_b2t_redshift[pclass_sz->n_z_psi_b2t-1];

 if (l<pclass_sz->array_psi_b2t_multipole[0])
    return 0.;//l = pclass_sz->array_psi_b2t_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l>pclass_sz->array_psi_b2t_multipole[pclass_sz->n_l_psi_b2t-1])
    return 0.;//l =  pclass_sz->array_psi_b2t_multipole[pclass_sz->n_l_psi_b2t-1];


// printf("l=%.3e\n",l);

 return exp(pwl_interp_2d(

                          pclass_sz->n_z_psi_b2t,
                          pclass_sz->n_l_psi_b2t,

                          pclass_sz->array_psi_b2t_redshift,
                          pclass_sz->array_psi_b2t_multipole,
                          pclass_sz->array_psi_b2t_psi,
                          1,
                          &z,
                          &l))-1.;
}



double get_psi_b2g_at_k_and_z(double l_asked, double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double l = log(l_asked);

 if (z<pclass_sz->array_psi_b2g_redshift[0])
    return 0.;//z = pclass_sz->array_psi_b2g_redshift[0];
 if (z>pclass_sz->array_psi_b2g_redshift[pclass_sz->n_z_psi_b2g-1])
    return 0.;//z = pclass_sz->array_psi_b2g_redshift[pclass_sz->n_z_psi_b2g-1];

 if (l<pclass_sz->array_psi_b2g_multipole[0])
    return 0.;//l = pclass_sz->array_psi_b2g_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l>pclass_sz->array_psi_b2g_multipole[pclass_sz->n_l_psi_b2g-1])
    return 0.;//l =  pclass_sz->array_psi_b2g_multipole[pclass_sz->n_l_psi_b2g-1];


// printf("l=%.3e\n",l);

 return exp(pwl_interp_2d(

                          pclass_sz->n_z_psi_b2g,
                          pclass_sz->n_l_psi_b2g,

                          pclass_sz->array_psi_b2g_redshift,
                          pclass_sz->array_psi_b2g_multipole,
                          pclass_sz->array_psi_b2g_psi,
                          1,
                          &z,
                          &l))-1.;
}


double get_psi_b2kg_at_k_and_z(double l_asked, double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double l = log(l_asked);

 if (z<pclass_sz->array_psi_b2kg_redshift[0])
    return 0.;//z = pclass_sz->array_psi_b2kg_redshift[0];
 if (z>pclass_sz->array_psi_b2kg_redshift[pclass_sz->n_z_psi_b2kg-1])
    return 0.;//z = pclass_sz->array_psi_b2kg_redshift[pclass_sz->n_z_psi_b2kg-1];

 if (l<pclass_sz->array_psi_b2kg_multipole[0])
    return 0.;//l = pclass_sz->array_psi_b2kg_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l>pclass_sz->array_psi_b2kg_multipole[pclass_sz->n_l_psi_b2kg-1])
    return 0.;//l =  pclass_sz->array_psi_b2kg_multipole[pclass_sz->n_l_psi_b2kg-1];


// printf("l=%.3e\n",l);

 return exp(pwl_interp_2d(

                          pclass_sz->n_z_psi_b2kg,
                          pclass_sz->n_l_psi_b2kg,

                          pclass_sz->array_psi_b2kg_redshift,
                          pclass_sz->array_psi_b2kg_multipole,
                          pclass_sz->array_psi_b2kg_psi,
                          1,
                          &z,
                          &l))-1.;
}




double get_psi_b1g_at_k_and_z(double l_asked, double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double l = log(l_asked);

 if (z<pclass_sz->array_psi_b1g_redshift[0])
    return 0.;//z = pclass_sz->array_psi_b1g_redshift[0];
 if (z>pclass_sz->array_psi_b1g_redshift[pclass_sz->n_z_psi_b1g-1])
    return 0.;//z = pclass_sz->array_psi_b1g_redshift[pclass_sz->n_z_psi_b1g-1];

 if (l<pclass_sz->array_psi_b1g_multipole[0])
    return 0.;//l = pclass_sz->array_psi_b1g_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l>pclass_sz->array_psi_b1g_multipole[pclass_sz->n_l_psi_b1g-1])
    return 0.;//l =  pclass_sz->array_psi_b1g_multipole[pclass_sz->n_l_psi_b1g-1];


// printf("l=%.3e\n",l);

 return exp(pwl_interp_2d(

                          pclass_sz->n_z_psi_b1g,
                          pclass_sz->n_l_psi_b1g,

                          pclass_sz->array_psi_b1g_redshift,
                          pclass_sz->array_psi_b1g_multipole,
                          pclass_sz->array_psi_b1g_psi,
                          1,
                          &z,
                          &l));
}



double get_psi_b1kg_at_k_and_z(double l_asked, double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double l = log(l_asked);

 if (z<pclass_sz->array_psi_b1kg_redshift[0])
    return 0.;//z = pclass_sz->array_psi_b1kg_redshift[0];
 if (z>pclass_sz->array_psi_b1kg_redshift[pclass_sz->n_z_psi_b1kg-1])
    return 0.;//z = pclass_sz->array_psi_b1kg_redshift[pclass_sz->n_z_psi_b1kg-1];

 if (l<pclass_sz->array_psi_b1kg_multipole[0])
    return 0.;//l = pclass_sz->array_psi_b1kg_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l>pclass_sz->array_psi_b1kg_multipole[pclass_sz->n_l_psi_b1kg-1])
    return 0.;//l =  pclass_sz->array_psi_b1kg_multipole[pclass_sz->n_l_psi_b1kg-1];


// printf("l=%.3e\n",l);

 return exp(pwl_interp_2d(

                          pclass_sz->n_z_psi_b1kg,
                          pclass_sz->n_l_psi_b1kg,

                          pclass_sz->array_psi_b1kg_redshift,
                          pclass_sz->array_psi_b1kg_multipole,
                          pclass_sz->array_psi_b1kg_psi,
                          1,
                          &z,
                          &l));
}





double get_psi_b1t_at_k_and_z(double l_asked, double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double l = log(l_asked);

 if (z<pclass_sz->array_psi_b1t_redshift[0])
    return 0.;//z = pclass_sz->array_psi_b1t_redshift[0];
 if (z>pclass_sz->array_psi_b1t_redshift[pclass_sz->n_z_psi_b1t-1])
    return 0.;//z = pclass_sz->array_psi_b1t_redshift[pclass_sz->n_z_psi_b1t-1];

 if (l<pclass_sz->array_psi_b1t_multipole[0])
    return 0.;//l = pclass_sz->array_psi_b1t_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l>pclass_sz->array_psi_b1t_multipole[pclass_sz->n_l_psi_b1t-1])
    return 0.;//l =  pclass_sz->array_psi_b1t_multipole[pclass_sz->n_l_psi_b1t-1];


 return exp(pwl_interp_2d(

                          pclass_sz->n_z_psi_b1t,
                          pclass_sz->n_l_psi_b1t,

                          pclass_sz->array_psi_b1t_redshift,
                          pclass_sz->array_psi_b1t_multipole,
                          pclass_sz->array_psi_b1t_psi,
                          1,
                          &z,
                          &l));
}




double get_psi_b1gt_at_k1_k2_and_z(double l_asked,double l_asked2, double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double l1 = log(l_asked);
  double l2 = log(l_asked2);
// printf("z=%.8e\n",z);
 if (z<pclass_sz->array_psi_b1gt_redshift[0])
    return 0.;//z = pclass_sz->array_psi_b1gt_redshift[0];
 if (z>pclass_sz->array_psi_b1gt_redshift[pclass_sz->n_z_psi_b1gt-1])
    return 0.;//z = pclass_sz->array_psi_b1gt_redshift[pclass_sz->n_z_psi_b1gt-1];

 if (l1<pclass_sz->array_psi_b1gt_multipole[0])
    return 0.;//l1 = pclass_sz->array_psi_b1gt_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l1>pclass_sz->array_psi_b1gt_multipole[pclass_sz->n_l_psi_b1gt-1])
    return 0.;//l1 =  pclass_sz->array_psi_b1gt_multipole[pclass_sz->n_l_psi_b1gt-1];

 if (l2<pclass_sz->array_psi_b1gt_multipole[0])
    return 0.;//l2 = pclass_sz->array_psi_b1gt_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l2>pclass_sz->array_psi_b1gt_multipole[pclass_sz->n_l_psi_b1gt-1])
    return 0.;//l2 =  pclass_sz->array_psi_b1gt_multipole[pclass_sz->n_l_psi_b1gt-1];

  // find the closest z's in the grid:
  int id_z_low;
  int id_z_up;
  r8vec_bracket(pclass_sz->n_z_psi_b1gt,pclass_sz->array_psi_b1gt_redshift,z,&id_z_low,&id_z_up);


 double ln_rho_low = pwl_interp_2d(pclass_sz->n_l_psi_b1gt,
                                  pclass_sz->n_l_psi_b1gt,
                                  pclass_sz->array_psi_b1gt_multipole,
                                  pclass_sz->array_psi_b1gt_multipole,
                                  pclass_sz->array_psi_b1gt_psi[id_z_low-1],
                                  1,
                                  &l1,
                                  &l2);

 double ln_rho_up = pwl_interp_2d(pclass_sz->n_l_psi_b1gt,
                                  pclass_sz->n_l_psi_b1gt,
                                  pclass_sz->array_psi_b1gt_multipole,
                                  pclass_sz->array_psi_b1gt_multipole,
                                  pclass_sz->array_psi_b1gt_psi[id_z_up-1],
                                  1,
                                  &l1,
                                  &l2);
 double ln_l_low = pclass_sz->array_psi_b1gt_redshift[id_z_low-1];
 double ln_l_up = pclass_sz->array_psi_b1gt_redshift[id_z_up-1];
 double result =  exp(ln_rho_low + ((z - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low));
 if (isnan(result)||isinf(result)){
   printf("get b1gt : z %.3e l_asked %.4e k1 %.4e k2 %.4e ln_rho_low %.4e ln_rho_up %.4e\n",z,l_asked,exp(l1),exp(l2),ln_rho_low,ln_rho_up);
   exit(0);
 }
 return result;


}




double get_psi_b1kgt_at_k1_k2_and_z(double l_asked,double l_asked2, double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double l1 = log(l_asked);
  double l2 = log(l_asked2);
// printf("z=%.8e\n",z);
 if (z<pclass_sz->array_psi_b1kgt_redshift[0])
    return 0.;//z = pclass_sz->array_psi_b1kgt_redshift[0];
 if (z>pclass_sz->array_psi_b1kgt_redshift[pclass_sz->n_z_psi_b1kgt-1])
    return 0.;//z = pclass_sz->array_psi_b1kgt_redshift[pclass_sz->n_z_psi_b1kgt-1];

 if (l1<pclass_sz->array_psi_b1kgt_multipole[0])
    return 0.;//l1 = pclass_sz->array_psi_b1kgt_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l1>pclass_sz->array_psi_b1kgt_multipole[pclass_sz->n_l_psi_b1kgt-1])
    return 0.;//l1 =  pclass_sz->array_psi_b1kgt_multipole[pclass_sz->n_l_psi_b1kgt-1];

 if (l2<pclass_sz->array_psi_b1kgt_multipole[0])
    return 0.;//l2 = pclass_sz->array_psi_b1kgt_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l2>pclass_sz->array_psi_b1kgt_multipole[pclass_sz->n_l_psi_b1kgt-1])
    return 0.;//l2 =  pclass_sz->array_psi_b1kgt_multipole[pclass_sz->n_l_psi_b1kgt-1];

  // find the closest z's in the grid:
  int id_z_low;
  int id_z_up;
  r8vec_bracket(pclass_sz->n_z_psi_b1kgt,pclass_sz->array_psi_b1kgt_redshift,z,&id_z_low,&id_z_up);


 double ln_rho_low = pwl_interp_2d(pclass_sz->n_l_psi_b1kgt,
                                  pclass_sz->n_l_psi_b1kgt,
                                  pclass_sz->array_psi_b1kgt_multipole,
                                  pclass_sz->array_psi_b1kgt_multipole,
                                  pclass_sz->array_psi_b1kgt_psi[id_z_low-1],
                                  1,
                                  &l1,
                                  &l2);

 double ln_rho_up = pwl_interp_2d(pclass_sz->n_l_psi_b1kgt,
                                  pclass_sz->n_l_psi_b1kgt,
                                  pclass_sz->array_psi_b1kgt_multipole,
                                  pclass_sz->array_psi_b1kgt_multipole,
                                  pclass_sz->array_psi_b1kgt_psi[id_z_up-1],
                                  1,
                                  &l1,
                                  &l2);
 double ln_l_low = pclass_sz->array_psi_b1kgt_redshift[id_z_low-1];
 double ln_l_up = pclass_sz->array_psi_b1kgt_redshift[id_z_up-1];
 double result =  exp(ln_rho_low + ((z - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low));
 if (isnan(result)||isinf(result)){
   printf("get b1kgt : z %.3e l_asked %.4e k1 %.4e k2 %.4e ln_rho_low %.4e ln_rho_up %.4e\n",z,l_asked,exp(l1),exp(l2),ln_rho_low,ln_rho_up);
   exit(0);
 }
 return result;


}



double get_psi_b1kgg_at_k1_k2_and_z(double l_asked,double l_asked2, double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double l1 = log(l_asked);
  double l2 = log(l_asked2);
// printf("z=%.8e\n",z);
 if (z<pclass_sz->array_psi_b1kgg_redshift[0])
    return 0.;//z = pclass_sz->array_psi_b1kgg_redshift[0];
 if (z>pclass_sz->array_psi_b1kgg_redshift[pclass_sz->n_z_psi_b1kgg-1])
    return 0.;//z = pclass_sz->array_psi_b1kgg_redshift[pclass_sz->n_z_psi_b1kgg-1];

 if (l1<pclass_sz->array_psi_b1kgg_multipole[0])
    return 0.;//l1 = pclass_sz->array_psi_b1kgg_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l1>pclass_sz->array_psi_b1kgg_multipole[pclass_sz->n_l_psi_b1kgg-1])
    return 0.;//l1 =  pclass_sz->array_psi_b1kgg_multipole[pclass_sz->n_l_psi_b1kgg-1];

 if (l2<pclass_sz->array_psi_b1kgg_multipole[0])
    return 0.;//l2 = pclass_sz->array_psi_b1kgg_multipole[0];
      // printf("dealing with mass conversion in hmf3\n");
 if (l2>pclass_sz->array_psi_b1kgg_multipole[pclass_sz->n_l_psi_b1kgg-1])
    return 0.;//l2 =  pclass_sz->array_psi_b1kgg_multipole[pclass_sz->n_l_psi_b1kgg-1];

  // find the closest z's in the grid:
  int id_z_low;
  int id_z_up;
  r8vec_bracket(pclass_sz->n_z_psi_b1kgg,pclass_sz->array_psi_b1kgg_redshift,z,&id_z_low,&id_z_up);


 double ln_rho_low = pwl_interp_2d(pclass_sz->n_l_psi_b1kgg,
                                  pclass_sz->n_l_psi_b1kgg,
                                  pclass_sz->array_psi_b1kgg_multipole,
                                  pclass_sz->array_psi_b1kgg_multipole,
                                  pclass_sz->array_psi_b1kgg_psi[id_z_low-1],
                                  1,
                                  &l1,
                                  &l2);

 double ln_rho_up = pwl_interp_2d(pclass_sz->n_l_psi_b1kgg,
                                  pclass_sz->n_l_psi_b1kgg,
                                  pclass_sz->array_psi_b1kgg_multipole,
                                  pclass_sz->array_psi_b1kgg_multipole,
                                  pclass_sz->array_psi_b1kgg_psi[id_z_up-1],
                                  1,
                                  &l1,
                                  &l2);
 double ln_l_low = pclass_sz->array_psi_b1kgg_redshift[id_z_low-1];
 double ln_l_up = pclass_sz->array_psi_b1kgg_redshift[id_z_up-1];
 double result =  exp(ln_rho_low + ((z - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low));
 if (isnan(result)||isinf(result)){
   printf("get b1kgt : z %.3e l_asked %.4e k1 %.4e k2 %.4e ln_rho_low %.4e ln_rho_up %.4e\n",z,l_asked,exp(l1),exp(l2),ln_rho_low,ln_rho_up);
   exit(0);
 }
 return result;


}



double get_dydz_at_z(double z_asked, struct class_sz_structure * pclass_sz)
{
  double z = log(1.+z_asked);


 if (z<pclass_sz->array_dydz_redshift[0])
    return 0.;
    // z = pclass_sz->array_dydz_redshift[0];
 if (z>pclass_sz->array_dydz_redshift[pclass_sz->n_z_dydz-1])
    return 0.;
    // z = pclass_sz->array_dydz_redshift[pclass_sz->n_z_dydz-1];




 return exp(pwl_value_1d(pclass_sz->n_z_dydz,
                         pclass_sz->array_dydz_redshift,
                         pclass_sz->array_dydz_at_z,
                         z));

}

double get_dNdlny_at_z_and_y(double z_asked, double y_asked, struct background * pba, struct class_sz_structure * pclass_sz){

  double z = log(1.+z_asked);
  double y = log(y_asked);

 if (z<pclass_sz->array_y_to_m_redshift[0])
    return 0.;//z = pclass_sz->array_y_to_m_redshift[0];
 if (z>pclass_sz->array_y_to_m_redshift[pclass_sz->n_z_y_to_m-1])
    return 0.;//z = pclass_sz->array_y_to_m_redshift[pclass_sz->n_z_y_to_m-1];

 if (y<pclass_sz->array_y_to_m_y[0])
    return 0.;//y = pclass_sz->array_y_to_m_y[0];

 if (y>pclass_sz->array_y_to_m_y[pclass_sz->n_y_y_to_m-1])
    return 0.;//y =  pclass_sz->array_y_to_m_y[pclass_sz->n_y_y_to_m-1];



 double m = get_y_to_m_at_z_and_y(z_asked,y_asked,pclass_sz);
 double dNdlnm = 1.;//get_volume_at_z(z_asked,pba)*get_dndlnM_at_z_and_M(z_asked,m,pclass_sz);
 dNdlnm = get_dndlnM_at_z_and_M(z_asked,m,pclass_sz);
 dNdlnm *= get_volume_at_z(z_asked,pba);
 double dlnmdlny = get_dlnm_dlny(log(y_asked),z_asked,pclass_sz);

 double result = dNdlnm*dlnmdlny;
 if (isinf(result))
  printf("inf in dndlny %.5e %.5e %.5e %.5e\n",m,get_dndlnM_at_z_and_M(z_asked,m,pclass_sz),get_volume_at_z(z_asked,pba),dlnmdlny);
 return result;

}


double get_theta_at_y_and_z(double y,
                            double z,
                            struct class_sz_structure * pclass_sz,
                            struct background * pba)
                         {
 if (z<pclass_sz->array_y_to_m_redshift[0])
    return -1.;//z = pclass_sz->array_y_to_m_redshift[0];
 if (z>pclass_sz->array_y_to_m_redshift[pclass_sz->n_z_y_to_m-1])
    return -1.;//z = pclass_sz->array_y_to_m_redshift[pclass_sz->n_z_y_to_m-1];

 if (y<pclass_sz->array_y_to_m_y[0])
    return -1.;//y = pclass_sz->array_y_to_m_y[0];

 if (y>pclass_sz->array_y_to_m_y[pclass_sz->n_y_y_to_m-1])
    return -1.;//y =  pclass_sz->array_y_to_m_y[pclass_sz->n_y_y_to_m-1];

double m = get_y_to_m_at_z_and_y(z,y,pclass_sz);
double theta = get_theta_at_m_and_z(m,z,pclass_sz,pba);
return theta;
                         }



double get_szcountsz_sigma_at_theta_in_patch(double theta,int index_patches,struct class_sz_structure *pclass_sz){

double y;
if (pclass_sz->use_skyaveraged_noise == 0){

    if (theta < pclass_sz->thetas[0]){
        int l1 = 0;
        int l2 = 1;
        double th1 = pclass_sz->thetas[l1];
        double th2 = pclass_sz->thetas[l2];
        double y1 = pclass_sz->ylims[index_patches][l1];
        double y2 = pclass_sz->ylims[index_patches][l2];
        y = y1 + (y2-y1)/(th2-th1)*(theta-th1);
        if (pclass_sz->use_edge_noise_values == 1){
           y = y1;
           }
        }
    else if (theta > pclass_sz->thetas[pclass_sz->nthetas-1]){
        int l1 = pclass_sz->nthetas - 1;
        int l2 = pclass_sz->nthetas - 2;
        double th1 = pclass_sz->thetas[l1];
        double th2 = pclass_sz->thetas[l2];
        double y1 = pclass_sz->ylims[index_patches][l1];
        double y2 = pclass_sz->ylims[index_patches][l2];
        y = y1 + (y2-y1)/(th2-th1)*(theta-th1);
        if (pclass_sz->use_edge_noise_values == 1){
           y = y2;
           }
        }

    else{
    y = pwl_value_1d(pclass_sz->nthetas,
                        pclass_sz->thetas,
                        pclass_sz->ylims[index_patches],
                        theta);
                      }
}
else{

    if (theta < pclass_sz->thetas[0]){
      int l1 = 0;
      int l2 = 1;
      double th1 = pclass_sz->thetas[l1];
      double th2 = pclass_sz->thetas[l2];
      double y1 = pclass_sz->sky_averaged_ylims[l1];
      double y2 = pclass_sz->sky_averaged_ylims[l2];
      y = y1 + (y2-y1)/(th2-th1)*(theta-th1);

      }
    else if (theta > pclass_sz->thetas[pclass_sz->nthetas-1]){
      int l1 = pclass_sz->nthetas - 1;
      int l2 = pclass_sz->nthetas - 2;
      double th1 = pclass_sz->thetas[l1];
      double th2 = pclass_sz->thetas[l2];
      double y1 = pclass_sz->sky_averaged_ylims[l1];
      double y2 = pclass_sz->sky_averaged_ylims[l2];
      y = y1 + (y2-y1)/(th2-th1)*(theta-th1);
    }
  else{
  // for sky averaged sigma's:
  y =  pwl_value_1d(pclass_sz->nthetas,
                       pclass_sz->thetas,
                       pclass_sz->sky_averaged_ylims,
                       theta); // ~5% difference
                     }
}

if (y<0) y = 0.;
if (isnan(y) || isinf(y) || (y==0)){
  // printf("nan or inf in get_szcountsz_sigma_at_theta_in_patch at theta  = %.5e and idpatch = %d\n",
  //       theta, index_patches);
  // printf("in this patch:\n");
  // int idth;
  // for (idth = 0; idth<pclass_sz->nthetas; idth++){
  //   printf("th = %.5e ylim = %.5e\n",pclass_sz->thetas[idth],pclass_sz->ylims[index_patches][idth]);
  // }
  y = 1e300;
}
return y;

}


// void tabulate_dlnm_dlnq(double * lnq_tab,
//                         );

double get_dlnm_dlnq(double lnq,
                     double z,
                     int idpatch,
                     struct class_sz_structure * pclass_sz,
                     struct background * pba)
                         {


// first we tabulate q as a function of m.
int ntab = 500;
double lnq_tab[ntab];
double lnm_tab[ntab];

int itab;
double lnm_tab_mmin = log(0.5*pclass_sz->M1SZ); // add some padding
double lnm_tab_mmax = log(2.*pclass_sz->M2SZ); // add some padding
double dlnm_tab = (lnm_tab_mmax-lnm_tab_mmin)/(ntab-1.);
for (itab = 0;itab<ntab;itab++){
lnm_tab[itab] = lnm_tab_mmin+itab*dlnm_tab;
double mtab = exp(lnm_tab[itab]);
double ytab = get_y_at_m_and_z(mtab,z,pclass_sz,pba);
double thetatab = get_theta_at_m_and_z(mtab,z,pclass_sz,pba);
double sigtab = get_szcountsz_sigma_at_theta_in_patch(thetatab,idpatch,pclass_sz);
lnq_tab[itab] = log(ytab/sigtab);
}

// now given the arrays we can compute the derivative by
// finite difference

double tol = 1e-3;
// double dlnq = 5*fabs(lnq_tab[ntab-1]-lnq_tab[0])/ntab;
double lnqp = lnq+tol;
double lnqm = lnq-tol;

double result;
double lnmp = pwl_value_1d(ntab,
                           lnq_tab,
                           lnm_tab,
                           lnqp);
double lnmm = pwl_value_1d(ntab,
                          lnq_tab,
                          lnm_tab,
                          lnqp);

result = (lnmp-lnmm)/2./tol;


// printf("dlnM\n");
//! JCH edit: I think Komatsu has forgotten the Jacobian factor dlnMdel/dlnM
//! as discussed in Eq. (5) of Komatsu-Seljak (2002)
//! Approximate via standard three-point finite difference
//! (checked w/ Mathematica implementation -- agrees very well)


return result;
}


double get_dlnm_dlny(double lny,
                     double z,
                     struct class_sz_structure * pclass_sz)
                         {
  // printf("dlnM\n");
//! JCH edit: I think Komatsu has forgotten the Jacobian factor dlnMdel/dlnM
//! as discussed in Eq. (5) of Komatsu-Seljak (2002)
//! Approximate via standard three-point finite difference
//! (checked w/ Mathematica implementation -- agrees very well)
double result;
double tol= 2.*(pclass_sz->array_y_to_m_y[1]-pclass_sz->array_y_to_m_y[0]);
//double mvir;
double fp,fm;
double lnyp = lny+tol;
double lnym = lny-tol;

fp = log(get_y_to_m_at_z_and_y(z,exp(lnyp),pclass_sz));
fm = log(get_y_to_m_at_z_and_y(z,exp(lnym),pclass_sz));

if (lny-tol<pclass_sz->array_y_to_m_y[0])
  return 0.;
if (lny+tol>pclass_sz->array_y_to_m_y[pclass_sz->n_y_y_to_m-1])
  return 0.;


result = (fp-fm)/2./tol;
// if (isinf(result)){}


// result = dlnm200ddlnm;


return result;
}



double get_y_to_m_at_z_and_y(double z_asked, double y_asked, struct class_sz_structure * pclass_sz)
{
  double z = log(1.+z_asked);
  double y = log(y_asked);

 if (z<pclass_sz->array_y_to_m_redshift[0])
    return 0.;//z = pclass_sz->array_y_to_m_redshift[0];
 if (z>pclass_sz->array_y_to_m_redshift[pclass_sz->n_z_y_to_m-1])
    return 0.;//z = pclass_sz->array_y_to_m_redshift[pclass_sz->n_z_y_to_m-1];

 if (y<pclass_sz->array_y_to_m_y[0])
    return 0.;//y = pclass_sz->array_y_to_m_y[0];

 if (y>pclass_sz->array_y_to_m_y[pclass_sz->n_y_y_to_m-1])
    return 0.;//y =  pclass_sz->array_y_to_m_y[pclass_sz->n_y_y_to_m-1];




 return pwl_interp_2d(pclass_sz->n_z_y_to_m,
                      pclass_sz->n_y_y_to_m,
                      pclass_sz->array_y_to_m_redshift,
                      pclass_sz->array_y_to_m_y,
                      pclass_sz->array_y_to_m_at_z_y,
                      1,
                      &z,
                      &y);

}


double get_m_to_xout_at_z_and_m(double z_asked, double m_asked, struct class_sz_structure * pclass_sz)
{
  double z = log(1.+z_asked);
  double m = log(m_asked);

 if (z<pclass_sz->array_m_to_xout_redshift[0])
    z = pclass_sz->array_m_to_xout_redshift[0];
 if (z>pclass_sz->array_m_to_xout_redshift[pclass_sz->n_z_m_to_xout-1])
    z = pclass_sz->array_m_to_xout_redshift[pclass_sz->n_z_m_to_xout-1];

 if (m<pclass_sz->array_m_to_xout_mass[0])
    m = pclass_sz->array_m_to_xout_mass[0];

 if (m>pclass_sz->array_m_to_xout_mass[pclass_sz->n_mass_m_to_xout-1])
    m =  pclass_sz->array_m_to_xout_mass[pclass_sz->n_mass_m_to_xout-1];




 return pwl_interp_2d(pclass_sz->n_z_m_to_xout,
                          pclass_sz->n_mass_m_to_xout,
                          pclass_sz->array_m_to_xout_redshift,
                          pclass_sz->array_m_to_xout_mass,
                          pclass_sz->array_m_to_xout_at_z_m,
                          1,
                          &z,
                          &m);

}


double get_dcib0dz_at_z_and_nu(double z_asked, double nu_asked, struct class_sz_structure * pclass_sz)
{
  double z = log(1.+z_asked);
  double m = log(nu_asked);

 if (z<pclass_sz->array_dcib0dz_redshift[0])
    z = pclass_sz->array_dcib0dz_redshift[0];
 if (z>pclass_sz->array_dcib0dz_redshift[pclass_sz->n_z_dcib0dz-1])
    z = pclass_sz->array_dcib0dz_redshift[pclass_sz->n_z_dcib0dz-1];

 if (m<pclass_sz->array_dcib0dz_nu[0])
    m = pclass_sz->array_dcib0dz_nu[0];

 if (m>pclass_sz->array_dcib0dz_nu[pclass_sz->n_nu_dcib0dz-1])
    m =  pclass_sz->array_dcib0dz_nu[pclass_sz->n_nu_dcib0dz-1];




 return exp(pwl_interp_2d(pclass_sz->n_z_dcib0dz,
                          pclass_sz->n_nu_dcib0dz,
                          pclass_sz->array_dcib0dz_redshift,
                          pclass_sz->array_dcib0dz_nu,
                          pclass_sz->array_dcib0dz_at_z_nu,
                          1,
                          &z,
                          &m));

}


double get_dndlnM_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);

 if (z<pclass_sz->array_z_dndlnM[0])
      return 1e-100;//z = pclass_sz->array_z_dndlnM[0];
 if (z>pclass_sz->array_z_dndlnM[pclass_sz->n_z_dndlnM-1])
      return 1e-100;//z = pclass_sz->array_z_dndlnM[pclass_sz->n_z_dndlnM-1];

 if (m<pclass_sz->array_m_dndlnM[0])
    return 1e-100;//m = pclass_sz->array_m_dndlnM[0];

 if (m>pclass_sz->array_m_dndlnM[pclass_sz->n_m_dndlnM-1])
      return 1e-100;//m =  pclass_sz->array_m_dndlnM[pclass_sz->n_m_dndlnM-1];


 return exp(pwl_interp_2d(pclass_sz->n_z_dndlnM,
                          pclass_sz->n_m_dndlnM,
                          pclass_sz->array_z_dndlnM,
                          pclass_sz->array_m_dndlnM,
                          pclass_sz->array_dndlnM_at_z_and_M,
                          1,
                          &z,
                          &m));
}

double get_m200m_to_m500c_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
  // printf("ok\n");


 return exp(pwl_interp_2d(pclass_sz->n_z_dndlnM,
                          pclass_sz->n_m_dndlnM,
                          pclass_sz->array_ln_1pz_m200m_to_m500c,
                          pclass_sz->array_m_m200m_to_m500c,
                          pclass_sz->array_m200m_to_m500c_at_z_and_M,
                          1,
                          &z,
                          &m));
// return 0.;
}


double get_m200m_to_m200c_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
 return exp(pwl_interp_2d(pclass_sz->n_z_dndlnM,
                          pclass_sz->n_m_dndlnM,
                          pclass_sz->array_ln_1pz_m200m_to_m200c,
                          pclass_sz->array_m_m200m_to_m200c,
                          pclass_sz->array_m200m_to_m200c_at_z_and_M,
                          1,
                          &z,
                          &m));
}


double get_m200c_to_m200m_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);

 if (z<pclass_sz->array_ln_1pz_m200c_to_m200m[0])
    z = pclass_sz->array_ln_1pz_m200c_to_m200m[0];
 if (z>pclass_sz->array_ln_1pz_m200c_to_m200m[pclass_sz->n_z_dndlnM-1])
    z = pclass_sz->array_ln_1pz_m200c_to_m200m[pclass_sz->n_z_dndlnM-1];

 if (m<pclass_sz->array_m_m200c_to_m200m[0])
    m = pclass_sz->array_m_m200c_to_m200m[0];

 if (m>pclass_sz->array_m_m200c_to_m200m[pclass_sz->n_m_dndlnM-1])
    m =  pclass_sz->array_m_m200c_to_m200m[pclass_sz->n_m_dndlnM-1];



 return exp(pwl_interp_2d(pclass_sz->n_z_dndlnM,
                          pclass_sz->n_m_dndlnM,
                          pclass_sz->array_ln_1pz_m200c_to_m200m,
                          pclass_sz->array_m_m200c_to_m200m,
                          pclass_sz->array_m200c_to_m200m_at_z_and_M,
                          1,
                          &z,
                          &m));
}


double get_m200c_to_m500c_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
 return exp(pwl_interp_2d(pclass_sz->n_z_dndlnM,
                          pclass_sz->n_m_dndlnM,
                          pclass_sz->array_ln_1pz_m200c_to_m500c,
                          pclass_sz->array_m_m200c_to_m500c,
                          pclass_sz->array_m200c_to_m500c_at_z_and_M,
                          1,
                          &z,
                          &m));
}


double get_m500c_to_m200c_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
 return exp(pwl_interp_2d(pclass_sz->n_z_dndlnM,
                          pclass_sz->n_m_dndlnM,
                          pclass_sz->array_ln_1pz_m500c_to_m200c,
                          pclass_sz->array_m_m500c_to_m200c,
                          pclass_sz->array_m500c_to_m200c_at_z_and_M,
                          1,
                          &z,
                          &m));
}


double get_m500c_to_m200m_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
 return exp(pwl_interp_2d(pclass_sz->n_z_dndlnM,
                          pclass_sz->n_m_dndlnM,
                          pclass_sz->array_ln_1pz_m500c_to_m200m,
                          pclass_sz->array_m_m500c_to_m200m,
                          pclass_sz->array_m500c_to_m200m_at_z_and_M,
                          1,
                          &z,
                          &m));
}


double get_m200c_to_mvir_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
 return exp(pwl_interp_2d(pclass_sz->n_z_dndlnM,
                          pclass_sz->n_m_dndlnM,
                          pclass_sz->array_ln_1pz_m200c_to_mvir,
                          pclass_sz->array_m_m200c_to_mvir,
                          pclass_sz->array_m200c_to_mvir_at_z_and_M,
                          1,
                          &z,
                          &m));
}



double get_m200m_to_mvir_at_z_and_M(double z_asked, double m_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
 return exp(pwl_interp_2d(pclass_sz->n_z_dndlnM,
                          pclass_sz->n_m_dndlnM,
                          pclass_sz->array_ln_1pz_m200m_to_mvir,
                          pclass_sz->array_m_m200m_to_mvir,
                          pclass_sz->array_m200m_to_mvir_at_z_and_M,
                          1,
                          &z,
                          &m));
}


double get_T10_alpha_at_z(double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  if (z<pclass_sz->T10_ln1pz[0])
   z = pclass_sz->T10_ln1pz[0];
  else if (z>pclass_sz->T10_ln1pz[pclass_sz->T10_lnalpha_size-1])
   z = pclass_sz->T10_ln1pz[pclass_sz->T10_lnalpha_size-1];

double result = exp(pwl_value_1d(pclass_sz->T10_lnalpha_size,
                          pclass_sz->T10_ln1pz,
                          pclass_sz->T10_lnalpha,
                          z));
// printf("z = %.3e  alpha = %.3e\n",z_asked,result);
return result;
}

double get_hmf_counter_term_nmin_at_z(double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  if (z<pclass_sz->array_redshift_hmf_counter_terms[0])
   z = pclass_sz->array_redshift_hmf_counter_terms[0];
  else if (z>pclass_sz->array_redshift_hmf_counter_terms[pclass_sz->n_z_hmf_counter_terms-1])
   z = pclass_sz->array_redshift_hmf_counter_terms[pclass_sz->n_z_hmf_counter_terms-1];

double result = pwl_value_1d( pclass_sz->n_z_hmf_counter_terms,
                              pclass_sz->array_redshift_hmf_counter_terms,
                              pclass_sz->array_hmf_counter_terms_nmin,
                              z);
//printf("z = %.3e  alpha = %.3e\n",z_asked,result);
return result;
}



double get_hmf_counter_term_b1min_at_z(double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  if (z<pclass_sz->array_redshift_hmf_counter_terms[0])
   z = pclass_sz->array_redshift_hmf_counter_terms[0];
  else if (z>pclass_sz->array_redshift_hmf_counter_terms[pclass_sz->n_z_hmf_counter_terms-1])
   z = pclass_sz->array_redshift_hmf_counter_terms[pclass_sz->n_z_hmf_counter_terms-1];

double result = pwl_value_1d( pclass_sz->n_z_hmf_counter_terms,
                              pclass_sz->array_redshift_hmf_counter_terms,
                              pclass_sz->array_hmf_counter_terms_b1min,
                              z);
//printf("z = %.3e  alpha = %.3e\n",z_asked,result);
return result;
}


double get_hmf_counter_term_b2min_at_z(double z_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  if (z<pclass_sz->array_redshift_hmf_counter_terms[0])
   z = pclass_sz->array_redshift_hmf_counter_terms[0];
  else if (z>pclass_sz->array_redshift_hmf_counter_terms[pclass_sz->n_z_hmf_counter_terms-1])
   z = pclass_sz->array_redshift_hmf_counter_terms[pclass_sz->n_z_hmf_counter_terms-1];

double result = pwl_value_1d( pclass_sz->n_z_hmf_counter_terms,
                              pclass_sz->array_redshift_hmf_counter_terms,
                              pclass_sz->array_hmf_counter_terms_b2min,
                              z);
//printf("z = %.3e  alpha = %.3e\n",z_asked,result);
return result;
}


double get_theta_at_m_and_z(double m, double z, struct class_sz_structure * pclass_sz, struct background * pba){

  double tau;
  int first_index_back = 0;
  double * pvecback;
  if (z==0)
    z += 1e-5; // avoid division by zero in theta

  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              pclass_sz->error_message);

  class_call(background_tau_of_z(pba,z,&tau),
             pclass_sz->error_message,
             pclass_sz->error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &first_index_back,
                               pvecback),
             pclass_sz->error_message,
             pclass_sz->error_message);

// Eq. 9 of https://arxiv.org/pdf/1303.5080.pdf
double H0 = pba->h*100.;
double Eh = pvecback[pba->index_bg_H]/pba->H0;
double d_A = pvecback[pba->index_bg_ang_distance]*pba->h;
double mp_bias = m/pclass_sz->HSEbias;
double thetastar2 = pclass_sz->thetastar * pow(H0/70.,-2./3.);
// below the pivot mass should always be 3e14*h, that fixes the theta_star normalisation
double theta500_for_mp_at_zp =  thetastar2 * pow(mp_bias/3.e14* (100./H0),pclass_sz->alpha_theta); // alpha_theta = 1/3
theta500_for_mp_at_zp *=    pow(Eh,-2./3)*pow(100.*d_A/(500.0*H0),-1.);
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


// double H0 = pba->h*100.;
double Eh = pvecback[pba->index_bg_H]/pba->H0;
double d_A = pvecback[pba->index_bg_ang_distance]*pba->h;
double rz = d_A*(1.+z);
double volume = _c_/1.0e5*rz*rz/Eh;
free(pvecback);

return volume;
}



double get_y_at_m_and_z(double m, double z, struct class_sz_structure * pclass_sz, struct background * pba){

  double tau;
  int first_index_back = 0;
  double * pvecback;
  if (z==0)
    z += 1e-5; // to avoid 1/0 division due to dA
  class_alloc(pvecback,
              pba->bg_size*sizeof(double),
              pclass_sz->error_message);

  class_call(background_tau_of_z(pba,z,&tau),
             pclass_sz->error_message,
             pclass_sz->error_message);

  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &first_index_back,
                               pvecback),
             pclass_sz->error_message,
             pclass_sz->error_message);


double H0 = pba->h*100.;
double Eh = pvecback[pba->index_bg_H]/pba->H0;
double d_A = pvecback[pba->index_bg_ang_distance]*pba->h;
double mp_bias = m/pclass_sz->HSEbias;
double yp;

if (pclass_sz->y_m_relation == 1){
        double A = pclass_sz->A_ym;
        double B = pclass_sz->B_ym;

        double f_rel;
        if (pclass_sz->apply_relativistic_correction_to_y_m == 0){
          f_rel = 1.;
        }
        else{
          double t = -0.00848*pow(mp_bias/(pclass_sz->m_pivot_ym*70./(pba->h*100.))*Eh,-0.585);
          f_rel = 1. + 3.79*t -28.2*t*t;
        }

        double h70 = 1.;//pba->h/0.70;

        yp = A*pow(Eh,2.)*pow(mp_bias/(pclass_sz->m_pivot_ym*pba->h)*h70,1.+B)*f_rel;

        // double a = -1.29389e-01;
        // double b = 9.991387e-01;
        // double c = -3.403211e+01;
        // double d = 1.279992e-01;
        // double e = 8.441768e-01;
        // double M200m_ws = mp_bias/pba->h;
        // yp = 1e-4*exp(a*exp(-pow(log(M200m_ws/3e14/e),2.))
        // +(b*(log(M200m_ws))+c)*(1.-exp(-pow(log(M200m_ws/3e14/d),2.))));
        // """Return optimization bias correction factor - multiply true y0 by this to get what the cluster finder recovers """



      }
else if (pclass_sz->y_m_relation == 0){

        double ystar2 = pow(10.,pclass_sz->ystar_ym)/pow(2., pclass_sz->alpha_ym)*0.00472724; // this factor is : 1./(5**2*1e8)*(np.pi/60/180)**-2 = 0.004727241144016912
        //  that last factor is because in the actual code they decided to normalzie D_A by 500 Mpc, and there is another 1e-4 normalization in Y500, both things combined lead to this 1/(500)**2/1e-4
        // justification: historical from https://www.aanda.org/articles/aa/pdf/2011/01/aa13999-10.pdf eq 5 ( arXiv:1001.0871)


        ystar2 *=  pow(H0/70.,-2.+pclass_sz->alpha_ym);
        double y500_for_mp_at_zp =  ystar2 * pow(mp_bias/pclass_sz->m_pivot_ym* (100./H0),pclass_sz->alpha_ym);
        y500_for_mp_at_zp *=   pow(Eh,pclass_sz->beta_ym) *pow(100.*d_A/(500.0*H0),-2.);


        yp = y500_for_mp_at_zp;


        // if (isinf(yp)){
          // printf("yp = %.5e m = %.5e z = %.5e Eh = %.5e d_A =%.5e\n",yp,mp_bias,z,Eh,d_A);
        //   // exit(0);
        // }

}

else if (pclass_sz->y_m_relation == 2){

        // double ystar2 = pow(10.,pclass_sz->ystar_ym)/pow(2., pclass_sz->alpha_ym)*0.00472724; ////8.9138435358806980e-004;
        //
        // ystar2 *=  pow(H0/70.,-2.+pclass_sz->alpha_ym);
        // double y500_for_mp_at_zp =  ystar2 * pow(mp_bias/pclass_sz->m_pivot_ym* (100./H0),pclass_sz->alpha_ym);
        // y500_for_mp_at_zp *=   pow(Eh,pclass_sz->beta_ym) *pow(100.*d_A/(500.0*H0),-2.);
        //
        //
        // yp = y500_for_mp_at_zp;
        // NIKA2



}

free(pvecback);

return yp;
}




double  get_L_sat_at_z_M_nu(double z_asked, double m_asked, double nu_asked, struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
  double nu = log(nu_asked);

  // printf("nu asked = %.3e\n",nu_asked);
  // exit(0);


  // double z = log(1.+z_asked);
  // double m = log(m_asked);
   if (z<pclass_sz->array_z_L_sat[0]){
      // z = pclass_sz->array_z_L_sat[0];
      printf("redshift min out of range in Lsat asked %.3e bound %.3e.\n",z,pclass_sz->array_z_L_sat[0]);
      exit(0);
    }
        // printf("dealing with mass conversion in hmf\n");
   if (z>pclass_sz->array_z_L_sat[pclass_sz->n_z_L_sat-1]){
      // z =  pclass_sz->array_z_L_sat[pclass_sz->n_z_L_sat-1];

      printf("redshift max out of range in Lsat asked %.3e bound %.3e.\n",z,pclass_sz->array_z_L_sat[pclass_sz->n_z_L_sat-1]);
      exit(0);
    }

   if (m<pclass_sz->array_m_L_sat[0]){
    // m = pclass_sz->array_m_L_sat[0];
      printf("mass min out of range in Lsat asked %.3e bound %.3e.\n",m,pclass_sz->array_m_L_sat[0]);
      exit(0);
  }
      // printf("dealing with mass conversion in hmf\n");
   if (m>pclass_sz->array_m_L_sat[pclass_sz->n_m_L_sat-1]){
      // m =  pclass_sz->array_m_L_sat[pclass_sz->n_m_L_sat-1];
      printf("mass max out of range in Lsat asked %.3e bound %.3e.\n",m,pclass_sz->array_m_L_sat[pclass_sz->n_m_L_sat-1]);
      exit(0);
    }

   if (nu<pclass_sz->array_nu_L_sat[0]){
    // nu = pclass_sz->array_nu_L_sat[0];
      printf("freq min out of range in Lsat asked %.8e bound %.8e.\n",exp(nu),exp(pclass_sz->array_nu_L_sat[0]));
      exit(0);
  }
      // printf("dealing with mass conversion in hmf\n");
   if (nu>pclass_sz->array_nu_L_sat[pclass_sz->n_nu_L_sat-1]){
      // nu =  pclass_sz->array_nu_L_sat[pclass_sz->n_nu_L_sat-1];
      printf("freq max out of range in Lsat asked %.3e bound %.3e.\n",exp(nu),exp(pclass_sz->array_nu_L_sat[pclass_sz->n_nu_L_sat-1]));
      exit(0);
    }

  // if (pclass_sz->tau_profile == 1){
  // find the closest l's in the grid:
  int id_l_low;
  int id_l_up;
  int n_nu = pclass_sz->n_nu_L_sat;
  int n_m = pclass_sz->n_m_L_sat;
  int n_z = pclass_sz->n_z_L_sat;
  r8vec_bracket(n_nu,pclass_sz->array_nu_L_sat,nu,&id_l_low,&id_l_up);

  // interpolate 2d at l_low:

 double ln_rho_low = pwl_interp_2d(n_m,
                                n_z,
                                pclass_sz->array_m_L_sat,
                                pclass_sz->array_z_L_sat,
                                pclass_sz->array_L_sat_at_M_z_nu[id_l_low-1],
                                1,
                                &m,
                                &z);

 double ln_rho_up = pwl_interp_2d(n_m,
                                n_z,
                                pclass_sz->array_m_L_sat,
                                pclass_sz->array_z_L_sat,
                                pclass_sz->array_L_sat_at_M_z_nu[id_l_up-1],
                                1,
                                &m,
                                &z);
 double ln_l_low = pclass_sz->array_nu_L_sat[id_l_low-1];
 double ln_l_up = pclass_sz->array_nu_L_sat[id_l_up-1];

//  printf("lnrh %.5e %.5e %d %d\n",ln_rho_low,ln_rho_up,id_l_low,id_l_up);
//  printf("nu %.5e %.5e %d %d %.5e\n",ln_l_low,ln_l_up,id_l_low,id_l_up,nu);
//  printf("expnu %.5e %.5e %d %d %.5e\n",exp(ln_l_low),exp(ln_l_up),id_l_low,id_l_up,exp(nu));

 return exp(ln_rho_low + ((nu - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low))-1.;
 // return ln_rho_low + ((l - ln_l_low) / (ln_l_up - ln_l_low)) * (ln_rho_up - ln_rho_low);


}



// // not used :
// double get_L_sat_at_z_and_M_at_nu(double z_asked,
//                                   double m_asked,
//                                   int index_nu,
//                                   struct background * pba,
//                                   struct class_sz_structure * pclass_sz){
//   double z = log(1.+z_asked);
//   double m = log(m_asked);
//    if (z<pclass_sz->array_z_L_sat[0])
//       z = pclass_sz->array_z_L_sat[0];
//         // printf("dealing with mass conversion in hmf\n");
//    if (z>pclass_sz->array_z_L_sat[pclass_sz->n_z_L_sat-1])
//       z =  pclass_sz->array_z_L_sat[pclass_sz->n_z_L_sat-1];
//    if (m<pclass_sz->array_m_L_sat[0])
//     m = pclass_sz->array_m_L_sat[0];
//       // printf("dealing with mass conversion in hmf\n");
//    if (m>pclass_sz->array_m_L_sat[pclass_sz->n_m_L_sat-1])
//       m =  pclass_sz->array_m_L_sat[pclass_sz->n_m_L_sat-1];
//  // printf("index_nu = %d\n",index_nu);
//  return exp(pwl_interp_2d(pclass_sz->n_z_L_sat,
//                           pclass_sz->n_m_L_sat,
//                           pclass_sz->array_z_L_sat,
//                           pclass_sz->array_m_L_sat,
//                           pclass_sz->array_L_sat_at_z_and_M_at_nu[index_nu],
//                           1,
//                           &z,
//                           &m))-1.;
// }

// double get_L_sat_at_z_and_M_at_nu_prime(double z_asked,
//                                   double m_asked,
//                                   struct background * pba,
//                                   struct class_sz_structure * pclass_sz){
//   double z = log(1.+z_asked);
//   double m = log(m_asked);
//  return exp(pwl_interp_2d(pclass_sz->n_z_L_sat,
//                           pclass_sz->n_m_L_sat,
//                           pclass_sz->array_z_L_sat,
//                           pclass_sz->array_m_L_sat,
//                           pclass_sz->array_L_sat_at_z_and_M_at_nu_prime,
//                           1,
//                           &z,
//                           &m))-1.;
// }
