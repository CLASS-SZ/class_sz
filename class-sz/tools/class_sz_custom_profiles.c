# include "class_sz.h"
# include "class_sz_tools.h"
# include "class_sz_custom_profiles.h"
# include "Patterson.h"
# include "r8lib.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "fft.h"


double get_radial_kernel_W_custom1_at_z(double z,
                                        struct class_sz_structure * pclass_sz){

double ln1pz = log(1.+z);
return exp(pwl_value_1d(pclass_sz->array_custom1_redshift_kernel_n_z,
                        pclass_sz->array_custom1_redshift_kernel_ln1pz,
                        pclass_sz->array_custom1_redshift_kernel_W,
                        ln1pz));
                                           }


double get_custom1_profile_at_k_m_z(double k_asked,
                                    double m_asked,
                                    double z_asked,
                                    struct class_sz_structure * pclass_sz){
  double z = log(1.+z_asked);
  double m = log(m_asked);
  double k = log(k_asked);

  int id_k_low;
  int id_k_up;
  int n_k = pclass_sz->n_k_custom1_profile;
  int n_m = pclass_sz->n_m_custom1_profile;
  int n_z = pclass_sz->n_z_custom1_profile;
  r8vec_bracket(n_k,pclass_sz->array_custom1_profile_ln_k,k,&id_k_low,&id_k_up);

  if (id_k_low == id_k_up){
    printf("bug in get_custom1_profile_at_k_m_z");
    exit(0);
  }

  if (m<pclass_sz->array_custom1_profile_ln_m[0])
    return 0.;
  if (m>pclass_sz->array_custom1_profile_ln_m[n_m-1])
    return 0.;

  if (k<pclass_sz->array_custom1_profile_ln_k[0])
    return 0.;
  if (k>pclass_sz->array_custom1_profile_ln_k[n_k-1])
    return 0.;

  // interpolate 2d at l_low:

 double ln_rho_low = pwl_interp_2d(n_m,
                                n_z,
                                pclass_sz->array_custom1_profile_ln_m,
                                pclass_sz->array_custom1_profile_ln_1pz,
                                pclass_sz->array_custom1_profile_u_at_lnk_lnm_ln1pz[id_k_low-1],
                                1,
                                &m,
                                &z);

 double ln_rho_up = pwl_interp_2d(n_m,
                                n_z,
                                pclass_sz->array_custom1_profile_ln_m,
                                pclass_sz->array_custom1_profile_ln_1pz,
                                pclass_sz->array_custom1_profile_u_at_lnk_lnm_ln1pz[id_k_up-1],
                                1,
                                &m,
                                &z);
 double ln_k_low = pclass_sz->array_custom1_profile_ln_k[id_k_low-1];
 double ln_k_up = pclass_sz->array_custom1_profile_ln_k[id_k_up-1];

 return exp(ln_rho_low + ((k - ln_k_low) / (ln_k_up - ln_k_low)) * (ln_rho_up - ln_rho_low));
}

double get_custom1_profile_at_x_m_z(double x_asked,
                                    double m_asked,
                                    double z_asked,
                                    struct class_sz_structure * pclass_sz)
{
  double z = log(1.+z_asked);
  double m = log(m_asked);
  double k = log(x_asked);

  int id_k_low;
  int id_k_up;
  int n_k = pclass_sz->n_k_custom1_profile;
  int n_m = pclass_sz->n_m_custom1_profile;
  int n_z = pclass_sz->n_z_custom1_profile;
  r8vec_bracket(n_k,pclass_sz->array_custom1_profile_ln_x,k,&id_k_low,&id_k_up);

  if (id_k_low == id_k_up){
    printf("bug in get_custom1_profile_at_k_m_z");
    exit(0);
  }

  if (m<pclass_sz->array_custom1_profile_ln_m[0])
    return 0.;
  if (m>pclass_sz->array_custom1_profile_ln_m[n_m-1])
    return 0.;
  if (k<pclass_sz->array_custom1_profile_ln_x[0])
    return 0.;
  if (k>pclass_sz->array_custom1_profile_ln_x[n_k-1])
    return 0.;
  // interpolate 2d at l_low:

 double ln_rho_low = pwl_interp_2d(n_m,
                                n_z,
                                pclass_sz->array_custom1_profile_ln_m,
                                pclass_sz->array_custom1_profile_ln_1pz,
                                pclass_sz->array_custom1_profile_u_at_lnx_lnm_ln1pz[id_k_low-1],
                                1,
                                &m,
                                &z);

 double ln_rho_up = pwl_interp_2d(n_m,
                                n_z,
                                pclass_sz->array_custom1_profile_ln_m,
                                pclass_sz->array_custom1_profile_ln_1pz,
                                pclass_sz->array_custom1_profile_u_at_lnx_lnm_ln1pz[id_k_up-1],
                                1,
                                &m,
                                &z);
 double ln_k_low = pclass_sz->array_custom1_profile_ln_x[id_k_low-1];
 double ln_k_up = pclass_sz->array_custom1_profile_ln_x[id_k_up-1];

 return exp(ln_rho_low + ((k - ln_k_low) / (ln_k_up - ln_k_low)) * (ln_rho_up - ln_rho_low));
}
