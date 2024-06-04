# include "class_sz.h"
# include "class_sz_tools.h"
# include "class_sz_custom_bias.h"
# include "Patterson.h"
# include "r8lib.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "fft.h"


double get_b_custom1_at_z(double z,
                          struct class_sz_structure * pclass_sz){

double ln1pz = log(1.+z);


return exp(pwl_value_1d(pclass_sz->array_b_custom1_n_z,
                        pclass_sz->array_b_custom1_ln1pz,
                        pclass_sz->array_b_custom1_bias,
                        ln1pz));
                                           }
