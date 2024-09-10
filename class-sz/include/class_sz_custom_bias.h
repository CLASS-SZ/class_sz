
#ifndef __CUSTOM_BIAS__
#define __CUSTOM_BIAS__

#include "common.h"
#include "r8lib.h"
#include <time.h>
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf_bessel.h"

double get_b_custom1_at_z(double z,
                          struct class_sz_structure * pclass_sz);


#endif
