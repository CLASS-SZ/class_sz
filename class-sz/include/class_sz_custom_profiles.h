
#ifndef __CUSTOM_PROFILES__
#define __CUSTOM_PROFILES__

#include "common.h"
#include "r8lib.h"
#include <time.h>
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf_bessel.h"

double get_radial_kernel_W_custom1_at_z(double z,
                                        struct class_sz_structure * pclass_sz);

double get_custom1_profile_at_k_m_z(double k_asked,
                                    double m_asked,
                                    double z_asked,
                                    struct class_sz_structure * pclass_sz);

double get_custom1_profile_at_x_m_z(double x_asked,
                                    double m_asked,
                                    double z_asked,
                                    struct class_sz_structure * pclass_sz);

// double get_b_custom1_at_z(double z,
//                           struct class_sz_structure * pclass_sz);


#endif
