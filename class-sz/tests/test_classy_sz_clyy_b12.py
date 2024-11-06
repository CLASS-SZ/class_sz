import pytest
import scipy.integrate as integrate
import numpy as np
from classy_sz import Class as Class_sz
import os 


h = 0.674
omega_b = 0.0224
Omega_m=0.315
sigma8=0.811
N_ur = 3.046
n_s = 0.965
tau_reio=0.054
N_ncdm=0
m_ncdm=0


def test_classy_sz_clyy_b12():
    
    Mclass_sz = Class_sz()
    print('PATH_TO_CLASS_SZ_DATA:', os.environ['PATH_TO_CLASS_SZ_DATA'])

    Mclass_sz.set({

    'h': h,
    'sigma8': sigma8,
    'n_s': n_s,
    'tau_reio': tau_reio,
    'omega_b': omega_b,
    'omega_cdm':  0.1207,#Omega_m*h**2-omega_b,  
    'N_ur': 3.046,


    'output': 'tSZ_1h',

    'pressure_profile':'B12',
    'delta_for_electron_pressure':'200c',
    "concentration_parameter":"D08",
    "ell_min" : 125,
    "ell_max" : 9725,
    'dell': 200,
    'dlogell': 0.,
        
    'M_min' : 1e11*h, 
    'M_max' : 1e16*h,

    'z_min': 0.005,
    'z_max': 6.,
        
        
    'n_z_pressure_profile': 500,
    'n_m_pressure_profile' : 500,
    'n_l_pressure_profile' : 500,
        
    'l_min_gas_pressure_profile' :  1.e-2,
    'l_max_gas_pressure_profile' :  5.e4,    

    'pressure_profile_epsrel':1e-4,
    'pressure_profile_epsabs':1e-100,
        

        
        
    'hm_consistency' : 0,
        

    'use_fft_for_profiles_transform' : 1,
    'x_min_gas_pressure_fftw' : 1e-5,
    'x_max_gas_pressure_fftw' : 1e5,
    'N_samp_fftw' : 8192,
        
        
    # 'ndim_masses' : 500,
    # 'ndim_redshifts' :100,
    'redshift_epsrel': 1e-6,
    'redshift_epsabs': 1e-100,
    'mass_epsrel':1e-6,
    'mass_epsabs':1e-100,    
        

    'truncate_gas_pressure_wrt_rvir' : 1,
    'x_outSZ': 2.,
    'mass_function' : 'T10M200m',
    'T10_alpha_fixed' : 1,
        
        
    'P_k_max_h/Mpc': 10.,
    'k_per_decade_class_sz':80.,
    'k_min_for_pk_class_sz':1e-4,
    'k_max_for_pk_class_sz': 10*h, 
        

    })
    Mclass_sz.compute()


    result_sum = np.asarray(Mclass_sz.cl_sz()['1h']).sum()

    assert result_sum == pytest.approx(51.6, rel=1e-2), f"Sum of '1h' values is {result_sum}, expected close to 51.6"
    # Print message if the test passes
    print(f"Test passed: sum is {result_sum}, close to 51.6")



