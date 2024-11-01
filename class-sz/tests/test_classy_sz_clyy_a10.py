import pytest
import scipy.integrate as integrate
import numpy as np
from classy_sz import Class as Class_sz
import os 


cosmo_params= {
'omega_b': 0.02242,
'omega_cdm':  0.11933,
'H0': 67.66, 
'tau_reio': 0.0561,
'ln10^{10}A_s': 3.047,
'n_s': 0.9665, 
"cosmo_model": 1, # use mnu-lcdm emulators
}

precision_params = {
'x_outSZ': 4., # truncate profile beyond x_outSZ*r_s

'n_m_pressure_profile' :50, # default: 100, decrease for faster
'n_z_pressure_profile' :50, # default: 100, decrease for faster
    

'use_fft_for_profiles_transform' : 1, # use fft's or not. 
# only used if use_fft_for_profiles_transform set to 1
'N_samp_fftw' : 512,
'x_min_gas_pressure_fftw' : 1e-4,
'x_max_gas_pressure_fftw' : 1e4,
    
    
'ndim_redshifts' :30,

    
'redshift_epsabs': 1.0e-40,
'redshift_epsrel': 0.0001,    

    
'mass_epsabs': 1.0e-40,
'mass_epsrel': 0.0001
}

classy_sz = Class_sz()
classy_sz.set(cosmo_params)
classy_sz.set(precision_params)
classy_sz.set({

'output': 'tSZ_tSZ_1h,tSZ_tSZ_2h',
    
"ell_min" : 2,
"ell_max" : 8000,
'dell': 0,
'dlogell': 0.2,
    
'z_min' : 0.005,
'z_max' : 3.0,
'M_min' : 1.0e10, 
'M_max' : 3.5e15,
 

'mass_function' : 'T08M500c',



'pressure profile':'custom_gnfw', # can be Battaglia, Arnaud, etc
    
"P0GNFW": 8.130,
"c500": 1.156,
"gammaGNFW": 0.3292,
"alphaGNFW": 1.0620,
"betaGNFW":5.4807,
    



})
classy_sz.compute_class_szfast()

l = np.asarray(classy_sz.cl_sz()['ell'])
cl_yy_1h = np.asarray(classy_sz.cl_sz()['1h'])
cl_yy_2h = np.asarray(classy_sz.cl_sz()['2h'])

classy_sz_A10 = Class_sz()
classy_sz_A10.set(cosmo_params)
classy_sz_A10.set(precision_params)
classy_sz_A10.set({

'output': 'tSZ_1h,tSZ_2h',
    
"ell_min" : 2,
"ell_max" : 8000,
'dell': 0,
'dlogell': 0.2,
    
'z_min' : 0.005,
'z_max' : 3.0,
'M_min' : 1.0e10, 
'M_max' : 3.5e15,
 

'mass function' : 'T08M500c',

'pressure profile':'A10', # can be Battaglia, Arnaud, etc

})
classy_sz_A10.compute_class_szfast()

l_A10 = np.asarray(classy_sz_A10.cl_sz()['ell'])
cl_yy_1h_A10 = np.asarray(classy_sz_A10.cl_sz()['1h'])
cl_yy_2h_A10 = np.asarray(classy_sz_A10.cl_sz()['2h'])


def test_classy_sz_clyy_a10_1h():
    assert cl_yy_1h_A10.sum() == pytest.approx(33.73, rel=1e-2), f"Sum of '1h' values is {cl_yy_1h_A10.sum()}, expected close to 33.73"
    print(f"Test passed for a10 tabulated 1h: sum is {cl_yy_1h_A10.sum()}, close to 33.73")


def test_classy_sz_clyy_a10_2h():
    assert cl_yy_2h_A10.sum() == pytest.approx(0.550, rel=1e-2), f"Sum of '2h' values is {cl_yy_2h_A10.sum()}, expected close to 0.55"
    print(f"Test passed for a10 tabulated 2h: sum is {cl_yy_2h_A10.sum()}, close to 0.55")

def test_classy_sz_clyy_a10_ffts_1h():
    assert cl_yy_1h.sum() == pytest.approx(35.32, rel=1e-2), f"Sum of '1h' values is {cl_yy_1h.sum()}, expected close to 35.32"
    print(f"Test passed for a10 ffts 1h: sum is {cl_yy_1h.sum()}, close to 35.32")

def test_classy_sz_clyy_a10_ffts_2h():
    assert cl_yy_2h.sum() == pytest.approx(0.559, rel=1e-2), f"Sum of '2h' values is {cl_yy_2h.sum()}, expected close to 0.559"
    print(f"Test passed for a10 ffts 2h: sum is {cl_yy_2h.sum()}, close to 0.559")