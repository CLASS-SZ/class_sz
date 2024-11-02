import pytest
import scipy.integrate as integrate
import numpy as np
from classy_sz import Class as Class_sz
import os 


cosmo_params = {
'omega_b': 0.02242,
'omega_cdm':  0.11933,
'H0': 67.66, # use H0 because this is what is used by the emulators and to avoid any ambiguity when comparing with camb. 
'tau_reio': 0.0561,
'ln10^{10}A_s': 3.047,
'n_s': 0.9665   
}


classy_sz = Class_sz()
classy_sz.set(cosmo_params)
classy_sz.set({
'output':' ',
'skip_hubble':0,
'jax': 0,
})

print(classy_sz.pars)

classy_sz.compute_class_szfast()

print("classy_sz.jax_mode",classy_sz.jax_mode)



def test_classy_sz_hz_fast():
    z = np.linspace(0.,20,1000)
    assert classy_sz.get_hubble_at_z(z,params_values_dict = cosmo_params)[0] == pytest.approx(0.00022571238269609295,rel=1e-9)


test_classy_sz_hz_fast()