import pytest
import scipy.integrate as integrate
import jax.numpy as jnp
import jax.scipy as jscipy
from classy_sz import Class as Class_sz
import os 
from jax import grad
import jax
import numpy as np


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
'skip_chi':0,
'jax': 1,
})



z = jnp.linspace(0.,10,10)


def test_classy_sz_daz_fast():

    
    classy_sz.compute_class_szfast()
    # print(classy_sz.get_angular_distance_at_z(z,params_values_dict = cosmo_params))

    expected_values = jnp.array([   0.          , 3663.16238962, 5619.19639645, 6811.47532722, 7626.15118295,
                               8226.71965179, 8692.09953795, 9066.6049671, 9375.83447402, 9637.42733764])
    
    result = classy_sz.get_angular_distance_at_z(z, params_values_dict=cosmo_params)*(1.+z)
    print("cosmo model default")
    print(result)

    
    np.testing.assert_allclose(result, expected_values, rtol=1e-5)

def test_classy_sz_daz_fast_cosmomodel_0():

    classy_sz.set({'cosmo_model': 0})
    classy_sz.compute_class_szfast()

    expected_values = jnp.array([0.          , 3663.04149234, 5618.94079371, 6811.03765429, 7625.75439281,
                               8226.01526502, 8691.41376217, 9065.71293446, 9375.23339903, 9636.58188782])
    
    result = classy_sz.get_angular_distance_at_z(z, params_values_dict=cosmo_params)*(1.+z)

    print("cosmo model 0")
    print(result)

    
    np.testing.assert_allclose(result, expected_values, rtol=1e-5)



