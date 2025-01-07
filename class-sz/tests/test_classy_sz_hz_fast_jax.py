import pytest
import scipy.integrate as integrate
import jax.numpy as jnp
import jax.scipy as jscipy
from classy_sz import Class as Class_sz
import os 
from jax import grad
import jax


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
'jax': 1,
})

print(classy_sz.pars)

classy_sz.compute_class_szfast()

print("classy_sz.jax_mode",classy_sz.jax_mode)


z = 1.
def Hubble(H0):
    cosmo_params.update({'H0':H0})
    hz = classy_sz.get_hubble_at_z(z,params_values_dict = cosmo_params)
    return hz



def test_classy_sz_hz_fast():
    z = jnp.linspace(0.,20,1000)
    assert classy_sz.get_hubble_at_z(z,params_values_dict = cosmo_params)[0] == pytest.approx(0.00022571238269609295,rel=1e-9)
    assert classy_sz.get_hubble_at_z(1.,params_values_dict = cosmo_params) == pytest.approx(0.0004023684284735257,rel=1e-9)
    print("%.25f"%classy_sz.Hubble(1.)[0])
    assert classy_sz.Hubble(1.)[0] == pytest.approx(0.0004023684284735257101623,rel=1e-9)


def test_jaxification():

    z = jnp.linspace(1., 20, 1000)
    hubble_values = classy_sz.Hubble(z)

    # Test if it's a JAX array
    print("isinstance(hubble_values, jnp.ndarray)",isinstance(hubble_values, jnp.ndarray))
    assert isinstance(hubble_values, jnp.ndarray), "Hubble(z) should return a JAX array"

    # Test if it supports JAX transformations
    try:
        jitted_hubble = jax.jit(classy_sz.Hubble)(z)
        supports_jit = True
        print("Does Hubble(z) support JAX jit?", supports_jit)
        assert jnp.allclose(hubble_values, jitted_hubble), "Jitted Hubble should match original"
    except Exception as e:
        assert False, f"Hubble(z) should support JAX jit, but got error: {str(e)}"

def test_classy_sz_hz_fast_grad():
    dHubble = grad(Hubble)
    print(dHubble(76.))