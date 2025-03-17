import matplotlib
import matplotlib.pyplot as plt
import os 
os.environ["JAX_PLATFORM_NAME"] = "cpu"
from classy_sz import Class as Class_sz

import jax
import jax.numpy as jnp


cosmo_params = {
'omega_b': 0.02242,
'omega_cdm':  0.11933,
'H0': 67.66, # use H0 because this is what is used by the emulators and to avoid any ambiguity when comparing with camb.
'tau_reio': 0.0561,
'ln10^{10}A_s': 3.047,
'n_s': 0.9665
}

# initialize computation
classy_sz = Class_sz()
classy_sz.set(cosmo_params)
classy_sz.set({
'jax': 1
})
classy_sz.initialize_classy_szfast()


z = 0.3
pks,ks = classy_sz.get_pkl_at_z(z,params_values_dict = cosmo_params)

plt.plot(ks,pks)
cosmo_params.update({'H0':82,'n_s':0.75})
pks,ks = classy_sz.get_pkl_at_z(z,params_values_dict = cosmo_params)
plt.plot(ks,pks)
plt.loglog()
plt.show()


def get_pkl_at_z_and_k(z, k, cosmo_params):
    print('getting k,pk')
    pks, ks = classy_sz.get_pkl_at_z(z, params_values_dict=cosmo_params)

    print('getting lnk,lnpk')
    # Convert to JAX arrays if not already
    pks = jnp.array(pks)
    ks = jnp.array(ks)

    # Interpolate in log-log space
    log_ks = jnp.log(ks)
    log_pks = jnp.log(pks)

    # Convert k to an array with shape (1,) to avoid issues with scalar interpolation
    log_k = jnp.log(jnp.array([k]))

    print('interp lnk,lnpk', log_k)
    print(log_ks.shape, log_pks.shape, log_k.shape)
    log_pk = jnp.interp(log_k, log_ks, log_pks)

    # print('return lnk,lnpk')
    # Extract scalar from the result if needed
    return jnp.exp(log_pk)[0]

# Example usage
k = 0.1  # Example k value
z = 1.5
H0 = 69.
cosmo_params.update({'H0': H0})
pkl = get_pkl_at_z_and_k(z, k, cosmo_params)
print('pkl',pkl)
