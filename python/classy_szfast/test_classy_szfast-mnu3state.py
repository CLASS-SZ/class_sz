# import classy_szfast as cszfast
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# from classy_szfast import tabulate_gas_pressure_profile_k,get_gas_pressure_profile_x_parallel
from classy_sz import Class
from classy_szfast import classy_szfast
import time

M = Class()


params = {
'output': 'tCl pCl lCl mPk',
'lensing': 'yes',
    
'omega_b': 0.02242,
'omega_cdm': 0.11933,
'H0': 67.66,
'tau_reio': 0.0561,
'ln10^{10}A_s': 3.047,
'n_s': 0.9665,
'k_pivot': 0.05,

'm_ncdm': '0.02,0.2,0.2',
'N_ncdm': 3,
'N_ur': 0.00641,



'l_max_scalars': 8000, #specified in the ACT likelihood Cobaya interface itself as well, but repeat here just in case
'perturb_sampling_stepsize': 0.05,

          }

# def compute_class_sz(M):
M.set(params)
M.set({
'output':'tCl,lCl,pCl',
'skip_background_and_thermo': 1,
'skip_pkl': 1,
'skip_pknl': 1,
'skip_sigma8_and_der': 1,
'skip_sigma8_at_z': 1,
'skip_chi': 1,
'skip_hubble': 1,
'cosmo_model' : 5,
    
'm_ncdm': '0.02,0.02,0.02',
'N_ncdm': 3,
'N_ur': 0.00641,
'modes' : 's',
})
M.compute_class_szfast()

cls_mnu_3states = M.lensed_cl(lmax=8000)
