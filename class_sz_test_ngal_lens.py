import numpy as np
import matplotlib.pyplot as plt
# import healpy as hp
from scipy.interpolate import interp1d
import sys
from classy_sz import Class
from mcfit.transforms import *
from scipy import interpolate
import time


#matplotlib.use('pdf')
font = {'size'   : 16, 'family':'STIXGeneral'}
plt.rcParams.update({
     "text.usetex": True,
     "font.family": "serif",
     "font.sans-serif": ['Computer Modern']})
plt.rc_context({'axes.autolimit_mode': 'round_numbers'})


def l_to_dl(lp):
    return lp*(lp+1.)/2./np.pi

# path_data = "/Users/aleksandra/Desktop/DES_tSZ/data/"
path_to_class_sz = "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/"

maglim_pdict={}
maglim_pdict['galaxy_sample']= 'custom'

# maglim_pdict['delta for galaxies'] = "200c"
# maglim_pdict['delta for matter density'] = "200c"
# maglim_pdict['mass function'] = 'T08M200c'
# maglim_pdict['concentration parameter'] = 'D08'

maglim_pdict['dlogell']= 0.3
#maglim_pdict['dell']= 50.0
maglim_pdict['ell_max']= 5000.0
maglim_pdict['ell_min']= 2.0
maglim_pdict['z_min']= 1.0e-8
maglim_pdict['z_max']= 2.0

### Precision
maglim_pdict['redshift_epsabs']= 1.0e-40
maglim_pdict['redshift_epsrel']= 0.0005
maglim_pdict['mass_epsabs']= 1.0e-40
maglim_pdict['mass_epsrel']= 0.0005
maglim_pdict['ndim_masses']= 150
maglim_pdict['ndim_redshifts']= 50
maglim_pdict['class_sz_verbose']= 10
maglim_pdict['nonlinear_verbose']= 0

maglim_pdict['use_fft_for_profiles_transform']= 1
maglim_pdict['N_samp_fftw'] = 1024 #precision parameter for the bessel transform to theta space
maglim_pdict['l_min_samp_fftw'] = 1e-12
maglim_pdict['l_max_samp_fftw']= 1e12
maglim_pdict['x_min_gas_pressure_fftw']= 1e-4
maglim_pdict['x_max_gas_pressure_fftw']= 1e3

maglim_pdict['P_k_max_h/Mpc']=100.0
maglim_pdict['k_min_for_pk_class_sz']= 0.0001
maglim_pdict['k_max_for_pk_class_sz']= 50.0
maglim_pdict['k_per_decade_class_sz']= 20.0

cosmo = {
    'Omega_b': 0.0486,
    'Omega_cdm': (0.315-0.0486),
    'h' : 0.6774,
    'tau_reio': 0.0543,
    #'ln10^{10}A_s': 2.7177747638974306,
    'sigma8': 0.8159,
    'n_s':  0.9649,
}

M = Class()

M.set(maglim_pdict)
M.set(cosmo)
M.set({
# 'output':'gal_lens_hf, mPk',
# 'output':'mPk',
'output':'gal_lens_hf',

# 'hm_consistency' : 1,
#
# 'M_min': 1e+10, # *cosmo['h'],
# 'M_max':  1e+15, # *cosmo['h'],

# [HOD]
#     logMmin = 11.74
#     sig_logM = 0.27
#     logM0 = 11.74
#     logM1 = 13.32
#     alpha_g = 1.66
#     fcen = 1.0
#     rmax_r200c = 1.0
#     rmax_rvir = 1.0
#     rsg_rs = 2.43
#     mass_def_for_rmax = '200c'
#HOD
# 'M0_HOD': 10.**11.74, # *cosmo['h'],
# 'M_min_HOD':10.**11.74, # *cosmo['h'], #Msun/h
# 'M1_prime_HOD':10.**13.32, #*cosmo['h'], #Msun/h
# 'sigma_log10M_HOD':0.27,
# 'alpha_s_HOD':1.66,
# 'f_cen_HOD': 1.,
# 'x_out_truncated_nfw_profile_satellite_galaxies':1., # so corresponds to 1xr200c
# 'csat_over_cdm' : 1.0,
'full_path_to_dndz_gal':  path_to_class_sz+'class_sz_auxiliary_files/nz_maglim_bin2.txt', # lens galaxies
# 'full_path_to_source_dndz_gal': path_to_class_sz+'class_sz_auxiliary_files/nz_source_normalized_bin4.txt', # source galaxies
# 'Delta_z_lens':0.00,
# 'Delta_z_source':0.00,


# 'use_pknl_in_2hterms': 0,
# 'P_k_max_h/Mpc':5e1,
# 'non_linear':'halofit',
'effective_galaxy_bias' : 1.3,
# 'skip_pkl':0,
# 'skip_pknl':0
'ndim_redshifts': 30,
'skip_background_and_thermo': 0,
'skip_chi': 1,
'skip_hubble': 1,
'skip_pkl': 0,
'skip_pknl': 0,
'skip_sigma8_and_der': 1,
'skip_sigma8_at_z': 1,
'cosmo_model': 0,
})

# M.compute_class_szfast()


t0 = time.time()
M.compute_class_szfast()
# M.compute()
t1 = time.time()

total = t1-t0

cl_gal_lens = M.cl_kg()
print('in time:',total,cl_gal_lens)

k1_a = np.geomspace(1e-3,10.,500)
z = 0.1
pk1 = np.vectorize(M.pk_lin)(k1_a,z)
print(pk1[:10])
