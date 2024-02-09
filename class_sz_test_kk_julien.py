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
cosmo_params = {
        'omega_b': 0.02242,
        'omega_cdm':  0.11933,
        'H0': 67.66, # use H0 because this is what is used by the emulators.
        'tau_reio': 0.0561,
        # 'sigma8': 0.81,
        # 'sigma8': 0.81,
        'ln10^{10}A_s': 3.0980,
        'n_s': 0.9665,

        'k_pivot': 0.05,
        'N_ncdm': 1,
        'N_ur': 2.0328,
        'm_ncdm': 0.06,
        # 'z_max_pk':20.

        'output': 'lens_lens_hf',#,lens_lens_1h,lens_lens_2h',
        # 'ndim_redshifts':30, # this is the number of redshift points used to tabulate Pk from the emulators. If you dont need Pk, set this to 4. 
        # the parameter ndim_redshifts is critical ! the more z points you add, the more calls to emulators
        # 'skip_background_and_thermo':1,
        # 'skip_chi': 0,
        # 'skip_hubble':0,

        'ell_max': 10000.0,
        'ell_min': 2.0,
        'dlogell': 0.1,
        'dell': 0,
        'redshift_epsrel': 0.00001,

            
        'z_min':1e-5,
        'z_max': 25.,
            
        # 'k_min_for_pk_class_sz' : 1e-4,
        # 'k_max_for_pk_class_sz' : 5e1,
        # 'k_per_decade_class_sz' : 20.,
        # 'P_k_max_h/Mpc' : 200.0,

            
        # 'ndim_masses' : 150, # important 128 is default ccl value
        'ndim_redshifts' : 100,
        'non_linear':'hmcode',
        # 'perturb_sampling_stepsize' : 0.005,
        # 'k_max_tau0_over_l_max':5.,
        # 'lensing':'yes',
        # 'modes':'s',
        # 'skip_background_and_thermo': 0,
        # 'skip_cmb': 1,

}


# M = Class()
def compute_class_fast(M):
    
    M.set(cosmo_params)
    M.set({

    

    })
    M.compute_class_szfast()

# M = Class()
def compute_class_slow(M):
    
    M.set(cosmo_params)
    M.set({
    
    'output': 'lCl,tCl,mPk',
    # 'perturb_sampling_stepsize' : 0.005,
    'k_max_tau0_over_l_max':0.01,
    'lensing':'yes',
    'modes':'s',
    'l_max_scalars':10000.,
    'l_logstep': 0.3,


    'overwrite_clpp_with_limber' : 0,
    'z_max_pk': 30.,
    'P_k_max_h/Mpc': 100.,
    'z_min': 1e-5,
    'z_max': 30.,
    'class_sz_verbose':0,
    'ell_min': 2.,
    'ell_max': 20000,
    'dell': 0.,
    'dlogell': 0.1
    })
    M.compute()

M_fast = Class()
start = time.perf_counter()
compute_class_fast(M_fast)
end = time.perf_counter()
print('>>> class_szfast took %.3f s'%(end-start))


cl_kk_fast = M_fast.cl_kk
ell_fast = np.asarray(cl_kk_fast()['ell'])
fac_fast = ell_fast*(ell_fast+1.)/2./np.pi
cl_kk_hf_fast = np.asarray(cl_kk_fast()['hf'])


l_class = M_fast.lensed_cl()['ell']
cl_class = M_fast.lensed_cl()['pp']*(M_fast.lensed_cl()['ell']*(M_fast.lensed_cl()['ell']+1.)/2.)**2.
# print(cl_kk_1h_fast,cl_kk_2h_fast)
# exit(0)

M_slow = Class()
start = time.perf_counter()
compute_class_slow(M_slow)
end = time.perf_counter()
print('>>> class_szslow took %.3f s'%(end-start))

# cl_kk_slow = M_slow.cl_kk
# ell_slow = np.asarray(cl_kk_slow()['ell'])
# fac_slow = ell_slow*(ell_slow+1.)/2./np.pi
# cl_kk_1h_slow = np.asarray(cl_kk_slow()['1h'])
# cl_kk_2h_slow = np.asarray(cl_kk_slow()['2h'])
l_class_full = M_slow.lensed_cl()['ell']
cl_class_full = M_slow.lensed_cl()['pp']*(M_slow.lensed_cl()['ell']*(M_slow.lensed_cl()['ell']+1.)/2.)**2.





label_size = 17
title_size = 22
legend_size = 13
handle_length = 1.5
fig, (ax1) = plt.subplots(1,1,figsize=(10,5))

ax = ax1
ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=label_size)
plt.setp(ax.get_xticklabels(), fontsize=label_size)
ax.grid( visible=True, which="both", alpha=0.1, linestyle='--')

ax.set_ylabel(r'$C_\ell$',size=title_size)
ax.set_xlabel(r'$\ell$',size=title_size)



ax.plot(ell_fast,cl_kk_hf_fast/fac_fast,ls='-',c='b',label=r'class hf')
ax.plot(l_class,cl_class,ls='--',c='k',label=r'class emu')
ax.plot(l_class_full,cl_class_full,ls='-.',c='r',label=r'class full')

# l,cl = np.loadtxt("/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/output/test_cpp.txt",unpack=True)
# cl *= (l*(l+1.)/2.)**2.
# ax.plot(l,cl,ls=':',c='g',label=r'class printed')




ax.loglog()

# ax.set_ylim(1e1,1e5)
# ax.set_xlim(2e-3,1e1)

ax.legend(loc=1,frameon=True,framealpha=1,fontsize=11)

# ax.set_title(r'$z=%f$'%z)




fig.tight_layout()
plt.show()
# plt.savefig('figures/pkz.pdf')
