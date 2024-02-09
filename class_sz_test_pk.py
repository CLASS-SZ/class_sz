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
'z_max_pk':20.

}

M_hmcode = Class()
M_hmcode.set(cosmo_params)
M_hmcode.set({
    'output': 'mPk',
    'non_linear':'hmcode'
})
start = time.perf_counter()
M_hmcode.compute()
end = time.perf_counter()
print('>>> class took %.3f s'%(end-start))

# M = Class()
Amod = -0.3
def compute_class_fast(M):
    
    M.set(cosmo_params)
    M.set({
        'output': 'mPk',#,lens_lens_1h,lens_lens_2h',
        # 'ndim_redshifts':30, # this is the number of redshift points used to tabulate Pk from the emulators. If you dont need Pk, set this to 4. 
        # the parameter ndim_redshifts is critical ! the more z points you add, the more calls to emulators
        # 'skip_background_and_thermo':1,
        # 'skip_chi': 0,
        # 'skip_hubble':0,

        # 'ell_max': 20000.0,
        # 'ell_min': 2.0,
        # 'dlogell': 0.1,
        # 'dell': 0,
        # 'redshift_epsrel': 0.0001,
        # 'mass_epsabs': 1e-40,
        # 'mass_epsrel': 0.0001,
            
        # 'M_min':1e7,
        # 'M_max':1e+17,
        # 'z_min':1e-5,
        # 'z_max': 5.,
            
        # 'k_min_for_pk_class_sz' : 1e-4,
        # 'k_max_for_pk_class_sz' : 5e1,
        # 'k_per_decade_class_sz' : 20.,
        # 'P_k_max_h/Mpc' : 200.0,

            
        # 'ndim_masses' : 150, # important 128 is default ccl value
        # 'ndim_redshifts' : 150,
        # 'non_linear':'hmcode',
        # # 'perturb_sampling_stepsize' : 0.005,
        # # 'k_max_tau0_over_l_max':5.,
            
        # 'hm_consistency': 1,

        'use_Amod': 1.,
        'Amod' : Amod,

    })
    M.compute_class_szfast()

M_fast = Class()
start = time.perf_counter()
compute_class_fast(M_fast)
end = time.perf_counter()
print('>>> class_szfast took %.3f s'%(end-start))


k1_a = np.geomspace(1e-3,10.,500)
h = M_fast.h()



label_size = 17
title_size = 22
legend_size = 13
handle_length = 1.5
fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(20,5))
ax = ax1
ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=label_size)
plt.setp(ax.get_xticklabels(), fontsize=label_size)
ax.grid( visible=True, which="both", alpha=0.1, linestyle='--')

ax.set_ylabel(r'$P(k)\quad\mathrm{[Mpc}/h]^3$',size=title_size)
ax.set_xlabel(r'$k\quad[h\mathrm{/Mpc]}$',size=title_size)

z = 0.
pk1_lin_fast = np.vectorize(M_fast.pk_lin)(k1_a,z)
pk1_lin_hmcode = np.vectorize(M_hmcode.pk_lin)(k1_a,z)
pk1_nonlin_fast = np.vectorize(M_fast.pk)(k1_a,z)
pk1_nonlin_hmcode = np.vectorize(M_hmcode.pk)(k1_a,z)


ax.plot(k1_a,pk1_lin_fast,label=r'fast (linear)',alpha=1.,c='k')
ax.plot(k1_a,pk1_nonlin_fast,label=r'fast (non-linear)',alpha=1.,c='k',lw=0.3)
ax.plot(k1_a,pk1_lin_hmcode,label=r'slow (linear)',alpha=1.,c='b',ls='-')
ax.plot(k1_a,pk1_nonlin_hmcode,label=r'hmcode',alpha=1.,c='b',ls='--')
ax.loglog()

ax.set_ylim(1e1,1e5)
ax.set_xlim(2e-3,1e1)

ax.legend(loc=1,frameon=True,framealpha=1,fontsize=11)

ax.set_title(r'$z=0$ Amod = %.3f'%Amod)

ax = ax2
ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=label_size)
plt.setp(ax.get_xticklabels(), fontsize=label_size)
ax.grid( visible=True, which="both", alpha=0.1, linestyle='--')

ax.set_ylabel(r'$P(k)\quad\mathrm{[Mpc}/h]^3$',size=title_size)
ax.set_xlabel(r'$k\quad[h\mathrm{/Mpc]}$',size=title_size)

z = 1.
pk1_lin_fast = np.vectorize(M_fast.pk_lin)(k1_a,z)
pk1_lin_hmcode = np.vectorize(M_hmcode.pk_lin)(k1_a,z)
pk1_nonlin_fast = np.vectorize(M_fast.pk)(k1_a,z)
pk1_nonlin_hmcode = np.vectorize(M_hmcode.pk)(k1_a,z)


ax.plot(k1_a,pk1_lin_fast,label=r'fast (linear)',alpha=1.,c='k')
ax.plot(k1_a,pk1_nonlin_fast,label=r'fast (non-linear)',alpha=1.,c='k',lw=0.3)
ax.plot(k1_a,pk1_lin_hmcode,label=r'slow (linear)',alpha=1.,c='b',ls='-')
ax.plot(k1_a,pk1_nonlin_hmcode,label=r'hmcode',alpha=1.,c='b',ls='--')
ax.loglog()

ax.set_ylim(1e1,1e5)
ax.set_xlim(2e-3,1e1)

ax.legend(loc=1,frameon=True,framealpha=1,fontsize=11)

ax.set_title(r'$z=1$ Amod = %.3f'%Amod)

ax = ax3
ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=label_size)
plt.setp(ax.get_xticklabels(), fontsize=label_size)
ax.grid( visible=True, which="both", alpha=0.1, linestyle='--')

ax.set_ylabel(r'$P(k)\quad\mathrm{[Mpc}/h]^3$',size=title_size)
ax.set_xlabel(r'$k\quad[h\mathrm{/Mpc]}$',size=title_size)

z = 3.
pk1_lin_fast = np.vectorize(M_fast.pk_lin)(k1_a,z)
pk1_lin_hmcode = np.vectorize(M_hmcode.pk_lin)(k1_a,z)
# if z>5.:
#     pk1_nonlin_fast = pk1_lin_fast
# else:
pk1_nonlin_fast = np.vectorize(M_fast.pk)(k1_a,z)
pk1_nonlin_hmcode = np.vectorize(M_hmcode.pk)(k1_a,z)


ax.plot(k1_a,pk1_lin_fast,label=r'fast (linear)',alpha=0.1,c='k',lw=3.3)
ax.plot(k1_a,pk1_lin_hmcode,label=r'slow (linear)',alpha=1.,c='b',ls='-')
ax.plot(k1_a,pk1_nonlin_fast,label=r'fast (non-linear)',alpha=0.1,c='r',lw=2.3)
ax.plot(k1_a,pk1_nonlin_hmcode,label=r'hmcode',alpha=1.,c='b',ls='--')
ax.loglog()

# ax.set_ylim(1e1,1e5)
ax.set_xlim(2e-3,1e1)

ax.legend(loc=1,frameon=True,framealpha=1,fontsize=11)

ax.set_title(r'$z=%f$ Amod = %.3f'%(z,Amod))




fig.tight_layout()
plt.show()
# plt.savefig('figures/pkz.pdf')
