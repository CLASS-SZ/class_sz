from classy_sz import Class
import numpy as np

Omegam0 = 0.3075
H0 = 67.74
Omegab = 0.0486
Omegac = 0.2589
Omegac + Omegab
omega_b = Omegab*(H0/100.)**2
omega_c = Omegac*(H0/100.)**2
hparam = H0/100.
maniyar_cosmo = {
'omega_b': omega_b,
'omega_cdm':  omega_c,
'h': H0/100.,
# 'tau_reio': 0.0561,
'ln10^{10}A_s': 3.048,
'n_s': 0.9665,
# 'sigma8':0.830,
'k_pivot': 0.05,
'N_ncdm': 1,
'N_ur': 2.0328,
'm_ncdm': 0.0
}

M = Class()
M.set({'output':'dndlnM,cib_cib_1h,cib_cib_2h,cib_monopole'})
M.set(maniyar_cosmo)
M.set({

'mass function' : 'T08M200c',
'use_maniyar_cib_model':1,

'maniyar_cib_etamax' : 5.12572945e-01,

'maniyar_cib_zc' : 1.5,
'maniyar_cib_tau' : 8.25475287e-01,
'maniyar_cib_fsub' : 0.134*np.log(10.),
'Most_efficient_halo_mass_in_Msun' : 5.34372069e+12,
'Size_of_halo_masses_sourcing_CIB_emission' :  1.5583436676980493,
#for the Lsat tabulation:
'freq_min': 9e1,
'freq_max': 1.5e3,
'dlogfreq' : 0.1,

'concentration parameter':'fixed', # this sets it to 5

'n_z_L_sat' :100,
'n_m_L_sat' :100,
'n_nu_L_sat':500,

'use_nc_1_for_all_halos_cib_HOD': 1,

'sub_halo_mass_function' : 'TW10',#'JvdB14',
'M_min_subhalo_in_Msun' : 1e5, # 1e5 see https://github.com/abhimaniyar/halomodel_cib_tsz_cibxtsz/blob/master/Cell_cib.py
'use_redshift_dependent_M_min': 0,
'M_min' : 1e8*hparam,
'M_max' : 1e15*hparam,
'z_min' : 0.012,
'z_max' : 10.,
'ell_min': 10.,
'ell_max':5e4,
'dlogell':0.3,


'ndim_redshifts': 210,
'ndim_masses':150,

'has_cib_flux_cut': 0,
'hm_consistency':0,

'epsabs_L_sat': 1e-40,
'epsrel_L_sat': 1e-9,
    
'damping_1h_term':0,

# "P_k_max_1/Mpc": 50.,
# 'k_max_for_pk_class_sz':50.
})

M.set({
       'cib_frequency_list_num' : 2,
       'cib_frequency_list_in_GHz' : '100,240',
       'class_sz_verbose': 10,
      })
M.compute()

