from classy_sz import Class
M = Class()
cosmo_params = {
'omega_b': 0.02242,
'omega_cdm':  0.11933,
'H0': 67.66, # use H0 because this is what is used by the emulators.
'tau_reio': 0.0561,
'ln10^{10}A_s': 3.047,
'n_s': 0.9665,

# 'k_pivot': 0.05,
# 'N_ncdm': 1,
# 'N_ur': 2.0328,
# 'm_ncdm': 0.06    

}

M = Class()
M.set(cosmo_params)
M.set({
'output' : 'gal_lens_1h,gal_lens_2h',#,galn_lens_1h,galn_lens_2h',
# 'output' : 'gal_lens_1h,gal_lens_2h',
# 'galaxy_samples_list_num' : 1, # the number of galaxy samples
# 'galaxy_samples_list' : '1', # the id string of each sample, can be any integer
# 'full_path_and_prefix_to_dndz_ngal':'/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/class_sz_auxiliary_files/WISC_bin3_ngal_example',
'full_path_to_dndz_gal':'/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/class-sz/class_sz_auxiliary_files/includes/WISC_bin3_ngal_example1.txt',
'galaxy_sample' : 'custom'
})
M.set({# class_sz parameters:

'mass function': 'T08M200c',
'concentration parameter' : 'D08',
'hm_consistency': 1,

'delta for galaxies' : '200c',
'delta for matter density' : '200c',
    
'x_out_truncated_nfw_profile' : 1,

# HOD parameters
'sigma_log10M_HOD' : 0.69,
'alpha_s_HOD' :  0.,
'M1_prime_HOD' : 5.03e12, # Msun/h
'M_min_HOD' : 6.25e11, # Msun/h
'M0_HOD' : 6.25e110,  # Msun/h ## set to very high value so Ns always 0.
'x_out_truncated_nfw_profile_satellite_galaxies' : 1.09,    
'f_cen_HOD' : 1., 
    
    
# # # HOD parameters
# 'sigma_log10M_HOD_ngal_0' : 0.69,
# 'alpha_s_HOD_ngal_0' :  0.,
# 'M1_prime_HOD_ngal_0' : 5.03e12, # Msun/h
# 'M_min_HOD_ngal_0' : 6.25e11, # Msun/h
# 'M0_HOD_ngal_0' : 6.25e11,  # Msun/h
# 'x_out_truncated_nfw_profile_satellite_galaxies_ngal_0' : 1.09,    
# 'f_cen_HOD_ngal_0' : 1., 
    

'M0 equal M_min (HOD)'  : 'no',



'M_min' : 1e11,
'M_max' : 5e15,
'ndim_masses' : 100,
    
    
'z_min' : 1e-3,
'z_max' : 3.,
'ndim_redshifts' : 100,
    
    
'dlogell' : 0.3,
'ell_max' : 20000.0,
'ell_min' : 2.0,


    
# precisions params:
'non_linear' : 'halofit',
'k_min_for_pk_class_sz' :  0.001,
'k_max_for_pk_class_sz' :  60.0,
'k_per_decade_class_sz' :  50,
'P_k_max_h/Mpc' :  50.0,

    
'redshift_epsabs' : 1e-40,
'redshift_epsrel' : 0.0001,
'mass_epsabs' : 1e-40,
'mass_epsrel' : 0.0001,    
'class_sz_verbose' : 10,

'cosmo_model' : 6,
        })
        
M.compute_class_szfast()
# M.compute()
cl_galn_lens = M.cl_galn_lens()
cl_gal_lens = M.cl_kg()

print(cl_gal_lens)
