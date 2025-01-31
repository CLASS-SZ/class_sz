#include "class_sz_macros.h"
class_sz_string_parameter(root,"root","root")
class_sz_string_parameter(A10_file,"/class_sz_auxiliary_files/includes/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_A10.txt","A10 file")
class_sz_string_parameter(P13_file,"/class_sz_auxiliary_files/includes/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_P13.txt","P13 file")
class_sz_string_parameter(Tinker_et_al_10_alpha_consistency_msyriac_file,"/class_sz_auxiliary_files/includes/Tinker_et_al_10_alpha_consistency_msyriac.txt","Tinker_et_al_10_alpha_consistency_msyriac_file")

class_sz_pclass_sz_parameter(jax,int,0)


class_sz_string_parameter(append_name_trispectrum_ref,"append_name_trispectrum_ref","append_name_trispectrum_ref")
class_sz_string_parameter(path_to_ref_trispectrum_for_cobaya,"path_to_ref_trispectrum_for_cobaya","path_to_ref_trispectrum_for_cobaya")
class_sz_string_parameter(full_path_to_noise_curve_for_y_y,"full_path_to_noise_curve_for_y_y","full_path_to_noise_curve_for_y_y")
class_sz_string_parameter(full_path_to_noise_curve_for_t_t,"full_path_to_noise_curve_for_t_t","full_path_to_noise_curve_for_t_t")
class_sz_string_parameter(full_path_to_dndz_gal,"/class_sz_auxiliary_files/includes/WISC_bin3_ngal_example1.txt","custom dndz file")
class_sz_string_parameter(full_path_and_prefix_to_dndz_ngal,"/class_sz_auxiliary_files/includes/WISC_bin3_ngal_example","custom dndz file")
class_sz_string_parameter(full_path_to_source_dndz_gal,"/class_sz_auxiliary_files/includes/WISC_bin3.txt","custom dndz file")
class_sz_string_parameter(full_path_to_redshift_dependent_M_min,"/class_sz_auxiliary_files/websky_halo_mass_completion_z_Mmin_in_Msun_over_h.txt","custom M_min file")


class_sz_string_parameter(full_path_to_n5k_gg_chi_K0,"/class_sz_auxiliary_files/excludes/galaxies/n5k_gg_chi_K0.txt","n5k_gg_chi_K0")
class_sz_string_parameter(full_path_to_n5k_z_chi,"/class_sz_auxiliary_files/excludes/galaxies/n5k_z_chi.txt","n5k_z_chi")
class_sz_string_parameter(full_path_to_n5k_z,"/class_sz_auxiliary_files/excludes/galaxies/n5k_z.txt","n5k_z")
class_sz_string_parameter(full_path_to_n5k_k,"/class_sz_auxiliary_files/excludes/galaxies/n5k_k.txt","n5k_k")
class_sz_string_parameter(full_path_to_n5k_pk_nl,"/class_sz_auxiliary_files/excludes/galaxies/n5k_pk_nl.txt","n5k_pk_nl")


class_sz_string_parameter(cmb_cls_filename,"/class_sz_auxiliary_files/cmb_cls_test.pickle","cmb_cls_filename")


//class_sz_string_parameter(UNWISE_dndz_file,"/class_sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz.txt","unWISE dndz file")
class_sz_string_parameter(UNWISE_dndz_file,"/class_sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz_cosmos.txt","unWISE dndz file")
class_sz_string_parameter(UNWISE_fdndz_file,"/class_sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_fdndz.txt","unWISE fdndz file")
class_sz_string_parameter(UNWISE_cosmos_dndz_file,"/class_sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz_cosmos.txt","unWISE fdndz file")
class_sz_string_parameter(WISC3_dndz_file,"/class_sz_auxiliary_files/WISC_bin3.txt","WISC bin 3 file")
class_sz_string_parameter(cib_Snu_file_snu,"/class_sz_auxiliary_files/includes/filtered_snu_planck_fine.txt","cib_Snu_file_snu")
class_sz_string_parameter(cib_Snu_file_z,"/class_sz_auxiliary_files/includes/filtered_snu_planck_z_fine.txt","cib_Snu_file_z")
class_sz_string_parameter(cib_Snu_file_nu,"/class_sz_auxiliary_files/includes/filtered_snu_planck_nu_fine.txt","cib_Snu_file_nu")


class_sz_string_parameter(ksz_filter_file,"/class_sz_auxiliary_files/UNWISE_galaxy_distributions/unwise_filter_functions_l_fl.txt","ksz filter file")
class_sz_string_parameter(ksz_template_file,"/class_sz_auxiliary_files/cl_ksz_bat.dat","ksz template file")
class_sz_string_parameter(ksz_reio_template_file,"/class_sz_auxiliary_files/FBN_kSZ_PS_patchy.d.txt","ksz template file, reio contribution")


class_sz_string_parameter(cmb_lensing_noise_file,"cmb_lensing_noise_file","cmb_lensing_noise_file")
//in this file:
//from Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/class_sz_auxiliary_files/noise_curves/nlkk_v3_1_0_deproj0_SENS2_fsky0p4_it_lT30-3000_lP30-5000.dat
// column ells and N_lensing_MV (all)
// columns in the original file
//[ells, N_lensing_TT, N_lensing_TE, N_lensing_EE, N_lensing_TB, N_lensing_EB, N_lensing_Pol (EE+EB), N_lensing_MV (all), N_curl_TT, N_curl_TE, N_curl_EE, N_curl_TB, N_curl_EB, N_curl_Pol (EE+EB), N_curl_MV (all)]

// other option for s4:
// from https://cmb-s4.uchicago.edu/wiki/index.php/Survey_Performance_Expectations
// kappa_deproj0_sens0_16000_lT30-3000_lP30-5000.dat
// columns in the original file are :
// ell, TT, TE, EE, TB, EB, EE+EB, TT+TE+EE+EB, TT (curl), TE (curl), EE (curl), TB (curl), EB (curl), EE+EB (curl), TT+TE+EE+EB (curl)

class_sz_string_parameter(SZ_cat_file,"/class_sz_auxiliary_files/includes/SZ_cat.txt","SZ_cat_file")


class_sz_string_parameter(classy_sz_verbose,"none","classy_sz_verbose")


class_sz_string_parameter(Planck_thetas_file,"/class_sz_auxiliary_files/includes/SZ_thetas.txt","sz_selection_function_thetas_file")
class_sz_string_parameter(SO_thetas_file,"/class_sz_auxiliary_files/nemo_sim_thetas_030722_50bins.txt","sz_selection_function_thetas_file")

class_sz_string_parameter(Planck_skyfracs_file,"/class_sz_auxiliary_files/includes/SZ_skyfracs.txt","sz_selection_function_skyfracs_file")
class_sz_string_parameter(SO_skyfracs_file,"/class_sz_auxiliary_files/nemo_sims_skyfracs_030722_50bins.txt","sz_selection_function_skyfracs_file")

class_sz_string_parameter(Planck_ylims_file,"/class_sz_auxiliary_files/includes/SZ_ylims.txt","sz_selection_function_ylims_file")
class_sz_string_parameter(SO_ylims_file,"/class_sz_auxiliary_files/nemo_sim_ylims_030722_50bins.txt","sz_selection_function_ylims_file")

class_sz_pclass_sz_parameter(no_spline_in_tinker,int,0)

class_sz_pclass_sz_parameter(sigma_derivative,int,0) // 0 is gradient, and 1 is mcfit
class_sz_pclass_sz_parameter(use_pknl_in_2hterms,int,0)
class_sz_pclass_sz_parameter(use_pknl_in_2hterms_IA_only,int,0)

class_sz_pclass_sz_parameter(use_pkl_in_linbias_calc,int,0)

class_sz_pclass_sz_parameter(use_pk_z_bins,int,0)
class_sz_pclass_sz_parameter(use_pknl_z_bins,int,0)
class_sz_pclass_sz_parameter(pk_z_bins_z1,double,0.)
class_sz_pclass_sz_parameter(pk_z_bins_z2,double,0.)
class_sz_pclass_sz_parameter(pk_z_bins_A0,double,0.)
class_sz_pclass_sz_parameter(pk_z_bins_A1,double,0.)
class_sz_pclass_sz_parameter(pk_z_bins_A2,double,0.)

class_sz_pclass_sz_parameter(n_k_density_profile,int,100)
class_sz_pclass_sz_parameter(n_m_density_profile,int,100)
class_sz_pclass_sz_parameter(n_z_density_profile,int,100)

class_sz_pclass_sz_parameter(fixed_c200m,double,7.)
class_sz_pclass_sz_parameter(fixed_c200c,double,5.)


class_sz_pclass_sz_parameter(n_m_matter_density_profile,int,100)
class_sz_pclass_sz_parameter(n_z_matter_density_profile,int,100)

class_sz_pclass_sz_parameter(szcc_qtrunc,double,1.)
class_sz_pclass_sz_parameter(szcounts_obsscatter,double,1.)

class_sz_pclass_sz_parameter(szcounts_qmin_fft_padded,double,-100.)
class_sz_pclass_sz_parameter(szcounts_qmax_fft_padded,double,100.)
class_sz_pclass_sz_parameter(szcounts_lnqmin_fft,double,-5.)
class_sz_pclass_sz_parameter(szcounts_lnqmax_fft,double,5.)
class_sz_pclass_sz_parameter(tol_dlnm_dlnq,double,0.5)
class_sz_pclass_sz_parameter(ntab_dlnm_dlnq,int,250)

// optimization bias correction for cluster counts
class_sz_pclass_sz_parameter(A_opt_bias_ym,double,0.)
class_sz_pclass_sz_parameter(B_opt_bias_ym,double,0.)

class_sz_pclass_sz_parameter(use_skyaveraged_noise,int,0)
class_sz_pclass_sz_parameter(use_edge_noise_values,int,0)
class_sz_pclass_sz_parameter(szcounts_fft_nqobs,int,50)
class_sz_pclass_sz_parameter(szcounts_fft_nexpected_qobs_n,int,50)
class_sz_pclass_sz_parameter(szcounts_fft_nexpected_qobs_min,double,5.)
class_sz_pclass_sz_parameter(szcounts_fft_nexpected_qobs_max,double,25.)
class_sz_pclass_sz_parameter(szcounts_fft_nz,int,50)
class_sz_pclass_sz_parameter(szcc_dof,double,0.)
class_sz_pclass_sz_parameter(szcounts_fft_nsigmayobs,int,20)
class_sz_pclass_sz_parameter(szcounts_fft_z_min,double,0.)
class_sz_pclass_sz_parameter(szcounts_fft_sigmayobs_min,double,5e-5)
class_sz_pclass_sz_parameter(szcounts_fft_z_max,double,2.)
class_sz_pclass_sz_parameter(szcounts_fft_sigmayobs_max,double,1e-2)

class_sz_pclass_sz_parameter(n_k_pressure_profile,int,50)
class_sz_pclass_sz_parameter(n_k_pressure_profile_2h,int,100)
class_sz_pclass_sz_parameter(n_m_pressure_profile,int,30)
class_sz_pclass_sz_parameter(n_l_pressure_profile,int,300)
class_sz_pclass_sz_parameter(n_z_pressure_profile,int,30)

class_sz_pclass_sz_parameter(n_k_custom1_profile,int,30)
class_sz_pclass_sz_parameter(n_m_custom1_profile,int,60)
class_sz_pclass_sz_parameter(n_z_custom1_profile,int,60)

class_sz_pclass_sz_parameter(array_custom1_redshift_kernel_n_z,int,1000)
class_sz_pclass_sz_parameter(array_b_custom1_n_z,int,1000)

class_sz_pclass_sz_parameter(n_z_psi_b1g,int,50)
class_sz_pclass_sz_parameter(n_l_psi_b1g,int,50)

class_sz_pclass_sz_parameter(n_z_psi_b1kg,int,50)
class_sz_pclass_sz_parameter(n_l_psi_b1kg,int,50)

class_sz_pclass_sz_parameter(n_nu_dcib0dz,int,80)
class_sz_pclass_sz_parameter(n_z_dcib0dz,int,68)

class_sz_pclass_sz_parameter(use_cmb_cls_from_file,int,0)
class_sz_pclass_sz_parameter(skip_cmb,int,0)
// class_sz_pclass_sz_parameter(cosmo_model,int,0)
class_sz_pclass_sz_parameter(skip_sigma8_at_z,int,0)
class_sz_pclass_sz_parameter(skip_sigma8_and_der,int,0)
class_sz_pclass_sz_parameter(skip_chi,int,0)
class_sz_pclass_sz_parameter(skip_hubble,int,0)
class_sz_pclass_sz_parameter(skip_pknl,int,0)
class_sz_pclass_sz_parameter(skip_pkl,int,0)
class_sz_pclass_sz_parameter(skip_pk,int,0)
class_sz_pclass_sz_parameter(want_pp,int,1)
class_sz_pclass_sz_parameter(skip_class_sz,int,0)
class_sz_pclass_sz_parameter(do_real_space_with_mcfit,int,0)
class_sz_pclass_sz_parameter(skip_background_and_thermo,int,0)
class_sz_pclass_sz_parameter(skip_input,int,0)

class_sz_pclass_sz_parameter(n_z_dydz,int,68)

class_sz_pclass_sz_parameter(n_z_psi_b2t,int,50)
class_sz_pclass_sz_parameter(n_l_psi_b2t,int,50)

class_sz_pclass_sz_parameter(n_z_psi_b2g,int,50)
class_sz_pclass_sz_parameter(n_l_psi_b2g,int,50)

class_sz_pclass_sz_parameter(n_z_psi_b2kg,int,50)
class_sz_pclass_sz_parameter(n_l_psi_b2kg,int,50)

class_sz_pclass_sz_parameter(n_z_psi_b1t,int,80)
class_sz_pclass_sz_parameter(n_l_psi_b1t,int,80)

class_sz_pclass_sz_parameter(n_k_n5k,int,50)
class_sz_pclass_sz_parameter(n_l_n5k,int,103)

class_sz_pclass_sz_parameter(k_min_n5k, double, 1.e-4);
class_sz_pclass_sz_parameter(k_max_n5k, double, 1.e2);

class_sz_pclass_sz_parameter(chi_min_n5k_samp_fftw, double, 1.e0);
class_sz_pclass_sz_parameter(chi_max_n5k_samp_fftw, double, 7.e3);

class_sz_pclass_sz_parameter(integrand_n5k_epsrel, double, 1.e-6);
class_sz_pclass_sz_parameter(integrand_n5k_epsabs, double, 7.e-40);

class_sz_pclass_sz_parameter(matter_density_norm_epsrel, double, 1.e-6);
class_sz_pclass_sz_parameter(matter_density_norm_epsabs, double, 1.e-100);


class_sz_pclass_sz_parameter(density_norm_epsrel, double, 1.e-6);
class_sz_pclass_sz_parameter(density_norm_epsabs, double, 1.e-100);

class_sz_pclass_sz_parameter(M1SZ_L_sat, double, 1.e9);
class_sz_pclass_sz_parameter(M2SZ_L_sat, double, 1.e17);
class_sz_pclass_sz_parameter(z1SZ_L_sat, double, 1.e-5);
class_sz_pclass_sz_parameter(z2SZ_L_sat, double, 6.);
class_sz_pclass_sz_parameter( epsabs_L_sat, double, 1e-15);
class_sz_pclass_sz_parameter( epsrel_L_sat, double, 1e-6);

class_sz_pclass_sz_parameter(n_z_L_sat , int,101);
class_sz_pclass_sz_parameter(n_m_L_sat , int,102);
class_sz_pclass_sz_parameter(n_nu_L_sat,int,103);

class_sz_pclass_sz_parameter(use_xout_in_density_profile_from_enclosed_mass , int,0);
class_sz_pclass_sz_parameter(tabulate_rhob_xout_at_m_and_z , int,0);
class_sz_pclass_sz_parameter(n_z_m_to_xout , int,101);
class_sz_pclass_sz_parameter(n_mass_m_to_xout , int,102);

class_sz_pclass_sz_parameter(use_bg_at_z_in_ksz2g_eff, int, 0);
class_sz_pclass_sz_parameter(use_fdndz_for_ksz2g_eff, int, 0);


class_sz_pclass_sz_parameter(n_y_y_to_m,int,50)
class_sz_pclass_sz_parameter(n_z_y_to_m,int,50)

class_sz_pclass_sz_parameter(n_z_psi_b1gt,int,50)
class_sz_pclass_sz_parameter(n_l_psi_b1gt,int,50)

class_sz_pclass_sz_parameter(n_z_psi_b1kgt,int,50)
class_sz_pclass_sz_parameter(n_l_psi_b1kgt,int,50)

class_sz_pclass_sz_parameter(n_z_psi_b1kgg,int,50)
class_sz_pclass_sz_parameter(n_l_psi_b1kgg,int,50)

class_sz_pclass_sz_parameter(N_samp_fftw,int,1024)
class_sz_pclass_sz_parameter(l_min_samp_fftw,double,1e-12)
class_sz_pclass_sz_parameter(l_max_samp_fftw,double,1e9)
class_sz_pclass_sz_parameter(k_min_samp_fftw,double,1e-12)
class_sz_pclass_sz_parameter(k_max_samp_fftw,double,1e9)
class_sz_pclass_sz_parameter(k_min_gas_density_profile,double,1e-3)
class_sz_pclass_sz_parameter(k_max_gas_density_profile,double,1e1)

class_sz_pclass_sz_parameter(k_min_gas_pressure_profile,double,1e-2)
class_sz_pclass_sz_parameter(k_max_gas_pressure_profile,double,1e2)


class_sz_pclass_sz_parameter(l_min_gas_pressure_profile,double,1e-2)
class_sz_pclass_sz_parameter(l_max_gas_pressure_profile,double,5e4)

class_sz_pclass_sz_parameter(fstar_ms,double,0.)

class_sz_pclass_sz_parameter(x_min_matter_density_fftw,double,1e-5)
class_sz_pclass_sz_parameter(x_max_matter_density_fftw,double,1e2)

class_sz_pclass_sz_parameter(x_min_custom1_fftw,double,1e-5)
class_sz_pclass_sz_parameter(x_max_custom1_fftw,double,1e2)

class_sz_pclass_sz_parameter(x_min_gas_density_fftw,double,1e-5)
class_sz_pclass_sz_parameter(x_max_gas_density_fftw,double,1e2)

class_sz_pclass_sz_parameter(x_min_gas_pressure_fftw,double,1e-5)
class_sz_pclass_sz_parameter(x_max_gas_pressure_fftw,double,1e5)


class_sz_pclass_sz_parameter(x_out_truncated_gas_density_profile_normalization,double,1.)

class_sz_pclass_sz_parameter(x_out_matter_density_profile_normalization,double,2.)
class_sz_pclass_sz_parameter(x_out_matter_density_profile,double,2.)

class_sz_pclass_sz_parameter(has_pk,int,0)
class_sz_pclass_sz_parameter(ngal_ngal_auto_only,int,0)


// class_sz_pclass_sz_parameter(has_b_custom1,int,0)

class_sz_pclass_sz_parameter(use_fft_for_profiles_transform,int,0)

class_sz_pclass_sz_parameter(has_tracer_bias_zdependence,int,0)

class_sz_pclass_sz_parameter(k_min_gas_pressure_profile_2h,double,1e-2)
class_sz_pclass_sz_parameter(k_max_gas_pressure_profile_2h,double,1e2)
class_sz_pclass_sz_parameter(normalize_gas_density_profile,int,0)

class_sz_pclass_sz_parameter(ell_min_kSZ2_gal_multipole_grid,double,2)
class_sz_pclass_sz_parameter(ell_max_kSZ2_gal_multipole_grid,double,1e5)


class_sz_pclass_sz_parameter(n_z_hmf_counter_terms,int,200)
class_sz_pclass_sz_parameter(array_profile_ln_PgNFW_at_lnl_over_ls_size,int,200)
// class_sz_pclass_sz_parameter(m_min_counter_terms,double,1e11) // this is then set to M1SZ ("m_min") in class_sz.c
// class_sz_pclass_sz_parameter(m_max_counter_terms,double,1e17)

class_sz_pclass_sz_parameter(hmf_apply_zthreshold_to_hmf_and_bias,int,1) // see sec 4 of https://arxiv.org/pdf/0803.2706.pdf -- CCL doesnt have this


class_sz_pclass_sz_parameter(mass_epsrel_ngbar,double,1e-6)
class_sz_pclass_sz_parameter(mass_epsabs_ngbar,double,1e-40)

class_sz_pclass_sz_parameter(m_to_xout_epsrel,double,1e-6)
class_sz_pclass_sz_parameter(m_to_xout_epsabs,double,1e-40)


// class_sz_pclass_sz_parameter(n_m_dndlnM,int,200)

//printf("-> File Name pr: %s\n",pclass_sz->WISC3_dndz_file);



#undef class_sz_pclass_sz_parameter
#undef class_sz_string_parameter
#undef class_sz_type_parameter
