#include "class_sz_macros.h"
class_sz_string_parameter(root,"root","root")
class_sz_string_parameter(append_name_trispectrum_ref,"append_name_trispectrum_ref","append_name_trispectrum_ref")
class_sz_string_parameter(path_to_ref_trispectrum_for_cobaya,"path_to_ref_trispectrum_for_cobaya","path_to_ref_trispectrum_for_cobaya")
class_sz_string_parameter(full_path_to_noise_curve_for_y_y,"full_path_to_noise_curve_for_y_y","full_path_to_noise_curve_for_y_y")
class_sz_string_parameter(full_path_to_noise_curve_for_t_t,"full_path_to_noise_curve_for_t_t","full_path_to_noise_curve_for_t_t")
class_sz_string_parameter(full_path_to_dndz_gal,"/sz_auxiliary_files/WISC_bin3.txt","custom dndz file")
class_sz_string_parameter(full_path_to_source_dndz_gal,"/sz_auxiliary_files/WISC_bin3.txt","custom dndz file")

//class_sz_string_parameter(UNWISE_dndz_file,"/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz.txt","unWISE dndz file")
class_sz_string_parameter(UNWISE_dndz_file,"/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz_cosmos.txt","unWISE dndz file")
class_sz_string_parameter(UNWISE_fdndz_file,"/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_fdndz.txt","unWISE fdndz file")
class_sz_string_parameter(UNWISE_cosmos_dndz_file,"/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz_cosmos.txt","unWISE fdndz file")
class_sz_string_parameter(WISC3_dndz_file,"/sz_auxiliary_files/WISC_bin3.txt","WISC bin 3 file")
class_sz_string_parameter(A10_file,"/sz_auxiliary_files/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_A10.txt","A10 file")
class_sz_string_parameter(P13_file,"/sz_auxiliary_files/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_P13.txt","P13 file")
class_sz_string_parameter(Tinker_et_al_10_alpha_consistency_msyriac_file,"/sz_auxiliary_files/Tinker_et_al_10_alpha_consistency_msyriac.txt","Tinker_et_al_10_alpha_consistency_msyriac_file")

class_sz_string_parameter(ksz_filter_file,"/sz_auxiliary_files/UNWISE_galaxy_distributions/unwise_filter_functions_l_fl.txt","ksz filter file")
class_sz_string_parameter(ksz_template_file,"/sz_auxiliary_files/cl_ksz_bat.dat","ksz template file")
class_sz_string_parameter(ksz_reio_template_file,"/sz_auxiliary_files/FBN_kSZ_PS_patchy.d.txt","ksz template file, reio contribution")


class_sz_string_parameter(Planck_thetas_file,"/sz_auxiliary_files/SZ_thetas.txt","Planck_thetas_file")
class_sz_string_parameter(SO_thetas_file,"/sz_auxiliary_files/so_3freqs_191121_thetas.txt","SO_thetas_file")
class_sz_string_parameter(Planck_skyfracs_file,"/sz_auxiliary_files/SZ_skyfracs.txt","Planck_skyfracs_file")
class_sz_string_parameter(SO_skyfracs_file,"/sz_auxiliary_files/so_3freqs_191121_skyfracs.txt","SO_skyfracs_file")
class_sz_string_parameter(Planck_ylims_file,"/sz_auxiliary_files/SZ_ylims.txt","Planck_ylims_file")
class_sz_string_parameter(SO_ylims_file,"/sz_auxiliary_files/so_3freqs_191121_ylims.txt","SO_ylims_file")


class_sz_ptsz_parameter(n_ell_density_profile,int,100)
class_sz_ptsz_parameter(n_m_density_profile,int,100)
class_sz_ptsz_parameter(n_z_density_profile,int,100)

class_sz_ptsz_parameter(n_ell_pressure_profile,int,50)
class_sz_ptsz_parameter(n_m_pressure_profile,int,300)
class_sz_ptsz_parameter(n_z_pressure_profile,int,300)

class_sz_ptsz_parameter(n_z_psi_b1g,int,50)
class_sz_ptsz_parameter(n_l_psi_b1g,int,50)

class_sz_ptsz_parameter(n_z_psi_b1kg,int,50)
class_sz_ptsz_parameter(n_l_psi_b1kg,int,50)

class_sz_ptsz_parameter(n_nu_dcib0dz,int,80)
class_sz_ptsz_parameter(n_z_dcib0dz,int,68)

class_sz_ptsz_parameter(n_z_dydz,int,68)

class_sz_ptsz_parameter(n_z_psi_b2t,int,50)
class_sz_ptsz_parameter(n_l_psi_b2t,int,50)

class_sz_ptsz_parameter(n_z_psi_b2g,int,50)
class_sz_ptsz_parameter(n_l_psi_b2g,int,50)

class_sz_ptsz_parameter(n_z_psi_b2kg,int,50)
class_sz_ptsz_parameter(n_l_psi_b2kg,int,50)

class_sz_ptsz_parameter(n_z_psi_b1t,int,80)
class_sz_ptsz_parameter(n_l_psi_b1t,int,80)


class_sz_ptsz_parameter(M1SZ_L_sat, double, 1.e9);
class_sz_ptsz_parameter(M2SZ_L_sat, double, 1.e17);
class_sz_ptsz_parameter(z1SZ_L_sat, double, 1.e-5);
class_sz_ptsz_parameter(z2SZ_L_sat, double, 6.);
class_sz_ptsz_parameter( epsabs_L_sat, double, 1e-15);
class_sz_ptsz_parameter( epsrel_L_sat, double, 1e-6);
class_sz_ptsz_parameter(n_z_L_sat , int,101);
class_sz_ptsz_parameter(n_m_L_sat , int,102);
class_sz_ptsz_parameter(n_nu_L_sat, int, 103);

class_sz_ptsz_parameter(use_bg_at_z_in_ksz2g_eff, int, 0);
class_sz_ptsz_parameter(use_fdndz_for_ksz2g_eff, int, 0);

class_sz_ptsz_parameter(n_z_psi_b1gt,int,50)
class_sz_ptsz_parameter(n_l_psi_b1gt,int,50)

class_sz_ptsz_parameter(n_z_psi_b1kgt,int,50)
class_sz_ptsz_parameter(n_l_psi_b1kgt,int,50)

class_sz_ptsz_parameter(N_samp_fftw,int,150)
class_sz_ptsz_parameter(l_min_samp_fftw,double,1e-12)
class_sz_ptsz_parameter(l_max_samp_fftw,double,1e9)
class_sz_ptsz_parameter(k_min_gas_density_profile,double,1e-4)
class_sz_ptsz_parameter(normalize_gas_density_profile,int,0)

class_sz_ptsz_parameter(ell_min_kSZ2_gal_multipole_grid,double,2)
class_sz_ptsz_parameter(ell_max_kSZ2_gal_multipole_grid,double,1e5)


class_sz_ptsz_parameter(n_z_hmf_counter_terms,int,200)
class_sz_ptsz_parameter(array_profile_ln_PgNFW_at_lnl_over_ls_size,int,200)
class_sz_ptsz_parameter(m_min_counter_terms,double,1e11)
class_sz_ptsz_parameter(m_max_counter_terms,double,1e17)


class_sz_ptsz_parameter(mass_epsrel_ngbar,double,1e-6)
class_sz_ptsz_parameter(mass_epsabs_ngbar,double,1e-40)


// class_sz_ptsz_parameter(n_m_dndlnM,int,200)

//printf("-> File Name pr: %s\n",ptsz->WISC3_dndz_file);



#undef class_sz_ptsz_parameter
#undef class_sz_string_parameter
#undef class_sz_type_parameter
