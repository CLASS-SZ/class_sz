#include "class_sz_macros.h"
class_sz_string_parameter(root,"root","root")
class_sz_string_parameter(append_name_trispectrum_ref,"append_name_trispectrum_ref","append_name_trispectrum_ref")
class_sz_string_parameter(path_to_ref_trispectrum_for_cobaya,"path_to_ref_trispectrum_for_cobaya","path_to_ref_trispectrum_for_cobaya")
class_sz_string_parameter(full_path_to_noise_curve_for_y_y,"full_path_to_noise_curve_for_y_y","full_path_to_noise_curve_for_y_y")
class_sz_string_parameter(full_path_to_dndz_gal,"/sz_auxiliary_files/WISC_bin3.txt","custom dndz file")
//class_sz_string_parameter(UNWISE_dndz_file,"/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz.txt","unWISE dndz file")
class_sz_string_parameter(UNWISE_dndz_file,"/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz_cosmos.txt","unWISE dndz file")
class_sz_string_parameter(UNWISE_fdndz_file,"/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_fdndz.txt","unWISE fdndz file")
class_sz_string_parameter(UNWISE_cosmos_dndz_file,"/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz_cosmos.txt","unWISE fdndz file")
class_sz_string_parameter(WISC3_dndz_file,"/sz_auxiliary_files/WISC_bin3.txt","WISC bin 3 file")
class_sz_string_parameter(A10_file,"/sz_auxiliary_files/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_A10.txt","A10 file")
class_sz_string_parameter(P13_file,"/sz_auxiliary_files/class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_P13.txt","P13 file")
class_sz_string_parameter(Tinker_et_al_10_alpha_consistency_msyriac_file,"/sz_auxiliary_files/Tinker_et_al_10_alpha_consistency_msyriac.txt","Tinker_et_al_10_alpha_consistency_msyriac_file")

class_sz_string_parameter(ksz_filter_file,"/sz_auxiliary_files/UNWISE_galaxy_distributions/unwise_filter_functions_l_fl.txt","ksz filter file")


class_sz_string_parameter(Planck_thetas_file,"/sz_auxiliary_files/SZ_thetas.txt","Planck_thetas_file")
class_sz_string_parameter(SO_thetas_file,"/sz_auxiliary_files/so_3freqs_120921_pass2_thetas.txt","SO_thetas_file")
class_sz_string_parameter(Planck_skyfracs_file,"/sz_auxiliary_files/SZ_skyfracs.txt","Planck_skyfracs_file")
class_sz_string_parameter(SO_skyfracs_file,"/sz_auxiliary_files/so_3freqs_120921_pass2_skyfracs.txt","SO_skyfracs_file")
class_sz_string_parameter(Planck_ylims_file,"/sz_auxiliary_files/SZ_ylims.txt","Planck_ylims_file")
class_sz_string_parameter(SO_ylims_file,"/sz_auxiliary_files/so_3freqs_120921_pass2_ylims.txt","SO_ylims_file")


class_sz_ptsz_parameter(n_ell_density_profile,int,100)
class_sz_ptsz_parameter(n_m_density_profile,int,100)
class_sz_ptsz_parameter(n_z_density_profile,int,100)

class_sz_ptsz_parameter(n_ell_pressure_profile,int,70)
class_sz_ptsz_parameter(n_m_pressure_profile,int,70)
class_sz_ptsz_parameter(n_z_pressure_profile,int,70)

class_sz_ptsz_parameter(n_z_psi_b1g,int,50)
class_sz_ptsz_parameter(n_l_psi_b1g,int,50)

class_sz_ptsz_parameter(n_nu_dcib0dz,int,80)
class_sz_ptsz_parameter(n_z_dcib0dz,int,68)

class_sz_ptsz_parameter(n_z_dydz,int,68)

class_sz_ptsz_parameter(n_z_psi_b2t,int,50)
class_sz_ptsz_parameter(n_l_psi_b2t,int,50)

class_sz_ptsz_parameter(n_z_psi_b2g,int,50)
class_sz_ptsz_parameter(n_l_psi_b2g,int,50)

class_sz_ptsz_parameter(n_z_psi_b1t,int,80)
class_sz_ptsz_parameter(n_l_psi_b1t,int,80)

class_sz_ptsz_parameter(n_z_psi_b1gt,int,50)
class_sz_ptsz_parameter(n_l_psi_b1gt,int,50)

class_sz_ptsz_parameter(ell_min_kSZ2_gal_multipole_grid,double,2)
class_sz_ptsz_parameter(ell_max_kSZ2_gal_multipole_grid,double,1e4)


class_sz_ptsz_parameter(n_z_hmf_counter_terms,int,200)
class_sz_ptsz_parameter(array_profile_ln_PgNFW_at_lnl_over_ls_size,int,200)
class_sz_ptsz_parameter(m_min_counter_terms,double,1e11)
class_sz_ptsz_parameter(m_max_counter_terms,double,1e17)

// class_sz_ptsz_parameter(n_m_dndlnM,int,200)

//printf("-> File Name pr: %s\n",ptsz->WISC3_dndz_file);



#undef class_sz_ptsz_parameter
#undef class_sz_string_parameter
#undef class_sz_type_parameter
