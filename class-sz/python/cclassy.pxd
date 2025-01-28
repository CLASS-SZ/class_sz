# Bunch of declarations from C to python. The idea here is to define only the
# quantities that will be used, for input, output or intermediate manipulation,
# by the python wrapper. For instance, in the precision structure, the only
# item used here is its error message. That is why nothing more is defined from
# this structure. The rest is internal in Class.
# If, for whatever reason, you need an other, existing parameter from Class,
# remember to add it inside this cdef.

DEF _MAX_NUMBER_OF_K_FILES_ = 30
DEF _MAXTITLESTRINGLENGTH_ = 8000
DEF _FILENAMESIZE_ = 256
DEF _LINE_LENGTH_MAX_ = 8192

cdef extern from "class.h":

    cdef char[10] _VERSION_

    ctypedef char FileArg[40]

    ctypedef char* ErrorMsg

    ctypedef char FileName[_FILENAMESIZE_]

    cdef enum linear_or_logarithmic:
        linear
        logarithmic

    cdef enum file_format:
        class_format
        camb_format

    cdef enum non_linear_method:
        nl_none
        nl_halofit
        nl_HMcode

    cdef enum pk_outputs:
        pk_linear
        pk_nonlinear

    cdef enum out_sigmas:
        out_sigma
        out_sigma_prime
        out_sigma_disp

    cdef struct precision:
        ErrorMsg error_message

    cdef struct background:
        ErrorMsg error_message
        int bg_size
        int index_bg_ang_distance
        int index_bg_lum_distance
        int index_bg_conf_distance
        int index_bg_a
        int index_bg_H
        int index_bg_D
        int index_bg_f
        int index_bg_Omega_m
        short long_info
        short inter_normal
        short  has_ncdm
        double T_cmb
        double h
        double H0
        double age
        double conformal_age
        double * m_ncdm_in_eV
        double Neff
        double Omega0_g
        double Omega0_b
        double Omega0_idr
        double T_idr
        double Omega0_idm_dr
        double Omega0_cdm
        double Omega0_dcdm
        double Omega0_ncdm_tot
        double Omega0_lambda
        double Omega0_fld
        double w0_fld
        double wa_fld
        double cs2_fld
        double Omega0_ur
        double Omega0_dcdmdr
        double Omega0_dr
        double Omega0_scf
        double Omega0_k
        int bt_size
        double Omega0_m
        double Omega0_r
        double Omega0_de
        double a_eq
        double H_eq
        double z_eq
        double tau_eq

    cdef struct thermo:
        ErrorMsg error_message
        int th_size
        int index_th_xe
        int index_th_Tb
        short inter_normal
        double tau_reio
        double z_reio
        double z_rec
        double tau_rec
        double rs_rec
        double ds_rec
        double da_rec
        double z_star
        double tau_star
        double rs_star
        double ds_star
        double ra_star
        double da_star
        double rd_star
        double z_d
        double tau_d
        double ds_d
        double rs_d
        double YHe
        double n_e
        double a_idm_dr
        double b_idr
        double nindex_idm_dr
        double m_idm

        int tt_size

    cdef struct perturbs:
        ErrorMsg error_message
        short has_scalars
        short has_vectors
        short has_tensors

        short has_density_transfers
        short has_velocity_transfers

        int has_pk_matter
        int l_lss_max

        int store_perturbations
        int k_output_values_num
        double k_output_values[_MAX_NUMBER_OF_K_FILES_]
        double k_max_for_pk
        int index_k_output_values[_MAX_NUMBER_OF_K_FILES_]
        char scalar_titles[_MAXTITLESTRINGLENGTH_]
        char vector_titles[_MAXTITLESTRINGLENGTH_]
        char tensor_titles[_MAXTITLESTRINGLENGTH_]
        int number_of_scalar_titles
        int number_of_vector_titles
        int number_of_tensor_titles


        double * scalar_perturbations_data[_MAX_NUMBER_OF_K_FILES_]
        double * vector_perturbations_data[_MAX_NUMBER_OF_K_FILES_]
        double * tensor_perturbations_data[_MAX_NUMBER_OF_K_FILES_]
        int size_scalar_perturbation_data[_MAX_NUMBER_OF_K_FILES_]
        int size_vector_perturbation_data[_MAX_NUMBER_OF_K_FILES_]
        int size_tensor_perturbation_data[_MAX_NUMBER_OF_K_FILES_]

        double * alpha_idm_dr
        double * beta_idr

        int * k_size
        int * ic_size
        int index_md_scalars

    cdef struct transfers:
        ErrorMsg error_message

    cdef struct primordial:
        ErrorMsg error_message
        double k_pivot
        double A_s
        double n_s
        double alpha_s
        double beta_s
        double r
        double n_t
        double alpha_t
        double V0
        double V1
        double V2
        double V3
        double V4
        double f_cdi
        double n_cdi
        double c_ad_cdi
        double n_ad_cdi
        double f_nid
        double n_nid
        double c_ad_nid
        double n_ad_nid
        double f_niv
        double n_niv
        double c_ad_niv
        double n_ad_niv
        double phi_min
        double phi_max
        int lnk_size

    cdef struct spectra:
        ErrorMsg error_message
        int has_tt
        int has_te
        int has_ee
        int has_bb
        int has_pp
        int has_tp
        int has_dd
        int has_td
        int has_ll
        int has_dl
        int has_tl
        int l_max_tot
        int ** l_max_ct
        int ct_size
        int * ic_size
        int * ic_ic_size
        int md_size
        int d_size
        int non_diag
        int index_ct_tt
        int index_ct_te
        int index_ct_ee
        int index_ct_bb
        int index_ct_pp
        int index_ct_tp
        int index_ct_dd
        int index_ct_td
        int index_ct_pd
        int index_ct_ll
        int index_ct_dl
        int index_ct_tl
        int * l_size
        int index_md_scalars

    cdef struct output:
        ErrorMsg error_message

    cdef struct class_sz_structure:
        ErrorMsg error_message
        double A_cib
        double A_cn
        double A_ir
        double A_rs
        double szunbinned_loglike
        double Sigma8OmegaM_SZ
        int nlSZ
        double chi_star
        double Rho_crit_0
        double z1SZ
        double z2SZ
        int ndim_redshifts
        int ndim_redshifts_for_integral
        double M1SZ
        double M2SZ
        int do_real_space_with_mcfit
        int jax
        int skip_class_sz
        int skip_pk
        int skip_cmb
        int cosmo_model
        int use_cmb_cls_from_file
        int skip_chi
        int skip_hubble
        int want_pp
        int skip_pkl
        int skip_pknl
        int skip_sigma8_at_z
        int skip_sigma8_and_der
        int skip_background_and_thermo
        int delta_def_matter_density
        double cl_gal_gal_A_sn 
        double * array_m_dndlnM
        int n_m_dndlnM
        double alpha_cib
        double T0_cib
        double beta_cib
        double gamma_cib
        double delta_cib
        double L0_cib
        double sigma2_LM_cib
        double m_eff_cib
        double z_for_pk_hm
        double P0GNFW
        double c500
        double gammaGNFW
        double alphaGNFW
        double betaGNFW
        double alpha_s_HOD
        double sigma_log10M_HOD
        double M1_prime_HOD
        double M_min_HOD
        double * alpha_s_HOD_ngal
        double * sigma_log10M_HOD_ngal
        double * M1_prime_HOD_ngal
        double * M_min_HOD_ngal
        double * M0_HOD_ngal
        double * dndz_stretch_ngal
        double * dndz_shift_ngal
        double * f_cen_HOD_ngal
        double * x_out_truncated_nfw_profile_satellite_galaxies_ngal
        double * photo_z_params_ngal
        double * centrals_only_ngal
        double dndz_shift_source_gal
        double dndz_stretch_source_gal
        double dndz_shift_gal
        double dndz_stretch_gal
        double shear_calibration_m
        double A_IA
        double eta_IA
        double x_inSZ
        double x_outSZ
        double HSEbias
        double mu_e
        double Omega_m_0
        double Omega_r_0
        double Omega0_cdm
        double Omega0_b
        double f_free
        int delta_def_custom1
        double x_out_custom1
        int  ndim_masses
        int nbins_M
        int n_k_for_pk_hm
        double logR1SZ
        double logR2SZ
        double delta_cSZ
        int szcounts_ntot
        double szcounts_fft_z_min
        double szcounts_fft_z_max
        double szcounts_qmax_fft_padded
        int szcounts_fft_nz
        double sn_cutoff
        double * szrate
        int szcat_size
        double * szcat_snr
        double * szcat_z
        int need_sigma
        int has_vrms2
        int has_pk
        double effective_galaxy_bias
        double alphaSZ
        double beta0SZ
        double gamma0SZ
        double phi0SZ
        double eta0SZ
        double alpha_break_pressure
        double M_break_pressure
        double sigma8_Pcb
        double * frequencies_for_cib
        double * cib_monopole
        double * cib_shotnoise
        double y_monopole
        int n_frequencies_for_cib
        double * ell
        double * cl_sz_1h
        double * cl_sz_2h
        double * cl_isw_tsz
        double * cl_te_y_y
        double * cl_tSZ_gallens_1h
        double * cl_tSZ_gallens_2h
        double * cl_gallens_lensmag_1h
        double * cl_gallens_lensmag_2h
        double * cl_tSZ_gal_1h
        double * cl_tSZ_gal_2h
        double * cl_IA_gal_2h
        double * cl_kSZ_kSZ_gal_1h
        double * cl_kSZ_kSZ_gal_2h
        double * cl_kSZ_kSZ_gal_3h
        double * cl_kSZ_kSZ_gal_1h_fft
        double * cl_kSZ_kSZ_gal_2h_fft
        double * cl_kSZ_kSZ_gal_3h_fft
        double * cl_kSZ_kSZ_gal_hf
        double * cl_kSZ_kSZ_gal_lensing_term
        double * cov_ll_kSZ_kSZ_gal
        double * cl_t2t2f
        double * cl_kSZ_kSZ_gallens_1h_fft
        double * cl_kSZ_kSZ_gallens_2h_fft
        double * cl_kSZ_kSZ_gallens_3h_fft
        double * cl_kSZ_kSZ_gallens_hf
        double * cl_kSZ_kSZ_gallens_lensing_term
        double * cov_ll_kSZ_kSZ_gallens
        double * cl_kSZ_kSZ_lens_1h_fft
        double * cl_kSZ_kSZ_lens_2h_fft
        double * cl_kSZ_kSZ_lens_3h_fft
        double * cl_gal_gal_lens_1h_fft
        double * cl_gal_gal_lens_2h_fft
        double * cl_gal_gal_lens_3h_fft
        double * cl_kSZ_kSZ_lens_hf
        double * cl_kSZ_kSZ_lens_lensing_term
        double * cov_ll_kSZ_kSZ_lens
        double * cl_tSZ_lensmag_1h
        double * cl_tSZ_lensmag_2h
        double * cl_tSZ_lens_1h
        double * cl_tSZ_lens_2h
        double * cl_gal_gallens_1h
        double * cl_gal_gallens_2h
        double * cl_gallens_gallens_1h
        double * cl_gallens_gallens_2h
        double * cl_gallens_lens_1h
        double * cl_gallens_lens_2h
        double * thetas_arcmin
        double * gamma_gal_gallens_1h
        double * gamma_gal_gallens_2h
        double * cl_kSZ_kSZ_1h
        double * cl_kSZ_kSZ_2h
        double * cl_gal_gal_1h
        double * cl_gal_gal_2h
        double * cl_gal_gal_hf
        double * cl_gal_lens_hf
        double * cl_tau_gal_1h
        double * cl_tau_gal_2h
        double * cl_tau_tau_1h
        double * cl_tau_tau_2h
        double * cl_gal_lens_1h
        double * cl_gal_lens_2h
        double * cl_lens_lens_1h
        double * cl_custom1_custom1_1h
        double * cl_custom1_custom1_2h
        double * cl_custom1_lens_1h
        double * cl_custom1_lens_2h
        double * cl_custom1_tSZ_1h
        double * cl_custom1_tSZ_2h
        double * cl_custom1_gal_1h
        double * cl_custom1_gal_2h
        double * cl_custom1_gallens_1h
        double * cl_custom1_gallens_2h
        double ** cl_custom1_cib_1h
        double ** cl_custom1_cib_2h
        double * cl_lens_lens_2h
        double * cl_lens_lens_hf
        double *** cl_cib_cib_1h
        double *** cl_ngal_ngal_1h
        double *** cl_ngal_ngal_2h
        double *** cl_ngal_ngal_hf
        double ** cl_ngal_lens_1h
        double ** cl_ngal_lens_2h
        double ** cl_ngal_lens_hf
        double ** cl_ngal_gallens_1h
        double ** cl_ngal_gallens_2h
        double ** cl_ngal_tsz_1h
        double ** cl_ngal_tsz_2h
        double ** cl_ngal_IA_2h
        double ** cl_nlensmag_gallens_1h
        double ** cl_nlensmag_gallens_2h
        double ** cl_nlensmag_tsz_1h
        double ** cl_nlensmag_tsz_2h
        double *** cl_cib_cib_2h
        double **  cl_gallens_cib_1h
        double **  cl_gallens_cib_2h
        double **  cl_tSZ_cib_1h
        double **  cl_tSZ_cib_2h
        double **  cl_gal_cib_1h
        double **  cl_gal_cib_2h
        double **  cl_lens_cib_1h
        double **  cl_lens_cib_2h
        double * cib_frequency_list
        int cib_frequency_list_num
        int * galaxy_samples_list
        int galaxy_samples_list_num
        double * pk_at_z_1h
        double * pk_at_z_2h
        double * bk_at_z_1h
        double * bk_at_z_2h
        double * bk_at_z_3h
        double * bk_ttg_at_z_1h
        double * bk_ttg_at_z_2h
        double * bk_ttg_at_z_3h
        double * b_tSZ_tSZ_tSZ_1halo
        double * b_tSZ_tSZ_tSZ_2h
        double * b_tSZ_tSZ_tSZ_3h
        double * b_kSZ_kSZ_tSZ_1h
        double * b_kSZ_kSZ_tSZ_2h
        double * b_kSZ_kSZ_tSZ_3h
        double * k_for_pk_hm
        double * pk_gg_at_z_1h
        double * pk_gg_at_z_2h
        double * pk_bb_at_z_1h
        double * pk_bb_at_z_2h
        double * pk_b_at_z_2h
        double * pk_em_at_z_1h
        double * pk_em_at_z_2h
        double * cl_gal_lensmag_1h
        double * cl_gal_lensmag_2h
        double * cl_gal_lensmag_hf
        double * cl_lens_lensmag_1h
        double * cl_lens_lensmag_2h
        double * cl_lens_lensmag_hf
        double * cl_lensmag_lensmag_1h
        double * cl_lensmag_lensmag_2h
        double * cl_lensmag_lensmag_hf
        double ** tllprime_sz
        double * cov_Y_N_mass_bin_edges
        double * cov_N_N
        double ** cov_N_N_hsv
        double * array_sigma_at_z_and_R
        double * array_dsigma2dR_at_z_and_R
        double * array_radius
        double * array_redshift
        double * array_lnk
        double * array_pknl_at_z_and_k
        double * array_pkl_at_z_and_k
        double * array_vrms2_at_z
        int array_b_custom1_n_z
        double * array_b_custom1_ln1pz
        double * array_b_custom1_bias
        double * array_custom1_redshift_kernel_W
        double * array_custom1_redshift_kernel_ln1pz
        int array_custom1_redshift_kernel_n_z
        int n_z_custom1_profile
        double * array_custom1_profile_ln_1pz
        int n_m_custom1_profile
        double * array_custom1_profile_ln_m
        int n_k_custom1_profile
        double * array_custom1_profile_ln_x
        double ** array_custom1_profile_u_at_lnx_lnm_ln1pz
        int has_custom1
        int has_b_custom1
        short has_class_sz_structure
        short sz_verbose
        double bin_dlog10_snr_last_bin
        double fNL
        int use_class_sz_fast_mode
        double P0_B12
        double beta_B12
        double alpha_B12
        double gamma_B12
        double xc_B12
        double alpha_m_P0_B12
        double alpha_m_xc_B12
        double alpha_m_beta_B12
        double alpha_z_P0_B12
        double alpha_z_xc_B12
        double alpha_z_beta_B12
        double mcut_B12
        double c_B12
        double alphap_m_P0_B12
        double alphap_m_xc_B12
        double alphap_m_beta_B12
        double alpha_c_P0_B12
        double alpha_c_xc_B12
        double alpha_c_beta_B12
        double A_rho0
        double A_alpha
        double A_beta
        double alpha_m_rho0
        double alpha_m_alpha
        double alpha_m_beta
        double alpha_z_rho0
        double alpha_z_alpha
        double alpha_z_beta
        double mcut
        double c_B16
        double alphap_m_rho0
        double alphap_m_alpha
        double alphap_m_beta
        double alpha_c_rho0
        double alpha_c_alpha
        double alpha_c_beta
        double xc_B16
        double gamma_B16
        


    cdef struct szcount:
        double ** dNdzdy_theoretical
        double ystar
        double alpha
        double sigmaM
        int has_completeness
        int nzSZ
        int size_logM
        int Nbins_z
        int Nbins_y
        double rho_m_at_z
        double dlogy
        double dz
        double * z_center
        double * logy


    cdef struct lensing:
        int has_tt
        int has_ee
        int has_te
        int has_bb
        int has_pp
        int has_tp
        int has_dd
        int has_td
        int has_ll
        int has_dl
        int has_tl
        int index_lt_tt
        int index_lt_te
        int index_lt_ee
        int index_lt_bb
        int index_lt_pp
        int index_lt_tp
        int index_lt_dd
        int index_lt_td
        int index_lt_ll
        int index_lt_dl
        int index_lt_tl
        int * l_max_lt
        int lt_size
        int has_lensed_cls
        int l_lensed_max
        int l_unlensed_max
        ErrorMsg error_message

    cdef struct nonlinear:
        short has_pk_matter
        int method
        int ic_size
        int ic_ic_size
        int k_size
        int ln_tau_size
        int tau_size
        int index_tau_min_nl
        double * k
        double * ln_tau
        double * tau
        double ** ln_pk_l
        double ** ln_pk_nl
        double * sigma8
        int has_pk_m
        int has_pk_cb
        int index_pk_m
        int index_pk_cb
        int index_pk_total
        int index_pk_cluster
        ErrorMsg error_message

    cdef struct file_content:
        char * filename
        int size
        FileArg * name
        FileArg * value
        short * read

    void lensing_free(void*)
    void spectra_free(void*)
    void transfer_free(void*)
    void primordial_free(void*)
    void perturb_free(void*)
    void thermodynamics_free(void*)
    void background_free(void*)
    void nonlinear_free(void*)
    void class_sz_free(void*)
    void szcounts_free(void*,void*)

    cdef int _FAILURE_
    cdef int _FALSE_
    cdef int _TRUE_

    int input_init(void*, void*, void*, void*, void*, void*, void*, void*, void*,
        void*, void*, void*, void*, char*)
    int background_init(void*,void*)
    int thermodynamics_init(void*,void*,void*)
    int perturb_init(void*,void*,void*,void*)
    int primordial_init(void*,void*,void*)
    int nonlinear_init(void*,void*,void*,void*,void*,void*)
    int transfer_init(void*,void*,void*,void*,void*,void*)
    int spectra_init(void*,void*,void*,void*,void*,void*,void*,void*,void*,void*)
    int lensing_init(void*,void*,void*,void*,void*)
    int class_sz_integrate_init(void*,void*,void*,void*,void*,void*,void*,void*,void*)
    int class_sz_tabulate_init(void*,void*,void*,void*,void*,void*,void*,void*,void*)
    int class_sz_cosmo_init(void*,void*,void*,void*,void*,void*,void*,void*,void*)
    int szcount_init(void*,void*,void*,void*,void*)

    int background_tau_of_z(void* pba, double z,double* tau)
    int background_at_tau(void* pba, double tau, short return_format, short inter_mode, int * last_index, double *pvecback)
    int background_output_titles(void * pba, char titles[_MAXTITLESTRINGLENGTH_])
    int background_output_data(void *pba, int number_of_titles, double *data)

    int thermodynamics_at_z(void * pba, void * pth, double z, short inter_mode, int * last_index, double *pvecback, double *pvecthermo)
    int thermodynamics_output_titles(void * pba, void *pth, char titles[_MAXTITLESTRINGLENGTH_])
    int thermodynamics_output_data(void *pba, void *pth, int number_of_titles, double *data)

    int perturb_output_data(void *pba,void *ppt, file_format output_format, double z, int number_of_titles, double *data)
    int perturb_output_firstline_and_ic_suffix(void *ppt, int index_ic, char first_line[_LINE_LENGTH_MAX_], FileName ic_suffix)
    int perturb_output_titles(void *pba, void *ppt,  file_format output_format, char titles[_MAXTITLESTRINGLENGTH_])

    int primordial_output_titles(void * ppt, void *ppm, char titles[_MAXTITLESTRINGLENGTH_])
    int primordial_output_data(void *ppt, void *ppm, int number_of_titles, double *data)

    int spectra_cl_at_l(void* psp,double l,double * cl,double * * cl_md,double * * cl_md_ic)
    int lensing_cl_at_l(void * ple,int l,double * cl_lensed)

    int spectra_pk_at_z(
        void * pba,
        void * psp,
        int mode,
        double z,
        double * output_tot,
        double * output_ic,
        double * output_cb_tot,
        double * output_cb_ic
        )

    int spectra_pk_at_k_and_z(
        void* pba,
        void * ppm,
        void * psp,
        double k,
        double z,
        double * pk,
        double * pk_ic,
        double * pk_cb,
        double * pk_cb_ic)

    int spectra_pk_nl_at_k_and_z(
        void* pba,
        void * ppm,
        void * psp,
        double k,
        double z,
        double * pk,
        double * pk_cb)

    int spectra_pk_nl_at_z(
        void * pba,
        void * psp,
        int mode,
        double z,
        double * output_tot,
        double * output_cb_tot)

    int nonlinear_pk_at_k_and_z(
        void * pba,
        void * ppm,
        void * pnl,
        int pk_output,
        double k,
        double z,
        int index_pk,
        double * out_pk,
        double * out_pk_ic)

    int nonlinear_pk_tilt_at_k_and_z(
        void * pba,
        void * ppm,
        void * pnl,
        int pk_output,
        double k,
        double z,
        int index_pk,
        double * pk_tilt)

    int nonlinear_sigmas_at_z(
        void * ppr,
        void * pba,
        void * pnl,
        double R,
        double z,
        int index_pk,
        int sigma_output,
        double * result)

    int nonlinear_pks_at_kvec_and_zvec(
        void * pba,
        void * pnl,
        int pk_output,
        double * kvec,
        int kvec_size,
        double * zvec,
        int zvec_size,
        double * out_pk,
        double * out_pk_cb)

    int nonlinear_hmcode_sigma8_at_z(void* pba, void* pnl, double z, double* sigma_8, double* sigma_8_cb)
    int nonlinear_hmcode_sigmadisp_at_z(void* pba, void* pnl, double z, double* sigma_disp, double* sigma_disp_cb)
    int nonlinear_hmcode_sigmadisp100_at_z(void* pba, void* pnl, double z, double* sigma_disp_100, double* sigma_disp_100_cb)
    int nonlinear_hmcode_sigmaprime_at_z(void* pba, void* pnl, double z, double* sigma_prime, double* sigma_prime_cb)
    int nonlinear_hmcode_window_nfw(void* pnl, double k, double rv, double c, double* window_nfw)

    int nonlinear_k_nl_at_z(void* pba, void* pnl, double z, double* k_nl, double* k_nl_cb)

    int spectra_firstline_and_ic_suffix(void *ppt, int index_ic, char first_line[_LINE_LENGTH_MAX_], FileName ic_suffix)

    int spectra_sigma(
                  void * pba,
                  void * ppm,
                  void * psp,
                  double R,
                  double z,
                  double * sigma)

    int spectra_sigma_cb(
                  void * pba,
                  void * ppm,
                  void * psp,
                  double R,
                  double z,
                  double * sigma_cb)

    int spectra_fast_pk_at_kvec_and_zvec(
                  void * pba,
                  void * psp,
                  double * kvec,
                  int kvec_size,
                  double * zvec,
                  int zvec_size,
                  double * pk_tot_out,
                  double * pk_cb_tot_out,
                  int nonlinear)

    double get_z_of_chi(double chi, void * pclass_sz, void * pba)

    double get_scale_dependent_bias_at_z_and_k(double z_asked,
                                               double k_asked,
                                               double bh,
                                               void * pclass_sz)
    double get_pk_nonlin_at_k_and_z(double k,
                                    double z,
                                    void * pba,
                                    void * ppm,
                                    void * pnl,
                                    void * pclass_sz)
    double get_pk_lin_at_k_and_z(double k,
                                    double z,
                                    void * pba,
                                    void * ppm,
                                    void * pnl,
                                    void * pclass_sz)

    double get_pk_nonlin_at_k_and_z_fast(double k,
                                    double z,
                                    void * pba,
                                    void * ppm,
                                    void * pnl,
                                    void * pclass_sz)
    double get_pk_lin_at_k_and_z_fast(double k,
                                    double z,
                                    void * pba,
                                    void * ppm,
                                    void * pnl,
                                    void * pclass_sz)
    double get_P_delta_at_m_and_z_b12(double m_asked,
                                       double z_asked,
                                       void * pclass_sz,
                                       void * pba)

    double get_rho_2h_at_r_and_m_and_z(double r_asked,
                                       double m_asked,
                                       double z_asked,
                                       void * pclass_sz,
                                       void * pba)

    double get_gas_pressure_2h_at_r_and_m_and_z(double r_asked,
                                                double m_asked,
                                                double z_asked,
                                                void * pclass_sz,
                                                void * pba)

    double get_mass_profile_at_x_M_z_nfw_200m(double x_asked,
                                             double m_asked,
                                             double z_asked,
                                             void * pba,
                                             void * pclass_sz)


    double get_gas_profile_at_x_M_z_nfw_200m(double x_asked,
                                             double m_asked,
                                             double z_asked,
                                             void * pba,
                                             void * pclass_sz)

    double get_gas_profile_at_x_M_z_bcm_200c(double r_asked,
                                      double m_asked,
                                      double z_asked,
                                      void * pba,
                                      void * pclass_sz)
    double get_planck_sigma_at_theta500(double theta500, void * pclass_sz)
    double get_szcountsz_sigma_at_theta_in_patch(double theta500, int id, void * pclass_sz)

    double get_rvir_of_m200c_at_z(double m_asked,
                                  double z_asked,
                                  void * pba,
                                  void * pclass_sz)
    double get_gas_profile_at_x_M_z_nfw_200c(double x_asked,
                                             double m_asked,
                                             double z_asked,
                                             void * pba,
                                             void * pclass_sz)

    double get_lensing_noise_at_ell(double l,
                                    void * pclass_sz)

    double evaluate_truncated_nfw_profile(double z,
                                          double k,
                                          double r_delta,
                                          double c_delta,
                                          double xout)

    double get_mean_galaxy_bias_at_z(double z, void * pclass_sz)



    double  get_szcounts_dndzdq_at_z_q(double z_asked, double qobs_asked, void * pclass_sz)


    double get_f_tinker10_at_nu_and_z(double nu, double z,void * pclass_sz)
    double get_T10_alpha_at_z(double z,void * pclass_sz)
    double get_f_tinker08_at_nu_and_z(double nu, double z, void * pclass_sz)

    double get_f_of_sigma_at_m_and_z(double m,
                                     double z,
                                     void * pba,
                                     void * pnl,
                                     void * pclass_sz)

    double get_gas_profile_at_x_M_z_b16_200c(double x_asked,
                                             double m_asked,
                                             double z_asked,
                                             double c_asked,
                                             double A_rho0,
                                             double A_alpha,
                                             double A_beta,
                                             double alpha_m_rho0,
                                             double alpha_m_alpha,
                                             double alpha_m_beta,
                                             double alpha_z_rho0,
                                             double alpha_z_alpha,
                                             double alpha_z_beta,
                                						 double mcut,
                                						 double alphap_m_rho0,
                                						 double alphap_m_alpha,
                                						 double alphap_m_beta,
                                						 double alpha_c_rho0,
                                						 double alpha_c_alpha,
                                						 double alpha_c_beta,
                                             double gamma,
                                             double xc,
                                             void * pba,
                                             void * pclass_sz)

    double get_f_b()
    double get_mu_e()
    double get_f_free()

    double get_rho_crit_at_z(double z_asked,
                             void * pba,
                             void * pclass_sz)

    double Delta_c_of_Omega_m(double Omega_m)

    double evaluate_rvir_of_mvir(double mvir,
                                 double delc,
                                 double rhoc,
                                 void * pclass_sz)

    double evaluate_cvir_of_mvir(double mvir,
                                 double z,
                                 void * pclass_sz,
                                 void * pba)

    double radial_kernel_W_galaxy_lensing_at_z(double z,
                                               void * pclass_sz)

    double get_source_galaxy_number_counts(double z,
                                           void * pclass_sz)
    double get_dkappacmbdz_at_l_and_z(double l,
                                      double z,
                                      void * pba,
                                      void * ppm,
                                      void * pnl,
                                      void * pclass_sz)

    double get_dkappacmbdz_pklin_at_l_and_z(double l,
                                            double z,
                                            void * pba,
                                            void * ppm,
                                            void * pnl,
                                            void * pclass_sz)
                                            
    double get_dyl2dzdlnm_at_l_z_and_m(double l,
                                       double z,
                                       double m,
                                       void * pba,
                                       void * pnl,
                                       void * pclass_sz)


    double get_dyldzdlnm_at_l_z_and_m(double l,
                                      double z,
                                      double m,
                                      void * pba,
                                      void * pnl,
                                      void * pclass_sz)

    double get_dygldzdlnm_at_l_z_and_m(double l,
                                      double z,
                                      double m,
                                      void * pba,
                                      void * pnl,
                                      void * pclass_sz)


    double get_m200m_to_m200c_at_z_and_M(double z_asked,
                                         double m_asked,
                                         void * pclass_sz)

    double get_normalization_gas_density_profile(double z_asked,
                                                 double m_asked,
                                                 void * pclass_sz)
    double get_m_to_xout_at_z_and_m(double z_asked,
                                    double m_asked,
                                    void * pclass_sz)

    double get_c200m_at_m_and_z_D08(double M, double z)
    double get_c200c_at_m_and_z_D08(double M, double z)
    double get_c200c_at_m_and_z_B13(double M, double z, void * ba, void * tsz)
    double get_c200m_at_m_and_z(double M, double z, void * ba, void * tsz)
    double get_c200c_at_m_and_z(double M, double z, void * ba, void * tsz)
    double get_c500c_at_m_and_z(double M, double z, void * ba, void * tsz)

    double radial_kernel_W_galaxy_ngal_at_z(int index_g,
                                            double * pvecback,
                                            double z,
                                            void * ba,
                                            void * tsz)

    double m_nfw(double x)

    double get_vrms2_at_z(double z,
                          void * pclass_sz)

    double get_hmf_counter_term_nmin_at_z(double z_asked,
                                          void * tsz)

    double get_hmf_counter_term_b1min_at_z(double z_asked,
                                          void * tsz)

    double get_hmf_counter_term_b2min_at_z(double z_asked,
                                          void * tsz)

    double get_volume_at_z(double z, void * pba)
    double get_IA_of_z(double z, void * pba, void * pclass_sz)
    double get_delta_mean_from_delta_crit_at_z(double delta_crit,double z,void * pclass_sz)
    double get_dndlnM_at_z_and_M(double z_asked,
                                 double m_asked,
                                 void * tsz)
    double gnu_tsz_of_nu_in_ghz(double nu_in_ghz,double Tcmb)
    double get_dcib0dz_at_z_and_nu(double z_asked,
                                   double nu_asked,
                                   void * tsz)
    double get_galaxy_number_counts(double z,void * tsz)

    double get_dydz_at_z(double z_asked,
                         void * tsz)
    double get_m200m_to_m500c_at_z_and_M(double z_asked,
                                         double m_asked,
                                         void * tsz)

    double get_m200c_to_m500c_at_z_and_M(double z_asked,
                                         double m_asked,
                                         void * tsz)
    double get_m200m_to_m200c_at_z_and_M(double z_asked,
                                         double m_asked,
                                         void * tsz)
    double get_m200c_to_m200m_at_z_and_M(double z_asked,
                                         double m_asked,
                                         void * tsz)
    double get_m500c_to_m200c_at_z_and_M(double z_asked,
                                         double m_asked,
                                         void * tsz)
    double get_m500c_to_m200m_at_z_and_M(double z_asked,
                                         double m_asked,
                                         void * tsz)

    double get_gas_density_profile_at_k_M_z(double l_asked,
                                        double m_asked,
                                        double z_asked,
                                        int normalize,
                                        void * tsz)

    double get_te_of_m500c_at_z_arnaud(double m, double z, void * pba, void * pclass_sz)
    double get_te_of_m500c_at_z_lee(double m, double z, void * pba, void * pclass_sz)


    double get_1e6xdy_from_battaglia_pressure_at_x_z_and_m200c(double z,
                                                               double m,
                                                               double x,
                                                               void * pba,
                                                               void * pclass_sz)

    double get_upp_from_gnfw_pressure_at_x_z_and_m500c(double z,
                                                           double m,
                                                           double x,
                                                           double d,
                                                           void * pba,
                                                           void * pclass_sz)




    double get_1e6xdy_from_gnfw_pressure_at_x_z_and_m500c(double z,
                                                           double m,
                                                           double x,
                                                           double d,
                                                           void * pba,
                                                           void * pclass_sz)


    double get_pressure_P_over_P_delta_at_x_M_z_b12_200c(double x_asked,
                                                         double m_asked,
                                                         double z_asked,
                                                         double c_asked,
                                                         double A_P0,
                                                         double A_xc,
                                                         double A_beta,
                                                         double alpha_m_P0,
                                                         double alpha_m_xc,
                                                         double alpha_m_beta,
                                                         double alpha_z_P0,
                                                         double alpha_z_xc,
                                                         double alpha_z_beta,
                                      							     double mcut,
                                      							     double alphap_m_P0,
                                      							     double alphap_m_xc,
                                      							     double alphap_m_beta,
                                      							     double alpha_c_P0,
                                      							     double alpha_c_xc,
                                      							     double alpha_c_beta,
                                                         double alpha,
                                                         double gamma,
                                                         void * pba,
                                                         void * tsz)

    double get_pressure_P_over_P_delta_at_x_gnfw_500c(double x_asked,
                                                          double P0GNFW,
                                                          double alphaGNFW,
                                                          double betaGNFW,
                                                          double gammaGNFW,
                                                          double c500,
                                                          void * pba,
                                                          void * tsz)

    double get_second_order_bias_at_z_and_nu(double z_asked,
                                             double nu_asked,
                                             void * tsz,
                                             void * pba)


    double get_cib_Snu_z_and_nu(double z_asked,
                                double nu_asked,
                                void * tsz)

    double get_first_order_bias_at_z_and_nu(double z_asked,
                                             double nu_asked,
                                             void * tsz)

    double get_y_at_m_and_z(double m, double z, void * pclass_sz, void * pba)
    double get_theta_at_m_and_z(double m, double z, void * pclass_sz, void * pba)
    double get_sigma_at_z_and_m(double z_asked,
                                double m_asked,
                                void * tsz,
                                void * pba)
    double get_dlnsigma_dlnR_at_z_and_m(double z_asked,
                                double m_asked,
                                void * tsz,
                                void * pba)
    double get_sigma8_at_z(double z_asked,
                          void * tsz,
                          void * pba)
    double get_nl_index_at_z_and_k(double z_asked,
                                    double k,
                                    void * tsz,
                                    void * nl)
    double get_nl_index_at_z_and_k_no_wiggles(double z_asked,
                                    double k,
                                    void * tsz,
                                    void * nl)
    double bispectrum_f2_kernel(double k1, double k2, double k3)

    double HOD_mean_number_of_satellite_galaxies(double z,
                                                 double M_halo,
                                                 double Nc_mean,
                                                 double M_min,
                                                 double alpha_s,
                                                 double M1_prime,
                                                 void * pclass_sz,
                                                 void * pba)

    double HOD_mean_number_of_central_galaxies(double z,
                                               double M_halo,
                                               double M_min,
                                               double sigma_log10M,
                                               double fcen,
                                               void * pclass_sz,
                                               void * pba)


    double get_fstar_of_m_at_z(double m,
                               double z,
                               void * pclass_sz)


    double get_matter_bispectrum_at_z_tree_level_PT(double k1_in_h_over_Mpc,
                                                         double k2_in_h_over_Mpc,
                                                         double k3_in_h_over_Mpc,
                                                         double z,
                                                         void * pclass_sz,
                                                         void * pba,
                                                         void * pnl,
                                                         void * ppm)


    double get_ttg_bispectrum_at_z_tree_level_PT(double k1_in_h_over_Mpc,
                                                         double k2_in_h_over_Mpc,
                                                         double k3_in_h_over_Mpc,
                                                         double z,
                                                         void * pclass_sz,
                                                         void * pba,
                                                         void * pnl,
                                                         void * ppm)
    double get_ttg_bispectrum_at_z_effective_approach(double k1_in_h_over_Mpc,
                                                         double k2_in_h_over_Mpc,
                                                         double k3_in_h_over_Mpc,
                                                         double z,
                                                         void * pclass_sz,
                                                         void * pba,
                                                         void * pnl,
                                                         void * ppm)

    double get_matter_bispectrum_at_z_effective_approach_smoothed(double k1_in_h_over_Mpc,
                                                         double k2_in_h_over_Mpc,
                                                         double k3_in_h_over_Mpc,
                                                         double z,
                                                         void * pclass_sz,
                                                         void * pba,
                                                         void * pnl,
                                                         void * ppm)

    double get_matter_bispectrum_at_z_effective_approach(double k1_in_h_over_Mpc,
                                                         double k2_in_h_over_Mpc,
                                                         double k3_in_h_over_Mpc,
                                                         double z,
                                                         void * pclass_sz,
                                                         void * pba,
                                                         void * pnl,
                                                         void * ppm)
    double get_matter_bispectrum_at_z_effective_approach_SC(double k1_in_h_over_Mpc,
                                                         double k2_in_h_over_Mpc,
                                                         double k3_in_h_over_Mpc,
                                                         double z,
                                                         void * pclass_sz,
                                                         void * pba,
                                                         void * pnl,
                                                         void * ppm)

    double evaluate_mean_galaxy_number_density_at_z(double z,
                                                    void * pclass_sz)

    double evaluate_unwise_m_min_cut(double z,
                                     int sample_id,
                                     void * pclass_sz)


    double get_nu_at_z_and_m(double z_asked,
                             double m_asked,
                             void * tsz,
                             void * pba)

    double get_radial_kernel_W_custom1_at_z(double z,
                                            void * pclass_sz)

    double get_custom1_profile_at_k_m_z(double k_asked,
                                        double m_asked,
                                        double z_asked,
                                        void * pclass_sz)

    double get_custom1_profile_at_x_m_z(double x_asked,
                                        double m_asked,
                                        double z_asked,
                                        void * pclass_sz)
                                        
    double get_b_custom1_at_z(double z,
                              void * pclass_sz)


    void fftlog_ComputeXiLMsloz(int l, 
                                 int m, 
                                 int N, 
                                 double * k, 
                                 double * pk, 
                                 double * r, 
                                 double * xi, 
                                 void * pclass_sz)


    double get_n5k_z_of_chi(double chi,
                            void * pclass_sz)

    double get_n5k_pk_at_z_and_k(double z_asked, 
                                double k_asked, 
                                void * pclass_sz)

    double get_n5k_cl_K1_at_chi(double chi,
                                void * pclass_sz)   

    int load_n5k_cl_K1(void * pclass_sz)
    int load_n5k_pk_zk(void * pclass_sz)
    int load_n5k_z_of_chi(void * pclass_sz)
