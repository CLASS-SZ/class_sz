
# $ python3 compute_and_plot_class_sz.py -param_name 'h' -p_val '[0.6732]'  -show_legend yes -show_error_bars no -output 'kSZ_kSZ_gal_3h,kSZ_kSZ_gal fft (3h)' -plot_kSZ_kSZ_gal yes  -save_figure yes -plot_redshift_dependent_functions no -mode run -x_min 1e1 -x_max 5e3
import argparse
couleur = 'green'
use_hod = 'yes'

import argparse
import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from scipy import interpolate
from scipy import stats
from datetime import datetime
import ast
import itertools



hplanck=6.62607004e-34 #m2 kg / s
kb = 1.38064852e-23 #m2 kg s-2 K-1
Tcmb =  2.725 #K

def g_nu(nu_in_GHz):
    nu_in_Hz = nu_in_GHz*1e9
    x = (hplanck*nu_in_Hz/kb/Tcmb)
    return (x/np.tanh(x/2.)-4.)


path_to_class_external_data = '/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts'
path_to_class = '/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/'
# path_to_class = '/Users/boris/Downloads/class_sz-6c5cb2cde2445ecf3e308fa8e90e10af4adf245f/'
# path_to_class = '/Users/boris/Downloads/class_sz/'
FIG_DIR = '/Users/boris/Work/CLASS-SZ/SO-SZ/figures'
# path_to_class = '/Users/boris/Downloads/class_sz-5bec681f30dad5f1b3529db7ffb4504764d77ffb/'
# path_to_class = '/Users/boris/Downloads/class_sz-9abe94d3f153d4a5c11599dd133f7bbf8abaaea7/'

# path_to_class = '/Users/boris/Downloads/class_sz-48ad2b333660ffeef8eb71bef208f59674aeef79/'
# path_to_class = '/Users/boris/Downloads/class_sz-d563960cb26a384462256a0a0d366829ab95ee04/'
# path_to_class = '/Users/boris/Downloads/class_sz-97ef4ef8feb782552a4897276b376e35e7f467c6/' # old version paper kusiak et al
# path_to_class = '/Users/boris/Downloads/class_sz-b008e6910a5dad0b1b0bdf8b9aa690dfdeb26de1/' # different

# path_to_class = '/Users/boris/Downloads/class_sz-837d50e0431028daf87b1d1badf51b60be4efac8/' # same as kusiak

# path_to_class = '/Users/boris/Downloads/class_sz-2067e49bf3eeca6b69a1abb42881c69e9993471f/' # slow and same as kusiak --> made fast
# path_to_class = '/Users/boris/Downloads/class_sz-cf17bc6bae8a666e3ec4173dd905ccc3d4dfe9ba/' # slow
# path_to_class = '/Users/boris/Downloads/class_sz-ad509c73e07b0652eee25716848348bb3d826be2/' # slow
# path_to_class = '/Users/boris/Downloads/class_sz-759fd6d3ab7054576f6cc025eacd6034e67739ff/' # fast




def scientific_notation(p_value,digit=2):
    str_xinj_asked = str("%.3e"%p_value)
    text_gamma_str1 = ''
    if p_value>1.:
        num = float(str_xinj_asked.split('e+')[0])
        exp = int(str_xinj_asked.split('e+')[1])
        if (round(num,digit)==10.):
            num = 1.
            exp = exp + 1
        if digit == 1:
            text_gamma_str1 = r'$%.1f \times 10^{%d}$'% (num,exp)
        if digit == 2:
            text_gamma_str1 = r'$%.2f \times 10^{%d}$'% (num,exp)
        if digit == 0:
            text_gamma_str1 = r'$%.0f \times 10^{%d}$'% (num,exp)
        if num == 1.:
            text_gamma_str1 = r'$10^{%d}$'% (exp)
    if p_value==0:
        text_gamma_str1 = r'$0$'
    elif p_value<1.:
        num = float(str_xinj_asked.split('e-')[0])
        exp = int(str_xinj_asked.split('e-')[1])
        if (round(num,digit)==10.):
            num = 1.
            exp = exp - 1
        if digit == 1:
            text_gamma_str1 = r'$%.1f \times 10^{-%d}$'% (num,exp)
        if digit == 2:
            text_gamma_str1 = r'$%.2f \times 10^{-%d}$'% (num,exp)
        if digit == 0:
            text_gamma_str1 = r'$%.0f \times 10^{-%d}$'% (num,exp)
        if num == 1.:
            text_gamma_str1 = r'$10^{-%d}$'% (exp)
    if p_value==1.:
        text_gamma_str1 = r'$1$'

    return text_gamma_str1



r_dict = {} #dictionary of results
def run(args):
    # if(args.plot_ref_data == 'yes'):
    #     parameter_file = 'class-sz_chill_B12_parameters.ini'
    # else:
    #     parameter_file = 'class-sz_parameters.ini'

    # important parameters are re-ajusted later, here we just load a template file:
    # parameter_file = 'class-sz_parameters_KFSW20.ini'
    #'class-sz_chill_B12_parameters.ini'
    #parameter_file = 'class-sz_parameters_rotti++20.ini'
    #parameter_file ='tSZ_params_ref_resolved-rotti++20_snr12_step_2.ini'
    # parameter_file = 'class-sz_chill_B12_parameters.ini'
    #parameter_file = 'class-sz_parameters_for_sz_ps_completeness_analysis_rotti++20.ini'


    os.chdir(path_to_class)
    #collect arguments
    param_name = args.param_name
    print(param_name)

    if (param_name == 'sigma8'):
        label_key = r'$\sigma_8$'
    elif (param_name == 'B'):
        label_key = r'$B$'
    elif (param_name == 'h'):
        label_key = r'$h$'
    elif (param_name == 'Omega_cdm'):
        label_key = r'$\Omega_\mathrm{m}$'
    elif (param_name == 'M_max'):
        label_key = r'$M_\mathrm{max}$'
    elif (param_name == 'z_max'):
        label_key = r'$z_\mathrm{max}$'
    elif (param_name == 'M_min'):
        label_key = r'$M_\mathrm{min}$'
    elif (param_name == 'm_ncdm'):
        label_key = r'$\Sigma m_\mathrm{\nu}$'
    elif (param_name == 'signal-to-noise cut-off for ps completeness analysis'):
        label_key = r'$q_\mathrm{cut}$'
    else:
        label_key = 'val'


    if args.p_val:
        p = np.asarray(ast.literal_eval(args.p_val))
        args.N = len(p)
        N = args.N
    else:
        p_min = float(args.p_min)
        p_max = float(args.p_max)
        N = int(args.N)
        spacing = args.spacing
        #define array of parameter values
        if(spacing=='log'):
            p = np.logspace(np.log10(p_min),np.log10(p_max),N)
        elif(spacing=='lin'):
            p = np.linspace(p_min,p_max,N)


    # p_min = float(args.p_min)
    # p_max = float(args.p_max)
    # N = int(args.N)
    # spacing = args.spacing


    #array to store results
    cl_1h = []
    trispectrum = []
    multipoles = []
    col = []
    val_label = []
    y_err = []
    cl_2h = []
    te_y_y = []
    kSZ_kSZ_gal_1h = []
    kSZ_kSZ_gal_1h_fft = []
    kSZ_kSZ_gal_2h_fft = []
    kSZ_kSZ_gal_3h_fft = []
    kSZ_kSZ_tSZ_1h = []
    kSZ_kSZ_tSZ_2h = []
    kSZ_kSZ_tSZ_3h = []
    kSZ_kSZ_gal_2h = []
    kSZ_kSZ_gal_3h = []
    kSZ_kSZ_gal_hf = []
    kSZ_kSZ_lensmag_1h = []
    tSZ_tSZ_tSZ_1h = []
    tSZ_gal_1h = []
    tSZ_gal_2h = []
    gal_gal_1h = []
    gal_gal_2h = []
    gal_gal_hf = []
    gal_lens_1h = []
    gal_lens_2h = []
    gal_lens_hf = []
    tSZ_lens_1h = []
    tSZ_lens_2h = []
    isw_lens = []
    isw_tsz = []
    isw_auto = []
    lens_lens_1h = []
    lens_lens_2h = []
    gal_lensmag_1h = []
    gal_lensmag_2h = []
    gal_gallens_1h = []
    gal_gallens_2h = []
    gal_lensmag_hf = []
    lensmag_lensmag_1h = []
    lensmag_lensmag_2h = []
    lensmag_lensmag_hf = []
    lens_lensmag_1h = []
    lens_lensmag_2h = []
    lens_lensmag_hf = []
    redshift_dependent_functions_z = []
    redshift_dependent_functions_Q = []

    #
    # #define array of parameter values
    # if(spacing=='log'):
    #     p = np.logspace(np.log10(p_min),np.log10(p_max),N)
    # elif(spacing=='lin'):
    #     p = np.linspace(p_min,p_max,N)

    # #load template parameter file into dictionnary
    p_dict = {}
    # with open(path_to_class_external_data+"/run_scripts/" + parameter_file) as f:
    #     for line in f:
    #         x = line.strip()
    #         if not x.startswith("#"):
    #             print(x.split('='))
    #             (key, val) = x.split('=')
    #             (key, val) = (key.strip(), val.strip())
    #             p_dict[key] = val

    # p_dict['path_to_class'] = "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz"
    # P18 1st column table 1 of https://arxiv.org/pdf/1807.06209.pdf
    p_dict['multipoles_sz'] = 'ell_mock'
    p_dict['omega_b'] = 0.022383
    p_dict['omega_cdm'] = 0.12011
    p_dict['h'] = 0.6732
    p_dict['tau_reio'] = 0.0543
    p_dict['ln10^{10}A_s'] = 3.0448
    p_dict['n_s'] = 0.96605
    p_dict['k_pivot'] = 0.05
    p_dict['N_ncdm'] = 1
    p_dict['N_ur'] = 2.0328
    p_dict['m_ncdm'] = 0.06
    p_dict['f_free'] = 1.


    p_dict['n_ell_density_profile'] =100
    p_dict['n_m_density_profile'] = 100 # 80
    p_dict['n_z_density_profile'] = 100 # 80

    p_dict['n_z_psi_b1g'] = 100
    p_dict['n_l_psi_b1g'] = 400

    p_dict['n_z_psi_b2g'] = 100
    p_dict['n_l_psi_b2g'] = 400

    p_dict['n_z_psi_b2t'] = 100
    p_dict['n_l_psi_b2t'] = 400

    p_dict['n_z_psi_b1t'] = 1000
    p_dict['n_l_psi_b1t'] = 100

    p_dict['n_z_psi_b1gt'] = 100
    p_dict['n_l_psi_b1gt'] = 100

    p_dict['N_samp_fftw'] = 2000
    p_dict['l_min_samp_fftw'] = 1e-9
    p_dict['l_max_samp_fftw'] = 1e9
    # p_dict['P0_B12'] = 30.
    #p_dict['halofit_k_per_decade'] = 3000.
    # p_dict['P_k_max_h/Mpc'] = 50.

    # set the k grid for compuation of sigma and dsigma used in HMF
    p_dict['k_min_for_pk_class_sz'] = 1E-4
    p_dict['k_max_for_pk_class_sz'] = 1E2
    p_dict['k_per_decade_class_sz'] = 20.
    p_dict['P_k_max_h/Mpc'] = p_dict['k_max_for_pk_class_sz']/p_dict['h']

    # p_dict["dlnk_for_pk_hm"] = 0.1
    # p_dict["z_for_pk_hm"] = 2.
    # p_dict["z_pk"] = p_dict["z_for_pk_hm"]
    # p_dict["k_min_for_pk_hm"] = 1e-4
    # p_dict["k_max_for_pk_hm"] = 2e1

    p_dict['non linear'] = 'halofit'

    p_dict['z_min'] = 1e-2
    p_dict['hm_consistency'] = 1
    p_dict['check_consistency_conditions'] = 1
    p_dict['M_min'] = 1e10#*p_dict['h']
    p_dict['M_max'] = 5e15#*p_dict['h']
    # p_dict['HMF_prescription_NCDM'] = 'No-pres'
    # p_dict['mass function'] = 'T10'  #fiducial  T10
    p_dict['mass function'] = 'T08M200c'  #fiducial  T10
    p_dict['galaxy_sample'] = "unwise"
    #p_dict['full path to dndz (normalized galaxy dist.)'] = "/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/run_scripts/yxg/data/dndz/unwise_"+couleur+".txt"

    p_dict['galaxy_sample_id'] = couleur

    p_dict['galaxy_sample'] = 'custom'
    p_dict['full_path_to_dndz_gal'] = '/Users/boris/Work/DES/nz_maglim_z_bin1.txt'
    p_dict['full_path_to_source_dndz_gal'] = '/Users/boris/Work/DES/nz_maglim_z_bin3.txt'
    p_dict['convert_cls_to_gamma'] = 1

    p_dict['concentration parameter'] = 'D08'


    # p_dict['M_min_HOD_mass_factor_unwise'] = 1.03
    # p_dict['use_analytical_truncated_nfw'] = 'yes'

    p_dict['x_out_truncated_nfw_profile_satellite_galaxies'] = 1

    p_dict['M0 equal M_min (HOD)'] = 'no'


    p_dict['gas profile'] = 'nfw' # 'nfw' or 'B16'
    p_dict['gas profile mode'] = 'agn'



    p_dict['sigma_log10M_HOD'] = 0.76
    p_dict['alpha_s_HOD'] = 2.08
    p_dict['M_min_HOD'] = 1.01e13
    p_dict['M1_prime_HOD'] = 1.18e14
    p_dict['M0_HOD'] = 0.

    # p_dict['x_out_nfw_profile'] = 2. # default 1
    p_dict['x_out_truncated_nfw_profile'] = 1.
    p_dict['pk_nonlinear_for_vrms2'] = 1
    # p_dict['hm_consistency'] =  0: nothing 1: counter terms 2: alpha(z)

    p_dict['delta for electron density'] ='200c'
    p_dict['delta for galaxies'] ='200m'
    p_dict['delta for matter density'] ='200m' # not used

    p_dict['sz_verbose'] = 3
    p_dict['root'] = 'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_'
    p_dict['write sz results to files'] = 'yes' # this writes  PS and f(z)

    p_dict['nfw_profile_epsabs'] = 1.e-33
    p_dict['nfw_profile_epsrel'] = 1.e-2

    p_dict['redshift_epsabs'] = 1.e-50 # fiducial 1e-30
    p_dict['redshift_epsrel'] = 1.e-8 # fiducial value 1e-8

    p_dict['mass_epsabs'] = 1.e-40 # fiducial 1e-30
    p_dict['mass_epsrel'] = 1e-5

    p_dict['dlogell'] = 0.1 # 0.1
    p_dict['ell_min'] = 2.
    p_dict['ell_max'] = 5000.


    # p_dict['N_kSZ2_gal_multipole_grid'] =  20 # fiducial 10

    p_dict['non linear'] = 'halofit'
    p_dict['nonlinear_verbose']= 1
    p_dict['ndim_masses'] = 150 # important 128 is default ccl value
    p_dict['ndim_redshifts'] = 150


    # p_dict['perturb_sampling_stepsize'] = 0.01
    # p_dict['k_max_tau0_over_l_max']=10.
    # p_dict['l_max_scalars'] = 5000




    if(args.output):
        print(args.output)
        p_dict['output'] = args.output
    if ("lens_lens_" in p_dict['output']):
        p_dict['z_max'] = 20.
    else:
        p_dict['z_max'] = 4.

    # if ("kSZ" in p_dict['output']):
    #     # p_dict['nfw_profile_epsabs'] = 1.e-6
    #     # p_dict['nfw_profile_epsrel'] = 1.e-10
    #
    #     p_dict['redshift_epsabs'] = 1.e-33 # fiducial 1e-30
    #     p_dict['redshift_epsrel'] = 5.e-4 # fiducial value 1e-8
    #
    #     p_dict['mass_epsabs'] = 1.e-33 # fiducial 1e-30
    #     p_dict['mass_epsrel'] = 5e-4
    #     #p_dict['halo occupation distribution'] = 'KFSW20'
    #
    #     p_dict['dell'] = 100
    #     p_dict['ell_max_mock'] = 4000.
    #     p_dict['ell_min_mock'] = 4300.
    #
    #     # p_dict['m_min_counter_terms'] = 1e8
    #     p_dict['ksz_filter_file'] = path_to_class+'/sz_auxiliary_files/UNWISE_galaxy_distributions/unwise_filter_functions_l_fl.txt'
        # print('setting ksz file')
        # p_dict['ksz_filter_file'] = '/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/sz_auxiliary_files/UNWISE_galaxy_distributions/AdvACT_kSZfilt_ellmax8000_smoothed_tapered_nosqrt_w1p5arcminbeam.txt'

        #p_dict['A10_file'] = "class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_A10.txt"

        # exit(0)
        # p_dict['nfw_profile_epsabs'] = 1.e-50
        # p_dict['nfw_profile_epsrel'] = 1.e-4
    if ("kSZ_kSZ_tSZ" in p_dict['output']):
            p_dict['redshift_epsabs'] = 1.e-30 # fiducial 1e-30
            p_dict['redshift_epsrel'] = 1.e-2 # fiducial value 1e-8
            p_dict['mass_epsabs'] = 1.e-30 # fiducial 1e-30
            p_dict['mass_epsrel'] = 1.e-3

            p_dict['dlogell'] = 0.05 # 0.1
            p_dict['ell_min'] = 2.
            p_dict['ell_max'] = 5000.
            p_dict['n_z_hmf_counter_terms'] = 300
            p_dict['ndim_masses'] = 200 # important 128 is default ccl value
            p_dict['ndim_redshifts'] = 120
            p_dict['P_k_max_h/Mpc'] = 100.
            # set the k grid for compuation of sigma and dsigma used in HMF
            p_dict['k_min_for_pk_class_sz'] = 1E-3
            p_dict['k_max_for_pk_class_sz'] = 60.
            p_dict['k_per_decade_class_sz'] = 50 # at least 40?
            p_dict['Frequency for y-distortion in GHz'] = 143.
            p_dict['bispectrum_lambda_2'] = 1.
            p_dict['bispectrum_lambda_3'] = 0.01

    if ("kSZ_kSZ_gal" in p_dict['output']):
            # p_dict['nfw_profile_epsabs'] = 1.e-6
            # p_dict['nfw_profile_epsrel'] = 1.e-10
            # p_dict['tol_background_integration'] = 1e-10
            # p_dict['back_integration_stepsize'] = 1e-5


            p_dict['redshift_epsabs'] = 1.e-30 # fiducial 1e-30
            p_dict['redshift_epsrel'] = 1.e-2 # fiducial value 1e-8

            p_dict['N_kSZ2_gal_multipole_grid'] =  70#70 # fiducial 70
            p_dict['N_kSZ2_gal_theta_grid'] =  70#70 # fiducial 70
            p_dict['ell_min_kSZ2_gal_multipole_grid'] = 2.
            p_dict['ell_max_kSZ2_gal_multipole_grid'] = 2e5

            p_dict['mass_epsabs'] = 1.e-30 # fiducial 1e-30
            p_dict['mass_epsrel'] = 1.e-3
            #p_dict['halo occupation distribution'] = 'KFSW20'

            p_dict['dlogell'] = 0.05 # 0.1
            # p_dict['dell'] = 100.
            # p_dict['dell'] = 100.
            p_dict['ell_min'] = 10.
            p_dict['ell_max'] = 5000.





            p_dict['n_z_hmf_counter_terms'] = 300

            # p_dict['m_min_counter_terms'] = 1e8
            p_dict['ksz_filter_file'] = path_to_class+'/sz_auxiliary_files/UNWISE_galaxy_distributions/unwise_filter_functions_l_fl.txt'
            # print('setting ksz file')
            # p_dict['ksz_filter_file'] = '/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/sz_auxiliary_files/UNWISE_galaxy_distributions/AdvACT_kSZfilt_ellmax8000_smoothed_tapered_nosqrt_w1p5arcminbeam.txt'

            #p_dict['A10_file'] = "class_sz_lnIgnfw-and-d2lnIgnfw-vs-lnell-over-ell500_A10.txt"

            # exit(0)
            # p_dict['nfw_profile_epsabs'] = 1.e-50
            # p_dict['nfw_profile_epsrel'] = 1.e-4

            p_dict['ndim_masses'] = 200 # important 128 is default ccl value
            p_dict['ndim_redshifts'] = 120
            #
            p_dict['P_k_max_h/Mpc'] = 100.


            # set the k grid for compuation of sigma and dsigma used in HMF
            p_dict['k_min_for_pk_class_sz'] = 1E-3
            p_dict['k_max_for_pk_class_sz'] = 60.
            p_dict['k_per_decade_class_sz'] = 50 # at least 40?


            p_dict['write sz results to files'] = 'yes'

    if ("pk" or "bk" in p_dict['output']):
        # exit(0)
        p_dict['z_max_pk'] = p_dict['z_max']
        p_dict["dlnk_for_pk_hm"] = 0.01
        p_dict["z_for_pk_hm"] = 1.
        p_dict["z_pk"] = p_dict["z_for_pk_hm"]
        p_dict["k_min_for_pk_hm"] = 1e-3
        p_dict["k_max_for_pk_hm"] = 1e1
        p_dict['nfw_profile_epsabs'] = 1.e-40
        p_dict['nfw_profile_epsrel'] = 1.e-5
    if (args.plot_bk == 'yes' or args.plot_bk_ttg == 'yes'):
        # p_dict['nfw_profile_epsabs'] = 1.e-6
        # p_dict['nfw_profile_epsrel'] = 1.e-10
        p_dict['omega_b'] = 0.0226
        p_dict['omega_cdm'] = 0.11
        p_dict['h'] = 0.71
        p_dict['tau_reio'] = 0.088
        p_dict.pop('ln10^{10}A_s')
        p_dict['A_s'] = 2.43e-9
        p_dict['n_s'] = 0.963
        p_dict['k_pivot'] = 0.002

    if "ttg" in p_dict['output']:
        p_dict['n_ell_density_profile'] = 200
        p_dict['n_m_density_profile'] = 200 # 80
        p_dict['n_z_density_profile'] = 200 # 80





    p_dict['z_max_pk'] = p_dict['z_max']


    bg_blue = 1.56
    bg_green = 2.23
    bg_red = 3.29

    smag_red = 0.842
    smag_green = 0.648
    smag_blue = 0.455

    # if(args.plot_ref_data == 'yes'):
    #     #L = np.loadtxt('/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/sz_auxiliary_files/chill_cltsz_data.txt')
    # L_ref = np.loadtxt('/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/sz_auxiliary_files/run_scripts/class-sz_szpowerspectrum_kSZ2_gal_ref.txt')
    # ell_ref = L_ref[:,0]
    # kSZ_kSZ_gal_1h_ref = L_ref[:,8]

    if couleur == 'red':
        L_ref_sf_july = np.loadtxt(path_to_class_external_data+'/kSZ2_gal_unwise/07-28-2020-l2Cl_filtered_planck_gilmarin_b1_YESJeans_kSZ2xgals-red_16.2.dat')
        L_ref_sf_august = np.loadtxt(path_to_class_external_data+'/kSZ2_gal_unwise/08-18-2020-l2Cl_filtered_planck_NLveldisp_b1_NOJeans_kSZ2xgals-red_16.2.dat')
        bgeff = bg_red
        smag = smag_red
        shot_noise = 29.6e-2



    elif couleur == 'green':
        L_ref_sf_july = np.loadtxt(path_to_class_external_data+'/kSZ2_gal_unwise/07-28-2020-l2Cl_filtered_planck_gilmarin_b1_YESJeans_kSZ2xgals-green.dat')
        L_ref_sf_august = np.loadtxt(path_to_class_external_data+'/kSZ2_gal_unwise/08-18-2020-l2Cl_filtered_planck_NLveldisp_b1_NOJeans_kSZ2xgals-green.dat')
        bgeff = bg_green
        smag = smag_green
        shot_noise = 1.81e-2

    elif couleur == 'blue':
        L_ref_sf_july = np.loadtxt(path_to_class_external_data+'/kSZ2_gal_unwise/07-28-2020-l2Cl_filtered_planck_gilmarin_b1_YESJeans_kSZ2xgals-blue.dat')
        L_ref_sf_august = np.loadtxt(path_to_class_external_data+'/kSZ2_gal_unwise/08-18-2020-l2Cl_filtered_planck_NLveldisp_b1_NOJeans_kSZ2xgals-blue.dat')
        bgeff = bg_blue
        smag = smag_blue
        shot_noise = 0.92e-2
        #print('bg_blue=%.3e'%bg_blue)

    t2g_1h_ref = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/class-sz_szpowerspectrum_1.0.txt')
    kSZ_kSZ_gal_1h_ref = t2g_1h_ref[:,8]
    kSZ_kSZ_gal_1h_ells_ref = t2g_1h_ref[:,0]

    t2g_2h_ref = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/class-sz_szpowerspectrum_1.0_2h.txt')
    kSZ_kSZ_gal_2h_ref = t2g_2h_ref[:,42]
    kSZ_kSZ_gal_2h_ells_ref = t2g_2h_ref[:,0]

    ell_ref_sf_august= L_ref_sf_august[0,:]
    kSZ_kSZ_gal_1h_ref_sf_august = L_ref_sf_august[1,:]

    ell_ref_sf_july = L_ref_sf_july[0,:]
    kSZ_kSZ_gal_1h_ref_sf_july = L_ref_sf_july[1,:]


    if couleur == 'red':
        L_ref_sf_july = np.loadtxt('/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/kSZ2_gal_unwise/07-28-2020-kappa_mag_l2Cl_filtered_planck_gilmarin_YESJeans_kSZ2xgals-red_16.2.dat')
        L_ref_sf_august = np.loadtxt('/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/kSZ2_gal_unwise/08-18-2020-kappa_mag_l2Cl_filtered_planck_NLveldisp_NOJeans_kSZ2xgals-red_16.2.dat')
        smageff = smag_red

    elif couleur == 'green':
        L_ref_sf_july = np.loadtxt(path_to_class_external_data+'/kSZ2_gal_unwise/07-28-2020-kappa_mag_l2Cl_filtered_planck_gilmarin_YESJeans_kSZ2xgals-green.dat')
        L_ref_sf_august = np.loadtxt(path_to_class_external_data+'/kSZ2_gal_unwise/08-18-2020-kappa_mag_l2Cl_filtered_planck_NLveldisp_NOJeans_kSZ2xgals-green.dat')
        smageff = smag_green

    elif couleur == 'blue':
        L_ref_sf_july = np.loadtxt(path_to_class_external_data+'/kSZ2_gal_unwise/07-28-2020-kappa_mag_l2Cl_filtered_planck_gilmarin_YESJeans_kSZ2xgals-blue.dat')
        L_ref_sf_august = np.loadtxt(path_to_class_external_data+'/kSZ2_gal_unwise/08-18-2020-kappa_mag_l2Cl_filtered_planck_NLveldisp_NOJeans_kSZ2xgals-blue.dat')
        smageff = smag_blue


    ell_ref_lensmag_sf_august= L_ref_sf_august[0,:]
    kSZ_kSZ_lensmag_1h_ref_sf_august = L_ref_sf_august[1,:]

    ell_ref_lensmag_sf_july = L_ref_sf_july[0,:]
    kSZ_kSZ_lensmag_1h_ref_sf_july = L_ref_sf_july[1,:]


    #     L_ref = np.loadtxt('/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/output/class-sz_tmp_szpowerspectrum_mnu_0d02_ref.txt')
    #     multipoles_ref = L_ref[:,0]
    #     cl_1h_ref = L_ref[:,1]


    #prepare the figure
    label_size = 12
    title_size = 15
    legend_size = 13
    handle_length = 1.5
    if (args.print_rel_diff == 'yes'):
        fig, (ax1,ax2) = plt.subplots(2,1,figsize=(7,7))
    else:
        fig, ax1 = plt.subplots(1,1,figsize=(7,5))
    ax = ax1
    ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
    ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=label_size)
    plt.setp(ax.get_xticklabels(), fontsize=label_size)
    ax.grid( b=True, which="both", alpha=0.3, linestyle='--')



    if (args.plot_redshift_dependent_functions == 'yes'):
        ax.set_xlabel(r'$z$',size=title_size)
        #ax.set_ylabel(r'$Q(z)$',size=title_size)
        ax.set_ylabel(r'$\bar{n}_g=\int_{M_{min}}^{M_{max}}dM dn/dM$',size=title_size)
        #ax.set_ylabel(r'$v_\mathrm{rms}\quad [\mathrm{km/s}]$',size=title_size)
        ax.set_xscale('log')
        ax.set_yscale('linear')

        #ax.set_xlim(0.,6.)
        # ax.set_ylim(0.,600.)
    else:
        ax.set_xscale('log')
        ax.set_yscale('linear')
        ax.set_xlabel(r'$\ell$',size=title_size)
        if (args.plot_trispectrum == 'yes'):
            ax.set_ylabel(r'$T_{\ell,\ell}$',size=title_size)
        elif (args.plot_te_y_y == 'yes'):
            ax.set_ylabel(r'$\mathrm{T_e^{tSZ}} \quad [\mathrm{keV}]$',size=title_size)
        elif (args.plot_kSZ_kSZ_gal == 'yes'):
            ax.set_ylabel(r'$\ell(\ell+1)C_\ell^\mathrm{^{kSZ^2-g}}/2\pi\quad [\mathrm{\mu K^2}]$',size=title_size)
        # elif (args.plot_kSZ_kSZ_gal_2h == 'yes'):
        #     ax.set_ylabel(r'$\ell(\ell+1)C_\ell^\mathrm{^{kSZ^2-g\,(2h)}}/2\pi\quad [\mathrm{\mu K^2}]$',size=title_size)
        # elif (args.plot_kSZ_kSZ_gal_3h == 'yes'):
        #     ax.set_ylabel(r'$\ell(\ell+1)C_\ell^\mathrm{^{kSZ^2-g\,(3h)}}/2\pi\quad [\mathrm{\mu K^2}]$',size=title_size)
        elif (args.plot_kSZ_kSZ_lensmag_1h == 'yes'):
            ax.set_ylabel(r'$\ell(\ell+1)C_\ell^\mathrm{^{kSZ^2-\mu}}/2\pi\quad [\mathrm{\mu K^2}]$',size=title_size)
        elif (args.plot_tSZ_tSZ_tSZ_1h == 'yes'):
            ax.set_ylabel(r'$\mathrm{b^{y-y-y}_{\ell_1,\ell_2,\ell_3}}$',size=title_size)
        elif (args.plot_tSZ_gal == 'yes'):
            ax.set_ylabel(r'$10^9 \ell(\ell+1)\mathrm{C^{yg}_\ell}/2\pi$',size=title_size)
        elif (args.plot_gal_gal == 'yes'):
            ax.set_ylabel(r'$10^{5}\times\mathrm{C^{gg}_\ell}$',size=title_size)
        elif (args.plot_gal_lens == 'yes'):
            ax.set_ylabel(r'$10^{5}\times\mathrm{C^{g\kappa}_\ell}$',size=title_size)
        elif (args.plot_gal_gallens == 'yes'):
            ax.set_ylabel(r'$10^{5}\times\mathrm{C^{g\gamma}_\ell}$',size=title_size)
        elif (args.plot_lens_lens == 'yes'):
            ax.set_ylabel(r'$[\ell(\ell+1)]^2\mathrm{C^{\phi\phi}_\ell}/2\pi$',size=title_size)
        elif (args.plot_tSZ_lens == 'yes'):
            ax.set_ylabel(r'$\ell^2(\ell+1)\mathrm{C^{y\phi}_\ell/2\pi}$ [$\mu$K]',size=title_size)
        elif (args.plot_isw_lens == 'yes'):
            ax.set_ylabel(r'$\ell(\ell+1)\mathrm{C^{ISW\times\phi}_\ell/2\pi}$',size=title_size)
        elif (args.plot_isw_tsz == 'yes'):
            ax.set_ylabel(r'$\ell(\ell+1)\mathrm{C^{ISW\times y}_\ell/2\pi}$',size=title_size)
        elif (args.plot_isw_auto == 'yes'):
            ax.set_ylabel(r'$\ell(\ell+1)\mathrm{C^{ISW\times ISW}_\ell/2\pi}$',size=title_size)
        elif (args.plot_pk == 'yes'):
            ax.set_xlabel(r'$k\,\, [h/\mathrm{Mpc}]$',size=title_size)
            ax.set_ylabel(r'$P(k)\,\, [(\mathrm{Mpc}/h)^3]$',size=title_size)
        elif (args.plot_bk == 'yes'):
            ax.set_xlabel(r'$k\,\, [h/\mathrm{Mpc}]$',size=title_size)
            ax.set_ylabel(r'$B(k,k,k)\,\, [(\mathrm{Mpc}/h)^6]$',size=title_size)
        elif (args.plot_b_kSZ_kSZ_tSZ == 'yes'):
            ax.set_xlabel(r'$\ell$',size=title_size)
            ax.set_ylabel(r'$B(\ell,\lambda_2 \ell,\lambda_3 \ell)\,\,TTY\,\, [(\mathrm{Mpc}/h)^6]$',size=title_size)
        elif (args.plot_bk_ttg == 'yes'):
            ax.set_xlabel(r'$k\,\, [h/\mathrm{Mpc}]$',size=title_size)
            ax.set_ylabel(r'$B(k,k,k) [ttg]\,\, [(\mathrm{Mpc}/h)^6]$',size=title_size)
        else:
            ax.set_ylabel(r'$10^{12}\ell(\ell+1)\mathrm{C^{yy}_\ell/2\pi}$',size=title_size)

    if(args.x_min):
        ax.set_xlim(left=float(args.x_min))
    if(args.x_max):
        ax.set_xlim(right=float(args.x_max))

    if(args.y_min):
        ax.set_ylim(bottom=float(args.y_min))
    if(args.y_max):
        ax.set_ylim(top=float(args.y_max))
    if(args.f_sky):
        f_sky = float(args.f_sky)
    else:
        f_sky = 1.

    #colors = iter(['b','green','orange'])
    colors = iter(cm.viridis(np.linspace(0., 0.9, N)))



    if (args.compute_scaling_with_param == 'yes'):
        #declare list to store cl at l=100
        cl_1h_100 = []
    if(int(args.N)>0):
        #loop over parameter values
        id_p = 0
        for p_val in p:
            #update dictionnary with current value
            p_dict[param_name] = p_val
            if args.mode == 'run':
                if id_p == 0:
                    subprocess.call(['rm','-r','-f',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])
                    subprocess.call(['mkdir','-p',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])
                with open(path_to_class+'sz_auxiliary_files/run_scripts/tmp/tmp.ini', 'w') as f:
                    for k, v in p_dict.items():
                        # print(str(k)+' = '+str(v))
                        f.write(str(k) + ' = '+ str(v) + '\n')
                startTime = datetime.now()
                subprocess.call(['./class',path_to_class+'sz_auxiliary_files/run_scripts/tmp/tmp.ini'])
                print("time in class -> " + str((datetime.now() - startTime)))
            #print(datetime.now() - startTime)


            if (param_name == 'Omega_cdm'):
                p_val = p_val + float(p_dict['Omega_b'])
            elif (param_name == 'm_ncdm'):
                p_val = 3.*p_val
            #val_label.append(label_key + ' = %.2e'%(p_val))
            if (param_name == 'm_ncdm'):
                val_label.append(label_key + ' = %3.0f'%(1000*p_val)+ ' meV')
            elif (param_name == 'use_central_hod'):
                val_label.append(label_key + ' = %d'%(p_val))
            elif (param_name == 'unwise_galaxy_sample_id'):
                val_label.append('%s'%(p_val))
            elif (param_name == 'use_hod'):
                val_label.append('%s'%(p_val))
            elif (param_name == 'M_min_HOD_mass_factor_unwise'):
                val_label.append(label_key + ' = ' + '%s'%(p_val))
            else:
                val_label.append(label_key + ' = ' + scientific_notation(p_val))
                #val_label.append(label_key + ' = ' + "%.2f"%(p_val))

            col.append(next(colors))


            if (args.plot_redshift_dependent_functions == 'yes'):
                R = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_redshift_dependent_functions.txt')
                redshift_dependent_functions_z.append(R[:,0])
                redshift_dependent_functions_Q.append(R[:,9]) # 6: v2rms, 7: sigma2_hsv, 8: m200m/m200c @ m200m = 10^{13.5} Msun/h
                if (val_label[id_p] == 'gr_shallow'):
                    ax.plot(redshift_dependent_functions_z[id_p],redshift_dependent_functions_Q[id_p],color='forestgreen',ls='--',alpha = 1.,label = val_label[id_p])
                elif (val_label[id_p] == 'green'):
                    ax.plot(redshift_dependent_functions_z[id_p],redshift_dependent_functions_Q[id_p],color='green',ls='-',alpha = 1.,label = val_label[id_p])
                elif (val_label[id_p] == 'blue'):
                    ax.plot(redshift_dependent_functions_z[id_p],redshift_dependent_functions_Q[id_p],color='blue',ls='-',alpha = 1.,label = val_label[id_p])
                elif (val_label[id_p] == 'red'):
                    ax.plot(redshift_dependent_functions_z[id_p],redshift_dependent_functions_Q[id_p],color='red',ls='-',alpha = 1.,label = val_label[id_p])
                else:
                    ax.plot(redshift_dependent_functions_z[id_p],redshift_dependent_functions_Q[id_p],color=col[id_p],ls='-',alpha = 1.,label = val_label[id_p])

            # if ('pk_at_z_1h' in p_dict['output']
            # or 'pk_at_z_2h' in p_dict['output']):
            if (args.plot_pk == 'yes'):
                C = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_halo_model_pk_at_z_k_pk1h_pk2h.txt')
                k_class_sz = C[:,0]
                pk_1h_class_sz = C[:,1]
                pk_2h_class_sz = C[:,2]
                pk_gg_1h_class_sz = C[:,3]
                pk_gg_2h_class_sz = C[:,4]
                C = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_pk_cb.dat')
                k_class = C[:,0]
                pk_linear_class = C[:,1]

                C = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_pk_cb_nl.dat')
                k_class = C[:,0]
                pk_nonlinear_class = C[:,1]
                ax.plot(k_class,pk_linear_class,color='r',ls=':',alpha = 1.,label = 'pk_linear')
                ax.plot(k_class,pk_nonlinear_class,color='r',ls='-.',alpha = 1.,label = 'pk_halofit')
                ax.plot(k_class_sz,pk_1h_class_sz,color='k',ls='--',alpha = 1.,label = 'pk_1h_class_sz')
                #ax.plot(k_class_sz,pk_2h_class_sz,color='k',ls='-.',alpha = 1.,label = 'pk_2h_class_sz')
                ax.plot(k_class_sz,pk_1h_class_sz+pk_2h_class_sz,color='k',ls='-',alpha = 1.,label = 'pk_tot_class_sz')
                ax.plot(k_class_sz,pk_gg_1h_class_sz,color='b',ls='--',alpha = 1.,label = 'pk_gg_1h_class_sz')
                #ax.plot(k_class_sz,pk_2h_class_sz,color='k',ls='-.',alpha = 1.,label = 'pk_2h_class_sz')
                ax.plot(k_class_sz,pk_gg_1h_class_sz+pk_gg_2h_class_sz,color='b',ls='-',alpha = 1.,label = 'pk_gg_tot_class_sz')

            if (args.plot_bk == 'yes'):
                C = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_halo_model_bk_at_z_k_bk1h_bk2h_bk3h.txt')
                k_class_sz = C[:,0]
                bk_1h_class_sz = C[:,1]
                bk_2h_class_sz = C[:,2]
                bk_3h_class_sz = C[:,3]
                print(bk_1h_class_sz)
                print(bk_2h_class_sz)
                print(bk_3h_class_sz)

                ax.plot(k_class_sz,bk_1h_class_sz,color='r',ls='--',alpha = 1.,label = 'bk_1h_class_sz')
                ax.plot(k_class_sz,bk_2h_class_sz,color='b',ls='-.',alpha = 1.,label = 'bk_2h_class_sz')
                ax.plot(k_class_sz,bk_3h_class_sz,color='g',ls=':',alpha = 1.,label = 'bk_3h_class_sz')
                ax.plot(k_class_sz,bk_1h_class_sz+bk_2h_class_sz+bk_3h_class_sz,color='k',ls='-',alpha = 1.,label = 'bk_tot_class_sz')

            if (args.plot_bk_ttg== 'yes'):
                C = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_halo_model_bk_ttg_at_z_k_bk1h_bk2h_bk3h.txt')
                k_class_sz = C[:,0]
                bk_ttg_1h_class_sz = C[:,1]
                bk_ttg_2h_class_sz = C[:,2]
                bk_ttg_3h_class_sz = C[:,3]
                print(bk_ttg_1h_class_sz)
                print(bk_ttg_2h_class_sz)
                print(bk_ttg_3h_class_sz)

                ax.plot(k_class_sz,bk_ttg_1h_class_sz,color='r',ls='--',alpha = 1.,label = 'bk_ttg_1h_class_sz')
                ax.plot(k_class_sz,bk_ttg_2h_class_sz,color='b',ls='-.',alpha = 1.,label = 'bk_ttg_2h_class_sz')
                ax.plot(k_class_sz,bk_ttg_3h_class_sz,color='g',ls=':',alpha = 1.,label = 'bk_ttg_3h_class_sz')
                ax.plot(k_class_sz,bk_ttg_1h_class_sz+bk_ttg_2h_class_sz+bk_ttg_3h_class_sz,color='k',ls='-',alpha = 1.,label = 'bk_ttg_tot_class_sz')


            elif ('tSZ_1h' in p_dict['output']
            or 'kSZ_kSZ_gal_1h' in p_dict['output']
            or 'kSZ_kSZ_gal fft (1h)' in p_dict['output']
            or 'kSZ_kSZ_gal fft (2h)' in p_dict['output']
            or 'kSZ_kSZ_gal fft (3h)' in p_dict['output']
            or 'kSZ_kSZ_gal_2h' in p_dict['output']
            or 'kSZ_kSZ_gal_3h' in p_dict['output']
            or 'kSZ_kSZ_gal_hf' in p_dict['output']
            or 'kSZ_kSZ_tSZ_1h' in p_dict['output']
            or 'kSZ_kSZ_tSZ_2h' in p_dict['output']
            or 'kSZ_kSZ_tSZ_3h' in p_dict['output']
            or 'kSZ_kSZ_lensmag_1h' in p_dict['output']
            or 'tSZ_lens_1h' in p_dict['output']
            or 'tSZ_lens_2h' in p_dict['output']
            or 'tSZ_gal' in p_dict['output']
            or 'lens_lens' in p_dict['output']
            or 'isw_lens' in p_dict['output']
            or 'gal_gal' in p_dict['output']
            or 'gal_lens' in p_dict['output']
            or 'gal_gallens' in p_dict['output']
            or 'lens_lensmag' in p_dict['output']
            or 'lensmag_lensmag' in p_dict['output']
            or 'isw_tsz' in p_dict['output']
            or 'isw_auto' in p_dict['output']):
            # or 'pk_at_z_1h' in p_dict['output']
            # or 'pk_at_z_2h' in p_dict['output']):
                if args.mode == 'run':
                    R = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum.txt')

                elif args.mode == 'plot':
                    R = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_szpowerspectrum_' + str(p_val) + '.txt')



                multipoles.append(R[:,0])
                cl_1h.append(R[:,1])
                y_err.append(R[:,5]/np.sqrt(f_sky))
                cl_2h.append(R[:,6])
                #print(cl_2h)
                trispectrum.append(R[:,3])
                te_y_y.append(R[:,7])
                kSZ_kSZ_gal_1h.append(R[:,8])
                kSZ_kSZ_gal_1h_fft.append(R[:,54])
                kSZ_kSZ_gal_2h_fft.append(R[:,55])
                kSZ_kSZ_gal_3h_fft.append(R[:,56])
                tSZ_lens_1h.append(R[:,9])
                tSZ_lens_2h.append(R[:,10])
                isw_lens.append(R[:,11])
                isw_tsz.append(R[:,12])
                isw_auto.append(R[:,13])
                tSZ_gal_1h.append(R[:,20])
                tSZ_tSZ_tSZ_1h.append(R[:,21])
                gal_gal_1h.append(R[:,22])
                gal_gal_2h.append(R[:,23])
                gal_lens_1h.append(R[:,24])
                gal_lens_2h.append(R[:,25])
                tSZ_gal_2h.append(R[:,26])
                lens_lens_1h.append(R[:,27])
                lens_lens_2h.append(R[:,28])                # L = [multipole,cl_1h]
                kSZ_kSZ_lensmag_1h.append(R[:,33])
                gal_lensmag_1h.append(R[:,36])
                gal_lensmag_2h.append(R[:,37])
                lensmag_lensmag_1h.append(R[:,38])
                lensmag_lensmag_2h.append(R[:,39])
                lens_lensmag_1h.append(R[:,40])
                lens_lensmag_2h.append(R[:,41])
                kSZ_kSZ_gal_2h.append(R[:,42])
                kSZ_kSZ_gal_3h.append(R[:,43])
                kSZ_kSZ_gal_hf.append(R[:,46])
                gal_gal_hf.append(R[:,47])
                gal_lens_hf.append(R[:,48])
                gal_lensmag_hf.append(R[:,49])
                lens_lensmag_hf.append(R[:,50])
                lensmag_lensmag_hf.append(R[:,51])
                gal_gallens_1h.append(R[:,57])
                gal_gallens_2h.append(R[:,58])
                kSZ_kSZ_tSZ_1h.append(R[:,59])
                kSZ_kSZ_tSZ_2h.append(R[:,60])
                kSZ_kSZ_tSZ_3h.append(R[:,61])

                # r_dict[p_val] = L

                #store value of Cl at ell=100
                if (args.compute_scaling_with_param == 'yes'):
                    x = np.log(multipoles[id_p])
                    y = np.log(cl_1h[id_p])
                    f = interpolate.interp1d(x, y)
                    cl_1h_100.append(np.exp(f(np.log(100.))))





                #ax.plot(x_axis,y_axis,color=next(colors),ls='-',alpha = 1.,label = val_label)
                verr = y_err[id_p]
                if (args.plot_trispectrum == 'yes'):
                    ax.plot(multipoles[id_p],trispectrum[id_p],color=col[id_p],ls='-.',alpha = 1.,label = val_label[id_p],marker='o')

                elif (args.plot_te_y_y == 'yes'):
                    print(te_y_y[id_p])
                    ax.plot(multipoles[id_p],te_y_y[id_p],color=col[id_p],ls='-.',alpha = 1.,label = val_label[id_p],marker='o')
                elif (args.plot_b_kSZ_kSZ_tSZ == 'yes'):
                    print(multipoles[id_p])
                    print(kSZ_kSZ_tSZ_1h[id_p])
                    print(kSZ_kSZ_tSZ_2h[id_p])
                    print(kSZ_kSZ_tSZ_3h[id_p])
                    fac =  1.#(2.726e6)**2*multipoles[id_p]*(multipoles[id_p]+1.)/2./np.pi
                    ax.plot(multipoles[id_p],kSZ_kSZ_tSZ_1h[id_p]*fac,color='k',
                            ls=':',alpha = 1.,
                            # label = val_label[id_p] + ' (1h)',
                            label = '1-halo')#,
                            # markersize = 3,
                            # marker='o')
                    ax.plot(multipoles[id_p],kSZ_kSZ_tSZ_2h[id_p]*fac,color='r',
                            ls='--',alpha = 1.,
                            # label = val_label[id_p] + ' (1h)',
                            label = '2-halo')#,
                            # markersize = 3,
                            # marker='o')
                    ax.plot(multipoles[id_p],kSZ_kSZ_tSZ_3h[id_p]*fac,color='b',
                            ls='-.',alpha = 1.,
                            # label = val_label[id_p] + ' (1h)',
                            label = '3-halo')#,
                            # markersize = 3,
                            # marker='o')
                    ax.plot(multipoles[id_p],(kSZ_kSZ_tSZ_1h[id_p]+kSZ_kSZ_tSZ_2h[id_p]+kSZ_kSZ_tSZ_3h[id_p])*fac,color='k',
                            ls='-',alpha = 0.12,
                            # label = val_label[id_p] + ' (1h)',
                            label = '1+2+3-halo')#,
                            # markersize = 3,
                            # marker='o')
                elif (args.plot_kSZ_kSZ_gal == 'yes'):
                    print(multipoles[id_p])
                    print(kSZ_kSZ_gal_1h_fft[id_p])
                    print(kSZ_kSZ_gal_2h_fft[id_p])
                    print(kSZ_kSZ_gal_3h_fft[id_p])
                    print(kSZ_kSZ_gal_1h[id_p])
                    print(kSZ_kSZ_gal_2h[id_p])
                    print(kSZ_kSZ_gal_3h[id_p])
                    print(kSZ_kSZ_gal_hf[id_p])
                    fac =  (2.726e6)**2*multipoles[id_p]*(multipoles[id_p]+1.)/2./np.pi
                    if kSZ_kSZ_gal_1h[id_p].all() != 0:
                        ax.plot(multipoles[id_p],kSZ_kSZ_gal_1h[id_p]*fac,color='k',
                                ls='-',alpha = 1.,
                                # label = val_label[id_p] + ' (1h)',
                                label = '1-halo',
                                markersize = 3,
                                marker='o')
                    # ax.plot(kSZ_kSZ_gal_1h_ells_ref,kSZ_kSZ_gal_1h_ref*fac,color='k',
                    #         ls='-',alpha = 1.,
                    #         # label = val_label[id_p] + ' (1h)',
                    #         label = '1-halo',
                    #         markersize = 3,
                    #         marker='o')

                    ax.plot(multipoles[id_p],kSZ_kSZ_gal_1h_fft[id_p]*fac,color='r',
                            ls='-',alpha = 1.,
                            # label = val_label[id_p] + ' (1h)',
                            label = '1-halo (FFTs)',
                            markersize = 3,
                            marker='o')

                    # ax.plot(kSZ_kSZ_gal_2h_ells_ref,kSZ_kSZ_gal_2h_ref*fac,color='green',
                    #         ls='-',alpha = 1.,
                    #         # label = val_label[id_p] + ' (1h)',
                    #         label = '2-halo',
                    #         markersize = 3,
                    #         marker='o')
                    ax.plot(multipoles[id_p],kSZ_kSZ_gal_2h_fft[id_p]*fac,color='b',
                            ls='-',alpha = 1.,
                            # label = val_label[id_p] + ' (1h)',
                            label = '2-halo (FFTs)',
                            markersize = 3,
                            marker='o')

                    ax.plot(multipoles[id_p],np.abs(kSZ_kSZ_gal_3h_fft[id_p])*fac,color='green',
                            ls='-',alpha = 1.,
                            # label = val_label[id_p] + ' (1h)',
                            label = '3-halo (FFTs)',
                            markersize = 3,
                            marker='o')
                    if kSZ_kSZ_gal_2h[id_p].all() != 0:
                        ax.plot(multipoles[id_p],kSZ_kSZ_gal_2h[id_p]*fac,color='b',
                                ls='--',alpha = 1.,
                                # label = val_label[id_p] + ' (2h)',
                                label = '2-halo',
                                markersize = 3,
                                marker='o')
                    if kSZ_kSZ_gal_3h[id_p].all() != 0:
                        ax.plot(multipoles[id_p],np.abs(kSZ_kSZ_gal_3h[id_p]*fac),color='orange',
                                ls='-.',alpha = 1.,
                                # label = val_label[id_p] + ' (3h)',
                                label = '3-halo',
                                markersize = 3,
                            marker='o')
                        # ax.plot(multipoles[id_p],-kSZ_kSZ_gal_3h[id_p]*fac,color='orange',
                        #         ls='--',alpha = 1.,
                        #         # label = val_label[id_p] + ' (3h)',
                        #         # label = '3-halo',
                        #         markersize = 3,
                        #         marker='o')
                    ax.plot(multipoles[id_p],np.abs(kSZ_kSZ_gal_1h_fft[id_p]*fac+kSZ_kSZ_gal_2h_fft[id_p]*fac+kSZ_kSZ_gal_3h_fft[id_p]*fac),
                            color='grey',
                            ls='-.',alpha = 1.,
                            # label = val_label[id_p] + ' (3h)',
                            label = '1+2+3-halo')#,
                            # markersize = 3,
                            # marker='o')
                    if kSZ_kSZ_gal_hf[id_p].all() != 0:
                        bgeff = 1.
                        ax.plot(multipoles[id_p],kSZ_kSZ_gal_hf[id_p]*fac*bgeff,
                                color='pink',
                                ls='-',alpha = 1.,
                                # label = val_label[id_p] + ' (3h)',
                                label = 'new computation - effective approach (class_sz)',
                                markersize = 2,
                                markerfacecolor = 'pink',
                                markeredgecolor = 'k',
                                marker='o')
                    # # total = (kSZ_kSZ_gal_1h[id_p]+kSZ_kSZ_gal_2h[id_p]+kSZ_kSZ_gal_3h[id_p])*fac
                    # if (kSZ_kSZ_gal_1h[id_p]+kSZ_kSZ_gal_2h[id_p]+np.nan_to_num(kSZ_kSZ_gal_3h[id_p])).all() != 0:
                    #     ax.plot(multipoles[id_p],(kSZ_kSZ_gal_1h[id_p]+kSZ_kSZ_gal_2h[id_p]+np.nan_to_num(kSZ_kSZ_gal_3h[id_p]))*fac,color='red',
                    #             ls='-',alpha = 1.,
                    #             # label = val_label[id_p] + ' (3h)',
                    #             label = '1+2+3-halo',
                    #             markersize = 2,
                    #             marker='o')
                    # # np.savetxt("/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/ksz2xunwise_hm_data/data_hm_kSZ2g_class_sz_in_muk2_l_dllin.txt",
                    # #            np.c_[multipoles[id_p],kSZ_kSZ_gal_1h[id_p]*fac])
                    # #fac =  (2.726e6)**2*ell_ref*(ell_ref+1.)/2./np.pi
                    # #ax.plot(ell_ref,kSZ_kSZ_gal_1h_ref*fac,c='k',ls='--',label='ref',alpha=0.7)
                    # if id_p == N-1:
                    #     fac = bgeff
                    #     # ax.plot(ell_ref_sf_july,kSZ_kSZ_gal_1h_ref_sf_july*fac,c='r',ls='--',marker='o',label='Simone x bg_eff (July)',alpha=0.7)
                    #     ax.plot(ell_ref_sf_august,kSZ_kSZ_gal_1h_ref_sf_august*fac,
                    #     c='k',
                    #     ls=':',
                    #     marker='o',
                    #     markersize = 2,
                    #     markerfacecolor = 'r',
                    #     label='previous computation (1605.02722)',
                    #     alpha=0.7)
                    #ax.plot(multipoles[id_p],-kSZ_kSZ_gal_1h[id_p],color=col[id_p],ls='--',alpha = 1.,marker='o')
                        # np.savetxt("/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/ksz2xunwise_hm_data/data_hm_bg_times_kSZ2g_sf_in_muk2_l_dllin_dlnonlin.txt",
                        #            np.c_[ell_ref_sf_july,kSZ_kSZ_gal_1h_ref_sf_july*fac,kSZ_kSZ_gal_1h_ref_sf_august*fac])

                elif (args.plot_kSZ_kSZ_lensmag_1h == 'yes'):
                    print(kSZ_kSZ_lensmag_1h[id_p])
                    fac =  (2.726e6)**2*multipoles[id_p]*(multipoles[id_p]+1.)/2./np.pi
                    ax.plot(multipoles[id_p],kSZ_kSZ_lensmag_1h[id_p]*fac,color=col[id_p],ls='-',alpha = 1.,label = 'class_sz hod 1-halo',marker='o')
                    #fac =  (2.726e6)**2*ell_ref*(ell_ref+1.)/2./np.pi
                    #ax.plot(ell_ref,kSZ_kSZ_gal_1h_ref*fac,c='k',ls='--',label='ref',alpha=0.7)
                    np.savetxt(path_to_class_external_data+"/ksz2xunwise_hm_data/data_hm_kSZ2mu_class_sz_in_muk2_l_dlnonlin_cvirtau_2.txt",
                               np.c_[multipoles[id_p],kSZ_kSZ_lensmag_1h[id_p]*fac])
                    fac = 1.
                    ax.plot(ell_ref_lensmag_sf_july,kSZ_kSZ_lensmag_1h_ref_sf_july*fac,c='r',ls='--',label='Simone  (July)',alpha=0.7)
                    ax.plot(ell_ref_lensmag_sf_august,kSZ_kSZ_lensmag_1h_ref_sf_august*fac,c='r',ls=':',label='Simone (August)',alpha=0.7)

                    np.savetxt(path_to_class_external_data+"/ksz2xunwise_hm_data/data_hm_kSZ2mu_sf_in_muk2_l_dllin_dlnonlin.txt",
                               np.c_[ell_ref_lensmag_sf_july,kSZ_kSZ_lensmag_1h_ref_sf_july*fac,kSZ_kSZ_lensmag_1h_ref_sf_august*fac])


                elif (args.plot_tSZ_tSZ_tSZ_1h == 'yes'):
                    print(tSZ_tSZ_tSZ_1h[id_p])
                    ax.plot(multipoles[id_p],tSZ_tSZ_tSZ_1h[id_p],color=col[id_p],ls='-',alpha = 1.,label = val_label[id_p],marker='o')
                    ax.plot(multipoles[id_p],-tSZ_tSZ_tSZ_1h[id_p],color=col[id_p],ls='--',alpha = 1.,marker='o')
                elif (args.plot_tSZ_lens == 'yes'):
                    print(tSZ_lens_1h[id_p])
                    print(tSZ_lens_2h[id_p])
                    #for (nu,colg) in zip((100,143,353),('k','r','b')):
                    for (nu,colg) in zip((100,353),('k','r')):
                        g = g_nu(nu)
                        print(g)
                        if g>0:
                            ax.plot(multipoles[id_p],g*tSZ_lens_1h[id_p]*Tcmb,color=colg,ls='-',alpha = 1.,label = '1-halo @ %.2f GHz'%nu)
                            ax.plot(multipoles[id_p],g*tSZ_lens_2h[id_p]*Tcmb,color='g',ls='-',alpha = 1.,label = '2-halo')
                            ax.plot(multipoles[id_p],g*(tSZ_lens_2h[id_p]+tSZ_lens_1h[id_p])*Tcmb,color=colg,ls=':',alpha = 1.,label = val_label[id_p])
                        else:
                            ax.plot(multipoles[id_p],-g*tSZ_lens_1h[id_p]*Tcmb,color=colg,ls='--',alpha = 1.,label='1-halo @ %.2f GHz'%nu)
                            ax.plot(multipoles[id_p],-g*tSZ_lens_2h[id_p]*Tcmb,color='g',ls='--',alpha = 1.,label='2-halo')
                            ax.plot(multipoles[id_p],-g*(tSZ_lens_2h[id_p]+tSZ_lens_1h[id_p])*Tcmb,color=colg,ls=':',alpha = 1.,label = val_label[id_p])
                elif (args.plot_tSZ_gal == 'yes'):
                    print(tSZ_gal_1h[id_p])
                    print(tSZ_gal_2h[id_p])
                    #for (nu,colg) in zip((100,143,353),('k','r','b')):
                    # b_ell = exp(-0.5 * ell * (ell + 1) * sigma^2) where sigma = fwhm / sqrt(8 * log(2))
                    beam_fwhm = 0.00290888 # 10' in radians
                    beam_sigma = beam_fwhm / np.sqrt(8 * np.log(2))
                    beam_factor = 1.#np.exp(-0.5 * multipoles[id_p] * (multipoles[id_p] + 1.) * beam_sigma**2.)
                    fac = 1e6/1e9/beam_factor
                    if (val_label[id_p] == 'gr_shallow'):
                        ax.plot(multipoles[id_p],(tSZ_gal_1h[id_p])/fac,color='forestgreen',ls='--',alpha = 1.,label = val_label[id_p])
                    elif (val_label[id_p] == 'green'):
                        ax.plot(multipoles[id_p],(tSZ_gal_1h[id_p])/fac,color='green',ls='-',alpha = 1.,label = val_label[id_p])
                    elif (val_label[id_p] == 'blue'):
                        ax.plot(multipoles[id_p],(tSZ_gal_1h[id_p])/fac,color='blue',ls='-',alpha = 1.,label = val_label[id_p])
                    elif (val_label[id_p] == 'red'):
                        ax.plot(multipoles[id_p],(tSZ_gal_1h[id_p])/fac,color='red',ls='-',alpha = 1.,label = val_label[id_p])
                    else:
                        ax.plot(multipoles[id_p],(tSZ_gal_1h[id_p])/fac,ls='--',alpha = 1.,label = 'class_sz hod: 1-halo '+val_label[id_p])
                        ax.plot(multipoles[id_p],(tSZ_gal_2h[id_p])/fac,ls='-',alpha = 1.,label = 'class_sz hod: 2-halo '+val_label[id_p])
                    #ax.plot(ell_KA20,cl_yg_1h_KA20,label='KA20')
                elif (args.plot_gal_gallens == 'yes'):
                    print(gal_gallens_1h[id_p])
                    print(gal_gallens_2h[id_p])
                    fac = multipoles[id_p]*(multipoles[id_p]+1.)/2./np.pi/1e5
                    if ('gal_gallens_1h' in p_dict['output']):
                        ax.plot(multipoles[id_p],(gal_gallens_1h[id_p])/fac,color=col[id_p],ls=':',alpha = 1.,label =  r'g$\gamma$ (1h)')
                    if ('gal_gallens_2h' in p_dict['output']):
                        ax.plot(multipoles[id_p],(gal_gallens_2h[id_p])/fac,color=col[id_p],ls='--',alpha = 1.,label =  r'g$\gamma$ (2h)')
                    if (('gal_gallens_1h' or  'gal_gallens_2h') in p_dict['output']):
                        ax.plot(multipoles[id_p],(gal_gallens_1h[id_p])/fac+(gal_gallens_2h[id_p])/fac,color=col[id_p],ls='-.',alpha = 1.,label =  r'g$\gamma$ (1+2h)')

                elif (args.plot_gal_gal == 'yes'):
                    print('plotting gxg')
                    fac = multipoles[id_p]*(multipoles[id_p]+1.)/2./np.pi/1e5
                    print(gal_gal_1h[id_p]/(fac*1e5))
                    print(gal_gal_2h[id_p])
                    print(gal_gal_hf[id_p])
                    print('gal_lensmag')
                    print(gal_lensmag_1h[id_p])
                    print(gal_lensmag_2h[id_p])
                    print(gal_lensmag_hf[id_p])
                    print('lensmag_lensmag')
                    print(lensmag_lensmag_1h[id_p])
                    print(lensmag_lensmag_2h[id_p])
                    print(lensmag_lensmag_hf[id_p])
                    # #ax.plot(multipoles[id_p],(gal_gal_1h[id_p])/fac,color=col[id_p],ls='-',alpha = 1.,label = 'class_sz gg-1h')
                    # #ax.plot(ell_KA20,cl_gg_1h_KA20,label='KA20-1h')
                    # #ax.plot(multipoles[id_p],(gal_gal_2h[id_p])/fac,color=col[id_p],ls='--',alpha = 1.,label = 'class_sz gg-2h')
                    # # for shot_noise see Table 1 of KFSW20
                    # if (val_label[id_p] == 'gr_shallow'):
                    #     ax.plot(multipoles[id_p],(gal_gal_2h[id_p])/fac,color='forestgreen',ls='--',alpha = 1.,label = val_label[id_p])
                    # elif (val_label[id_p] == 'green'):
                    #     shot_noise = (1846/3.046174198e-4)**-1*1e5
                    #     ax.plot(multipoles[id_p],(gal_gal_2h[id_p])/fac+shot_noise,color='green',ls='-',alpha = 1.,label = val_label[id_p])
                    #     ax.plot(multipoles[id_p],shot_noise*multipoles[id_p]/multipoles[id_p], c='green',ls='--')
                    # elif (val_label[id_p] == 'blue'):
                    #     shot_noise = (3409/3.046174198e-4)**-1*1e5
                    #     ax.plot(multipoles[id_p],(gal_gal_2h[id_p])/fac+shot_noise,color='blue',ls='-',alpha = 1.,label = val_label[id_p])
                    #     ax.plot(multipoles[id_p],shot_noise*multipoles[id_p]/multipoles[id_p], c='blue',ls='--')
                    # elif (val_label[id_p] == 'red'):
                    #     shot_noise = (144/3.046174198e-4)**-1*1e5
                    #     ax.plot(multipoles[id_p],(gal_gal_2h[id_p])/fac+shot_noise,color='red',ls='-',alpha = 1.,label = val_label[id_p])
                    #     ax.plot(multipoles[id_p],shot_noise*multipoles[id_p]/multipoles[id_p], c='red',ls='--')
                    #     #ax.plot(multipoles[id_p],(gal_gal_1h[id_p])/fac,color='red',ls='--',alpha = 1.,label = val_label[id_p])
                    # else:
                    #     #ax.plot(multipoles[id_p],(gal_gal_2h[id_p])/fac,color=col[id_p],ls='-',alpha = 1.,label = 'class_sz')
                    #
                    #     #couleur = 'forestgreen'
                    #     if (val_label[id_p] == 'yes'):
                    #         ax.plot(multipoles[id_p],(gal_gal_2h[id_p])/fac,color=couleur,ls='-',alpha = 1.,
                    #         marker =  'o',markersize = 5,markerfacecolor='None',
                    #         label = 'class_sz hod (2-halo)')
                    #         # ax.plot(multipoles[id_p],(gal_gal_1h[id_p])/fac,color=couleur,ls='--',alpha = 1.,
                    #         #  marker =  'o',markersize = 5,markerfacecolor='None',
                    #         #  label = 'simplified hod: ' + val_label[id_p] + ' (1-halo)')
                    #         #ax.plot(multipoles[id_p],(gal_gal_1h[id_p])/fac+(gal_gal_2h[id_p])/fac,color=couleur,ls=':',alpha = 1.,
                    #         # marker =  'o',markersize = 5,markerfacecolor='None',
                    #         # label = 'simplified hod: ' + val_label[id_p] + ' (1+2-halo)')
                    #     else:
                    #ax.plot(multipoles[id_p],(gal_gal_1h[id_p])/fac,color=col[id_p],ls='-.',alpha = 0.7)
                    #ax.plot(multipoles[id_p],(gal_gal_2h[id_p])/fac,color=col[id_p],ls='-',alpha = 0.4)

                    # if (p_dict['use_hod'] == "yes"):
                    # halofit = np.loadtxt(path_to_class_external_data+"/ksz2xunwise_hm_data/data_hm_1e5cl_gg_halofit_l_gg_gmu_mumu_%s.txt"%(couleur))
                    # ax.plot(halofit[:,0],halofit[:,1],color='k',ls='-',alpha = 1.,label = r'halofit $gg$')
                    # ax.plot(halofit[:,0],halofit[:,2],color='green',ls='--',alpha = 1.,label = r'halofit $g\mu$')
                    # ax.plot(halofit[:,0],halofit[:,3],color='orange',ls='-.',alpha = 1.,label = r'halofit $\mu\mu$')
                    if ('gal_gal_1h' in p_dict['output']):
                        ax.plot(multipoles[id_p],(gal_gal_1h[id_p])/fac,color=col[id_p],ls=':',alpha = 1.,label =  'gg (1h)')
                    if ('gal_gal_2h' in p_dict['output']):
                        ax.plot(multipoles[id_p],(gal_gal_2h[id_p])/fac,color=col[id_p],ls='--',alpha = 1.,label =  'gg (2h)')
                    if (('gal_gal_1h' or  'gal_gal_2h') in p_dict['output']):
                        ax.plot(multipoles[id_p],(gal_gal_1h[id_p])/fac+(gal_gal_2h[id_p])/fac,color=col[id_p],ls='-.',alpha = 1.,label =  'gg (1+2h)')
                    # np.savetxt(path_to_class_external_data+"/ksz2xunwise_hm_data/data_hm_1e5cl_gg_gamma_%s_l_1h_1p2h.txt"%(str(p_val)),
                    #            np.c_[multipoles[id_p],(gal_gal_1h[id_p])/fac,(gal_gal_1h[id_p]+gal_gal_2h[id_p])/fac])
                # elif (p_dict['use_hod'] == "no"):
                    if ('gal_gal_hf' in p_dict['output']):
                        ax.plot(multipoles[id_p],(gal_gal_hf[id_p])/fac,color='r',ls='-.',alpha = 1.,
                        #marker =  'o',markersize = 3,
                        label = r'$gg (hf)$')
                    # ax.plot(multipoles[id_p],(gal_gal_1h[id_p]+gal_gal_2h[id_p])/fac,color='k',ls='-',alpha = 1.,
                    # #marker =  'o',markersize = 3,
                    # label = r'$gg$')
                    if (('gal_lensmag_1h' or 'gal_lensmag_2h') in p_dict['output']):
                        ax.plot(multipoles[id_p],(5.*smag-2.)*(gal_lensmag_1h[id_p]+gal_lensmag_2h[id_p])/fac,color='green',ls='--',alpha = 1.,
                        #marker =  'o',markersize = 3,
                        label = r'$g\mu$')
                    if ('gal_lensmag_hf' in p_dict['output']):
                        ax.plot(multipoles[id_p],(5.*smag-2.)*(gal_lensmag_hf[id_p])/fac,color='green',ls='--',alpha = 1.,
                        #marker =  'o',markersize = 3,
                        label = r'$g\mu$ (hf)')

                    if (('lensmag_lensmag_1h' or 'lensmag_lensmag_2h') in p_dict['output']):
                        ax.plot(multipoles[id_p],(5.*smag-2.)**2.*(lensmag_lensmag_1h[id_p]+lensmag_lensmag_2h[id_p])/fac,color='orange',ls='-.',alpha = 1.,
                        #marker =  'o',markersize = 3,
                        label = r'$\mu\mu$')

                    if ('lensmag_lensmag_hf' in p_dict['output']):
                        ax.plot(multipoles[id_p],(5.*smag-2.)**2.*(lensmag_lensmag_hf[id_p])/fac,color='orange',ls='-.',alpha = 1.,
                        #marker =  'o',markersize = 3,
                        label = r'$\mu\mu$ (hf)')
                    # ax.plot(multipoles[id_p],(gal_gal_1h[id_p]+gal_gal_2h[id_p])/fac+(5.*smag-2.)*(gal_lensmag_1h[id_p]+gal_lensmag_2h[id_p])/fac,color='red',ls='-',alpha = 1.,
                    # #marker =  'o',markersize = 3,
                    # label = r'$gg+g\mu$')
                    cl_tot = (gal_gal_1h[id_p]+gal_gal_2h[id_p])/fac\
                                             +(5.*smag-2.)*(gal_lensmag_1h[id_p]+gal_lensmag_2h[id_p])/fac\
                                             +(5.*smag-2.)**2.*(lensmag_lensmag_1h[id_p]+lensmag_lensmag_2h[id_p])/fac
                    if not np.all(cl_tot==0):
                        ax.plot(multipoles[id_p],cl_tot,color='orange',ls='-',alpha = 1.,
                        #marker =  'o',markersize = 3,
                        label = r'$gg+g\mu+\mu\mu$')
                    cl_tot = (gal_gal_hf[id_p])/fac \
                             +(5.*smag-2.)*(gal_lensmag_hf[id_p])/fac \
                                             +(5.*smag-2.)**2.*(lensmag_lensmag_hf[id_p])/fac
                    if not np.all(cl_tot==0):
                        ax.plot(multipoles[id_p],cl_tot,color='red',ls='-',alpha = 1.,
                        #marker =  'o',markersize = 3,
                        label = r'$gg+g\mu+\mu\mu$ (hf)')
                        # np.savetxt(path_to_class_external_data+"/ksz2xunwise_hm_data/data_hm_1e5cl_gg_halofit_l_gg_gmu_mumu_%s.txt"%(couleur),
                        #            np.c_[multipoles[id_p],(gal_gal_1h[id_p]+gal_gal_2h[id_p])/fac,(5.*smag-2.)*(gal_lensmag_1h[id_p]+gal_lensmag_2h[id_p])/fac,(5.*smag-2.)**2.*(lensmag_lensmag_1h[id_p]+lensmag_lensmag_2h[id_p])/fac])
                    #shot_noise = (144/3.046174198e-4)**-1*1e5 # red
                    #shot_noise = (3409/3.046174198e-4)**-1*1e5 # blue
                    if id_p == N-1:
                        #shot_noise = (1846/3.046174198e-4)**-1*1e5 # green
                        #shot_noise = 0.296#(144/3.046174198e-4)**-1*1e5 # red
                        #shot_noise = (144/3.046174198e-4)**-1*1e5 # red
                        # #ax.plot(multipoles[id_p],(multipoles[id_p])/(multipoles[id_p])*shot_noise,color='k',ls='--',alpha = .5,
                        # #label = 'shot noise')
                        # halofit_approach = np.loadtxt(path_to_class_external_data+'/class_sz_KFSW_green_halofit/class-sz_szpowerspectrum_no.txt')
                        # halofit_approach_multipoles = halofit_approach[:,0]
                        # halofit_approach_gal_lens_2h = halofit_approach[:,23]
                        # print('halofit')
                        # print(halofit_approach_gal_lens_2h)
                        #
                        #
                        #
                        # fac = halofit_approach_multipoles*(halofit_approach_multipoles+1.)/2./np.pi/1e5
                        # ax.plot(halofit_approach_multipoles,halofit_approach_gal_lens_2h/fac,ls='--',c='r',label = 'Halofit')
                        #
                        # np.savetxt("/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/ksz2xunwise_hm_data/data_hm_1e5cl_gg_halofit_l_1e5cl.txt",
                        #            np.c_[halofit_approach_multipoles,halofit_approach_gal_lens_2h/fac])
                        if couleur == 'red':
                            alex = np.loadtxt(path_to_class_external_data+'/alex_measurements/Bandpowers_Auto_Sample7.dat')
                        if couleur == 'blue':
                            alex = np.loadtxt(path_to_class_external_data+'/alex_measurements/Bandpowers_Auto_Sample1.dat')
                        if couleur == 'green':
                            alex = np.loadtxt(path_to_class_external_data+'/alex_measurements/Bandpowers_Auto_Sample2.dat')
                        ax.errorbar(alex[0,:],1e5*alex[1,:]-shot_noise,yerr=[1e5*alex[2,:],1e5*alex[2,:]],ls='none',
                                    marker='o',markersize=3,capsize=5,c='k',label='KFSW20')
                        print(len(alex[0,:]))
                        print(len(alex[1,:]))
                        print(len(alex[2,:]))
                        # np.savetxt(path_to_class_external_data+"/ksz2xunwise_hm_data/data_hm_1e5cl_gg_kfsw20_l_1e5cl_1e5errcl.txt",
                        #            np.c_[alex[0,:],1e5*alex[1,:]-shot_noise,1e5*alex[2,:]])

                        #print(alex[:,0])
                    # ax.plot(multipoles[id_p],(gal_gal_1h[id_p])/fac,color=couleur,ls='--',alpha = 1.,
                    # marker =  '*',markersize = 7,
                    # label = 'simplified hod: ' + val_label[id_p] + ' (1-halo)')

                elif (args.plot_gal_lens == 'yes'):
                    print('plotting gxkappa')

                    fac = multipoles[id_p]*(multipoles[id_p]+1.)/2./np.pi/1e5
                    #if id_p == 0 :
                    #    fac /= 2
                    print('gal_lens')
                    print(gal_lens_1h[id_p])
                    print(gal_lens_2h[id_p])
                    print(gal_lens_hf[id_p])
                    print('lens_lensmag')
                    print(lens_lensmag_1h[id_p])
                    print(lens_lensmag_2h[id_p])
                    print(lens_lensmag_hf[id_p])
                    # if (p_dict['use_hod'] == "yes"):
                    # halofit = np.loadtxt(path_to_class_external_data+"/ksz2xunwise_hm_data/data_hm_1e5cl_gkappa_halofit_l_gk_muk_%s.txt"%(couleur))
                    # ax.plot(halofit[:,0],halofit[:,1],color='k',ls='-',alpha = 1.,label = r'halofit $\kappa g$')
                    # ax.plot(halofit[:,0],halofit[:,2],color='green',ls='--',alpha = 1.,label = r'halofit $\kappa \mu$')
                    if ('gal_lens_1h' in p_dict['output']):
                        ax.plot(multipoles[id_p],(gal_lens_1h[id_p])/fac,color=col[id_p],ls=':',alpha = 1.,label =  'gk (1h)')
                    if ('gal_lens_2h' in p_dict['output']):
                        ax.plot(multipoles[id_p],(gal_lens_2h[id_p])/fac,color=col[id_p],ls='--',alpha = 1.,label =  'gk (2h)')
                    if (('gal_lens_1h' or 'gal_lens_2h') in p_dict['output']):
                        ax.plot(multipoles[id_p],(gal_lens_1h[id_p])/fac+(gal_lens_2h[id_p])/fac,color=col[id_p],ls='-.',alpha = 1.,label =  'gk (1+2h)')
                    if ('gal_lens_hf' in p_dict['output']):
                        ax.plot(multipoles[id_p],(gal_lens_hf[id_p])/fac,color=col[id_p],ls='--',alpha = 1.,label =  'gk (hf)')

                    # np.savetxt("/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/ksz2xunwise_hm_data/data_hm_1e5cl_gkappa_gamma_%s_l_1h_1p2h.txt"%(str(p_val)),
                    #            np.c_[multipoles[id_p],(gal_lens_1h[id_p])/fac,(gal_lens_1h[id_p]+gal_lens_2h[id_p])/fac])

                #ax.plot(multipoles[id_p],(gal_lens_2h[id_p])/fac,color=col[id_p],ls='-',alpha = 0.2)
                #ax.plot(multipoles[id_p],(gal_lens_1h[id_p])/fac,color=col[id_p],ls='-.',alpha = 1.,label='1H')
                # elif (p_dict['use_hod'] == "no"):
                    # ax.plot(multipoles[id_p],(gal_lens_1h[id_p]+gal_lens_2h[id_p])/fac,color='k',ls='-',alpha = 0.5,
                    # #marker =  'o',markersize = 3,
                    # label = r'$\kappa g$')
                    if ('lens_lensmag_1h' in p_dict['output']):
                        ax.plot(multipoles[id_p],(5.*smag-2.)*(lens_lensmag_1h[id_p]+lens_lensmag_2h[id_p])/fac,color='green',ls='--',alpha = 1.,
                        #marker =  'o',markersize = 3,
                        label = r'$\kappa \mu$')

                    cl_tot = (gal_lens_1h[id_p]+gal_lens_2h[id_p])/fac+(5.*smag-2.)*(lens_lensmag_1h[id_p]+lens_lensmag_2h[id_p])/fac
                    if not np.all(cl_tot == 0):
                        ax.plot(multipoles[id_p],cl_tot,color='red',ls='-',alpha = 1.,
                        #marker =  'o',markersize = 3,
                        label = r'$\kappa g+\kappa\mu$')

                    cl_tot = (5.*smag-2.)*(lens_lensmag_hf[id_p])/fac
                    if not np.all(cl_tot == 0):
                        ax.plot(multipoles[id_p],cl_tot,color='green',ls='--',alpha = 1.,
                        #marker =  'o',markersize = 3,
                        label = r'$\kappa \mu$ (hf)')
                    cl_tot = (gal_lens_hf[id_p])/fac+(5.*smag-2.)*(gal_lensmag_hf[id_p])/fac
                    if not np.all(cl_tot == 0):
                        ax.plot(multipoles[id_p],cl_tot,color='red',ls='-',alpha = 1.,
                        #marker =  'o',markersize = 3,
                        label = r'$\kappa g+\kappa\mu$ (hf)')

                    # np.savetxt("/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/ksz2xunwise_hm_data/data_hm_1e5cl_gkappa_halofit_l_gk_muk_%s.txt"%(couleur),
                    #            np.c_[multipoles[id_p],(gal_lens_1h[id_p]+gal_lens_2h[id_p])/fac,(5.*smag-2.)*(lens_lensmag_1h[id_p]+lens_lensmag_2h[id_p])/fac])

                    ## plot data:
                    # if id_p == N-1:
                    #
                    #     if couleur == 'red':
                    #         #alex = np.loadtxt(path_to_class_external_data+'/alex_measurements/Bandpowers_Auto_Sample7.dat')
                    #         alex = np.loadtxt(path_to_class_external_data+'/alex_measurements/Bandpowers_Cross_Sample7_Sample5.dat')
                    #         OK = np.loadtxt(path_to_class_external_data+'/FINALred_c_ell_kappa_g.txt')
                    #         alex_kg = np.loadtxt(path_to_class_external_data+'/theory_kg_plus_kmu_low_ind_1.txt')
                    #         alex_kg = alex_kg[:]
                    #         alex_kg_ell = np.arange(0,len(alex_kg))
                    #         ax.plot(alex_kg_ell,1e5*alex_kg_ell*alex_kg,label='alex kg+km',ls=':',c='r')
                    #         alex_km = np.loadtxt(path_to_class_external_data+'/theory_kmu_low_ind_1.txt')
                    #         alex_km = alex_km[:]
                    #         alex_km_ell = np.arange(0,len(alex_km))
                    #         #ax.plot(alex_km_ell,1e5*alex_km_ell*alex_km,label='alex km',ls=':')
                    #     if couleur == 'blue':
                    #         #alex = np.loadtxt(path_to_class_external_data+'/alex_measurements/Bandpowers_Auto_Sample1.dat')
                    #         alex = np.loadtxt(path_to_class_external_data+'/alex_measurements/Bandpowers_Cross_Sample1_Sample5.dat')
                    #         OK = np.loadtxt(path_to_class_external_data+'/FINALblue_c_ell_kappa_g.txt')
                    #     if couleur == 'green':
                    #         alex = np.loadtxt(path_to_class_external_data+'/alex_measurements/Bandpowers_Cross_Sample2_Sample5.dat')
                    #         OK = np.loadtxt(path_to_class_external_data+'/FINALgreen_c_ell_kappa_g.txt')
                    #         #alex = np.loadtxt(path_to_class_external_data+'/alex_measurements/Bandpowers_Auto_Sample2.dat')
                    #
                    #     # ax.errorbar(alex[0,:],1e5*alex[1,:]*alex[0,:],yerr=[1e5*alex[2,:]*alex[0,:],1e5*alex[2,:]*alex[0,:]],ls='none',
                    #     #             marker='o',markersize=3,capsize=5,c='k',label='KFSW20')
                    #     ax.errorbar(alex[0,:],1e5*alex[1,:],yerr=[1e5*alex[2,:],1e5*alex[2,:]],ls='none',
                    #                 marker='o',markersize=3,capsize=5,c='k',label='KFSW20')
                    #
                    #     # np.savetxt("/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz_external_data_and_scripts/ksz2xunwise_hm_data/data_hm_1e5lcl_gkappa_kfsw20_l_1e5lcl_1e5errlcl.txt",
                    #     #            np.c_[alex[0,:],1e5*alex[1,:]*alex[0,:],1e5*alex[2,:]*alex[0,:]])

                elif (args.plot_lens_lens == 'yes'):
                    print(lens_lens_1h[id_p])
                    print(lens_lens_2h[id_p])

                    fac = 4./(multipoles[id_p]*(multipoles[id_p]+1.))
                    if val_label[id_p] == 'no':
                        ax.plot(multipoles[id_p],fac*lens_lens_2h[id_p],color='k',ls='-',alpha = 0.5,label = 'halofit approach')
                    else:
                        ax.plot(multipoles[id_p],fac*lens_lens_1h[id_p],color=col[id_p],ls=':',alpha = 1.,label = '1-h')
                        ax.plot(multipoles[id_p],fac*lens_lens_2h[id_p],color=col[id_p],ls='--',alpha = 1.,label = '2-h')
                        ax.plot(multipoles[id_p],fac*lens_lens_2h[id_p]+fac*lens_lens_1h[id_p],color=col[id_p],ls='-',alpha = 1.,label = '1+2-h')
                    if id_p == 0:
                        try:
                            class_cls = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_cl.dat')
                            fac = (class_cls[:,0]*(class_cls[:,0]+1.))
                            ax.plot(class_cls[:,0],fac*class_cls[:,1],ls=':',c='orange',label = 'Cl^lens-lens class')
                        except OSError:
                            print("ask for lcl if you want to see class's limber phiphi")



            plt.draw()
            plt.pause(0.05)
            if args.mode == 'run':
                #save some results and remove the temporary files
                try:
                    subprocess.call(['mv',path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum.txt', path_to_class+'sz_auxiliary_files/run_scripts/tmp/' + 'class-sz_szpowerspectrum_' + str(p_val) + '.txt'])
                except:
                    print("no szpowerspectrum file")
            id_p += 1
        if (args.plot_lens_lens == 'yes'):
            ax.legend(loc=1,ncol = 1,frameon=True,fontsize=11)
            ax.set_xscale('log')
            ax.set_yscale('log')
        elif (args.plot_kSZ_kSZ_gal == 'yes'):
            # ax.set_xscale('log')
            # ax.set_yscale('log')
            ax.set_xscale('linear')
            ax.set_yscale('linear')
            ax.legend(loc=2,ncol = 1,frameon=True,fontsize=8)
        elif (args.plot_kSZ_kSZ_lensmag_1h == 'yes'):
            ax.set_xscale('log')
            ax.legend(loc=4,ncol = 1,frameon=True,fontsize=11)
        elif (args.plot_gal_lens == 'yes'):
            ax.set_xscale('linear')
            ax.set_yscale('linear')
            ax.legend(loc=1,ncol = 1,frameon=True,fontsize=11)
        elif (args.plot_gal_gal == 'yes'):
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend(loc=3,ncol = 1,frameon=True,fontsize=11)
        elif (args.plot_b_kSZ_kSZ_tSZ == 'yes'):
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend(loc=3,ncol = 1,frameon=True,fontsize=11)
        elif (args.plot_pk == 'yes'):
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend(loc=4,ncol = 1,frameon=True,fontsize=11)
        elif (args.plot_bk == 'yes'):
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_ylim(1e2,1e11)
            ax.set_xlim(1e-3,1e1)
            ax.legend(loc=1,ncol = 3,frameon=True,fontsize=11)
        elif (args.plot_bk_ttg == 'yes'):
            ax.set_xscale('log')
            ax.set_yscale('log')
            # ax.set_ylim(1e2,1e11)
            ax.set_xlim(1e-3,1e1)
            ax.legend(loc=1,ncol = 3,frameon=True,fontsize=11)
        else:
            ax.set_xscale('linear')
            ax.legend(loc=1,ncol = 1,frameon=True,fontsize=11)
        # ax.set_yscale('log')
        # ax.set_xscale('log')

    # ax.set_xscale('log')
    # ax.set_yscale('log')


    #
    # if (args.show_legend == 'yes'):
    #     if (args.plot_isw_tsz == 'yes' or args.plot_isw_auto == 'yes'):
    #         ax1.legend(loc=1)
    #     elif (args.plot_tSZ_lens == 'yes' or args.plot_isw_lens == 'yes'  or args.plot_tSZ_gal == 'yes'):
    #         ax1.legend(loc=3,ncol = 2)
    #     else:
    #         ax1.legend(loc=2)
    #     if (args.print_rel_diff == 'yes'):
    #         ax2.legend(loc=1,ncol = 1)
    #     else:
    #         ax1.legend(loc=2,ncol = 1)
    print("done")

    fig.tight_layout()


    if (args.save_fig == 'yes'):
        if (args.plot_pk == 'yes'):
            print('saving fig')
            FIG_NAME = '/pk'
            plt.savefig(FIG_DIR + FIG_NAME +".pdf")
        if (args.plot_bk == 'yes'):
            FIG_NAME = '/bk'
            plt.savefig(FIG_DIR + FIG_NAME +".pdf")
        if (args.plot_bk_ttg == 'yes'):
            FIG_NAME = '/bk_ttg'
            plt.savefig(FIG_DIR + FIG_NAME +".pdf")
        if (args.plot_gal_lens == 'yes'):
            FIG_NAME = '/cl_gkappa_'
            plt.savefig(FIG_DIR + FIG_NAME +"_"+couleur+".pdf")
        elif (args.plot_gal_gallens == 'yes'):
            FIG_NAME = '/cl_ggamma_'
            plt.savefig(FIG_DIR + FIG_NAME +"_"+".pdf")
        elif (args.plot_gal_gal == 'yes'):
            FIG_NAME = '/cl_gg_'
            plt.savefig(FIG_DIR + FIG_NAME +"_"+couleur+".pdf")
        elif (args.plot_lens_lens == 'yes'):
            FIG_NAME = '/cl_phiphi_'
            plt.savefig(FIG_DIR + FIG_NAME +".pdf")
        elif (args.plot_kSZ_kSZ_gal == 'yes'):
            FIG_NAME = '/cl_ksz2g_'
            plt.savefig(FIG_DIR + FIG_NAME +"_"+couleur+".pdf")
        elif (args.plot_b_kSZ_kSZ_tSZ == 'yes'):
            FIG_NAME = '/b_ksz_kSZ_tSZ_'
            plt.savefig(FIG_DIR + FIG_NAME +"_"+couleur+".pdf")
        elif (args.plot_kSZ_kSZ_lensmag_1h == 'yes'):
            FIG_NAME = '/cl_ksz2mu_'
            plt.savefig(FIG_DIR + FIG_NAME +"_"+couleur+".pdf")
        elif (args.plot_tSZ_gal == 'yes'):
            FIG_NAME = '/cl_yxg_'
            plt.savefig(FIG_DIR + FIG_NAME +"_"+couleur+".pdf")
    plt.show(block=True)






        #subprocess.call(['rm','-r',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])









def main():
    parser=argparse.ArgumentParser(description="Plot cosmotherm spectra")
    parser.add_argument("-mode",help="plot or run" ,dest="mode", type=str, required=True)
    parser.add_argument("-param_name",help="name of varying parameter" ,dest="param_name", type=str, required=True)
    parser.add_argument("-min",help="minimum value of parameter" ,dest="p_min", type=str, required=False)
    parser.add_argument("-max",help="maximum value of parameter" ,dest="p_max", type=str, required=False)
    parser.add_argument("-N",help="number of evaluations" ,dest="N", type=int, required=False)
    parser.add_argument("-p_val",help="list of param values" ,dest="p_val", type=str, required=False)
    parser.add_argument("-spacing",help="linear (lin) or log spacing (log)" ,dest="spacing", type=str, required=False)
    parser.add_argument("-show_legend",help="show legend on figure? ('yes' or 'no')" ,dest="show_legend", type=str, required=True)
    parser.add_argument("-show_error_bars",help="show error bars on figure? ('yes' or 'no')" ,dest="show_error_bars", type=str, required=True)
    parser.add_argument("-y_min",help="ylim for y-axis" ,dest="y_min", type=str, required=False)
    parser.add_argument("-y_max",help="ylim for y-axis" ,dest="y_max", type=str, required=False)
    parser.add_argument("-x_min",help="xlim for x-axis" ,dest="x_min", type=str, required=False)
    parser.add_argument("-x_max",help="xlim for x-axis" ,dest="x_max", type=str, required=False)
    parser.add_argument("-f_sky",help="sky fraction f_sky" ,dest="f_sky", type=str, required=False)
    parser.add_argument("-output",help="what quantities to plot" ,dest="output", type=str, required=False)
    parser.add_argument("-plot_ref_data",help="some other spectra" ,dest="plot_ref_data", type=str, required=False)
    parser.add_argument("-compute_scaling_with_param",help="Compute alpha in C_l ~ p^alpha at l=100" ,dest="compute_scaling_with_param", type=str, required=False)
    parser.add_argument("-save_tsz_ps",help="save file with tsz power spectrum in output directory" ,dest="save_tsz_ps", type=str, required=False)
    parser.add_argument("-save_figure",help="save figure" ,dest="save_fig", type=str, required=False)
    parser.add_argument("-print_rel_diff",help="[cl-cl_ref]/cl_ref" ,dest="print_rel_diff", type=str, required=False)
    parser.add_argument("-plot_trispectrum",help="Tll" ,dest="plot_trispectrum", type=str, required=False)
    parser.add_argument("-plot_pk",help="pk" ,dest="plot_pk", type=str, required=False)
    parser.add_argument("-plot_bk",help="bk" ,dest="plot_bk", type=str, required=False)
    parser.add_argument("-plot_b_kSZ_kSZ_tSZ",help="bk" ,dest="plot_b_kSZ_kSZ_tSZ", type=str, required=False)
    parser.add_argument("-plot_bk_ttg",help="bk" ,dest="plot_bk_ttg", type=str, required=False)
    parser.add_argument("-plot_te_y_y",help="Tll" ,dest="plot_te_y_y", type=str, required=False)
    parser.add_argument("-plot_kSZ_kSZ_gal",help="plot_kSZ_kSZ_gal" ,dest="plot_kSZ_kSZ_gal", type=str, required=False)
    parser.add_argument("-plot_kSZ_kSZ_lensmag_1h",help="kSZ_kSZ_lensmag_1h" ,dest="plot_kSZ_kSZ_lensmag_1h", type=str, required=False)
    parser.add_argument("-plot_tSZ_tSZ_tSZ_1h",help="tSZ_tSZ_tSZ_1h" ,dest="plot_tSZ_tSZ_tSZ_1h", type=str, required=False)
    parser.add_argument("-plot_tSZ_lens",help="tSZ_lens" ,dest="plot_tSZ_lens", type=str, required=False)
    parser.add_argument("-plot_tSZ_gal",help="tSZ_gal" ,dest="plot_tSZ_gal", type=str, required=False)
    parser.add_argument("-plot_gal_gal",help="gal_gal" ,dest="plot_gal_gal", type=str, required=False)
    parser.add_argument("-plot_gal_lens",help="gal_lens" ,dest="plot_gal_lens", type=str, required=False)
    parser.add_argument("-plot_gal_gallens",help="gal_gallens" ,dest="plot_gal_gallens", type=str, required=False)
    parser.add_argument("-plot_isw_lens",help="isw_lens" ,dest="plot_isw_lens", type=str, required=False)
    parser.add_argument("-plot_lens_lens",help="lens_lens" ,dest="plot_lens_lens", type=str, required=False)
    parser.add_argument("-plot_isw_tsz",help="isw_tsz" ,dest="plot_isw_tsz", type=str, required=False)
    parser.add_argument("-plot_isw_auto",help="isw_auto" ,dest="plot_isw_auto", type=str, required=False)
    parser.add_argument("-plot_redshift_dependent_functions",help="redshift dependent functions" ,dest="plot_redshift_dependent_functions", type=str, required=False)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
	main()
