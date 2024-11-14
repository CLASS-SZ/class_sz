import pytest
import scipy.integrate as integrate
import numpy as np
from classy_sz import Class as Class_sz
import os 
import time


print('PATH_TO_CLASS_SZ_DATA:', os.environ['PATH_TO_CLASS_SZ_DATA'])

# Dictionary of parameters
pdict = {
    'A_ym': 4.949940743969354e-05,
    'B': 1.25,
    'B_ym': 0.12329999999999997,
    'C_ym': 0.0,
    'HMF_prescription_NCDM': 1,
    'M_max': 8400000000000000.0,
    'M_min': 5600000000000.0,
    'N_samp_fftw': 6400,
    'concentration_parameter': 'B13',
    'cosmo_model': 1,
    'experiment': 1,
    'h': 0.7,
    'has_selection_function': 1,
    'm_ncdm': 0.06,
    'm_pivot_ym_[Msun]': 300000000000000.0,
    'mass_function': 'T08M500c',
    'mass_epsabs': 1e-100,
    'mass_epsrel': 1e-06,
    'n_m_dndlnM': 8000,
    'n_s': 0.96,
    'ndim_masses': 100,
    'ndim_redshifts': 5000,
    'no_spline_in_tinker': 1,
    'ntab_dlnm_dlnq': 1000,
    'omega_b': 0.023995299999999997,
    'omega_cdm': 0.1303547,
    'output': 'sz_cluster_counts_fft',
    'redshift_epsabs': 1e-100,
    'redshift_epsrel': 1e-06,
    'sigma8': 0.811,
    'sigmaM_ym': 0.173,
    'signal-to-noise_cut-off_for_survey_cluster_completeness': 5.0,
    'sz_selection_function_skyfracs_file': os.environ['PATH_TO_CLASS_SZ_DATA'] + '/class_sz/class-sz/class_sz_auxiliary_files/includes/so_sim_sz_mf_noise_skyfracs.txt',
    'sz_selection_function_thetas_file': os.environ['PATH_TO_CLASS_SZ_DATA'] + '/class_sz/class-sz/class_sz_auxiliary_files/includes/so_sim_sz_mf_noise_thetas.txt',
    'sz_selection_function_ylims_file': os.environ['PATH_TO_CLASS_SZ_DATA'] + '/class_sz/class-sz/class_sz_auxiliary_files/includes/so_sim_sz_mf_noise_sigmas.txt',
    'szcc_dof': 0.0,
    'szcc_qtrunc': -1.0,
    'szcounts_fft_nz': 550,
    'szcounts_fft_z_max': 3.5999999999999996,
    'szcounts_fft_z_min': 0.008,
    'szcounts_qmax_fft_padded': 500.0,
    'tau_reio': 0.06,
    'tol_dlnm_dlnq': 0.001,
    'use_m200c_in_ym_relation': 0,
    'use_m500c_in_ym_relation': 1,
    'use_skyaveraged_noise': 0,
    'y_m_relation': 1,
    'z_max': 3.5999999999999996,
    'z_min': 0.008
}

pdict_binned = {'A_ym': 4.949940743969354e-05,
 'B': 1.25,
 'B_ym': 0.12329999999999997,
 'C_ym': 0.0,
 'HMF_prescription_NCDM': 1,
 'M_max': 8400000000000000.0,
 'M_min': 5600000000000.0,
 'N_ncdm': 1,
 'N_samp_fftw': 6400,
 'N_ur': 2.0328,
 'T_ncdm': 0.71611,
 'bin_dlog10_snr': 0.25,
 'bin_dz_cluster_counts': 0.1,
 'bin_z_max_cluster_counts': 2.01,
 'bin_z_min_cluster_counts': 0.01,
 'cluster_count_completeness_grid_z_cutoff_low': 0.4,
 'cluster_count_completeness_grid_z_cutoff_mid': 1.0,
 'concentration parameter': 'B13',
 'cosmo_model': 0,
 'dlnM_cluster_count_completeness_grid': 0.01,
 'dlny': 0.025,
 'dz_cluster_count_completeness_grid_high_z': 0.01,
 'dz_cluster_count_completeness_grid_low_z': 0.001,
 'dz_cluster_count_completeness_grid_mid_z': 0.001,
 'experiment': 1,
 'h': 0.7,
 'has_selection_function': 1,
 'lnymax': 10.0,
 'lnymin': -20.0,
 'log10_snr_max': 2.0,
 'log10_snr_min': 0.6,
 'm_ncdm': 0.06,
 'm_pivot_ym_[Msun]': 300000000000000.0,
 'mass function': 'T08M500c',
 'mass_epsabs': 1e-100,
 'mass_epsrel': 1e-06,
 'n_m_dndlnM': 8000,
 'n_s': 0.96,
 'ndim_masses': 100,
 'ndim_redshifts': 5000,
 'no_spline_in_tinker': 1,
 'ntab_dlnm_dlnq': 1000,
 'omega_b': 0.023995299999999997,
 'omega_cdm': 0.1303547,
 'output': 'sz_cluster_counts',
 'redshift_epsabs': 1e-100,
 'redshift_epsrel': 1e-06,
 'sigma8': 0.811,
 'sigmaM_ym': 0.173,
 'sigma_derivative': 0,
 'signal-to-noise_cut-off_for_survey_cluster_completeness': 5.0,
 'sz_selection_function_skyfracs_file': os.environ['PATH_TO_CLASS_SZ_DATA'] + '/class_sz/class-sz/class_sz_auxiliary_files/includes/so_sim_sz_mf_noise_skyfracs.txt',
 'sz_selection_function_thetas_file': os.environ['PATH_TO_CLASS_SZ_DATA'] + '/class_sz/class-sz/class_sz_auxiliary_files/includes/so_sim_sz_mf_noise_thetas.txt',
 'sz_selection_function_ylims_file': os.environ['PATH_TO_CLASS_SZ_DATA'] + '/class_sz/class-sz/class_sz_auxiliary_files/includes/so_sim_sz_mf_noise_sigmas.txt',
 'szcc_dof': 0.0,
 'szcc_qtrunc': -1.0,
 'szcounts_fft_nz': 550,
 'szcounts_fft_z_max': 3.5999999999999996,
 'szcounts_fft_z_min': 0.008,
 'szcounts_qmax_fft_padded': 500.0,
 'tau_reio': 0.06,
 'tol_dlnm_dlnq': 0.001,
 'use_m200c_in_ym_relation': 0,
 'use_m500c_in_ym_relation': 1,
 'use_skyaveraged_noise': 0,
 'y_m_relation': 1,
 'z_max': 3.5999999999999996,
 'z_min': 0.008}



def get_Ntot(classy_sz):
    # Define parameters for redshift and q arrays
    zmin = 0.01
    zmax = 3.0
    nz = 585
    z_arr = np.linspace(zmin, zmax, nz)

    q_threshold = 5.0
    q_max = 200.0
    nq = 5000
    q_arr = np.geomspace(q_threshold, q_max, nq)

    # Define the function to integrate
    get_dndzdq = np.vectorize(classy_sz.get_szcounts_dndzdq_at_z_q)

    # Integrate over z
    Nz = []
    for zp in z_arr:
        Nz.append(integrate.simpson(get_dndzdq(zp, q_arr), x=q_arr, axis=0))
    Nz = np.asarray(Nz)

    # Integrate over q
    Nq = []
    for qp in q_arr:
        Nq.append(integrate.simpson(get_dndzdq(z_arr, qp), x=z_arr, axis=0))
    Nq = np.asarray(Nq)

    # Total integral over redshift and q
    Ntot = np.trapz(Nz, x=z_arr)

    print("get Ntot:", Ntot)
    return Ntot


def get_Ntot_binned(classy_sz):
    dNdzdy_theoretical = classy_sz.dndzdy_theoretical()['dndzdy']
    z_center = classy_sz.dndzdy_theoretical()['z_center']
    z_edges = classy_sz.dndzdy_theoretical()['z_edges']
    log10y_center = classy_sz.dndzdy_theoretical()['log10y_center']
    log10y_edges = classy_sz.dndzdy_theoretical()['log10y_edges']

    N_z,N_y = np.shape(dNdzdy_theoretical)
    N_clusters_z_theory = []
    for iz in range(N_z):
        N_clusters_z_theory.append(np.sum(dNdzdy_theoretical[iz][0:]))
    N_clusters_y_theory = []

    for iy in range(N_y):
        N_clusters_y_theory.append(np.sum(np.asarray(dNdzdy_theoretical)[:,iy]))
        
    # print total number of clusters: 
    Ntot = np.sum(N_clusters_z_theory)

    print("get Ntot binned:", Ntot)
    return Ntot


# Function to run the test
def test_sz_cluster_counts_cosmo_model_1():
    start_time = time.time()
    # Initialize the Class object
    classy_sz = Class_sz()


    # Set the parameters in the Class object
    classy_sz.set(pdict)

    # Compute the results
    classy_sz.compute_class_szfast()

    Ntot = get_Ntot(classy_sz)


    # Check the result is close to the expected value (15017)
    assert Ntot == pytest.approx(15017, rel=0.01), f"Ntot is {Ntot}, expected close to 15017"
    # Print message if the test passes
    print(f"Test passed: Ntot is {Ntot}, close to 15017")
    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Test executed in {elapsed_time:.2f} seconds")



# Function to run the test
def test_sz_cluster_counts_cosmo_model_0():
    start_time = time.time()
    # Initialize the Class object
    classy_sz = Class_sz()


    # Set the parameters in the Class object
    pdict['cosmo_model'] = 0
    classy_sz.set(pdict)

    # Compute the results
    classy_sz.compute_class_szfast()

    Ntot = get_Ntot(classy_sz)


    # Check the result is close to the expected value (15017)
    assert Ntot == pytest.approx(15003, rel=0.01), f"Ntot is {Ntot}, expected close to 15003"
    # Print message if the test passes
    print(f"Test passed: Ntot is {Ntot}, close to 15003")

    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Test executed in {elapsed_time:.2f} seconds")



# Function to run the test
def test_sz_cluster_counts_binned():
    start_time = time.time()
    # Initialize the Class object
    classy_sz = Class_sz()

    # Set the parameters in the Class object
    classy_sz.set(pdict_binned)

    # Compute the results
    classy_sz.compute_class_szfast()

    Ntot = get_Ntot_binned(classy_sz)

    # Check the result is close to the expected value (15017)
    assert Ntot == pytest.approx(14983.433909, rel=0.01), f"Ntot is {Ntot}, expected close to 14983.433909"
    # Print message if the test passes
    print(f"Test passed for binned: Ntot is {Ntot}, close to 14983.433909")

    # End timing
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Test executed in {elapsed_time:.2f} seconds")