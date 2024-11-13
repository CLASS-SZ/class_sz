import pytest
import scipy.integrate as integrate
import numpy as np
from classy_sz import Class as Class_sz
import os 


h = 0.6774
Omega_b = 0.0486
Omega_c= (0.315-0.0486)
sigma8=0.8159
N_ur = 3.046
n_s = 0.9649
tau_reio=0.0543
N_ncdm=0
m_ncdm=0

def test_classy_sz_clgyappag():

   Mclass_sz = Class_sz()
   print('PATH_TO_CLASS_SZ_DATA:', os.environ['PATH_TO_CLASS_SZ_DATA'])

   Mclass_sz.set({

   'h': h,
   'Omega_b': Omega_b,
   'Omega_cdm':  Omega_c,
   'sigma8': sigma8,
   'n_s': n_s,
   'tau_reio': tau_reio,
   'N_ur': N_ur,
   'cosmo_model':0, 


   'output': 'tSZ_gal_1h,tSZ_gal_2h',

   'hm_consistency' : 1,  
   'M_min': 1e+10, # *cosmo['h'],
   'M_max':  1e+15, # *cosmo['h'],

   'pressure profile': 'B12',
   'x_outSZ': 4.,
   'truncate_gas_pressure_wrt_rvir': 1.,
   'delta for electron pressure':'200c',


   'galaxy_sample': 'custom',
   'delta_for_galaxies': "200c",
   'delta_for_matter_density': "200c",
   'mass_function': 'T08M200c',
   'concentration parameter': 'D08',

   'dell':200.0, 
   'multipoles_sz': 'ell_mock',
   'ell_max': 3000.0,
   'ell_min': 200.0,
   'z_min': 1.0e-8,
   'z_max': 2.0,

   'M0_HOD': 0, 
   'M_min_HOD':10.**11.74, # *cosmo['h'], #Msun/h
   'M1_prime_HOD':10.**13.24846, #*cosmo['h'], #Msun/h
   'sigma_log10M_HOD':0.3,
   'alpha_s_HOD':1.85,
   'f_cen_HOD': 1.,
   'x_out_truncated_nfw_profile_satellite_galaxies':1., # so corresponds to 1xr200c
   'csat_over_cdm' : 1.0,
   'full_path_to_dndz_gal':  os.environ['PATH_TO_CLASS_SZ_DATA'] + '/class_sz/class-sz/class_sz_auxiliary_files/includes/nz_maglim_bin2.txt', # lens galaxies
   'Delta_z_lens':0.00,
   'Delta_z_source':0.00,     

   'redshift_epsabs': 1.0e-40,
   'redshift_epsrel': 0.0005,
   'mass_epsabs': 1.0e-40,
   'mass_epsrel': 0.0005,
   'ndim_masses': 150,
   'ndim_redshifts': 150,
   'class_sz_verbose': 0,
   'nonlinear_verbose': 0,

   'use_fft_for_profiles_transform': 1,
   'N_samp_fftw': 1024, #precision parameter for the bessel transform to theta space
   'l_min_samp_fftw': 1e-12,
   'l_max_samp_fftw': 1e12,
   'x_min_gas_pressure_fftw': 1e-4,
   'x_max_gas_pressure_fftw': 1e3,

   'P_k_max_h/Mpc':100.0,
   'k_min_for_pk_class_sz': 0.0001,
   'k_max_for_pk_class_sz': 50.0,
   'k_per_decade_class_sz': 20.0,   

    })
   Mclass_sz.compute_class_szfast()


   gy_1h = 1e-6*np.asarray(Mclass_sz.cl_yg()['1h'])
   gy_2h = 1e-6*np.asarray(Mclass_sz.cl_yg()['2h'])
   ell_list = np.asarray(Mclass_sz.cl_yg()['ell'])


   #values computed with Shivam Pandey's independent code 
   ell_shivam = np.array([ 200.,  400.,  600.,  800., 1000., 1200., 1400., 1600., 1800.,
      2000., 2200., 2400., 2600., 2800., 3000., 3200., 3400., 3600.,
      3800., 4000., 4200., 4400., 4600., 4800., 5000., 5200., 5400.,
      5600., 5800., 6000., 6200., 6400., 6600., 6800., 7000., 7200.,
      7400., 7600., 7800., 8000., 8200., 8400., 8600., 8800., 9000.,
      9200., 9400., 9600., 9800.])
   gy_1h_shivam = np.array([2.40489236e-12, 2.22051755e-12, 2.01025631e-12, 1.79322157e-12,
      1.58424355e-12, 1.39074152e-12, 1.21531649e-12, 1.05823098e-12,
      9.18825628e-13, 7.96094622e-13, 6.88823232e-13, 5.95634219e-13,
      5.15060728e-13, 4.45640098e-13, 3.85995867e-13, 3.34877405e-13,
      2.91165373e-13, 2.53857132e-13, 2.22051555e-13, 1.94939584e-13,
      1.71804315e-13, 1.52021991e-13, 1.35059684e-13, 1.20466655e-13,
      1.07864084e-13, 9.69356641e-14, 8.74187120e-14, 7.90960695e-14,
      7.17889666e-14, 6.53488681e-14, 5.96515045e-14, 5.45918310e-14,
      5.00807465e-14, 4.60429027e-14, 4.24150149e-14, 3.91442105e-14,
      3.61865740e-14, 3.35052926e-14, 3.10693259e-14, 2.88520623e-14,
      2.68304747e-14, 2.49844586e-14, 2.32962518e-14, 2.17502818e-14,
      2.03330274e-14, 1.90324510e-14, 1.78381045e-14, 1.67406600e-14,
      1.57313969e-14])
   gy_2h_shivam = np.array([1.87088117e-12, 5.11641281e-13, 2.14550366e-13, 1.09884223e-13,
      6.32605917e-14, 3.93470652e-14, 2.58700200e-14, 1.77360121e-14,
      1.25657606e-14, 9.14253820e-15, 6.80016787e-15, 5.15327094e-15,
      3.96858662e-15, 3.09954177e-15, 2.45110577e-15, 1.96000365e-15,
      1.58311871e-15, 1.29045896e-15, 1.06077100e-15, 8.78742483e-16,
      7.33181576e-16, 6.15805993e-16, 5.20419212e-16, 4.42337937e-16,
      3.77988401e-16, 3.24618810e-16, 2.80092470e-16, 2.42737384e-16,
      2.11235494e-16, 1.84540043e-16, 1.61813663e-16, 1.42381825e-16,
      1.25697744e-16, 1.11315776e-16, 9.88710325e-17, 8.80634676e-17,
      7.86453576e-17, 7.04112952e-17, 6.31901791e-17, 5.68388478e-17,
      5.12369821e-17, 4.62830171e-17, 4.18908971e-17, 3.79874956e-17,
      3.45104892e-17, 3.14065931e-17, 2.86301043e-17, 2.61416663e-17,
      2.39072489e-17])

   for i, l in enumerate(ell_list[:14]):
      f = l*(l+1)/2/np.pi
      diff_1h = np.abs(gy_1h[i]-f*gy_1h_shivam[i])
      diff_2h = np.abs(gy_2h[i]-f*gy_2h_shivam[i])
      print(f" Relative error in '1h' term of g x y at ell={l} is {diff_1h/gy_1h[i]}, expected < 2e-2")
      print(f" Relative error in '2h' term of g x y at ell={l} is {diff_2h/gy_2h[i]}, expected  < 7e-2")
      assert f*gy_1h_shivam[i]  == pytest.approx(gy_1h[i], rel=2.5e-2), f" Relative error in '1h' term of g x y at ell={l} is {diff_1h/gy_1h[i]}, expected < 2e-2"
      assert f*gy_2h_shivam[i]  == pytest.approx(gy_2h[i], rel=7e-2), f" Relative error in '2h' term of g x y at ell={l} is {diff_2h/gy_2h[i]}, expected  < 7e-2"

