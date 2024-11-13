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

def test_classy_sz_clgkappag():
    
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
    # 'cosmo_model':0, 


    'output': 'gamma_gal_gallens_1h,gamma_gal_gallens_2h,',

    'hm_consistency' : 1,  
    'M_min': 1e+10, # *cosmo['h'],
    'M_max':  1e+15, # *cosmo['h'],
        

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
    'full_path_to_source_dndz_gal':os.environ['PATH_TO_CLASS_SZ_DATA'] + '/class_sz/class-sz/class_sz_auxiliary_files/includes/nz_source_normalized_bin4.txt', # source galaxies
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


    'P_k_max_h/Mpc':100.0,
    'k_min_for_pk_class_sz': 0.0001,
    'k_max_for_pk_class_sz': 50.0,
    'k_per_decade_class_sz': 20.0,   

    })
    Mclass_sz.compute()


    gk_1h = np.asarray(Mclass_sz.cl_ggamma()['1h'])
    gk_2h = np.asarray(Mclass_sz.cl_ggamma()['2h'])
    ell_list = np.asarray(Mclass_sz.cl_ggamma()['ell'])

   #  print(gk_1h)
   #  print(gk_2h)


    #values computed with Shivam Pandey's independent code 
    ell_shivam = np.array([ 200.,  400.,  600.,  800., 1000., 1200., 1400., 1600., 1800.,
       2000., 2200., 2400., 2600., 2800., 3000., 3200., 3400., 3600.,
       3800., 4000., 4200., 4400., 4600., 4800., 5000., 5200., 5400.,
       5600., 5800., 6000., 6200., 6400., 6600., 6800., 7000., 7200.,
       7400., 7600., 7800., 8000., 8200., 8400., 8600., 8800., 9000.,
       9200., 9400., 9600., 9800.])
    gk_1h_shivam = np.array([1.01896624e-08, 9.89574719e-09, 9.48411630e-09, 8.97669740e-09,
       8.40269021e-09, 7.79091071e-09, 7.16646760e-09, 6.54951241e-09,
       5.95513837e-09, 5.39386474e-09, 4.87233987e-09, 4.39409168e-09,
       3.96023157e-09, 3.57008010e-09, 3.22170296e-09, 2.91235138e-09,
       2.63881329e-09, 2.39767996e-09, 2.18553872e-09, 1.99910440e-09,
       1.83530160e-09, 1.69130986e-09, 1.56458250e-09, 1.45284783e-09,
       1.35409889e-09, 1.26657684e-09, 1.18875095e-09, 1.11929720e-09,
       1.05707646e-09, 1.00111329e-09, 9.50575638e-10, 9.04755664e-10,
       8.63052102e-10, 8.24954324e-10, 7.90028179e-10, 7.57903636e-10,
       7.28264216e-10, 7.00838089e-10, 6.75390653e-10, 6.51718405e-10,
       6.29643909e-10, 6.09011683e-10, 5.89684866e-10, 5.71542495e-10,
       5.54477277e-10, 5.38393771e-10, 5.23206900e-10, 5.08840706e-10,
       4.95227290e-10])
    gk_2h_shivam = np.array([6.03185563e-08, 1.68868641e-08, 7.32773562e-09, 3.91512365e-09,
       2.36352932e-09, 1.54622644e-09, 1.07096197e-09, 7.74026174e-10,
       5.78222181e-10, 4.43536906e-10, 3.47699869e-10, 2.77585814e-10,
       2.25087403e-10, 1.84997969e-10, 1.53862070e-10, 1.29321083e-10,
       1.09725430e-10, 9.38967242e-11, 8.09773673e-11, 7.03328757e-11,
       6.14869287e-11, 5.40772448e-11, 4.78250291e-11, 4.25134523e-11,
       3.79722511e-11, 3.40665505e-11, 3.06886443e-11, 2.77518783e-11,
       2.51860494e-11, 2.29339100e-11, 2.09484888e-11, 1.91910187e-11,
       1.76293243e-11, 1.62365563e-11, 1.49901938e-11, 1.38712516e-11,
       1.28636481e-11, 1.19536977e-11, 1.11297021e-11, 1.03816191e-11,
       9.70079337e-12, 9.07973571e-12, 8.51194274e-12, 7.99174751e-12,
       7.51419592e-12, 7.07494383e-12, 6.67017093e-12, 6.29650842e-12,
       5.95097789e-12])

    for i, l in enumerate(ell_list[:14]):
        f = l*(l+1)/2/np.pi
        diff_1h = gk_1h[i]-f*gk_1h_shivam[i]
        diff_2h = gk_2h[i]-f*gk_2h_shivam[i]
        print("")
        assert f*gk_1h_shivam[i]  == pytest.approx(gk_1h[i], rel=2e-2), f" Relative error in '1h' term of g x kappa_g at ell={l} is {diff_1h/gk_1h[i]}, expected < 2e-2"
        assert f*gk_2h_shivam[i]  == pytest.approx(gk_2h[i], rel=7e-2), f" Relative error in '2h' term of g x kappa_g at ell={l} is {diff_2h/gk_2h[i]}, expected  < 7e-2"

