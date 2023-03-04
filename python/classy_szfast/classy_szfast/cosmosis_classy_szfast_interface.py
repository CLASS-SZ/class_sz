from builtins import str
import os
from cosmosis.datablock import names, option_section
import sys
import traceback
import warnings

# add class directory to the path
dirname = os.path.split(__file__)[0]
# enable debugging from the same directory
if not dirname.strip():
    dirname = '.'

import numpy as np

# These are pre-defined strings we use as datablock
# section names
cosmo = names.cosmological_parameters
distances = names.distances
cmb_cl = names.cmb_cl


def setup(options):
    # class_version = options.get_string(option_section, "version", "3.2.0")
    # pyversion = f"{sys.version_info.major}.{sys.version_info.minor}"
    # install_dir = dirname + f"/class_v{class_version}/classy_install/lib/python{pyversion}/site-packages/"
    # with open(f"{install_dir}/easy-install.pth") as f:
    #     pth = f.read().strip()
    #     install_dir = install_dir + pth
    # install_dir = '/Users/boris/opt/miniconda3/lib/python3.9/site-packages/classy_sz-2.9.4-py3.9-macosx-11.0-arm64.egg/classy_sz.cpython-39-darwin.so'
    # sys.path.insert(0, install_dir)

    import classy_sz as classy
    print(f"Loaded classy_sz from {classy.__file__}")

    # Read options from the ini file which are fixed across
    # the length of the chain
    config = {
        'lmax': options.get_int(option_section, 'lmax', default=2000),
        'zmax': options.get_double(option_section, 'zmax', default=3.0),
        'kmax': options.get_double(option_section, 'kmax', default=50.0),
        'debug': options.get_bool(option_section, 'debug', default=False),
        'lensing': options.get_bool(option_section, 'lensing', default=True),
        'cmb': options.get_bool(option_section, 'cmb', default=True),
        'mpk': options.get_bool(option_section, 'mpk', default=True),
        'save_matter_power_lin': options.get_bool(option_section, 'save_matter_power_lin', default=True),
        'save_cdm_baryon_power_lin': options.get_bool(option_section, 'save_cdm_baryon_power_lin', default=False),
    }


    for _, key in options.keys(option_section):
        if key.startswith('class_'):
            config[key] = options[option_section, key]


    # Create the object that connects to Class
    config['cosmo'] = classy.Class()

    # Return all this config information
    return config

def choose_outputs(config):
    outputs = []
    if config['cmb']:
        outputs.append("tCl pCl")
    if config['lensing']:
        outputs.append("lCl")
    if config["mpk"]:
        outputs.append("mPk")
    return " ".join(outputs)

def get_class_inputs(block, config):

    # Get parameters from block and give them the
    # names and form that class expects
    nnu = block.get_double(cosmo, 'nnu', 3.046)
    nmassive = block.get_int(cosmo, 'num_massive_neutrinos', default=0)

    # print(nmassive)
    # print(block.get_double(cosmo, 'mnu', default=0.06))
    # exit(0)

    params = {
        'output': choose_outputs(config),
        'lensing':   'yes' if config['lensing'] else 'no',
        'A_s':       block[cosmo, 'A_s'],
        'n_s':       block[cosmo, 'n_s'],
        'H0':        100 * block[cosmo, 'h0'],
        'omega_b':   block[cosmo, 'ombh2'],
        'omega_cdm': block[cosmo, 'omch2'],
        'tau_reio':  block[cosmo, 'tau'],
        'T_cmb':     block.get_double(cosmo, 'TCMB', default=2.726),
        # 'N_ur':      nnu - nmassive,
        # 'N_ncdm':    nmassive,
        # 'm_ncdm':    block.get_double(cosmo, 'mnu', default=0.06)

     # CLASS SETTINGS FOR COSMOPOWER
      'N_ncdm': 1,
      'N_ur': 2.0328,
      'm_ncdm': 0.06,

     # class_sz fast :
     'ndim_redshifts' : 25 # number of z's in our pk interpolator

    }
    # print(params)
    # exit(0)

    if config["cmb"] or config["lensing"]:
        params.update({
          'l_max_scalars': config["lmax"],
        })


    if config["mpk"]:
        params.update({
            'P_k_max_h/Mpc':  config["kmax"],
            'z_pk': ', '.join(str(z) for z in np.arange(0.0, config['zmax'], 0.01)),
            'z_max_pk': config['zmax'],
        })



    if block.has_value(cosmo, "massless_nu"):
        warnings.warn("Parameter massless_nu is being ignored. Set nnu, the effective number of relativistic species in the early Universe.")

    if (block.has_value(cosmo, "omega_nu") or block.has_value(cosmo, "omnuh2")) and not (block.has_value(cosmo, "mnu")):
        warnings.warn("Parameter omega_nu and omnuh2 are being ignored. Set mnu and num_massive_neutrinos instead.")


    for key,val in config.items():
        if key.startswith('class_'):
            params[key[6:]] = val


    return params




def get_class_outputs(block, c, config):
    ##
    # Derived cosmological parameters
    ##

    h0 = block[cosmo, 'h0']

    ##
    # Matter power spectrum
    ##

    # Ranges of the redshift and matter power
    dz = 0.01
    kmin = 1e-4
    kmax = config['kmax'] * h0
    nk = 100

    # Define k,z we want to sample
    z = np.arange(0.0, config["zmax"] + dz, dz)
    k = np.logspace(np.log10(kmin), np.log10(kmax), nk)
    nz = len(z)

    # Extract (interpolate) P(k,z) at the requested
    # sample points.
    # print(c.pars['output'])
    # exit(0)
    # if 'mPk' in c.pars['output']:
    if 1==1:
        #print('trying to get mPk')
        #print('first getting sigma8')
        # exit(0)
        block[cosmo, 'sigma_8'] = c.sigma8()
        #print('got sigma8:',c.sigma8())
        #exit(0)

        # Total matter power spectrum (saved as grid)
        if config['save_matter_power_lin']:
            P = np.zeros((k.size, z.size))
            for i, ki in enumerate(k):
                for j, zi in enumerate(z):
                    P[i, j] = c.pk_lin(ki, zi)
            # print(np.shape(P * h0**3))
            # exit(0)
            # print('matter Pk:',P)
            # exit(0)
            block.put_grid("matter_power_lin", "k_h", k / h0, "z", z, "P_k", P * h0**3)

        # CDM+baryons power spectrum
        if config['save_cdm_baryon_power_lin']:
            # P = np.zeros((k.size, z.size))
            # for i, ki in enumerate(k):
            #     for j, zi in enumerate(z):
            #         P[i, j] = c.pk_cb_lin(ki, zi)
            # block.put_grid('cdm_baryon_power_lin', 'k_h', k/h0, 'z', z, 'P_k', P*h0**3)\
            print('save_cdm_baryon_power_lin not available with class_szfast yet')
            exit(0)

        # Get growth rates and sigma_8
        D = [c.scale_independent_growth_factor(zi) for zi in z] # original class routine
        # print('got D:',D)
        f = [c.scale_independent_growth_factor_f(zi) for zi in z] # original class routine
        # print('got f:',f)
        # fsigma = [c.effective_f_sigma8(zi) for zi in z] # get this from the emulators
        fsigma = [c.get_effective_f_sigma8(zi) for zi in z]
        # print('got fsigma8:',fsigma)
        # sigma_8_z = [c.sigma(8.0/h0, zi, h_units=True) for zi in z] # get this from the emulators, had to add h0 by hand in class v2.9.4 the hunits thing is not treated properly
        sigma_8_z = [c.get_sigma8_at_z(zi) for zi in z] # get this from the emulators, had to add h0 by hand in class v2.9.4 the hunits thing is not treated properly
        # print('sigma8z:',sigma_8_z)
        # print('z:',z)
        # exit(0)
        block[names.growth_parameters, "z"] = z
        block[names.growth_parameters, "sigma_8"] = np.array(sigma_8_z)
        block[names.growth_parameters, "fsigma_8"] = np.array(fsigma)
        block[names.growth_parameters, "d_z"] = np.array(D)
        block[names.growth_parameters, "f_z"] = np.array(f)
        block[names.growth_parameters, "a"] = 1/(1+z)

        #print('getting nonlinear')
        #print(c.nonlinear_method)
        #print(config)
        tmp_nonlinear_method = c.nonlinear_method
        if 'class_non_linear' in config:
            tmp_nonlinear_method = 1
            #print(config['class_non_linear'])

        if tmp_nonlinear_method != 0:
            #print('entering nonlinear block')
            for i, ki in enumerate(k):
                for j, zi in enumerate(z):
                    P[i, j] = c.pk(ki, zi)
            # print('got Pknl:',P)
            block.put_grid("matter_power_nl", "k_h", k / h0, "z", z, "P_k", P * h0**3)
        #     exit(0)
        # exit(0)


    ##
    # Distances and related quantities
    ##

    # print("moving on")
    # save redshifts of samples
    block[distances, 'z'] = z
    block[distances, 'nz'] = nz
    block[distances, 'a'] = 1/(1+z)
    # Save distance samples
    # print("moving on")
    d_a = np.array([c.angular_distance(zi) for zi in z]) # original class routine
    block[distances, 'd_a'] = d_a # original class routine
    # block[distances, 'd_l'] = np.array([c.luminosity_distance(zi) for zi in z])
    # print("moving on after ang")
    block[distances, 'd_l'] = d_a * (1 + z)**2 # original class routine
    block[distances, 'd_m'] = d_a * (1 + z) # original class routine
    # print("moving on after distances")
    # Save some auxiliary related parameters
    block[distances, 'age'] = c.age() # original class routine
    # print("moving on after age:",c.age())
    block[distances, 'rs_zdrag'] = c.rs_drag() # original class routine
    # print("moving on done after drag:",c.rs_drag())
    # exit(0)
    ##
    # Now the CMB C_ell
    ##
    if config["cmb"]:
        c_ell_data = c.lensed_cl() if config['lensing'] else c.raw_cl()
        ell = c_ell_data['ell']
        ell = ell[2:]

        # Save the ell range
        block[cmb_cl, "ell"] = ell

        # t_cmb is in K, convert to mu_K, and add ell(ell+1) factor
        tcmb_muk = block[cosmo, 'tcmb'] * 1e6
        f = ell * (ell + 1.0) / 2 / np.pi * tcmb_muk**2

        # Save each of the four spectra
        for s in ['tt', 'ee', 'te', 'bb']:
            block[cmb_cl, s] = c_ell_data[s][2:] * f


def execute(block, config):
    import classy_sz as classy
    c = config['cosmo']

    try:
        # Set input parameters
        params = get_class_inputs(block, config)
        params_values = params.copy()
        A_s = params_values['A_s']
        lnAs = np.log(1e10*A_s)
        params_values['ln10^{10}A_s'] = lnAs
        params_values.pop('A_s')
        params_values['output'] = '' # remove all outputs cause things are taken care of by the emulators
        params_values.pop('non_linear') # also pretend we dont compute mPk...
        params_values.pop('z_max_pk')
        params_values.pop('P_k_max_h/Mpc')
        params_values.pop('z_pk')
        # params_values['ndim_redshifts'] =
        # print(params_values)
        # exit(0)
        c.set(params_values)
        # exit(0)


        # Run calculations
        # c.compute()
        # run the fast class_sz calculation using emulators
        # print('computing class_sz fast:')
        c.compute_class_szfast()
        # print('class_sz fast computed')
        # exit(0)
        # restore the original parameter file:
        # c.set(params)

        # Extract outputs
        get_class_outputs(block, c, config)
    except classy.CosmoError as error:
        if config['debug']:
            sys.stderr.write("Error in class or class_sz. You set debug=T so here is more debug info:\n")
            traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write("Error in class or class_sz. Set debug=T for info: {}\n".format(error))
        return 1
    finally:
        # Reset for re-use next time
        c.struct_cleanup()
    return 0


def cleanup(config):
    config['cosmo'].empty()
