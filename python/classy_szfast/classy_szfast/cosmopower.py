from .utils import *
from .config import *

cosmopower_derived_params_names = ['100*theta_s',
                                   'sigma8',
                                   'YHe',
                                   'z_reio',
                                   'Neff',
                                   'tau_rec',
                                   'z_rec',
                                   'rs_rec',
                                   'ra_rec',
                                   'tau_star',
                                   'z_star',
                                   'rs_star',
                                   'ra_star',
                                   'rs_drag']

cp_l_max_scalars = 11000 # max multipole of train ing data



emulator_dict = {}
emulator_dict['lcdm'] = {}
emulator_dict['mnu'] = {}
emulator_dict['neff'] = {}
emulator_dict['wcdm'] = {}




emulator_dict['lcdm']['TT'] = 'TT_v1'
emulator_dict['lcdm']['TE'] = 'TE_v1'
emulator_dict['lcdm']['EE'] = 'EE_v1'
emulator_dict['lcdm']['PP'] = 'PP_v1'
emulator_dict['lcdm']['PKNL'] = 'PKNL_v1'
emulator_dict['lcdm']['PKL'] = 'PKL_v1'
emulator_dict['lcdm']['DER'] = 'DER_v1'
emulator_dict['lcdm']['DAZ'] = 'DAZ_v1'
emulator_dict['lcdm']['HZ'] = 'HZ_v1'
emulator_dict['lcdm']['S8Z'] = 'S8Z_v1'

emulator_dict['mnu']['TT'] = 'TT_mnu_v1'
emulator_dict['mnu']['TE'] = 'TE_mnu_v1'
emulator_dict['mnu']['EE'] = 'EE_mnu_v1'
emulator_dict['mnu']['PP'] = 'PP_mnu_v1'
emulator_dict['mnu']['PKNL'] = 'PKNL_mnu_v1'
emulator_dict['mnu']['PKL'] = 'PKL_mnu_v1'
emulator_dict['mnu']['DER'] = 'DER_mnu_v1'
emulator_dict['mnu']['DAZ'] = 'DAZ_mnu_v1'
emulator_dict['mnu']['HZ'] = 'HZ_mnu_v1'
emulator_dict['mnu']['S8Z'] = 'S8Z_mnu_v1'


emulator_dict['neff']['TT'] = 'TT_neff_v1'
emulator_dict['neff']['TE'] = 'TE_neff_v1'
emulator_dict['neff']['EE'] = 'EE_neff_v1'
emulator_dict['neff']['PP'] = 'PP_neff_v1'
emulator_dict['neff']['PKNL'] = 'PKNL_neff_v1'
emulator_dict['neff']['PKL'] = 'PKL_neff_v1'
emulator_dict['neff']['DER'] = 'DER_neff_v1'
emulator_dict['neff']['DAZ'] = 'DAZ_neff_v1'
emulator_dict['neff']['HZ'] = 'HZ_neff_v1'
emulator_dict['neff']['S8Z'] = 'S8Z_neff_v1'



emulator_dict['wcdm']['TT'] = 'TT_w_v1'
emulator_dict['wcdm']['TE'] = 'TE_w_v1'
emulator_dict['wcdm']['EE'] = 'EE_w_v1'
emulator_dict['wcdm']['PP'] = 'PP_w_v1'
emulator_dict['wcdm']['PKNL'] = 'PKNL_w_v1'
emulator_dict['wcdm']['PKL'] = 'PKL_w_v1'
emulator_dict['wcdm']['DER'] = 'DER_w_v1'
emulator_dict['wcdm']['DAZ'] = 'DAZ_w_v1'
emulator_dict['wcdm']['HZ'] = 'HZ_w_v1'
emulator_dict['wcdm']['S8Z'] = 'S8Z_w_v1'



cp_tt_nn = {}
cp_te_nn = {}
cp_ee_nn = {}
cp_pp_nn = {}
cp_pknl_nn = {}
cp_pkl_nn = {}
cp_der_nn = {}
cp_da_nn = {}
cp_h_nn = {}
cp_s8_nn = {}

for mp in ['lcdm','mnu','neff','wcdm']:
    path_to_emulators = path_to_cosmopower_organization + mp +'/'

    cp_tt_nn[mp] = cosmopower_NN(restore=True,
                             restore_filename=path_to_emulators + 'TTTEEE/' + emulator_dict[mp]['TT'])

    cp_te_nn[mp] = cosmopower_PCAplusNN(restore=True,
                                    restore_filename=path_to_emulators + 'TTTEEE/' + emulator_dict[mp]['TE'])

    cp_ee_nn[mp] = cosmopower_NN(restore=True,
                             restore_filename=path_to_emulators + 'TTTEEE/' + emulator_dict[mp]['EE'])

    cp_pp_nn[mp] = cosmopower_NN(restore=True,
                             restore_filename=path_to_emulators + 'PP/' + emulator_dict[mp]['PP'])

    cp_pknl_nn[mp] = cosmopower_NN(restore=True,
                               restore_filename=path_to_emulators + 'PK/' + emulator_dict[mp]['PKNL'])

    cp_pkl_nn[mp] = cosmopower_NN(restore=True,
                              restore_filename=path_to_emulators + 'PK/' + emulator_dict[mp]['PKL'])

    cp_der_nn[mp] = cosmopower_NN(restore=True,
                              restore_filename=path_to_emulators + 'derived-parameters/' + emulator_dict[mp]['DER'])

    cp_da_nn[mp] = cosmopower_NN(restore=True,
                             restore_filename=path_to_emulators + 'growth-and-distances/' + emulator_dict[mp]['DAZ'])

    cp_h_nn[mp] = cosmopower_NN(restore=True,
                            restore_filename=path_to_emulators + 'growth-and-distances/' + emulator_dict[mp]['HZ'])

    cp_s8_nn[mp] = cosmopower_NN(restore=True,
                             restore_filename=path_to_emulators + 'growth-and-distances/' + emulator_dict[mp]['S8Z'])
