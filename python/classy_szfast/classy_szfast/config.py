import subprocess
import os


# set some paths, e.g.:
# This path needs to be adjsuted (TBD: set this automatically)
#path_to_sd_projects = "/scratch/nas_chluba/specdist/bolliet/specdist_ml/"


#the path to save the output from cosmotherm
#path_to_recfast_results =  path_to_sd_projects + "specdist/specdist/recfast_outputs"

# make a directory
#subprocess.call(['mkdir','-p',path_to_recfast_results])



root_path = os.path.abspath("")
# assuming cosmopower organization codes are one level up:
# path_to_cosmopower_organization = root_path + '/../cosmopower-organization/'
# path_to_cosmopower_organization = '/Users/boris/Work/CLASS-SZ/SO-SZ/cosmopower-organization/'
path_to_cosmopower_organization = 'cosmopower-organization/'
# print(path_to_cosmopower_organization)

# path_to_emulators = path_to_cosmopower_organization + 'lcdm/'
#
# # print(path_to_emulators)
#
# path_to_emulators = path_to_cosmopower_organization + 'lcdm/'
# str_cmd_subprocess = ["ls",path_to_emulators]
# print('inside lcdm:')
# print(subprocess.call(str_cmd_subprocess))
#
# path_to_emulators = path_to_cosmopower_organization + 'mnu/'
# str_cmd_subprocess = ["ls",path_to_emulators]
# print('\ninside mnu:')
# subprocess.call(str_cmd_subprocess)
#
# path_to_emulators = path_to_cosmopower_organization + 'wcdm/'
# str_cmd_subprocess = ["ls",path_to_emulators]
# print('\ninside wcdm:')
# subprocess.call(str_cmd_subprocess)
#
# path_to_emulators = path_to_cosmopower_organization + 'neff/'
# str_cmd_subprocess = ["ls",path_to_emulators]
# print('\ninside neff:')
# subprocess.call(str_cmd_subprocess)
