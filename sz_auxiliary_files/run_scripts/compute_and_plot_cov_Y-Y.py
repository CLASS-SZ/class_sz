#! python3
import argparse
import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
path_to_class = '/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/'
FIG_DIR = '/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/sz_auxiliary_files/run_scripts/figures'
r_dict = {} #dictionary of results


def run(args):
    os.chdir(path_to_class)

    #load template parameter file into dictionnary
    p_dict = {}
    with open(path_to_class+"sz_auxiliary_files/run_scripts/class-sz_parameters.ini") as f:
        for line in f:
            x = line.strip()
            if not x.startswith("#"):
                print(x.split('='))
                (key, val) = x.split('=')
                (key, val) = (key.strip(), val.strip())
                p_dict[key] = val
    f_sky = 0.5
    param_name = 'f_sky'
    p_dict[param_name] = f_sky
    param_name = 'output'
    p_dict[param_name] = 'tSZ_1h,tSZ_Trispectrum'
    #p_dict[param_name] = 'tSZ_1h'
    param_name = 'multipoles_sz'
    p_dict[param_name] = 'ell_mock'
    param_name = 'ell_min_mock'
    p_dict[param_name] = 2.
    param_name = 'ell_max_mock'
    p_dict[param_name] = 1.e4
    param_name = 'redshift_epsabs'
    p_dict[param_name] = 1.e-30
    param_name = 'mass_epsabs'
    p_dict[param_name] = 1.e-30
    p_dict['mass_epsrel'] = 1.e-4
    p_dict['redshift_epsrel'] = 1.e-4
    p_dict['include_ssc'] = 'no'

    table = []
    for i in np.arange(0,4):
        if i == 0:
            for j in np.arange(2,10):
                table.append(j*np.power(10,i))
        else :
            for j in np.arange(1,10):
                table.append(j*np.power(10,i))
    table.append(1e4)
    ell_ticks_table = np.asarray(table)

    ell_array_major = np.asarray([10,1e2,1e3,1e4])
    ell_array_label = [r'$10$',r'$10^{2}$',r'$10^{3}$',r'$10^{4}$']


    def func_where_is_ell(ellp):
        return np.log(ellp/float(p_dict['ell_min_mock']))*1e2/np.log(float(p_dict['ell_max_mock'])/float(p_dict['ell_min_mock']))


    ell_array_major = func_where_is_ell(ell_array_major)
    ell_ticks_table = func_where_is_ell(ell_ticks_table)


    #create a temporary ini file with parameter value
    subprocess.call(['mkdir',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])
    with open(path_to_class+'sz_auxiliary_files/run_scripts/tmp/tmp.ini', 'w') as f:
        for k, v in p_dict.items():
            f.write(str(k) + ' = '+ str(v) + '\n')
    subprocess.call(['./class',path_to_class+'sz_auxiliary_files/run_scripts/tmp/tmp.ini'])

    #store the results
    R = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum.txt')
    C = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum_r_cl_clp.txt')
    x_axis = R[:,0]
    #remove the temporary files
    subprocess.call(['rm','-r',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])

    SZ_ps_cov = C
    ell_array = x_axis
    bin_ell = np.arange(0,101)
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter

    label_size = 13
    title_size = 15
    legend_size = 15
    handle_length = 1.5

    fig = plt.figure(figsize=(8,6))

    ax = fig.add_subplot(111)
    fig.subplots_adjust(right=0.75)

    Y = bin_ell
    X = bin_ell

    Xm,Ym = np.meshgrid(X,Y)
    Z = np.log10(SZ_ps_cov)
    print(Z)

    surf = ax.imshow(Z, interpolation='nearest', origin='lower', aspect='auto',
            extent=[X[0], X[-1], Y[0], Y[-1]])

    ax.tick_params(axis = 'x',which='both',length=5,direction='out', pad=7)
    ax.tick_params(axis = 'y',which='both',length=5,direction='out', pad=5)


    y_label_list = ell_array_label
    ax.set_xticks(ell_ticks_table, minor=True)
    ax.set_xticks(ell_array_major)
    ax.set_xticklabels(y_label_list)

    y_label_list = ell_array_label
    ax.set_yticks(ell_ticks_table, minor=True)
    ax.set_yticks(ell_array_major)
    ax.set_yticklabels(y_label_list)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


    plt.setp(ax.get_yticklabels(), rotation=0,
             fontsize=label_size)
    plt.setp(ax.get_xticklabels(),  rotation=0,
             fontsize=label_size)

    ax.set_ylabel(r'$\mathrm{multipole\,}\,\,\ell$',size=title_size,labelpad = 5)
    ax.set_xlabel(r'$\mathrm{multipole\,}\,\,\ell^\prime$',size=title_size,labelpad = -1)


    cbaxes = fig.add_axes([0.81, 0.23, 0.05, 0.6])

    axcb = fig.colorbar(surf,cax=cbaxes,pad=0.)
    plt.subplots_adjust(wspace = .1)
    axcb.set_label(r'$\mathrm{log}_{10}[\frac{\mathrm{cov(C_\ell,C_{\ell^\prime})^{SSC}}}{\sqrt{\mathrm{cov(C_\ell,C_\ell)}\mathrm{cov(C_{\ell^\prime},C_{\ell^\prime})}}}]$',rotation=-90,labelpad = 40,y=0.5,size=15)

    #fig.tight_layout()
    FIG_NAME = '/cov_Y-Y'
    plt.savefig(FIG_DIR + FIG_NAME +".pdf")
    plt.show()



def main():
	parser=argparse.ArgumentParser(description="Plot cosmotherm spectra")

	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)

if __name__=="__main__":
	main()
