#! python3
import argparse
import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
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



    param_name = 'output'
    #p_dict[param_name] = 'tSZ_cov_Y_N'
    #p_dict[param_name] = 'tSZ_cov_hsv'
    p_dict[param_name] = 'tSZ_cov_N_N'

    f_sky = 0.5
    param_name = 'f_sky'
    p_dict[param_name] = f_sky
    param_name = 'redshift_epsabs'
    p_dict[param_name] = 1.e-30
    param_name = 'mass_epsabs'
    p_dict[param_name] = 1.e-30
    param_name = 'M1SZ'
    p_dict[param_name] = 1.e11
    param_name = 'M2SZ'
    p_dict[param_name] = 1.e16
    p_dict['include_ssc'] = 'yes'
    param_name = 'number of mass bins for cov(Y,N)'
    p_dict[param_name] = float(args.nbins_mass)


    nbins_M = float(args.nbins_mass)

    # M_bins = []
    M1SZ = float(p_dict['M1SZ'])
    M2SZ = float(p_dict['M2SZ'])


    #create a temporary ini file with parameter value
    subprocess.call(['rm','-rf',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])
    subprocess.call(['mkdir',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])
    with open(path_to_class+'sz_auxiliary_files/run_scripts/tmp/tmp.ini', 'w') as f:
        for k, v in p_dict.items():
            f.write(str(k) + ' = '+ str(v) + '\n')
    subprocess.call(['./class',path_to_class+'sz_auxiliary_files/run_scripts/tmp/tmp.ini'])

    #store the results

    C = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum_cov_N_N_diagonal.txt')
    m_bins_lower_edges = np.asarray(C[:,0])
    m_bins_upper_edges = np.asarray(C[:,1])
    m_bins_center = np.exp(0.5*(np.log(m_bins_lower_edges)+np.log(m_bins_upper_edges)))
    x_axis = m_bins_center

    N_counts = np.asarray(C[:,2])
    shot_error = np.sqrt(np.asarray(C[:,2]))/N_counts
    sample_error = np.sqrt(np.asarray(C[:,3]))/N_counts


    #prepare the figure
    label_size = 12
    title_size = 15
    legend_size = 13
    handle_length = 1.5
    fig, (ax1) = plt.subplots(1,1,figsize=(7,5))
    ax = ax1
    ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
    ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_yticklabels(), rotation='vertical', fontsize=label_size)
    plt.setp(ax.get_xticklabels(), fontsize=label_size)
    plt.grid( b=True, which="both", alpha=0.3, linestyle='--')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(r'$M\,\, [\mathrm{M_\odot}/h]$',size=title_size)
    ax.set_ylabel(r'$\sigma_N/N$',size=title_size)

    y_axis = shot_error
    ax.plot(x_axis,y_axis,color='k',ls='--',alpha = 1.,label = "Poisson")
    y_axis = sample_error
    ax.plot(x_axis,y_axis,color='k',ls='-',alpha = 1.,label = "SSC")

    ax2 = ax1.twinx()
    ax = ax2
    y_axis = N_counts
    ax.plot(x_axis,y_axis,color='r',ls='-',alpha = 1.,label = "counts",marker = 'o',markerfacecolor='None')
    ax.set_ylabel(r'$N$',size=title_size,color='r')
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax1.legend(loc=1)
    plt.draw()

    fig.tight_layout()
    FIG_NAME = '/cov_N-N_diagonal'
    fig.savefig(FIG_DIR + FIG_NAME +".pdf")
    plt.show(block=False)
    #plt.close()

    M1SZ = float(p_dict['M1SZ'])
    M2SZ = float(p_dict['M2SZ'])

    # for i in range(0,50):
    #     M_bins.append(np.exp(np.log(M1SZ)+i*(np.log(M2SZ)-np.log(M1SZ))/(nbins_M-1.)))
    # M_bins = np.asarray(M_bins)
    # M_bins

    def func_where_is_mp(mp):
        return np.log(mp/M1SZ)*1e2/np.log(M2SZ/M1SZ)

    m_array_major = np.asarray([1e11,1e12,1e13,1e14,1e15,1e16])
    m_array_label = [r'$10^{11}$',r'$10^{12}$',r'$10^{13}$',r'$10^{14}$',r'$10^{15}$',r'$10^{16}$']
    m_array_major = func_where_is_mp(m_array_major)
    one_to_nine = np.arange(1,10)
    one_to_five = np.arange(1,6)
    thirteen_to_fifteen = np.arange(11,16)
    table = []
    for i in thirteen_to_fifteen:
        if i != 16:
            for j in one_to_nine:
                table.append(j*np.power(10,i))
        else :
            for j in one_to_five:
                table.append(j*np.power(10,i))

    m_ticks_table = np.asarray(table)
    m_ticks_table = func_where_is_mp(m_ticks_table)

    C = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum_r_N_N.txt')
    x_axis = m_bins_center
    #remove the temporary files
    #subprocess.call(['rm','-r',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])

    SZ_ps_cov = C
    m_array = x_axis
    bin_m = np.arange(0,101)
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter

    label_size = 13
    title_size = 15
    legend_size = 15
    handle_length = 1.5

    fig2 = plt.figure(figsize=(8,6))
    fig = fig2

    ax = fig.add_subplot(111)
    fig.subplots_adjust(right=0.75)

    Y = bin_m
    X = bin_m

    Xm,Ym = np.meshgrid(X,Y)
    Z = np.log10(SZ_ps_cov)
    print(Z)

    if (p_dict['include_ssc']!='yes'):
        colors=['yellow']
        cmap = mpl.colors.ListedColormap(colors)
    else:
        cmap = plt.get_cmap('viridis')

    surf = ax.imshow(Z, interpolation='nearest', origin='lower', aspect='auto',
            extent=[X[0], X[-1], Y[0], Y[-1]],cmap=cmap)

    ax.tick_params(axis = 'x',which='both',length=5,direction='out', pad=7)
    ax.tick_params(axis = 'y',which='both',length=5,direction='out', pad=5)


    y_label_list = m_array_label
    ax.set_xticks(m_ticks_table, minor=True)
    ax.set_xticks(m_array_major)
    ax.set_xticklabels(y_label_list)

    y_label_list = m_array_label
    ax.set_yticks(m_ticks_table, minor=True)
    ax.set_yticks(m_array_major)
    ax.set_yticklabels(y_label_list)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


    plt.setp(ax.get_yticklabels(), rotation=0,
             fontsize=label_size)
    plt.setp(ax.get_xticklabels(),  rotation=0,
             fontsize=label_size)

    ax.set_xlabel(r'$\mathrm{mass\,\,\, M}\quad\,\,[\mathrm{M_{\odot}}/h]$',size=title_size,labelpad = 5)
    ax.set_ylabel(r'$\mathrm{mass\,\,\, M}\quad\,\,[\mathrm{M_{\odot}}/h]$',size=title_size,labelpad = -1)




    if (p_dict['include_ssc']=='yes'):
        cbaxes = fig.add_axes([0.81, 0.23, 0.05, 0.6])
        axcb = fig.colorbar(surf,cax=cbaxes,pad=0.)
        plt.subplots_adjust(wspace = .1)
        axcb.set_label(r'$\mathrm{log}_{10}[\frac{\mathrm{cov(N_i,N_j)}}{\sqrt{\mathrm{cov(N_i,N_i)}\mathrm{cov(N_j,N_j)}}}]$',rotation=-90,labelpad = 40,y=0.5,size=15)

    #fig.tight_layout()
    FIG_NAME = '/cov_N-N'
    fig.savefig(FIG_DIR + FIG_NAME +".pdf")
    plt.show()


def main():
	parser=argparse.ArgumentParser(description="Plot covariance matrix N-N")
	parser.add_argument("-nbins_mass",help="number of mass bins" ,dest="nbins_mass", type=str, required=True)
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)

if __name__=="__main__":
	main()
