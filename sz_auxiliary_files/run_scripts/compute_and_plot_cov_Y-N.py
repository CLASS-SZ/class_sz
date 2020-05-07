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
    #collect arguments
    # param_name = args.param_name
    #
    # if (param_name == 'sigma8'):
    #     label_key = r'$\sigma_8$'
    # if (param_name == 'b'):
    #     param_name = 'HSEbias'
    #     label_key = r'$b$'
    # if (param_name == 'h'):
    #     label_key = r'$h$'
    # if (param_name == 'Omega_cdm'):
    #     label_key = r'$\Omega_\mathrm{m}$'
    #
    # else:
    #     label_key = 'val'


    # p_min = float(args.p_min)
    # p_max = float(args.p_max)
    # N = int(args.N)
    # spacing = args.spacing

    #define array of parameter values
    # if(spacing=='log'):
    #     p = np.logspace(np.log10(p_min),np.log10(p_max),N)
    # elif(spacing=='lin'):
    #     p = np.linspace(p_min,p_max,N)

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
    p_dict[param_name] = 'tSZ_cov_Y_N'

    param_name = 'number of mass bins for cov(Y,N)'
    p_dict[param_name] = float(args.nbins_mass)


    nbins_M = float(args.nbins_mass)

    M_bins = []
    M1SZ = float(p_dict['M1SZ'])
    M2SZ = float(p_dict['M2SZ'])

    for i in range(0,50):
        M_bins.append(np.exp(np.log(M1SZ)+i*(np.log(M2SZ)-np.log(M1SZ))/(nbins_M-1.)))
    M_bins = np.asarray(M_bins)
    M_bins

    def func_where_is_mp(mp):
        return np.log(mp/M1SZ)*1e2/np.log(M2SZ/M1SZ)

    mp_array_major = np.asarray([1e11,1e12,1e13,1e14,1e15])
    mp_array_label_major = [r'$10^{11}$',r'$10^{12}$',r'$10^{13}$',r'$10^{14}$',r'$10^{15}$']
    mp_array_major = func_where_is_mp(mp_array_major)
    one_to_nine = np.arange(1,10)
    one_to_five = np.arange(1,6)
    thirteen_to_fifteen = np.arange(13,16)
    table = []
    for i in thirteen_to_fifteen:
        if i != 15:
            for j in one_to_nine:
                table.append(j*np.power(10,i))
        else :
            for j in one_to_five:
                table.append(j*np.power(10,i))

    m_ticks_table = np.asarray(table)
    m_ticks_table = func_where_is_mp(m_ticks_table)


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
    subprocess.call(['rm','-rf',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])
    subprocess.call(['mkdir',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])
    with open(path_to_class+'sz_auxiliary_files/run_scripts/tmp/tmp.ini', 'w') as f:
        for k, v in p_dict.items():
            f.write(str(k) + ' = '+ str(v) + '\n')
    subprocess.call(['./class',path_to_class+'sz_auxiliary_files/run_scripts/tmp/tmp.ini'])

    #store the results
    R = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum.txt')
    C = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum_r_Y_N.txt')
    x_axis = R[:,0]
    # y_axis = R[:,7]
    # L = [x_axis,y_axis]
    # r_dict[p_val] = L
    #remove the temporary files
    #subprocess.call(['rm','-r',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])


    #print(p_dict)
    SZ_ps_cov = C
    ell_array = x_axis
    bin_M = np.arange(0,101)
    bin_ell = np.arange(0,101)
    from mpl_toolkits.mplot3d import Axes3D
    #import matplotlib.pyplot as plt
    #from matplotlib import cm
    #from matplotlib.ticker import LinearLocator, FormatStrFormatter
    from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter

    #import numpy as np
    label_size = 13
    title_size = 15
    legend_size = 15
    handle_length = 1.5

    #pylab.rc('text', usetex = False)
    fig = plt.figure(figsize=(8,6))
    #ax = Axes3D(fig)
    #ax = fig.gca(projection='3d')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(right=0.75)
    #plt.subplots_adjust(left=0, bottom=0, right=10,top =10)




    Y = bin_ell
    X = bin_M
    #X = p_array
    #Y = w_array


    Xm,Ym = np.meshgrid(X,Y)
    Z = np.log10(SZ_ps_cov)
    #print(Z)

    #extent = np.min(Xm), np.max(Xm), np.min(Ym), np.max(Ym)
    #surf = ax.imshow(Z, cmap='viridis', interpolation='bilinear')

    #ax.plot_surface(Ym,Xm,Z,color='blue',alpha=0.5)

    #surf=ax.contourf(X, Y, Z, 50, cmap='viridis',interpolation='bilinear')

    #surf = plt.pcolormesh(Xm,Ym,Z, cmap='viridis',shading='gouraud')
    surf = ax.imshow(Z, interpolation='nearest', origin='lower', aspect='auto',
            extent=[X[0], X[-1], Y[0], Y[-1]])

    #ax.colorbar()
    #ax.clabel(surf, inline=1, fontsize=10)
    #ax.clabel(CS, inline=0, fontsize=10,fmt=fmt)



    #ax.zaxis.set_rotate_label(True)
    #ax.set_yscale('log')

    ax.tick_params(axis = 'x',which='both',length=5,direction='out', pad=7)
    ax.tick_params(axis = 'y',which='both',length=5,direction='out', pad=5)


    #x_label_list = mp_array_label
    ax.set_xticks(m_ticks_table, minor=True)
    ax.set_xticks(mp_array_major)
    ax.set_xticklabels(mp_array_label_major)


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
    ax.set_xlabel('$\mathrm{mass\,\,\, M}\quad\,\,[\mathrm{M_{\odot}}/h]$',size=title_size,labelpad = -1)



    #sel='_by_bh'
    #ax.view_init(55,65)
    #ax.set_xlim3d(0.98*np.min(Ym), np.max(Ym))
    #ax.set_ylim(0.98*np.min(Ym),3.)
    #ax.set_xlim(2e11, 3e14)


    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.)





    cbaxes = fig.add_axes([0.81, 0.23, 0.05, 0.6])
    #cax = plt.axes([0.9, 0.1, 0.075, 0.8])
    #cb = fix.colorbar(ax1, cax = cbaxes)

    # axins = inset_axes(ax,
    #                width="5%",  # width = 5% of parent_bbox width
    #                height="50%",  # height : 50%
    #                loc='lower left',
    #                bbox_to_anchor=(1.05, 0., 1, 1),
    #                bbox_transform=ax.transAxes,
    #                borderpad=2,
    #                )
    axcb = fig.colorbar(surf,cax=cbaxes,pad=0.)
    plt.subplots_adjust(wspace = .1)
    #axcb = fig.colorbar(lc,cax=axins)
    #axcb.set_clim(-8.5, 0)
    axcb.set_label(r'$\mathrm{log}_{10}[\frac{\mathrm{cov(N,C_\ell)}}{\mathrm{cov(N,N)}\mathrm{cov(C_\ell,C_\ell)}}]$',rotation=-90,labelpad = 40,y=0.5,size=15)

    #fig.subplots_adjust(left=0.01,right=0.1)


    #fig.tight_layout()
    FIG_NAME = '/cov_Y-N'
    #plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    plt.savefig(FIG_DIR + FIG_NAME +".pdf")
    plt.show()



def main():
	parser=argparse.ArgumentParser(description="Plot covariance matrix Y-N")
	parser.add_argument("-nbins_mass",help="number of mass bins" ,dest="nbins_mass", type=str, required=True)
	#parser.add_argument("-param_name",help="name of parameter" ,dest="param_name", type=str, required=False)
	#parser.add_argument("-param_value",help="value of parameter" ,dest="param_value", type=str, required=False)
	# parser.add_argument("-min",help="minimum value of parameter" ,dest="p_min", type=str, required=True)
	# parser.add_argument("-max",help="maximum value of parameter" ,dest="p_max", type=str, required=True)
	# parser.add_argument("-N",help="number of evaluations" ,dest="N", type=int, required=True)
	# parser.add_argument("-spacing",help="linear (lin) or log spacing (log)" ,dest="spacing", type=str, required=True)
	# parser.add_argument("-show_legend",help="show legend on figure? ('yes' or 'no')" ,dest="show_legend", type=str, required=True)
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)

if __name__=="__main__":
	main()
