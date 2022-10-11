# To compute cib x cib:
# $ python compute_and_plot_cib_simple.py -param_name 'h' -p_val '[0.6711]' -output 'cib_cib_1h,cib_cib_2h' -plot_cib_cib yes -save_figure yes -mode run
# To compute tSZ x CIB
# $ python compute_and_plot_cib_simple.py -param_name 'h' -p_val '[0.6711]' -output 'tSZ_cib_1h,tSZ_cib_2h' -plot_tSZ_cib yes -save_figure yes -mode run
# To compute lens x CIB
# $ python compute_and_plot_cib_simple.py -param_name 'h' -p_val '[0.6711]' -output 'lens_cib_1h,lens_cib_2h' -plot_lens_cib yes -save_figure yes -mode run

# Note the last option "-mode run" if you want to run class_sz, and "-mode plot" if you just want to plot
# what was last computed.
import matplotlib.pyplot as plt
import argparse
import numpy as np
import os
import subprocess

from matplotlib.pyplot import cm
from scipy import interpolate
from scipy import stats
from datetime import datetime
import ast
import itertools

import pyccl as ccl

# do_ccl_comparison = 'yes'
do_ccl_comparison = 'no'
# path_to_ccl = '/Users/boris/Work/CCL/'
freq_cib_1 = 545.
freq_cib_2 = 545.
# table 1 of https://arxiv.org/pdf/1309.0382.pdf
#1: freq GHz 2: Flux cut mJy
# 100 - 400
# 143 - 350
# 217 - 225
# 353 - 315
# 545 - 350
# 857 - 710
# 3000  - 1000
cib_Snu_1 = 710.
cib_Snu_2 = 710.



# set path to the class_sz code
path_to_class = '/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/'
# path_to_class = '/Users/boris/Downloads/class_sz-df37ebccd13a2e266cecd9bf553b729e93101e7e/'
# set path to to the repository containing the folder 'cib_files' with your cib's in it
# set path to where to save the figures
FIG_DIR = '/Users/boris/Work/CLASS-SZ/SO-SZ/figures'
subprocess.call(['mkdir','-p',path_to_class+'sz_auxiliary_files/run_scripts'])

# usefule function for formatting numbers
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
    if p_value<1.:
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

def run(args):
    # important parameters are re-ajusted later, here we just load a template file:
    # parameter_file = 'class-sz_parameters_MM20.ini'

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
    elif (param_name == 'M2SZ'):
        label_key = r'$M_\mathrm{max}$'
    elif (param_name == 'z2SZ'):
        label_key = r'$z_\mathrm{max}$'
    elif (param_name == 'M1SZ'):
        label_key = r'$M_\mathrm{min}$'
    elif (param_name == 'm_ncdm'):
        label_key = r'$\Sigma m_\mathrm{\nu}$'
    elif (param_name == 'signal-to-noise cut-off for ps completeness analysis'):
        label_key = r'$q_\mathrm{cut}$'
    else:
        label_key = 'val'


    if args.p_val:
        try:
            p = np.asarray(ast.literal_eval(args.p_val))
        except ValueError:
            s = str(args.p_val)

            s = s.replace("[", "")
            s = s.replace("]", "")
            s = s.split(",")
            print(s[0])
            p = s

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

    #array to store results
    cl_1h = []
    trispectrum = []
    multipoles = []
    col = []
    val_label = []
    y_err = []
    cl_2h = []
    te_y_y = []
    kSZ_kSZ_gal_1halo = []
    tSZ_tSZ_tSZ_1h = []
    tSZ_gal_1h = []
    tSZ_gal_2h = []
    gal_gal_1h = []
    gal_gal_2h = []
    gal_lens_1h = []
    gal_lens_2h = []
    tSZ_lens_1h = []
    tSZ_lens_2h = []
    isw_lens = []
    isw_tsz = []
    isw_auto = []
    lens_lens_1h = []
    lens_lens_2h = []
    cib_cib_1h = []
    cib_cib_2h = []
    tSZ_cib_1h = []
    tSZ_cib_2h = []
    lens_cib_1h = []
    lens_cib_2h = []


    #build parameter file into dictionnary
    p_dict = {}


    # 'M500' is Tinker et al 2008 @ M500
    # 'T10' is Tinker et al 2010 @ m200_mean

    p_dict['mass function'] = 'T08M200c'
    p_dict['concentration parameter'] = 'D08'
    p_dict['delta for cib'] = '200m'
    p_dict['hm_consistency'] = 1
    p_dict['damping_1h_term'] = 0




    # parameters for Cosmology
    p_dict['Omega_cdm'] = 0.3175-0.022068/0.6711/0.6711
    p_dict['omega_b'] = 0.022068
    p_dict['h'] = 0.6711
    p_dict['A_s'] = 2.2e-9
    p_dict['n_s'] = .9624
    p_dict['k_pivot'] = 0.05


    p_dict['N_ncdm'] = 1
    p_dict['N_ur'] = 0.00641
    p_dict['deg_ncdm'] = 3
    p_dict['m_ncdm'] = 0.02
    p_dict['T_ncdm'] = 0.71611


    # mass bounds
    p_dict['M_min'] = 1e10*p_dict['h']
    p_dict['M_max'] = 1e16*p_dict['h']

    # redshift bounds
    p_dict['z_min'] = 1e-5
    p_dict['z_max'] = 6. # fiducial for MM20 : 6

    p_dict['freq_min'] = 10.
    p_dict['freq_max'] = 5e4 # fiducial for MM20 : 6
    p_dict['n_nu_L_sat'] = 200

    # HOD parameters for CIB
    p_dict['M_min_HOD'] = pow(10.,10)
    # p_dict['M1_prime_HOD'] =pow(10.,125.1536196)*p_dict['h']





    # CIB parametes see McCarthy & Madhavacheril 2020
    p_dict['Redshift evolution of dust temperature'] =  0.36
    p_dict['Dust temperature today in Kelvins'] = 24.4
    p_dict['Emissivity index of sed'] = 1.75
    p_dict['Power law index of SED at high frequency'] = 1.7
    p_dict['Redshift evolution of L − M normalisation'] = 3.6
    p_dict['Most efficient halo mass in Msun'] = pow(10.,12.6)
    p_dict['Normalisation of L − M relation in [Jy MPc2/Msun]'] = 6.4e-8
    p_dict['Size of of halo masses sourcing CIB emission'] = 0.5
    p_dict['has_cib_flux_cut'] = 0

    # List of frequency bands for cib
    if freq_cib_1 == freq_cib_2:
        p_dict['cib_frequency_list_num'] = 1
        p_dict['cib_frequency_list_in_GHz'] = str(freq_cib_1)
        p_dict['cib_Snu_cutoff_list [mJy]'] = str(cib_Snu_1)
        # p_dict['cib_frequency_list_num'] = 5
        # p_dict['cib_frequency_list_in_GHz'] = '217,353,545,857,3000'
        p_dict["Frequency_id nu for cib in GHz (to save in file)"] = 0
        p_dict["Frequency_id nu^prime for cib in GHz (to save in file)"] = 0

    else:
        p_dict['cib_frequency_list_num'] = 2
        p_dict['cib_frequency_list_in_GHz'] = str(freq_cib_1)+','+str(freq_cib_2)
        p_dict['cib_Snu_cutoff_list [mJy]'] = str(cib_Snu_1)+','+str(cib_Snu_2)
        # p_dict['cib_frequency_list_num'] = 5
        # p_dict['cib_frequency_list_in_GHz'] = '217,353,545,857,3000'
        p_dict["Frequency_id nu for cib in GHz (to save in file)"] = 0
        p_dict["Frequency_id nu^prime for cib in GHz (to save in file)"] = 1
    # verbose paramete of class_sz
    p_dict['class_sz_verbose'] = 10
    p_dict['root'] = 'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_'
    p_dict['write sz results to files'] = 'yes' # this writes  PS and f(z)
    # precision parameters
    p_dict['pressure_profile_epsabs'] = 1.e-8
    p_dict['pressure_profile_epsrel'] = 1.e-3
    # precision for redshift integal
    p_dict['redshift_epsabs'] = 1e-40#1.e-40
    p_dict['redshift_epsrel'] = 1e-4#1.e-10 # fiducial value 1e-8
    # precision for mass integal
    p_dict['mass_epsabs'] = 1e-40 #1.e-40
    p_dict['mass_epsrel'] = 1e-4#1e-10
    # precision for Luminosity integral (sub-halo mass function)
    p_dict['L_sat_epsabs'] = 1e-40 #1.e-40
    p_dict['L_sat_epsrel'] = 1e-3#1e-10


    p_dict['z_max_pk'] = p_dict['z_max']

    # multipole array
    p_dict['dlogell'] = 0.4
    p_dict['ell_max'] = 5000.
    p_dict['ell_min'] = 2.



    p_dict['pressure profile'] = 'Custom. GNFW'
    p_dict['P0GNFW']  = 6.41
    p_dict['c500']  = 1.81
    p_dict['gammaGNFW']  = 0.31
    p_dict['alphaGNFW']  = 1.33
    p_dict['betaGNFW']  = 4.13
    p_dict['B'] = 1.25
    p_dict['x_outSZ'] = 4.



    p_dict['output'] = args.output
    p_dict['root'] = 'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_'
    # p_dict['path_to_class'] = path_to_class


    L_ref = np.loadtxt(path_to_class + 'sz_auxiliary_files/cib_files/cib_1h_'+str(int(freq_cib_1))+'x'+str(int(freq_cib_2))+'.txt')
    ell_MM20 = L_ref[:,0]
    cl_cib_cib_1h_MM20 = L_ref[:,1]
    L_ref = np.loadtxt(path_to_class + 'sz_auxiliary_files/cib_files/cib_2h_'+str(int(freq_cib_1))+'x'+str(int(freq_cib_2))+'.txt')
    cl_cib_cib_2h_MM20 = L_ref[:,1]


    if do_ccl_comparison == 'yes':
        # Initialize cosmology and halo model ingredients
        # bm = np.loadtxt(path_to_ccl + "/benchmarks/data/cib_class_sz_szpowerspectrum.txt",
        #          unpack=True)
        cosmo = ccl.Cosmology(Omega_b=p_dict['omega_b']/p_dict['h']**2.,
                           Omega_c=p_dict['Omega_cdm'],
                           h=p_dict['h'],
                           n_s=p_dict['n_s'],
                           A_s=p_dict['A_s'],
                           Neff=3.046)
        mdef = ccl.halos.MassDef200m()
        cM = ccl.halos.ConcentrationDuffy08(mdef)
        nM = ccl.halos.MassFuncTinker10(cosmo, mdef, norm_all_z=True)
        bM = ccl.halos.HaloBiasTinker10(cosmo, mdef)
        hmc = ccl.halos.HMCalculator(cosmo, nM, bM, mdef)
        if freq_cib_1 == freq_cib_2:
            pr = ccl.halos.HaloProfileCIBShang12(cM, freq_cib_1, Mmin=1E10)
            pr.update_parameters(nu_GHz=freq_cib_1,
                              alpha=0.36,
                              T0=24.4,
                              beta=1.75,
                              gamma=1.7,
                              s_z=3.6,
                              log10meff=12.6,
                              sigLM=np.sqrt(0.5),
                              Mmin=1E10,
                              L0=6.4E-8)

            pr2pt = ccl.halos.Profile2ptCIB()

            # CIB tracer
            tr = ccl.CIBTracer(cosmo, z_min=p_dict['z_min'])

            # 3D power spectrum
            pk = ccl.halos.halomod_Pk2D(cosmo, hmc, pr, prof_2pt=pr2pt)

        # Angular power spectrum
        # ls = bm[0]
        ls = np.geomspace(p_dict['ell_min'],p_dict['ell_max'],100)

        cl = ccl.angular_cl(cosmo, tr, tr, ls, p_of_k_a=pk)


        dl = cl*ls*(ls+1)/(2*np.pi)
        l_ccl = ls
        dl_ccl = dl






    #prepare the figure
    label_size = 12
    title_size = 15
    legend_size = 13
    handle_length = 1.5
    fig, ax1 = plt.subplots(1,1,figsize=(7,5))
    ax = ax1
    ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
    ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=label_size)
    plt.setp(ax.get_xticklabels(), fontsize=label_size)
    ax.grid( b=True, which="both", alpha=0.3, linestyle='--')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\ell$',size=title_size)
    if (args.plot_tSZ_cib == 'yes'):
            ax.set_ylabel(r'$\ell(\ell+1)\mathrm{C^{tSZ\times cib}_\ell/2\pi}$',size=title_size)
    elif (args.plot_lens_cib == 'yes'):
            ax.set_ylabel(r'$\ell(\ell+1)\mathrm{C^{lens\times cib}_\ell/2\pi}$',size=title_size)
    elif (args.plot_cib_cib == 'yes'):
        ax.set_ylabel(r'$\ell(\ell+1)\mathrm{C^{cib\times cib}_\ell/2\pi}$',size=title_size)
    else:
        print("you need to specify a quantity to plot, e.g., '-plot_cib_cib yes'")

    if(args.x_min):
        ax.set_xlim(left=float(args.x_min))
    if(args.x_max):
        ax.set_xlim(right=float(args.x_max))

    if(args.y_min):
        ax.set_ylim(bottom=float(args.y_min))
    if(args.y_max):
        ax.set_ylim(top=float(args.y_max))
    if do_ccl_comparison == 'yes':
        ax.plot(l_ccl,dl_ccl,label='ccl')


    colors = iter(cm.viridis(np.linspace(0, 1, N)))

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
                        f.write(str(k) + ' = '+ str(v) + '\n')
                startTime = datetime.now()
                subprocess.call(['./class',path_to_class+'sz_auxiliary_files/run_scripts/tmp/tmp.ini'])
                print("time in class -> " + str((datetime.now() - startTime)))
            #print(datetime.now() - startTime)

            # try:
            # val_label.append(label_key + ' = ' + scientific_notation(p_val))
                #val_label.append(label_key + ' = ' + "%.2f"%(p_val))

            col.append(next(colors))

            if args.mode == 'run':
                R = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum.txt')
            elif args.mode == 'plot':
                R = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_szpowerspectrum_' + str(p_val) + '.txt')

            multipoles.append(R[:,0])
            cl_1h.append(R[:,1])
            cl_2h.append(R[:,6])
            #print(cl_2h)
            trispectrum.append(R[:,3])
            te_y_y.append(R[:,7])
            kSZ_kSZ_gal_1halo.append(R[:,8])
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
            lens_lens_2h.append(R[:,28])
            tSZ_cib_1h.append(R[:,29])
            tSZ_cib_2h.append(R[:,30])
            cib_cib_1h.append(R[:,31])
            cib_cib_2h.append(R[:,32])
            lens_cib_1h.append(R[:,34])
            lens_cib_2h.append(R[:,35])


            if (args.plot_cib_cib == 'yes'):
                print('plotting cibxcib')
                ax.plot(ell_MM20,cl_cib_cib_1h_MM20,label='MM20-1h')
                ax.plot(ell_MM20,cl_cib_cib_2h_MM20,ls='-',label='MM20-2h')
                print(cib_cib_1h[id_p])
                print(cib_cib_2h[id_p])
                # print('ok')
                print(multipoles[id_p])
                ax.plot(multipoles[id_p],(cib_cib_2h[id_p]),color='r',ls='--',alpha = 1.,
                marker =  'o',markersize = 2,
                label = 'class_sz (2-halo)')
                # print('ok')
                ax.plot(multipoles[id_p],(cib_cib_1h[id_p]),color='k',ls='--',alpha = 1.,
                marker =  '*',markersize = 1,
                label = 'class_sz (1-halo)')
                # print('ok')
                ax.plot(multipoles[id_p],(cib_cib_1h[id_p]+cib_cib_2h[id_p]),color='grey',ls='-',alpha = 1.,
                marker =  '*',markersize = 1,
                label = 'class_sz (1+2-halo)')

            elif (args.plot_tSZ_cib == 'yes'):
                print(tSZ_cib_1h[id_p])
                print(tSZ_cib_2h[id_p])
                ax.plot(multipoles[id_p],(tSZ_cib_2h[id_p]),color='r',ls='--',alpha = 1.,
                marker =  'o',markersize = 2,
                label = 'class_sz (2-halo)')
                ax.plot(multipoles[id_p],(tSZ_cib_1h[id_p]),color='k',ls='--',alpha = 1.,
                marker =  '*',markersize = 1,
                label = 'class_sz (1-halo)')
                ax.plot(multipoles[id_p],(tSZ_cib_1h[id_p]+tSZ_cib_2h[id_p]),color='grey',ls='-',alpha = 1.,
                marker =  '*',markersize = 1,
                label = 'class_sz (1+2-halo)')

            elif (args.plot_lens_cib == 'yes'):
                print(lens_cib_1h[id_p])
                print(lens_cib_2h[id_p])
                ax.plot(multipoles[id_p],multipoles[id_p]*(lens_cib_2h[id_p]),color='r',ls='--',alpha = 1.,
                marker =  'o',markersize = 2,
                label = 'class_sz (2-halo)')
                ax.plot(multipoles[id_p],multipoles[id_p]*(lens_cib_1h[id_p]),color='k',ls='--',alpha = 1.,
                marker =  '*',markersize = 1,
                label = 'class_sz (1-halo)')
                ax.plot(multipoles[id_p],multipoles[id_p]*(lens_cib_1h[id_p]+lens_cib_2h[id_p]),color='grey',ls='-',alpha = 1.,
                marker =  '*',markersize = 1,
                label = 'class_sz hod (1+2-halo)')
            # print('ok3')
            ax.legend(loc=4,ncol = 1)
            # print('ok4')
            plt.draw()
            # print('ok5')
            plt.pause(0.05)
            # print('ok6')
            if args.mode == 'run':
                #save some results and remove the temporary files
                # print('ok3')
                subprocess.call(['mv',path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum.txt', path_to_class+'sz_auxiliary_files/run_scripts/tmp/' + 'class-sz_szpowerspectrum_' + str(p_val) + '.txt'])
                # print('ok3')
            id_p += 1


        #end loop on p over param values

    fig.tight_layout()

    if (args.save_fig == 'yes'):
        FIG_NAME = '/my_figure'
        plt.savefig(FIG_DIR + FIG_NAME +".pdf")

    plt.show(block=True)







def main():
    parser=argparse.ArgumentParser(description="Plot cosmotherm spectra")
    parser.add_argument("-mode",help="plot or run" ,dest="mode", type=str, required=True)
    parser.add_argument("-param_name",help="name of varying parameter" ,dest="param_name", type=str, required=True)
    parser.add_argument("-min",help="minimum value of parameter" ,dest="p_min", type=str, required=False)
    parser.add_argument("-max",help="maximum value of parameter" ,dest="p_max", type=str, required=False)
    parser.add_argument("-N",help="number of evaluations" ,dest="N", type=int, required=False)
    parser.add_argument("-p_val",help="list of param values" ,dest="p_val", type=str, required=False)
    parser.add_argument("-spacing",help="linear (lin) or log spacing (log)" ,dest="spacing", type=str, required=False)
    parser.add_argument("-y_min",help="ylim for y-axis" ,dest="y_min", type=str, required=False)
    parser.add_argument("-y_max",help="ylim for y-axis" ,dest="y_max", type=str, required=False)
    parser.add_argument("-x_min",help="xlim for x-axis" ,dest="x_min", type=str, required=False)
    parser.add_argument("-x_max",help="xlim for x-axis" ,dest="x_max", type=str, required=False)
    parser.add_argument("-output",help="what quantities to compute" ,dest="output", type=str, required=False)
    parser.add_argument("-save_figure",help="save figure" ,dest="save_fig", type=str, required=False)
    parser.add_argument("-plot_cib_cib",help="cib_cib" ,dest="plot_cib_cib", type=str, required=False)
    parser.add_argument("-plot_tSZ_cib",help="tSZ_cib" ,dest="plot_tSZ_cib", type=str, required=False)
    parser.add_argument("-plot_lens_cib",help="lens_cib" ,dest="plot_lens_cib", type=str, required=False)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
	main()
