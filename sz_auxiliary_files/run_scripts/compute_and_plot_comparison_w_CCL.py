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
import pyccl as ccl


def run(args):
    os.chdir(path_to_class)
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
    p_dict[param_name] = 'tSZ_1h,dndlnM'

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
    param_name = 'sz_verbose'
    p_dict[param_name]  = 3
    param_name = 'm_ncdm'
    p_dict[param_name]  = 0.
    param_name = 'HMF_prescription_NCDM'
    p_dict[param_name]  = 'No-pres'
    param_name = 'B'
    p_dict[param_name]  = 1.41
    param_name = 'pressure profile'
    p_dict[param_name]  = 'Custom. GNFW'
    p_dict['P0GNFW'] = 6.41
    p_dict['c500'] = 1.81
    p_dict['gammaGNFW'] = 0.31
    p_dict['alphaGNFW'] = 1.33
    p_dict['betaGNFW'] = 4.13
    p_dict['pressure_profile_epsabs'] = 1.e-9
    p_dict['pressure_profile_epsrel'] = 1.e-2
    p_dict['ell_max_mock'] = 15000.
    p_dict['ell_min_mock'] = 2.

    M1SZ = float(p_dict['M1SZ'])
    M2SZ = float(p_dict['M2SZ'])


    #create a temporary ini file with parameter value
    subprocess.call(['rm','-rf',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])
    subprocess.call(['mkdir',path_to_class+'sz_auxiliary_files/run_scripts/tmp'])
    with open(path_to_class+'sz_auxiliary_files/run_scripts/tmp/tmp.ini', 'w') as f:
        for k, v in p_dict.items():
            f.write(str(k) + ' = '+ str(v) + '\n')
    subprocess.call(['./class',path_to_class+'sz_auxiliary_files/run_scripts/tmp/tmp.ini'])

    C = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum_dndlnM_redshifts.txt')
    dndlnM_redshifts = C
    C = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum_dndlnM_masses.txt')
    dndlnM_masses = C

    dndlnM = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum_dndlnM.txt')
    cl_1h = []
    R = np.loadtxt(path_to_class+'sz_auxiliary_files/run_scripts/tmp/class-sz_tmp_szpowerspectrum.txt')
    multipoles = R[:,0]
    cl_1h = R[:,1]
    x_axis_class_sz = dndlnM_masses
    y_axis_class_sz = dndlnM[:,0]



    cosmo = ccl.Cosmology(
        Omega_b=float(p_dict['Omega_b']),
        Omega_c=float(p_dict['Omega_cdm']),
        h=float(p_dict['h']),
        n_s=float(p_dict['n_s']),
        A_s=float(p_dict['A_s']),
        w0=-1,
        wa=0,
        m_nu=float(p_dict['m_ncdm']),
        m_nu_type='normal',
        Neff=3.046,
        Omega_k=0)
    mass_def = ccl.halos.MassDef(500, 'critical')
    hmf = ccl.halos.MassFuncTinker08(cosmo, mass_def=mass_def)
    hbf = ccl.halos.HaloBiasTinker10(cosmo, mass_def=mass_def)
    hmc = ccl.halos.HMCalculator(cosmo, hmf, hbf, mass_def)
    prf = ccl.halos.HaloProfilePressureArnaud(mass_bias=1./float(p_dict['B']))
    pk = ccl.halos.halomod_Pk2D(cosmo, hmc, prf, get_2h=False)
    tr_y = ccl.tSZTracer(cosmo, z_max=4.)


    dndlnM_ccl = []
    for zz in dndlnM_redshifts:
        dndlnM_ccl_zz = []
        for mm in dndlnM_masses:
            dndlnM_ccl_p = hmf.get_mass_function(cosmo, mm/float(p_dict['h']), 1/(1+zz))*float(p_dict['h'])**-3.*(1./np.log(10.))
            dndlnM_ccl_zz.append(dndlnM_ccl_p)
        dndlnM_ccl_zz = np.asarray(dndlnM_ccl_zz)
        dndlnM_ccl.append(dndlnM_ccl_zz)

    c_ell_yy_ccl = []
    for ell in multipoles:
        c_ell_yy_ccl.append(ccl.angular_cl(cosmo, tr_y, tr_y, ell, p_of_k_a=pk)*(ell)*(ell+1.)/2./np.pi*1e12)
    c_ell_yy_ccl =  np.asarray(c_ell_yy_ccl)


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
    ax.set_yscale('linear')

    ax.set_xlabel(r'$M\,\, [\mathrm{M_\odot}/h]$',size=title_size)
    ax.set_ylabel(r'$\Delta dndlnM/dndlnM$',size=title_size)

    x_axis_class_sz = dndlnM_masses
    color=iter(cm.viridis(np.linspace(1,0,len(dndlnM_redshifts))))
    for i in range(len(dndlnM_redshifts)):
        col =next(color)
        y_axis_class_sz = dndlnM[:,i]
        y_axis_ccl = dndlnM_ccl[i]
        #ax.plot(x_axis_class_sz,y_axis_class_sz,color='k',ls='-',alpha = 1.,label = "class_sz")
        label = "z=%.2e"%dndlnM_redshifts[i]
        ax.plot(x_axis_class_sz,y_axis_ccl/y_axis_class_sz-1.,color=col,ls='-',alpha = 1.,label = label)

    ax1.legend(loc=2)
    plt.draw()
    plt.ylim(-2.,2.)

    fig.tight_layout()
    FIG_NAME = '/dndlnM'
    fig.savefig(FIG_DIR + FIG_NAME +".pdf")
    plt.show(block=False)


    #cl^yy
    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(7,10))
    ax = ax1
    ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
    ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_yticklabels(), rotation='vertical', fontsize=label_size)
    plt.setp(ax.get_xticklabels(), fontsize=label_size)
    ax.grid( b=True, which="both", alpha=0.3, linestyle='--')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(r'$\ell$',size=title_size)
    ax.set_ylabel(r'$10^{12}\ell(\ell+1)C_\ell^{yy}/2\pi$',size=title_size)

    x_axis_class_sz = multipoles
    color='k'

    col = color
    y_axis_class_sz = cl_1h
    y_axis_ccl = c_ell_yy_ccl

    ax.plot(x_axis_class_sz,y_axis_class_sz,c='k',ls='-',alpha = 1.,label = 'class_sz')
    ax.plot(x_axis_class_sz,y_axis_ccl,c='r',ls='-',alpha = 1,label='ccl')

    ax1.legend(loc=2)
    plt.draw()

    ax = ax2
    ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
    ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_yticklabels(), rotation='vertical', fontsize=label_size)
    plt.setp(ax.get_xticklabels(), fontsize=label_size)
    ax.grid( b=True, which="both", alpha=0.3, linestyle='--')

    ax.set_xscale('log')
    ax.set_yscale('linear')

    ax.set_xlabel(r'$\ell$',size=title_size)
    ax.set_ylabel(r'$C_\ell^{yy,\mathrm{ccl}}/C_\ell^{yy,\mathrm{class\_sz}}-1$'+'   [%]',size=title_size)

    x_axis_class_sz = multipoles
    color='k'
    col = color
    y_axis_class_sz = cl_1h
    y_axis_ccl = c_ell_yy_ccl
    ax.plot(x_axis_class_sz,100.*(y_axis_ccl-y_axis_class_sz)/y_axis_class_sz,c='r',ls='-',alpha = 1)

    ax1.legend(loc=2)
    plt.draw()

    fig.tight_layout()
    FIG_NAME = '/cl_1h'
    fig.savefig(FIG_DIR + FIG_NAME +".pdf")
    plt.show(block=True)




def main():
	parser=argparse.ArgumentParser(description="compute and compare n(z) with CCL")
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)

if __name__=="__main__":
	main()
