from headers_constants import *
from Cell_cib import *
from Cell_tSZ import *
from Cell_CIBxtSZ import *
from plot_cell import *
from input_var import *

"""
Although we are calculating the halo mass function, halo bias, and the
Fourier transform of the NFW profile here, the computation can be speeded up
by precomputing them before and storing them in a file and then reading
them here.
"""

Planck = {'name': 'Planck',
          'do_cib': 1, 'do_tsz': 1, 'do_cibxtsz': 1,
          'freq_cib': [100., 143., 217., 353., 545., 857.],
          'cc': np.array([1.076, 1.017, 1.119, 1.097, 1.068, 0.995, 0.960]),
          'cc_cibmean': np.array([1.076, 1.017, 1.119, 1.097, 1.068, 0.995, 0.960]),
          'freq_cibmean': np.array([100., 143., 217., 353., 545., 857.]),
          'fc': np.ones(7),
          }


# Planck = {'name': 'Planck',
#           'do_cib': 1, 'do_tsz': 0, 'do_cibxtsz': 0,
#           'freq_cib': [ 545.],
#           'cc': np.array([1.068]),
#           'cc_cibmean': np.array([1.068]),
#           'freq_cibmean': np.array([545.]),
#           'fc': np.ones(1),
#           }

# Herschel = {'name': 'Herschel-spire',
#             'do_cib': 1, 'do_tsz': 0, 'do_cibxtsz': 0,
#             'freq_cib': [600., 857., 1200.],
#             'cc': np.array([0.974, 0.989, 0.988]),
#             'cc_cibmean': np.array([0.974, 0.989, 0.988]),
#             'freq_cibmean': np.array([600., 857., 1200.]),
#             'fc': np.ones(3),
#             }
#
# CCAT = {'name': 'CCAT-p',
#         'do_cib': 1, 'do_tsz': 0, 'do_cibxtsz': 0,
#         'freq_cib': [220., 280., 350., 410., 850.],
#         'cc': np.ones(5),
#         'cc_cibmean': np.ones(5),
#         'freq_cibmean': np.array([220., 280., 350., 410., 850.]),
#         'fc': np.ones(5),
#          }

exp = Planck

# ############### planck cib data #########################

ell = np.geomspace(10., 5e4, 100)
redshifts = np.loadtxt('data_files/redshifts.txt')
z1 = np.linspace(min(redshifts), max(redshifts), 200)

z = redshifts  # z1  # redshifts

logmass = np.arange(8., 15.00, 0.05)

print('len logmass',len(logmass))
mass = 10**logmass

driver = data_var(exp, mass, z, ell)

# exit(0)

# if exp['do_cib'] == 1:
clcib = cl_cib(driver)

cl1h_cib = clcib.onehalo_int()
cl2h_cib = clcib.twohalo_int()

# plotting the CIB power spectra for freq[nu1]xfreq[nu2] GHz

# fam = "serif"
# plt.rcParams["font.family"] = fam



#### do class_sz computation:


Omegam0 = 0.3075
H0 = 67.74
Omegab = 0.0486
Omegac = 0.2589
Omegac + Omegab
omega_b = Omegab*(H0/100.)**2
omega_c = Omegac*(H0/100.)**2
hparam = H0/100.
maniyar_cosmo = {
'omega_b': omega_b,
'omega_cdm':  omega_c,
'h': H0/100.,
# 'tau_reio': 0.0561,
'ln10^{10}A_s': 3.048,
'n_s': 0.9665,
# 'sigma8':0.830,
'k_pivot': 0.05,
'N_ncdm': 1,
'N_ur': 2.0328,
'm_ncdm': 0.0
}
print('computing class_sz:')
from classy_sz import Class
M = Class()
M.set({'output':'dndlnM,cib_cib_1h,cib_cib_2h'})
# M.set({'output':'cib_cib_1h'})
# M.set(common_settings)
M.set(maniyar_cosmo)
# M.set(websky_cib_params)
M.set({

'mass function' : 'T08M200c',
'use_maniyar_cib_model':1,

'maniyar_cib_etamax' : 5.12572945e-01,#0.42,

'maniyar_cib_zc' : 1.5,
'maniyar_cib_tau' : 8.25475287e-01,#1.17,
'maniyar_cib_fsub' : 0.134*np.log(10.),
'Most efficient halo mass in Msun' : 5.34372069e+12,#10.**12.94,
'Size of of halo masses sourcing CIB emission' :  1.5583436676980493,#1.75**2.,
#for the Lsat tabulation:
'freq_min': 9e1,
'freq_max': 8.57e2,
'dlogfreq' : 0.1,

'concentration parameter':'fixed',

'n_z_L_sat' :200,
'n_m_L_sat' :200,
'n_nu_L_sat':200,


# 'Redshift evolution of dust temperature' :  0.2,
# 'Dust temperature today in Kelvins' : 20.7,
# 'Emissivity index of sed' : 1.6,
# 'Power law index of SED at high frequency' : 1.7, # not given in WebSky paper, actually not relevant since we dont use high freqs in websky.
# 'Redshift evolution of L − M normalisation' : 1.28, # try 2.4 see slack.
# 'Most efficient halo mass in Msun' : 10.**12.3,
# 'Normalisation of L − M relation in [Jy MPc2/Msun]' : 1e-7,  # not given in WebSky paper
# 'Size of of halo masses sourcing CIB emission' : 0.3,
# 'z_plateau_cib' : 2.,

# M_min_HOD is the threshold above which nc = 1:
# 'M_min_HOD' : 10.**10.1, # not used here
'use_nc_1_for_all_halos_cib_HOD': 1,

'sub_halo_mass_function' : 'TW10',#'JvdB14',
'M_min_subhalo_in_Msun' : 1e5, # 1e5 see https://github.com/abhimaniyar/halomodel_cib_tsz_cibxtsz/blob/master/Cell_cib.py
'use_redshift_dependent_M_min': 0,
#'full_path_to_redshift_dependent_M_min':'/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/sz_auxiliary_files/websky_halo_mass_completion_z_Mmin_in_Msun_over_h.txt',
'M_min' : 1e7*hparam,#mabhi.min()*maniyar_cosmo['h'], # not used
'M_max' : 1e15*hparam,#mabhi.max()*maniyar_cosmo['h'],
'z_min' : 0.012,
'z_max' : 10.,
'ell_min': 10.,
'ell_max':5e4,
'dlogell':0.6,


'ndim_redshifts': 210,
'ndim_masses':150,
# table 1 of https://arxiv.org/pdf/1309.0382.pdf
#1: freq GHz 2: Flux cut mJy
# 100 - 400
# 143 - 350
# 217 - 225
# 353 - 315
# 545 - 350
# 857 - 710
# 3000  - 1000
#cib_Snu_1 = 315.
#cib_Snu_2 = 315.
#'cib_Snu_cutoff_list [mJy]':'315',
'has_cib_flux_cut': 0,
'hm_consistency':0,

'epsabs_L_sat': 1e-40,
'epsrel_L_sat': 1e-9,

# "P_k_max_1/Mpc": 50.,
# 'k_max_for_pk_class_sz':50.
})

M.set({
       'cib_frequency_list_num' : 6,
       'cib_frequency_list_in_GHz' : '100,143,217,353,545,857',
      })
M.compute()
cl_cib_cib = M.cl_cib_cib()
print(cl_cib_cib)

freqs = '545'
Planck_cib_dict = {'name': 'Planck',
          'do_cib': 1, 'do_tsz': 1, 'do_cibxtsz': 1,
          'freq_cib': [100., 143., 217., 353., 545., 857.],
          'cc': np.array([1.076, 1.017, 1.119, 1.097, 1.068, 0.995, 0.960]),
          'cc_cibmean': np.array([1.076, 1.017, 1.119, 1.097, 1.068, 0.995, 0.960]),
          'freq_cibmean': np.array([100., 143., 217., 353., 545., 857.]),
          'fc': np.ones(7),
          }
idfreqabh = 4
faccib = 1.#Planck_cib_dict['cc'][idfreqabh]**2

# a simple conversion from cl's to dl's
def l_to_dl(lp):
    return lp*(lp+1.)/2./np.pi

freq = ['100', '143', '217', '353', '545', '857']
# freq = ['545']
# nu1, nu2 = 4, 4
# plot_Cell(ell, cl1h_cib, cl2h_cib, nu1, nu2, freq, 'CIB')
two_halo = cl2h_cib
one_halo = cl1h_cib

# plot_dim = int(1/2*(-1 + np.sqrt(1 + 8*len(freq))))
plot_dim = len(freq)
f, axes = plt.subplots(figsize=(8, 8),
                         sharex=True,
                         #sharey=True,
                         ncols=plot_dim,
                         nrows=plot_dim)

plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

for i in range(plot_dim):
    for j in range(plot_dim):
        if i<j:
            axes[i, j].axis('off')
        else:
            ax = axes[i, j]
            freqs = freq[i]
            freqsp = freq[j]
            nu1 = i
            nu2 = j
            ax.plot(ell, np.abs(one_halo[nu1, nu2, :]), 'b-.', label='Maniyar et al 1h',lw=0.6)
            ax.plot(ell, np.abs(two_halo[nu1, nu2, :]), 'b-', label='Maniyar et al 2h',lw=0.6)
            # ax.plot(ell, np.abs(total[nu1, nu2, :]), 'b', label='total')
            l = np.asarray(cl_cib_cib[freqs+'x'+freqs]['ell'])
            ax.plot(l,cl_cib_cib[freqs+'x'+freqsp]['1h']/l_to_dl(l)*faccib+cl_cib_cib[freqs+'x'+freqsp]['2h']/l_to_dl(l)*faccib*0.,ls='-.',c='k',label='class_sz 1h ')
            ax.plot(l,cl_cib_cib[freqs+'x'+freqsp]['1h']/l_to_dl(l)*faccib*0.+cl_cib_cib[freqs+'x'+freqsp]['2h']/l_to_dl(l)*faccib,ls='-',c='k',label='class_sz 2h ')
            ax.set_xscale("log")
            ax.set_yscale("log")
            # ax.legend(loc='upper right', prop={'size': 8}, frameon=False)
            # ax.set_xticks([100, 500, 1000])
            # ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            # ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.set_ylabel(r'$\mathrm{C_l}\: [\mathrm{Jy}^2\: \mathrm{sr}^{-1}]$', fontsize=7)
            ax.set_xlabel(r'$\;\ell$', fontsize=7)
# ax.set_title(r''+mod+' %s x %s GHz' % (freq[nu1], freq[nu2]))
plt.show(block=True)


    # plt.figure()

# ############################### tSZ params ############################
#
# if exp['do_tsz'] == 1:
#     cltsz = cl_tsz(driver)
#     cl1h_tsz = cltsz.C_ell_1h()
#     cl2h_tsz = cltsz.C_ell_2h()
#     # self.B = 1.41
#
#     # plotting the tSZ power spectra for freq[nu1]xfreq[nu2] GHz
#     freq = ['100', '143', '217', '353', '545', '857']
#     nu1, nu2 = 0, 0
#     plot_Cell(ell, cl1h_tsz, cl2h_tsz, nu1, nu2, freq, 'tSZ')
#     # print (cl1h_tsz[0, 0])
#     # print (cl2h_tsz[0, 0])
#
# # ################################ cib x tSZ ########################
#
# if exp['do_cibxtsz'] == 1:
#     # """
#     cib_cls = cl_cib(driver)
#     tsz_cls = cl_tsz(driver)
#     cibtsz = cl_cibxtsz(cib_cls, tsz_cls)
#     cl1h_cibtsz = cibtsz.onehalo()  # *Kcmb_MJy*1e6
#     cl2h_cibtsz = cibtsz.twohalo()
#
#     # plotting the CIBxtSZ power spectra for freq[nu1]xfreq[nu2] GHz
#     freq = ['100', '143', '217', '353', '545', '857']
#     nu1, nu2 = 0, 0
#     plot_Cell(ell, cl1h_cibtsz, cl2h_cibtsz, nu1, nu2, freq, 'CIB x tSZ')
    # """
