from headers_constants import *


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
M.set({'output':'dndlnM'})
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
'maniyar_cib_fsub' : 0.,#0.134*np.log(10.),
'Most_efficient_halo_mass_in_Msun' : 5.34372069e+12,#10.**12.94,
'Size_of_halo_masses_sourcing_CIB_emission' :  1.5583436676980493,#1.75**2.,
#for the monopole computation:
'freq_min': 9e1,
'freq_max': 8e2,
'dlogfreq' : 0.1,

'concentration parameter':'fixed',




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
'M_min' : 1e8*hparam,#mabhi.min()*maniyar_cosmo['h'], # not used
'M_max' : 1e15*hparam,#mabhi.max()*maniyar_cosmo['h'],
'z_min' : 0.012,
'z_max' : 11.,
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

"P_k_max_1/Mpc": 100.,
'k_max_for_pk_class_sz':100.
})

# M.set({
#        'cib_frequency_list_num' : 1,
#        'cib_frequency_list_in_GHz' : '545',
#       })
M.compute()
# cl_cib_cib = M.cl_cib_cib()
print('finnished class_sz hmf computation')
zp = 1.3

logmass = np.arange(8, 15.00, 0.05)
print('len logmasses in msun:',len(logmass))
mabhi = 10.**logmass
dndlogm_class_sz = np.vectorize(M.get_dndlnM_at_z_and_M)(zp,mabhi*hparam)*np.log(10.)*hparam**3
print(dndlogm_class_sz)
# exit(0)



class data_var(object):

    def __init__(self, exp, mass, z, ell):
        # ############### cib data #########################
        self.exp = exp
        name = self.exp['name']
        deltah_cib = 200
        self.z_c = 1.5
        self.freqcib = self.exp['freq_cib']
        self.cc = self.exp['cc']
        self.fc = self.exp['fc']

        self.cc_cibmean = self.exp['cc_cibmean']
        self.freq_cibmean = self.exp['freq_cibmean']

        self.mass = mass
        self.z = z
        self.ell = ell

        nm = len(self.mass)
        nz = len(self.z)

        # ########## reading in the matter power spectra #############
        redshifts = np.loadtxt('data_files/redshifts.txt')

        if min(self.z) < min(redshifts) or max(self.z) > max(redshifts):
            print ("If the redshift range is outside of [%s to %s], then " +
                   "values of the matter power spectrum and effective " +
                   "CIB SEDs are extrapolated and might be incorrect.") % (min(redshifts), max(redshifts))
        ll = [str(x) for x in range(1, 211)]
        """
        please note that the matter power spectrum files here are arranged in
        a reversed order i.e. redshift decreases as you go from _1 file to _210.
        Perhaps better to generate your own power spectra for the redshift
        array you want to consider and read that here. Also, note that I am
        getting rid of the reduced Hubble constant here from units of k and Pk.
        """
        addr = 'data_files/matter_power_spectra'
        pkarray = np.loadtxt('%s/test_class_sz_highk_lin_matterpower_210.dat' % (addr))
        k = pkarray[:, 0]*cosmo.h
        Pk = np.zeros((len(k), len(redshifts)))
        for i in range(len(redshifts)):
            pkarray = np.loadtxt("%s/test_class_sz_highk_lin_matterpower_%s.dat" % (addr, ll[209-i]))
            Pk[:, i] = pkarray[:, 1]/cosmo.h**3

        pkinterpz = interp1d(redshifts, Pk, kind='linear', bounds_error=False, fill_value="extrapolate")

        self.k_array = np.zeros((len(self.ell), len(self.z)))
        self.Pk_int = np.zeros(self.k_array.shape)
        """
        Pk_int 2-d array for corresponding redshifts and
        given ell range such that k = ell/chi i.e. for every
        redshift
        """
        chiz = cosmo.comoving_distance(self.z).value
        for i in range(len(self.ell)):
            self.k_array[i, :] = self.ell[i]/chiz
            for j in range(len(self.z)):
                pkz = pkinterpz(self.z[j])
                self.Pk_int[i, j] = np.interp(self.k_array[i, j], k, pkz)


        if self.exp['do_cib'] == 1 or self.exp['do_cibxtsz'] == 1:
            # ######### reading and interpolating the SEDs
            """
            The effective SEDs for the CIB for Planck (100, 143, 217, 353, 545,
            857) and
            IRAS (3000) GHz frequencies.
            Here we are shwoing the CIB power spectra corressponding to the
            Planck
            frequency channels. If you want to calculate the Hershel/Spire
            power spectra, use corresponding files in the data folder.
            """
            if name == 'Planck':
                snuaddr = 'data_files/filtered_snu_planck.fits'
                hdulist = fits.open(snuaddr)
                redshifts = hdulist[1].data
                snu_eff = hdulist[0].data  # in Jy/Lsun
                hdulist.close()
                fsnu_eff = interp1d(redshifts, snu_eff, kind='linear',
                                    bounds_error=False, fill_value="extrapolate")
                self.snu = fsnu_eff(self.z)
            elif name == 'Herschel-spire':
                snuaddr = 'data_files/filtered_snu_spire.fits'
                hdulist = fits.open(snuaddr)
                redshifts = hdulist[1].data
                snu_eff = hdulist[0].data  # in Jy/Lsun
                hdulist.close()
                fsnu_eff = interp1d(redshifts, snu_eff, kind='linear',
                                    bounds_error=False, fill_value="extrapolate")
                self.snu = fsnu_eff(self.z)
            else:
                # ######### unfiltered SEDs ###########################

                list_of_files = sorted(glob.glob('data_files/TXT_TABLES_2015/./*.txt'))
                a = list_of_files[95]
                b = list_of_files[96]
                for i in range(95, 208):
                    list_of_files[i] = list_of_files[i+2]
                list_of_files[208] = a
                list_of_files[209] = b

                wavelengths = np.loadtxt('data_files/TXT_TABLES_2015/EffectiveSED_B15_z0.012.txt')[:, [0]]
                # wavelengths are in microns
                freq = c_light/wavelengths
                # c_light is in Km/s, wavelength is in microns and we would like to
                # have frequency in GHz. So gotta multiply by the following
                # numerical factor which comes out to be 1
                # numerical_fac = 1e3*1e6/1e9
                numerical_fac = 1.
                freqhz = freq*1e3*1e6
                freq *= numerical_fac
                freq_rest = freqhz*(1+redshifts)

                n = np.size(wavelengths)

                snu_unfiltered = np.zeros([n, len(redshifts)])
                for i in range(len(list_of_files)):
                    snu_unfiltered[:, i] = np.loadtxt(list_of_files[i])[:, 1]
                L_IR15 = self.L_IR(snu_unfiltered, freq_rest, redshifts)
                # print (L_IR15)

                for i in range(len(list_of_files)):
                    snu_unfiltered[:, i] = snu_unfiltered[:, i]*L_sun/L_IR15[i]

                # Currently unfiltered snus are ordered in increasing wavelengths,
                # we re-arrange them in increasing frequencies i.e. invert it

                freq = freq[::-1]
                snu_unfiltered = snu_unfiltered[::-1]
                fsnu_unfiltered = RectBivariateSpline(freq, redshifts,
                                                      snu_unfiltered)

                nuinp = self.freqcib
                self.snu = fsnu_unfiltered(nuinp, self.z)


        if self.exp['do_cib'] == 1:
            # ######### CIB halo model parameters ###################
            cibparresaddr = 'data_files/one_halo_bestfit_allcomponents_lognormal_sigevol_1p5zcutoff_nospire_fcpl_onlyautoshotpar_no3000_gaussian600n857n1200_planck_spire_hmflog10.txt'
            self.Meffmax, self.etamax, self.sigmaMh, self.tau = np.loadtxt(cibparresaddr)[:4, 0]
            # self.Meffmax, self.etamax, self.sigmaMh, self.tau = 8753289339381.791, 0.4028353504978569, 1.807080723258688, 1.2040244128818796

            print(
            'self.Meffmax, self.etamax, self.sigmaMh, self.tau',self.Meffmax,
            self.etamax,
            self.sigmaMh,
            self.tau
            )

            # ######## hmf, bias, nfw ###########
            print ("Calculating the halo mass function, halo bias, nfw " +
                   "profile " +
                   "for given mass and redshift for CIB calculations.")

            self.hmf = np.zeros((nm, nz))
            dndlogm_class_sz = np.zeros((nm, nz))
            self.u_nfw = np.ones((nm, len(self.k_array[:, 0]), nz))
            self.bias_m_z = np.zeros((nm, nz))
            delta_h = deltah_cib

            for r in range(nz):
                print('---------  computing stuuff at z = ',self.z[r])
                pkz = pkinterpz(self.z[r])
                instance = hmf_unfw_bias.h_u_b(k, pkz, self.z[r],
                                               cosmo, delta_h, self.mass)
                # self.hmf[:, r] = instance.dn_dlogm()
                # nfw_u[:, :, r] = instance.nfwfourier_u()
                # self.bias_m_z[:, r] = instance.b_nu()
                nucsz = np.vectorize(M.get_nu_at_z_and_m)(self.z[r],mabhi*hparam)
                bcsz = np.vectorize(M.get_first_order_bias_at_z_and_nu)(self.z[r],nucsz)
                self.bias_m_z[:, r] = bcsz
                dndlogm_class_sz[:,r] = np.vectorize(M.get_dndlnM_at_z_and_M)(self.z[r],mabhi*hparam)*np.log(10.)*hparam**3

                # print('hmf_shape:',np.shape(self.hmf))
                # print(self.hmf[:,r]/dndlogm_class_sz[:,r])
                self.hmf[:,r] = dndlogm_class_sz[:,r]
                # print()
                # exit(0)
                # print('bias_shape:',np.shape(self.bias_m_z))
                # instance2 = hmf_unfw_bias.h_u_b(self.k_array[:, r],
                #                                 self.Pk_int[:, r], self.z[r],
                #                                 cosmo, delta_h, self.mass)
                # dim(len(m), len(k))
                rd = np.vectorize(M.get_r_delta_of_m_delta_at_z)(200.,mabhi*hparam,self.z[r])
                kp = self.k_array[:, r]#self.k_array[i, :] = self.ell[i]/chiz
                ucsz =  np.ones((nm, len(self.k_array[:, 0]), nz))
                for ikpp,kpp in enumerate(kp):
                    ucsz[:,ikpp,r] = np.vectorize(M.get_truncated_nfw_profile_at_z_k_rd_cd_xout)(self.z[r],kpp/hparam,rd,5.,1.)
                print(np.shape(ucsz))
                self.u_nfw[:, :, r] = ucsz[:,:,r]
                # self.u_nfw[:, :, r] = instance2.nfwfourier_u()#/instance2.nfwfourier_u()

            # exit(0)

        if self.exp['do_tsz'] == 1 or self.exp['do_cibxtsz'] == 1:
            # ############################### tSZ params #####################
            xstep = 50
            lnx = np.linspace(-6, 1, xstep)
            self.x = 10**lnx
            # self.nutsz = np.array([100., 143., 217., 353., 545., 857.])*ghz
            self.nutsz = np.array(self.freqcib)*ghz
            #nus = ['100', '143', '217', '353', '545', '857']
            self.delta_h_tsz = 500  # 500 # 200
            self.B = 1.5  # 1.41
            self.m500 = np.repeat(self.mass[..., np.newaxis], len(self.z),
                                  axis=1)

            print ("Calculating the halo mass function, halo bias, nfw " +
                   "profile " +
                   "for given mass and redshift for tSZ calculations.")
            # self.hmf = np.zeros((len(self.m500), nz))
            # self.bias_m_z = np.zeros((len(self.m500[:, 0]), nz))
            # self.u_nfw = np.zeros((nm, len(self.k_array[:, 0]), nz))
            # delta_h = self.delta_h_tsz
            #
            # for r in range(nz):
            #     pkz = pkinterpz(self.z[r])
            #     instance = hmf_unfw_bias.h_u_b(k, pkz, self.z[r],
            #                                    cosmo, delta_h, self.m500[:, 0])
            #     self.hmf[:, r] = instance.dn_dlogm()
            #     # nfw_u[:, :, r] = instance.nfwfourier_u()
            #     self.bias_m_z[:, r] = instance.b_nu()
            #     instance2 = hmf_unfw_bias.h_u_b(self.k_array[:, r],
            #                                     self.Pk_int[:, r], self.z[r],
            #                                     cosmo, delta_h, self.mass)
            #     self.u_nfw[:, :, r] = instance2.nfwfourier_u()#/instance2.nfwfourier_u()
            # self.hmf = np.zeros((nm, nz))
            # self.u_nfw = np.zeros((nm, len(self.k_array[:, 0]), nz))
            # self.bias_m_z = np.zeros((nm, nz))
            # delta_h = deltah_cib
            #
            # for r in range(nz):
            #     pkz = pkinterpz(self.z[r])
            #     instance = hmf_unfw_bias.h_u_b(k, pkz, self.z[r],
            #                                    cosmo, delta_h, self.mass)
            #     self.hmf[:, r] = instance.dn_dlogm()
            #     # nfw_u[:, :, r] = instance.nfwfourier_u()
            #     self.bias_m_z[:, r] = instance.b_nu()
            #     instance2 = hmf_unfw_bias.h_u_b(self.k_array[:, r],
            #                                     self.Pk_int[:, r], self.z[r],
            #                                     cosmo, delta_h, self.mass)
            #     self.u_nfw[:, :, r] = instance2.nfwfourier_u()#



        if self.exp['do_cibxtsz'] == 1:
            # mass definition here is m500 for both the CIB and tSZ. The best fit values
            # for CIB parameters change due to this.

            self.snu = self.snu[:-1, :]
            self.cc = self.exp['cc'][:-1]
            self.fc = self.exp['fc'][:-1]
            # ################################ cib x tSZ #####################
            # cibxtszparresaddr = 'data_files/one_halo_bestfit_allcomponents_lognormal_sigevol_highk_deltah500_onlyautoshotpar_no3000_gaussian600n857n1200_planck_spire_hmflog10.txt'
            # self.Meffmax, self.etamax, self.sigmaMh, self.tau = np.loadtxt(cibxtszparresaddr)[:4, 0]
            # self.Meffmax_cross, self.etamax_cross, self.sigmaMh_cross = 6962523672799.227, 0.4967291547804018, 1.8074450009861387
            # self.tau_cross = 1.2016980179374213
            cibparresaddr = 'data_files/one_halo_bestfit_allcomponents_lognormal_sigevol_1p5zcutoff_nospire_fcpl_onlyautoshotpar_no3000_gaussian600n857n1200_planck_spire_hmflog10.txt'
            self.Meffmax, self.etamax, self.sigmaMh, self.tau = np.loadtxt(cibparresaddr)[:4, 0]
