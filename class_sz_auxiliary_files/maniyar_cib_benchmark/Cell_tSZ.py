from headers_constants import *

# In this code, the mass input is m500.


class cl_tsz(object):

    def __init__(self, data_var):
        self.dv = data_var
        self.nu = self.dv.nutsz  # please note that this is in Hz instead of GHz
        self.m = self.dv.m500
        self.z = self.dv.z
        self.cosmo = cosmo
        self.delta_h = self.dv.delta_h_tsz  # 500 for tSZ
        self.x = self.dv.x
        self.ell = self.dv.ell
        self.B = self.dv.B  # mass bias
        self.hmf = self.dv.hmf
        self.power = self.dv.Pk_int
        self.biasmz = self.dv.bias_m_z
        self.M_tilde = self.m/self.B

    def f_nu(self):
        """
        Normally it is calculated as:
        x = h_p*self.nu / k_B / T_cmb
        gnu = x*((np.exp(x) + 1)/(np.exp(x) - 1)) - 4
        However, we have to convolve it with the bandpass filter at a given
        frequency for a given experiment. Here we directly use the values for
        f_nu which are convolved with the Planck bandpass filters at
        100, 143, 217, 353, 545, and 857 GHz.
        """
        if self.dv.exp['name'] == 'Planck':
            result = np.array([-4.031, -2.785, 0.187, 6.205, 14.455, 26.335])
        else:
            print ("bandpassed values for f_nu factor not available for %s."+
                   "Using standard non-bandpassed formula" % (self.dv.exp['name']))
            x = h_p*self.nu / k_B / T_cmb
            gnu = x*((np.exp(x) + 1)/(np.exp(x) - 1)) - 4
            result = gnu
        return result

    def r_delta(self):
        """
        radius of the halo containing amount of matter corresponding to delta
        times the critical density of the universe inside that halo
        """
        rho_crit = (self.cosmo.critical_density(self.z)).to(u.Msun/u.Mpc**3).value
        r3 = 3*self.M_tilde/(4*np.pi*self.delta_h*rho_crit)
        return r3**(1./3.)  # dim m,z in Mpc

    def ell_delta(self):
        return self.cosmo.angular_diameter_distance(self.z).value/self.r_delta()
    # dim m,z unitless

    def E_z(self):
        return np.sqrt(self.cosmo.Om0*(1+self.z)**3 + self.cosmo.Ode0)

    def C(self):
        M_tilde = self.M_tilde
        a = 1.65*(self.cosmo.h/0.7)**2*self.E_z()**(8./3)
        b = ((self.cosmo.h/0.7)*M_tilde/(3e14))**(2./3 + 0.12)
        return a*b  # dim m,z final units are eV*cm**-3

    def P_e(self):
        # pressure profile
        # values of the constants taken from
        # https://www.aanda.org/articles/aa/pdf/2013/02/aa20040-12.pdf
        C_t = self.C()*eV_to_J/cm_to_m**3  # converting to SI units
        gamma = 0.31
        alpha = 1.33
        beta = 4.13
        P_0 = 6.41
        c_500 = 1.81
        a = C_t[:, :, None]*P_0*(c_500*self.x)**-gamma
        b = (1+(c_500*self.x)**alpha)**((gamma-beta)/alpha)
        return a*b  # fin dim=m,z,x  J/m^3

    """
    def y_ell(self):
        r500 = self.r_delta()*Mpc_to_m  # m,z
        l500 = self.ell_delta()  # m,z
        a = (sig_T/(m_e*(c_light*Km_to_m)**2))*(4*np.pi*r500/l500**2)  # m,z

        integral = np.zeros((len(self.ell), len(self.m[:, 0]), len(self.z), len(self.x)))
        Pe = self.P_e()  # m,z,x
        Pex2 = Pe*self.x**2  # m,z,x
        x_ls = self.x/l500[:, :, None]  # m,z,x
        for i in range(len(self.ell)):
            integral[i, :, :, :] = Pex2*np.sin(self.ell[i]*x_ls)/(self.ell[i] *
                                                                  x_ls)
        # dim of integral m,z,ell,x
        intgn = intg.simps(integral, x=self.x, axis=-1, even='avg')  # ell,m,z
        return intgn*a  # dim ell,m,z  unitless
    """

    def y_ell_tab(self):
        """
        Integrated pressure profile
        to incresase the calculation speed, we have already tabulated the
        $y_\ell$ values.
        """
        r500 = self.r_delta()*Mpc_to_m  # m,z
        l500 = self.ell_delta()  # m,z
        fact = (sig_T/(m_e*(c_light*Km_to_m)**2))*(4*np.pi*r500/l500**2)  # m,z

        C_t = self.C()*eV_to_J/cm_to_m**3  # converting to SI units
        P_0 = 6.41
        yl = np.loadtxt('data_files/y_ell_integration.txt')
        intgn = np.zeros((len(self.ell), len(self.m[:, 0]), len(self.z)))

        for i in range(len(self.ell)):
            for j in range(len(self.m[:, 0])):
                l_l500 = self.ell[i]/l500[j, :]  # z
                y_int = np.interp(np.log(l_l500), yl[:, 0], yl[:, 1])
                # if min(np.log(l_l500)) < min(yl[:, 0]) or max(np.log(l_l500)) > max(yl[:, 0]):
                intgn[i, j, :] = np.exp(y_int)
        return P_0*intgn*fact*C_t  # dim ell,m,z  unitless

    def dVc_dz(self):  # dim z
        return c_light*self.cosmo.comoving_distance(self.z).value**2/(self.cosmo.H0.value*self.E_z())

    def C_ell_1h(self):
        if self.dv.exp['name'] == 'Planck':
            """
            Kcmb_MJy are factors for Planck frequency channels to convert
            units from Kcmb to MJy
            """
            Kcmb_MJy = np.array([244.1, 371.74, 483.69, 287.45, 58.04, 2.27])
        else:
            print ("factors to convert units from Kcmb to MJy for %s experiment are not provided." +
                   "So the final units here will be Kcmb^2" % (self.dv.exp['name']))
            Kcmb_MJy = np.ones(len(self.nu))

        cl = np.zeros((len(self.nu), len(self.nu), len(self.ell)))
        a_z = self.dVc_dz()  # z
        # y_l = self.y_ell()  # ell,m,z  hmf=m,z
        y_l = self.y_ell_tab()  # ell,m,z  hmf=m,z

        y_l2 = y_l**2

        intgral1 = self.hmf*y_l2  # ell,m,z

        dlogm = np.log10(self.m[1, 0] / self.m[0, 0])
        # intgn1 = intg.simps(intgral1, dx=dlogm, axis=1, even='avg')  # ell,z
        intgn1 = intg.simps(intgral1, x=np.log10(self.m[:, 0]), axis=1, even='avg')  # ell,z

        intgral2 = a_z*intgn1

        res = intg.simps(intgral2, x=self.z, axis=1, even='avg')  # ell

        fnu = self.f_nu()*1e6*Kcmb_MJy
        for f in range(len(self.nu)):
            cl[f, :, :] = np.outer(fnu, res)*fnu[f]  # *T_cmb**2

        return cl

    def tsz_hmf_bias(self):
        # y_l = self.y_ell()  # ell,m,z
        y_l = self.y_ell_tab()  # ell,m,z
        intgrl = self.hmf*self.biasmz*y_l  # # ell,m,z  hmf,bias=m,z
        dlogm = np.log10(self.m[1, 0] / self.m[0, 0])
        # fin = intg.simps(intgrl, dx=dlogm, axis=1, even='avg')  # ell,z
        res = intg.simps(intgrl, x=np.log10(self.m[:, 0]), axis=1, even='avg')  # ell,z
        return res**2  # ell,z

    def C_ell_2h(self):
        if self.dv.exp['name'] == 'Planck':
            """
            Kcmb_MJy are factors for Planck frequency channels to convert
            units from Kcmb to MJy
            """
            Kcmb_MJy = np.array([244.1, 371.74, 483.69, 287.45, 58.04, 2.27])
        else:
            print ("factors to convert units from Kcmb to MJy for %s experiment are not provided." +
                   "So the final units here will be Kcmb^2" % (self.dv.exp['name']))
            Kcmb_MJy = np.ones(len(self.nu))

        cl = np.zeros((len(self.nu), len(self.nu), len(self.ell)))
        V_z = self.dVc_dz()  # z
        fnu = self.f_nu()*1e6*Kcmb_MJy  # nu  # Jy units
        ylhmfbias2 = self.tsz_hmf_bias()  # ell,z
        intgrl = V_z*self.power*ylhmfbias2  # ell,z
        res = intg.simps(intgrl, x=self.z, axis=1, even='avg')  # ell
        for f in range(len(self.nu)):
            cl[f, :, :] = np.outer(fnu, res)*fnu[f]  # *T_cmb**2
        return cl

    def cltot(self):
        cl1h = self.C_ell_1h()
        cl2h = self.C_ell_2h()
        return cl1h+cl2h
