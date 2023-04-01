from headers_constants import *
from scipy.interpolate import InterpolatedUnivariateSpline as _spline


# code to calculate the halo mass function, halo bias, and the Fourier
# transform of the NFW profile


class h_u_b:

    def __init__(self, kk, power, z, cosmo, delta_h, mh):
        self.kk = kk
        self.power = power
        self.z = z
        self.delta_h = delta_h
        self.cosmo = cosmo
        self.mh = mh
        self.delta_wrt = 'crit'  # ''mean' or 'crit' mass can be calculated either with
        # respect to the mean background density or critical background density.
        # say 'mean' if we want to calculate it wrt to mean background density

    def mean_density0(self):
        cosmo = self.cosmo
        return (cosmo.Om0*cosmo.critical_density0).to(u.Msun/u.Mpc**3).value

    def mass_to_radius(self):
        """
        Lagrangian radius of a dark matter halo
        """
        rho_mean = self.mean_density0()
        r3 = 3*self.mh/(4*np.pi*rho_mean)
        return r3**(1./3.)

    def r_delta(self):
        """
        radius of the halo containing amount of matter corresponding to delta times
        the critical density of the universe inside that halo
        """
        cosmo = self.cosmo
        rho_crit = (cosmo.critical_density(self.z)).to(u.Msun/u.Mpc**3).value
        r3 = 3*self.mh/(4*np.pi*self.delta_h*rho_crit)
        return r3**(1./3.)

    def W(self, rk):  # Fourier transform of top hat window function
        """
        limit of 1.4e-6 put as it doesn't add much to the final answer and
        helps for faster convergence
        """
        return np.where(rk > 1.4e-6, (3*(sin(rk)-rk*cos(rk)) / rk**3), 1)

    def sigma(self, rad):
        """
        matter variance for given power spectrum, wavenumbers k and radius of
        the halo. radius is calculated below from mass of the halo. It has to
        be noted that there's no factor
        of \Delta = 200 to calculate the radius.
        """
        rk = np.outer(rad, self.kk)
        rest = self.power * self.kk**3
        # dlnk = np.log(self.kk[1] / self.kk[0])
        lnk = np.log(self.kk)
        integ = rest*self.W(rk)**2
        sigm = (0.5/np.pi**2) * intg.simps(integ, x=lnk, axis=-1)
        return np.sqrt(sigm)

    def nu_delta(self):  # peak heights
        rad = self.r_delta()
        """
        to calculate the peak heights: nu_delta,  we use r_delta rather than
        the simple Lagrangian radius calculated using mass_to_radius function.
        This will be used in c-M relation to calculate the NFW profile.
        """
        delta_c = 1.686  # critical density of the universe. Redshift evolution
        # is small and neglected
        sig = self.sigma(rad)
        return delta_c / sig  # length of mass array

    def nu(self):  # peak heights
        rad = self.mass_to_radius()
        """
        to calculate the peak heights: nu, we use the
        simple Lagrangian radius calculated using mass_to_radius function.
        """
        delta_c = 1.686  # critical density of the universe. Redshift evolution
        # is small and neglected
        sig = self.sigma(rad)
        return delta_c / sig  # length of mass array

    def dw_dlnkr(self, rk):
        return np.where(rk > 1e-3, (9*rk*cos(rk)+3*sin(rk)*(rk**2-3))/rk**3, 0)

    def dlns2_dlnr(self, rad):
        rk = np.outer(rad, self.kk)
        rest = self.power * self.kk**3
        w = self.W(rk)
        dw = self.dw_dlnkr(rk)
        inte = w*dw*rest
        lnk = np.log(self.kk)
        s = self.sigma(rad)
        return intg.simps(inte, x=lnk, axis=-1, even='avg')/(np.pi**2*s**2)

    def dlnr_dlnm(self):  # The derivative of log radius with log mass.
        return 1./3.

    def dlns2_dlnm(self, rad):
        return self.dlns2_dlnr(rad) * self.dlnr_dlnm()

    def dlns_dlnm(self, rad):
        return 0.5 * self.dlns2_dlnm(rad)

    def delta_halo(self):
        """ Overdensity of a halo w.r.t mean density or ctirical density"""
        if self.delta_wrt == 'mean':
            return self.delta_h

        elif self.delta_wrt == 'crit':
            return self.delta_h / self.cosmo.Om(self.z)

    def fsigma(self, rad):
        # https://arxiv.org/pdf/0803.2706.pdf
        """
        this is the function giving f(sigma)
        All the paprameters mentioned here are taken for \Delta = 200
        z is the redshift and redshift evolution has to be considered here
        It has to be noted that values provided in Tinker 2008 paper are
        wrt to mean background density. If we want to calculate the best fit
        parameter values wrt to critical background density, we have to change
        the corresponding delta_halo value and interpolate the values for the
        corresponding delta_halo values.
        """
        z = self.z
        a = {  # -- A
             "A_200": 1.858659e-01,
             "A_300": 1.995973e-01,
             "A_400": 2.115659e-01,
             "A_600": 2.184113e-01,
             "A_800": 2.480968e-01,
             "A_1200": 2.546053e-01,
             "A_1600": 2.600000e-01,
             "A_2400": 2.600000e-01,
             "A_3200": 2.600000e-01,
             # -- a
             "a_200": 1.466904,
             "a_300": 1.521782,
             "a_400": 1.559186,
             "a_600": 1.614585,
             "a_800": 1.869936,
             "a_1200": 2.128056,
             "a_1600": 2.301275,
             "a_2400": 2.529241,
             "a_3200": 2.661983,
             # --- b
             "b_200": 2.571104,
             "b_300": 2.254217,
             "b_400": 2.048674,
             "b_600": 1.869559,
             "b_800": 1.588649,
             "b_1200": 1.507134,
             "b_1600": 1.464374,
             "b_2400": 1.436827,
             "b_3200": 1.405210,
             # --- c
             "c_200": 1.193958,
             "c_300": 1.270316,
             "c_400": 1.335191,
             "c_600": 1.446266,
             "c_800": 1.581345,
             "c_1200": 1.795050,
             "c_1600": 1.965613,
             "c_2400": 2.237466,
             "c_3200": 2.439729}
        A_exp = -0.14
        a_exp = 0.06
        delta_virs = np.array([200, 300, 400, 600, 800, 1200, 1600, 2400, 3200])
        dhalo = self.delta_halo()
        if dhalo not in delta_virs:
            A_array = np.array([a["A_%s" % d] for d in delta_virs])
            a_array = np.array([a["a_%s" % d] for d in delta_virs])
            b_array = np.array([a["b_%s" % d] for d in delta_virs])
            c_array = np.array([a["c_%s" % d] for d in delta_virs])

            A_intfunc = _spline(delta_virs, A_array)
            a_intfunc = _spline(delta_virs, a_array)
            b_intfunc = _spline(delta_virs, b_array)
            c_intfunc = _spline(delta_virs, c_array)

            A_0 = A_intfunc(dhalo)
            a_0 = a_intfunc(dhalo)
            b_0 = b_intfunc(dhalo)
            c_0 = c_intfunc(dhalo)
        else:
            A_0 = a["A_%s" % (int(dhalo))]
            a_0 = a["a_%s" % (int(dhalo))]
            b_0 = a["b_%s" % (int(dhalo))]
            c_0 = a["c_%s" % (int(dhalo))]

        s = self.sigma(rad)
        A = A_0 * (1+z)**A_exp
        a = a_0 * (1+z)**a_exp
        alpha = 10**(-(0.75/np.log10(dhalo/75.))**1.2)
        b = b_0 * (1+z)**alpha
        return A * ((s/b)**-a + 1)*np.exp(-c_0/s**2)

    def dn_dm(self):
        rad = self.mass_to_radius()
        return self.fsigma(rad) * self.mean_density0() * np.abs(self.dlns_dlnm(rad)) / self.mh**2

    def dn_dlnm(self):
        return self.mh * self.dn_dm()  # *np.log(10)

    def dn_dlogm(self):
        return self.mh * self.dn_dm() * np.log(10)

# ################ NFW profile calculation #######################

    """
    Code to calculate the Fourier transform of the NFW profile. The analytical
    formula has been taken from arXiv:1206.6890v1 (or arXiv:astro-ph/0006319v2)
    """

    def sine_cosine_int(self, x):
        r"""
        sine and cosine integrals required to calculate the Fourier transform
        of the NFW profile.
        $ si(x) = \int_0^x \frac{\sin(t)}{t} dt \\
        ci(x) = - \int_x^\infty \frac{\cos(t)}{t}dt$
        """
        si, ci = sici(x)
        return si, ci

    def k_R(self):
        """
        we need a c-M relation to calculate the fourier transform of the NFW
        profile. We use the relation from https://arxiv.org/pdf/1407.4730.pdf
        where they use the slope of the power spectrum with respect to the
        wavenumber in their formalism. This power spectrum slope has to be
        evaluated at a certain value of k such that k = kappa*2*pi/rad
        This kappa value comes out to be 0.69 according to their calculations.
        """
        rad = self.mass_to_radius()
        kappa = 0.69
        return kappa * 2 * np.pi / rad

    def dlnpk_dlnk(self):
        """
        When the power spectrum is obtained from CAMB, slope of the ps wrt k
        shows wiggles at lower k which corresponds to the BAO features. Also
        at high k, there's a small dip in the slope of the ps which is due to
        the effect of the baryons (this is not very important for the current
        calculations though). We are using the analysis from the paper
        https://arxiv.org/pdf/1407.4730.pdf where they have used the power
        spectrum from Eisenstein and Hu 1998 formalism where these effects have
        been negelected and therefore they don't have the wiggles and the bump.
        In order to acheive this, we have to smooth out the
        slope of the ps at lower k. But we have checked that the results
        do not vary significantly with the bump.
        """
        grad = np.zeros(self.power.shape, np.float)
        grad[0:-1] = np.diff(np.log(self.power)) / np.diff(np.log(self.kk))
        grad[-1] = (np.log(self.power[-1]) - np.log(self.power[-2]))/(np.log(self.kk[-1]) - np.log(self.kk[-2]))
        # grad_smoothed = savgol_filter(grad, 51, 2)
        grad_smoothed = grad
        kr = self.k_R()
        return np.interp(kr, self.kk, grad_smoothed)

    def nu_c(self):
        nu_delta = self.nu_delta()
        return nu_delta

    def nu_to_c200c(self):  # length of mass array
        use_mean = False  # 2 relations provided. mean and median.
        phi0_median, phi0_mean = 6.58, 7.14
        phi1_median, phi1_mean = 1.37, 1.60
        eta0_median, eta0_mean = 6.82, 4.10
        eta1_median, eta1_mean = 1.42, 0.75
        alpha_median, alpha_mean = 1.12, 1.40
        beta_median, beta_mean = 1.69, 0.67
        _nu = self.nu_c()
        n_k = self.dlnpk_dlnk()
        if use_mean:
            c_min = phi0_mean + phi1_mean * n_k
            nu_min = eta0_mean + eta1_mean * n_k
            return c_min * \
                ((_nu/nu_min)**-alpha_mean + (_nu/nu_min)**beta_mean)/2
        else:
            c_min = phi0_median + phi1_median * n_k
            nu_min = eta0_median + eta1_median * n_k
            return c_min * \
                ((_nu/nu_min)**-alpha_median + (_nu/nu_min)**beta_median)/2

    def r_star(self):
        r"""
        characteristic radius also calld r_s in other literature.
        $ c \equiv \frac{r_{200}{r_s}$
        """
        c_200c = self.nu_to_c200c()
        c_200c = c_200c/c_200c*5.
        r200 = self.r_delta()
        return r200/c_200c  # length of mass array

    def ampl_nfw(self, c):
        r"""
        Dimensionless amplitude of the NFW profile.
        Gives:
            $\frac{1}{\log(1+c) - \frac{c}{1+c}}$
        """
        return 1. / (np.log(1.+c) - c/(1.+c))

    def nfwfourier_u(self):
        rs = self.r_star()*(1.+self.z)  ## should be times (1+z)
        c = self.nu_to_c200c()*5.
        c = c/c*5.
        a = self.ampl_nfw(c)
        mu = np.outer(self.kk, rs)
        Si1, Ci1 = self.sine_cosine_int(mu + mu * c)
        Si2, Ci2 = self.sine_cosine_int(mu)
        unfw = a*(cos(mu)*(Ci1-Ci2) + sin(mu)*(Si1-Si2)-sin(mu*c) / (mu+mu*c))
        return unfw.transpose()  # dim(len(m), len(k))

    """
    def mh_to_c(self):
        halo mass to concentration parameter c relation.
        Here we take the relation from Duffy et al. Arxiv 0804.2486
        There are 3 definitions given: r_virial, r_200 and r_mean.
        We here take the r_mean definition where r_min contains the halo mass
        inside r_mean with the mean internal density 200 times the mean
        background density. We take the same definition of the radius
        while calculating the halo mass function.
        IT HAS TO BE NOTED THAT THIS RELATION WAS CALCULATED USING WMAP5
        COSMOLOGY WITH
        [Omega_M,Omega_b,Omega_Lambda,h,sigma_8,n_s] given by
        [0.258, 0.0441, 0.742, 0.719, 0.796, 0.963].
        But we are using a new cosmology here in our analysis (Planck15).
        So we have to find a C-M relation which takes into account
        the change in cosmology.
        ****** this is not being used as it's cosmology dependent. We are using
        the relation from this paper https://arxiv.org/pdf/1407.4730.pdf
        ******
        m_pivot = 2e12
        A_mean = 10.14
        B_mean = -0.081
        C_mean = -1.01
        return A_mean * (self.mh/m_pivot)**B_mean * (1 + self.z)**C_mean
    """

# ######################## halo bias b(m, z) ###########################
# Tinker at al. 2010

    def b_nu(self):
        delta = self.delta_halo()
        y = np.log10(delta)
        A = 1.0 + 0.24*y*np.exp(-(4./y)**4)
        aa = 0.44*y - 0.88
        B = 0.183
        b = 1.5
        C = 0.019 + 0.107*y + 0.19*np.exp(-(4./y)**4)
        c = 2.4
        nuu = self.nu()
        dc = 1.686  # neglecting the redshift evolution
        return 1 - (A*nuu**aa/(nuu**aa + dc**aa)) + B*nuu**b + C*nuu**c
