from .utils import *
from .config import *
import numpy as np
from .cosmopower import *
from .pks_and_sigmas import *
import scipy
import time
from multiprocessing import Process
from mcfit import TophatVar
from scipy.interpolate import CubicSpline
import pickle

H_units_conv_factor = {"1/Mpc": 1, "km/s/Mpc": Const.c_km_s}

class Class_szfast(object):
    def __init__(self,
                #lowring=False,  some options if needed
                 **kwargs):
        # some parameters
        # self.xy = xy
        # self.lowring = lowring

        # cosmopower emulators
        self.cp_path_to_cosmopower_organization = path_to_cosmopower_organization + '/'
        self.cp_tt_nn = cp_tt_nn
        self.cp_te_nn = cp_te_nn
        self.cp_ee_nn = cp_ee_nn
        self.cp_pp_nn = cp_pp_nn
        self.cp_pknl_nn  = cp_pknl_nn
        self.cp_pkl_nn = cp_pkl_nn
        self.cp_der_nn = cp_der_nn
        self.cp_da_nn = cp_da_nn
        self.cp_h_nn = cp_h_nn
        self.cp_s8_nn = cp_s8_nn

        if dofftlog_alphas == True:
            self.cp_pkl_fftlog_alphas_nus = cp_pkl_fftlog_alphas_nus
            self.cp_pkl_fftlog_alphas_real_nn  = cp_pkl_fftlog_alphas_real_nn
            self.cp_pkl_fftlog_alphas_imag_nn = cp_pkl_fftlog_alphas_imag_nn

        self.cosmo_model = 'lcdm'
        self.use_Amod = 0
        self.Amod = 0 

        self.cp_lmax = cp_l_max_scalars
        self.cp_ls = np.arange(2,self.cp_lmax+1)

        self.cp_kmax = 50.
        self.cp_kmin = 1e-4
        self.cp_nk = 5000
        self.cp_ndspl_k = 10

        self.cp_predicted_tt_spectrum =np.zeros(self.cp_lmax)
        self.cp_predicted_te_spectrum =np.zeros(self.cp_lmax)
        self.cp_predicted_ee_spectrum =np.zeros(self.cp_lmax)
        self.cp_predicted_pp_spectrum =np.zeros(self.cp_lmax)


        self.cszfast_ldim = 20000 # used for the cls arrays

        self.cszfast_pk_grid_nz = 100 # has to be same as narraySZ, i.e., ndim_redshifts; it is setup hereafter if ndim_redshifts is passed
        self.cszfast_pk_grid_zmax = 5. # current max z of our pk emulators (sept 23)
        self.cszfast_pk_grid_z = np.linspace(0.,self.cszfast_pk_grid_zmax,self.cszfast_pk_grid_nz)
        self.cszfast_pk_grid_ln1pz = np.log(1.+self.cszfast_pk_grid_z)

        self.cszfast_pk_grid_kmin = 1e-4
        self.cszfast_pk_grid_kmax = 50.
        self.cszfast_pk_grid_k = np.geomspace(self.cp_kmin,self.cp_kmax,self.cp_nk)[::self.cp_ndspl_k]
        self.cszfast_pk_grid_lnk = np.log(self.cszfast_pk_grid_k)
        self.cszfast_pk_grid_nk = len(np.geomspace(self.cp_kmin,self.cp_kmax,self.cp_nk)[::self.cp_ndspl_k]) # has to be same as ndimSZ, and the same as dimension of cosmopower pk emulators
        for k,v in kwargs.items():
            # if k == 'ndim_masses':
            #     self.cszfast_pk_grid_nk = v
            #     self.cszfast_pk_grid_k = np.geomspace(self.cszfast_pk_grid_kmin,self.cszfast_pk_grid_kmax,self.cszfast_pk_grid_nk)
            if k == 'ndim_redshifts':
                self.cszfast_pk_grid_nz = v
                self.cszfast_pk_grid_z = np.linspace(0.,self.cszfast_pk_grid_zmax,self.cszfast_pk_grid_nz)
                self.cszfast_pk_grid_ln1pz = np.log(1.+self.cszfast_pk_grid_z)
                self.cszfast_pk_grid_pknl_flat = np.zeros(self.cszfast_pk_grid_nz*self.cszfast_pk_grid_nk)
                self.cszfast_pk_grid_pkl_flat = np.zeros(self.cszfast_pk_grid_nz*self.cszfast_pk_grid_nk)

            if k == 'cosmo_model':
                # print('updating cosmo model')
                cosmo_model_dict = {0: 'lcdm',
                                    1: 'mnu',
                                    2: 'neff',
                                    3: 'wcdm',
                                    4: 'ede',
                                    5: 'mnu-3states'}
                self.cosmo_model = cosmo_model_dict[v]

            if k == 'use_Amod':
                self.use_Amod = v
                self.Amod  = kwargs['Amod']
                # print('self.cosmo_model',self.cosmo_model)
        # print('self.cszfast_pk_grid_nk',self.cszfast_pk_grid_nk)
        # print('self.cszfast_pk_grid_nz',self.cszfast_pk_grid_nz)

        # exit(0)
        self.cp_z_interp = np.linspace(0.,20.,5000)

        self.csz_base = None
        # self.csz_base.compute()


        # z_arr = np.linspace(0.,zmax,nz) # z-array of redshift data [21oct22] oct 26 22: nz = 1000, zmax = 20
        #
        # nk = self.cp_nk
        # ndspl = self.cp_ndspl_k
        # k_arr = np.geomspace(self.cp_kmin,self.cp_kmax,nk)[::ndspl]  # oct 26 22 : (1e-4,50.,5000), jan 10: ndspl


        self.cszfast_zgrid_zmin = 0.
        self.cszfast_zgrid_zmax = 4.
        self.cszfast_zgrid_nz = 250
        self.cszfast_zgrid = np.linspace(self.cszfast_zgrid_zmin,
                                         self.cszfast_zgrid_zmax,
                                         self.cszfast_zgrid_nz)


        self.cszfast_mgrid_mmin = 1e10
        self.cszfast_mgrid_mmax = 1e15
        self.cszfast_mgrid_nm = 50
        self.cszfast_mgrid = np.geomspace(self.cszfast_mgrid_mmin,
                                          self.cszfast_mgrid_mmax,
                                          self.cszfast_mgrid_nm)

        self.cszfast_gas_pressure_xgrid_xmin = 1e-2
        self.cszfast_gas_pressure_xgrid_xmax = 1e2
        self.cszfast_gas_pressure_xgrid_nx = 100
        self.cszfast_gas_pressure_xgrid = np.geomspace(self.cszfast_gas_pressure_xgrid_xmin,
                                                       self.cszfast_gas_pressure_xgrid_xmax,
                                                       self.cszfast_gas_pressure_xgrid_nx)


    def find_As(self,params_cp):
        # params_cp = self.params_cp
        t0 = time.time()

        sigma_8_asked = params_cp["sigma8"]
        # print(params_cp)
        def to_root(ln10_10_As_goal):
            params_cp["ln10^{10}A_s"] = ln10_10_As_goal[0]
            params_dict = {}
            for k,v in params_cp.items():
                params_dict[k]=[v]
            return self.cp_der_nn[self.cosmo_model].ten_to_predictions_np(params_dict)[0][1]-sigma_8_asked

        lnA_s = optimize.root(to_root,
                              x0=3.046,
                              #tol = 1e-10,
                              method="hybr")
        params_cp['ln10^{10}A_s'] = lnA_s.x[0]# .x[0]
        params_cp.pop('sigma8')
        # print("T total in find As",time.time()-t0)#self.t_total)
        # print(params_cp)
        return 1


    def get_H0_from_thetas(self,params_values):
        # print(params_values)
        theta_s_asked = params_values['100*theta_s']
        def fzero(H0_goal):
          params_values['H0'] = H0_goal[0]
          params_dict = {}
          for k,v in params_values.items():
              params_dict[k]=[v]
          # print(params_dict)
          predicted_der_params = self.cp_der_nn[self.cosmo_model].ten_to_predictions_np(params_dict)
          return predicted_der_params[0][0]-theta_s_asked
        sol = optimize.root(fzero,
                          #[40., 99.],
                          x0 = 100.*(3.54*theta_s_asked**2-5.455*theta_s_asked+2.548),
                          #jac=jac,
                          tol = 1e-10,
                          method='hybr')

        params_values.pop('100*theta_s')
        params_values['H0'] = sol.x[0]
        return 1


    def calculate_pkl_fftlog_alphas(self,zpk = 0.,**params_values_dict):
        params_values = params_values_dict.copy()
        params_dict = {}
        for k,v in params_values.items():
            params_dict[k]=[v]
        params_dict['z_pk_save_nonclass'] = [zpk]

        predicted_testing_alphas_creal = self.cp_pkl_fftlog_alphas_real_nn[self.cosmo_model].predictions_np(params_dict)[0]
        predicted_testing_alphas_cimag = self.cp_pkl_fftlog_alphas_imag_nn[self.cosmo_model].predictions_np(params_dict)[0]
        predicted_testing_alphas_cimag = np.append(predicted_testing_alphas_cimag,0.)
        creal = predicted_testing_alphas_creal
        cimag = predicted_testing_alphas_cimag
        Nmax = len(self.cszfast_pk_grid_k)
        cnew = np.zeros(Nmax+1,dtype=complex)
        for i in range(Nmax+1):
            if i<int(Nmax/2):
                cnew[i] = complex(creal[i],cimag[i])
            elif i==int(Nmax/2):
                cnew[i] = complex(creal[i],cimag[i])
            else:
                j = i-int(Nmax/2)
                cnew[i] = complex(creal[::-1][j],-cimag[::-1][j])
        # self.predicted_fftlog_pkl_alphas = cnew
        return cnew

    def get_pkl_reconstructed_from_fftlog(self,zpk = 0.,**params_values_dict):
        #c_n_math = self.predicted_fftlog_pkl_alphas
        c_n_math = self.calculate_pkl_fftlog_alphas(zpk = zpk,**params_values_dict)
        nu_n_math = self.cp_pkl_fftlog_alphas_nus[self.cosmo_model]['arr_0']
        Nmax = int(len(c_n_math)-1)
        term1 = c_n_math[int(Nmax/2)]*self.cszfast_pk_grid_k**(nu_n_math[int(Nmax/2)])
        term2_array = [c_n_math[int(Nmax/2)+i]*self.cszfast_pk_grid_k**(nu_n_math[int(Nmax/2)+i]) for i in range(1, int(Nmax/2)+1)]
        pk_reconstructed = (term1 + 2*np.sum(term2_array,axis=0)).real
        return self.cszfast_pk_grid_k,pk_reconstructed

    def calculate_cmb(self,
                      want_tt=True,
                      want_te=True,
                      want_ee=True,
                      want_pp=1,
                      **params_values_dict):
        params_values = params_values_dict.copy()
        # params_values['ln10^{10}A_s'] = params_values.pop("logA")
        # print('in cmb:',params_values)
        params_dict = {}
        # print('cosmo_model',self.cosmo_model)
        for k,v in params_values.items():
            params_dict[k]=[v]

        if 'm_ncdm' in params_dict.keys():
            if isinstance(params_dict['m_ncdm'][0],str): 
                params_dict['m_ncdm'] =  [float(params_dict['m_ncdm'][0].split(',')[0])]

        # print('params_dict',params_dict)

        if want_tt:
            self.cp_predicted_tt_spectrum = self.cp_tt_nn[self.cosmo_model].ten_to_predictions_np(params_dict)[0]
        if want_te:
            self.cp_predicted_te_spectrum = self.cp_te_nn[self.cosmo_model].predictions_np(params_dict)[0]
        if want_ee:
            self.cp_predicted_ee_spectrum = self.cp_ee_nn[self.cosmo_model].ten_to_predictions_np(params_dict)[0]
        if want_pp:
            self.cp_predicted_pp_spectrum = self.cp_pp_nn[self.cosmo_model].ten_to_predictions_np(params_dict)[0]
        # print('>>> clssy_szfast.py cmb computed')

    def load_cmb_cls_from_file(self,**params_values_dict):
        cls_filename = params_values_dict['cmb_cls_filename']
        with open(cls_filename, 'rb') as handle:
            cmb_cls_loaded = pickle.load(handle)
        nl_cls_file = len(cmb_cls_loaded['ell'])
        cls_ls = cmb_cls_loaded['ell']
        dlfac = cls_ls*(cls_ls+1.)/2./np.pi
        cls_tt = cmb_cls_loaded['tt']*dlfac
        cls_te = cmb_cls_loaded['te']*dlfac
        cls_ee = cmb_cls_loaded['ee']*dlfac
        # cls_pp = cmb_cls_loaded['pp']*dlfac # check the normalization.
        nl_cp_cmb = len(self.cp_predicted_tt_spectrum)
        nl_req = min(nl_cls_file,nl_cp_cmb)

        self.cp_predicted_tt_spectrum[:nl_req-2] = cls_tt[2:nl_req]
        self.cp_predicted_te_spectrum[:nl_req-2] = cls_te[2:nl_req]
        self.cp_predicted_ee_spectrum[:nl_req-2] = cls_ee[2:nl_req]
        # self.cp_predicted_pp_spectrum[:nl_req-2] = cls_pp[2:nl_req]


    def calculate_pkl(self,
                      # cosmo_model = self.cosmo_model,
                      **params_values_dict):
        nz = self.cszfast_pk_grid_nz # number of z-points in redshift data [21oct22] --> set to 80
        zmax = self.cszfast_pk_grid_zmax # max redshift of redshift data [21oct22] --> set to 4 because boltzmannbase.py wants to extrapolate
        z_arr = np.linspace(0.,zmax,nz) # z-array of redshift data [21oct22] oct 26 22: nz = 1000, zmax = 20

        nk = self.cp_nk
        ndspl = self.cp_ndspl_k
        k_arr = np.geomspace(self.cp_kmin,self.cp_kmax,nk)[::ndspl]  # oct 26 22 : (1e-4,50.,5000), jan 10: ndspl


        params_values = params_values_dict.copy()
        # print('in pkl:',params_values)

        params_dict = {}
        for k,v in zip(params_values.keys(),params_values.values()):
            params_dict[k]=[v]

        if 'm_ncdm' in params_dict.keys():
            if isinstance(params_dict['m_ncdm'][0],str): 
                params_dict['m_ncdm'] =  [float(params_dict['m_ncdm'][0].split(',')[0])]



        predicted_pk_spectrum_z = []

        if self.use_Amod:
            # print(">>> using Amod :",self.Amod)
            for zp in z_arr:
                params_dict_pp = params_dict.copy()
                params_dict_pp['z_pk_save_nonclass'] = [zp]
                pkl_p = self.cp_pkl_nn[self.cosmo_model].predictions_np(params_dict_pp)[0]
                pknl_p = self.cp_pknl_nn[self.cosmo_model].predictions_np(params_dict_pp)[0]
                pk_ae  = pkl_p + self.Amod*(pknl_p-pkl_p)
                predicted_pk_spectrum_z.append(pk_ae)
        else:
            for zp in z_arr:
                params_dict_pp = params_dict.copy()
                params_dict_pp['z_pk_save_nonclass'] = [zp]
                predicted_pk_spectrum_z.append(self.cp_pkl_nn[self.cosmo_model].predictions_np(params_dict_pp)[0])

        predicted_pk_spectrum = np.asarray(predicted_pk_spectrum_z)

        # weird scaling to get rid off
        # scaling factor for the pk emulator:
        ls = np.arange(2,self.cp_nk+2)[::ndspl] # jan 10 ndspl
        dls = ls*(ls+1.)/2./np.pi
        pk = 10.**predicted_pk_spectrum
        pk_re =  ((dls)**-1*pk)
        pk_re = np.transpose(pk_re)

        # print(z_arr,np.log(k_arr),np.log(pk_re))
        # print(np.log(pk_re).min())
        # print(np.log(pk_re).max())
        # print(np.shape(np.log(pk_re)))
        # print('coordinates')
        # # z_coords, k_coords = np.meshgrid(z_arr, np.log(k_arr), indexing='ij')
        # # z_coords, k_coords = z_coords.flatten(), k_coords.flatten()
        # z_coords, k_coords = z_arr,np.log(k_arr)
        # print(len(z_coords),len(k_coords))
        # #
        # # print(np.shape(list(zip(z_arr,np.log(k_arr)))))
        # pk_values = np.log(pk_re).T
        # # pk_values = pk_values.ravel()
        # print(np.shape(pk_values),pk_values)
        #
        # self.pkl_linearnd_interp = LinearNDInterpolator(np.asarray(z_coords),np.asarray(k_coords), pk_values)
        # exit(0)
        # # self.pkl_cloughtocher_interp = CloughTocher2DInterpolator(list(zip(z_arr,np.log(k_arr))), np.log(pk_re).T)
        # exit(0)
        # self.pkl_cloughtocher_interp = CloughTocher2DInterpolator(list(zip(z_arr,np.log(k_arr))), pk_value)
        # is this line making things slower?
        self.pkl_interp = PowerSpectrumInterpolator(z_arr,k_arr,np.log(pk_re).T,logP=True)
        # self.pkl_interp = None
        self.cszfast_pk_grid_pk = pk_re
        self.cszfast_pk_grid_pkl_flat = pk_re.flatten()
        return pk_re, k_arr, z_arr


    def calculate_sigma(self,
                        # cosmo_model = self.cosmo_model,
                        # z_asked = None,
                        # r_asked = None,
                        **params_values_dict):
        params_values = params_values_dict.copy()
        k = self.cszfast_pk_grid_k
        # self.cszfast_pk_grid_z
        # print(self.cszfast_pk_grid_pk,np.shape(self.cszfast_pk_grid_pk))
        P = self.cszfast_pk_grid_pk
        var = P.copy()
        dvar = P.copy()
        for iz,zp in enumerate(self.cszfast_pk_grid_z):
            R, var[:,iz] = TophatVar(k, lowring=True)(P[:,iz], extrap=True)
            # dvar[:,iz] = np.gradient(var[:,iz], np.log(R))
            # if params_values_dict['sigma_derivative'] == 1: ## need more points here !!
            #     rds,dvar[:,iz] =  TophatVar(k,lowring=True,deriv=1)(P[:,iz]*k,extrap=True)
            # else:
            #     dvar[:,iz] = np.gradient(var[:,iz], R)

            dvar[:,iz] = np.gradient(var[:,iz], R)
            # from inigo: R_vec,self.dvar = TophatVar(k,lowring=True,deriv=1)(P[:,iz]*k,extrap=True)
            # from inigo: self.dsigma_vec = self.dvar/(2.*self.sigma_vec)
            # print('R:',R,np.shape(R))
            # varR = CubicSpline(R, var[:,iz])
            # print(zp,np.sqrt(varR(8)))
        # print(params_values)
        # print('in sigma:',params_values)
        # h = params_values['H0']/100.
        # var = var.T
        # dvar = dvar.T

        self.cszfast_pk_grid_lnr = np.log(R)
        self.cszfast_pk_grid_sigma2 = var
        self.cszfast_pk_grid_sigma2_flat = var.flatten()
        self.cszfast_pk_grid_lnsigma2_flat = 0.5*np.log(var.flatten())
        # self.cszfast_pk_grid_lnsigma2_flat = self.cszfast_pk_grid_lnsigma2_flat.T
        self.cszfast_pk_grid_dsigma2 = dvar
        self.cszfast_pk_grid_dsigma2_flat = dvar.flatten()
        # if z_asked != None and r_asked != None:
        #     print(z_asked[0],r_asked[0])
        #     return z_asked, r_asked
        # else:
        #     return 0
        return 0


    def calculate_sigma8_and_der(self,
                         # cosmo_model = self.cosmo_model,
                         **params_values_dict):
        params_values = params_values_dict.copy()
        # print('in pkl:',params_values)

        params_dict = {}
        for k,v in zip(params_values.keys(),params_values.values()):
            params_dict[k]=[v]

        if 'm_ncdm' in params_dict.keys():
            if isinstance(params_dict['m_ncdm'][0],str): 
                params_dict['m_ncdm'] =  [float(params_dict['m_ncdm'][0].split(',')[0])]

        self.cp_predicted_der = self.cp_der_nn[self.cosmo_model].ten_to_predictions_np(params_dict)[0]
        self.sigma8 = self.cp_predicted_der[1]
        return 0


    def calculate_sigma8_at_z(self,
                             # cosmo_model = self.cosmo_model,
                             **params_values_dict):
        params_values = params_values_dict.copy()
        # print('in pkl:',params_values)

        params_dict = {}
        for k,v in zip(params_values.keys(),params_values.values()):
            params_dict[k]=[v]

        if 'm_ncdm' in params_dict.keys():
            if isinstance(params_dict['m_ncdm'][0],str): 
                params_dict['m_ncdm'] =  [float(params_dict['m_ncdm'][0].split(',')[0])]


        s8z  = self.cp_s8_nn[self.cosmo_model].predictions_np(params_dict)
        # print(self.s8z)
        self.s8z_interp = scipy.interpolate.interp1d(
                                                    np.linspace(0.,20.,5000),
                                                    s8z[0],
                                                    kind='linear',
                                                    axis=-1,
                                                    copy=True,
                                                    bounds_error=None,
                                                    fill_value=np.nan,
                                                    assume_sorted=False)

    def calculate_pknl(self,
                      # cosmo_model = self.cosmo_model,
                      **params_values_dict):
        nz = self.cszfast_pk_grid_nz # number of z-points in redshift data [21oct22] --> set to 80
        zmax = self.cszfast_pk_grid_zmax # max redshift of redshift data [21oct22] --> set to 4 because boltzmannbase.py wants to extrapolate
        z_arr = np.linspace(0.,zmax,nz) # z-array of redshift data [21oct22] oct 26 22: nz = 1000, zmax = 20

        nk = self.cp_nk
        ndspl = self.cp_ndspl_k
        k_arr = np.geomspace(self.cp_kmin,self.cp_kmax,nk)[::ndspl]  # oct 26 22 : (1e-4,50.,5000), jan 10: ndspl


        params_values = params_values_dict.copy()
        # print('in pknl:',params_values)

        params_dict = {}
        for k,v in zip(params_values.keys(),params_values.values()):
            params_dict[k]=[v]

        if 'm_ncdm' in params_dict.keys():
            if isinstance(params_dict['m_ncdm'][0],str): 
                params_dict['m_ncdm'] =  [float(params_dict['m_ncdm'][0].split(',')[0])]


        predicted_pk_spectrum_z = []

        for zp in z_arr:
            params_dict_pp = params_dict.copy()
            params_dict_pp['z_pk_save_nonclass'] = [zp]
            predicted_pk_spectrum_z.append(self.cp_pknl_nn[self.cosmo_model].predictions_np(params_dict_pp)[0])

        predicted_pk_spectrum = np.asarray(predicted_pk_spectrum_z)

        # weird scaling to get rid off
        # scaling factor for the pk emulator:
        ls = np.arange(2,self.cp_nk+2)[::ndspl] # jan 10 ndspl
        dls = ls*(ls+1.)/2./np.pi
        pk = 10.**predicted_pk_spectrum
        pk_re =  ((dls)**-1*pk)
        pk_re = np.transpose(pk_re)

        # print(z_arr,np.log(k_arr),np.log(pk_re))
        # exit(0)
        # self.pknl_linearnd_interp = LinearNDInterpolator(list(zip(z_arr,np.log(k_arr))), pk_re.T)
        # self.pknl_cloughtocher_interp = CloughTocher2DInterpolator(list(zip(z_arr,np.log(k_arr))), pk_re.T)

        # is this line making things slower?
        self.pknl_interp = PowerSpectrumInterpolator(z_arr,k_arr,np.log(pk_re).T,logP=True)
        # self.pknl_interp = None
        self.cszfast_pk_grid_pknl = pk_re
        self.cszfast_pk_grid_pknl_flat = pk_re.flatten()
        return pk_re, k_arr, z_arr


    def calculate_hubble(self,
                                 # cosmo_model = self.cosmo_model,
                                 **params_values_dict):
        params_values = params_values_dict.copy()

        params_dict = {}
        for k,v in zip(params_values.keys(),params_values.values()):
            params_dict[k]=[v]

        if 'm_ncdm' in params_dict.keys():
            if isinstance(params_dict['m_ncdm'][0],str): 
                params_dict['m_ncdm'] =  [float(params_dict['m_ncdm'][0].split(',')[0])]



        self.cp_predicted_hubble = self.cp_h_nn[self.cosmo_model].ten_to_predictions_np(params_dict)[0]
        # print(self.cp_predicted_hubble)
        # z_interp =
        self.hz_interp = scipy.interpolate.interp1d(
                                            self.cp_z_interp,
                                            self.cp_predicted_hubble,
                                            kind='linear',
                                            axis=-1,
                                            copy=True,
                                            bounds_error=None,
                                            fill_value=np.nan,
                                            assume_sorted=False)

    def calculate_chi(self,
                      # cosmo_model = self.cosmo_model,
                      **params_values_dict):
        # def test_integrand_func(x, alpha, beta, i, j, k, l):
        #     return x * alpha * beta + i * j * k
        #
        #
        # grid = np.random.random((5000, 2))
        #
        # res, err = pyquad.quad_grid(test_integrand_func, 0, 1, grid, (1.0, 1.0, 1.0, 1.0))
        #
        # print(res,err)
        # def integrand_chi(z,alpha, beta, i, j, k, l):
        #     z = z-1.
        #     return 1./self.get_Hubble(z)
        # zmax = 1.
        # grid = np.random.random((10000000, 2))
        # chiz,err = pyquad.quad_grid(integrand_chi, 1., 1.+zmax, grid, (1.0, 1.0, 1.0, 1.0))
        # print(chiz)
        params_values = params_values_dict.copy()

        params_dict = {}
        for k,v in zip(params_values.keys(),params_values.values()):
            params_dict[k]=[v]

        if 'm_ncdm' in params_dict.keys():
            if isinstance(params_dict['m_ncdm'][0],str): 
                params_dict['m_ncdm'] =  [float(params_dict['m_ncdm'][0].split(',')[0])]


        self.cp_predicted_da  = self.cp_da_nn[self.cosmo_model].predictions_np(params_dict)[0]
        self.chi_interp = scipy.interpolate.interp1d(
                                                    self.cp_z_interp,
                                                    self.cp_predicted_da*(1.+self.cp_z_interp),
                                                    kind='linear',
                                                    axis=-1,
                                                    copy=True,
                                                    bounds_error=None,
                                                    fill_value=np.nan,
                                                    assume_sorted=False)

    def get_cmb_cls(self,ell_factor=True,Tcmb_uk = Tcmb_uk):
        cls = {}
        cls['ell'] = np.arange(self.cszfast_ldim)
        cls['tt'] = np.zeros(self.cszfast_ldim)
        cls['te'] = np.zeros(self.cszfast_ldim)
        cls['ee'] = np.zeros(self.cszfast_ldim)
        cls['pp'] = np.zeros(self.cszfast_ldim)
        cls['bb'] = np.zeros(self.cszfast_ldim)
        cls['tt'][2:self.cp_lmax+1] = (Tcmb_uk)**2.*self.cp_predicted_tt_spectrum.copy()
        cls['te'][2:self.cp_lmax+1] = (Tcmb_uk)**2.*self.cp_predicted_te_spectrum.copy()
        cls['ee'][2:self.cp_lmax+1] = (Tcmb_uk)**2.*self.cp_predicted_ee_spectrum.copy()
        cls['pp'][2:self.cp_lmax+1] = self.cp_predicted_pp_spectrum.copy()/4. ## this is clkk... works for so lensinglite lkl
        # cls['bb'][2:self.cp_lmax+1] = self.cp_predicted_pp_spectrum.copy()/4. ## this is clkk... works for so lensinglite lkl
        # print('doing gets')
        # For planck likelihood:
        # lcp = np.asarray(cls['ell'][2:nl+2])
        # cls['pp'][2:nl+2] = self.pp_spectra[0].copy()/(lcp*(lcp+1.))**2.
        # cls['pp'][2:nl+2] *= (lcp*(lcp+1.))**2./2./np.pi



        if ell_factor==False:
            fac_l = np.zeros(self.cszfast_ldim)
            fac_l[2:self.cp_lmax+1] = 1./(cls['ell'][2:self.cp_lmax+1]*(cls['ell'][2:self.cp_lmax+1]+1.)/2./np.pi)
            cls['tt'][2:self.cp_lmax+1] *= fac_l[2:self.cp_lmax+1]
            cls['te'][2:self.cp_lmax+1] *= fac_l[2:self.cp_lmax+1]
            cls['ee'][2:self.cp_lmax+1] *= fac_l[2:self.cp_lmax+1]
            # cls['bb'] = np.zeros(self.cszfast_ldim)
        return cls


    def get_pknl_at_k_and_z(self,k_asked,z_asked):
    # def get_pkl_at_k_and_z(self,k_asked,z_asked,method = 'cloughtocher'):
        # if method == 'linear':
        #     pk = self.pkl_linearnd_interp(z_asked,np.log(k_asked))
        # elif method == 'cloughtocher':
        #     pk = self.pkl_cloughtocher_interp(z_asked,np.log(k_asked))
        # return np.exp(pk)
        return self.pknl_interp.P(z_asked,k_asked)


    # def get_pkl_at_k_and_z(self,k_asked,z_asked,method = 'cloughtocher'):
    def get_pkl_at_k_and_z(self,k_asked,z_asked):
        # if method == 'linear':
        #     pk = self.pknl_linearnd_interp(z_asked,np.log(k_asked))
        # elif method == 'cloughtocher':
        #     pk = self.pknl_cloughtocher_interp(z_asked,np.log(k_asked))
        # return np.exp(pk)
        return self.pkl_interp.P(z_asked,k_asked)

    # function used to overwrite the classy function in fast mode.
    def get_sigma_at_r_and_z(self,r_asked,z_asked):
        # if method == 'linear':
        #     pk = self.pknl_linearnd_interp(z_asked,np.log(k_asked))
        # elif method == 'cloughtocher':
        #     pk = self.pknl_cloughtocher_interp(z_asked,np.log(k_asked))
        # return np.exp(pk)
        k = self.cszfast_pk_grid_k
        P_at_z = self.get_pkl_at_k_and_z(k,z_asked)
        R,var = TophatVar(k, lowring=True)(P_at_z, extrap=True)
        varR = CubicSpline(R, var)
        sigma_at_r_and_z = np.sqrt(varR(r_asked))
        return sigma_at_r_and_z



    def get_hubble(self, z,units="1/Mpc"):
        return np.array(self.hz_interp(z)*H_units_conv_factor[units])

    def get_chi(self, z):
        return np.array(self.chi_interp(z))

    def get_gas_pressure_profile_x(self,z,m,x):
        return 0#np.vectorize(self.csz_base.get_pressure_P_over_P_delta_at_x_M_z_b12_200c)(x,m,z)


    # def get_gas_pressure_profile_x_parallel(self,index_z,param_values_array,**kwargs):
    #     # zp = self.cszfast_zgrid[iz]
    #     # x = self.get_gas_pressure_profile_x(zp,3e13,self.cszfast_gas_pressure_xgrid)
    #     x = 1
    #     return x


    def tabulate_gas_pressure_profile_k(self):
        z_asked,m_asked,x_asked = 0.2,3e14,np.geomspace(1e-3,1e2,500)
        start = time.time()
        px = self.get_gas_pressure_profile_x(z_asked,m_asked,x_asked)
        end = time.time()
        # print(px)
        # print('end tabuulate pressure profile:',end-start)
        # print('grid tabulate pressure profile')
        start = time.time()
        # px = self.get_gas_pressure_profile_x(self.cszfast_zgrid,self.cszfast_mgrid,self.cszfast_gas_pressure_xgrid)
        # for zp in self.cszfast_zgrid:
        #     for mp in self.cszfast_mgrid:
        #         zp = mp
        #         # px = self.get_gas_pressure_profile_x(zp,mp,self.cszfast_gas_pressure_xgrid)
        # zp,mp = 0.1,2e14
        px = self.get_gas_pressure_profile_x(self.cszfast_zgrid[:,None,None],
                                             self.cszfast_mgrid[None,:,None],
                                             self.cszfast_gas_pressure_xgrid[None,None,:])

        # fn=functools.partial(get_gas_pressure_profile_x_parallel,
        #                      param_values_array=None)
        # pool = multiprocessing.Pool()
        # results = pool.map(fn,range(self.cszfast_zgrid_nz))
        # pool.close()
        #
        # def task(iz):
        #     zp = self.cszfast_zgrid[iz]
        #     x = self.get_gas_pressure_profile_x(zp,m_asked,self.cszfast_gas_pressure_xgrid)
        # processes = [Process(target=task, args=(i,)) for i in range(self.cszfast_zgrid_nz)]
        #         # start all processes
        # for process in processes:
        #     process.start()
        # # wait for all processes to complete
        # for process in processes:
        #     process.join()
        # # report that all tasks are completed
        # print('Done', flush=True)
        end = time.time()

        # print(px)
        # print('end grid tabuulate pressure profile:',end-start)
    #
    # def get_gas_pressure_profile_x_parallel(index_z,param_values_array,**kwargs):
    #     # zp = self.cszfast_zgrid[iz]
    #     # x = self.get_gas_pressure_profile_x(zp,3e13,self.cszfast_gas_pressure_xgrid)
    #     x = 1
    #     return x


    def get_sigma8_at_z(self,z):
        return self.s8z_interp(z)


    def rs_drag(self):
        try:
            return self.cp_predicted_der[13]
        except AttributeError:
            return 0
    #################################
    # gives an estimation of f(z)*sigma8(z) at the scale of 8 h/Mpc, computed as (d sigma8/d ln a)
    # Now used by cobaya wrapper.
    def get_effective_f_sigma8(self, z, z_step=0.1,params_values_dict={}):
        """
        effective_f_sigma8(z)

        Returns the time derivative of sigma8(z) computed as (d sigma8/d ln a)

        Parameters
        ----------
        z : float
                Desired redshift
        z_step : float
                Default step used for the numerical two-sided derivative. For z < z_step the step is reduced progressively down to z_step/10 while sticking to a double-sided derivative. For z< z_step/10 a single-sided derivative is used instead.

        Returns
        -------
        (d ln sigma8/d ln a)(z) (dimensionless)
        """

        s8z_interp =  self.s8z_interp
        # we need d sigma8/d ln a = - (d sigma8/dz)*(1+z)

        # if possible, use two-sided derivative with default value of z_step
        if z >= z_step:
            result = (s8z_interp(z-z_step)-s8z_interp(z+z_step))/(2.*z_step)*(1+z)
            # return (s8z_interp(z-z_step)-s8z_interp(z+z_step))/(2.*z_step)*(1+z)
        else:
            # if z is between z_step/10 and z_step, reduce z_step to z, and then stick to two-sided derivative
            if (z > z_step/10.):
                z_step = z
                result = (s8z_interp(z-z_step)-s8z_interp(z+z_step))/(2.*z_step)*(1+z)
                # return (s8z_interp(z-z_step)-s8z_interp(z+z_step))/(2.*z_step)*(1+z)
            # if z is between 0 and z_step/10, use single-sided derivative with z_step/10
            else:
                z_step /=10
                result = (s8z_interp(z)-s8z_interp(z+z_step))/z_step*(1+z)
                # return (s8z_interp(z)-s8z_interp(z+z_step))/z_step*(1+z)
        # print('fsigma8 result : ',result)
        return result
