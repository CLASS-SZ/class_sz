import numpy as np
import os
from scipy.ndimage.interpolation import shift





SZ_data_directory  = '/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/sz_auxiliary_files/'
planck_SZ_datafile = 'SZ_and_fg_models.txt'

#fixed trispectrum files
#we use a fixed trispectrum given in the files:
fiducial_file_covmat	 = 'cobaya_class_sz_likelihoods/cobaya_reference_trispectrum/tSZ_trispectrum_ref.txt'
fiducial_file_cell       = 'cobaya_class_sz_likelihoods/cobaya_reference_trispectrum/tSZ_c_ell_ref.txt'
fiducial_file_params = 'cobaya_class_sz_likelihoods/cobaya_reference_trispectrum/tSZ_param_ref.txt'


#By default we use fixed trispectrum for the Planck likelihood
class SOSZPowerSpectrumLikelihood(object):

    def __init__(self):
        #Planck 2015 y-map 2pt likelihood
        self.datafile = planck_SZ_datafile
        self.data_directory = SZ_data_directory

        D = np.loadtxt(os.path.join(self.data_directory, self.datafile))

        #fill arrays with Planck tSZ and FG data
        self.ell_plc = D[:,0]
        self.y2AndFg_plc = D[:,1]
        self.sigma_tot_plc = D[:,2]
        self.A_CIB_MODEL = D[:,3]
        self.A_RS_MODEL = D[:,4]
        self.A_IR_MODEL = D[:,5]
        self.A_CN_MODEL = D[:,6]

        #sky_fraction of Planck y-map (communicated by Barbara Commis)
        self.f_sky = 0.47

        #number of data points (multipole bins)
        self.num_points = np.shape(self.ell_plc)[0]


        self.fid_values_exist = False
        self.fiducial_file_params = fiducial_file_params
        self.fiducial_file_cell = fiducial_file_cell
        self.fiducial_file_covmat = fiducial_file_covmat
        #If reference trispectrun is available -> read in the values
        if os.path.exists(os.path.join(self.data_directory, self.fiducial_file_params)):
            self.fid_values_exist = True
            #read-in the fiducial data
            D = np.loadtxt(os.path.join(self.data_directory, self.fiducial_file_cell))
            self.ell = D[:,0]
            self.data = D[:,1]

            #check that the reference tSZ data is computed at same multipoles as planck data
            try:
                if(np.any(self.ell-self.ell_plc)):
                    print("[reading ref tSZ files] the reference multipoles do not match the planck data.")
                    exit(0)
            except ValueError:
                print("[reading ref tSZ files] the reference multipoles do not match the planck data.")
                exit(0)

            #compute number of multipoles in each bin
            #[useful when we need to comput ethe gaussian sample variance]
            #[not used in the original planck likelihood as sigma_G is tabulated]
            multipoles_bin_center = self.ell
            mp_shift = shift(multipoles_bin_center, -1, cval=np.NaN)
            diff_ell_shift =  mp_shift - multipoles_bin_center
            diff_ell_up = diff_ell_shift
            diff_ell_up[len(diff_ell_shift)-1] = diff_ell_shift[len(diff_ell_shift)-2]
            mp_shift = shift(multipoles_bin_center, 1, cval=np.NaN)
            diff_ell_shift =  -mp_shift + multipoles_bin_center
            diff_ell_down = diff_ell_shift
            diff_ell_down[0] = diff_ell_shift[1]
            nell =  (diff_ell_up + diff_ell_down)/2.
            #Compute sigma_G^2
            #self.cvg =1./self.f_sky*(np.asarray(self.data)*np.sqrt(2./(2.*np.asarray(self.ell)+1.)))**2./nell

            #For the planck likelihood we read-in the values of sigma_G:
            self.cvg = np.asarray(self.sigma_tot_plc)**2.
            self.cvg = np.diag(self.cvg)
            self.covmat = np.asarray(np.loadtxt(os.path.join(self.data_directory, self.fiducial_file_covmat)))
            #Add sigm_NG (trispectrum strored in self.covmat) to sigma_G (stored in self.cvg)
            self.covmat = self.covmat + self.cvg
            self.inv_covmat = np.linalg.inv(self.covmat)
            self.det_covmat = np.linalg.det(self.covmat)
            print('[reading ref tSZ files] read-in completed.')
        else:
            print('[reading ref tSZ files] reference trispectrum unavailable -> creating reference trispectrum.')



    def __call__(self, _theory={}):
        #collect current SZ quantities/info
        D = _theory.get_Cl_sz()
        params = D[0]
        ell = np.asarray(list(D[1].values()))
        cl = np.asarray(list(D[2].values()))
        trispectrum = list(D[3].values())
        output_string = params.get('output')

        print(cl)


        if self.fid_values_exist is True:
            if 'tSZ_Trispectrum' in output_string:
                print("Make sure you do not require 'tSZ_Trispectrum' in the 'output' field of the input param file, as this will speed-up the calculations.")
                exit(0)
            inv_covmat = self.inv_covmat
            det_covmat = self.det_covmat

            print("det = ")
            print(det_covmat)

            #foreground residua amplitudes [TBD]
            A_CIB = 0.
            A_RS = 0.
            A_IR = 0.
            # A_CN amplitude is set by looking at
            # SZ_and_fg_models-high_ell.txt at high multipole (ell = 2742),
            # where the correlated noise largely dominates over the other
            # components (see Bolliet et al. 1712.00788).
            A_CN = 0.9033

            # define array for difference between theory and observations
            difference = np.ndarray(self.num_points, 'float64')
            for i in range(self.num_points):
                difference[i] = self.y2AndFg_plc[i] - (cl[i] + A_CIB*self.A_CIB_MODEL[i] + A_RS*self.A_RS_MODEL[i] + A_IR*self.A_IR_MODEL[i] + A_CN*self.A_CN_MODEL[i])

            # chisquare
            chi2 = np.dot(difference, np.dot(inv_covmat, difference))
            #since trispectrum and sigma_G are fixed we do not need the determinant term as it is a constant.
            return  - 0.5*chi2 #-0.5*np.log(det_covmat)



        if self.fid_values_exist is False:
            if 'tSZ_Trispectrum' not in output_string:
                print("[writting ref tSZ files] In order to compute the reference trispectrum, You need to require 'tSZ_Trispectrum' in  'output' field of the input param file")
                exit(0)


            try:
                if(np.any(ell-self.ell_plc)):
                    print("[writting ref tSZ files] the reference multipoles do not match the planck data.")
                    exit(0)
            except ValueError:
                print("[writting ref tSZ files] the reference multipoles do not match the planck data.")
                exit(0)

            #write parameters for fixed reference trispectrum
            fid_file = open(os.path.join(self.data_directory, self.fiducial_file_params), 'w')
            fid_file.write('# Parameters for the reference trispectrum\n\n')
            for k, v in params.items():
                fid_file.write(str(k) + ' = '+ str(v) + '\n\n')

            #write fixed reference trispectrum
            T_l_l_prime = np.zeros((self.num_points, self.num_points), 'float64')
            for j in range (self.num_points):
                for i in range (self.num_points):
                    if i<=j:
                        t_l_l_prime = (ell[i]*(ell[i]+1.)/(2.*np.pi))*(ell[j]*(ell[j]+1.)/(2.*np.pi))*trispectrum[j][i]
                        T_l_l_prime[j][i] = t_l_l_prime/self.f_sky
                        T_l_l_prime[i][j] = T_l_l_prime[j][i]
            fid_file = open(os.path.join(self.data_directory, self.fiducial_file_covmat), 'w')
            np.savetxt(os.path.join(self.data_directory, self.fiducial_file_covmat),T_l_l_prime,fmt='%1.4e')

            #write cl's corresponding to fixed reference trispectrum
            fid_file = open(os.path.join(self.data_directory, self.fiducial_file_cell), 'w')
            np.savetxt(os.path.join(self.data_directory, self.fiducial_file_cell),
                                    np.c_[ell,cl],
                       fmt='%1.4e')
            print("[writting ref tSZ files] A reference trispectrum has been written.")
            print("[writting ref tSZ files] You can now run chains using this reference trispectrum for non-gaussian errors.")
            print("[writting ref tSZ files] To do so, make sure you do not require 'tSZ_Trispectrum' in the 'output' field input param file, as this will speed-up the calculations.")
            exit(0)
