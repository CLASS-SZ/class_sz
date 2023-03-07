"""
.. module:: classy
    :synopsis: Python wrapper around CLASS
.. moduleauthor:: Karim Benabed <benabed@iap.fr>
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
.. moduleauthor:: Julien Lesgourgues <lesgourg@cern.ch>

This module defines a class called Class. It is used with Monte Python to
extract cosmological parameters.

# JL 14.06.2017: TODO: check whether we should free somewhere the allocated fc.filename and titles, data (4 times)

"""
from math import exp,log
from collections import defaultdict #BB: added for class_sz
import numpy as np
cimport numpy as np
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
import cython
cimport cython
from scipy import interpolate
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline as _spline

# testing stuffsL: BB
# import tensorflow as tf
# import cosmopower
# path_to_cosmopower_organization = '/Users/boris/Work/CLASS-SZ/SO-SZ/cosmopower-organization'
# path_to_emulators = path_to_cosmopower_organization + 'lcdm/'
# from mcfit import P2xi
import time
from classy_szfast import classy_szfast



# Nils : Added for python 3.x and python 2.x compatibility
import sys
def viewdictitems(d):
    if sys.version_info >= (3,0):
        return d.items()
    else:
        return d.viewitems()

ctypedef np.float_t DTYPE_t
ctypedef np.int_t DTYPE_i

# Import the .pxd containing definitions
from cclassy cimport *

DEF _MAXTITLESTRINGLENGTH_ = 8000

__version__ = _VERSION_.decode("utf-8")

# Implement a specific Exception (this might not be optimally designed, nor
# even acceptable for python standards. It, however, does the job).
# The idea is to raise either an AttributeError if the problem happened while
# reading the parameters (in the normal Class, this would just return a line in
# the unused_parameters file), or a NameError in other cases. This allows
# MontePython to handle things differently.
class CosmoError(Exception):
    def __init__(self, message=""):
        self.message = message.decode() if isinstance(message,bytes) else message

    def __str__(self):
        return '\n\nError in Class: ' + self.message


class CosmoSevereError(CosmoError):
    """
    Raised when Class failed to understand one or more input parameters.

    This case would not raise any problem in Class default behaviour. However,
    for parameter extraction, one has to be sure that all input parameters were
    understood, otherwise the wrong cosmological model would be selected.
    """
    pass


class CosmoComputationError(CosmoError):
    """
    Raised when Class could not compute the cosmology at this point.

    This will be caught by the parameter extraction code to give an extremely
    unlikely value to this point
    """
    pass


cdef class Class:
    """
    Class wrapping, creates the glue between C and python

    The actual Class wrapping, the only class we will call from MontePython
    (indeed the only one we will import, with the command:
    from classy import Class

    """
    # List of used structures, defined in the header file. They have to be
    # "cdefined", because they correspond to C structures
    cdef precision pr
    cdef background ba
    cdef thermo th
    cdef perturbs pt
    cdef primordial pm
    cdef nonlinear nl
    cdef transfers tr
    cdef spectra sp
    cdef output op
    cdef lensing le
    cdef tszspectrum tsz  #BB: added for class_sz
    cdef szcount csz  #BB: added for class_sz
    cdef file_content fc

    cpdef int computed # Flag to see if classy has already computed with the given pars
    cpdef int allocated # Flag to see if classy structs are allocated already
    cpdef object _pars # Dictionary of the parameters
    cpdef object ncp   # Keeps track of the structures initialized, in view of cleaning.


    cpdef double sigma8_fast
    cpdef object class_szfast

    # Defining two new properties to recover, respectively, the parameters used
    # or the age (set after computation). Follow this syntax if you want to
    # access other quantities. Alternatively, you can also define a method, and
    # call it (see _T_cmb method, at the very bottom).
    property pars:
        def __get__(self):
            return self._pars
    property state:
        def __get__(self):
            return True
    property Omega_nu:
        def __get__(self):
            return self.ba.Omega0_ncdm_tot
    property nonlinear_method:
        def __get__(self):
            return self.nl.method

#    def set_default(self):
#        _pars = {
#            "output":"tCl mPk",}
#        self.set(**_pars)
# BB: modified for class_sz
    def set_default(self):
        _pars = {
            "output":"tSZ",}
        self.set(**_pars)

    def __cinit__(self, default=False):
        cpdef char* dumc
        self.allocated = False
        self.computed = False
        self._pars = {}
        self.fc.size=0
        self.fc.filename = <char*>malloc(sizeof(char)*30)
        assert(self.fc.filename!=NULL)
        dumc = "NOFILE"
        sprintf(self.fc.filename,"%s",dumc)
        self.ncp = set()
        if default: self.set_default()

    def __dealloc__(self):
        if self.allocated:
          self.struct_cleanup()
        self.empty()
        # Reset all the fc to zero if its not already done
        if self.fc.size !=0:
            self.fc.size=0
            free(self.fc.name)
            free(self.fc.value)
            free(self.fc.read)
            free(self.fc.filename)

    # Set up the dictionary
    def set(self,*pars,**kars):
        oldpars = self._pars.copy()
        if len(pars)==1:
            self._pars.update(dict(pars[0]))
        elif len(pars)!=0:
            raise CosmoSevereError("bad call")
        self._pars.update(kars)
        if viewdictitems(self._pars) <= viewdictitems(oldpars):
          return # Don't change the computed states, if the new dict was already contained in the previous dict
        self.computed=False
        return True

    def empty(self):
        self._pars = {}
        self.computed = False

    # Create an equivalent of the parameter file. Non specified values will be
    # taken at their default (in Class)
    def _fillparfile(self):
        cdef char* dumc

        if self.fc.size!=0:
            free(self.fc.name)
            free(self.fc.value)
            free(self.fc.read)
        self.fc.size = len(self._pars)
        self.fc.name = <FileArg*> malloc(sizeof(FileArg)*len(self._pars))
        assert(self.fc.name!=NULL)

        self.fc.value = <FileArg*> malloc(sizeof(FileArg)*len(self._pars))
        assert(self.fc.value!=NULL)

        self.fc.read = <short*> malloc(sizeof(short)*len(self._pars))
        assert(self.fc.read!=NULL)

        # fill parameter file
        i = 0
        for kk in self._pars:

            dumcp = kk.encode()
            dumc = dumcp
            sprintf(self.fc.name[i],"%s",dumc)
            dumcp = str(self._pars[kk]).encode()
            dumc = dumcp
            sprintf(self.fc.value[i],"%s",dumc)
            self.fc.read[i] = _FALSE_
            i+=1

    # Called at the end of a run, to free memory
    def struct_cleanup(self):
        # print('cleaning up memoryyyy!!!')
        # print('level:',self.ncp)
        if(self.allocated != True):
          return
        if self.tsz.use_class_sz_fast_mode == 1:
            szcounts_free(&self.csz,&self.tsz)
            szpowerspectrum_free(&self.tsz)
        if "szcount" in self.ncp: #BB: added for class_sz
            szcounts_free(&self.csz,&self.tsz)
        if "szpowerspectrum" in self.ncp:  #BB: added for class_sz
            szpowerspectrum_free(&self.tsz)
        if "lensing" in self.ncp:
            lensing_free(&self.le)
        if "spectra" in self.ncp:
            spectra_free(&self.sp)
        if "transfer" in self.ncp:
            transfer_free(&self.tr)
        if "nonlinear" in self.ncp:
            nonlinear_free(&self.nl)
        if "primordial" in self.ncp:
            primordial_free(&self.pm)
        if "perturb" in self.ncp:
            perturb_free(&self.pt)
        if "thermodynamics" in self.ncp:
            thermodynamics_free(&self.th)
        if "background" in self.ncp:
            background_free(&self.ba)
        self.allocated = False
        self.computed = False

    def _check_task_dependency(self, level):
        """
        Fill the level list with all the needed modules

        .. warning::

            the ordering of modules is obviously dependent on CLASS module order
            in the main.c file. This has to be updated in case of a change to
            this file.

        Parameters
        ----------

        level : list
            list of strings, containing initially only the last module required.
            For instance, to recover all the modules, the input should be
            ['lensing']

        """
        if "szcount" in level:  #BB: added for class_sz
            if "szpowerspectrum" not in level:
              level.append("szpowerspectrum")
        if "szpowerspectrum" in level:  #BB: added for class_sz
            if "lensing" not in level:
              level.append("class_sz_cosmo")
        if "class_sz_cosmo" in level:  #BB: added for class_sz
            if "lensing" not in level:
              level.append("lensing")
        #if "distortions" in level:
        #    if "lensing" not in level:
        #        level.append("lensing")
        if "lensing" in level:
            if "spectra" not in level:
                level.append("spectra")
        if "spectra" in level:
            if "transfer" not in level:
                level.append("transfer")
        if "transfer" in level:
            if "nonlinear" not in level:
                level.append("nonlinear")
        if "nonlinear" in level:
            if "primordial" not in level:
                level.append("primordial")
        if "primordial" in level:
            if "perturb" not in level:
                level.append("perturb")
        if "perturb" in level:
            if "thermodynamics" not in level:
                level.append("thermodynamics")
        if "thermodynamics" in level:
            if "background" not in level:
                level.append("background")
        if len(level)!=0 :
            if "input" not in level:
                level.append("input")
        return level

    def _pars_check(self, key, value, contains=False, add=""):
        val = ""
        if key in self._pars:
            val = self._pars[key]
            if contains:
                if value in val:
                    return True
            else:
                if value==val:
                    return True
        if add:
            sep = " "
            if isinstance(add,str):
                sep = add

            if contains and val:
                    self.set({key:val+sep+value})
            else:
                self.set({key:value})
            return True
        return False

    def compute(self, level=["szcount"]):
        """
        #compute(level=["lensing"])
        compute(level=["szcount"])

        Main function, execute all the _init methods for all desired modules.
        This is called in MontePython, and this ensures that the Class instance
        of this class contains all the relevant quantities. Then, one can deduce
        Pk, Cl, etc...

        Parameters
        ----------
        level : list
                list of the last module desired. The internal function
                _check_task_dependency will then add to this list all the
                necessary modules to compute in order to initialize this last
                one. The default last module is "lensing".

        .. warning::

            level default value should be left as an array (it was creating
            problem when casting as a set later on, in _check_task_dependency)

        """
        cdef ErrorMsg errmsg

        # Append to the list level all the modules necessary to compute.
        level = self._check_task_dependency(level)

        # Check if this function ran before (self.computed should be true), and
        # if no other modules were requested, i.e. if self.ncp contains (or is
        # equivalent to) level. If it is the case, simply stop the execution of
        # the function.
        if self.computed and self.ncp.issuperset(level):
            return

        # Check if already allocated to prevent memory leaks
        if self.allocated:
            self.struct_cleanup()

        # Otherwise, proceed with the normal computation.
        self.computed = False

        # Equivalent of writing a parameter file
        self._fillparfile()

        # self.ncp will contain the list of computed modules (under the form of
        # a set, instead of a python list)
        self.ncp=set()
        # Up until the empty set, all modules are allocated
        # (And then we successively keep track of the ones we allocate additionally)
        self.allocated = True

        # BB playground for emulators
        # cszfast = classy_szfast()
        # print(self._pars)
        # params_settings = self._pars


        # print('calculating cmb')
        # start = time.time()
        # cszfast.calculate_cmb(**params_settings)
        # end = time.time()
        # print('cmb calculation took:',end-start)


        # print('calculating pkl')
        # start = time.time()
        # cszfast.calculate_pkl(**params_settings)
        # end = time.time()
        # print('pk calculation took:',end-start)

        # print('calculating pknl')
        # start = time.time()
        # cszfast.calculate_pknl(**params_settings)
        # end = time.time()
        # print('pknl calculation took:',end-start)

        # print('tabulate sigma')
        # start = time.time()
        # cszfast.calculate_sigma(**params_settings)
        # end = time.time()
        # print('end tabulate sigma:',end-start)
        # print(cszfast.cszfast_pk_grid_dsigma2_flat[:10],np.shape(cszfast.cszfast_pk_grid_dsigma2_flat))
        # print(cszfast.cszfast_pk_grid_sigma2_flat[:10],np.shape(cszfast.cszfast_pk_grid_sigma2_flat))


        # --------------------------------------------------------------------
        # Check the presence for all CLASS modules in the list 'level'. If a
        # module is found in level, executure its "_init" method.
        # --------------------------------------------------------------------
        # The input module should raise a CosmoSevereError, because
        # non-understood parameters asked to the wrapper is a problematic
        # situation.
        if "input" in level:
            if input_init(&self.fc, &self.pr, &self.ba, &self.th,
                          &self.pt, &self.tr, &self.pm, &self.sp,
                          &self.nl, &self.le, &self.tsz, &self.csz, &self.op, errmsg) == _FAILURE_:
                raise CosmoSevereError(errmsg)
            self.ncp.add("input")
            # This part is done to list all the unread parameters, for debugging
            problem_flag = False
            problematic_parameters = []
            for i in range(self.fc.size):
                if self.fc.read[i] == _FALSE_:
                    problem_flag = True
                    problematic_parameters.append(self.fc.name[i].decode())
            if problem_flag:
                raise CosmoSevereError(
                    "Class did not read input parameter(s): %s\n" % ', '.join(
                    problematic_parameters))

        # The following list of computation is straightforward. If the "_init"
        # methods fail, call `struct_cleanup` and raise a CosmoComputationError
        # with the error message from the faulty module of CLASS.
        if "background" in level:
            if background_init(&(self.pr), &(self.ba)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.ba.error_message)
            self.ncp.add("background")

        if "thermodynamics" in level:
            if thermodynamics_init(&(self.pr), &(self.ba),
                                   &(self.th)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.th.error_message)
            self.ncp.add("thermodynamics")

        if "perturb" in level:
            if perturb_init(&(self.pr), &(self.ba),
                            &(self.th), &(self.pt)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pt.error_message)
            self.ncp.add("perturb")

        if "primordial" in level:
            if primordial_init(&(self.pr), &(self.pt),
                               &(self.pm)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.pm.error_message)
            self.ncp.add("primordial")

        if "nonlinear" in level:
            if nonlinear_init(&self.pr, &self.ba, &self.th,
                              &self.pt, &self.pm, &self.nl) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.nl.error_message)
            self.ncp.add("nonlinear")

        if "transfer" in level:
            if transfer_init(&(self.pr), &(self.ba), &(self.th),
                             &(self.pt), &(self.nl), &(self.tr)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.tr.error_message)
            self.ncp.add("transfer")

        if "spectra" in level:
            if spectra_init(&(self.pr), &(self.ba), &(self.pt),
                            &(self.pm), &(self.nl), &(self.tr),
                            &(self.sp)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.sp.error_message)
            self.ncp.add("spectra")

        if "lensing" in level:
            if lensing_init(&(self.pr), &(self.pt), &(self.sp),
                            &(self.nl), &(self.le)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.le.error_message)
            self.ncp.add("lensing")

        if "class_sz_cosmo" in level:
            if class_sz_cosmo_init(&(self.ba), &(self.th), &(self.pt), &(self.nl), &(self.pm),
            &(self.sp),&(self.le),&(self.tsz),&(self.pr)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.tsz.error_message)
            self.ncp.add("class_sz_cosmo")

        if "szpowerspectrum" in level:
            if szpowerspectrum_init(&(self.ba), &(self.th), &(self.pt), &(self.nl), &(self.pm),
            &(self.sp),&(self.le),&(self.tsz),&(self.pr)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.tsz.error_message)
            self.ncp.add("szpowerspectrum")


        if "szcount" in level:
            if szcount_init(&(self.ba), &(self.nl), &(self.pm),
            &(self.tsz),&(self.csz)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.tsz.error_message)
            self.ncp.add("szcount")


        self.computed = True

        # At this point, the cosmological instance contains everything needed. The
        # following functions are only to output the desired numbers
        return

    def compute_class_szfast(self):
        # print("print parameters:")
        # print(self._pars)
        self.compute(level=["thermodynamics"])
        # print(self._pars)
        params_settings = self._pars


        # BB playground for emulators
        # cszfast = classy_szfast()
        cszfast = classy_szfast(**params_settings)
        # print(self._pars)
        # params_settings = self._pars


        # print('calculating cmb')
        start = time.time()
        cszfast.calculate_cmb(**params_settings)
        end = time.time()
        # print('cmb calculation took:',end-start)


        # print('calculating pkl')
        start = time.time()
        cszfast.calculate_pkl(**params_settings)
        end = time.time()
        # print('pk calculation took:',end-start)

        # print('calculating pknl')
        start = time.time()
        cszfast.calculate_pknl(**params_settings)
        end = time.time()
        # print('pknl calculation took:',end-start)

        # print('calculating sigma8')
        start = time.time()
        cszfast.calculate_sigma8_and_der(**params_settings)
        end = time.time()
        # print('der calculation took:',end-start)
        # print('sigma8:',cszfast.sigma8)
        self.sigma8_fast = cszfast.sigma8

        cszfast.calculate_sigma8_at_z(**params_settings)


        self.tsz.use_class_sz_fast_mode = 1

        if class_sz_cosmo_init(&(self.ba), &(self.th), &(self.pt), &(self.nl), &(self.pm),
        &(self.sp),&(self.le),&(self.tsz),&(self.pr)) == _FAILURE_:
            self.struct_cleanup()
            raise CosmoComputationError(self.tsz.error_message)

        # print('tabulate sigma')
        start = time.time()
        cszfast.calculate_sigma(**params_settings)
        # print('lnsigma2',cszfast.cszfast_pk_grid_lnsigma2_flat)
        # print('dsigma2',cszfast.cszfast_pk_grid_dsigma2_flat)
        # print('ln1pz',cszfast.cszfast_pk_grid_ln1pz)
        # print('lnr',cszfast.cszfast_pk_grid_lnr)
        if self.tsz.need_sigma == 1:
          index_z_r = 0
          for index_z in range(self.tsz.n_arraySZ):
            for index_r in range(self.tsz.ndimSZ):
                  self.tsz.array_radius[index_r] = cszfast.cszfast_pk_grid_lnr[index_r]
                  self.tsz.array_redshift[index_z] = cszfast.cszfast_pk_grid_ln1pz[index_z]
                  self.tsz.array_sigma_at_z_and_R[index_z_r] = cszfast.cszfast_pk_grid_lnsigma2_flat[index_z_r]
                  self.tsz.array_dsigma2dR_at_z_and_R[index_z_r] = cszfast.cszfast_pk_grid_dsigma2_flat[index_z_r]
                  self.tsz.array_pkl_at_z_and_k[index_z_r] = cszfast.cszfast_pk_grid_pk_flat[index_z_r]
                  self.tsz.array_pknl_at_z_and_k[index_z_r] = cszfast.cszfast_pk_grid_pknl_flat[index_z_r]
                  self.tsz.array_lnk[index_r] = cszfast.cszfast_pk_grid_lnk[index_r]
                  index_z_r += 1
        end = time.time()
        # print('end tabulate sigma:',end-start)

        self.class_szfast = cszfast

        if szpowerspectrum_init(&(self.ba), &(self.th), &(self.pt), &(self.nl), &(self.pm),
        &(self.sp),&(self.le),&(self.tsz),&(self.pr)) == _FAILURE_:
            self.struct_cleanup()
            raise CosmoComputationError(self.tsz.error_message)
        if szcount_init(&(self.ba), &(self.nl), &(self.pm),
        &(self.tsz),&(self.csz)) == _FAILURE_:
            self.struct_cleanup()
            raise CosmoComputationError(self.tsz.error_message)

        self.computed = True
        return




    def compute_class_sz(self,pdict_to_update):
        self._fillparfile()
        for k,v in pdict_to_update.items():
          if k == 'M_min':
            self.tsz.M1SZ = pdict_to_update[k]
          if k == 'M_max':
            self.tsz.M2SZ = pdict_to_update[k]
          if k == 'z_min':
            self.tsz.z1SZ = pdict_to_update[k]
          if k == 'z_max':
            self.tsz.z2SZ = pdict_to_update[k]
          if k == 'fNL':
            self.tsz.fNL = pdict_to_update['fNL']
          if k == 'A_rho0':
            self.tsz.A_rho0 = pdict_to_update['A_rho0']
          if k == 'A_alpha':
            self.tsz.A_alpha = pdict_to_update['A_alpha']
          if k == 'A_beta':
              self.tsz.A_beta = pdict_to_update['A_beta']
          if k == 'alpha_m_rho0':
            self.tsz.alpha_m_rho0 = pdict_to_update['alpha_m_rho0']
          if k == 'alpha_m_alpha':
            self.tsz.alpha_m_alpha = pdict_to_update['alpha_m_alpha']
          if k == 'alpha_m_beta':
              self.tsz.alpha_m_beta = pdict_to_update['alpha_m_beta']
          if k == 'alpha_z_rho0':
            self.tsz.alpha_z_rho0 = pdict_to_update['alpha_z_rho0']
          if k == 'alpha_z_alpha':
            self.tsz.alpha_z_alpha = pdict_to_update['alpha_z_alpha']
          if k == 'alpha_z_beta':
              self.tsz.alpha_z_beta = pdict_to_update['alpha_z_beta']
          if k == 'c_B16':
              self.tsz.c_B16 = pdict_to_update['cp_B16']
          if k == 'mcut':
              self.tsz.mcut = pdict_to_update['mcut']
          if k == 'alphap_m_rho0':
              self.tsz.alphap_m_rho0 = pdict_to_update['alphap_m_rho0']
          if k == 'alphap_m_alpha':
             self.tsz.alphap_m_alpha = pdict_to_update['alphap_m_alpha']
          if k == 'alphap_m_beta':
              self.tsz.alphap_m_beta = pdict_to_update['alphap_m_beta']
          if k == 'alpha_c_rho0':
            self.tsz.alpha_c_rho0 = pdict_to_update['alpha_c_rho0']
          if k == 'alpha_c_alpha':
            self.tsz.alpha_c_alpha = pdict_to_update['alpha_c_alpha']
          if k == 'alpha_c_beta':
              self.tsz.alpha_c_beta = pdict_to_update['alpha_c_beta']
          if k == 'gamma_B16':
              self.tsz.gamma_B16 = pdict_to_update['gamma_B16']
          if k == 'xc_B16':
              self.tsz.xc_B16 = pdict_to_update['xc_B16']
          if k == 'P0_B12':
              self.tsz.P0_B12 = pdict_to_update['P0_B12']
          if k == 'beta_B12':
              self.tsz.beta_B12 = pdict_to_update['beta_B12']
          if k == 'alpha_B12':
              self.tsz.alpha_B12 = pdict_to_update['alpha_B12']
          if k == 'gamma_B12':
              self.tsz.gamma_B12 = pdict_to_update['gamma_B12']
          if k == 'xc_B12':
              self.tsz.xc_B12 = pdict_to_update['xc_B12']
          if k == 'alpha_m_P0_B12':
              self.tsz.alpha_m_P0_B12 = pdict_to_update['alpha_m_P0_B12']
          if k == 'alpha_m_xc_B12':
              self.tsz.alpha_m_xc_B12 = pdict_to_update['alpha_m_xc_B12']
          if k == 'alpha_m_beta_B12':
              self.tsz.alpha_m_beta_B12 = pdict_to_update['alpha_m_beta_B12']
          if k == 'alpha_z_P0_B12':
              self.tsz.alpha_z_P0_B12 = pdict_to_update['alpha_z_P0_B12']
          if k == 'alpha_z_xc_B12':
              self.tsz.alpha_z_xc_B12 = pdict_to_update['alpha_z_xc_B12']
          if k == 'alpha_z_beta_B12':
              self.tsz.alpha_z_beta_B12 = pdict_to_update['alpha_z_beta_B12']
          if k == 'c_B12':
              self.tsz.c_B12 = pdict_to_update['cp_B12']
          if k == 'mcut_B12':
              self.tsz.mcut_B12 = pdict_to_update['mcut_B12']
          if k == 'alphap_m_P0_B12':
              self.tsz.alphap_m_P0_B12 = pdict_to_update['alphap_m_P0_B12']
          if k == 'alphap_m_xc_B12':
              self.tsz.alphap_m_xc_B12 = pdict_to_update['alphap_m_xc_B12']
          if k == 'alphap_m_beta_B12':
              self.tsz.alphap_m_beta_B12 = pdict_to_update['alphap_m_beta_B12']
          if k == 'alpha_c_P0_B12':
              self.tsz.alpha_c_P0_B12 = pdict_to_update['alpha_c_P0_B12']
          if k == 'alpha_c_xc_B12':
              self.tsz.alpha_c_xc_B12 = pdict_to_update['alpha_c_xc_B12']
          if k == 'alpha_c_beta_B12':
              self.tsz.alpha_c_beta_B12 = pdict_to_update['alpha_c_beta_B12']
          if k == 'x_outSZ':
              self.tsz.x_outSZ = pdict_to_update['x_outSZ']
        # print('array_redshift:',
        #       self.tsz.array_redshift[0],
        #       self.tsz.array_redshift[1],
        #       self.tsz.array_redshift[self.tsz.n_arraySZ-1])
        # print('array_radius:',
        #       self.tsz.array_radius[0],
        #       self.tsz.array_radius[1],
        #       self.tsz.array_radius[self.tsz.ndimSZ-1])
        # print('array_dsigma2dR_at_z_and_R:',
        #       self.tsz.array_dsigma2dR_at_z_and_R[0],
        #       self.tsz.array_dsigma2dR_at_z_and_R[1],
        #       self.tsz.array_dsigma2dR_at_z_and_R[self.tsz.n_arraySZ*self.tsz.ndimSZ-1])
        # print('array_sigma_at_z_and_R:',
        #       self.tsz.array_sigma_at_z_and_R[0],
        #       self.tsz.array_sigma_at_z_and_R[1],
        #       self.tsz.array_sigma_at_z_and_R[self.tsz.n_arraySZ*self.tsz.ndimSZ-1])

        # # BB playground for emulators

        # print(self._pars)
        # params_settings = self._pars
        # cszfast = classy_szfast(**params_settings)


        # print('calculating cmb')
        # start = time.time()
        # cszfast.calculate_cmb(**params_settings)
        # end = time.time()
        # print('cmb calculation took:',end-start)


        # print('calculating pkl')
        # start = time.time()
        # cszfast.calculate_pkl(**params_settings)
        # end = time.time()
        # print('pk calculation took:',end-start)

        # print('calculating pknl')
        # start = time.time()
        # cszfast.calculate_pknl(**params_settings)
        # end = time.time()
        # print('pknl calculation took:',end-start)


        # print('tabulate sigma')
        # start = time.time()
        # cszfast.calculate_sigma(**params_settings)
        # print('lnsigma2',cszfast.cszfast_pk_grid_lnsigma2_flat)
        # print('dsigma2',cszfast.cszfast_pk_grid_dsigma2_flat)
        # print('ln1pz',cszfast.cszfast_pk_grid_ln1pz)
        # print('lnr',cszfast.cszfast_pk_grid_lnr)
        # index_z_r = 0
        # for index_z in range(self.tsz.n_arraySZ):
        #   for index_r in range(self.tsz.ndimSZ):
        #         self.tsz.array_radius[index_r] = cszfast.cszfast_pk_grid_lnr[index_r]
        #         self.tsz.array_redshift[index_z] = cszfast.cszfast_pk_grid_ln1pz[index_z]
        #         self.tsz.array_sigma_at_z_and_R[index_z_r] = cszfast.cszfast_pk_grid_lnsigma2_flat[index_z_r]
        #         self.tsz.array_dsigma2dR_at_z_and_R[index_z_r] = cszfast.cszfast_pk_grid_dsigma2_flat[index_z_r]
        #         self.tsz.array_pkl_at_z_and_k[index_z_r] = cszfast.cszfast_pk_grid_pk_flat[index_z_r]
        #         self.tsz.array_pknl_at_z_and_k[index_z_r] = cszfast.cszfast_pk_grid_pknl_flat[index_z_r]
        #         self.tsz.array_lnk[index_r] = cszfast.cszfast_pk_grid_lnk[index_r]
        #         index_z_r += 1
        # end = time.time()
        # print('end tabulate sigma:',end-start)
        if szpowerspectrum_init(&(self.ba), &(self.th), &(self.pt), &(self.nl), &(self.pm),
        &(self.sp),&(self.le),&(self.tsz),&(self.pr)) == _FAILURE_:
            self.struct_cleanup()
            raise CosmoComputationError(self.tsz.error_message)
        self.computed = True
        return




    def raw_cl(self, lmax=-1, nofail=False):
        """
        raw_cl(lmax=-1, nofail=False)

        Return a dictionary of the primary C_l

        Parameters
        ----------
        lmax : int, optional
                Define the maximum l for which the C_l will be returned
                (inclusively). This number will be checked against the maximum l
                at which they were actually computed by CLASS, and an error will
                be raised if the desired lmax is bigger than what CLASS can
                give.
        nofail: bool, optional
                Check and enforce the computation of the spectra module
                beforehand, with the desired lmax.

        Returns
        -------
        cl : dict
                Dictionary that contains the power spectrum for 'tt', 'te', etc... The
                index associated with each is defined wrt. Class convention, and are non
                important from the python point of view. It also returns now the
                ell array.
        """
        cdef int lmaxR
        cdef double *rcl = <double*> calloc(self.sp.ct_size,sizeof(double))

        # Quantities for tensor modes
        cdef double **cl_md = <double**> calloc(self.sp.md_size, sizeof(double*))
        for index_md in range(self.sp.md_size):
            cl_md[index_md] = <double*> calloc(self.sp.ct_size, sizeof(double))

        # Quantities for isocurvature modes
        cdef double **cl_md_ic = <double**> calloc(self.sp.md_size, sizeof(double*))
        for index_md in range(self.sp.md_size):
            cl_md_ic[index_md] = <double*> calloc(self.sp.ct_size*self.sp.ic_ic_size[index_md], sizeof(double))

        # Define a list of integers, refering to the flags and indices of each
        # possible output Cl. It allows for a clear and concise way of looping
        # over them, checking if they are defined or not.
        has_flags = [
            (self.sp.has_tt, self.sp.index_ct_tt, 'tt'),
            (self.sp.has_ee, self.sp.index_ct_ee, 'ee'),
            (self.sp.has_te, self.sp.index_ct_te, 'te'),
            (self.sp.has_bb, self.sp.index_ct_bb, 'bb'),
            (self.sp.has_pp, self.sp.index_ct_pp, 'pp'),
            (self.sp.has_tp, self.sp.index_ct_tp, 'tp'),]
        spectra = []

        for flag, index, name in has_flags:
            if flag:
                spectra.append(name)

        if not spectra:
            raise CosmoSevereError("No Cl computed")
        lmaxR = self.sp.l_max_tot
        if lmax == -1:
            lmax = lmaxR
        if lmax > lmaxR:
            if nofail:
                self._pars_check("l_max_scalars",lmax)
                self.compute(["lensing"])
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        # Initialise all the needed Cls arrays
        cl = {}
        for elem in spectra:
            cl[elem] = np.zeros(lmax+1, dtype=np.double)

        # Recover for each ell the information from CLASS
        for ell from 2<=ell<lmax+1:
            if spectra_cl_at_l(&self.sp, ell, rcl, cl_md, cl_md_ic) == _FAILURE_:
                raise CosmoSevereError(self.sp.error_message)
            for flag, index, name in has_flags:
                if name in spectra:
                    cl[name][ell] = rcl[index]
        cl['ell'] = np.arange(lmax+1)

        free(rcl)
        for index_md in range(self.sp.md_size):
            free(cl_md[index_md])
            free(cl_md_ic[index_md])
        free(cl_md)
        free(cl_md_ic)

        return cl

    def lensed_cl(self, lmax=-1,nofail=False):
        """
        lensed_cl(lmax=-1, nofail=False)

        Return a dictionary of the lensed C_l, computed by CLASS, without the
        density C_ls. They must be asked separately with the function aptly
        named density_cl

        Parameters
        ----------
        lmax : int, optional
                Define the maximum l for which the C_l will be returned (inclusively)
        nofail: bool, optional
                Check and enforce the computation of the lensing module beforehand

        Returns
        -------
        cl : dict
                Dictionary that contains the power spectrum for 'tt', 'te', etc... The
                index associated with each is defined wrt. Class convention, and are non
                important from the python point of view.
        """
        cdef int lmaxR
        cdef double *lcl = <double*> calloc(self.le.lt_size,sizeof(double))

        if self.tsz.use_class_sz_fast_mode == 1:
          # bypass all together and replace by cosmopower call.
          cls = {}
          cls['ell'] = np.arange(20000)
          cls['tt'] = np.zeros(20000)
          cls['te'] = np.zeros(20000)
          cls['ee'] = np.zeros(20000)
          cls['pp'] = np.zeros(20000)

          nl = len(self.class_szfast.cp_predicted_tt_spectrum)
          lcp = np.asarray(cls['ell'][2:nl+2])
          # print('nl:',nl)
          cls['tt'][2:nl+2] = self.class_szfast.cp_predicted_tt_spectrum
          cls['tt'][2:nl+2] *= 1./(lcp*(lcp+1.)/2./np.pi)
          cls['te'][2:nl+2] = self.class_szfast.cp_predicted_te_spectrum
          cls['te'][2:nl+2] *= 1./(lcp*(lcp+1.)/2./np.pi)
          cls['ee'][2:nl+2] = self.class_szfast.cp_predicted_ee_spectrum
          cls['ee'][2:nl+2] *= 1./(lcp*(lcp+1.)/2./np.pi)
          # cls['pp'][2:nl+2] = self.class_szfast.cp_predicted_pp_spectrum/4. ## this is clkk... works for so lensinglite lkl
          cls['pp'][2:nl+2] = self.class_szfast.cp_predicted_pp_spectrum/(lcp*(lcp+1.))**2.
          # # here for the planck lensing lkl, using lfactor option gives:
          # cls['pp'][2:nl+2] = self.cp_predicted_pp_spectrum/(lcp*(lcp+1.))**2.
          # cls['pp'][2:nl+2] *= (lcp*(lcp+1.))**2./2./np.pi
          free(lcl)
          return cls
        else:

          # Define a list of integers, refering to the flags and indices of each
          # possible output Cl. It allows for a clear and concise way of looping
          # over them, checking if they are defined or not.
          has_flags = [
              (self.le.has_tt, self.le.index_lt_tt, 'tt'),
              (self.le.has_ee, self.le.index_lt_ee, 'ee'),
              (self.le.has_te, self.le.index_lt_te, 'te'),
              (self.le.has_bb, self.le.index_lt_bb, 'bb'),
              (self.le.has_pp, self.le.index_lt_pp, 'pp'),
              (self.le.has_tp, self.le.index_lt_tp, 'tp'),]
          spectra = []

          for flag, index, name in has_flags:
              if flag:
                  spectra.append(name)

          if not spectra:
              raise CosmoSevereError("No lensed Cl computed")
          lmaxR = self.le.l_lensed_max

          if lmax == -1:
              lmax = lmaxR
          if lmax > lmaxR:
              if nofail:
                  self._pars_check("l_max_scalars",lmax)
                  self.compute(["lensing"])
              else:
                  raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

          cl = {}
          # Simple Cls, for temperature and polarisation, are not so big in size
          for elem in spectra:
              cl[elem] = np.zeros(lmax+1, dtype=np.double)
          for ell from 2<=ell<lmax+1:
              if lensing_cl_at_l(&self.le,ell,lcl) == _FAILURE_:
                  raise CosmoSevereError(self.le.error_message)
              for flag, index, name in has_flags:
                  if name in spectra:
                      cl[name][ell] = lcl[index]
          cl['ell'] = np.arange(lmax+1)

          free(lcl)
          return cl

    def density_cl(self, lmax=-1, nofail=False):
        """
        density_cl(lmax=-1, nofail=False)

        Return a dictionary of the primary C_l for the matter

        Parameters
        ----------
        lmax : int, optional
            Define the maximum l for which the C_l will be returned (inclusively)
        nofail: bool, optional
            Check and enforce the computation of the lensing module beforehand

        Returns
        -------
        cl : numpy array of numpy.ndarrays
            Array that contains the list (in this order) of self correlation of
            1st bin, then successive correlations (set by non_diagonal) to the
            following bins, then self correlation of 2nd bin, etc. The array
            starts at index_ct_dd.
        """
        cdef int lmaxR
        cdef double *dcl = <double*> calloc(self.sp.ct_size,sizeof(double))

        # Quantities for tensor modes
        cdef double **cl_md = <double**> calloc(self.sp.md_size, sizeof(double*))
        for index_md in range(self.sp.md_size):
            cl_md[index_md] = <double*> calloc(self.sp.ct_size, sizeof(double))

        # Quantities for isocurvature modes
        cdef double **cl_md_ic = <double**> calloc(self.sp.md_size, sizeof(double*))
        for index_md in range(self.sp.md_size):
            cl_md_ic[index_md] = <double*> calloc(self.sp.ct_size*self.sp.ic_ic_size[index_md], sizeof(double))

        lmaxR = self.pt.l_lss_max
        has_flags = [
            (self.sp.has_dd, self.sp.index_ct_dd, 'dd'),
            (self.sp.has_td, self.sp.index_ct_td, 'td'),
            (self.sp.has_ll, self.sp.index_ct_ll, 'll'),
            (self.sp.has_dl, self.sp.index_ct_dl, 'dl'),
            (self.sp.has_tl, self.sp.index_ct_tl, 'tl')]
        spectra = []

        for flag, index, name in has_flags:
            if flag:
                spectra.append(name)
                l_max_flag = self.sp.l_max_ct[self.sp.index_md_scalars][index]
                if l_max_flag < lmax and lmax > 0:
                    raise CosmoSevereError(
                        "the %s spectrum was computed until l=%i " % (
                            name.upper(), l_max_flag) +
                        "but you asked a l=%i" % lmax)

        if not spectra:
            raise CosmoSevereError("No density Cl computed")
        if lmax == -1:
            lmax = lmaxR
        if lmax > lmaxR:
            if nofail:
                self._pars_check("l_max_lss",lmax)
                self._pars_check("output",'nCl')
                self.compute()
            else:
                raise CosmoSevereError("Can only compute up to lmax=%d"%lmaxR)

        cl = {}

        # For density Cls, the size is bigger (different redshfit bins)
        # computes the size, given the number of correlations needed to be computed
        size = int((self.sp.d_size*(self.sp.d_size+1)-(self.sp.d_size-self.sp.non_diag)*
                (self.sp.d_size-1-self.sp.non_diag))/2);
        for elem in ['dd', 'll', 'dl']:
            if elem in spectra:
                cl[elem] = {}
                for index in range(size):
                    cl[elem][index] = np.zeros(
                        lmax+1, dtype=np.double)
        for elem in ['td', 'tl']:
            if elem in spectra:
                cl[elem] = np.zeros(lmax+1, dtype=np.double)

        for ell from 2<=ell<lmax+1:
            if spectra_cl_at_l(&self.sp, ell, dcl, cl_md, cl_md_ic) == _FAILURE_:
                raise CosmoSevereError(self.sp.error_message)
            if 'dd' in spectra:
                for index in range(size):
                    cl['dd'][index][ell] = dcl[self.sp.index_ct_dd+index]
            if 'll' in spectra:
                for index in range(size):
                    cl['ll'][index][ell] = dcl[self.sp.index_ct_ll+index]
            if 'dl' in spectra:
                for index in range(size):
                    cl['dl'][index][ell] = dcl[self.sp.index_ct_dl+index]
            if 'td' in spectra:
                cl['td'][ell] = dcl[self.sp.index_ct_td]
            if 'tl' in spectra:
                cl['tl'][ell] = dcl[self.sp.index_ct_tl]
        cl['ell'] = np.arange(lmax+1)

        free(dcl)
        for index_md in range(self.sp.md_size):
            free(cl_md[index_md])
            free(cl_md_ic[index_md])
        free(cl_md)
        free(cl_md_ic)

        return cl

    def z_of_r (self,z_array):
        cdef double tau=0.0
        cdef int last_index=0 #junk
        cdef double * pvecback
        r = np.zeros(len(z_array),'float64')
        dzdr = np.zeros(len(z_array),'float64')

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        i = 0
        for redshift in z_array:
            if background_tau_of_z(&self.ba,redshift,&tau)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
                raise CosmoSevereError(self.ba.error_message)

            # store r
            r[i] = pvecback[self.ba.index_bg_conf_distance]
            # store dz/dr = H
            dzdr[i] = pvecback[self.ba.index_bg_H]

            i += 1

        free(pvecback)
        return r[:],dzdr[:]

    def luminosity_distance(self, z):
        """
        luminosity_distance(z)
        """
        cdef double tau=0.0
        cdef int last_index = 0  # junk
        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba, z, &tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba, tau, self.ba.long_info,
                self.ba.inter_normal, &last_index, pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)
        lum_distance = pvecback[self.ba.index_bg_lum_distance]
        free(pvecback)
        return lum_distance

    # Gives the total matter pk for a given (k,z)
    def pk(self,double k,double z):
        """
        Gives the total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk
        if self.tsz.use_class_sz_fast_mode == 1:
          pk = self.class_szfast.get_pknl_at_k_and_z(k,z)
        else:
          if (self.pt.has_pk_matter == _FALSE_):
              raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")

          if (self.nl.method == nl_none):
              if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_linear,k,z,self.nl.index_pk_m,&pk,NULL)==_FAILURE_:
                  raise CosmoSevereError(self.nl.error_message)
          else:
              if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_nonlinear,k,z,self.nl.index_pk_m,&pk,NULL)==_FAILURE_:
                  raise CosmoSevereError(self.nl.error_message)

        return pk

    # Gives the cdm+b pk for a given (k,z)
    def pk_cb(self,double k,double z):
        """
        Gives the cdm+b pk (in Mpc**3) for a given k (in 1/Mpc) and z (will be non linear if requested to Class, linear otherwise)

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk_cb

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")
        if (self.nl.has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed (probably because there are no massive neutrinos) so you cannot ask for it")

        if (self.nl.method == nl_none):
            if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_linear,k,z,self.nl.index_pk_cb,&pk_cb,NULL)==_FAILURE_:
                raise CosmoSevereError(self.nl.error_message)
        else:
            if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_nonlinear,k,z,self.nl.index_pk_cb,&pk_cb,NULL)==_FAILURE_:
                raise CosmoSevereError(self.nl.error_message)

        return pk_cb

    # Gives the total matter pk for a given (k,z)
    def pk_lin(self,double k,double z):
        """
        Gives the linear total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk_lin

        if self.tsz.use_class_sz_fast_mode == 1:
          pk_lin = self.class_szfast.get_pkl_at_k_and_z(k,z)
        else:
          if (self.pt.has_pk_matter == _FALSE_):
              raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")

          if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_linear,k,z,self.nl.index_pk_m,&pk_lin,NULL)==_FAILURE_:
              raise CosmoSevereError(self.nl.error_message)

        return pk_lin

    # Gives the total matter pk for a given (k,z)
    def pk_nonlin(self,double k,double z):
        """
        Gives the linear total matter pk (in Mpc**3) for a given k (in 1/Mpc) and z

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk_lin

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")

        if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_nonlinear,k,z,self.nl.index_pk_m,&pk_lin,NULL)==_FAILURE_:
            raise CosmoSevereError(self.nl.error_message)

        return pk_lin

    # Gives the cdm+b pk for a given (k,z)
    def pk_cb_lin(self,double k,double z):
        """
        Gives the linear cdm+b pk (in Mpc**3) for a given k (in 1/Mpc) and z

        .. note::

            there is an additional check that output contains `mPk`,
            because otherwise a segfault will occur

        """
        cdef double pk_cb_lin

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. You must add mPk to the list of outputs.")

        if (self.nl.has_pk_cb == _FALSE_):
            raise CosmoSevereError("P_cb not computed by CLASS (probably because there are no massive neutrinos)")

        if nonlinear_pk_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_linear,k,z,self.nl.index_pk_cb,&pk_cb_lin,NULL)==_FAILURE_:
            raise CosmoSevereError(self.nl.error_message)

        return pk_cb_lin

    def get_pk(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=3] pk = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk[index_k,index_z,index_mu] = self.pk(k[index_k,index_z,index_mu],z[index_z])
        return pk

    def get_pk_cb(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=3] pk_cb = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk_cb[index_k,index_z,index_mu] = self.pk_cb(k[index_k,index_z,index_mu],z[index_z])
        return pk_cb

    def get_pk_lin(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the linear power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=3] pk = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk[index_k,index_z,index_mu] = self.pk_lin(k[index_k,index_z,index_mu],z[index_z])
        return pk

    def get_pk_cb_lin(self, np.ndarray[DTYPE_t,ndim=3] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, int mu_size):
        """ Fast function to get the linear power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=3] pk_cb = np.zeros((k_size,z_size,mu_size),'float64')
        cdef int index_k, index_z, index_mu

        for index_k in xrange(k_size):
            for index_z in xrange(z_size):
                for index_mu in xrange(mu_size):
                    pk_cb[index_k,index_z,index_mu] = self.pk_cb_lin(k[index_k,index_z,index_mu],z[index_z])
        return pk_cb

    def get_pk_and_k_and_z(self, nonlinear=True, only_clustering_species = False):
        """
        Returns a grid of matter power spectrum values and the z and k
        at which it has been fully computed. Useful for creating interpolators.

        Parameters
        ----------
        nonlinear : bool
                Whether the returned power spectrum values are linear or non-linear (default)
        nonlinear : bool
                Whether the returned power spectrum is for galaxy clustering and excludes massive neutrinos, or always includes evrything (default)
        """
        cdef np.ndarray[DTYPE_t,ndim=2] pk_at_k_z = np.zeros((self.nl.k_size, self.nl.ln_tau_size),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] k = np.zeros((self.nl.k_size),'float64')
        cdef np.ndarray[DTYPE_t,ndim=1] z = np.zeros((self.nl.ln_tau_size),'float64')
        cdef int index_k, index_tau, index_pk
        cdef double z_max_nonlinear, z_max_requested

        # consistency checks

        if self.nl.has_pk_matter == False:
            raise CosmoSevereError("You ask classy to return an array of P(k,z) values, but the input parameters sent to CLASS did not require any P(k,z) calculations; add 'mPk' in 'output'")

        if nonlinear == True and self.nl.method == nl_none:
            raise CosmoSevereError("You ask classy to return an array of nonlinear P(k,z) values, but the input parameters sent to CLASS did not require any non-linear P(k,z) calculations; add e.g. 'halofit' or 'HMcode' in 'nonlinear'")

        # check wich type of P(k) to return (total or clustering only, i.e. without massive neutrino contribution)
        if (only_clustering_species == True):
            index_pk = self.nl.index_pk_cluster
        else:
            index_pk = self.nl.index_pk_total

        # get list of redshfits

        if self.nl.ln_tau_size == 1:
            raise CosmoSevereError("You ask classy to return an array of P(k,z) values, but the input parameters sent to CLASS did not require any P(k,z) calculations for z>0; pass either a list of z in 'z_pk' or one non-zero value in 'z_max_pk'")
        else:
            for index_tau in xrange(self.nl.ln_tau_size):
                if index_tau == self.nl.ln_tau_size-1:
                    z[index_tau] = 0.
                else:
                    z[index_tau] = self.z_of_tau(np.exp(self.nl.ln_tau[index_tau]))

        # check consitency of the list of redshifts

        if nonlinear == True:
            z_max_nonlinear = self.z_of_tau(self.nl.tau[self.nl.index_tau_min_nl])
            z_max_requested = z[0]
            if ((self.nl.tau_size - self.nl.ln_tau_size) < self.nl.index_tau_min_nl):
                raise CosmoSevereError("get_pk_and_k_and_z() is trying to return P(k,z) up to z_max=%e (to encompass your requested maximum value of z); but the input parameters sent to CLASS were such that the non-linear P(k,z) could only be consistently computed up to z=%e; increase the input parameter 'P_k_max_h/Mpc' or 'P_k_max_1/Mpc', or increase the precision parameters 'halofit_min_k_max' and/or 'hmcode_min_k_max', or decrease your requested z_max"%(z_max_requested,z_max_nonlinear))

        # get list of k

        for index_k in xrange(self.nl.k_size):
            k[index_k] = self.nl.k[index_k]

        # get P(k,z) array

        for index_tau in xrange(self.nl.ln_tau_size):
            for index_k in xrange(self.nl.k_size):
                if nonlinear == True:
                    pk_at_k_z[index_k, index_tau] = np.exp(self.nl.ln_pk_nl[index_pk][index_tau * self.nl.k_size + index_k])
                else:
                    pk_at_k_z[index_k, index_tau] = np.exp(self.nl.ln_pk_l[index_pk][index_tau * self.nl.k_size + index_k])

        return pk_at_k_z, k, z

    # Gives sigma(R,z) for a given (R,z)
    def sigma(self,double R,double z, h_units = False):
        """
        Gives sigma (total matter) for a given R and z
        (R is the radius in units of Mpc, so if R=8/h this will be the usual sigma8(z)

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        cdef double sigma
        if self.tsz.use_class_sz_fast_mode == 1:
          sigma = self.class_szfast.get_sigma_at_r_and_z(R,z)
          # if h_units == True:
          #   sigma = self.class_szfast.get_sigma_at_r_and_z(R*self.h(),z)
          # else:
          #   sigma = self.class_szfast.get_sigma_at_r_and_z(R,z)
        else:
          if (self.pt.has_pk_matter == _FALSE_):
              raise CosmoSevereError("No power spectrum computed. In order to get sigma(R,z) you must add mPk to the list of outputs.")

          if (self.pt.k_max_for_pk < self.ba.h):
              raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

          if nonlinear_sigmas_at_z(&self.pr,&self.ba,&self.nl,R,z,self.nl.index_pk_m,out_sigma,&sigma)==_FAILURE_:
              raise CosmoSevereError(self.nl.error_message)

        return sigma

    # Gives sigma_cb(R,z) for a given (R,z)
    def sigma_cb(self,double R,double z):
        """
        Gives sigma (cdm+b) for a given R and z
        (R is the radius in units of Mpc, so if R=8/h this will be the usual sigma8(z)

        .. note::

            there is an additional check to verify whether output contains `mPk`,
            and whether k_max > ...
            because otherwise a segfault will occur

        """
        cdef double sigma_cb

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get sigma(R,z) you must add mPk to the list of outputs.")

        if (self.nl.has_pk_cb == _FALSE_):
            raise CosmoSevereError("sigma_cb not computed by CLASS (probably because there are no massive neutrinos)")

        if (self.pt.k_max_for_pk < self.ba.h):
            raise CosmoSevereError("In order to get sigma(R,z) you must set 'P_k_max_h/Mpc' to 1 or bigger, in order to have k_max > 1 h/Mpc.")

        if nonlinear_sigmas_at_z(&self.pr,&self.ba,&self.nl,R,z,self.nl.index_pk_cb,out_sigma,&sigma_cb)==_FAILURE_:
            raise CosmoSevereError(self.nl.error_message)

        return sigma_cb

    # Gives effective logarithmic slope of P_L(k,z) (total matter) for a given (k,z)
    def pk_tilt(self,double k,double z):
        """
        Gives effective logarithmic slope of P_L(k,z) (total matter) for a given k and z
        (k is the wavenumber in units of 1/Mpc, z is the redshift, the output is dimensionless)

        .. note::

            there is an additional check to verify whether output contains `mPk` and whether k is in the right range

        """
        cdef double pk_tilt

        if (self.pt.has_pk_matter == _FALSE_):
            raise CosmoSevereError("No power spectrum computed. In order to get pk_tilt(k,z) you must add mPk to the list of outputs.")

        if (k < self.nl.k[1] or k > self.nl.k[self.nl.k_size-2]):
            raise CosmoSevereError("In order to get pk_tilt at k=%e 1/Mpc, you should compute P(k,z) in a wider range of k's"%k)

        if nonlinear_pk_tilt_at_k_and_z(&self.ba,&self.pm,&self.nl,pk_linear,k,z,self.nl.index_pk_total,&pk_tilt)==_FAILURE_:
            raise CosmoSevereError(self.nl.error_message)

        return pk_tilt

    #calculates the hmcode window_function of the Navarrow Frenk White Profile
    def nonlinear_hmcode_window_nfw(self,double k,double rv,double c):
        """
        Gives window_nfw for a given wavevector k, virial radius rv and concentration c

        """
        cdef double window_nfw


        if nonlinear_hmcode_window_nfw(&self.nl,k,rv,c,&window_nfw)==_FAILURE_:
                 raise CosmoSevereError(self.sp.error_message)

        return window_nfw


    #################################
    # gives an estimation of f(z)*sigma8(z) at the scale of 8 h/Mpc, computed as (d sigma8/d ln a)
    def effective_f_sigma8(self, z, z_step=0.1):
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

        # we need d sigma8/d ln a = - (d sigma8/dz)*(1+z)

        # if possible, use two-sided derivative with default value of z_step
        if z >= z_step:
            return (self.sigma(8/self.h(),z-z_step,h_units=True)-self.sigma(8/self.h(),z+z_step,h_units=True))/(2.*z_step)*(1+z)
        else:
            # if z is between z_step/10 and z_step, reduce z_step to z, and then stick to two-sided derivative
            if (z > z_step/10.):
                z_step = z
                return (self.sigma(8/self.h(),z-z_step,h_units=True)-self.sigma(8/self.h(),z+z_step,h_units=True))/(2.*z_step)*(1+z)
            # if z is between 0 and z_step/10, use single-sided derivative with z_step/10
            else:
                z_step /=10
                return (self.sigma(8/self.h(),z,h_units=True)-self.sigma(8/self.h(),z+z_step,h_units=True))/z_step*(1+z)



    def age(self):
        if self.tsz.use_class_sz_fast_mode == 1:
          return self.ba.age
        else:
          self.compute(["background"])
          return self.ba.age

    def h(self):
        return self.ba.h

    def n_s(self):
        return self.pm.n_s

    def tau_reio(self):
        return self.th.tau_reio

    def Omega_m(self):
        return self.ba.Omega0_m

    def Omega_r(self):
        return self.ba.Omega0_r

    def theta_s_100(self):
        return 100.*self.th.rs_rec/self.th.da_rec/(1.+self.th.z_rec)

    def theta_star_100(self):
        return 100.*self.th.rs_star/self.th.da_star/(1.+self.th.z_star)

    def Omega_Lambda(self):
        return self.ba.Omega0_lambda

    def Omega_g(self):
        return self.ba.Omega0_g

    def Omega_b(self):
        return self.ba.Omega0_b

    def omega_b(self):
        return self.ba.Omega0_b * self.ba.h * self.ba.h

    def Neff(self):
        return self.ba.Neff

    def k_eq(self):
        self.compute(["background"])
        return self.ba.a_eq*self.ba.H_eq

    def sigma8(self):
        if self.tsz.use_class_sz_fast_mode == 1:
          return self.sigma8_fast
        else:
          self.compute(["nonlinear"])
          return self.nl.sigma8[self.nl.index_pk_m]

    def get_sigma8_at_z(self,z):
        return self.class_szfast.get_sigma8_at_z(z)

    def get_effective_f_sigma8(self,z):
        return self.class_szfast.get_effective_f_sigma8(z)


    #def neff(self):
    #    self.compute(["spectra"])
    #    return self.sp.neff

    def sigma8_cb(self):
        self.compute(["nonlinear"])
        return self.nl.sigma8[self.nl.index_pk_cb]

    def rs_drag(self):
        if self.tsz.use_class_sz_fast_mode == 1:
          return self.th.rs_d
        else:
          self.compute(["thermodynamics"])
          return self.th.rs_d

    def z_reio(self):
        self.compute(["thermodynamics"])
        return self.th.z_reio

    def angular_distance(self, z):
        """
        angular_distance(z)

        Return the angular diameter distance (exactly, the quantity defined by Class
        as index_bg_ang_distance in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D_A = pvecback[self.ba.index_bg_ang_distance]

        free(pvecback)

        return D_A


    def scale_independent_growth_factor(self, z):
        """
        scale_independent_growth_factor(z)

        Return the scale invariant growth factor D(a) for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_D in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D = pvecback[self.ba.index_bg_D]

        free(pvecback)

        return D

    def scale_independent_growth_factor_f(self, z):
        """
        scale_independent_growth_factor_f(z)

        Return the scale invariant growth factor f(z)=d ln D / d ln a for CDM perturbations
        (exactly, the quantity defined by Class as index_bg_f in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        f = pvecback[self.ba.index_bg_f]

        free(pvecback)

        return f

    def z_of_tau(self, tau):
        """
        Redshift corresponding to a given conformal time.

        Parameters
        ----------
        tau : float
                Conformal time
        """
        cdef double z
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        z = 1./pvecback[self.ba.index_bg_a]-1.

        free(pvecback)

        return z

    def Hubble(self, z):
        """
        Hubble(z)

        Return the Hubble rate (exactly, the quantity defined by Class as index_bg_H
        in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        H = pvecback[self.ba.index_bg_H]

        free(pvecback)

        return H

    def Om_m(self, z):
        """
        Omega_m(z)

        Return the matter density fraction (exactly, the quantity defined by Class as index_bg_Omega_m
        in the background module)

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        Om_m = pvecback[self.ba.index_bg_Omega_m]

        free(pvecback)

        return Om_m


    def ionization_fraction(self, z):
        """
        ionization_fraction(z)

        Return the ionization fraction for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback
        cdef double * pvecthermo

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
        pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if thermodynamics_at_z(&self.ba,&self.th,z,self.th.inter_normal,&last_index,pvecback,pvecthermo) == _FAILURE_:
            raise CosmoSevereError(self.th.error_message)

        xe = pvecthermo[self.th.index_th_xe]

        free(pvecback)
        free(pvecthermo)

        return xe

    def baryon_temperature(self, z):
        """
        baryon_temperature(z)

        Give the baryon temperature for a given redshift z

        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback
        cdef double * pvecthermo

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))
        pvecthermo = <double*> calloc(self.th.th_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if thermodynamics_at_z(&self.ba,&self.th,z,self.th.inter_normal,&last_index,pvecback,pvecthermo) == _FAILURE_:
            raise CosmoSevereError(self.th.error_message)

        Tb = pvecthermo[self.th.index_th_Tb]

        free(pvecback)
        free(pvecthermo)

        return Tb

    def T_cmb(self):
        """
        Return the CMB temperature
        """
        return self.ba.T_cmb

    # redundent with a previous Omega_m() funciton,
    # but we leave it not to break compatibility
    def Omega0_m(self):
        """
        Return the sum of Omega0 for all non-relativistic components
        """
        return self.ba.Omega0_m


    def cl_sz(self):
        """
        (SZ) Return the 1-halo and 2-halo terms of tSZ power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_sz_1h[index])
            cl['2h'].append(self.tsz.cl_sz_2h[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl

    def cl_sz_at_nu_in_GHz_in_microK2(self,nu_in_GHz):
        frequency_in_Hz = nu_in_GHz*1e9
        T_cmb = self.T_cmb()
        _h_P_=6.62606896e-34
        _k_B_=1.3806504e-23
        Tcmb_gNU = T_cmb*((_h_P_*frequency_in_Hz/(_k_B_*T_cmb))*(1./np.tanh((_h_P_*frequency_in_Hz/(_k_B_*T_cmb))/2.))-4.)
        r = {}
        r['ell'] = np.asarray(self.cl_sz()['ell'])
        r['1h'] = np.asarray(self.cl_sz()['1h'])*Tcmb_gNU**2.
        r['2h'] = np.asarray(self.cl_sz()['2h'])*Tcmb_gNU**2.
        return r

    def get_volume_dVdzdOmega_at_z(self,z):
        return get_volume_at_z(z,&self.ba)

    def get_galaxy_number_counts(self,z):
        return get_galaxy_number_counts(z,&self.tsz)


    def get_cov_N_N(self):
        """
        (class_sz) Return the covariance of cluster number counts
        """
        cl = {}
        cl['poisson'] = []
        cl['hsv'] = []
        cl['r'] = []
        cl['mass_bin_edges'] = []
        for index_M_bins in range(self.tsz.nbins_M):
            cl['poisson'].append(self.tsz.cov_N_N[index_M_bins])
            cl['hsv'].append([])
            cl['r'].append([])
            cl['mass_bin_edges'].append(self.tsz.cov_Y_N_mass_bin_edges[index_M_bins])
            for index_M_bins_prime in range(self.tsz.nbins_M):
                cl['hsv'][index_M_bins].append(self.tsz.cov_N_N_hsv[index_M_bins][index_M_bins_prime])
                cl['r'][index_M_bins].append(self.tsz.cov_N_N_hsv[index_M_bins][index_M_bins_prime]/np.sqrt(self.tsz.cov_N_N[index_M_bins]+self.tsz.cov_N_N_hsv[index_M_bins][index_M_bins])/np.sqrt(self.tsz.cov_N_N[index_M_bins_prime]+self.tsz.cov_N_N_hsv[index_M_bins_prime][index_M_bins_prime]))
        cl['mass_bin_edges'].append(self.tsz.cov_Y_N_mass_bin_edges[self.tsz.nbins_M])

        cl['poisson'] = np.array(cl['poisson'])
        cl['hsv'] = np.array(cl['hsv'])
        cl['r'] = np.array(cl['r'])
        cl['mass_bin_edges'] = np.array(cl['mass_bin_edges'])
        return cl

    def cl_kg_kg(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of galaxy-lensing x galaxy-lensing power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_gallens_gallens_1h[index])
            cl['2h'].append(self.tsz.cl_gallens_gallens_2h[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl

    def cl_yg(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of y x g power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_tSZ_gal_1h[index])
            cl['2h'].append(self.tsz.cl_tSZ_gal_2h[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl


    def cl_kSZ_kSZ_g(self):
        """
        (class_sz) Return the 1-halo, 2-halo and 3-halo terms of kSZ x kSZ x g power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        cl['3h'] = []
        cl['1h (slow)'] = []
        cl['2h (slow)'] = []
        cl['3h (slow)'] = []
        cl['hf'] = []
        cl['covmat'] = []
        cl['lensing term'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h (slow)'].append(self.tsz.cl_kSZ_kSZ_gal_1h[index])
            cl['2h (slow)'].append(self.tsz.cl_kSZ_kSZ_gal_2h[index])
            cl['3h (slow)'].append(self.tsz.cl_kSZ_kSZ_gal_3h[index])
            cl['hf'].append(self.tsz.cl_kSZ_kSZ_gal_hf[index])
            cl['1h'].append(self.tsz.cl_kSZ_kSZ_gal_1h_fft[index])
            cl['2h'].append(self.tsz.cl_kSZ_kSZ_gal_2h_fft[index])
            cl['3h'].append(self.tsz.cl_kSZ_kSZ_gal_3h_fft[index])
            cl['lensing term'].append(self.tsz.cl_kSZ_kSZ_gal_lensing_term[index])
            cl['covmat'].append(self.tsz.cov_ll_kSZ_kSZ_gal[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl

    def cl_kSZ_kSZ_kg(self):
        """
        (class_sz) Return the 1-halo, 2-halo and 3-halo terms of kSZ x kSZ x kg power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        cl['3h'] = []
        cl['hf'] = []
        cl['covmat'] = []
        cl['lensing term'] = []
        for index in range(self.tsz.nlSZ):
            cl['hf'].append(self.tsz.cl_kSZ_kSZ_gallens_hf[index])
            cl['1h'].append(self.tsz.cl_kSZ_kSZ_gallens_1h_fft[index])
            cl['2h'].append(self.tsz.cl_kSZ_kSZ_gallens_2h_fft[index])
            cl['3h'].append(self.tsz.cl_kSZ_kSZ_gallens_3h_fft[index])
            cl['lensing term'].append(self.tsz.cl_kSZ_kSZ_gallens_lensing_term[index])
            cl['covmat'].append(self.tsz.cov_ll_kSZ_kSZ_gallens[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl

    def cl_kSZ_kSZ_kcmb(self):
        """
        (class_sz) Return the 1-halo, 2-halo and 3-halo terms of kSZ x kSZ x kg power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        cl['3h'] = []
        cl['hf'] = []
        cl['covmat'] = []
        cl['lensing term'] = []
        for index in range(self.tsz.nlSZ):
            cl['hf'].append(self.tsz.cl_kSZ_kSZ_lens_hf[index])
            cl['1h'].append(self.tsz.cl_kSZ_kSZ_lens_1h_fft[index])
            cl['2h'].append(self.tsz.cl_kSZ_kSZ_lens_2h_fft[index])
            cl['3h'].append(self.tsz.cl_kSZ_kSZ_lens_3h_fft[index])
            cl['lensing term'].append(self.tsz.cl_kSZ_kSZ_lens_lensing_term[index])
            cl['covmat'].append(self.tsz.cov_ll_kSZ_kSZ_lens[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl



    def cl_te_y_y(self):
        """
        (class_sz) Return Te x y x y power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['teyy'] = []
        for index in range(self.tsz.nlSZ):
            cl['teyy'].append(self.tsz.cl_te_y_y[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl

    def cl_ym(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of y x mu (lensing magnification) power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_tSZ_lensmag_1h[index])
            cl['2h'].append(self.tsz.cl_tSZ_lensmag_2h[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl


    def cl_y_kcmb(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of y x k (cmb lensing) power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_tSZ_lens_1h[index])
            cl['2h'].append(self.tsz.cl_tSZ_lens_2h[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl

    def cl_t2t2f(self):
        cl = {}
        cl['ell'] = []
        cl['t2t2f'] = []
        for index in range(self.tsz.nlSZ):
            cl['ell'].append(self.tsz.ell[index])
            cl['t2t2f'].append(self.tsz.cl_t2t2f[index])
        return cl



    def cl_gg(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of gxg power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        cl['hf'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_gal_gal_1h[index])
            cl['2h'].append(self.tsz.cl_gal_gal_2h[index])
            cl['hf'].append(self.tsz.cl_gal_gal_hf[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl


    def cl_ksz(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of kszxksz power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_kSZ_kSZ_1h[index])
            cl['2h'].append(self.tsz.cl_kSZ_kSZ_2h[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl


    def cl_ggamma(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of galaxy x galaxy lensing power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_gal_gallens_1h[index])
            cl['2h'].append(self.tsz.cl_gal_gallens_2h[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl

    def gamma_ggamma(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of tangential shear
        """
        cl = {}
        cl['thetas'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.gamma_gal_gallens_1h[index])
            cl['2h'].append(self.tsz.gamma_gal_gallens_2h[index])
            cl['thetas'].append(self.tsz.thetas_arcmin[index])
        return cl

    def cl_kg(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of kappa (lensing) x galaxy power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        cl['hf'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_gal_lens_1h[index])
            cl['2h'].append(self.tsz.cl_gal_lens_2h[index])
            cl['hf'].append(self.tsz.cl_gal_lens_hf[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl

    def pk_at_z_hm(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of 3d P(k) matter power spectrum
        """
        cl = {}
        cl['k'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.n_k_for_pk_hm):
            cl['1h'].append(self.tsz.pk_at_z_1h[index])
            cl['2h'].append(self.tsz.pk_at_z_2h[index])
            cl['k'].append(self.tsz.k_for_pk_hm[index])
        return cl


    def bk_at_z_hm(self):
        """
        (class_sz) Return the 1-halo, 2-halo and 3-halo terms of 3d B(k) matter bispectrum
        """
        cl = {}
        cl['k'] = []
        cl['1h'] = []
        cl['2h'] = []
        cl['3h'] = []
        for index in range(self.tsz.n_k_for_pk_hm):
            cl['1h'].append(self.tsz.bk_at_z_1h[index])
            cl['2h'].append(self.tsz.bk_at_z_2h[index])
            cl['3h'].append(self.tsz.bk_at_z_3h[index])
            cl['k'].append(self.tsz.k_for_pk_hm[index])
        return cl


    def bk_ttg_at_z_hm(self):
        """
        (class_sz) Return the 1-halo, 2-halo and 3-halo terms of 3d B(k) ttg bispectrum
        """
        cl = {}
        cl['k'] = []
        cl['1h'] = []
        cl['2h'] = []
        cl['3h'] = []
        for index in range(self.tsz.n_k_for_pk_hm):
            cl['1h'].append(self.tsz.bk_ttg_at_z_1h[index])
            cl['2h'].append(self.tsz.bk_ttg_at_z_2h[index])
            cl['3h'].append(self.tsz.bk_ttg_at_z_3h[index])
            cl['k'].append(self.tsz.k_for_pk_hm[index])
        return cl

    def b_yyy(self):
        """
        (class_sz) Return the 1-halo, 2-halo and 3-halo terms of tsz bispectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        cl['3h'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.b_tSZ_tSZ_tSZ_1halo[index])
            cl['2h'].append(self.tsz.b_tSZ_tSZ_tSZ_2h[index])
            cl['3h'].append(self.tsz.b_tSZ_tSZ_tSZ_3h[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl

    def pk_gg_at_z_hm(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of 3d P(k) gg power spectrum
        """
        cl = {}
        cl['k'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.n_k_for_pk_hm):
            cl['1h'].append(self.tsz.pk_gg_at_z_1h[index])
            cl['2h'].append(self.tsz.pk_gg_at_z_2h[index])
            cl['k'].append(self.tsz.k_for_pk_hm[index])
        return cl

    def pk_b_at_z_2h(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of 3d P(k) bb power spectrum
        """
        cl = {}
        cl['k'] = []
        cl['2h'] = []
        for index in range(self.tsz.n_k_for_pk_hm):
            cl['2h'].append(self.tsz.pk_b_at_z_2h[index])
            cl['k'].append(self.tsz.k_for_pk_hm[index])
        return cl

    def pk_bb_at_z_hm(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of 3d P(k) bb power spectrum
        """
        cl = {}
        cl['k'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.n_k_for_pk_hm):
            cl['1h'].append(self.tsz.pk_bb_at_z_1h[index])
            cl['2h'].append(self.tsz.pk_bb_at_z_2h[index])
            cl['k'].append(self.tsz.k_for_pk_hm[index])
        return cl

    def pk_em_at_z_hm(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of 3d P(k) electron-matter power spectrum
        """
        cl = {}
        cl['k'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.n_k_for_pk_hm):
            cl['1h'].append(self.tsz.pk_em_at_z_1h[index])
            cl['2h'].append(self.tsz.pk_em_at_z_2h[index])
            cl['k'].append(self.tsz.k_for_pk_hm[index])
        return cl

    def cl_kk(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of kappa x kappa (lensing) power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        cl['hf'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_lens_lens_1h[index])
            cl['2h'].append(self.tsz.cl_lens_lens_2h[index])
            cl['hf'].append(self.tsz.cl_lens_lens_hf[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl

    def cib_monopole(self):
        """
        (class_sz) Return the cib monopole as a function of frequency
        """
        cl = {}
        cl['nu'] = []
        cl['I0'] = []
        for index in range(self.tsz.n_frequencies_for_cib):
            cl['nu'].append(self.tsz.frequencies_for_cib[index])
            cl['I0'].append(self.tsz.cib_monopole[index])
        return cl

    def cib_shotnoise(self):
        """
        (class_sz) Return the cib shotnoise as a function of frequency
        """
        cl = {}
        cl['nu'] = []
        cl['shotnoise'] = []
        for id_nu1 in range(self.tsz.cib_frequency_list_num):
            nu1 = self.tsz.cib_frequency_list[id_nu1]
            cl['nu'].append(nu1)
            cl['shotnoise'].append(self.tsz.cib_shotnoise[id_nu1])
        return cl

    def cl_cib_cib(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of cib x cib power spectrum
        """
        cl_cib = {}

        for id_nu1 in range(self.tsz.cib_frequency_list_num):
            for id_nu2 in range(0,id_nu1+1):
                nu1 = self.tsz.cib_frequency_list[id_nu1]
                nu2 = self.tsz.cib_frequency_list[id_nu2]

                cl = {}
                cl['ell'] = []
                cl['1h'] = []
                cl['2h'] = []
                for index in range(self.tsz.nlSZ):
                    cl['1h'].append(self.tsz.cl_cib_cib_1h[id_nu1][id_nu2][index])
                    cl['2h'].append(self.tsz.cl_cib_cib_2h[id_nu1][id_nu2][index])
                    cl['ell'].append(self.tsz.ell[index])
                cl_cib[str(int(nu1))+'x'+str(int(nu2))] = cl
        return cl_cib

    def cl_galn_galn(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of gal x gal power spectrum
        """
        cl_gg = {}

        for id_nu1 in range(self.tsz.galaxy_samples_list_num):
            for id_nu2 in range(0,id_nu1+1):
                nu1 = self.tsz.galaxy_samples_list[id_nu1]
                nu2 = self.tsz.galaxy_samples_list[id_nu2]

                cl = {}
                cl['ell'] = []
                cl['1h'] = []
                cl['2h'] = []
                cl['hf'] = []
                for index in range(self.tsz.nlSZ):
                    cl['1h'].append(self.tsz.cl_ngal_ngal_1h[id_nu1][id_nu2][index])
                    cl['2h'].append(self.tsz.cl_ngal_ngal_2h[id_nu1][id_nu2][index])
                    cl['hf'].append(self.tsz.cl_ngal_ngal_hf[id_nu1][id_nu2][index])
                    cl['ell'].append(self.tsz.ell[index])
                cl_gg[str(int(nu1))+'x'+str(int(nu2))] = cl
        return cl_gg

    def cl_galn_lens(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of gal x lens power spectrum
        """
        cl_gk = {}

        for id_nu1 in range(self.tsz.galaxy_samples_list_num):
            nu1 = self.tsz.galaxy_samples_list[id_nu1]
            cl = {}
            cl['ell'] = []
            cl['1h'] = []
            cl['2h'] = []
            cl['hf'] = []
            for index in range(self.tsz.nlSZ):
                cl['1h'].append(self.tsz.cl_ngal_lens_1h[id_nu1][index])
                cl['2h'].append(self.tsz.cl_ngal_lens_2h[id_nu1][index])
                cl['hf'].append(self.tsz.cl_ngal_lens_hf[id_nu1][index])
                cl['ell'].append(self.tsz.ell[index])
            cl_gk[str(int(nu1))] = cl
        return cl_gk


    def cl_tSZ_cib(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of tSZ x cib power spectrum
        """
        cl_y_cib = {}

        for id_nu1 in range(self.tsz.cib_frequency_list_num):
            nu1 = self.tsz.cib_frequency_list[id_nu1]
            cl = {}
            cl['ell'] = []
            cl['1h'] = []
            cl['2h'] = []
            for index in range(self.tsz.nlSZ):
                cl['1h'].append(self.tsz.cl_tSZ_cib_1h[id_nu1][index])
                cl['2h'].append(self.tsz.cl_tSZ_cib_2h[id_nu1][index])
                cl['ell'].append(self.tsz.ell[index])
            cl_y_cib[str(int(nu1))] = cl
        return cl_y_cib

    def cl_gal_cib(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of galaxy x cib power spectrum
        """
        cl_g_cib = {}

        for id_nu1 in range(self.tsz.cib_frequency_list_num):
            nu1 = self.tsz.cib_frequency_list[id_nu1]
            cl = {}
            cl['ell'] = []
            cl['1h'] = []
            cl['2h'] = []
            for index in range(self.tsz.nlSZ):
                cl['1h'].append(self.tsz.cl_gal_cib_1h[id_nu1][index])
                cl['2h'].append(self.tsz.cl_gal_cib_2h[id_nu1][index])
                cl['ell'].append(self.tsz.ell[index])
            cl_g_cib[str(int(nu1))] = cl
        return cl_g_cib

    def cl_lens_cib(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of lens x cib power spectrum
        """
        cl_phi_cib = {}

        for id_nu1 in range(self.tsz.cib_frequency_list_num):
            nu1 = self.tsz.cib_frequency_list[id_nu1]
            cl = {}
            cl['ell'] = []
            cl['1h'] = []
            cl['2h'] = []
            for index in range(self.tsz.nlSZ):
                cl['1h'].append(self.tsz.cl_lens_cib_1h[id_nu1][index])
                cl['2h'].append(self.tsz.cl_lens_cib_2h[id_nu1][index])
                cl['ell'].append(self.tsz.ell[index])
            cl_phi_cib[str(int(nu1))] = cl
        return cl_phi_cib


    def cl_km(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of kappa (lensing) x mu (lensing magnification) power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_lens_lensmag_1h[index])
            cl['2h'].append(self.tsz.cl_lens_lensmag_2h[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl

    def cl_mm(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of mu x mu (lensing magnification) power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        cl['hf'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_lensmag_lensmag_1h[index])
            cl['2h'].append(self.tsz.cl_lensmag_lensmag_2h[index])
            cl['hf'].append(self.tsz.cl_lensmag_lensmag_hf[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl

    def cl_gm(self):
        """
        (class_sz) Return the 1-halo and 2-halo terms of galaxy x mu (lensing magnification) power spectrum
        """
        cl = {}
        cl['ell'] = []
        cl['1h'] = []
        cl['2h'] = []
        cl['hf'] = []
        for index in range(self.tsz.nlSZ):
            cl['1h'].append(self.tsz.cl_gal_lensmag_1h[index])
            cl['2h'].append(self.tsz.cl_gal_lensmag_2h[index])
            cl['hf'].append(self.tsz.cl_gal_lensmag_hf[index])
            cl['ell'].append(self.tsz.ell[index])
        return cl


    def get_scale_dependent_bias_at_z_and_k(self,z_asked,k_asked,bh):
        return get_scale_dependent_bias_at_z_and_k(z_asked,k_asked,bh,&self.tsz)

    def get_szcounts_dndzdq_at_z_q(self,z_asked,qobs_asked):
        return get_szcounts_dndzdq_at_z_q(z_asked, qobs_asked, &self.tsz)



    def get_params_sz(self):
        """
        (SZ) Return the current parameters
        """
        return self._pars


    def tllprime_sz(self):
        """
        (SZ) Return the trispectrum
        """
        T_ll = defaultdict(list)
        cdef int index_l,index_l_prime
        for index_l in range(self.tsz.nlSZ):
            for index_l_prime in range(index_l+1):
                T_ll[index_l].append(self.tsz.tllprime_sz[index_l][index_l_prime])
        return T_ll


    def A_rs(self):
        """
        (SZ) Return the foreground A_rs coefficient
        """
        return self.tsz.A_rs

    def B_sz(self):
        """
        (SZ) Return the mass bias
        """
        return self.tsz.HSEbias

    def M1SZ(self):
        """
        (SZ) Return the lower mass bound
        """
        return self.tsz.M1SZ


    def ystar(self):
        """
        (SZ) Return ystar
        """
        return self.csz.ystar


    def alpha(self):
        """
        (SZ) Return alpha
        """
        return self.csz.alpha

    def sigmaM(self):
        """
        (SZ) Return sigmaM
        """
        return self.csz.sigmaM


    def sigma8_Pcb(self):
        return self.tsz.sigma8_Pcb

    def A_cib(self):
        """
        (SZ) Return the foreground A_cib coefficient
        """
        return self.tsz.A_cib

    def A_sn(self):
        """
        (SZ) Return the foreground A_cib coefficient
        """
        return self.tsz.cl_gal_gal_A_sn


    def A_ir(self):
        """
        (SZ) Return the foreground A_ir coefficient
        """
        return self.tsz.A_ir


    def A_cn(self):
        """
        (SZ) Return the correlated noise A_cn coefficient
        """
        return self.tsz.A_cn

    def get_te_of_m500c_at_z_arnaud(self,m,z):
        return get_te_of_m500c_at_z_arnaud(m,z,&self.ba,&self.tsz)

    def get_lensing_noise_at_ell(self,l):
        return get_lensing_noise_at_ell(l,&self.tsz)

    def get_te_of_m500c_at_z_lee(self,m,z):
        return get_te_of_m500c_at_z_lee(m,z,&self.ba,&self.tsz)

    def get_f_tinker10_at_nu_and_z(self,nu,z):
        return get_f_tinker10_at_nu_and_z(nu,z,&self.tsz)

    def get_f_tinker08_at_nu_and_z(self,nu,z):
        return get_f_tinker08_at_nu_and_z(nu,z,&self.tsz)


    def get_T10_alpha_at_z(self,z):
        return get_T10_alpha_at_z(z,&self.tsz)

    def get_dndlnM_at_z_and_M(self,z,m):
        #try:
        r = get_dndlnM_at_z_and_M(z,m,&self.tsz)
        #except:
        #print("You are asking for a quantity that hasn't been computed properly.")
        #print("Check your settings and try again?")
        #r = 0.
        return r

    def get_dcib0dz_at_z_and_nu(self,z,nu):
        r = get_dcib0dz_at_z_and_nu(z,nu,&self.tsz)
        return r

    def get_dydz_at_z(self,z):
        r = get_dydz_at_z(z,&self.tsz)
        return r

    def get_mean_y(self):
        r = self.tsz.y_monopole
        return r

    def get_mean_galaxy_bias_at_z(self,z):
        return get_mean_galaxy_bias_at_z(z,&self.tsz)

    def get_gnu_tsz_of_nu_in_ghz(nu_in_ghz,Tcmb):
        r = gnu_tsz_of_nu_in_ghz(nu_in_ghz,Tcmb)
        return r


    def get_hmf_counter_term_nmin_at_z(self,z):
        return get_hmf_counter_term_nmin_at_z(z,&self.tsz)

    def get_hmf_counter_term_b1min_at_z(self,z):
        return get_hmf_counter_term_b1min_at_z(z,&self.tsz)

    def get_hmf_counter_term_b2min_at_z(self,z):
        return get_hmf_counter_term_b2min_at_z(z,&self.tsz)

    def get_gas_density_profile_at_k_M_z(self,l_asked,m_asked,z_asked, include_normalization = 'no'):
        tau_normalization = 1.
        #if (include_normalization == 'yes'):
        #    if (self.tsz.tau_profile == 0): # nfw case
        #        rho0 =  mvir
        #        tau_normalization =  self.ba.Omega0_b/self.tsz.Omega_m_0/self.tsz.mu_e*self.tsz.f_free
        #    elif (self.tsz.tau_profile == 1): # nfw case
        return tau_normalization*get_gas_density_profile_at_k_M_z(l_asked,m_asked,z_asked,&self.tsz)

    def get_rho_2h_at_r_and_m_and_z(self,r_asked,m_asked,z_asked):
        return get_rho_2h_at_r_and_m_and_z(r_asked,m_asked,z_asked,&self.tsz,&self.ba)


    def get_P_delta_at_m_and_z_b12(self,m_asked,z_asked):
        # this is in ev/cm3, see https://arxiv.org/pdf/2202.02275.pdf
        # return get_P_delta_at_m_and_z_b12(r_asked,m_asked,z_asked,&self.tsz,&self.ba)
        r200c = self.get_r_delta_of_m_delta_at_z(200,m_asked,z_asked)
        f_b =  self.get_f_b()
        Eh = self.Hubble(z_asked)/self.Hubble(0)
        P200 = m_asked/r200c*f_b*2.61051e-18*(100.*self.ba.h*Eh)**2.
        return P200


    def get_gas_pressure_2h_at_r_and_m_and_z(self,r_asked,m_asked,z_asked):
        return get_gas_pressure_2h_at_r_and_m_and_z(r_asked,m_asked,z_asked,&self.tsz,&self.ba)


    def get_r_delta_of_m_delta_at_z(self,delta,m_delta,z):
        return (m_delta*3./4./np.pi/delta/self.get_rho_crit_at_z(z))**(1./3.)

    def get_m200m_to_m200c_at_z_and_M(self,z_asked,m_asked):
        return get_m200m_to_m200c_at_z_and_M(z_asked,m_asked,&self.tsz)

    def get_normalization_gas_density_profile(self,z_asked,m_asked):
        return get_normalization_gas_density_profile(z_asked,m_asked,&self.tsz)

    def get_m_to_xout_at_z_and_m(self,z_asked,m_asked):
        return get_m_to_xout_at_z_and_m(z_asked,m_asked,&self.tsz)

    def get_c200m_at_m_and_z_D08(self,M,z):
        return get_c200m_at_m_and_z_D08(M,z)

    def get_c200c_at_m_and_z_D08(self,M,z):
        return get_c200c_at_m_and_z_D08(M,z)

    def get_c200c_at_m_and_z_B13(self,M,z):
        return get_c200c_at_m_and_z_B13(M,z,&self.ba,&self.tsz)

    def get_f_b(self):
        return self.ba.Omega0_b/self.tsz.Omega_m_0

    def get_f_free(self):
        return self.tsz.f_free

    def get_mu_e(self):
        return self.tsz.mu_e

    def get_Omega_m_0(self):
        return self.tsz.Omega_m_0
    def get_Omega_r_0(self):
        return self.tsz.Omega_r_0

    def get_m_nfw(self,x):
        return m_nfw(x)

    def get_rho_crit_at_z(self,z_asked):
        return get_rho_crit_at_z(z_asked,&self.ba,&self.tsz)

    def get_pk_nonlin_at_k_and_z(self,k, z):
        return get_pk_nonlin_at_k_and_z(k,z,&self.ba,&self.pm,&self.nl,&self.tsz)

    def get_pk_lin_at_k_and_z(self,k, z):
        return get_pk_lin_at_k_and_z(k,z,&self.ba,&self.pm,&self.nl,&self.tsz)

    def get_pk_nonlin_at_k_and_z_fast(self,k, z):
        return get_pk_nonlin_at_k_and_z_fast(k,z,&self.ba,&self.pm,&self.nl,&self.tsz)

    def get_pk_lin_at_k_and_z_fast(self,k, z):
        return get_pk_lin_at_k_and_z_fast(k,z,&self.ba,&self.pm,&self.nl,&self.tsz)

    def get_gas_profile_at_x_M_z_b16_200c(self,
                                          r_asked,
                                          m_asked,
                                          z_asked,
                                          c_asked = 0.,
                                          A_rho0 = 4.e3,
                                          A_alpha = 0.88,
                                          A_beta = 3.83,
                                          alpha_m_rho0 = 0.29,
                                          alpha_m_alpha = -0.03,
                                          alpha_m_beta = 0.04,
                                          alpha_z_rho0 = -0.66,
                                          alpha_z_alpha = 0.19,
                                          alpha_z_beta = -0.025,
                                          mcut = 1e14,
                                          alphap_m_rho0 = 0.29,
                                          alphap_m_alpha = -0.03,
                                          alphap_m_beta = 0.04,
                                          alpha_c_rho0 = 0.,
                                          alpha_c_alpha = 0.,
                                          alpha_c_beta = 0.,
                                          gamma = -0.2,
                                          xc = 0.5
                                          ):
        return get_gas_profile_at_x_M_z_b16_200c(r_asked,
                                                 m_asked,
                                                 z_asked,
                                                 c_asked,
                                                 A_rho0,
                                                 A_alpha,
                                                 A_beta,
                                                 alpha_m_rho0,
                                                 alpha_m_alpha,
                                                 alpha_m_beta,
                                                 alpha_z_rho0,
                                                 alpha_z_alpha,
                                                 alpha_z_beta,
                                                 mcut,
                                                 alphap_m_rho0,
                                                 alphap_m_alpha,
                                                 alphap_m_beta,
                                                 alpha_c_rho0,
                                                 alpha_c_alpha,
                                                 alpha_c_beta,
                                                 gamma,
                                                 xc,
                                                 &self.ba,
                                                 &self.tsz)

    def get_pressure_P_over_P_delta_at_x_M_z_b12_200c(self,
                                                        x_asked,
                                                        m_asked,
                                                        z_asked,
                                                        c_asked = 0.,
                                                        A_P0 = 18.1,
                                                        A_xc = 0.497,
                                                        A_beta = 4.35,
                                                        alpha_m_P0 = 0.154,
                                                        alpha_m_xc = -0.00865,
                                                        alpha_m_beta = 0.0393,
                                                        alpha_z_P0 = -0.758,
                                                        alpha_z_xc = 0.731,
                                                        alpha_z_beta = 0.415,
                                                        mcut = 1e14,
                                                        alphap_m_P0 = 0.154,
                                                        alphap_m_xc = -0.00865,
                                                        alphap_m_beta = 0.0393,
                                                        alpha_c_P0 = 0.,
                                                        alpha_c_xc = 0.,
                                                        alpha_c_beta = 0.,
                                                        alpha = 1.,
                                                        gamma = -0.3):
        return  get_pressure_P_over_P_delta_at_x_M_z_b12_200c(x_asked,
                                                              m_asked,
                                                              z_asked,
                                                              c_asked,
                                                              A_P0,
                                                              A_xc,
                                                              A_beta,
                                                              alpha_m_P0,
                                                              alpha_m_xc,
                                                              alpha_m_beta,
                                                              alpha_z_P0,
                                                              alpha_z_xc,
                                                              alpha_z_beta,
                                                              mcut,
                                                              alphap_m_P0,
                                                              alphap_m_xc,
                                                              alphap_m_beta,
                                                              alpha_c_P0,
                                                              alpha_c_xc,
                                                              alpha_c_beta,
                                                              alpha,
                                                              gamma,
                                                              &self.ba,
                                                              &self.tsz)


    def get_pressure_P_over_P_delta_at_x_gnfw_500c(self,
                                                     x_asked,
                                                     P0GNFW = 8.130,
                                                     alphaGNFW = 1.0620,
                                                     betaGNFW = 5.4807,
                                                     gammaGNFW = 0.3292,
                                                     c500 = 1.156):
        return  get_pressure_P_over_P_delta_at_x_gnfw_500c(x_asked,
                                                           P0GNFW,
                                                           alphaGNFW,
                                                           betaGNFW,
                                                           gammaGNFW,
                                                           c500,
                                                           &self.ba,
                                                           &self.tsz)
    def get_dA(self,z):
        """
        angular_distance(z) in Mpc/h


        Parameters
        ----------
        z : float
                Desired redshift
        """
        cdef double tau
        cdef int last_index #junk
        cdef double * pvecback

        pvecback = <double*> calloc(self.ba.bg_size,sizeof(double))

        if background_tau_of_z(&self.ba,z,&tau)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        if background_at_tau(&self.ba,tau,self.ba.long_info,self.ba.inter_normal,&last_index,pvecback)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        D_A = pvecback[self.ba.index_bg_ang_distance]

        free(pvecback)

        return D_A*self.ba.h

    def get_rad_to_arcmin(self,theta_rad):
        return (60.*180.)/np.pi*theta_rad

    def get_truncated_nfw_profile_at_z_k_rd_cd_xout(self,z,k,r_delta,c_delta,xout):
        return evaluate_truncated_nfw_profile(z,k,r_delta,c_delta,xout)

    def get_arcmin_to_rad(self,theta_arcmin):
        return np.pi/(60.*180.)*theta_arcmin

    def get_gas_profile_at_x_M_z_nfw_200m(self,r_asked,m_asked,z_asked):
        return get_gas_profile_at_x_M_z_nfw_200m(r_asked,m_asked,z_asked,&self.ba,&self.tsz)

    def get_gas_profile_at_x_M_z_nfw_200c(self,r_asked,m_asked,z_asked):
        return get_gas_profile_at_x_M_z_nfw_200c(r_asked,m_asked,z_asked,&self.ba,&self.tsz)

    def get_planck_sigma_at_theta500(self, theta500):
        return get_planck_sigma_at_theta500(theta500, &self.tsz)

    def get_szcountsz_sigma_at_theta_in_patch(self, id,theta500):
        return get_szcountsz_sigma_at_theta_in_patch(theta500, id, &self.tsz)

    def get_second_order_bias_at_z_and_nu(self,z,nu):
        return get_second_order_bias_at_z_and_nu(z,nu,&self.tsz,&self.ba)

    def get_first_order_bias_at_z_and_nu(self,z,nu):
        return get_first_order_bias_at_z_and_nu(z,nu,&self.tsz)

    def get_sigma_at_z_and_m(self,z,m):
        return get_sigma_at_z_and_m(z,m,&self.tsz,&self.ba)

    def get_dlnsigma_dlnR_at_z_and_m(self,z,m):
        return get_dlnsigma_dlnR_at_z_and_m(z,m,&self.tsz,&self.ba)

    def get_sigma8_at_z(self,z):
        return get_sigma8_at_z(z,&self.tsz,&self.ba)
    def get_y_at_m_and_z(self,m,z):
        return get_y_at_m_and_z(m, z, &self.tsz, &self.ba)

    def get_theta_at_m_and_z(self,m,z):
        return get_theta_at_m_and_z(m, z, &self.tsz, &self.ba)

    def get_m200m_to_m500c_at_z_and_M(self,z,m):
        return get_m200m_to_m500c_at_z_and_M(z,m,&self.tsz)

    def get_m200c_to_m500c_at_z_and_M(self,z,m):
        return get_m200c_to_m500c_at_z_and_M(z,m,&self.tsz)


    def get_m200c_to_m200m_at_z_and_M(self,z,m):
        return get_m200c_to_m200m_at_z_and_M(z,m,&self.tsz)

    def get_m200m_to_m200c_at_z_and_M(self,z,m):
        return get_m200m_to_m200c_at_z_and_M(z,m,&self.tsz)



    def get_m500c_to_m200c_at_z_and_M(self,z,m):
        return get_m500c_to_m200c_at_z_and_M(z,m,&self.tsz)



    def calcFRel(self,z, M500, obsFreqGHz = 148.0):
        """Calculates relativistic correction to SZ effect at specified frequency, given z, M500 in MSun.

        This assumes the Arnaud et al. (2005) M-T relation, and applies formulae of Itoh et al. (1998)

        As for H13, we return fRel = 1 + delta_SZE (see also Marriage et al. 2011)

        """
        h = self.ba.h
        M500 = M500/h # in Msun
        Ez = Eh = self.Hubble(z)/self.Hubble(0)

        # NOTE: we should define constants somewhere else...
        hplanck=6.63e-34
        kB=1.38e-23
        sigmaT=6.6524586e-29
        me=9.11e-31
        e=1.6e-19
        c=3e8
        TCMB = self.T_cmb()

        # Using Arnaud et al. (2005) M-T to get temperature
        A=3.84e14
        B=1.71
        #TkeV=5.*np.power(((cosmoModel.efunc(z)*M500)/A), 1/B)   # HMF/Astropy
        #TkeV=5.*np.power(((cosmoModel.Ez(z)*M500)/A), 1/B)   # Colossus
        TkeV=5.*np.power(((Ez*M500)/A), 1/B)
        TKelvin=TkeV*((1000*e)/kB)

        # Itoh et al. (1998) eqns. 2.25 - 2.30
        thetae=(kB*TKelvin)/(me*c**2)
        X=(hplanck*obsFreqGHz*1e9)/(kB*TCMB)
        Xtw=X*(np.cosh(X/2.)/np.sinh(X/2.))
        Stw=X/np.sinh(X/2.)

        Y0=-4+Xtw

        Y1=-10. + (47/2.)*Xtw - (42/5.)*Xtw**2 + (7/10.)*Xtw**3 + np.power(Stw, 2)*(-(21/5.) + (7/5.)*Xtw)

        Y2=-(15/2.) +  (1023/8.)*Xtw - (868/5.)*Xtw**2 + (329/5.)*Xtw**3 - (44/5.)*Xtw**4 + (11/30.)*Xtw**5 \
            + np.power(Stw, 2)*(-(434/5.) + (658/5.)*Xtw - (242/5.)*Xtw**2 + (143/30.)*Xtw**3) \
            + np.power(Stw, 4)*(-(44/5.) + (187/60.)*Xtw)

        Y3=(15/2.) + (2505/8.)*Xtw - (7098/5.)*Xtw**2 + (14253/10.)*Xtw**3 - (18594/35.)*Xtw**4 + (12059/140.)*Xtw**5 - (128/21.)*Xtw**6 + (16/105.)*Xtw**7 \
            + np.power(Stw, 2)*(-(7098/10.) + (14253/5.)*Xtw - (102267/35.)*Xtw**2 + (156767/140.)*Xtw**3 - (1216/7.)*Xtw**4 + (64/7.)*Xtw**5) \
            + np.power(Stw, 4)*(-(18594/35.) + (205003/280.)*Xtw - (1920/7.)*Xtw**2 + (1024/35.)*Xtw**3) \
            + np.power(Stw, 6)*(-(544/21.) + (992/105.)*Xtw)

        Y4=-(135/32.) + (30375/128.)*Xtw - (62391/10.)*Xtw**2 + (614727/40.)*Xtw**3 - (124389/10.)*Xtw**4 \
            + (355703/80.)*Xtw**5 - (16568/21.)*Xtw**6 + (7516/105.)*Xtw**7 - (22/7.)*Xtw**8 + (11/210.)*Xtw**9 \
            + np.power(Stw, 2)*(-(62391/20.) + (614727/20.)*Xtw - (1368279/20.)*Xtw**2 + (4624139/80.)*Xtw**3 - (157396/7.)*Xtw**4 \
            + (30064/7.)*Xtw**5 - (2717/7.)*Xtw**6 + (2761/210.)*Xtw**7) \
            + np.power(Stw, 4)*(-(124389/10.) + (6046951/160.)*Xtw - (248520/7.)*Xtw**2 + (481024/35.)*Xtw**3 - (15972/7.)*Xtw**4 + (18689/140.)*Xtw**5) \
            + np.power(Stw, 6)*(-(70414/21.) + (465992/105.)*Xtw - (11792/7.)*Xtw**2 + (19778/105.)*Xtw**3) \
            + np.power(Stw, 8)*(-(682/7.) + (7601/210.)*Xtw)

        deltaSZE=((X**3)/(np.exp(X)-1)) * ((thetae*X*np.exp(X))/(np.exp(X)-1)) * (Y0 + Y1*thetae + Y2*thetae**2 + Y3*thetae**3 + Y4*thetae**4)

        fRel=1+deltaSZE

        return fRel

    def get_ng_bar_at_z(self,double z):
        return evaluate_mean_galaxy_number_density_at_z(z, &self.tsz)

    def get_N_satellites(self,double z,double M_halo,double Nc_mean,double M_min,double alpha_s,double M1_prime):
        Ns = HOD_mean_number_of_satellite_galaxies(z,M_halo,Nc_mean,M_min,alpha_s,M1_prime,&self.tsz,&self.ba)
        return Ns

    def get_N_centrals(self,double z,double M_halo,double M_min,double sigma_log10M,double fc):
        Nc = HOD_mean_number_of_central_galaxies(z,M_halo,M_min,sigma_log10M,fc,&self.tsz,&self.ba)
        return Nc


    def get_unwise_m_min_cut_at_z(self,double z,int sample_id):
        return evaluate_unwise_m_min_cut(z,sample_id,&self.tsz)

    def get_yc_at_m_and_z_H13(self,m,z,A,B):
        Eh = self.Hubble(z)/self.Hubble(0)
        f_rel = self.calcFRel(z, m)
        vec_theta = np.vectorize(self.get_theta_at_m_and_z)
        theta_500 = vec_theta(m,z)
        h = self.ba.h
        QF = np.loadtxt('/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz/sz_auxiliary_files/so_3freqs_020621_theta500arcmin_Q.txt')
        get_Q_of_theta500 = interpolate.interp1d(QF[:,0],QF[:,1],fill_value='extrapolate')
        yp = A*pow(Eh,2.)*pow(m/(3.e14*h),1.+B)*f_rel*get_Q_of_theta500(theta_500)
        return yp

    def _cumulativeNumberDensity(self, z):
        """
        Returns N > M (per cubic Mpc).

        """

        h= self.ba.h
        mm = np.zeros(self.tsz.n_m_dndlnM, dtype=np.double)
        for imp in range(len(mm)):
            mm[imp] = self.tsz.array_m_dndlnM[imp]
        M= np.exp(mm)
        dndlnM=np.vectorize(self.get_dndlnM_at_z_and_M)(z,M)
        dndM=dndlnM/M
        ngtm=integrate.cumtrapz(dndlnM[::-1], np.log(M), initial = 0)[::-1]

        MUpper=np.arange(np.log(M[-1]), np.log(10**18), np.log(M[1])-np.log(M[0]))
        extrapolator=_spline(np.log(M), np.log(dndlnM), k=1)
        MF_extr=extrapolator(MUpper)
        intUpper=integrate.simps(np.exp(MF_extr), dx=MUpper[2] - MUpper[1], even='first')
        ngtm=ngtm+intUpper

        return ngtm


    def getPLog10M(self, z):
        """Returns the log10(mass) probability distribution at the given z, for the logarithmic mass
        binning and mass definition set when the MockSurvey object was constructed.

        Args:
            z (:obj:`float`): Redshift at which to calculate P(log10(M)).

        Returns:
            Array corresponding to the log10(mass) probability distribution.

        """
        numberDensity=self._cumulativeNumberDensity(z)
        mm = np.zeros(self.tsz.n_m_dndlnM, dtype=np.double)
        for imp in range(len(mm)):
            mm[imp] = self.tsz.array_m_dndlnM[imp]
        M= np.exp(mm)
        PLog10M=numberDensity/np.trapz(numberDensity, M)
        return PLog10M

    def getM500_from_y0_H13(self,y0,y0Err,A,B,sigma_int,z,applyMFDebiasCorrection = True):
        # y0 is fixed_y_c*1e-4 in the catalogue
        # log_y0Err is fixed_err_y_c in the catalogue
        # returns:
        # M500*1e14, M500_errPlus*1e14, M500_errMinus*1e14 in Msun
        M500c_zk = np.geomspace(1.00000000e+12,1e16,300)
        log10M500c_zk = np.log10(M500c_zk)
        #log10M500c_zk
        log_y0 = np.log(y0)
        log_y0Err = y0Err/y0
        h = self.ba.h
        m = (10**log10M500c_zk)*h
        log_y0pred = np.log(self.get_yc_at_m_and_z_H13(m,z,A,B)) # feeds in m500c in Msun/h
        Py0GivenM=np.exp(-np.power(log_y0-log_y0pred, 2)/(2*(np.power(log_y0Err, 2)+np.power(sigma_int, 2))))
        Py0GivenM=Py0GivenM/np.trapz(Py0GivenM, log10M500c_zk)
        # Mass function de-bias
        if applyMFDebiasCorrection == True:
            PLog10M_full=self.getPLog10M(z)
            log10Ms = log10M500c_zk
            mm = np.zeros(self.tsz.n_m_dndlnM, dtype=np.double)
            for imp in range(len(mm)):
                mm[imp] = self.tsz.array_m_dndlnM[imp]
            PLog10M =  interpolate.interp1d(mm,PLog10M_full,fill_value='extrapolate')(np.log(10)*log10Ms)
            PLog10M=PLog10M/np.trapz(PLog10M, np.log(10)*log10Ms)/np.log(10)
        else:
            PLog10M=1.0
        Pz = 1.
        P=Py0GivenM*PLog10M*Pz
        PArr = P
        PArr=np.array(PArr)
        P = PArr
        P=P/np.trapz(P, log10M500c_zk)

        calcErrors = True
        log10M = log10M500c_zk
        # Find max likelihood and integrate to get error bars
        tckP=interpolate.splrep(log10M, P)
        fineLog10M=np.linspace(log10M.min(), log10M.max(), 10000)
        fineP=interpolate.splev(fineLog10M, tckP)
        fineP=fineP/np.trapz(fineP, fineLog10M)
        index=np.argmax(fineP)

        clusterLogM500=fineLog10M[index]
        clusterM500=np.power(10, clusterLogM500)/1e14


        for n in range(fineP.shape[0]):
            minIndex=index-n
            maxIndex=index+n
            if minIndex < 0 or maxIndex > fineP.shape[0]:
                # This shouldn't happen; if it does, probably y0 is in the wrong units
                # Previously we threw an exception here, but we can't if using this for forced photometry
                #print("WARNING: outside M500 range - check y0 units or for problem at cluster location in map (if not in forced photometry mode)")
                clusterM500MinusErr=0.
                clusterM500PlusErr=0.
                break
            p=np.trapz(fineP[minIndex:maxIndex], fineLog10M[minIndex:maxIndex])
            if p >= 0.6827:
                clusterLogM500Min=fineLog10M[minIndex]
                clusterLogM500Max=fineLog10M[maxIndex]
                clusterM500MinusErr=(np.power(10, clusterLogM500)-np.power(10, clusterLogM500Min))/1e14
                clusterM500PlusErr=(np.power(10, clusterLogM500Max)-np.power(10, clusterLogM500))/1e14
                break

        return clusterM500*1e14, clusterM500MinusErr*1e14, clusterM500PlusErr*1e14



    def get_bispectrum_f2_kernel(self, double k1, double k2, double k3):
        return bispectrum_f2_kernel(k1, k2, k3)

    def get_nu_at_z_and_m(self,z,m):
        # (delc/sigma)**2
        return get_nu_at_z_and_m(z,m,&self.tsz,&self.ba)
    def get_matter_bispectrum_at_z_effective_approach_smoothed(self,k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,z):
        return get_matter_bispectrum_at_z_effective_approach_smoothed(k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,z,&self.tsz,&self.ba,&self.nl,&self.pm)

    def get_ttg_bispectrum_at_z_effective_approach(self,k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,z):
        return get_ttg_bispectrum_at_z_effective_approach(k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,z,&self.tsz,&self.ba,&self.nl,&self.pm)

    def get_ttg_bispectrum_at_z_tree_level_PT(self,k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,z):
        return get_ttg_bispectrum_at_z_tree_level_PT(k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,z,&self.tsz,&self.ba,&self.nl,&self.pm)

    def get_matter_bispectrum_at_z_effective_approach(self,k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,z):
        return get_matter_bispectrum_at_z_effective_approach(k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,z,&self.tsz,&self.ba,&self.nl,&self.pm)
    def get_matter_bispectrum_at_z_effective_approach_SC(self,k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,z):
        return get_matter_bispectrum_at_z_effective_approach_SC(k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,z,&self.tsz,&self.ba,&self.nl,&self.pm)
    def get_matter_bispectrum_at_z_tree_level_PT(self,k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,z):
        return get_matter_bispectrum_at_z_tree_level_PT(k1_in_h_over_Mpc,k2_in_h_over_Mpc,k3_in_h_over_Mpc,z,&self.tsz,&self.ba,&self.nl,&self.pm)

    def get_nl_index_at_z_and_k(self,z,k1_in_h_over_Mpc):
        return get_nl_index_at_z_and_k(z,k1_in_h_over_Mpc,&self.tsz,&self.nl)
    def get_nl_index_at_z_and_k_no_wiggles(self,z,k1_in_h_over_Mpc):
        return get_nl_index_at_z_and_k_no_wiggles(z,k1_in_h_over_Mpc,&self.tsz,&self.nl)



    def get_vrms2_at_z(self,z):
        return get_vrms2_at_z(z,&self.tsz)


    def get_1e6xdy_from_battaglia_pressure_at_x_z_and_m200c(self,z,m,x):
        return get_1e6xdy_from_battaglia_pressure_at_x_z_and_m200c(z,m,x,&self.ba,&self.tsz)

    def get_1e6xdy_from_gnfw_pressure_at_x_z_and_m500c(self,z,m,x,d):
        return get_1e6xdy_from_gnfw_pressure_at_x_z_and_m500c(z,m,x,d,&self.ba,&self.tsz)

    def szunbinned_loglike(self):
        return self.tsz.szunbinned_loglike
    def szcounts_ntot(self):
        return self.tsz.szcounts_ntot
    def szcounts_ntot_rates_loglike(self, int nq = 500, double qmax = -1):
        # # lndndzdq =  []
        # # print(self.tsz.szcat_size)
        # # first we compute ntot
        zmin = self.tsz.szcounts_fft_z_min
        zmax = self.tsz.szcounts_fft_z_max
        nz = self.tsz.szcounts_fft_nz
        q_threshold = self.tsz.sn_cutoff
        if qmax == -1:
           q_max = self.tsz.szcounts_qmax_fft_padded
        q_arr = np.geomspace(q_threshold, q_max,nq)
        z_arr = np.linspace(zmin,zmax,nz)
        get_dndzdq = np.vectorize(self.get_szcounts_dndzdq_at_z_q)
        Nz = []
        for zp in z_arr:
           Nz.append(np.trapz(get_dndzdq(zp,q_arr)*q_arr,x=np.log(q_arr)))
        Nz = np.asarray(Nz)
        Ntot = np.trapz(Nz,x=z_arr)
        # # compute rates
        rates = []
        for index in range(self.tsz.szcat_size):
            rates.append(self.tsz.szrate[index])
        rates = np.asarray(rates)
        # # compute loglike
        loglike = - Ntot + np.sum(np.log(rates))
        return {'ntot':Ntot,'rates':rates,'loglike':loglike}
        # return Nz

    def dndzdy_theoretical(self):
        """
        (SZ) Return the dndzdy theoretical
        """

        dndzdy =  []
        z_center =  []
        z_edges =  []
        log10y_center =  []
        log10y_edges =  []
        cdef int index_y,index_z

        for index_z in range(self.csz.Nbins_z):
            z_center.append(self.csz.z_center[index_z])
            z_edges.append(self.csz.z_center[index_z]-0.5*self.csz.dz)
            dndzdy_index_z = []
            for index_y in range(self.csz.Nbins_y):
                dndzdy_index_z.append(self.csz.dNdzdy_theoretical[index_z][index_y])
            dndzdy.append(dndzdy_index_z)

        z_edges.append(self.csz.z_center[self.csz.Nbins_z-1]+0.5*self.csz.dz)

        for index_y in range(self.csz.Nbins_y):
            log10y_center.append(self.csz.logy[index_y])
            log10y_edges.append(self.csz.logy[index_y]-0.5*self.csz.dlogy)
        #log10y_center.append(self.csz.logy[self.csz.Nbins_y-1])
        #log10y_edges.append(self.csz.logy[self.csz.Nbins_y-1]-0.5*self.tsz.bin_dlog10_snr_last_bin)
        log10y_edges.append(self.csz.logy[self.csz.Nbins_y-1]+0.5*self.csz.dlogy)
        return {'dndzdy':dndzdy,'z_center':z_center,'z_edges':z_edges,'log10y_center':log10y_center,'log10y_edges':log10y_edges}


    def get_background(self):
        """
        Return an array of the background quantities at all times.

        Parameters
        ----------

        Returns
        -------
        background : dictionary containing background.
        """
        cdef char *titles
        cdef double* data
        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if background_output_titles(&self.ba, titles)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.ba.bt_size

        data = <double*>malloc(sizeof(double)*timesteps*number_of_titles)

        if background_output_data(&self.ba, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.ba.error_message)

        background = {}

        for i in range(number_of_titles):
            background[names[i]] = np.zeros(timesteps, dtype=np.double)
            for index in range(timesteps):
                background[names[i]][index] = data[index*number_of_titles+i]

        free(titles)
        free(data)
        return background

    def get_thermodynamics(self):
        """
        Return the thermodynamics quantities.

        Returns
        -------
        thermodynamics : dictionary containing thermodynamics.
        """
        cdef char *titles
        cdef double* data

        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if thermodynamics_output_titles(&self.ba, &self.th, titles)==_FAILURE_:
            raise CosmoSevereError(self.th.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.th.tt_size

        data = <double*>malloc(sizeof(double)*timesteps*number_of_titles)

        if thermodynamics_output_data(&self.ba, &self.th, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.th.error_message)

        thermodynamics = {}

        for i in range(number_of_titles):
            thermodynamics[names[i]] = np.zeros(timesteps, dtype=np.double)
            for index in range(timesteps):
                thermodynamics[names[i]][index] = data[index*number_of_titles+i]

        free(titles)
        free(data)
        return thermodynamics

    def get_primordial(self):
        """
        Return the primordial scalar and/or tensor spectrum depending on 'modes'.
        'output' must be set to something, e.g. 'tCl'.

        Returns
        -------
        primordial : dictionary containing k-vector and primordial scalar and tensor P(k).
        """
        cdef char *titles
        cdef double* data

        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if primordial_output_titles(&self.pt, &self.pm, titles)==_FAILURE_:
            raise CosmoSevereError(self.pm.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.pm.lnk_size

        data = <double*>malloc(sizeof(double)*timesteps*number_of_titles)

        if primordial_output_data(&self.pt, &self.pm, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.pm.error_message)

        primordial = {}

        for i in range(number_of_titles):
            primordial[names[i]] = np.zeros(timesteps, dtype=np.double)
            for index in range(timesteps):
                primordial[names[i]][index] = data[index*number_of_titles+i]

        free(titles)
        free(data)
        return primordial


    @cython.returns(dict)
    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    @cython.cdivision(True)
    @cython.ccall
    def get_perturbations(self):
        """
        Return scalar, vector and/or tensor perturbations as arrays for requested
        k-values.

        .. note::

            you need to specify both 'k_output_values', and have some
            perturbations computed, for instance by setting 'output' to 'tCl'.

        Returns
        -------
        perturbations : dict of array of dicts
                perturbations['scalar'] is an array of length 'k_output_values' of
                dictionary containing scalar perturbations.
                Similar for perturbations['vector'] and perturbations['tensor'].
        """

        perturbations = {}

        if self.pt.k_output_values_num<1:
            return perturbations

        cdef:
            Py_ssize_t j
            Py_ssize_t i
            Py_ssize_t number_of_titles
            Py_ssize_t timesteps
            list names
            list tmparray
            dict tmpdict
            double[:,::1] data_mv
            double ** thedata
            int * thesizes

        # Doing the exact same thing 3 times, for scalar, vector and tensor. Sorry
        # for copy-and-paste here, but I don't know what else to do.
        for mode in ['scalar','vector','tensor']:
            if mode=='scalar' and self.pt.has_scalars:
                thetitles = <bytes> self.pt.scalar_titles
                thedata = self.pt.scalar_perturbations_data
                thesizes = self.pt.size_scalar_perturbation_data
            elif mode=='vector' and self.pt.has_vectors:
                thetitles = <bytes> self.pt.vector_titles
                thedata = self.pt.vector_perturbations_data
                thesizes = self.pt.size_vector_perturbation_data
            elif mode=='tensor' and self.pt.has_tensors:
                thetitles = <bytes> self.pt.tensor_titles
                thedata = self.pt.tensor_perturbations_data
                thesizes = self.pt.size_tensor_perturbation_data
            else:
                continue
            thetitles = str(thetitles.decode())
            names = thetitles.split("\t")[:-1]
            number_of_titles = len(names)
            tmparray = []
            if number_of_titles != 0:
                for j in range(self.pt.k_output_values_num):
                    timesteps = thesizes[j]//number_of_titles
                    tmpdict={}
                    data_mv = <double[:timesteps,:number_of_titles]> thedata[j]
                    for i in range(number_of_titles):
                        tmpdict[names[i]] = np.asarray(data_mv[:,i])
                    tmparray.append(tmpdict)
            perturbations[mode] = tmparray

        return perturbations

    def get_transfer(self, z=0., output_format='class'):
        """
        Return the density and/or velocity transfer functions for all initial
        conditions today. You must include 'mTk' and/or 'vCTk' in the list of
        'output'. The transfer functions can also be computed at higher redshift z
        provided that 'z_pk' has been set and that 0<z<z_pk.

        Parameters
        ----------
        z  : redshift (default = 0)
        output_format  : ('class' or 'camb') Format transfer functions according to
                         CLASS convention (default) or CAMB convention.

        Returns
        -------
        tk : dictionary containing transfer functions.
        """
        cdef char *titles
        cdef double* data
        cdef char ic_info[1024]
        cdef FileName ic_suffix
        cdef file_format outf

        if (not self.pt.has_density_transfers) and (not self.pt.has_velocity_transfers):
            return {}

        if output_format == 'camb':
            outf = camb_format
        else:
            outf = class_format

        index_md = self.pt.index_md_scalars;
        titles = <char*>calloc(_MAXTITLESTRINGLENGTH_,sizeof(char))

        if perturb_output_titles(&self.ba,&self.pt, outf, titles)==_FAILURE_:
            raise CosmoSevereError(self.pt.error_message)

        tmp = <bytes> titles
        tmp = str(tmp.decode())
        names = tmp.split("\t")[:-1]
        number_of_titles = len(names)
        timesteps = self.pt.k_size[index_md]

        size_ic_data = timesteps*number_of_titles;
        ic_num = self.pt.ic_size[index_md];

        data = <double*>malloc(sizeof(double)*size_ic_data*ic_num)

        if perturb_output_data(&self.ba, &self.pt, outf, <double> z, number_of_titles, data)==_FAILURE_:
            raise CosmoSevereError(self.pt.error_message)

        transfers = {}

        for index_ic in range(ic_num):
            if perturb_output_firstline_and_ic_suffix(&self.pt, index_ic, ic_info, ic_suffix)==_FAILURE_:
                raise CosmoSevereError(self.pt.error_message)
            ic_key = <bytes> ic_suffix

            tmpdict = {}
            for i in range(number_of_titles):
                tmpdict[names[i]] = np.zeros(timesteps, dtype=np.double)
                for index in range(timesteps):
                    tmpdict[names[i]][index] = data[index_ic*size_ic_data+index*number_of_titles+i]

            if ic_num==1:
                transfers = tmpdict
            else:
                transfers[ic_key] = tmpdict

        free(titles)
        free(data)

        return transfers

    def get_current_derived_parameters(self, names):
        """
        get_current_derived_parameters(names)

        Return a dictionary containing an entry for all the names defined in the
        input list.

        Parameters
        ----------
        names : list
                Derived parameters that can be asked from Monte Python, or
                elsewhere.

        Returns
        -------
        derived : dict

        .. warning::

            This method used to take as an argument directly the data class from
            Monte Python. To maintain compatibility with this old feature, a
            check is performed to verify that names is indeed a list. If not, it
            returns a TypeError. The old version of this function, when asked
            with the new argument, will raise an AttributeError.

        """
        if type(names) != type([]):
            raise TypeError("Deprecated")

        derived = {}
        for name in names:
            if name == 'h':
                value = self.ba.h
            elif name == 'H0':
                value = self.ba.h*100
            elif name == 'Omega0_lambda' or name == 'Omega_Lambda':
                value = self.ba.Omega0_lambda
            elif name == 'Omega0_fld':
                value = self.ba.Omega0_fld
            elif name == 'age':
                value = self.ba.age
            elif name == 'conformal_age':
                value = self.ba.conformal_age
            elif name == 'm_ncdm_in_eV':
                value = self.ba.m_ncdm_in_eV[0]
            elif name == 'm_ncdm_tot':
                value = self.ba.Omega0_ncdm_tot*self.ba.h*self.ba.h*93.14
            elif name == 'Neff':
                value = self.ba.Neff
            elif name == 'Omega_m':
                value = self.ba.Omega0_m
            elif name == 'omega_m':
                value = self.ba.Omega0_m*self.ba.h**2
            elif name == 'xi_idr':
                value = self.ba.T_idr/self.ba.T_cmb
            elif name == 'N_dg':
                value = self.ba.Omega0_idr/self.ba.Omega0_g*8./7.*pow(11./4.,4./3.)
            elif name == 'Gamma_0_nadm':
                value = self.th.a_idm_dr*(4./3.)*(self.ba.h*self.ba.h*self.ba.Omega0_idr)
            elif name == 'a_dark':
                value = self.th.a_idm_dr
            elif name == 'tau_reio':
                value = self.th.tau_reio
            elif name == 'z_reio':
                value = self.th.z_reio
            elif name == 'z_rec':
                value = self.th.z_rec
            elif name == 'tau_rec':
                value = self.th.tau_rec
            elif name == 'rs_rec':
                value = self.th.rs_rec
            elif name == 'rs_rec_h':
                value = self.th.rs_rec*self.ba.h
            elif name == 'ds_rec':
                value = self.th.ds_rec
            elif name == 'ds_rec_h':
                value = self.th.ds_rec*self.ba.h
            elif name == 'ra_rec':
                value = self.th.da_rec*(1.+self.th.z_rec)
            elif name == 'ra_rec_h':
                value = self.th.da_rec*(1.+self.th.z_rec)*self.ba.h
            elif name == 'da_rec':
                value = self.th.da_rec
            elif name == 'da_rec_h':
                value = self.th.da_rec*self.ba.h
            elif name == 'z_star':
                value = self.th.z_star
            elif name == 'tau_star':
                value = self.th.tau_star
            elif name == 'rs_star':
                value = self.th.rs_star
            elif name == 'ds_star':
                value = self.th.ds_star
            elif name == 'ra_star':
                value = self.th.ra_star
            elif name == 'da_star':
                value = self.th.da_star
            elif name == 'rd_star':
                value = self.th.rd_star
            elif name == 'z_d':
                value = self.th.z_d
            elif name == 'tau_d':
                value = self.th.tau_d
            elif name == 'ds_d':
                value = self.th.ds_d
            elif name == 'ds_d_h':
                value = self.th.ds_d*self.ba.h
            elif name == 'rs_d':
                value = self.th.rs_d
            elif name == 'rs_d_h':
                value = self.th.rs_d*self.ba.h
            elif name == '100*theta_s':
                value = 100.*self.th.rs_rec/self.th.da_rec/(1.+self.th.z_rec)
            elif name == '100*theta_star':
                value = 100.*self.th.rs_star/self.th.da_star/(1.+self.th.z_star)
            elif name == 'YHe':
                value = self.th.YHe
            elif name == 'n_e':
                value = self.th.n_e
            elif name == 'A_s':
                value = self.pm.A_s
            elif name == 'ln10^{10}A_s':
                value = log(1.e10*self.pm.A_s)
            elif name == 'n_s':
                value = self.pm.n_s
            elif name == 'alpha_s':
                value = self.pm.alpha_s
            elif name == 'beta_s':
                value = self.pm.beta_s
            elif name == 'r':
                # This is at the pivot scale
                value = self.pm.r
            elif name == 'r_0002':
                # at k_pivot = 0.002/Mpc
                value = self.pm.r*(0.002/self.pm.k_pivot)**(
                    self.pm.n_t-self.pm.n_s-1+0.5*self.pm.alpha_s*log(
                        0.002/self.pm.k_pivot))
            elif name == 'n_t':
                value = self.pm.n_t
            elif name == 'alpha_t':
                value = self.pm.alpha_t
            elif name == 'V_0':
                value = self.pm.V0
            elif name == 'V_1':
                value = self.pm.V1
            elif name == 'V_2':
                value = self.pm.V2
            elif name == 'V_3':
                value = self.pm.V3
            elif name == 'V_4':
                value = self.pm.V4
            elif name == 'epsilon_V':
                eps1 = self.pm.r*(1./16.-0.7296/16.*(self.pm.r/8.+self.pm.n_s-1.))
                eps2 = -self.pm.n_s+1.-0.7296*self.pm.alpha_s-self.pm.r*(1./8.+1./8.*(self.pm.n_s-1.)*(-0.7296-1.5))-(self.pm.r/8.)**2*(-0.7296-1.)
                value = eps1*((1.-eps1/3.+eps2/6.)/(1.-eps1/3.))**2
            elif name == 'eta_V':
                eps1 = self.pm.r*(1./16.-0.7296/16.*(self.pm.r/8.+self.pm.n_s-1.))
                eps2 = -self.pm.n_s+1.-0.7296*self.pm.alpha_s-self.pm.r*(1./8.+1./8.*(self.pm.n_s-1.)*(-0.7296-1.5))-(self.pm.r/8.)**2*(-0.7296-1.)
                eps23 = 1./8.*(self.pm.r**2/8.+(self.pm.n_s-1.)*self.pm.r-8.*self.pm.alpha_s)
                value = (2.*eps1-eps2/2.-2./3.*eps1**2+5./6.*eps1*eps2-eps2**2/12.-eps23/6.)/(1.-eps1/3.)
            elif name == 'ksi_V^2':
                eps1 = self.pm.r*(1./16.-0.7296/16.*(self.pm.r/8.+self.pm.n_s-1.))
                eps2 = -self.pm.n_s+1.-0.7296*self.pm.alpha_s-self.pm.r*(1./8.+1./8.*(self.pm.n_s-1.)*(-0.7296-1.5))-(self.pm.r/8.)**2*(-0.7296-1.)
                eps23 = 1./8.*(self.pm.r**2/8.+(self.pm.n_s-1.)*self.pm.r-8.*self.pm.alpha_s)
                value = 2.*(1.-eps1/3.+eps2/6.)*(2.*eps1**2-3./2.*eps1*eps2+eps23/4.)/(1.-eps1/3.)**2
            elif name == 'exp_m_2_tau_As':
                value = exp(-2.*self.th.tau_reio)*self.pm.A_s
            elif name == 'phi_min':
                value = self.pm.phi_min
            elif name == 'phi_max':
                value = self.pm.phi_max
            elif name == 'sigma8':
                value = self.nl.sigma8[self.nl.index_pk_m]
            elif name == 'sigma8_cb':
                value = self.nl.sigma8[self.nl.index_pk_cb]
            elif name == 'k_eq':
                value = self.ba.a_eq*self.ba.H_eq

            else:
                raise CosmoSevereError("%s was not recognized as a derived parameter" % name)
            derived[name] = value
        return derived

    def nonlinear_scale(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_scale(z, z_size)

        Return the nonlinear scale for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] k_nl = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] k_nl_cb = np.zeros(z_size,'float64')
        #cdef double *k_nl
        #k_nl = <double*> calloc(z_size,sizeof(double))
        for index_z in range(z_size):
            if nonlinear_k_nl_at_z(&self.ba,&self.nl,z[index_z],&k_nl[index_z],&k_nl_cb[index_z]) == _FAILURE_:
                raise CosmoSevereError(self.nl.error_message)

        return k_nl

    def nonlinear_scale_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """

make        nonlinear_scale_cb(z, z_size)

        Return the nonlinear scale for all the redshift specified in z, of size

        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] k_nl = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] k_nl_cb = np.zeros(z_size,'float64')
        #cdef double *k_nl
        #k_nl = <double*> calloc(z_size,sizeof(double))
        if (self.ba.Omega0_ncdm_tot == 0.):
            raise CosmoSevereError(
                "No massive neutrinos. You must use pk, rather than pk_cb."
                )
        for index_z in range(z_size):
            if nonlinear_k_nl_at_z(&self.ba,&self.nl,z[index_z],&k_nl[index_z],&k_nl_cb[index_z]) == _FAILURE_:
                raise CosmoSevereError(self.nl.error_message)

        return k_nl_cb

    def nonlinear_hmcode_sigma8(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigma8(z, z_size)

        Return sigma_8 for all the redshift specified in z, of size

        """
        cdef int index_z

        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if nonlinear_hmcode_sigma8_at_z(&self.ba,&self.nl,z[index_z],&sigma_8[index_z],&sigma_8_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_8

    def nonlinear_hmcode_sigma8_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigma8(z, z_size)

        Return sigma_8 for all the redshift specified in z, of size

        """
        cdef int index_z

        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_8_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if nonlinear_hmcode_sigma8_at_z(&self.ba,&self.nl,z[index_z],&sigma_8[index_z],&sigma_8_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_8_cb

    def nonlinear_hmcode_sigmadisp(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigmadisp(z, z_size)

        Return sigma_disp for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if nonlinear_hmcode_sigmadisp_at_z(&self.ba,&self.nl,z[index_z],&sigma_disp[index_z],&sigma_disp_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_disp

    def nonlinear_hmcode_sigmadisp_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigmadisp(z, z_size)

        Return sigma_disp for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if nonlinear_hmcode_sigmadisp_at_z(&self.ba,&self.nl,z[index_z],&sigma_disp[index_z],&sigma_disp_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_disp_cb

    def nonlinear_hmcode_sigmadisp100(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigmadisp100(z, z_size)

        Return sigma_disp_100 for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_100 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_100_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if nonlinear_hmcode_sigmadisp100_at_z(&self.ba,&self.nl,z[index_z],&sigma_disp_100[index_z],&sigma_disp_100_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_disp_100

    def nonlinear_hmcode_sigmadisp100_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigmadisp100(z, z_size)

        Return sigma_disp_100 for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_100 = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_disp_100_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if nonlinear_hmcode_sigmadisp100_at_z(&self.ba,&self.nl,z[index_z],&sigma_disp_100[index_z],&sigma_disp_100_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_disp_100_cb

    def nonlinear_hmcode_sigmaprime(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigmaprime(z, z_size)

        Return sigma_disp for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_prime = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_prime_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if nonlinear_hmcode_sigmaprime_at_z(&self.ba,&self.nl,z[index_z],&sigma_prime[index_z],&sigma_prime_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_prime

    def nonlinear_hmcode_sigmaprime_cb(self, np.ndarray[DTYPE_t,ndim=1] z, int z_size):
        """
        nonlinear_hmcode_sigmaprime(z, z_size)

        Return sigma_disp for all the redshift specified in z, of size
        z_size

        Parameters
        ----------
        z : numpy array
                Array of requested redshifts
        z_size : int
                Size of the redshift array
        """
        cdef int index_z
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_prime = np.zeros(z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] sigma_prime_cb = np.zeros(z_size,'float64')

#        for index_z in range(z_size):
#            if nonlinear_hmcode_sigmaprime_at_z(&self.ba,&self.nl,z[index_z],&sigma_prime[index_z],&sigma_prime_cb[index_z]) == _FAILURE_:
#                raise CosmoSevereError(self.nl.error_message)

        return sigma_prime_cb

    def __call__(self, ctx):
        """
        Function to interface with CosmoHammer

        Parameters
        ----------
        ctx : context
                Contains several dictionaries storing data and cosmological
                information

        """
        data = ctx.get('data')  # recover data from the context

        # If the module has already been called once, clean-up
        if self.state:
            self.struct_cleanup()

        # Set the module to the current values
        self.set(data.cosmo_arguments)
        self.compute(["lensing"])

        # Compute the derived paramter value and store them
        params = ctx.getData()
        self.get_current_derived_parameters(
            data.get_mcmc_parameters(['derived']))
        for elem in data.get_mcmc_parameters(['derived']):
            data.mcmc_parameters[elem]['current'] /= \
                data.mcmc_parameters[elem]['scale']
            params[elem] = data.mcmc_parameters[elem]['current']

        ctx.add('boundary', True)
        # Store itself into the context, to be accessed by the likelihoods
        ctx.add('cosmo', self)

    def get_pk_array(self, np.ndarray[DTYPE_t,ndim=1] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, nonlinear):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=1] pk = np.zeros(k_size*z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] pk_cb = np.zeros(k_size*z_size,'float64')

        if nonlinear == 0:
            nonlinear_pks_at_kvec_and_zvec(&self.ba, &self.nl, pk_linear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        else:
            nonlinear_pks_at_kvec_and_zvec(&self.ba, &self.nl, pk_nonlinear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        return pk

    def get_pk_cb_array(self, np.ndarray[DTYPE_t,ndim=1] k, np.ndarray[DTYPE_t,ndim=1] z, int k_size, int z_size, nonlinear):
        """ Fast function to get the power spectrum on a k and z array """
        cdef np.ndarray[DTYPE_t, ndim=1] pk = np.zeros(k_size*z_size,'float64')
        cdef np.ndarray[DTYPE_t, ndim=1] pk_cb = np.zeros(k_size*z_size,'float64')

        if nonlinear == 0:
            nonlinear_pks_at_kvec_and_zvec(&self.ba, &self.nl, pk_linear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        else:
            nonlinear_pks_at_kvec_and_zvec(&self.ba, &self.nl, pk_nonlinear, <double*> k.data, k_size, <double*> z.data, z_size, <double*> pk.data, <double*> pk_cb.data)

        return pk_cb

    def Omega0_k(self):
        """ Curvature contribution """
        return self.ba.Omega0_k

    def Omega0_cdm(self):
        return self.ba.Omega0_cdm
