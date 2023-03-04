from cobaya.theories.classy import classy
from copy import deepcopy
from typing import NamedTuple, Sequence, Union, Optional
from cobaya.tools import load_module
import logging
import os
import numpy as np

class classy_sz(classy):
    use_class_sz_fast_mode = 0 # this is passed in the yaml file
    def initialize(self):
        """Importing CLASS from the correct path, if given, and if not, globally."""
        self.classy_module = self.is_installed()
        if not self.classy_module:
            raise NotInstalledError(
                self.log, "Could not find CLASS_SZ. Check error message above.")
        from classy_sz import Class, CosmoSevereError, CosmoComputationError
        global CosmoComputationError, CosmoSevereError
        self.classy = Class()
        super(classy,self).initialize()
        # Add general CLASS stuff
        self.extra_args["output"] = self.extra_args.get("output", "")
        if "sBBN file" in self.extra_args:
            self.extra_args["sBBN file"] = (
                self.extra_args["sBBN file"].format(classy=self.path))
        # Derived parameters that may not have been requested, but will be necessary later
        self.derived_extra = []
        self.log.info("Initialized!")

        # # class_sz default params for lkl
        # self.extra_args["output"] = 'tSZ_1h'
        # self.extra_args["multipoles_sz"] = 'P15'
        # self.extra_args['nlSZ'] = 18


    # # here modify if you want to bypass stuff in the class computation
    # def calculate(self, state, want_derived=True, **params_values_dict):
    #     print("Bypassing class_sz")




    def must_provide(self, **requirements):
        if "Cl_sz" in requirements:
            # make sure cobaya still runs as it does for standard classy
            requirements.pop("Cl_sz")
            # specify the method to collect the new observable
            self.collectors["Cl_sz"] = Collector(
                    method="cl_sz", # name of the method in classy.pyx
                    args_names=[],
                    args=[])

        if "sz_binned_cluster_counts" in requirements:
            # make sure cobaya still runs as it does for standard classy
            requirements.pop("sz_binned_cluster_counts")
            # specify the method to collect the new observable
            self.collectors["sz_binned_cluster_counts"] = Collector(
                    method="dndzdy_theoretical", # name of the method in classy.pyx
                    args_names=[],
                    args=[])

        if "sz_unbinned_cluster_counts" in requirements:
            # make sure cobaya still runs as it does for standard classy
            requirements.pop("sz_unbinned_cluster_counts")
            # specify the method to collect the new observable
            self.collectors["sz_unbinned_cluster_counts"] = Collector(
                    method="szcounts_ntot_rates_loglike", # name of the method in classy.pyx
                    args_names=[],
                    args=[])

        super().must_provide(**requirements)

    # get the required new observable
    def get_Cl_sz(self):
        cls = {}
        cls = deepcopy(self._current_state["Cl_sz"])
        return cls

    # get the required new observable
    def get_sz_unbinned_cluster_counts(self):
        cls = deepcopy(self._current_state["sz_unbinned_cluster_counts"])
        # print(cls)
        return cls['loglike']

    # get the required new observable
    def get_sz_binned_cluster_counts(self):
        cls = {}
        cls = deepcopy(self._current_state["sz_binned_cluster_counts"])
        return cls

    def calculate(self, state, want_derived=True, **params_values_dict):
        # Set parameters
        params_values = params_values_dict.copy()
        params_values['ln10^{10}A_s'] = params_values.pop("logA")
        self.set(params_values)
        # Compute!
        try:
            if self.use_class_sz_fast_mode == 1:
                self.classy.compute_class_szfast()
            # self.classy.compute_class_szfast()
            else:
                self.classy.compute()
        # "Valid" failure of CLASS: parameters too extreme -> log and report
        except self.classy_module.CosmoComputationError as e:
            if self.stop_at_error:
                self.log.error(
                    "Computation error (see traceback below)! "
                    "Parameters sent to CLASS: %r and %r.\n"
                    "To ignore this kind of error, make 'stop_at_error: False'.",
                    state["params"], dict(self.extra_args))
                raise
            else:
                self.log.debug("Computation of cosmological products failed. "
                               "Assigning 0 likelihood and going on. "
                               "The output of the CLASS error was %s" % e)
            return False
        # CLASS not correctly initialized, or input parameters not correct
        except self.classy_module.CosmoSevereError:
            self.log.error("Serious error setting parameters or computing results. "
                           "The parameters passed were %r and %r. To see the original "
                           "CLASS' error traceback, make 'debug: True'.",
                           state["params"], self.extra_args)
            raise  # No LoggedError, so that CLASS traceback gets printed
        # Gather products
        for product, collector in self.collectors.items():
            # Special case: sigma8 needs H0, which cannot be known beforehand:
            if "sigma8" in self.collectors:
                self.collectors["sigma8"].args[0] = 8 / self.classy.h()
            method = getattr(self.classy, collector.method)
            arg_array = self.collectors[product].arg_array
            if isinstance(arg_array, int):
                arg_array = np.atleast_1d(arg_array)
            if arg_array is None:
                state[product] = method(
                    *self.collectors[product].args, **self.collectors[product].kwargs)
            elif isinstance(arg_array, Sequence) or isinstance(arg_array, np.ndarray):
                arg_array = np.array(arg_array)
                if len(arg_array.shape) == 1:
                    # if more than one vectorised arg, assume all vectorised in parallel
                    n_values = len(self.collectors[product].args[arg_array[0]])
                    state[product] = np.zeros(n_values)
                    args = deepcopy(list(self.collectors[product].args))
                    for i in range(n_values):
                        for arg_arr_index in arg_array:
                            args[arg_arr_index] = \
                                self.collectors[product].args[arg_arr_index][i]
                        state[product][i] = method(
                            *args, **self.collectors[product].kwargs)
                elif len(arg_array.shape) == 2:
                    if len(arg_array) > 2:
                        raise NotImplementedError("Only 2 array expanded vars so far.")
                    # Create outer combinations
                    x_and_y = np.array(np.meshgrid(
                        self.collectors[product].args[arg_array[0, 0]],
                        self.collectors[product].args[arg_array[1, 0]])).T
                    args = deepcopy(list(self.collectors[product].args))
                    result = np.empty(shape=x_and_y.shape[:2])
                    for i, row in enumerate(x_and_y):
                        for j, column_element in enumerate(x_and_y[i]):
                            args[arg_array[0, 0]] = column_element[0]
                            args[arg_array[1, 0]] = column_element[1]
                            result[i, j] = method(
                                *args, **self.collectors[product].kwargs)
                    state[product] = (
                        self.collectors[product].args[arg_array[0, 0]],
                        self.collectors[product].args[arg_array[1, 0]], result)
                else:
                    raise ValueError("arg_array not correctly formatted.")
            elif arg_array in self.collectors[product].kwargs:
                value = np.atleast_1d(self.collectors[product].kwargs[arg_array])
                state[product] = np.zeros(value.shape)
                for i, v in enumerate(value):
                    kwargs = deepcopy(self.collectors[product].kwargs)
                    kwargs[arg_array] = v
                    state[product][i] = method(
                        *self.collectors[product].args, **kwargs)
            else:
                raise LoggedError(self.log, "Variable over which to do an array call "
                                            f"not known: arg_array={arg_array}")
            if collector.post:
                state[product] = collector.post(*state[product])
        # Prepare derived parameters
        d, d_extra = self._get_derived_all(derived_requested=want_derived)
        if want_derived:
            state["derived"] = {p: d.get(p) for p in self.output_params}
            # Prepare necessary extra derived parameters
        state["derived_extra"] = deepcopy(d_extra)




    @classmethod
    def is_installed(cls, **kwargs):
        return load_module('classy_sz')


# this just need to be there as it's used to fill-in self.collectors in must_provide:
class Collector(NamedTuple):
    method: str
    args: Sequence = []
    args_names: Sequence = []
    kwargs: dict = {}
    arg_array: Union[int, Sequence] = None
    post: Optional[callable] = None
