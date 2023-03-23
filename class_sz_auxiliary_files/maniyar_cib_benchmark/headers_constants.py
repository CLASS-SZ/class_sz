import numpy as np
from numpy import sin, cos
from astropy.cosmology import Planck15
import scipy.integrate as intg
import astropy.units as u
from scipy.special import *
from numpy import log
import scipy.constants as con
from astropy.io import fits
import hmf_unfw_bias
from scipy.interpolate import UnivariateSpline, interp1d, interp2d, RectBivariateSpline


KC = 1.0e-10  # Kennicutt constant for Chabrier IMF
T_cmb = 2.725  # CMB temp
c_light = 299792458e-3  # Km/sec

m_e = con.electron_mass  # Kg  9.10938356e-31 Kg
h_p = con.h  # Plank's constant 6.62607004e-34 SI units
k_B = con.k  # Boltzmann constant 1.38064852e-23 SI
sig_T = con.physical_constants['Thomson cross section'][0]  # m^2 6.6524587158e-29
H_0 = Planck15.H0.value  # Hubble's constant today
M_s = 1.99e30  # Mass of the Sun Kg

eV_to_J = 1.6e-19
cm_to_m = 1e-2
Mpc_to_m = 3.086e22  # Mpc to m
Km_to_m = 1e3
ghz = 1e9

w_jy = 1e26  # Watt to Jy
nW = 1e9

cosmo = Planck15
