import numpy as np
from datetime import datetime
import multiprocessing
import time
import functools
import re
from pkg_resources import resource_filename
import os
from scipy import optimize
from scipy.integrate import quad
from scipy.interpolate import interp1d
import math
from numpy import linalg as LA
import mcfit
from mcfit import P2xi
import cosmopower
# import classy_sz as csz



from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import CloughTocher2DInterpolator

kb = 1.38064852e-23 #m2 kg s-2 K-1
clight = 299792458. #m/s
hplanck=6.62607004e-34 #m2 kg / s
firas_T0 = 2.728 #pivot temperature used in the Max Lkl Analysis
firas_T0_bf = 2.725 #best-fitting temperature

Tcmb_uk = 2.7255e6

G_newton = 6.674e-11
rho_crit_over_h2_in_GeV_per_cm3 = 1.0537e-5


nu_21_cm_in_GHz =  1./21.1*clight*1.e2/1.e9
x_21_cm = hplanck*nu_21_cm_in_GHz/kb/firas_T0_bf*1.e9

kappa_c = 2.1419 # 4M_2-3M_c see below eq. 9b of https://arxiv.org/pdf/1506.06582.pdf

beta_mu = 2.1923

G1 = np.pi**2./6
G2 = 2.4041
G3 = np.pi**4/15.
a_rho = G2/G3
alpha_mu = 2.*G1/3./G2 # = 1/beta_mu = π^2/18ζ(3) see eq. 4.15 CUSO lectures.

z_mu_era = 3e5
z_y_era = 5e4
z_reio_min = 6
z_reio_max = 25
z_recombination_min = 800
z_recombination_max = 1500

# Physical constants
# ------------------
# Light speed
class Const:
    c_km_s = 299792.458  # speed of light
    h_J_s = 6.626070040e-34  # Planck's constant
    kB_J_K = 1.38064852e-23  # Boltzmann constant
