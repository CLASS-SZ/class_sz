import numpy as np

def custom1_b(z,classy_sz,*params):
    # # e.g., we want (chi_star - chi)/chi
    # chi_star = classy_sz.chi_star()
    # chi = classy_sz.get_chi(z)
    # H = classy_sz.Hubble(z)
    # Rho_crit_0 = classy_sz.Rho_crit_0()

    # h = classy_sz.h()
    # H0 = classy_sz.Hubble(0)


    # # example: cmb lensing
    # if z==0.:
    #     w = 1e-100
    # else:
    #     w = 3./2.*(H0/h)**2/Rho_crit_0*(chi/(1.+z))**-1.*(chi_star-chi)/chi_star # 

    return 2.# + 2*z
