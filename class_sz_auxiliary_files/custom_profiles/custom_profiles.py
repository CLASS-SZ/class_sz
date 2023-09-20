import numpy as np
def custom1_ux(x,m,z,classy_sz,*params):
    # m is m_delta and delta is the definition passed for custom1 profile
    # x is r/r_delta

    # delta is collected as:
    delta_def = classy_sz.delta_def_custom1()
    delta = classy_sz.get_delta_from_delta_def_at_z(delta_def,z)

    # to get r_delta:
    r_delta = classy_sz.get_r_delta_of_m_delta_at_z(delta,m,z)

    # to get c_delta:
    c_delta = classy_sz.get_c_delta_at_m_and_z(m,z,delta_def)

    # to get x_out:
    x_out = classy_sz.x_out_custom1()


    # example: nfw profile/cmb lensing
    rs = r_delta/c_delta
    xs = x*r_delta/rs
    prof = xs**-1.*(1.+ xs)**-2.
    norm = m*c_delta**3/(np.log(1.+x_out*c_delta)-x_out*c_delta/(1.+x_out*c_delta))

    return norm*prof

def custom1_W(z,classy_sz,*params):
    # e.g., we want (chi_star - chi)/chi
    chi_star = classy_sz.chi_star()
    chi = classy_sz.get_chi(z)
    H = classy_sz.Hubble(z)
    Rho_crit_0 = classy_sz.Rho_crit_0()

    h = classy_sz.h()
    H0 = classy_sz.Hubble(0)


    # example: cmb lensing
    if z==0.:
        w = 1e-100
    else:
        w = 3./2.*(H0/h)**2/Rho_crit_0*(chi/(1.+z))**-1.*(chi_star-chi)/chi_star # 

    return w
