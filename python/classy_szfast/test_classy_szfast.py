# import classy_szfast as cszfast
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# from classy_szfast import tabulate_gas_pressure_profile_k,get_gas_pressure_profile_x_parallel
from classy_sz import Class
from classy_szfast import classy_szfast
import time


# cszfast = classy_szfast()


# exit()



params_settings = {# LambdaCDM parameters
                   # last column of Table 1 of https://arxiv.org/pdf/1807.06209.pdf:
                   'H0': 67.37,
                   'omega_b':0.02233,
                   'omega_cdm': 0.1198,
                   'ln10^{10}A_s':3.043,
                   'tau_reio': 0.0540,
                   'n_s': 0.9652,
                   'ndim_redshifts' : 150,
                   }

cszfast = classy_szfast(**params_settings)
print(cszfast.cp_path_to_cosmopower_organization)

l_max_scalars = 11000

mp = 'lcdm'

params_cp = {}
for key,value in params_settings.items():
    params_cp[key] = [value]



predicted_tt_spectrum = cszfast.cp_tt_nn[mp].ten_to_predictions_np(params_cp)
predicted_te_spectrum = cszfast.cp_te_nn[mp].predictions_np(params_cp)
predicted_ee_spectrum = cszfast.cp_ee_nn[mp].ten_to_predictions_np(params_cp)



l = np.arange(2,l_max_scalars+1)
print(len(l))
plt.figure()
plt.plot(l,predicted_tt_spectrum[0],ls='-',lw=1.,label='tt')
plt.plot(l,np.abs(predicted_te_spectrum[0]),ls='-',lw=1.,label='te')
plt.plot(l,predicted_ee_spectrum[0],ls='-',lw=1.,label='ee')

plt.legend(framealpha=1)
plt.xscale('log')
plt.yscale('log')

plt.ylabel('Dimensionless Dl\'s')
plt.xlabel('l')
plt.grid(which='both',alpha=0.1)
# plt.show(block=True)
plt.show(block=False)


print('starting computations')

print('calculating cmb')
start = time.time()
cszfast.calculate_cmb(**params_settings)
end = time.time()
print('cmb calculation took:',end-start)

print('get cmb cls')
start = time.time()
cmb_cls = cszfast.get_cmb_cls(ell_factor = False)
end = time.time()
print('get cmb cls took:',end-start)
print(cmb_cls)

plt.figure()
plt.plot(cmb_cls['ell'],cmb_cls['tt'],ls='-',lw=1.,label='tt')
plt.plot(cmb_cls['ell'],np.abs(cmb_cls['te']),ls='-',lw=1.,label='te')
plt.plot(cmb_cls['ell'],cmb_cls['ee'],ls='-',lw=1.,label='ee')
plt.plot(cmb_cls['ell'],cmb_cls['pp'],ls='-',lw=1.,label='pp')

plt.legend(framealpha=1)
plt.xscale('log')
plt.yscale('log')

plt.ylabel('cl\'s')
plt.xlabel('l')
plt.grid(which='both',alpha=0.1)
# plt.show(block=True)
plt.show(block=False)


print('calculating pkl')
start = time.time()
cszfast.calculate_pkl(**params_settings)
end = time.time()
print('pk calculation took:',end-start)

print('calculating pknl')
start = time.time()
cszfast.calculate_pknl(**params_settings)
end = time.time()
print('pknl calculation took:',end-start)

# exit(0)

print('getting pkl linear ')
start = time.time()
k_asked = 50.
z_asked = 0.1
pkl = cszfast.get_pkl_at_k_and_z(k_asked,z_asked)
end = time.time()
print(pkl)
print('get pkl linear took:',end-start)

print('getting pkl cloughtocher ')
start = time.time()
k_asked = 50.
z_asked = 0.1
pkl = cszfast.get_pkl_at_k_and_z(k_asked,z_asked)
end = time.time()
print(pkl)
print('get pkl cloughtocher took:',end-start)


print('getting pknl linear ')
start = time.time()
k_asked = 50.
z_asked = 0.1
pknl = cszfast.get_pknl_at_k_and_z(k_asked,z_asked)
end = time.time()
print(pknl)
print('get pknl linear took:',end-start)

print('getting pknl cloughtocher ')
start = time.time()
k_asked = 50.
z_asked = 0.1
pknl = cszfast.get_pknl_at_k_and_z(k_asked,z_asked)
end = time.time()
print(pknl)
print('get pknl cloughtocher took:',end-start)


z_arr = [0.,2.5,5.]
k_arr = np.geomspace(1e-3,1e1,100)
print('getting pks at 3 zs ')
start = time.time()
for zp in z_arr:
    pknl = cszfast.get_pknl_at_k_and_z(k_arr,zp)
    pkl = cszfast.get_pkl_at_k_and_z(k_arr,zp)
end = time.time()
print('get pks at 3 zs took:',end-start)

plt.figure()
for zp in z_arr:
    pknl = cszfast.get_pknl_at_k_and_z(k_arr,zp)
    pkl = cszfast.get_pkl_at_k_and_z(k_arr,zp)
    # print(pkl,pknl)

    plt.plot(k_arr,pkl,ls='-',lw=1.,label=zp)
    plt.plot(k_arr,pknl,ls='--',lw=1.,label=zp)

plt.legend(framealpha=1)
plt.xscale('log')
plt.yscale('log')

plt.ylabel('pk\'s')
plt.xlabel('k')
plt.grid(which='both',alpha=0.1)
# plt.show(block=True)
# plt.show(block=False)

print('calculate hubble')
cszfast.calculate_hubble(**params_settings)
z_arr = np.linspace(0.,2.,500)
start = time.time()
hz = cszfast.get_hubble(z_arr)
end = time.time()
print(hz)
print('finish calculate hubble:',end-start)
print('calculate chi')
start = time.time()
cszfast.calculate_chi(**params_settings)
end = time.time()
print('finnish calculate chi:',end-start)
chiz = cszfast.get_chi(z_arr)
print('chiz',chiz)

print('calculate pressure profile')
z_asked,m_asked,x_asked = 0.2,3e14,0.3
px = cszfast.get_gas_pressure_profile_x(z_asked,m_asked,x_asked)
print(px)
print('end calculate pressure profile')

print('calculate pressure profile')
z_asked,m_asked,x_asked = 0.2,3e14,np.geomspace(1e-3,1e2,500)
start = time.time()
px = cszfast.get_gas_pressure_profile_x(z_asked,m_asked,x_asked)
end = time.time()
print(px)
print('end calculate pressure profile:',end-start)

print('tabulate pressure profile')
# z = 1
# m = 1
# cszfast.tabulate_gas_pressure_profile_k()

print('end tabulate pressure profile:',end-start)

print('tabulate sigma')
start = time.time()
cszfast.calculate_sigma()
end = time.time()
print('end tabulate sigma:',end-start)
# print(cszfast.cszfast_pk_grid_sigma2,np.shape(cszfast.cszfast_pk_grid_sigma2))
# print(cszfast.cszfast_pk_grid_sigma2_flat,np.shape(cszfast.cszfast_pk_grid_sigma2_flat))
print(cszfast.cszfast_pk_grid_dsigma2_flat,np.shape(cszfast.cszfast_pk_grid_dsigma2_flat))


print('getting sigma ')
start = time.time()
k_asked = 50.
z_asked = 0.1
pknl = cszfast.get_pknl_at_k_and_z(k_asked,z_asked)
end = time.time()
print(pknl)
print('get pknl cloughtocher took:',end-start)
