import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def plot_Cell(ell, one_halo, two_halo, nu1, nu2, freq, mod):
    total = one_halo + two_halo
    f, ax = plt.subplots(figsize=(8, 8))
    ax.plot(ell, np.abs(one_halo[nu1, nu2, :]), 'b-.', label='1-halo')
    ax.plot(ell, np.abs(two_halo[nu1, nu2, :]), 'b--', label='2-halo')
    ax.plot(ell, np.abs(total[nu1, nu2, :]), 'b', label='total')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(loc='upper right', prop={'size': 12}, frameon=False)
    ax.set_xticks([100, 500, 1000])
    # ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_ylabel(r'$\mathrm{C_l}\: [\mathrm{Jy}^2\: \mathrm{sr}^{-1}]$', fontsize=14)
    ax.set_xlabel(r'Multipole' r'$\;\ell$', fontsize=14)
    ax.set_title(r''+mod+' %s x %s GHz' % (freq[nu1], freq[nu2]))
