##################################################### Import system libraries ######################################################
import matplotlib as mpl
mpl.rcdefaults()
mpl.rcParams.update(mpl.rc_params_from_file('python/meine-matplotlibrc'))
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (
    nominal_values as noms,
    std_devs as stds,
)
################################################ Finish importing system libraries #################################################

################################################ Adding subfolder to system's path #################################################
import os
import sys
import inspect
# realpath() will make your script run, even if you symlink it :)
cmd_folder = os.path.realpath(os.path.abspath(
    os.path.split(inspect.getfile(inspect.currentframe()))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

old_path = os.getcwd()
os.chdir("python")

############################################# Finish adding subfolder to system's path #############################################

################################################ Finish importing custom libraries #################################################


def plot_test(file_in, file_out, file_r_out, realteil, iterativ):
    plt.clf()
    markers = ['o', 'D', 'v', 's', '*', 'p']
    Kx, err, rate, Nx = np.genfromtxt(file_in, unpack=True)

    a = int(min(Nx))
    b = int(max(Nx))+1
    r = np.zeros(b-a)
    s = np.zeros(b-a)
    assert (b-a) == len(markers), 'b-a='+str(b-a)+'\t no_markers='+str(len(markers))
    f = open(file_r_out, 'w')

    for i, mark in zip(range(a, b), markers):
        map = Nx == i
        temp = rate[map]
        r[i] = round(np.mean(temp[1:-1]), 2)
        s[i] = round(np.std(temp[1:-1]), 2)
        if not(i == b-1):
            f.write(str(r[i])+'\pm'+str(s[i]) + ',\; ')
        else:
            f.write(str(r[i])+'\pm'+str(s[i]))
        plt.plot(Kx[map], err[map]*1e-24, '--', marker=mark, markersize=8,          # mit N_D
                 fillstyle='none', label=r'$N={}$'.format(i))
        # plt.plot(Kx[map], err[map], '--', marker=mark, markersize=8,              # ohne N_D
        # fillstyle='none', label=r'$N={}$'.format(i))

    f.close()
    # plt.text(0.5, 0.1, r'$r={}$'.format(r), rotation=-10.)

    plt.xscale('log', basex=2)
    plt.yscale('log')
    plt.xlabel(r'$K_x$')
    if realteil:
        # plt.ylabel(r'$\left\lVert \, \operatorname{Re}(u-u^h) \, \right\rVert_{L^2(\Omega)}$')
        if iterativ:
            plt.ylabel(r'$e_{\mathcal{T},r}^{\alpha} / N_D$')                       # mit N_D
            # plt.ylabel(r'$e_{\mathcal{T},r}^{\alpha}$')                           # ohne N_D
        else:
            plt.ylabel(r'$e_r^{\alpha} / N_D$')                                     # mit N_D
            # plt.ylabel(r'$e_r^{\alpha}$')                                         # ohne N_D
    else:
        # plt.ylabel(r'$\left\lVert \, \operatorname{Im}(u-u^h) \, \right\rVert_{L^2(\Omega)}$')
        if iterativ:
            plt.ylabel(r'$e_{\mathcal{T},i}^{\alpha} / N_D$')                       # mit N_D
            # plt.ylabel(r'$e_{\mathcal{T},i}^{\alpha}$')                           # ohne N_D
        else:
            plt.ylabel(r'$e_i^{\alpha} / N_D$')                                     # mit N_D
            # plt.ylabel(r'$e_i^{\alpha}$')                                         # ohne N_D
    plt.legend(loc='best', prop={'size': 9}, markerscale=0.8, handlelength=3)
    plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
    plt.savefig(file_out)


# file_out = '../plots/test1_r.pdf'
# file_in = 'data/test1_real.txt'
# file_r_out = 'test1_r_r.tex'
# plot_test(file_in, file_out, file_r_out, True, True)
# file_out = '../plots/test1_i.pdf'
# file_in = 'data/test1_imag.txt'
# file_r_out = 'test1_r_i.tex'
# plot_test(file_in, file_out, file_r_out, False, True)
# file_rel_out = '../plots/test1_r_rel.pdf'
# file_rel_in = 'data/test1_rel_real.txt'
# file_r_rel_out = 'test1_r_r_rel.tex'
# plot_test(file_rel_in, file_rel_out, file_r_rel_out, True, False)
# file_rel_out = '../plots/test1_i_rel.pdf'
# file_rel_in = 'data/test1_rel_imag.txt'
# file_r_rel_out = 'test1_r_i_rel.tex'
# plot_test(file_rel_in, file_rel_out, file_r_rel_out, False, False)

# file_out = '../plots/test5/convGL0.pdf'
# file_in = 'data/test5/test5_conv_GL0.txt'
# file_r_out = 'test5_convGL0_r.tex'
# plot_test(file_in, file_out, file_r_out, True, True)
# file_out = '../plots/test5/convGL1.pdf'
# file_in = 'data/test5/test5_conv_GL1.txt'
# file_r_out = 'test5_convGL1_r.tex'
# plot_test(file_in, file_out, file_r_out, True, True)
# file_out = '../plots/test5/nonConvGL0.pdf'
# file_in = 'data/test5/test5_nonConv_GL0.txt'
# file_r_out = 'test5_nonConvGL0_r.tex'
# plot_test(file_in, file_out, file_r_out, True, True)
# file_out = '../plots/test5/nonConvGL1.pdf'
# file_in = 'data/test5/test5_nonConv_GL1.txt'
# file_r_out = 'test5_nonConvGL1_r.tex'
# plot_test(file_in, file_out, file_r_out, True, True)
# file_out = '../plots/test5/convGL1_i.pdf'
# file_in = 'data/test5/test5_conv_GL1_i.txt'
# file_r_out = 'test5_convGL1_i.tex'
# plot_test(file_in, file_out, file_r_out, False, True)

# file_out = '../plots/test3/real_GL1.pdf'
# file_in = 'data/test3_real_GL1.txt'
# file_r_out = 'test3_real_r.tex'
# plot_test(file_in, file_out, file_r_out, True, False)
# file_out = '../plots/test3/imag_GL1.pdf'
# file_in = 'data/test3_imag_GL1.txt'
# file_r_out = 'test3_imag_r.tex'
# plot_test(file_in, file_out, file_r_out, False, False)


# convenience file writing for standard make Files
f = open('.pysuccess', 'w')
f.write('MarktkraM')
f.close()

os.chdir(old_path)
