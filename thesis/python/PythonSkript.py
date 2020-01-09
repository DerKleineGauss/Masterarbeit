##################################################### Import system libraries ######################################################
import matplotlib as mpl
# mpl.rcdefaults()
# mpl.rcParams.update(mpl.rc_params_from_file('python/meine-matplotlibrc'))
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


def plot_test(file_in, file_out, file_r_out, realteil):
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
            f.write(str(r[i])+'\pm'+str(s[i]) + ', ')
        else:
            f.write(str(r[i])+'\pm'+str(s[i]))
        plt.plot(Kx[map], err[map], '--', marker=mark, markersize=8,
                 fillstyle='none', label=r'$N_x={}$'.format(i))

    f.close()
    # plt.text(0.5, 0.1, r'$r={}$'.format(r), rotation=-10.)

    plt.xscale('log', basex=2)
    plt.yscale('log')
    plt.xlabel(r'$K_x$')
    if realteil:
        plt.ylabel(r'$\left\lVert \, \Re(u-u^h) \, \right\rVert_{L^2(\Omega)}$')
    else:
        plt.ylabel(r'$\left\lVert \, \Im(u-u^h) \, \right\rVert_{L^2(\Omega)}$')
    plt.legend(loc='best', prop={'size': 9}, markerscale=0.8, handlelength=3)
    plt.text(0.5, 0.5, 'AHHHH')
    plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
    plt.savefig(file_out)


# file_out = '../plots/test1_r.pdf'
# file_in = 'data/test1_real.txt'
# file_r_out = 'test1_r_r.tex'
# plot_test(file_in, file_out, file_r_out, True)
# file_out = '../plots/test1_i.pdf'
# file_in = 'data/test1_imag.txt'
# file_r_out = 'test1_r_i.tex'
# plot_test(file_in, file_out, file_r_out, False)
# file_rel_out = '../plots/test1_r_rel.pdf'
# file_rel_in = 'data/test1_rel_real.txt'
# file_r_rel_out = 'test1_r_r_rel.tex'
# plot_test(file_rel_in, file_rel_out, file_r_rel_out, True)
# file_rel_out = '../plots/test1_i_rel.pdf'
# file_rel_in = 'data/test1_rel_imag.txt'
# file_r_rel_out = 'test1_r_i_rel.tex'
# plot_test(file_rel_in, file_rel_out, file_r_rel_out, False)

from scipy.io import loadmat
annots = loadmat('data/test6Results_Lr66.mat')
x = annots['test6Results'][0][0][0]
V = annots['test6Results'][0][0][1]
I_over_x = annots['test6Results'][0][0][2]
I = annots['test6Results'][0][0][3]
plt.plot(V, I)
# plt.savefig('out.pdf')
print(V)
print(I)
plt.show()


################################ FREQUENTLY USED CODE ################################
#
########## IMPORT ##########
# t, U, U_err = np.genfromtxt('messdaten/data.txt', unpack=True)
# t *= 1e-3


########## ERRORS ##########
# R_unc = ufloat(R[0],R[2])
# U = 1e3 * unp.uarray(U, U_err)
# Rx_mean = np.mean(Rx)                 # Mittelwert und syst. Fehler
# Rx_mean_with_error = mean(Rx, 0)      # unp.uarray mit Fehler und Fehler des Mittelwertes, die 0 gibt an, dass in einem R^2 array jeweils die Zeilen gemittelt werden sollen
# Rx_mean_err = MeanError(noms(Rx))     # nur der Fehler des Mittelwertes
#
# Relative Fehler zum späteren Vergleich in der Diskussion
# RelFehler_G = (G_mess - G_lit) / G_lit
# RelFehler_B = (B_mess - B_lit) / B_lit
# write('build/RelFehler_G.tex', make_SI(RelFehler_G*100, r'\percent', figures=1))
# write('build/RelFehler_B.tex', make_SI(RelFehler_B*100, r'\percent', figures=1))


########## CURVE FIT ##########
# def f(t, a, b, c, d):
#     return a * np.sin(b * t + c) + d
#
# params = ucurve_fit(f, t, U, p0=[1, 1e3, 0, 0])   # p0 bezeichnet die Startwerte der zu fittenden Parameter
# params = ucurve_fit(reg_linear, x, y)             # linearer Fit
# params = ucurve_fit(reg_quadratic, x, y)          # quadratischer Fit
# params = ucurve_fit(reg_cubic, x, y)              # kubischer Fit
# a, b = params
# write('build/parameter_a.tex', make_SI(a * 1e-3, r'\kilo\volt', figures=1))       # type in Anz. signifikanter Stellen
# write('build/parameter_b.tex', make_SI(b * 1e-3, r'\kilo\hertz', figures=2))      # type in Anz. signifikanter Stellen

########## MAKE_SI ##########
# make_SI(m_eff_n*1e32, r'\kilo\gramm', figures=1, exp='e-32')
# ufloats, die einen exponenten besitzen, werden als (2.3+/-0.3)e-32 dargestellt, das Problem ist hier der Exponent,
# man sollte sich bei einem Fehler der Form:
# TypeError: non-empty format string passed to object.__format__
# den Wert über print() ausgeben lassen und dann obige Modifikation (links *1e32, und dann exp='e-32') nutzen.
# (hierzu gibt es momentan /Juni2017 ein issue)


########## PLOTTING ##########
# plt.clf()                   # clear actual plot before generating a new one
#
# automatically choosing limits with existing array T1
# t_plot = np.linspace(np.amin(T1), np.amax(T1), 100)
# plt.xlim(t_plot[0]-1/np.size(T1)*(t_plot[-1]-t_plot[0]), t_plot[-1]+1/np.size(T1)*(t_plot[-1]-t_plot[0]))
#
# hard coded limits
# t_plot = np.linspace(-0.5, 2 * np.pi + 0.5, 1000) * 1e-3
#
# standard plotting
# plt.plot(t_plot * 1e3, f(t_plot, *noms(params)) * 1e-3, 'b-', label='Fit')
# plt.plot(t * 1e3, U * 1e3, 'rx', label='Messdaten')
# plt.errorbar(B * 1e3, noms(y) * 1e5, fmt='rx', yerr=stds(y) * 1e5, label='Messdaten')        # mit Fehlerbalken
# plt.xscale('log')                                                                            # logarithmische x-Achse
# plt.xlim(t_plot[0] * 1e3, t_plot[-1] * 1e3)
# plt.xlabel(r'$t \:/\: \si{\milli\second}$')
# plt.ylabel(r'$U \:/\: \si{\kilo\volt}$')
# plt.legend(loc='best')
# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/aufgabenteil_a_plot.pdf')


########## WRITING TABLES ##########
# IF THERE IS ONLY ONE COLUMN IN A TABLE (workaround):
# a=np.array([Wert_d[0]])
# b=np.array([Rx_mean])
# c=np.array([Rx_mean_err])
# d=np.array([Lx_mean*1e3])
# e=np.array([Lx_mean_err*1e3])
#
# write('build/Tabelle_b.tex', make_table([a,b,c,d,e],[0, 1, 0, 1, 1]))     # fehlerbehaftete arrays (uarrays) sollten rechts eine 1 bekommen (1 signifikante Stelle)
# write('build/Tabelle_b_texformat.tex', make_full_table(
#     caption = 'Messdaten Kapazitätsmessbrücke.',
#     label = 'table:A2',
#     source_table = 'build/Tabelle_b.tex',
#     stacking = [1,2,3,4,5],              # Hier aufpassen: diese Zahlen bezeichnen diejenigen resultierenden Spaltennummern, die Multicolumns sein sollen (beginnend bei 0)
#     units = ['Wert',
#     r'$C_2 \:/\: \si{\nano\farad}$',
#     r'$R_2 \:/\: \si{\ohm}$',
#     r'$R_3 / R_4$', '$R_x \:/\: \si{\ohm}$',
#     r'$C_x \:/\: \si{\nano\farad}$'],
#     replaceNaN = True,                      # default = false
#     replaceNaNby = 'not a number'))         # default = '-'
#
# Aufsplitten von Tabellen, falls sie zu lang sind
# t1, t2 = np.array_split(t * 1e3, 2)
# U1, U2 = np.array_split(U * 1e-3, 2)
# write('build/loesung-table.tex', make_table([t1, U1, t2, U2], [3, None, 3, None]))  # type in Nachkommastellen
#
# Verschmelzen von Tabellen (nur Rohdaten, Anzahl der Spalten muss gleich sein)
# write('build/Tabelle_b_composed.tex', make_composed_table(['build/Tabelle_b_teil1.tex','build/Tabelle_b_teil2.tex']))


########## ARRAY FUNCTIONS ##########
# np.arange(2,10)                   # Erzeugt aufwärts zählendes Array von 2 bis 10
# np.zeros(15)                      # Erzeugt Array mit 15 Nullen
# np.ones(15)                       # Erzeugt Array mit 15 Einsen
#
# np.amin(array)                    # Liefert den kleinsten Wert innerhalb eines Arrays
# np.argmin(array)                  # Gibt mir den Index des Minimums eines Arrays zurück
# np.amax(array)                    # Liefert den größten Wert innerhalb eines Arrays
# np.argmax(array)                  # Gibt mir den Index des Maximums eines Arrays zurück
#
# a1,a2 = np.array_split(array, 2)  # Array in zwei Hälften teilen
# np.size(array)                    # Anzahl der Elemente eines Arrays ermitteln
#
# np.delte(A,3)                     # liefert ein Array, in dem der Eintrag mit Index 3 des arrays
# A gelöscht wurde (nachfolgende indices werden aufgerückt)


########## ARRAY INDEXING ##########
# y[n - 1::n]                       # liefert aus einem Array jeden n-ten Wert als Array


########## DIFFERENT STUFF ##########
# R = const.physical_constants["molar gas constant"]      # Array of value, unit, error, siehe https://docs.scipy.org/doc/scipy-0.19.0/reference/constants.html
# search_replace_within_file('build/Tabelle_test.tex','find me','found you')    # Selbsterklärend


# convenience file writing for standard make Files
f = open('.pysuccess', 'w')
f.write('MarktkraM')
f.close()

os.chdir(old_path)
