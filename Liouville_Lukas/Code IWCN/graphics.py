# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 13:25:33 2019

@author: lschulz
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 14:27:43 2018

@author: lschulz
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 13:55:59 2018

@author: lschulz
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import scipy.io 

def grayify_cmap(cmap):
    """Return a grayscale version of the colormap"""
    cmap = plt.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))
    
    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    colors[:, :3] = luminance[:, np.newaxis]
    
    return cmap.from_list(cmap.name + "_grayscale", colors, cmap.N)

def figsize(scale):
     fig_width_pt = 600/2                      # Get this from LaTeX using \the\textwidth
     inches_per_pt = 1.0/72.27                       # Convert pt to inch
     golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
     fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
     fig_height = fig_width*golden_mean*1.5              # height in inches
     fig_size = [fig_width,fig_height]
     return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
     "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
     "text.usetex": True,                # use LaTeX to write all text
     "font.family": "serif",
     "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
     "font.sans-serif": [],
     "font.monospace": [],
     "axes.labelsize": 11,               # LaTeX default is 10pt font.
     "font.size": 11,
     "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
     "xtick.labelsize": 11,
     "ytick.labelsize": 11,
     "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
     "pgf.preamble": [
         r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
         r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
         ]
     }

mpl.rcParams.update(pgf_with_latex)
plt.rcParams['figure.figsize']=figsize(1)


plt.close('all')

# all data inputs follow right now

x = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\x.mat')
x = x['x']

y = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\y.mat')
y = y['y']

k = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\k.mat')
k = k['k']

C = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\C50.mat')
C = C['C']

n10 = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\n10.mat')
n10 = n10['n']

n20 = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\n20.mat')
n20 = n20['n']

n30 = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\n30.mat')
n30 = n30['n']

n40 = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\n40.mat')
n40 = n40['n']

n50 = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\n50.mat')
n50 = n50['n']

n250 = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\n250.mat')
n250 = n250['n']

V_stat = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\V_stat.mat')
V_stat = V_stat['V_stat']

N = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\N.mat')
N = N['Nd_V']


plt.figure()
plt.pcolormesh(x[0,:]/1E-9, y[0,:]/1E-9,np.imag(C), cmap='seismic')
plt.xlabel('$\chi$ in nm')
plt.ylabel('$\\xi$ in nm')
plt.tight_layout()
plt.colorbar()
#plt.savefig('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\F_i50.png', dpi=600)



plt.figure()
plt.plot(x[0,:]/1E-9, np.real(n10[0,:]/1E+6), color='b', linewidth=2)
#plt.plot(x[0,:]/1E-9, np.real(n20[0,:]/1E+6))
#plt.plot(x[0,:]/1E-9, np.real(n30[0,:]/1E+6), color='r', linewidth=2, marker='d', markersize=6, markevery=(3,3))
#plt.plot(x[0,:]/1E-9, np.real(n40[0,:]/1E+6))
plt.plot(x[0,:]/1E-9, np.real(n50[0,:]/1E+6), color='k', linewidth=2, linestyle='-', marker='s', markersize=4, markevery=(1,3))

plt.plot(x[0,:]/1E-9, np.real(n250[0,:])/1E+6, color='m', linewidth=2, linestyle=':', marker='o', markersize=4, markevery=(2,3))
plt.xlabel('$\chi$ in nm')
plt.ylabel('carrier density $n(\chi)$ in $cm^{-3}$')
plt.legend(('$N=10$','$N=50$','$N=250$'))
plt.xlim([np.min(x)/1E-9, np.max(x)/1E-9])
plt.tight_layout()
plt.savefig('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\F2.png', dpi=600)

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.plot(x[0,:]/1E-9, V_stat[0,:], color='r', label='$V(\chi)$')
plt.xlim((np.min(x[0,:]/1E-9),np.max(x[0,:]/1E-9)))
plt.xlabel('$\chi$ in nm')
ax1.set_ylabel('potential $V(\chi)$ in eV')
ax1.set_ylim((-0.025,0.225))


ax2 = ax1.twinx()
ax2.plot(x[0,:]/1E-9, N[0,:]/1E+6, linestyle=':', color='b', label='$N(\chi)$')
ax2.set_ylabel('doping profil $N(\chi)$ in cm$^{-3}$')
fig.legend(loc=0, bbox_to_anchor=(1,0.95), bbox_transform=ax1.transAxes)
plt.tight_layout()

#plt.savefig('D:\\Lukas Schulz\\Quantum_Transport_Equations\\CONFERENCES\\ICWN-2019\\Code\\F3.png', dpi=600)

