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
     fig_height = fig_width*golden_mean*1.5             # height in inches
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

Lx = 66E-9
Nx = 132
x = np.linspace(-Lx/2, +Lx/2, Nx)

Ly = 130E-9
Ny = 131
y = np.linspace(-Ly/2, +Ly/2, Ny)
y = (y[0:Ny-1]+y[1:Ny])/2
Ny = Ny-1

x_part      = 0.33
m0          = 9.109E-31
q           = 1.602E-19
rel         = 0.57
h           = 6.626E-34
hb          = h/2/np.pi
kB          = 1.38E-23
T           = 300
me          = 0.063*m0

# all data inputs follow right now

af_fb = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\LVN_BASIS_CAP\\Code\\C_fb.mat')
af_fb = af_fb['C']

af_sc = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\LVN_BASIS_CAP\\Code\\C_sc.mat')
af_sc = af_sc['C']

V_fb = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\LVN_BASIS_CAP\\Code\\V_fb.mat')
V_fb = V_fb['V_tran']

V_sc = scipy.io.loadmat('D:\\Lukas Schulz\\Quantum_Transport_Equations\\LVN_BASIS_CAP\\Code\\V_sc.mat')
V_sc = V_sc['V_stat']

n_fb = (af_fb[int(Ny/2), :]+af_fb[int(Ny/2+1), :])/2
n_sc = (af_sc[int(Ny/2), :]+af_sc[int(Ny/2+1), :])/2

plt.figure()
plt.pcolormesh(x/1E-9, y/1E-9,np.real(af_fb), cmap='seismic')
plt.xlabel('$\chi$ in nm')
plt.ylabel('$\\xi$ in nm')
plt.tight_layout()
plt.colorbar()

#plt.savefig('D:\\Lukas Schulz\\Quantum_Transport_Equations\\LVN_BASIS_CAP\\Code\\dm_fb.png', dpi=1200)

plt.figure()
plt.pcolormesh(x/1E-9, y/1E-9,np.real(af_sc), cmap='seismic')
plt.xlabel('$\chi$ in nm')
plt.ylabel('$\\xi$ in nm')
plt.tight_layout()
plt.colorbar()

#plt.savefig('D:\\Lukas Schulz\\Quantum_Transport_Equations\\LVN_BASIS_CAP\\Code\\dm_sc.png', dpi=1200)


fig, ax1 = plt.subplots()
plt.plot(x/1E-9, V_fb[0,:], color='b', linewidth=2)
plt.plot(x/1E-9, V_sc[:,0], color='r', linewidth=2)

plt.ylabel('potential $V(\chi)$ in eV')
plt.xlabel('$\chi$ in nm')

ax2 = ax1.twinx()
plt.plot(x/1E-9, np.real(n_fb)/1E+6, color='b', linewidth=2)
plt.plot(x/1E-9, np.real(n_sc)/1E+6, color='r', linewidth=2)

plt.ylabel('carrier density $n(\chi)$ in $cm^{-3}$')
plt.xlim([np.min(x/1E-9), np.max(x/1E-9)])
plt.tight_layout()

plt.savefig('D:\\Lukas Schulz\\Quantum_Transport_Equations\\LVN_BASIS_CAP\\Code\\compo.png', dpi=1200)
















