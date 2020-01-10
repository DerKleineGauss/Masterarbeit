import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
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

annots = loadmat('data/test6Results_Lr66.mat')
x = annots['test6Results'][0][0][0]
x = x[0, :]
V = annots['test6Results'][0][0][1]
V = V[0, :]
I_over_x = annots['test6Results'][0][0][2]
I = annots['test6Results'][0][0][3]
I = I[0, :]
print(np.size(I))
plt.plot(V, I)
# plt.savefig('out.pdf')
print(V)
print(I)
plt.show()
