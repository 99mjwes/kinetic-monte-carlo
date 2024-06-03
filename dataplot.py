import numpy as np
import matplotlib
matplotlib.use('tkagg')  # or 'tkagg'
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 26})
from numpy import log, exp
import os
import pandas as pd


script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

print(os.getcwd())
data = pd.read_csv("results.csv", delimiter=',', header=0)
errors = pd.read_csv("std.csv", delimiter=',', header=None)
xy = data.to_numpy()

rel_sigma = errors.to_numpy()
rel_sigma[np.isnan(rel_sigma)] = 1e-8
rel_sigma[rel_sigma == 0] = 1e-8


i = xy[:,0]                     # index
t = xy[:,1]                     # time
a0 = xy[:,2]                    # reaction rate
N = np.array(xy[:,3])           # number of particles
R = np.array(xy[:,4:])          # species concentration

ReactantNames = list(data.columns.values)

N = 1/N
norm_R = N[:,np.newaxis] * R    # Normalized species concentration

print(norm_R.shape)
print(rel_sigma[:,2:].shape)



fig, ax = plt.subplots(1, 2, figsize=(32, 16))

ax[0].set_title("Species Concentration")
K = ax[0].loglog(t, norm_R, lw=3)
ax[0].set_xlabel("Time [s]")
ax[0].set_xlim([1e-3,1e0])
ax[0].set_ylim([1e-5,1.1])
for i, axes in enumerate(K):
    axes.set_label(ReactantNames[i+4])
    ax[0].fill_between(t, norm_R[:,i] - 2*norm_R[:,i]*rel_sigma[:,i+2], norm_R[:,i] +2*norm_R[:,i]*rel_sigma[:,i+2], alpha=0.3, color=axes.get_color())
ax[0].set_ylabel("Relative Concentration")
ax[0].legend()
ax[0].grid()

ax[1].set_title("Reaction Rate")
ax[1].loglog(t, a0, 'k-', label=r"Reaction Rate", lw=3)
ax[1].fill_between(t, a0 - 2*a0*rel_sigma[:,0], a0 +2*a0*rel_sigma[:,0], alpha=0.3, color='k')
ax[1].set_xlabel("Time [s]")
ax[1].set_xlim([1e-3,1e0])
#  ax[1].set_ylim([1e6,1e7])
ax[1].legend()
ax[1].set_ylabel(r"Reaction Rate [s$^{-1}$]")
ax[1].grid()

plt.show(block=True)

