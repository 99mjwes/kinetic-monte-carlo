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

xy = data.to_numpy()

i = xy[:,0]                     # index
t = xy[:,1]                     # time
a0 = xy[:,2]                    # reaction rate
a0_error = xy[:,3]              # reaction rate error
N = np.array(xy[:,4])           # number of particles
N_error = xy[:,5]               # number of particles error 
Raw = np.array(xy[:,6:])        # species concentration

print(xy.shape)

R = Raw[:, ::2]
Error = Raw[:, 1::2]

rel_sigma = np.hstack((a0_error[:,np.newaxis], N_error[:,np.newaxis], Error))

ReactantNames = list(data.columns.values)


N = 1/N
norm_R = N[:,np.newaxis] * R    # Normalized species concentration

print(norm_R.shape)
print(rel_sigma[:,2:].shape)



fig, ax = plt.subplots(1, 2, figsize=(32, 16))

ax[0].set_title("Species Concentration")
K = ax[0].loglog(t, norm_R, lw=3)
ax[0].set_xlabel("Time [s]")
ax[0].set_xlim([t[0], t[-1]])
# ax[0].set_ylim([1e-5,1.1])
for i, axes in enumerate(K):
    axes.set_label(ReactantNames[2*i+6])
    ax[0].fill_between(t, norm_R[:,i] - 2*norm_R[:,i]*rel_sigma[:,i+2], norm_R[:,i] +2*norm_R[:,i]*rel_sigma[:,i+2], alpha=0.3, color=axes.get_color())
ax[0].set_ylabel("Relative Concentration")
ax[0].legend()
ax[0].grid()

ax[1].set_title("Reaction Rate")
ax[1].loglog(t, a0, 'k-', label=r"Reaction Rate", lw=3)
ax[1].fill_between(t, a0 - 2*a0*rel_sigma[:,0], a0 +2*a0*rel_sigma[:,0], alpha=0.3, color='k')
ax[1].set_xlabel("Time [s]")
ax[1].set_xlim([t[0], t[-1]])
#  ax[1].set_ylim([1e6,1e7])
ax[1].legend()
ax[1].set_ylabel(r"Reaction Rate [s$^{-1}$]")
ax[1].grid()

plt.show(block=True)

