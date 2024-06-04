import numpy as np
import matplotlib
matplotlib.use('tkagg')  # or 'tkagg'
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 26})
from numpy import log, exp

import pandas as pd
from scipy.integrate import RK45
from numpy import float128
import time

from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

data = pd.read_csv("results.csv", delimiter=',', header=0)
inputdata = pd.read_csv("kinetic_test_5.csv", delimiter=';', header=1)
inputdata = inputdata.apply(pd.to_numeric, errors='coerce')
inputdata = inputdata.dropna(axis=0, how='all')

errors = pd.read_csv("std.csv", delimiter=',', header=None)
rel_sigma = errors.to_numpy()
rel_sigma[np.isnan(rel_sigma)] = 1e-8
rel_sigma[rel_sigma == 0] = 1e-8

mixvals = inputdata.to_numpy(dtype=float, na_value=np.nan)


init_conc = []

ii = 0
while (mixvals[0,0] != mixvals[0,0]):
    init_conc = mixvals[0,:]
    mixvals = mixvals[1:,:]
    ii += 1
    

init_conc = init_conc[~np.isnan(init_conc)]
print(f"Removed {ii} rows of NaN values")



Rn = init_conc.size
print(f"Number of species: {Rn}")

mixvals = mixvals[:, 4:2*Rn+4]

ii = jj = 0

nan_array = np.sum(np.isnan(mixvals), axis=1)
mixvals = mixvals[nan_array == 0]

print(f"Removed {sum(nan_array != 0)} rows of NaN values")

print(f"Number of reactions: {sum(nan_array == 0)}")

reaction_factor = mixvals

xy = data.to_numpy()

i = xy[:,0]                     # index
t = xy[:,1]                     # time
a0 = xy[:,2]                    # reaction rate
N = np.array(xy[:,3])           # number of particles
R = np.array(xy[:,4:])          # species concentration

ReactantNames = list(data.columns.values)[4:]
assert len(ReactantNames) <= Rn, f"Number of species in the input file ({Rn}) does not match the number of species in the results file ({len(ReactantNames)})"


coeffs = pd.read_csv("coeffs.csv", delimiter=',', header=0)
names = coeffs.columns.values
coeffs = coeffs.to_numpy(dtype=float)[0,:]

T = coeffs[0]
E = coeffs[1]
Te = coeffs[2]
ne = coeffs[3]
V = coeffs[4]
P = coeffs[5]
coeffs = coeffs[6:-1]

assert len(coeffs) == reaction_factor.shape[0], f"Number of coefficients ({len(coeffs)}) does not match the number of reactions ({reaction_factor.shape[0]})"

print(f"Assuming that the system is at {T} K and {P} Pa")
R0 = 8.314 # J/(mol K)
avogadros_number = 6.02214076e23


nn = np.sum(init_conc)

print(f"V = {V}")

Rm = len(coeffs)

RungeKutta = True


def reaction_rate(Conc):
    factor = np.zeros_like(coeffs, dtype=np.float128)
    for i in range(Rm):
        factor[i] = coeffs[i] * np.prod(np.power(Conc, reaction_factor[i,:Rn]))
    return factor

def ode_O(t, y):
    dydt = np.zeros_like(y, dtype=np.float128)
    factor = reaction_rate(y)
    rf =  factor[:, np.newaxis] * reaction_factor
    for j in range(Rn):
        dydt[j] = np.sum(rf[:,j + Rn] - rf[:,j])
    return dydt

y0 = init_conc/V
t0 = 0

start_time = time.time()


rkfun = RK45(fun=ode_O, t0=0, y0=float128(y0), t_bound=0.1, rtol=1e-6)

t_values = [t0]
res_values = [y0]

i_max = 2_000_000

for i in range(i_max):
    # get solution step state
    rkfun.step()
    t_values.append(rkfun.t)
    res_values.append(rkfun.y)
    
    if i % 1000 == 0:
        print(f"RK45 at t = {t_values[-1]}" + 10*" ", end='\r')
    # break loop after modeling is finished
    if rkfun.status == 'finished':
        break

sum_res = np.sum(res_values, axis=1)
A0 = np.zeros_like(t_values, dtype=float128)
for i in range(len(t_values)):
    A0[i] = np.sum(reaction_rate(res_values[i]))*V

res_values = np.array(res_values, dtype=float128)/sum_res[:,np.newaxis]
t_values = np.array(t_values, dtype = float128)

end_time = time.time()

print(35*" ")
print(f"RK45 finished at t = {t_values[-1]}")
print(f"Duration of sim {end_time - start_time:0.3f}s)")
print(res_values[-1])

N = 1/N
norm_R = N[:,np.newaxis] * R    # Normalized species concentration

fig, ax = plt.subplots(1, 2, figsize=(32, 16))

ax[0].set_title("Species Concentration")
K = ax[0].loglog(t, norm_R, lw=3)
H = ax[0].loglog(t_values, res_values, lw=3, ls='--')
ax[0].set_xlabel("Time [s]")
ax[0].set_xlim([1e-4,1e0])
ax[0].set_ylim([1e-5,1.1e0])

ax[0].legend([(K[i], H[i]) for i in range(len(ReactantNames))], list(ReactantNames),
              handler_map={tuple: HandlerTuple(ndivide=None)})

for i, axes in enumerate(K):
    ax[0].fill_between(t, norm_R[:,i] - 2*norm_R[:,i]*rel_sigma[:,i+2], norm_R[:,i] +2*norm_R[:,i]*rel_sigma[:,i+2], alpha=0.3, color=axes.get_color())


ax[0].set_ylabel("Relative Concentration")
ax[0].grid()


ax[1].set_title("Reaction Rate")
G = ax[1].loglog(t, a0, 'k-', t_values, A0, 'm--', lw=3)
G[0].set_label(r"Reaction Rate KMC")
G[1].set_label(r"Reaction Rate RK45")
ax[1].fill_between(t, a0 - 2*a0*rel_sigma[:,0], a0 +2*a0*rel_sigma[:,0], alpha=0.3, color='k')
ax[1].set_xlabel("Time [s]")
ax[1].set_xlim([1e-4,1e0])
ax[1].set_ylim([1e6,1e8])
ax[1].legend()
ax[1].set_ylabel(r"Reaction Rate [s$^{-1}$]")
ax[1].grid()

plt.show(block=True)
