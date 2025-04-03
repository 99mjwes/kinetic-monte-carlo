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

# ---------------------------

# Put tsv filename here:
input_file = "test5.tsv"

# Put latex reference for the default equation here:
default_equation_latex_reference = "eq:ABC"

# ---------------------------

inputdata = pd.read_csv(input_file, delimiter='\t', header=1)

mixvals = inputdata.to_numpy(dtype=str)
mixvals = mixvals[mixvals[:,0] != 'nan']

Rn = (mixvals.shape[1]-4)//2
M = mixvals.shape[0]

print(f"Writing LaTeX table to 'table.tex' with {M} reactions and {Rn} species.")

names = inputdata.columns[4:Rn+4]

ii = 0

latex_table = []
latex_table.append(r'\begin{tabular}{|cccccc|}')
latex_table.append(r'\hline')
latex_table.append(r' &Reaction & A & B & C & Equation \\')
latex_table.append(r'\hline')


for i in range(M):
        # First column in the LaTeX table
        input = []
        output = []
        coeffs = []
        eq = ''
        
        for j in range(Rn):
            input_value = int(float(mixvals[i, j+4]))
            if (input_value != 0 & ~np.isnan(input_value)) :
                if input_value == 1:
                    input.append(f'{names[j]} ')
                else:
                    input.append(f'{input_value} {names[j]} ')
                
            output_value = int(float(mixvals[i, j+Rn+4]))
            if (output_value != 0 & ~np.isnan(output_value)) :
                if output_value == 1:
                    output.append(f'{names[j]} ')
                else:
                    output.append(f'{output_value} {names[j]} ')
        
        for j in range(3):
            
            val = mixvals[i, j]
            if val == 'nan':
                coeffs.append(' ')
                break
            elif int(float(val)) - float(val) == 0.0 :
                val = str(int(float(val)))

            coeffs.append('\\num{' + val + '}')
            
        if (mixvals[i, 3] == 'nan' or mixvals[i, 3] == 'default'):
            eq = '\\eqref{' + default_equation_latex_reference + '}'
        else:
            eq = mixvals[i, 3]
            
        latex_table.append('$R_{' + f'{i+1}' + '}$ & ' + ' \\ce{' + '+ '.join(input) + '-> ' + '+ '.join(output) + '} & ' + ' & '.join(coeffs) + ' & $' + eq + '$ \\\\')
        
latex_table.append(r'\hline')

# Add the LaTeX table footer
latex_table.append(r'\end{tabular}')

with open('table.tex', 'w') as f:
    for line in latex_table:
        f.write(line + '\n')
