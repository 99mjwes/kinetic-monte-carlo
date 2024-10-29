# Kinetic Monte-Carlo Code for Simulation of Atmospheric Chemistry
A Kinetic Monte Carlo Code for Simulation of Atmospheric Chemistry in an AMReX Framework

This Kinetic Monte Carlo (KMC) model can accurately
track the evolution of chemical reactions with an arbitrary number of reaction and species. It can
account for the pressure, temperature and electric field strength of the system. 
Further, it outputs the results in a format that can easily be interpreted and analyzed.

## Installation and Compilation 
This programme requires AMReX to work. Get it at https://amrex-codes.github.io/amrex/.

Download the source code into an appropriate folder with

    git clone https://github.com/99mjwes/kinetic-monte-carlo/

and modify the *GNUmakefile* to link to the AMReX folder. Then compile the programme with

    make -j

Create an *input* plain text file for the parameters. Write one parameter per line in the format of 
*parameter = value*

Then, run the programme with

    ./main1d.gnu.ex inputs

## Parameters

- **T** - The temperature of the cell in Kelvin
- **P** - The pressure of the cell in Pascal
- **E** - The electric field strenght in the cell in Townsend
- **filename** - Name of reaction input file. Defaults to *reactions.tsv*
- **savename** - Name of the output file. Defaults to *results.csv*
- **runtime** - Simulation runtime in seconds.
- **first_save** - First time point for which to save a data point.
- **n_saves** - Number of data points in the saved results.
- **max_steps** - Maximum number of reactions to perform in a simulation. Liberal estimate used if omitted.
- **n_iter** - Number of simulation iterations. A higher value reduces variance in the result, but increases computational time. Defaults to 1. 
- **ne** - Initial electron density in $m^{-3}$. Used to compute the initial electron quantity if the latter is not provided.
- **num_workers** - Number of workers/threads to use. Defaults to 1.
- **seed** - Initial seed for the KMC. Random number is used if omitted.
- **export_coeffs** - Exports coefficients for the *dataRK45.py* comparison tool. Defaults to *true*


## Data Input

See the included *simpletest.tsv* for an
example of how to structure the reactions.

Data inputs takes its form of a tsv-file with tabulator separated values. The first lines may contain
any text or data for descriptive purposes. The data part of the file must be a sheet of size 4 + 2N
by M + 2, where N is the number of species and M is the number of reactions.

The first valid line must have the character entry "A", while the next three
entries should be characters corresponding to the coefficients "B" and "C" and the type of equation
to use. The next N entries contains the names of the species. Names with a
leading "*" will not be saved in the output files, this can help reduce the clutter of the output if
many species are present. The last N entries should be empty. From the example file:

| | | | | | | | | | |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **A** | **B** | **C** | **Equation** | O | O2 | O3 | - | - | - |

For the next line the first four entries should be left blank. The next N entries will contain the
initial quantity of the species, while, again, the last N entries should be empty.

| | | | | | | | | | |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| - | - | - | - | 100000 | 10 | 0 | - | - | - |

For every following line the first three entries will specify the coefficients A, B and C and the fourth entry will provide the coefficient equation to use,
where the default will be used if the entry is left blank. The next N entries will specify the
reaction inputs and the last N values will be the reaction outputs. Outputs should be listed in the same order as inputs. If the parameter A is left blank,
the reaction will not be defined.

| | | | | | | | | | |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
8.31E-033 | -0.63 | 0 | default | 2 | 1 | 0 | 0 | 2 | 0
2E-11 | 0 |-2300 | default | 1 | 0 | 1 | 0 | 2 | 0
6.9E-034 | -1.25 | 0 | default | 1 | 2 | 0 | 0 | 1 | 1

### Notes
- If the programme detects two identical reactions it will throw a warning. You may
desire to have two similar reactions describing different mechanisms, and the programme will allow
that, but unintended reaction duplications can skew the results significantly.

- The custom equations follow the AMReX Parser format which is very
versatile in its capabilities.

## Data Output

Data output happens in two files; results.csv and coeffs.csv. These are indeed
comma separated and provides various ways of interpreting the results of a simulation.

results.csv stores the processed results from the simulation. It contains 4 + N columns with
4 + 2*n_saves rows. The first column specifies the saved step in integer intervals, the second is the
simulation time, the third and fourth is the reaction rate and its standard deviation and then the total reactant quantity with standard deviation in the fith and sixth and so on. The last four rows
specifies the start and end reactant quantities respectively. In simulations with multiple iterations,
the resulting quantity is computed as the average of each run.

### Notes

- For the purposes of  data visualization and cross validation *dataplot.py* and *dataRK45.py* are provided with the fluid model. 

- These are sensitive to the filenames of the results, and the
latter needs access to the input data as well. Modify them accordingly.