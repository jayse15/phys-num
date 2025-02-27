import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os

plt.rcParams.update({
    'text.usetex': True,               # Use LaTeX for all text rendering
    'font.family': 'serif',            # Set font family to serif
    'font.serif': ['Computer Modern'], # Use Computer Modern, the default LaTeX font
    'axes.labelsize': 16,              # Font size for axis labels
    'axes.titlesize': 18,              # Font size for titles
    'legend.fontsize': 14,             # Font size for legends
    'xtick.labelsize': 14,             # Font size for x-tick labels
    'ytick.labelsize': 14,             # Font size for y-tick labels
    'figure.dpi': 150,                 # DPI for displaying figures
})

# Parameters
executable = 'a.out'  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"./"
os.chdir(repertoire)

input_filename = 'configuration.in'  # Name of the input file


nsteps = np.array([4000, 40e3])#, 6000, 10000, 14e3, 20e3])
nsimul = len(nsteps)  # Number of simulations to perform

tfin = 259200
dt = tfin / nsteps


paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

# Simulations
outputs = []  # List to store output file names
convergence_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

error = np.zeros(nsimul)

for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation
    t = data[:, 0]

    vx = data[-1, 1]  # final position, velocity, energy
    vy = data[-1, 2]
    xx = data[-1, 3]
    yy = data[-1, 4]
    En = data[-1, 5]
    convergence_list.append(xx)
    error[i] = np.abs(En - data[0, 5]) # We use energy since it is the only known quantity at the end

    lw = 1.5
    fs = 16

    fig, ax = plt.subplots(constrained_layout=True)
    ax.plot(data[:, 3], data[:, 4])
    ax.set_xlabel('x [m]', fontsize=fs)
    ax.set_ylabel('y [m]', fontsize=fs)
    plt.show()

    fig2, ax2 = plt.subplots(constrained_layout=True)
    ax2.plot(data[:, 0], data[:, 5])
    ax2.set_xlabel(r't [s]', fontsize=fs)
    ax2.set_ylabel(r'$E_{mec}$ [J]', fontsize=fs)
    plt.show()



# uncomment the following to debug
#import pdb
#pbd.set_trace()
plt.figure()
plt.loglog(dt, error, 'r+-', linewidth=lw)
plt.xlabel(r'$\Delta$ t [s]', fontsize=fs)
plt.ylabel('final energy error [J]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)

"""
Si on n'a pas la solution analytique: on représente la quantite voulue
(ci-dessous v_y, modifier selon vos besoins)
en fonction de (Delta t)^norder, ou norder est un entier.
"""
norder = 2  # Modify if needed

plt.figure()
plt.plot(dt**norder, convergence_list, 'k+-', linewidth=lw)
plt.xlabel(fr'$\Delta$ $t^{norder}$ [s]', fontsize=fs)
plt.ylabel('v_y [m/s]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)

plt.show()
