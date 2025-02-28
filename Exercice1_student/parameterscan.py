import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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


nsteps = np.array([4000, 40e3])
nsimul = len(nsteps)  # Number of simulations to perform

tfin = 259200 # Check tfin is the same in configuration.in !!
dt = tfin / nsteps
mt     = 5.972e24
ml     = 7.348e22
dist   = 385000000
xt = -dist*ml/(mt+ml)
xl = dist*mt/(mt+ml)
rl     = 1737100
rt     = 6378100


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
    convergence_list.append(vy)
    error[i] = np.abs(En - data[0, 5]) # We use energy since it is the only known quantity at the end

    lw = 1.5
    fs = 16

    fig, ax = plt.subplots(constrained_layout=True)
    ax.plot(data[:, 3], data[:, 4])
    ax.set_xlabel('x [m]', fontsize=fs)
    ax.set_ylabel('y [m]', fontsize=fs)
    moon = patches.Circle((xl, 0), rl, fill=True, color='gray', lw=2, label='Lune')
    #earth = patches.Circle((xt, 0), rt, fill=True, color='blue', lw=2, label='Earth')
    ax.add_patch(moon)
    #ax.add_patch(earth)
    ax.set_aspect('equal', adjustable='box')
    ax.legend()
    plt.show()

    fig2, ax2 = plt.subplots(constrained_layout=True)
    ax2.plot(data[:, 0], data[:, 5])
    ax2.set_xlabel(r't [s]', fontsize=fs)
    ax2.set_ylabel(r'$E_{mec}$ [J]', fontsize=fs)
    plt.show()



# uncomment the following to debug
#import pdb
#pbd.set_trace()

"""
On représente l'erreur sur l'energie en fonction de N_steps et
v_y en fonction de (Delta t)^norder, ou norder est un entier.
"""
norder = 1  # Modify if needed

plt.figure()
plt.loglog(nsteps, error, 'r+-', linewidth=lw)
plt.loglog(nsteps, nsteps**-norder, linewidth=lw, ls='--', c='green', label=rf'~$1/N_{{steps}}^{norder}$')
plt.xlabel(r'$N_{steps}$', fontsize=fs)
plt.ylabel(r'$|E_{mec}(t=0)-E_{mec}(t_f)|$ [J]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.legend()
plt.grid(True)

plt.show()
plt.figure()
plt.plot(dt**norder, convergence_list, 'k+-', linewidth=lw)
plt.xlabel(fr'$\Delta$ $t^{norder}$ [s]', fontsize=fs)
plt.ylabel(r'$v_y$ [m/s]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)

plt.show()
