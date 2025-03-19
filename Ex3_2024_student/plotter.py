import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pdb
import os

plt.rcParams.update({
    'text.usetex': True,               # Use LaTeX for all text rendering
    'font.family': 'serif',            # Set font family to serif
    'font.serif': ['Computer Modern'], # Use Computer Modern
    'figure.dpi': 150,                 # DPI for displaying figures
})

# Parameters
executable = 'a.out'  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"./"
os.chdir(repertoire)

input_filename = 'configuration.in'  # Name of the input file

parameters = {}
with open(input_filename, "r") as file:
    for line in file:
        # Ignore comments (marked by "//")
        line = line.split("//")[0].strip()
        if not line:
            continue  # Skip empty lines

        # Split by "=" to get key-value pairs
        key, value = line.split("=")
        key = key.strip()
        value = value.strip()

        # Convert to appropriate type
        if value.replace('.', '', 1).replace('e', '', 1).replace('-', '', 1).isdigit():  # Float check
            parameters[key] = float(value)
        else:  # Assume string (for file names)
            parameters[key] = value

tfin = parameters['tFin']
mS     = parameters['mS']
mJ     = parameters['mJ']
rS     = parameters['rS']
rJ     = parameters['rJ']
a     = parameters['a']
nsel_physics = parameters['nsel_physics']

if(nsel_physics==1):
    xS = 0
else:
    xS = -a*mJ/(mS+mJ)
    xJ = a*mS/(mS+mJ)

traj = True # Set to true if we want to generate trajectories

nsteps = [1e3, 2e3, 5e3, 10e3, 20e3, 30e3, 50e3, 60e3, 70e3, 100e3, 150e3, 200e3]
epsilon = [10e3, 8e3, 5e3, 4e3, 3e3, 2e3, 1e3, 800, 600, 500, 300, 200]

# Simulations
convergence_lists_x, convergence_lists_y = [], []
errors = []
nsteps_list, dts = [], []
for alpha in param_2:
    output = []  # List to store output file names
    convergence_list_x, convergence_list_y = [], []
    if traj == True:
        nsteps = np.array([4000, 40e3])
    else :
        nsteps = np.array([30e3, 40e3, 50e3, 75e3, 100e3, 150e3, 200e3]) # implicit and explicit too unstable for small N

    nsimul = len(nsteps)  # Number of simulations to perform
    dt = tfin / nsteps
    param_1 = nsteps
    nsteps_list.append(nsteps)
    dts.append(dt)


    for i in range(nsimul):
        output_file = f"{paramstr_1}={param_1[i]}{paramstr_2}={alpha}.out"
        output.append(output_file)
        cmd = f"{repertoire}{executable} {input_filename} {paramstr_1}={param_1[i]:.15g} {paramstr_2}={alpha:.15g} output={output_file}"
        print(cmd)
        subprocess.run(cmd, shell=True)
        print('Done.')

    error = np.zeros(nsimul)
    datas = []
    for i in range(nsimul):  # Iterate through the results of all simulations
        data = np.loadtxt(output[i])  # Load the output file of the i-th simulation
        t = data[:, 0]

        vx = data[-1, 1]  # final position, velocity, energy
        vy = data[-1, 2]
        xx = data[-1, 3]
        yy = data[-1, 4]
        En = data[-1, 5]
        convergence_list_x.append(xx)
        convergence_list_y.append(yy)
        error[i] = np.abs(En - data[0, 5]) # We use energy since it is the only known quantity at the end

        datas.append(data)

    convergence_lists_x.append(convergence_list_x)
    convergence_lists_y.append(convergence_list_y)
    errors.append(error)

    lw = 1.5
    fs = 16
    if traj == True :
        fig, ax = plt.subplots(constrained_layout=True)
        ax.plot(datas[0][:, 3], datas[0][:, 4], label = r"$N_{steps}=4000$", c='orange')
        ax.plot(datas[1][:, 3], datas[1][:, 4], label = r"$N_{steps}=40000$", c='purple')
        ax.set_xlabel(r"$x'$ [m]", fontsize=fs)
        ax.set_ylabel(r"$y'$ [m]", fontsize=fs)
        moon = patches.Circle((xl, 0), rl, fill=True, color='gray', lw=2, label='Lune')
        #earth = patches.Circle((xt, 0), rt, fill=True, color='blue', lw=2, label='Earth')
        ax.add_patch(moon)
        #ax.add_patch(earth)
        ax.set_aspect('equal', adjustable='box')
        ax.legend()
        plt.show()

        fig2, ax2 = plt.subplots(constrained_layout=True)
        ax2.plot(datas[0][:, 0], datas[0][:, 5], label = r"$N_{steps}=4000$", c='orange')
        ax2.plot(datas[1][:, 0], datas[1][:, 5], label = r"$N_{steps}=40000$", c='purple')
        ax2.set_xlabel(r't [s]', fontsize=fs)
        ax2.set_ylabel(r'$E_{mec}$ [J]', fontsize=fs)
        ax2.legend()
        plt.show()
