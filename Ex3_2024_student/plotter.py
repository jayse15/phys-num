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

nsteps = np.array([20e3, 30e3, 50e3, 60e3, 70e3, 100e3])
epsilon = np.array([10e3, 8e3, 5e3, 4e3, 3e3, 2e3, 1e3, 800, 600, 500, 300, 200])
param_list = nsteps
param = 'nsteps'

# Simulations
output = []  # List to store output file names
convergence_list_x, convergence_list_y = [], []
nsimul = len(param_list)  # Number of simulations to perform
dt = tfin / nsteps

for i in range(nsimul):
    output_file = f"data/{param}={param_list[i]}.out"
    output.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} {param}={param_list[i]} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

error=[]
datas=[]
for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(output[i])  # Load the output file of the i-th simulation
    t = data[:, 0]

    xx = data[-1, 3]
    yy = data[-1, 4]
    En = data[-1, 5]
    convergence_list_x.append(xx)
    convergence_list_y.append(yy)
    error.append(np.abs(En - data[0, 5])) # We use energy since it is the only known quantity at the end

    datas.append(data)

lw = 1.5
fs = 16
if traj == True :
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    for i in range(nsimul):
        ax1.plot(datas[i][:, 3], datas[i][:, 4], label = rf"${param}={param_list[i]}$")
        ax2.plot(datas[i][:, 0], datas[i][:, 5], label = f"{param}={param_list[i]}")

    sun = patches.Circle((xS, 0), rS, fill=True, color='orange', lw=2, label='Soleil')
    #earth = patches.Circle((xt, 0), rt, fill=True, color='blue', lw=2, label='Earth')
    ax1.add_patch(sun)
    #ax.add_patch(earth)
    ax1.set_aspect('equal', adjustable='box')
    ax1.legend()
    ax1.set_xlabel(r"$x'$ [m]", fontsize=fs)
    ax1.set_ylabel(r"$y'$ [m]", fontsize=fs)

    ax2.set_xlabel(r't [s]', fontsize=fs)
    ax2.set_ylabel(r'$E_{mec}$ [J]', fontsize=fs)
    ax2.legend()

    plt.show()
