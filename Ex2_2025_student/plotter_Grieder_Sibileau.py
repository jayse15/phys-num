import numpy as np
import subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pdb
import os
from scipy.stats import linregress

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

N_excit = parameters['N_excit']
Nperiod = parameters['Nperiod']
Omega = parameters['Omega']
mu = parameters['mu']
m = parameters['m']
L = parameters['L']
B0 = parameters['B0']
theta0 = parameters['theta0']

nsteps_per = np.array([50, 60, 70, 80, 100, 125, 150, 200, 250, 300])
om0 = np.sqrt(mu*B0/(m*L**2/12))
Omega = 2*om0
nsimul = len(nsteps_per)


if N_excit>0 :
      tFin = N_excit*2*np.pi/Omega
      dt   = 2*np.pi/(Omega*nsteps_per)
else:
      tFin = Nperiod*2*np.pi/om0
      dt   = 2*np.pi/(om0*nsteps_per)

def theta_a(t):
    return theta0*np.cos(om0*t)

def thetadot_a(t):
    return -om0*theta0*np.sin(om0*t)

import re
def modify_config(filename, variable, new_value):
    with open(filename, 'r') as file:
        lines = file.readlines()

    with open(filename, 'w') as file:
        for line in lines:
            if re.match(rf'^\s*{variable}\s*=', line):  # Handles spaces around '='
                file.write(f"{variable} = {new_value}\n")  # Keeps consistent formatting
            else:
                file.write(line)

# Set to true to see trajectories and Emec
traj=False 
# Set to true for corresponding question, e.g A=True for question a)
A = False
B = False
C = True
D = False
E = False

# Simulations
output = []
if C: 
    thetas = np.arange(0, 2*np.pi, np.pi/16)
    thetas_dot = np.arange(0, 1e3, 5e2)
    for O in thetas: 
        for O_dot in thetas_dot: 
            modify_config(input_filename, 'theta0', O)
            modify_config(input_filename, 'thetadot0', O_dot)
            output_file = f"{'nsteps'}={nsteps_per[-1]}.out"
            output.append(output_file)
            cmd = f"{repertoire}{executable} {input_filename} {'nsteps'}={nsteps_per[-1]} output={output_file}"
            print(cmd)
            subprocess.run(cmd, shell=True)
            print('Done.')
else : 
    for i in range(nsimul):
        output_file = f"{'nsteps'}={nsteps_per[i]}.out"
        output.append(output_file)
        cmd = f"{repertoire}{executable} {input_filename} {'nsteps'}={nsteps_per[i]} output={output_file}"
        print(cmd)
        subprocess.run(cmd, shell=True)
        print('Done.')

lw = 1.5
fs = 20

errors = np.zeros(nsimul)
convergence_list=[]
datas = []

for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(output[i])  # Load the output file of the i-th simulation
    t = data[:, 0]

    theta_f = data[-1, 1]  # final position, velocity, energy
    theta_dot_f = data[-1, 2]

    error = np.sqrt(om0**2*(theta_f-theta_a(tFin))**2 + (theta_dot_f-thetadot_a(tFin))**2)
    errors[i] = error
    convergence_list.append(theta_f)

    datas.append(data)

    if traj==True :
        # plot trajectories
        plt.figure()
        plt.plot(data[:,0], data[:,1], 'r-', linewidth=lw, ms=2)
        plt.xlabel(r'$t$ [s]', fontsize=fs)
        plt.ylabel(r'$\theta$', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(True)
        plt.show()

        # plot energy
        plt.figure()
        plt.plot(data[:,0], data[:,3], 'r-', linewidth=lw, ms=2)
        plt.xlabel(r'$t$ [s]', fontsize=fs)
        plt.ylabel(r'$E_{mec}$ [J]', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(True)
        plt.show()

    if i == nsimul-1 :
        if  B:
            # plot phase space
            phase = data[:, 1:3]
            plt.figure()
            plt.plot(phase[:,0], phase[:,1], 'o', ms=1)
            plt.xlabel(r'$\theta$', fontsize=fs)
            plt.ylabel(r'$\dot{\theta}$', fontsize=fs)
            plt.xticks(fontsize=fs)
            plt.yticks(fontsize=fs)
            plt.grid(True)
            plt.show()

            # plot dEmec vs Pnc
            Pnc = data[:, [0,4]]
            dEmec = data[:, [0,5]]
            plt.figure()
            plt.plot(Pnc[:,0], Pnc[:,1], linewidth=lw, c='purple', alpha=0.5)
            plt.xlabel(r'$t$  [s]', fontsize=fs)
            plt.ylabel(r'$P_{nc}(t)$ [J/s]', fontsize=fs)
            plt.xticks(fontsize=fs)
            plt.yticks(fontsize=fs)
            plt.grid(True)
            plt.show()

            plt.figure()
            plt.plot(dEmec[:,0], dEmec[:,1], linewidth=lw, c='orange', alpha=0.5)
            plt.xlabel(r'$t$  [s]', fontsize=fs)
            plt.ylabel(r'$\frac{dE_{mec}(t)}{dt}$ [J/s]', fontsize=fs)
            plt.xticks(fontsize=fs)
            plt.yticks(fontsize=fs)
            plt.grid(True)
            plt.show()

if C:
    # plot Poincarre section
    poincare_list = []
    for out in output:
        data = np.loadtxt(out) 
        t = data[:, 0]
        theta_f = data[-1, 1]  
        theta_dot_f = data[-1, 2]
        datas.append(data)
        times = np.arange(0, len(data), nsteps_per[-1])
        poincare = data[times, 1:3]
        poincare_list.append(poincare)
    poincare_array = np.vstack(poincare_list)  # Stack list into 2D NumPy array
    plt.figure()
    plt.plot(poincare_array[:, 0], poincare_array[:, 1], 'o', linewidth=lw)
    plt.xlabel(r'$\theta$', fontsize=fs)
    plt.ylabel(r'$\dot{\theta}$', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.grid(True)
    plt.show()


if A:
    # plot final error
    plt.figure()
    plt.loglog(dt, errors, 'r+', linewidth=lw, ms=10)

    # Perform linear regression for convergence order
    slope, intercept, r_value, p_value, std_err = linregress(np.log(dt), np.log(errors))

    y_fit = np.exp(intercept) * dt**slope
    plt.loglog(dt, y_fit, c='black', ls='-', label=rf"$y \sim \Delta t^{{{slope:.2f}}}$", linewidth=lw)

    plt.xlabel(r'$\Delta t$ [s]', fontsize=fs)
    plt.ylabel(r'$\delta (t_{\mathrm{fin}})$ [J]', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.grid(True)
    plt.legend()
    plt.show()

if B:
    # plot final positions
    plt.figure()
    plt.plot(dt**2, convergence_list, c='k', marker='+', markeredgecolor='red', linewidth=lw, ms=10)
    plt.xlabel(r'$\Delta t^2$ [s]', fontsize=fs)
    plt.ylabel(r'$\theta(t_{\mathrm{fin}})$ [J]', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.grid(True)
    plt.show()
