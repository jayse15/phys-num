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
    'figure.dpi': 300,                 # DPI for displaying figures
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

nsteps_per = np.array([100, 150, 200, 300, 500, 1000, 1500, 2000, 2500, 3000])
dt   = 2*np.pi/(Omega*nsteps_per)
nsimul = len(nsteps_per)
om0 = np.sqrt(mu*B0/(m*L**2/12))

if N_excit>0 :
      tFin = N_excit*2*np.pi/Omega
else:
      tFin = Nperiod*2*np.pi/Omega

dt   = 2*np.pi/(Omega*nsteps_per)

def theta_a(t):
    return theta0*np.cos(om0*t)

def thetadot_a(t):
    return -om0*theta0*np.sin(om0*t)


# Simulations
output = []
for i in range(nsimul):
    output_file = f"{'nsteps'}={nsteps_per[i]}.out"
    output.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} {'nsteps'}={nsteps_per[i]} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

lw = 1.5
fs = 16

errors = np.zeros(nsimul)
convergence_list=[]
datas = []
traj=False # Set to true to see trajectories and Emec
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
        plt.plot(data[:,0], data[:,1], 'r+-', linewidth=lw)
        plt.xlabel(r'$t$', fontsize=fs)
        plt.ylabel(r'$\theta$', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(True)
        plt.show()

        # plot energy
        plt.figure()
        plt.plot(data[:,0], data[:,3], 'r+-', linewidth=lw)
        plt.xlabel(r'$t$', fontsize=fs)
        plt.ylabel(r'$E_{mec}$', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(True)
        plt.show()

    if i == nsimul-1 :
        # plot Poincarre section
        times = np.arange(0, len(data), 20)
        poincare = data[times, 1:3]
        plt.figure()
        plt.plot(poincare[:,0], poincare[:,1], linewidth=lw)
        plt.xlabel(r'$\theta$', fontsize=fs)
        plt.ylabel(r'$\dot{\theta}$', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.grid(True)
        plt.show()


# plot final error
plt.figure()
plt.loglog(dt, errors, 'r+', linewidth=lw, ms=10)

# Perform linear regression for convergence order
slope, intercept, r_value, p_value, std_err = linregress(np.log(dt), np.log(errors))

y_fit = np.exp(intercept) * dt**slope
plt.loglog(dt, y_fit, c='black', ls='-', label=rf"$y \sim \Delta t^{{{slope:.2f}}}$", linewidth=lw)

plt.xlabel(r'$\Delta t$', fontsize=fs)
plt.ylabel(r'$\delta (t_f)$', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)
plt.legend()
plt.show()
