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

N_excit = parameters['N_excit']
Nperiod = parameters['Nperiod']
Omega = parameters['Omega']
mu = parameters['mu']
m = parameters['m']
L = parameters['L']
B0 = parameters['B0']
theta0 = parameters['theta0']

nsteps_per = np.array([20, 30, 50, 75, 100, 150, 200, 300, 500, 1000, 2000])
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
for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(output[i])  # Load the output file of the i-th simulation
    t = data[:, 0]

    theta = data[-1, 1]  # final position, velocity, energy
    theta_dot = data[-1, 2]

    error = np.sqrt(om0**2*(theta - -theta_a(tFin))**2 + (theta_dot-thetadot_a(tFin)))
    errors[i] = error
    convergence_list.append(theta)

    datas.append(data)

    plt.figure()
    plt.plot(data[:,0], data[:,1]%(2*np.pi), 'r+-', linewidth=lw)
    plt.xlabel(r'$t$', fontsize=fs)
    plt.ylabel(r'$\theta$', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.grid(True)
    plt.show()
    
    if i == nsimul-1 : 
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


# For alpha = 0.5
norder = 2
C = 10e6
plt.figure()
plt.loglog(dt, errors, 'r+-', linewidth=lw)
#plt.loglog(dt, C*dt**2, 'g--', linewidth=lw, label=rf'$\sim \Delta t^2$')
plt.xlabel(r'$\Delta t$', fontsize=fs)
plt.ylabel(r'$\delta (t_f)$', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)
plt.show()
