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

nsteps = np.array([30e3, 50e3, 60e3, 70e3, 100e3, 200e3, 300e3])
epsilon = np.array([0.001, 0.01, 0.02, 0.05, 0.1, 0.5, 1])

# Simulations
output_e = []
output_n = []
convergence_list_x_n, convergence_list_y_n = [], []
convergence_list_x_e, convergence_list_y_e = [], []
nsimul = len(nsteps)  # Number of simulations to perform
esimul = len(epsilon)

for i in range(nsimul):
    output_file = f"data/nsteps={nsteps[i]}.out"
    output_n.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} nsteps={nsteps[i]} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

for i in range(esimul):
    output_file = f"data/tol={epsilon[i]}.out"
    output_e.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} tol={epsilon[i]} adapt=true output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

error_n=[]
datas_n=[]
for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(output_n[i])  # Load the output file of the i-th simulation
    t = data[:, 0]

    xx = data[-1, 3]
    yy = data[-1, 4]
    En = data[-1, 5]
    convergence_list_x_n.append(xx)
    convergence_list_y_n.append(yy)
    error_n.append(np.abs(En - data[0, 5])) # We use energy since it is the only known quantity at the end

    datas_n.append(data)

error_e=[]
datas_e=[]
jsteps=[]
for i in range(esimul):  # Iterate through the results of all simulations
    data = np.loadtxt(output_e[i])  # Load the output file of the i-th simulation
    t = data[:, 0]

    xx = data[-1, 3]
    yy = data[-1, 4]
    En = data[-1, 5]
    convergence_list_x_e.append(xx)
    convergence_list_y_e.append(yy)
    error_e.append(np.abs(En - data[0, 5])) # We use energy since it is the only known quantity at the end
    jsteps.append(len(data))
    datas_e.append(data)


lw = 1.5
fs = 16
if traj == True :
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()

    rmin=[]
    rmax=[]
    vmin=[]
    vmax=[]

    param_list=epsilon #change for tolerance
    datas=datas_e #change for tolerance
    p=r'\epsilon' #change for tolerance

    for i in range(nsimul):
        label = f'{param_list[i]:.2e}'
        label = label.replace(
            'e+0', r'\times 10^{').replace(
                'e+', r'\times 10^{').replace(
                    'e-0', r'\times 10^{-').replace(
                        'e-', r'\times 10^{-') + '}'
        ax1.plot(datas[i][:, 3], datas[i][:, 4], label = rf"${p}={label}$")
        ax2.plot(datas[i][:, 0], datas[i][:, 5], label = rf"${p}={label}$")

        v = np.sqrt(datas[i][:, 1]**2 + datas[i][:, 2]**2)
        r = np.sqrt(datas[i][:, 3]**2 + datas[i][:, 4]**2)

        rmin.append(r.min())
        rmax.append(r.max())
        vmin.append(v.min())
        vmax.append(v.max())

    ax1.set_aspect('equal', adjustable='box')
    ax1.legend()
    ax1.set_xlabel(r"$x$ [m]", fontsize=fs)
    ax1.set_ylabel(r"$y$ [m]", fontsize=fs)

    ax2.set_xlabel(r't [s]', fontsize=fs)
    ax2.set_ylabel(r'$E_{mec}$ [J]', fontsize=fs)
    ax2.legend()

    plt.show()

    plt.figure()
    plt.plot(param_list, rmin, 'k+-', label=r'$r_{\mathrm{min}}$')
    plt.xlabel(rf"${p}$", fontsize=fs)
    plt.ylabel(r"$r$ [m]", fontsize=fs)
    plt.hlines(y=10, xmin=param_list[0], xmax=param_list[-1], colors='r', label=r'$r_{\mathrm{min, true}}$')
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(param_list, rmax, 'k+-', label=r'$r_{\mathrm{max}}$')
    plt.xlabel(rf"${p}$", fontsize=fs)
    plt.ylabel(r"$r$ [m]", fontsize=fs)
    plt.hlines(y=10, xmin=param_list[0], xmax=param_list[-1], colors='r', label=r'$r_{\mathrm{max, true}}$')
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(param_list, vmin, 'k+-', label=r'$v_{\mathrm{min}}$')
    plt.xlabel(rf"${p}$", fontsize=fs)
    plt.ylabel(r"$v$ [m/s]", fontsize=fs)
    plt.hlines(y=10, xmin=param_list[0], xmax=param_list[-1], colors='r', label=r'$v_{\mathrm{min, true}}$')
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(param_list, vmax, 'k+-', label=r'$v_{\mathrm{max}}$')
    plt.xlabel(rf"${p}$", fontsize=fs)
    plt.ylabel(r"$v$ [m/s]", fontsize=fs)
    plt.hlines(y=10, xmin=param_list[0], xmax=param_list[-1], colors='r', label=r'$v_{\mathrm{max, true}}$')
    plt.legend()
    plt.show()



plt.figure()
plt.plot(nsteps, convergence_list_x_n, 'k+-', label='Schéma fixe')
plt.plot(jsteps, convergence_list_x_e, 'r+-', label='Schéma adaptif')
plt.xlabel(r"$N_{\mathrm{steps}}$", fontsize=fs)
plt.ylabel(r"$x_f$ [m/s]", fontsize=fs)
plt.legend()
plt.show()

plt.figure()
plt.plot(nsteps, convergence_list_y_n, 'k+-', label='Schéma fixe')
plt.plot(jsteps, convergence_list_y_e, 'r+-', label='Schéma adaptif')
plt.xlabel(r"$N_{\mathrm{steps}}$", fontsize=fs)
plt.ylabel(r"$y_f$ [m/s]", fontsize=fs)
plt.legend()
plt.show()

plt.figure()
plt.loglog(nsteps, error_n, 'k+-', label='Schéma fixe')
plt.loglog(jsteps, error_e, 'r+-', label='Schéma adaptif')
plt.xlabel(r"$N_{\mathrm{steps}}$", fontsize=fs)
plt.ylabel(r"$\Delta E_{\mathrm{mec}}$ [J/kg]", fontsize=fs)
plt.legend()
plt.show()

print(jsteps)
