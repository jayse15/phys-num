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
vx0 = parameters['vx0']
vy0 = parameters['vy0']

s_per_year = 3.1536e7
GM=6.674e-11

L0= abs(2*a*vy0)
E0 = 0.5*(vx0**2 + vy0**2) + GM*mS/(2*a)
sqrt_term = np.sqrt(GM**2 * mS**2 + 2*E0*L0**2)
rmin_true = (GM*mS - sqrt_term)/(2*E0)
rmax_true = (GM*mS + sqrt_term)/(2*E0)
vmin_true = L0/rmax_true
vmax_true = L0/vmin_true

if(nsel_physics==1):
    xS = 0
else:
    xS = -a*mJ/(mS+mJ)
    xJ = a*mS/(mS+mJ)

def Rprim_to_R(x, y, t):
    Om = np.sqrt(GM*mS/(pow(a, 2)*xJ))
    x_new = x*np.cos(Om*t) - y*np.sin(Om*t)
    y_new = x*np.sin(Om*t) + y*np.cos(Om*t)
    return x_new, y_new

traj = True # Set to true if we want to generate trajectories
adapt_traj = True # Set to true to have adaptive trajectories
inertial_traj = True # Set to true for traj in inertial frame for jup

nsteps = np.array([20e3, 30e3, 50e3, 60e3, 80e3, 100e3, 150e3, 200e3, 250e3, 300e3])
epsilon = np.geomspace(1e-3, 1e6, num=10)

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
t_n=[]
for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(output_n[i])  # Load the output file of the i-th simulation
    t_n = data[:, 0]

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

    rmin=[]
    rmax=[]
    vmin=[]
    vmax=[]

    if adapt_traj:
        param_list=epsilon
        datas=datas_e
        p=r'\epsilon'
        n=esimul
    else :
        param_list=nsteps
        datas=datas_n
        p=r'N_{\mathrm{steps}}'
        n=nsimul

    for i in range(n):

        v = np.sqrt(datas[i][:, 1]**2 + datas[i][:, 2]**2)
        r = np.sqrt(datas[i][:, 3]**2 + datas[i][:, 4]**2)

        rmin.append(r.min())
        rmax.append(r.max())
        vmin.append(v.min())
        vmax.append(v.max())

    i=0
    if adapt_traj: i=0
    else : i=-1
    x, y, E, t = datas[i][:, 3], datas[i][:, 4], datas[i][:, 5], datas[i][:, 0]
    coord = ''
    if nsel_physics==2 : coord+="'"

    plt.figure()
    plt.plot(x, y, 'purple', label = rf"${p}={param_list[i]}$")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel(rf"$x{coord}$ [m]", fontsize=fs)
    plt.ylabel(rf"$y{coord}$ [m]", fontsize=fs)
    plt.legend()
    plt.show()

    if inertial_traj:
        x, y = Rprim_to_R(x, y, t)
        xS, yS = Rprim_to_R(xS, 0, t)
        xJ, yJ = Rprim_to_R(xJ, 0, t)

        plt.figure()
        plt.plot(xS, yS, 'r', label='Soleil')
        plt.plot(xJ, yJ, 'orange', label='Jupiter')
        plt.plot(x, y, 'purple', label = rf"${p}={param_list[i]}$")
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel(rf"$x$ [m]", fontsize=fs)
        plt.ylabel(rf"$y$ [m]", fontsize=fs)
        plt.legend()
        plt.show()


    plt.figure()
    plt.plot(t/s_per_year, E, label = rf"${p}={param_list[i]}$")
    plt.xlabel(r't [an]', fontsize=fs)
    plt.ylabel(r'$E_{mec}$ [J]', fontsize=fs)
    plt.legend()
    plt.show()

    if nsel_physics==1:
        plt.figure()
        plt.plot(param_list, rmin, 'k+-', label=r'$r_{\mathrm{min}}$')
        plt.xlabel(rf"${p}$", fontsize=fs)
        plt.ylabel(r"$r$ [m]", fontsize=fs)
        plt.hlines(y=rmin_true, xmin=param_list[0], xmax=param_list[-1], colors='r', label=r'$r_{\mathrm{min, true}}$')
        plt.legend()
        plt.show()

        plt.figure()
        plt.plot(param_list, rmax, 'k+-', label=r'$r_{\mathrm{max}}$')
        plt.xlabel(rf"${p}$", fontsize=fs)
        plt.ylabel(r"$r$ [m]", fontsize=fs)
        plt.hlines(y=rmax_true, xmin=param_list[0], xmax=param_list[-1], colors='r', label=r'$r_{\mathrm{max, true}}$')
        plt.legend()
        plt.show()

        plt.figure()
        plt.plot(param_list, vmin, 'k+-', label=r'$v_{\mathrm{min}}$')
        plt.xlabel(rf"${p}$", fontsize=fs)
        plt.ylabel(r"$v$ [m/s]", fontsize=fs)
        plt.hlines(y=vmin_true, xmin=param_list[0], xmax=param_list[-1], colors='r', label=r'$v_{\mathrm{min, true}}$')
        plt.legend()
        plt.show()

        plt.figure()
        plt.plot(param_list, vmax, 'k+-', label=r'$v_{\mathrm{max}}$')
        plt.xlabel(rf"${p}$", fontsize=fs)
        plt.ylabel(r"$v$ [m/s]", fontsize=fs)
        plt.hlines(y=vmax_true, xmin=param_list[0], xmax=param_list[-1], colors='r', label=r'$v_{\mathrm{max, true}}$')
        plt.legend()
        plt.show()



plt.figure()
plt.plot(1/np.array(jsteps)**4, convergence_list_x_e, 'r+-')
plt.xlabel(r"$1/N_{\mathrm{steps}}^4$", fontsize=fs)
plt.ylabel(r"$x_f$ [m]", fontsize=fs)
plt.show()

plt.figure()
plt.plot(1/np.array(jsteps)**4, convergence_list_y_e, 'r+-')
plt.xlabel(r"$1/N_{\mathrm{steps}}^4$", fontsize=fs)
plt.ylabel(r"$y_f$ [m]", fontsize=fs)
plt.show()

plt.figure()
plt.plot(1/nsteps**4, convergence_list_x_n, 'b+-')
plt.xlabel(r"$1/N_{\mathrm{steps}}^4$", fontsize=fs)
plt.ylabel(r"$x_f$ [m]", fontsize=fs)
plt.show()

plt.figure()
plt.plot(1/nsteps**4, convergence_list_y_n, 'b+-')
plt.xlabel(r"$1/N_{\mathrm{steps}}^4$", fontsize=fs)
plt.ylabel(r"$y_f$ [m]", fontsize=fs)
plt.show()

plt.figure()
plt.loglog(nsteps, error_n, 'b+-', label='Schéma fixe')
plt.loglog(jsteps, error_e, 'r+-', label='Schéma adaptif')
plt.xlabel(r"$N_{\mathrm{steps}}$", fontsize=fs)
plt.ylabel(r"$\Delta E_{\mathrm{mec}}$ [J/kg]", fontsize=fs)
plt.legend()
plt.show()

dt = np.diff(datas_e[-1][:, 0])
t_e = datas_e[-1][:, 0]
plt.figure()
plt.plot(t_e[:-1]/s_per_year, dt/s_per_year, 'r')
plt.xlabel(r'$t$ [an]', fontsize=fs)
plt.ylabel(r'$\Delta t$ [an]', fontsize=fs)
plt.show()
