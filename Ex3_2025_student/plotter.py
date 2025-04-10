
import numpy as np
import subprocess
import matplotlib.pyplot as plt
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

L0= 2*a*vy0
E0 = 0.5*(vx0**2 + vy0**2) - GM*mS/(2*a)
sqrt_term = np.sqrt((GM*mS)**2 + 2*E0*L0**2)
rmin_true = (-GM*mS + sqrt_term)/(2*E0)
rmax_true = (-GM*mS - sqrt_term)/(2*E0)
vmin_true = L0/rmax_true
vmax_true = L0/rmin_true

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

def to_latex_sci(num, precision=2):
    if num == 0:
        return r'$0$'

    sign = '-' if num < 0 else ''
    num = abs(num)

    exponent = int(np.floor(np.log10(num)))
    mantissa = num / 10**exponent

    if round(mantissa, precision)==1 :
        return rf"10^{{{exponent}}}"

    return rf'{sign}{mantissa:.{precision}f}\cdot 10^{{{exponent}}}'



traj = True # Set to true if we want to generate trajectories
adapt_traj = True # Set to true to have adaptive trajectories
inertial_traj = False # Set to true for traj in inertial frame for jup

nsteps = np.array([ 30e3, 40e3, 50e3, 60e3, 70e3, 80e3, 100e3, 200e3, 300e3])
epsilon = np.array([1e-3, 1e-2, 1e-1, 10, 100, 1e3, 1e4, 1e5])

# Simulations
output_e = []
output_n = []
convergence_list_x_n, convergence_list_y_n = [], []
convergence_list_x_e, convergence_list_y_e = [], []
nsimul = len(nsteps)  # Number of simulations to perform
esimul = len(epsilon)
param_list=[]

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

rmin_n=[]
rmax_n=[]
vmin_n=[]
vmax_n=[]
rmin_e=[]
rmax_e=[]
vmin_e=[]
vmax_e=[]

for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(output_n[i])  # Load the output file of the i-th simulation

    t_n = data[:, 0]
    xx = data[-1, 3]
    yy = data[-1, 4]
    En = data[-1, 5]
    convergence_list_x_n.append(xx)
    convergence_list_y_n.append(yy)
    error_n.append(np.abs(En - data[0, 5])) # We use energy since it is the only known quantity at the end
    v = np.sqrt(data[:, 1]**2 + data[:, 2]**2)
    r = np.sqrt(data[:, 3]**2 + data[:, 4]**2)

    rmin_n.append(r.min())
    rmax_n.append(r.max())
    vmin_n.append(v.min())
    vmax_n.append(v.max())


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
    v = np.sqrt(data[:, 1]**2 + data[:, 2]**2)
    r = np.sqrt(data[:, 3]**2 + data[:, 4]**2)

    rmin_e.append(r.min())
    rmax_e.append(r.max())
    vmin_e.append(v.min())
    vmax_e.append(v.max())

    jsteps.append(len(data))
    datas_e.append(data)


lw = 1.5
fs = 16
if traj == True :

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

    i=0
    if adapt_traj: i=0
    else : i=-1
    x, y, E, t = [datas[0][:, 3],
                  datas[-1][:,3]], [datas[0][:, 4],
                                    datas[-1][:,4]], [datas[0][:, 5],
                                                      datas[-1][:,5]], [datas[0][:, 0],
                                                                        datas[-1][:,0]]
    coord = ''
    if nsel_physics==2 : coord+="'"

    plt.figure()
    plt.plot(x[0], y[0], 'purple', label = rf"${p}={to_latex_sci(param_list[0], 1)}$")
    plt.plot(x[1], y[1], 'green', label = rf"${p}={to_latex_sci(param_list[-1], 1)}$")
    plt.plot(xS, 0, 'ro', ms=2, label='Soleil')

    if nsel_physics==2 :
        plt.plot(xJ, 0, 'o', c='orange', ms=1, label='Jupiter')

    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel(rf"$x{coord}$ [m]", fontsize=fs)
    plt.ylabel(rf"$y{coord}$ [m]", fontsize=fs)
    plt.legend()
    plt.show()

    if inertial_traj:
        x, y = Rprim_to_R(x[i], y[i], t[i])
        xS, yS = Rprim_to_R(xS, 0, t[i])
        xJ, yJ = Rprim_to_R(xJ, 0, t[i])

        plt.figure()
        plt.plot(xS, yS, 'r', label='Soleil')
        plt.plot(xJ, yJ, 'orange', label='Jupiter')
        plt.plot(x, y, 'purple', label = rf"${p}={to_latex_sci(param_list[i], 1)}$")
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel(rf"$x$ [m]", fontsize=fs)
        plt.ylabel(rf"$y$ [m]", fontsize=fs)
        plt.legend()
        plt.show()


    plt.figure()
    plt.plot(t[0]/s_per_year, E[0], 'purple', label = rf"${p}={to_latex_sci(param_list[0], 1)}$")
    plt.plot(t[1]/s_per_year, E[1], 'green', label = rf"${p}={to_latex_sci(param_list[-1], 1)}$")
    plt.xlabel(r't [an]', fontsize=fs)
    plt.ylabel(r'$E_{mec}$ [J]', fontsize=fs)
    plt.legend()
    plt.show()

    if nsel_physics==1:
        plt.figure()
        plt.loglog(nsteps, np.abs(rmin_n-rmin_true)/rmin_true, 'b+-', label=r'$\varepsilon_{r_{\mathrm{min}}}$')
        plt.loglog(nsteps, np.abs(rmax_n-rmax_true)/rmax_true, 'r+-', label=r'$\varepsilon_{r_{\mathrm{max}}}$')
        plt.loglog(nsteps, np.abs(vmin_n-vmin_true)/vmin_true, 'g+-', label=r'$\varepsilon_{v_{\mathrm{min}}}$')
        plt.loglog(nsteps, np.abs(vmax_n-vmax_true)/vmax_true, '+-', c='orange', label=r'$\varepsilon_{v_{\mathrm{max}}}$')
        plt.xlabel(r"$N_{\mathrm{steps}}$", fontsize=fs)
        plt.ylabel(r"$\varepsilon_{\mathrm{rel}}$", fontsize=fs)
        plt.legend()
        plt.show()

        plt.figure()
        plt.loglog(jsteps, np.abs(rmin_e-rmin_true)/rmin_true, 'b+-', label=r'$\varepsilon_{r_{\mathrm{min}}}$')
        plt.loglog(jsteps, np.abs(rmax_e-rmax_true)/rmax_true, 'r+-', label=r'$\varepsilon_{r_{\mathrm{max}}}$')
        plt.loglog(jsteps, np.abs(vmin_e-vmin_true)/vmin_true, 'g+-', label=r'$\varepsilon_{v_{\mathrm{min}}}$')
        plt.loglog(jsteps, np.abs(vmax_e-vmax_true)/vmax_true, '+-', c='orange', label=r'$\varepsilon_{v_{\mathrm{max}}}$')
        plt.xlabel(r"$\epsilon$", fontsize=fs)
        plt.ylabel(r"$\varepsilon_{\mathrm{rel}}$", fontsize=fs)
        plt.legend()
        plt.show()


jsteps = np.array(jsteps)

plt.figure()
plt.plot(1/jsteps**4, convergence_list_x_e, 'r+-')
plt.xlabel(r"$1/N_{\mathrm{steps}}^4$", fontsize=fs)
plt.ylabel(r"$x_f$ [m]", fontsize=fs)
plt.show()

plt.figure()
plt.plot(1/jsteps**4, convergence_list_y_e, 'r+-')
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

# Perform linear regression for convergence order
slope, intercept, r_value, p_value, std_err = linregress(np.log(nsteps), np.log(error_n))
y_fit = np.exp(intercept) * nsteps**slope
plt.loglog(nsteps, y_fit, c='b', ls='-', label=rf"$y = {to_latex_sci(np.exp(intercept),3)}/N^{{({-slope:.3f}\pm{std_err:.3f})}}_{{steps}}$", linewidth=lw)

slope, intercept, r_value, p_value, std_err = linregress(np.log(jsteps), np.log(error_e))
y_fit = np.exp(intercept) * jsteps**slope
plt.loglog(jsteps, y_fit, c='r', ls='-', label=rf"$y ={to_latex_sci(np.exp(intercept),3)}/N_{{steps}}^{{({-slope:.2f}\pm{std_err:.2f})}}$", linewidth=lw)

plt.loglog(nsteps, error_n, 'b+', label='Schéma fixe')
plt.loglog(jsteps, error_e, 'r+', label='Schéma adaptif')
plt.xlabel(r"$N_{\mathrm{steps}}$", fontsize=fs)
plt.ylabel(r"$\Delta E_{\mathrm{mec}}$ [J/kg]", fontsize=fs)
plt.legend(loc='upper left')
plt.show()

dt = np.diff(datas_e[-1][:, 0])
t_e = datas_e[-1][:, 0]
plt.figure()
plt.plot(t_e[:-1]/s_per_year, dt/s_per_year, 'g')
plt.xlabel(r'$t$ [an]', fontsize=fs)
plt.ylabel(r'$\Delta t$ [an]', fontsize=fs)
plt.show()