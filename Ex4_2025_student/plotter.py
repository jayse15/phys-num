
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from scipy.stats import linregress

plt.rcParams.update({
    'text.usetex': True,               # Use LaTeX for all text rendering
    'font.family': 'serif',            # Set font family to serif
    'font.serif': ['Computer Modern'], # Use Computer Modern
    'figure.dpi': 100,                 # DPI for displaying figures
    'font.size': 16,
    'lines.linewidth':2,
    'lines.markersize':5,
    'axes.formatter.limits':(-3, 3)
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

rho0 = parameters['rho0']
VR = parameters['VR']
r1 = parameters['r1']
R = parameters['R']

def phi_0(r):
    return (R**2-r**2)/4

def E_0(r):
    return r/2

N = np.array([5, 6, 8, 10, 15, 20, 30, 40, 50, 60, 70, 100])
outputs = []
datas=[]
nsimul = len(N)

for i in range(nsimul):
    output_file = f"data/N={N[i]}"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} N1={N[i]} N2={N[i]} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

E, D, Phi = [], [], []
conv_phi = []
for i in range(nsimul):  # Iterate through the results of all simulations
    e = np.loadtxt(outputs[i]+'_E.out')  # Load the output file of the i-th simulation
    d = np.loadtxt(outputs[i]+'_D.out')  # Load the output file of the i-th simulation
    p = np.loadtxt(outputs[i]+'_phi.out')  # Load the output file of the i-th simulation
    E.append(e)
    D.append(d)
    Phi.append(p)
    conv_phi.append(p[0, 1])

r = E[1][:, 0]
plt.figure()
plt.plot(r, E[1][:, 1], 'o',  c='orange', label='Data', markerfacecolor='none')
plt.plot(r, E_0(r), 'k--', label='True')
plt.xlabel(r'$r$[m]')
plt.ylabel(r'$E(r)$ [V/m]')
plt.legend()
plt.show()

plt.figure()
plt.plot(D[1][:, 0], D[1][:, 1], 'go', markerfacecolor='none')
plt.xlabel(r'$r$[m]')
plt.ylabel(r'$D(r)$ [C/m$^2$]')
plt.show()

r = Phi[1][:, 0]
plt.figure()
plt.plot(r, Phi[1][:, 1], 'ro', label='Data', markerfacecolor='none')
plt.plot(r, phi_0(r), 'k--', label='True')
plt.xlabel(r'$r$[m]')
plt.ylabel(r'$\phi(r)$ [V]')
plt.legend()
plt.show()

plt.figure()
plt.plot((1/N)**2, conv_phi, 'k+-')
plt.xlabel(r'$1/N^2_1(=1/N^2_2)$')
plt.ylabel(r'$\phi(0)$ [V]')
plt.show()
