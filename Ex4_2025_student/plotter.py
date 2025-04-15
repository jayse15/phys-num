
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from scipy.stats import linregress

plt.rcParams.update({
    'text.usetex': True,               # Use LaTeX for all text rendering
    'font.family': 'serif',            # Set font family to serif
    'font.serif': ['Computer Modern'], # Use Computer Modern
    'figure.dpi': 250,                 # DPI for displaying figures
    'font.size': 16,
    'lines.linewidth':2,
    'lines.markersize':7,
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
    """Turns number into Latex string in scientific notation."""
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
e0=8.85418781e-12
triv = False

def phi_0(r):
    return (R**2-r**2)/4

def E_0(r):
    return r/2

def rho_true(r):
    return np.where(r < r1, rho0 * np.sin(np.pi * r / r1), 0)

N2 = np.array([10])
A = np.array([2]) # proportionnality constants between N1 and N2
nsimul = len(N2)
datas=[]
outputs=[]
plt.figure()
for i in range(len(A)):
    for j in range(nsimul):
        output_file = f"data/N2={N2[j]}"
        outputs.append(output_file)
        cmd = f"{repertoire}{executable} {input_filename} N1={A[i]*N2[j]} N2={N2[j]} output={output_file}"
        print(cmd)
        subprocess.run(cmd, shell=True)
        print('Done.')

    E, D, Phi, rho = [], [], [], []
    conv_phi = []
    err_phi = []
    traj=False
    for k in range(nsimul):  # Iterate through the results of all simulations
        e = np.loadtxt(outputs[k]+'_E.out')
        d = np.loadtxt(outputs[k]+'_D.out')
        p = np.loadtxt(outputs[k]+'_phi.out')
        r = np.loadtxt(outputs[k]+'_rho.out')
        E.append(e)
        D.append(d)
        Phi.append(p)
        rho.append(r)
        if triv:
            err = abs(p[0,1]-phi_0(0)) #error on phi(0) for trivial case
            err_phi.append(err)
        else:
            pr1 = p[np.abs(p[:, 0] - r1).argmin()][1] # phi(r1) found as closest phi to r1
            conv_phi.append(pr1)

        if traj:
            # Plot E, D and Phi
            r_plot = np.linspace(0, R, 200)
            plt.figure()
            plt.plot(e[:,0], e[:, 1], '-',  c='orange')
            plt.plot(r_plot, E_0(r_plot), 'k--', label='True')
            plt.xlabel(r'$r$[m]')
            plt.ylabel(r'$E(r)$ [V/m]')
            plt.grid(alpha=0.8)
            plt.show()

            plt.figure()
            plt.plot(d[:, 0], d[:, 1], 'g-')
            plt.xlabel(r'$r$[m]')
            plt.ylabel(r'$D(r)$ [C/m$^2$]')
            plt.grid(alpha=0.8)
            plt.show()

            plt.figure()
            plt.plot(p[:,0], p[:, 1], 'ro', label='Data')
            plt.plot(r_plot, phi_0(r_plot), 'k--', label=r'$\phi_{th}$')
            plt.xlabel(r'$r$[m]')
            plt.ylabel(r'$\phi(r)$ [V]')
            plt.grid(alpha=0.8)
            plt.legend()
            plt.show()
    a=A[i]
    if a>1: a=int(a)
    elif a==1 : a=''
    plt.plot(1/N2**2, conv_phi, '+-', label=rf'$N_1={a}N_2$') # Plots convergence for ratio A[i]
plt.xlabel(r'$1/N_2^2$')
plt.ylabel(r'$\phi(r_1)$ [V]')
plt.grid(alpha=0.8)
plt.legend()
plt.show()

if triv:
    # Perform linear regression for convergence order
    slope, intercept, r_value, p_value, std_err = linregress(np.log(N2), np.log(err_phi))
    y_fit = np.exp(intercept) * N2**slope

    # Log-Log plot for convergence in trivial case
    plt.figure()
    plt.loglog(N2, y_fit, 'k--', label=rf"$y = {to_latex_sci(np.exp(intercept),3)}/N^{{({-slope:.3f}\pm{std_err:.3f})}}$")
    plt.loglog(N2, err_phi, 'rx')
    plt.xlabel(r'$N_1(=N_2)$')
    plt.ylabel(r'$\Delta\phi(r=0)$ [V]')
    plt.grid(alpha=0.8)
    plt.legend()
    plt.show()
else :
    # Plot of rho calculated with D vs rho_true
    plt.plot(rho[0][:, 0], rho[0][:, 1], 'o', lw=3, c='purple', label=r'$\nabla\cdot D/\epsilon_0$')
    plt.plot(rho[0][:, 0], rho_true(rho[0][:, 0]), 'g--', label=r'$\rho_{\mathrm{lib}}/\epsilon_0$')
    plt.xlabel(r'$r$ [m]')
    plt.ylabel(r'$\rho/\epsilon_0$ [V/m$^2$]')
    plt.grid(alpha=0.8)
    plt.legend()
    plt.show()

    # Q calculated with gauss (L=1)
    Q_lib = 2*np.pi*D[0][:, 0]*D[0][:, 1]
    Q_tot = 2*np.pi*e0*E[0][:, 0]*E[0][:, 1]

    # Q calculated with integral of rho over volume and middle point approximation
    f = rho[0][:, 0] * e0*rho[0][:, 1]
    dr = rho[0][1:, 0]- rho[0][:-1, 0]
    partial = 0.5 * (f[:-1] + f[1:]) * dr
    integral = np.cumsum(partial)
    Q_lib_rho = 2.0 * np.pi * integral

    # Different Q plotted together
    plt.plot(D[0][:, 0], Q_lib, lw=3, label=r'$Q_{\mathrm{lib}}$')
    plt.plot(E[0][:, 0], Q_tot, label=r'$Q_{\mathrm{tot}}$')
    plt.plot(rho[0][1:, 0], Q_lib_rho, '--k', label=r'$Q_{\mathrm{lib}}$ avec $\rho_{\mathrm{lib}}$')
    plt.plot(D[0][:, 0], Q_tot-Q_lib, label=r'$Q_{\mathrm{pol}}$')
    plt.xlabel(r'$r$ [m]')
    plt.ylabel(r'$Q$ [C]')
    plt.legend()
    plt.grid(alpha=0.8)
    plt.show()
