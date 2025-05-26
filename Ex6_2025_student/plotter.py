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
        # Ignore comments (marked by "%")
        line = line.split("%")[0].strip()
        if not line:
            continue  # Skip empty lines

        # Split by "=" to get key-value pairs
        key, value = line.split("=")
        key = key.strip()
        value = value.strip()

        # Convert to appropriate type
        if value.replace('.', '', 1).replace('e', '', 1).replace('-', '', 1).isdigit():  # Float check
            parameters[key] = float(value)
        elif (value=='false' or value=='true'):
            parameters[key] = (value=='true')
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

L = parameters['xR']-parameters['xL']
x0 = parameters['x0']
om0 = parameters['om0']
nx = parameters['Nintervals']
nsteps = parameters['Nsteps']
tfin = parameters['tfin']
k = 2*parameters['n']*np.pi/L
t_trans = 0.045

def x_class(t):
    return x0*np.cos(om0*t) + k/om0 * np.sin(om0*t)

def p_class(t):
    return -om0*x0*np.sin(om0*t) + k*np.cos(om0*t)

evolve = False
heat = True
classic = False
conv = False
nx_fixe = False
tunnel = False
tunnel2 = False
datas=[]
outputs=[]
conv_list = []
conv_values = []
nsimul=1
V0 = []

if conv:
    conv_values = 2*np.array([80, 140, 160, 200, 280, 350, 600, 800, 1000])
    nsimul=len(conv_values)

if tunnel:
    if tunnel2:
        V0 = np.array([500, 1300, 3000])
        heat = True
    else:
        V0 = np.geomspace(500, 15e3, 25)
    nsimul = len(V0)


for i in range(nsimul):
    output_file = f"data/simul"
    if conv:
        if nx_fixe:
            output_file+=f'_nsteps={conv_values[i]}'
        else:
            output_file+=f'_nx={conv_values[i]}'
    if tunnel :
        output_file+=f'_V0={V0[i]}'

    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} output={output_file}"
    if conv:
        if nx_fixe:
            cmd+=f' Nsteps={conv_values[i]}'
        else:
            cmd+=f' Nintervals={conv_values[i]}'

    if tunnel:
        cmd+=f' V0={V0[i]}'

    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

    psi = np.loadtxt(outputs[i]+'_psi2.out')
    pot = np.loadtxt(outputs[i]+'_pot.out')
    obs = np.loadtxt(outputs[i]+'_obs.out')

    psi = psi.reshape(psi.shape[0], -1, 3)
    abs = psi[:, :, 0]
    real = psi[:, :, 1]
    imag = psi[:, :, 2]
    x = pot[:, 0]
    t = obs[:, 0]
    p_neg = obs[:, 1]
    p_pos = obs[:, 2]
    E = obs[:, 3]
    x_moy = obs[:, 4]
    x2_moy = obs[:, 5]
    p_moy = obs[:, 6]
    p2_moy = obs[:, 7]
    delta_P = np.sqrt(p2_moy - p_moy**2)
    delta_x = np.sqrt(x2_moy - x_moy**2)

    if conv:
        conv_list.append(x_moy[-1])

    if tunnel:
        index = np.argmin(np.abs(t - t_trans))
        conv_list.append([E[index], p_pos[index]])



    if evolve :
        plt.ion()
        fig, ax = plt.subplots()
        line1, = ax.plot(x, abs[0], 'k', label=r'$|\psi(x, t)|$')
        line2, = ax.plot(x, real[0], 'r', label=r'$Re(\psi(x, t))$')
        line3, = ax.plot(x, imag[0], 'b', label=r'$Im(\psi(x, t))$')
        ax.legend()
        ax.set_ylim(real.min() - 0.1, abs.max() + 0.1)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$\psi(x, t)$')

        # Update loop
        for frame in range(1, len(abs)):
            line1.set_ydata(abs[frame])
            line2.set_ydata(real[frame])
            line3.set_ydata(imag[frame])
            ax.set_title(f"Wave at t = {t[frame]:.3f} s")
            plt.draw()
            plt.pause(0.0001)

        plt.ioff()
        plt.show()

    if heat :
        # Plot heatmap
        plt.figure()
        extent = [x.min(), x.max(), t.min(), t.max()]
        max_abs = abs.max()

        plt.imshow(abs, aspect='auto', extent=extent, origin='lower', cmap='seismic',
                   vmin=-max_abs, vmax=max_abs)
        plt.colorbar(label=r'$|\psi(x, t)|$')
        plt.xlabel(r"$x$")
        plt.ylabel(r"$t$ [s]")
        if tunnel2:
            signs = ['>', r'\approx', '<']
            plt.title(rf'$\langle E \rangle {signs[i]} V_0 = {V0[i]}$')
        plt.show()

        if not tunnel2:
            plt.figure()
            extent = [x.min(), x.max(), t.min(), t.max()]
            max_abs = np.abs(real).max()

            plt.imshow(real, aspect='auto', extent=extent, origin='lower', cmap='seismic',
                    vmin=-max_abs, vmax=max_abs)
            plt.colorbar(label=r'$Re(\psi(x, t))$')
            plt.xlabel(r"$x$")
            plt.ylabel(r"$t$ [s]")
            plt.show()


    if classic :
        plt.figure()
        plt.plot(t, x_moy, 'r', label=r'Quantique')
        plt.plot(t, x_class(t), 'g--', label=r'Classique')
        plt.xlabel(r"$t$ [s]")
        plt.ylabel(r"$\langle x \rangle$")
        plt.grid()
        plt.legend(loc='lower left')
        plt.show()

        plt.figure()
        plt.plot(t, p_moy, 'r', label=r'Quantique')
        plt.plot(t, p_class(t), 'g--', label=r'Classique')
        plt.xlabel(r"$t$ [s]")
        plt.ylabel(r"$\langle p \rangle$")
        plt.grid()
        plt.legend(loc='lower left')
        plt.show()

        plt.figure()
        plt.plot(t, p_neg+p_pos, 'b')
        plt.xlabel(r"$t$ [s]")
        plt.ylabel(r"$P_{\mathrm{tot}}$")
        plt.grid()
        plt.show()

        plt.figure()
        plt.plot(t, E, 'orange')
        plt.xlabel(r"$t$ [s]")
        plt.ylabel(r"$\langle E \rangle$")
        plt.grid()
        plt.show()

        plt.figure()
        plt.plot(t, delta_P*delta_x, 'g')
        plt.hlines(0.5, 0, tfin, colors='r', ls='--', label=r'$\frac{\hbar}{2}$')
        plt.xlabel(r"$t$ [s]")
        plt.ylabel(r"$\Delta x \cdot \Delta p$")
        plt.grid()
        plt.legend()
        plt.show()

    if tunnel2:
        plt.figure()
        plt.plot(t, p_pos, 'r', label=r'$P_{x>0}(t)$')
        plt.plot(t, p_neg, 'b', label=r'$P_{x<0}(t)$')
        plt.xlabel(r"$t$ [s]")
        plt.ylabel(r"$P$")
        signs = ['>', r'\approx', '<']
        plt.title(rf'$\langle E \rangle {signs[i]} V_0 = {V0[i]}$')
        plt.legend()
        plt.show()





if conv:
    plt.figure()
    value = 0
    if nx_fixe:
        plt.xlabel(r"$\Delta t^2$ [s]")
        value = (tfin/conv_values)
    else :
        plt.xlabel(r"$h_x^2$")
        value = (L/conv_values)

    plt.plot(value**2, conv_list, 'k+-')
    plt.ylabel(r"$\langle x \rangle (t_{\mathrm{fin}})$")
    plt.grid()
    plt.show()

if (tunnel and not tunnel2):
    values = np.array(conv_list)
    plt.figure()
    plt.plot(values[:, 0]/V0, values[:, 1], '+-', color='purple')
    plt.xlabel(r"$\langle E \rangle / V_0$")
    plt.ylabel(r"$P_{\mathrm{trans}}$")
    plt.grid()
    plt.title(r'$t_{\mathrm{trans}}=0.045$ s')
    plt.show()
