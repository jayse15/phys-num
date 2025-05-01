
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


datas=[]
outputs=[]
CFL = parameters['CFL']
nx = parameters['nx']
ninit = parameters['n_init']
h0 = parameters['h00']
L = parameters['L']
g=9.81

om_n = np.sqrt(g*h0)*(ninit+0.5)*np.pi / L
oms = np.linspace(om_n-10, om_n+10, 40)
print(om_n)

states = ['static']
nsimul = len(states)
#nsimul = len(oms)
evolve = False
heat = True

E = []
for i in range(nsimul):
    output_file = f"data/{states[i]}.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} output={output_file} initial_state={states[i]}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

    f = np.loadtxt(outputs[i]+'_f')
    e = np.loadtxt(outputs[i]+'_en')
    x = np.loadtxt(outputs[i]+'_x')
    v = np.loadtxt(outputs[i]+'_v')
    t = f[:, 0]
    f = f[:, 1:]
    e = e[:, 1:]

    E.append(e.max())


    if evolve :
        plt.ion()
        fig, ax = plt.subplots()
        line, = ax.plot(x, f[0], 'b+-')
        ax.set_ylim(f.min() - 0.1, f.max() + 0.1)

        # Update loop
        for frame in range(1, len(f)):
            line.set_ydata(f[frame])
            ax.set_title(f"Wave at t = {t[frame]:.3f} s")
            plt.draw()
            plt.pause(0.01)

        plt.ioff()
        plt.show()

    if heat :
        # Plot heatmap
        plt.figure()
        extent = [x.min(), x.max(), t.min(), t.max()]
        max_abs = np.abs(f).max()

        plt.imshow(f, aspect='auto', extent=extent, origin='lower', cmap='seismic',
                   vmin=-max_abs, vmax=max_abs)
        plt.colorbar(label=r'Amplitude $f(x, t)$')
        plt.xlabel(r"$x$ [m]")
        plt.ylabel(r"$t$ [s]")
        plt.title(rf"$\beta_{{CFL}}={CFL}$, $n_x={int(nx)}$")
        plt.show()


plt.figure()
plt.plot(oms, E, 'r+-')
plt.xlabel(r"$\omega$ [rad/s]")
plt.ylabel(r"$\hat{E}$")
plt.show()
