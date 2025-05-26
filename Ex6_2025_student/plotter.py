import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from scipy.stats import linregress

plt.rcParams.update({
    'text.usetex': True,               # Use LaTeX for all text rendering
    'font.family': 'serif',            # Set font family to serif
    'font.serif': ['Computer Modern'], # Use Computer Modern
    'figure.dpi': 275,                 # DPI for displaying figures
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

evolve = False
datas=[]
outputs=[]
for i in range(1):
    output_file = f"data/simul"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')



    if evolve :
        plt.ion()
        fig, ax = plt.subplots()
        line, = ax.plot(x, f[0], 'b')
        ax.set_ylim(f.min() - 0.1, f.max() + 0.1)

        # Update loop
        for frame in range(1, len(f)):
            line.set_ydata(f[frame])
            ax.set_title(f"Wave at t = {t[frame]:.3f} s")
            plt.draw()
            plt.pause(0.01)

        plt.ioff()
        plt.show()
