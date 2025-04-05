
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


nsteps = np.array([ 30e3, 40e3, 50e3, 60e3, 70e3, 80e3, 100e3, 200e3, 300e3])
epsilon = np.array([1e-3, 1e-2, 1e-1, 10, 100, 1e3, 1e4, 1e5])
outputs = []
datas=[]
nsimul = len(nsteps)

for i in range(nsimul):
    output_file = f"data/nsteps={nsteps[i]}.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} nsteps={nsteps[i]} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')


for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation

    t_n = data[:, 0]
    xx = data[-1, 3]
    yy = data[-1, 4]
    En = data[-1, 5]
    v = np.sqrt(data[:, 1]**2 + data[:, 2]**2)
    r = np.sqrt(data[:, 3]**2 + data[:, 4]**2)


    datas.append(data)
