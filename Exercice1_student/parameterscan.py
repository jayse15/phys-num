import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os

# Parameters
# TODO adapt to what you need (folder path executable input filename)
executable = 'exe.exe'  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"C:/Users/Administrator/physnum/2025/EX1/SOLUTION/"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Name of the input file


nsteps = np.array([4000, 6000, 10000, 14e3, 20e3]) # TODO change
nsimul = len(nsteps)  # Number of simulations to perform

tfin = 259200  # TODO: Verify that the value of tfin is EXACTLY the same as in the input file
with open(input_filename, 'r') as file:
    lines = file.readlines()
    for line in lines:
        if 'tfin' in line:
            # Assuming the line is like "tfin = 259200"
            try:
                line_content = line.split('!')[0]
                file_tfin = int(line_content.split('=')[1].strip())
                if file_tfin != tfin:
                    print(f"Warning: tfin in the file ({file_tfin}) does not match the expected value ({tfin}).")
            except ValueError:
                print("Error: Unable to parse tfin from the input file.")
            break

dt = tfin / nsteps


paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

# Simulations
outputs = []  # List to store output file names
convergence_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    cmd = f"{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

error = np.zeros(nsimul)

for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation
    t = data[:, 0]

    vx = data[-1, 1]  # final position, velocity, energy
    vy = data[-1, 2]
    xx = data[-1, 3]
    yy = data[-1, 4]
    En = data[-1, 5]
    convergence_list.append(xx)
    # TODO compute the error for each simulation
    error[i] =  np.abs(En - data[0, 5])

lw = 1.5
fs = 16

fig, ax = plt.subplots(constrained_layout=True)
ax.plot(data[:, 1], data[:, 2])
ax.set_xlabel('x [m]', fontsize=fs)
ax.set_ylabel('y [m]', fontsize=fs)



# uncomment the following if you want debug
#import pdb
#pbd.set_trace()
plt.figure()
plt.loglog(dt, error, 'r+-', linewidth=lw)
plt.xlabel('\Delta t [s]', fontsize=fs)
plt.ylabel('final position error [m]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)

"""
Si on n'a pas la solution analytique: on représente la quantite voulue
(ci-dessous v_y, modifier selon vos besoins)
en fonction de (Delta t)^norder, ou norder est un entier.
"""
norder = 1  # Modify if needed

plt.figure()
plt.plot(dt**norder, convergence_list, 'k+-', linewidth=lw)
plt.xlabel('\Delta t [s]', fontsize=fs)
plt.ylabel('v_y [m/s]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)

plt.show()
