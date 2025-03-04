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

tfin = 259200 # Check tfin is the same in configuration.in !!
mt     = 5.972e24
ml     = 7.348e22
dist   = 385000000
xt = -dist*ml/(mt+ml)
xl = dist*mt/(mt+ml)
rl     = 1737100
rt     = 6378100
traj = False # Set to true if we want to generate trajectories

paramstr_1 = 'nsteps'  # Parameter name to scan

paramstr_2 = 'alpha'
param_2 = [0, 0.5, 1] # Parameter values to scan

# Simulations
convergence_lists_x, convergence_lists_y = [], []
errors = []
nsteps_list, dts = [], []
for alpha in param_2:
    output = []  # List to store output file names
    convergence_list_x, convergence_list_y = [], []
    if traj == True:
        nsteps = np.array([4000, 40e3])
    else :
        nsteps = np.array([30e3, 40e3, 50e3, 75e3, 100e3, 150e3, 200e3]) # implicit and explicit too unstable for small N

    nsimul = len(nsteps)  # Number of simulations to perform
    dt = tfin / nsteps
    param_1 = nsteps
    nsteps_list.append(nsteps)
    dts.append(dt)


    for i in range(nsimul):
        output_file = f"{paramstr_1}={param_1[i]}{paramstr_2}={alpha}.out"
        output.append(output_file)
        cmd = f"{repertoire}{executable} {input_filename} {paramstr_1}={param_1[i]:.15g} {paramstr_2}={alpha:.15g} output={output_file}"
        print(cmd)
        subprocess.run(cmd, shell=True)
        print('Done.')

    error = np.zeros(nsimul)
    datas = []
    for i in range(nsimul):  # Iterate through the results of all simulations
        data = np.loadtxt(output[i])  # Load the output file of the i-th simulation
        t = data[:, 0]

        vx = data[-1, 1]  # final position, velocity, energy
        vy = data[-1, 2]
        xx = data[-1, 3]
        yy = data[-1, 4]
        En = data[-1, 5]
        convergence_list_x.append(xx)
        convergence_list_y.append(yy)
        error[i] = np.abs(En - data[0, 5]) # We use energy since it is the only known quantity at the end

        datas.append(data)

    convergence_lists_x.append(convergence_list_x)
    convergence_lists_y.append(convergence_list_y)
    errors.append(error)

    lw = 1.5
    fs = 16
    if traj == True :
        fig, ax = plt.subplots(constrained_layout=True)
        ax.plot(datas[0][:, 3], datas[0][:, 4], label = r"$N_{steps}=4000$", c='orange')
        ax.plot(datas[1][:, 3], datas[1][:, 4], label = r"$N_{steps}=40000$", c='purple')
        ax.set_xlabel(r"$x'$ [m]", fontsize=fs)
        ax.set_ylabel(r"$y'$ [m]", fontsize=fs)
        moon = patches.Circle((xl, 0), rl, fill=True, color='gray', lw=2, label='Lune')
        #earth = patches.Circle((xt, 0), rt, fill=True, color='blue', lw=2, label='Earth')
        ax.add_patch(moon)
        #ax.add_patch(earth)
        ax.set_aspect('equal', adjustable='box')
        ax.legend()
        plt.show()

        fig2, ax2 = plt.subplots(constrained_layout=True)
        ax2.plot(datas[0][:, 0], datas[0][:, 5], label = r"$N_{steps}=4000$", c='orange')
        ax2.plot(datas[1][:, 0], datas[1][:, 5], label = r"$N_{steps}=40000$", c='purple')
        ax2.set_xlabel(r't [s]', fontsize=fs)
        ax2.set_ylabel(r'$E_{mec}$ [J]', fontsize=fs)
        ax2.legend()
        plt.show()



    # uncomment the following to debug
    #import pdb
    #pbd.set_trace()


# For alpha = 0 and 1
norder = 1
C = 10e8 # Constant for log log graph comparasion

plt.figure()
plt.loglog(nsteps_list[0], errors[0], 'r+-', linewidth=lw, label=r'$\alpha =0$')
plt.loglog(nsteps_list[2], errors[2], 'b+-', linewidth=lw, label=r'$\alpha =1$')
plt.loglog(nsteps_list[0], C*nsteps_list[0]**-norder, 'g--', linewidth=lw, label=r'$\sim 1/N_{steps}$')
plt.xlabel(r'$N_{steps}$', fontsize=fs)
plt.ylabel(r'$\Delta E_{mec}(t_f)$ [J]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.legend()
plt.grid(True)
plt.show()

plt.figure()
plt.plot(dts[0]**norder, convergence_lists_x[0], 'r+-', linewidth=lw, label=r'$\alpha =0$')
plt.plot(dts[2]**norder, convergence_lists_x[2], 'b+-', linewidth=lw, label=r'$\alpha =1$')
plt.xlabel(fr'$\Delta$ $t$ [$s$]', fontsize=fs)
plt.ylabel(r"$x'(t_f)$ [m]", fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.legend()
plt.grid(True)
plt.show()

plt.figure()
plt.plot(dts[0]**norder, convergence_lists_y[0], 'r+-', linewidth=lw, label=r'$\alpha =0$')
plt.plot(dts[2]**norder, convergence_lists_y[2], 'b+-', linewidth=lw, label=r'$\alpha =1$')
plt.xlabel(fr'$\Delta$ $t$ [$s$]', fontsize=fs)
plt.ylabel(r"$y'(t_f)$ [m]", fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.legend()
plt.grid(True)
plt.show()

# For alpha = 0.5
norder = 2
C = 10e6
plt.figure()
plt.loglog(nsteps_list[1], errors[1], 'r+-', linewidth=lw)
plt.loglog(nsteps_list[1], C*nsteps_list[1]**-norder, 'g--', linewidth=lw, label=rf'$\sim 1/N_{{steps}}^{norder}$')
plt.xlabel(r'$N_{steps}$', fontsize=fs)
plt.ylabel(r'$\Delta E_{mec}(t_f)$ [J]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.legend()
plt.grid(True)
plt.show()

plt.figure()
plt.plot(dts[1]**norder, convergence_lists_x[1], 'k+-', linewidth=lw)
plt.xlabel(fr'$\Delta$ $t^{norder}$ [$s^{norder}$]', fontsize=fs)
plt.ylabel(r"$x'(t_f)$ [m]", fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)
plt.show()

plt.figure()
plt.plot(dts[1]**norder, convergence_lists_y[1], 'k+-', linewidth=lw)
plt.xlabel(fr'$\Delta$ $t^{norder}$ [$s^{norder}$]', fontsize=fs)
plt.ylabel(r"$y'(t_f)$ [m]", fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)
plt.show()
