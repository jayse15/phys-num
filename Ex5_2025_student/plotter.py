
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from scipy.stats import linregress

plt.rcParams.update({
    'text.usetex': True,               # Use LaTeX for all text rendering
    'font.family': 'serif',            # Set font family to serif
    'font.serif': ['Computer Modern'], # Use Computer Modern
    'figure.dpi': 300,                 # DPI for displaying figures
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

datas=[]
outputs=[]
CFL = parameters['CFL']
nx = parameters['nx']
ninit = parameters['n_init']
h0 = parameters['h00']
L = parameters['L']
impose_n = parameters['impose_nsteps']
f_hat = parameters['f_hat']
hL = parameters['hL']
hR = parameters['hR']
xa = parameters['xa']
xb = parameters['xb']
x1 = parameters['x1']
x2 = parameters['x2']

g=9.81
dx = L/nx

def f_a(x):
    """Fonction mode propre analytique"""
    return np.cos(np.pi*(ninit+0.5)/L * x)

def err(f, x, dx):
    """Erreur entre f et f_analytique"""
    a = np.abs(f-f_a(x))
    return ((a[1:] + a[:-1])/2 * dx).sum()

def v2(x):
    """h(x) non-uniform."""
    x = np.asarray(x)
    h_values = np.where(
        x <= xa,
        hL,
        np.where(
            x >= xb,
            hR,
            0.5 * (hL + hR) + 0.5 * (hL - hR) * np.cos(np.pi * (x - xa) / (xb - xa))
        )
    )
    return g*h_values

def WKB_A(x):
    """WKB for case A"""
    A0 = f_hat/(v2(x2/2 + x1/2))**(0.25)
    return A0*(v2(x))**(0.25)

def WKB_B(x):
    """WKB for case B"""
    A0 = f_hat*(v2(x2/2 + x1/2))**(0.25)
    return A0/(v2(x))**(0.25)

def WKB_C(x):
    """WKB for case C"""
    A0 = f_hat*(v2(x2/2 + x1/2))**(0.75)
    return A0/(v2(x))**(0.75)


om_n = np.sqrt(g*h0)*(ninit+0.5)*np.pi / L # Mode propre
nTn = 30 # Nombre de périodes de transit
tfin = nTn*2*np.pi/om_n # Periode
oms = [om_n - 3 + i * 6 / (100 - 1) for i in range(100)] # Omegas pour résonnance

states = ['right']
n = 30
nsteps=np.array([n]) #, 2*n, 3*n, 4*n, 8*n, 12*n, 16*n, 25*n, 32*n, 40*n, 50*n, 64*n])
Nx = np.array([nx]) #, 2*nx, 3*nx, 4*nx, 8*nx, 12*nx, 16*nx, 25*nx, 32*nx, 40*n, 50*nx, 64*nx])
#nsimul = len(states)
#nsimul = len(nsteps)
#nsimul = len(oms)
evolve = False # évolution continue de la vague
heat = False # Heatmap de l'amplitude, x et t
mode = False # Mode propres
conv = False
E_ = False
tsunami = True

if impose_n :
    dt=tfin/nsteps
    dx=L/Nx
    CFL = np.sqrt(g*h0)*dt/dx


E = []
error=[]
for i in range(1):
    output_file = f"data/xa={450}km.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} output={output_file} initial_state={states[0]}"
    if mode : cmd+=f" tfin={tfin} nsteps={nsteps[0]} nx={Nx[0]} om={oms[i]}"
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


    if conv:
        error.append(err(f[-1], x, dx[i]))

    E.append(e.max())


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

        if mode :
            x_plot=np.linspace(x.min(), x.max(), 300)
            plt.figure()
            plt.plot(x, f[-1], 'b+-', label=r'$f_{num}(x, t=T_n)$')
            plt.plot(x_plot, f_a(x_plot), 'r', label=r'$f_{ana}(x, t=T_n)$')
            plt.xlabel(r'$x$ [m]')
            plt.ylabel(r'$y$ [u.a.]')
            plt.title(rf"$\beta_{{CFL}}={CFL[i]:.3f}$, $n_x={int(Nx[i])}$, $n={int(ninit)}$")
            plt.legend(loc='upper right')
            plt.grid(alpha=0.8)
            plt.show()

    if tsunami:
        max_i = f.argmax(axis=1)
        maxes = []
        x_maxes = []

        num_cols = f.shape[1]

        for j, i_max in enumerate(max_i):
            if i_max <= 0:
                cols = [0, 1, 2]
            elif i_max >= num_cols - 1:
                cols = [num_cols - 3, num_cols - 2, num_cols - 1]
            else:
                cols = [i_max - 1, i_max, i_max + 1]

            maxes.append(f[j, cols])
            x_maxes.append(x[cols])

        X = np.array(x_maxes)
        Y = np.array(maxes)
        f_max = []
        x_max = []
        for k in range(len(X)):
            a, b, c = np.polyfit(X[k], Y[k], deg=2)
            crete = -b / (2 * a)
            x_max.append(crete)
            f_max.append(a * crete ** 2 + b * crete + c)

        x_max = np.array(x_max)
        f_max = np.array(f_max)



        k=1
        u_num = (x_max[2*k:] - x_max[:-2*k]) / (t[2*k:] - t[:-2*k])
        x_pos = x_max[k:-k]

        # Plotting
        plt.figure()
        plt.plot(x_max/1000, f_max, 'b', label=r'Solution numérique')
        plt.plot(x_max/1000, WKB_B(x_max), 'g--', label=r'Solution WKB')
        plt.xlabel(r"$x$ [km]")
        plt.ylabel(r"$f_{\mathrm{max}}(x,t)$ [m]")
        plt.grid(alpha=0.8)
        plt.title(rf"$\beta_{{CFL}}={CFL}$, $n_x={int(nx)}$")
        plt.legend()
        plt.show()

        plt.figure()
        plt.plot(x_pos/1000, u_num, 'r', label=r'Solution numérique')
        plt.plot(x_max/1000, np.sqrt(v2(x_max)), 'g--', label=r'Solution analytique')
        plt.xlabel(r"$x$ [km]")
        plt.ylabel(r"$u(x,t)$ [m/s]")
        plt.grid(alpha=0.8)
        plt.title(rf"$\beta_{{CFL}}={CFL}$, $n_x={int(nx)}$")
        plt.legend()
        plt.show()





    if heat :
        # Plot heatmap
        plt.figure()
        extent = [x.min()/1000, x.max()/1000, t.min()/3600, t.max()/3600]
        max_abs = np.abs(f).max()

        plt.imshow(f, aspect='auto', extent=extent, origin='lower', cmap='seismic',
                   vmin=-max_abs, vmax=max_abs)
        plt.colorbar(label=r'Amplitude $f(x, t)$ [m]')
        plt.xlabel(r"$x$ [km]")
        plt.ylabel(r"$t$ [h]")
        plt.title(rf"$\beta_{{CFL}}={CFL}$, $n_x={int(nx)}$")
        plt.show()

if conv:
    # Perform linear regression for convergence order
    slope, intercept, r_value, p_value, std_err = linregress(np.log(dt), np.log(error))
    y_fit = np.exp(intercept) * dt**slope

    plt.figure()
    plt.loglog(dt, error, 'rx')
    plt.loglog(dt, y_fit, 'k--', label=rf"$y = {np.exp(intercept):.3f}\Delta t^{{({slope:.3f}\pm{std_err:.3f})}}$")
    plt.xlabel(r"$\Delta t$ [s]")
    plt.ylabel(r"Erreur")
    plt.legend()
    plt.title(rf"$\beta_{{CFL}}={CFL[0]:.3f}$, $n=3$")
    plt.show()

if E_:
    E=np.array(E)
    plt.figure()
    plt.plot(oms, E, 'b+-', lw=1.5)
    plt.vlines(om_n, ymin=E.min(), ymax=E.max(), colors='r', lw=1, label=r'$\omega_n$')
    plt.xlabel(r"$\omega$ [rad/s]")
    plt.ylabel(r"$\hat{E}$ [u.a.]")
    plt.legend()
    plt.title(rf"$\beta_{{CFL}}={CFL}$, $n=3$, $n_x={int(nx)}$, $t_{{\mathrm{{fin}}}}={int(nTn)}T_n$")
    plt.show()
