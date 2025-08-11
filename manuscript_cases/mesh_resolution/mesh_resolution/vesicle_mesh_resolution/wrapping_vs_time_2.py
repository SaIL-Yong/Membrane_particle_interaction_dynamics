import os
import glob
import numpy as np
import igl
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# --------------------
# Publication-style settings
# --------------------
cm2in = lambda cm: cm / 2.54
fig_w_cm, fig_h_cm = 8.30, 6.0
fig_w_in, fig_h_in = cm2in(fig_w_cm), cm2in(fig_h_cm)

mpl.rcParams.update({
    'figure.figsize':    (fig_w_in, fig_h_in),
    'figure.dpi':        600,
    'font.family':       'Arial',
    'font.size':         12,
    'axes.labelsize':    12,
    'xtick.labelsize':   8,
    'ytick.labelsize':   12,
    'legend.fontsize':   10,
    'lines.linewidth':   1,
})

colors = ['#FFA500', 'green', 'blue']
mass_dirs = ["nv=2562_run2", "nv=10242"]
mass_vals = [2562, 10242]

def read_off(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    if lines[0] != 'OFF':
        raise ValueError("Missing OFF header")
    nv, nf = map(int, lines[1].split()[:2])
    V = np.array([list(map(float, l.split())) for l in lines[2:2+nv]])
    F = np.array([list(map(int, l.split()[1:4])) for l in lines[2+nv:2+nv+nf]])
    return V, F

def mesh_area(path):
    V, F = read_off(path)
    return np.sum(igl.doublearea(V, F)) / 2.0

def read_strength(path):
    with open(path) as f:
        for line in f:
            if "Particle adhesion strength:" in line:
                return float(line.split(":", 1)[1])
    raise ValueError("Strength line not found")

def read_log(path):
    header = None
    rows = []
    with open(path) as f:
        for L in f:
            L = L.strip()
            if L.startswith("Iteration"):
                header = L.split()
                continue
            if header:
                parts = L.split()
                if len(parts) == len(header):
                    rows.append(list(map(float, parts)))
    return header, np.array(rows)

def main():
    fig, ax = plt.subplots()
    all_times = []

    for color, dname, mval in zip(colors, mass_dirs, mass_vals):
        d = os.path.join(os.getcwd(), dname)
        logf = os.path.join(d, "logfile.txt")
        off_list = glob.glob(os.path.join(d, "hauser_cube.off"))
        if not (os.path.exists(logf) and off_list):
            continue

        area = mesh_area(off_list[0])
        strength = read_strength(logf)
        header, data = read_log(logf)
        if header is None or data.size == 0:
            continue

        t_idx = header.index("Time")
        e_idx = header.index("AdhesionEnergy")
        t = data[:, t_idx]
        E = data[:, e_idx]

        wrap = -E / (strength * area)
        ax.plot(t, wrap, '-', color=color, label=f"{mval}")
        all_times.append(t)

    ax.set_xlabel("Time")
    ax.set_ylabel("Wrapping fraction")

    # ----------- Dynamic tick spacing on time axis -------------
    all_t_concat = np.concatenate(all_times)
    x_range = all_t_concat.max() - all_t_concat.min()
    spacing = x_range / 6 if x_range > 0 else 1.0

    ax.xaxis.set_major_locator(MultipleLocator(spacing))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.tick_params(axis='both', direction='in', which='major', length=6)
    ax.tick_params(axis='both', direction='in', which='minor', length=4)

    ax.legend(frameon=False, loc='lower right')

    ax.text(0.56, 0.12, r"$N_v$", transform=ax.transAxes,
            fontsize=12, verticalalignment="bottom", horizontalalignment="left")

    plt.tight_layout(pad=0.0)
    plt.savefig("wrapping_vs_time_mesh_resolution.jpg", dpi=1200, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()
