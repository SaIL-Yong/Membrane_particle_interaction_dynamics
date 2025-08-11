import os
import glob
import numpy as np
import igl
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator

# --------------------
# Publication‐style settings
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
    'xtick.labelsize':   12,
    'ytick.labelsize':   12,
    'legend.fontsize':   10,
    'lines.linewidth':   1,
})

colors    = ['#FFA500', 'green', 'blue']
mass_dirs = ["nv=252", "nv=702", "nv=2802"]
mass_vals = [252, 702, 2802]

def read_off(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    if lines[0] != 'OFF':
        raise ValueError("Missing OFF header")
    nv, nf = map(int, lines[1].split()[:2])
    V = np.array([list(map(float, l.split())) for l in lines[2:2+nv]])
    F = np.array([list(map(int,   l.split()[1:4])) for l in lines[2+nv:2+nv+nf]])
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
    rows   = []
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
        d    = os.path.join(os.getcwd(), dname)
        logf = os.path.join(d, "logfile.txt")
        offf = glob.glob(os.path.join(d, "hauser_cube.off"))
        if not (os.path.exists(logf) and offf):
            continue

        area     = mesh_area(offf[0])
        strength = read_strength(logf)
        header, data = read_log(logf)
        if header is None or data.size == 0:
            continue

        t_idx = header.index("Time")
        e_idx = header.index("AdhesionEnergy")
        t     = data[:, t_idx]
        E     = data[:, e_idx]

        wrap = -E / (strength * area)
        ax.plot(t, wrap, '-', color=color, label=f"{mval}")
        all_times.append(t)

    ax.set_xlabel("Time")
    ax.set_ylabel("Wrapping fraction")

    # force x‐ticks at 0, 400, 800, 1200, ...
    ax.xaxis.set_major_locator(MultipleLocator(400))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    # clip x‐axis to the next multiple of 400
    all_t = np.concatenate(all_times)
    max_t = np.ceil(all_t.max() / 400) * 400
    #ax.set_xlim(0, max_t)
    # after you’ve plotted everything, before saving:
    ax.set_xlim(0, 1200)    
    # y‐ticks every 0.2
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    # tick styling
    ax.tick_params(axis='both', direction='in', which='major', length=6)
    ax.tick_params(axis='both', direction='in', which='minor', length=4)

    ax.legend(frameon=False, loc='lower right')

    ax.text(
        0.56, 0.17,
        r"$N_p$",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="bottom",
        horizontalalignment="left"
    )

    plt.tight_layout(pad=0.0)
    plt.savefig("wrapping_vs_time_mesh_resolution.jpg", dpi=1200, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()
