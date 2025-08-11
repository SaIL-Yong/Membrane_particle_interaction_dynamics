import os
import glob
import numpy as np
import igl
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter

# --------------------
# Publication-style settings
# --------------------
cm2in = lambda cm: cm / 2.54
fig_w_cm, fig_h_cm = 8.30, 6.0      # single-column width, shorter height
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

base_cmap = mpl.cm.plasma

# --------------------
# Directories & labels
# --------------------
mass_dirs = [
    "nv=252",
    "nv=702",
    "nv=2802",
    
]
mass_vals = [252,702, 2802]
#colors    = base_cmap(np.linspace(0.1, 0.9, len(mass_dirs)))
colors=['#FFA500','green','blue']
# --------------------
# File readers
# --------------------
def read_off(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    if lines[0] != 'OFF':
        raise ValueError("Missing OFF header")
    nv, nf = map(int, lines[1].split()[:2])
    V = np.array([list(map(float, l.split())) 
                  for l in lines[2:2+nv]])
    F = np.array([list(map(int, l.split()[1:4])) 
                  for l in lines[2+nv:2+nv+nf]])
    return V, F

def mesh_area(path):
    V, F = read_off(path)
    return np.sum(igl.doublearea(V, F)) / 2.0

def read_strength(path):
    with open(path) as f:
        for line in f:
            if "Particle adhesion strength:" in line:
                return float(line.split(":",1)[1])
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

# --------------------
# Main plotting
# --------------------
def main():
    fig, ax = plt.subplots()

    for color, dname, mval in zip(colors, mass_dirs, mass_vals):
        d     = os.path.join(os.getcwd(), dname)
        logf  = os.path.join(d, "logfile.txt")
        # pick the first .off file in the folder
        off_list = glob.glob(os.path.join(d, "hauser_cube.off"))
        if not (os.path.exists(logf) and off_list):
            continue
        area     = mesh_area(off_list[0])
        strength = read_strength(logf)
        header, data = read_log(logf)
        if header is None or data.size == 0:
            continue

        t_idx = header.index("Time")
        e_idx = header.index("AdhesionEnergy")
        t      = data[:, t_idx]
        E      = data[:, e_idx]

        wrap   = -E / (strength * area)
        t_norm = t / t.max()

        ax.plot(t_norm, wrap, '-', color=color,
                label=f"{mval}")

    ax.set_xlabel("Normalized time")
    ax.set_ylabel("Wrapping fraction")

    # x-axis ticks: 0.0,0.2,…,1.0
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    # y-axis: coarser ticks
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    # tick styling
    ax.tick_params(axis='both', direction='in',
                   which='major', length=6)
    ax.tick_params(axis='both', direction='in',
                   which='minor', length=4)

    ax.legend(frameon=False, loc='best')

    # draw legend at lower right
    legend = ax.legend(frameon=False, loc='lower right')

# now add a math textbox just to its left, down near the bottom
    textstr = r"$N_p$"
    props   = dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8)

# tune x,y so the box sits just left of the legend in axes coords:
    ax.text(
    0.58, 0.15,            # x=65% across, y=10% up from bottom
    textstr,
    transform=ax.transAxes,
    fontsize=12,
    verticalalignment="bottom",
    horizontalalignment="left",
    #bbox=props
    )

    plt.tight_layout(pad=0.0)
    plt.savefig("wrapping_vs_time_mass_ratio.jpg", dpi=1200, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()
