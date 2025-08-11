import os
import numpy as np
import igl
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter

# --------------------
# Publication-style settings
# --------------------
cm2in = lambda cm: cm / 2.54
fig_width_cm, fig_height_cm = 8.30, 10.0
fig_width_in, fig_height_in = cm2in(fig_width_cm), cm2in(fig_height_cm)

mpl.rcParams.update({
    'figure.figsize':    (fig_width_in, fig_height_in),
    'figure.dpi':        600,
    'font.family':       'Arial',
    'font.size':         12,
    'axes.labelsize':    12,
    'xtick.labelsize':   10,
    'ytick.labelsize':   10,
    'legend.fontsize':   10,
    'lines.linewidth':   1,
})

# --------------------
# Simulation parameters
# --------------------
u_vals = [2.0, 4.0, 6.0, 8.0, 10.0]
u_dirs = [f"u={v}" for v in u_vals]
colors = mpl.cm.plasma(np.linspace(0.1, 0.9, len(u_vals)))

# --------------------
# Helpers
# --------------------
def read_off(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
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
                return float(line.split(":",1)[1])
    raise ValueError("Strength line not found")

def read_log(path):
    header, rows = None, []
    with open(path) as f:
        for L in f:
            L = L.strip()
            if L.startswith("Iteration"):
                header = L.split()
            elif header:
                parts = L.split()
                if len(parts) == len(header):
                    rows.append(list(map(float, parts)))
    return header, np.array(rows)

# --------------------
# Plot helper
# --------------------
def plot_wrapping(ax, base_path, show_legend=False):
    for color, u_dir, u_mod in zip(colors, u_dirs, u_vals):
        d    = os.path.join(base_path, u_dir)
        logf = os.path.join(d, "logfile.txt")
        offp = os.path.join(d, "vertical_superspheroid.off")
        if not (os.path.exists(logf) and os.path.exists(offp)):
            continue

        strength = read_strength(logf)
        area     = mesh_area(offp)
        header, data = read_log(logf)
        if header is None or data.size == 0:
            continue

        t      = data[:, header.index("Time")]
        E      = data[:, header.index("AdhesionEnergy")]
        wrap   = -E / (strength*area)
        t_norm = t / t.max()

        ax.plot(t_norm, wrap, '-', color=color,
                label=rf"$u_{{mod}}={u_mod}$")

    ax.set_ylabel("Wrapping fraction")
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', direction='in', which='major', length=6)
    ax.tick_params(axis='both', direction='in', which='minor', length=4)

    if show_legend:
        ax.legend(frameon=False, loc='best')

# --------------------
# Main
# --------------------
def main():
    y_path = "/mnt/c/Users/didarula/Desktop/production_run/rod_flat_sperical"
    z_path = "/mnt/c/Users/didarula/Desktop/production_run/rod_spherical"

    fig, (ax_y, ax_z) = plt.subplots(
        2, 1,
        sharex=True,
        gridspec_kw={'hspace': 0},
        figsize=(fig_width_in, fig_height_in),
        dpi=600
    )

    plot_wrapping(ax_y, y_path, show_legend=False)
    plot_wrapping(ax_z, z_path, show_legend=True)

    ax_y.spines['bottom'].set_visible(False)
    ax_z.spines['top'].set_visible(False)

    ax_y.tick_params(labelbottom=False)
    ax_z.set_xlabel("Normalized time")

    plt.savefig("wrapping_fraction_vs_time.png", dpi=600, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()
