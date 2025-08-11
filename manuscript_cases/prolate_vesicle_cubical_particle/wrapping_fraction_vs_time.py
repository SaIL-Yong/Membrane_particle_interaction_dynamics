import os
import numpy as np
import igl
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter

# --------------------
# Publication‐style settings
# --------------------
cm2in = lambda cm: cm / 2.54
fig_width_cm  = 8.30             # single‐column width for Soft Matter
fig_height_cm = 10.0             # arbitrary height
fig_width_in  = cm2in(fig_width_cm)
fig_height_in = cm2in(fig_height_cm)

mpl.rcParams.update({
    'figure.figsize':       (fig_width_in, fig_height_in),
    'figure.dpi':           600,
    'font.family':          'Arial',
    'font.size':            12,
    'axes.labelsize':       12,
    'xtick.labelsize':      12,
    'ytick.labelsize':      12,
    'legend.fontsize':      10,
    'lines.linewidth':      1,
})

base_cmap = plt.get_cmap('plasma')

# numeric values of u_mod, and matching directory names "u=2.0", etc.
u_vals = [2.0, 4.0, 6.0, 8.0, 10.0]
u_dirs = [f"u={v}" for v in u_vals]

# sample five distinct colors from plasma (avoid the endpoints)
colors = base_cmap(np.linspace(0.1, 0.9, len(u_vals)))

# --------------------
# File readers
# --------------------
def read_off(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    if lines[0] != 'OFF':
        raise ValueError("Missing OFF header")
    nv, nf = map(int, lines[1].split()[:2])
    verts = np.array([list(map(float, l.split()))
                      for l in lines[2:2+nv]])
    faces = np.array([list(map(int, l.split()[1:4]))
                      for l in lines[2+nv:2+nv+nf]])
    return verts, faces

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
# Plot helper
# --------------------
def plot_axis(ax, base_path, show_legend=False):
    for col, u_dir, u_mod in zip(colors, u_dirs, u_vals):
        d    = os.path.join(base_path, u_dir)
        logf = os.path.join(d, "logfile.txt")
        offp = os.path.join(d,
            "shifted_cube_prolate_y_axis.off"
            if 'y_axis' in base_path
            else "shifted_cube_prolate_z_axis.off"
        )
        if not (os.path.exists(logf) and os.path.exists(offp)):
            continue

        strength = read_strength(logf)
        area     = mesh_area(offp)
        header, data = read_log(logf)
        if header is None:
            continue

        t_i = header.index("Time")
        e_i = header.index("AdhesionEnergy")
        t   = data[:, t_i]
        E   = data[:, e_i]

        wrap   = -E / (strength * area)
        t_norm = t / t.max()

        # Legend label in LaTeX math mode
        label = rf'$u_{{mod}}={u_mod}$'
        ax.plot(t_norm, wrap, '-', color=col, label=label)

    ax.set_ylabel("Wrapping fraction")

    # Major ticks every 0.2 and formatted to one decimal place
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # Minor ticks
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.tick_params(axis='both', direction='in', which='major', length=6)
    ax.tick_params(axis='both', direction='in', which='minor', length=4)

    if show_legend:
        ax.legend(
            frameon=False,
            loc='center left',
            bbox_to_anchor=(0.1, 0.6),
            borderaxespad=0
        )

# --------------------
# Main
# --------------------
def main():
    fig, (ax_y, ax_z) = plt.subplots(
        2, 1,
        sharex=True,
        constrained_layout=True,
        gridspec_kw={'hspace': 0}
    )

    y_path = "/mnt/c/Users/didarula/Desktop/production_run/hauser_cube_prolate/y_axis"
    z_path = "/mnt/c/Users/didarula/Desktop/production_run/hauser_cube_prolate/z_axis"

    plot_axis(ax_y, y_path, show_legend=False)
    plot_axis(ax_z, z_path, show_legend=True)

    # Eliminate vertical space and tighten margins
    plt.subplots_adjust(hspace=0.0)
    plt.tight_layout(pad=0.0)

    ax_z.set_xlabel("Normalized time")
    ax_y.tick_params(labelbottom=False)

    plt.savefig("wrapping_fraction_y_vs_z.jpg", dpi=600, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()