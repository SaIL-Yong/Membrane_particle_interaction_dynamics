import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter

# --------------------
# Publication-style settings
# --------------------
cm2in = lambda cm: cm / 2.54
fig_width_cm, fig_height_cm = 8.30, 10.0      # single-column width × height
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
# Log-file reader
# --------------------
def read_log(path):
    header = None
    rows = []
    with open(path) as f:
        for line in f:
            L = line.strip()
            if L.startswith("Iteration"):
                header = L.split()
                continue
            if header:
                parts = L.split()
                if len(parts) == len(header):
                    try:
                        rows.append(list(map(float, parts)))
                    except ValueError:
                        pass
    return header, np.array(rows)

# --------------------
# Plot helper
# --------------------
def plot_bending(ax, base_path, show_legend=False):
    for color, u_dir, u_mod in zip(colors, u_dirs, u_vals):
        logfile = os.path.join(base_path, u_dir, "logfile.txt")
        if not os.path.exists(logfile):
            continue
        header, data = read_log(logfile)
        if header is None or data.size == 0:
            continue

        t        = data[:, header.index("Time")]
        BE       = data[:, header.index("BendingEnergy")]
        t_norm   = t / t.max()
        BE_norm  = BE / (8*np.pi*0.01)

        ax.plot(t_norm, BE_norm, '-', color=color,
                label=rf"$u_{{mod}}={u_mod}$")

    ax.set_ylabel("Bending energy")
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_locator(MultipleLocator(0.6))
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

    plot_bending(ax_y, y_path, show_legend=False)
    plot_bending(ax_z, z_path, show_legend=True)

    # seamless join
    ax_y.spines['bottom'].set_visible(False)
    ax_z.spines['top'].set_visible(False)

    ax_y.tick_params(labelbottom=False)
    ax_z.set_xlabel("Normalized time")

    plt.savefig("bending_vs_time.png", dpi=600, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()
