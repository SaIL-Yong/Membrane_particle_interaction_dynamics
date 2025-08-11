import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter

# --------------------
# Publication‐style settings
# --------------------
cm2in = lambda cm: cm / 2.54
fig_w_cm, fig_h_cm = 8.30, 10.0
fig_w_in, fig_h_in = cm2in(fig_w_cm), cm2in(fig_h_cm)

mpl.rcParams.update({
    'figure.figsize':    (fig_w_in, fig_h_in),
    'figure.dpi':        600,
    'font.family':       'Arial',
    'font.size':         10,
    'axes.labelsize':    10,
    'xtick.labelsize':   10,
    'ytick.labelsize':   10,
    'legend.fontsize':   8,
    'lines.linewidth':   1,
})

# --------------------
# Simulation parameters
# --------------------
u_vals = [2.0, 4.0, 6.0, 8.0, 10.0]
u_dirs = [f"u={v}" for v in u_vals]
cmap   = mpl.cm.plasma
colors = cmap(np.linspace(0.1, 0.9, len(u_vals)))

# --------------------
# Log‐file reader
# --------------------
def read_log(path):
    header = None
    rows   = []
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
# Plotting helper
# --------------------
def plot_bending(ax, base_path, show_legend=False):
    for color, u_dir, u_mod in zip(colors, u_dirs, u_vals):
        logf = os.path.join(base_path, u_dir, "logfile.txt")
        if not os.path.exists(logf):
            continue

        header, data = read_log(logf)
        if header is None or data.size == 0:
            continue

        t_i  = header.index("Time")
        be_i = header.index("BendingEnergy")
        t    = data[:, t_i]
        BE   = data[:, be_i]

        # normalize
        t_norm  = t / t.max()
        BE_norm = BE / (8 * np.pi * 0.01)

        ax.plot(t_norm, BE_norm, '-', color=color,
                label=rf"$u_{{mod}}={u_mod}$")

    # only y‐axis label on both
    ax.set_ylabel("Bending energy")
    
    # x‐axis: exactly as before
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    # y‐axis ticks
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    # separate calls for major and minor tick styling
    ax.tick_params(axis='both', direction='in', which='major', length=6)
    ax.tick_params(axis='both', direction='in', which='minor', length=4)

    if show_legend:
        ax.legend(frameon=False, loc='best')

# --------------------
# Main routine
# --------------------
def main():
    fig, (ax_y, ax_z) = plt.subplots(
        2, 1,
        sharex=True,
        constrained_layout=True,
        gridspec_kw={'hspace': 0}
    )

    # adjust to your actual directories
    y_path = "/mnt/c/Users/didarula/Desktop/production_run/hauser_cube_prolate/y_axis"
    z_path = "/mnt/c/Users/didarula/Desktop/production_run/hauser_cube_prolate/z_axis"

    plot_bending(ax_y, y_path, show_legend=False)
    plot_bending(ax_z, z_path, show_legend=True)

    # only bottom plot gets the x‐label
    ax_y.tick_params(labelbottom=False)
    ax_z.set_xlabel("Normalized time")

    plt.subplots_adjust(hspace=0.0)
    plt.tight_layout(pad=0.0)
    plt.savefig("bending_energy_normalized_vs_time.jpg", dpi=600, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()
# This script reads bending energy data from simulation log files and plots it against normalized time.