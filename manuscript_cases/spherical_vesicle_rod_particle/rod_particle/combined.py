import os
import numpy as np
import igl
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter, MaxNLocator

# --------------------
# Publication-style settings
# --------------------
cm2in = lambda cm: cm / 2.54
fig_w_cm, fig_h_cm = 3*8.30, 10.0      # three single-column widths side by side
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
# Common log-file reader
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
# Wrapping-fraction helpers
# --------------------
def read_off(fname):
    with open(fname) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    if lines[0] != 'OFF':
        raise ValueError("Missing OFF header")
    nv, nf = map(int, lines[1].split()[:2])
    verts = np.array([list(map(float, l.split())) for l in lines[2:2+nv]])
    faces = np.array([list(map(int, l.split()[1:4])) for l in lines[2+nv:2+nv+nf]])
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

# --------------------
# Plotting routines
# --------------------
def plot_bending(ax, base_path, show_legend=False):
    for color, u_dir, u_mod in zip(colors, u_dirs, u_vals):
        logfile = os.path.join(base_path, u_dir, "logfile.txt")
        if not os.path.exists(logfile): continue
        header, data = read_log(logfile)
        if header is None or data.size==0: continue

        t = data[:, header.index("Time")]
        BE= data[:, header.index("BendingEnergy")]
        t_norm  = t / t.max()
        BE_norm = BE / (8*np.pi*0.01)

        ax.plot(t_norm, BE_norm, '-', color=color, label=rf"$u_{{mod}}={u_mod}$")

    ax.set_ylabel("Bending energy")
    ax.set_xlim(0,1)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_locator(MultipleLocator(0.6))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', direction='in', which='major', length=6)
    ax.tick_params(axis='both', direction='in', which='minor', length=4)
    if show_legend:
        ax.legend(frameon=False, loc='best')

def plot_reduced(ax, base_path, show_legend=False):
    for color, u_dir, u_mod in zip(colors, u_dirs, u_vals):
        logfile = os.path.join(base_path, u_dir, "logfile.txt")
        if not os.path.exists(logfile): continue
        header, data = read_log(logfile)
        if header is None or data.size==0: continue

        t = data[:, header.index("Time")]
        RV= data[:, header.index("ReducedVolume")]
        t_norm = t / t.max()

        ax.plot(t_norm, RV, '-', color=color, label=rf"$u_{{mod}}={u_mod}$")

    ax.set_ylabel("Reduced volume")
    ax.set_xlim(0,1)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.tick_params(axis='both', direction='in', which='major', length=6)
    ax.tick_params(axis='both', direction='in', which='minor', length=4)
    if show_legend:
        ax.legend(frameon=False, loc='best')

def plot_wrapping(ax, base_path, show_legend=False):
    for color, u_dir, u_mod in zip(colors, u_dirs, u_vals):
        d = os.path.join(base_path, u_dir)
        logfile = os.path.join(d, "logfile.txt")
        offp    = os.path.join(d, "vertical_superspheroid.off")
        if not (os.path.exists(logfile) and os.path.exists(offp)): continue

        strength = read_strength(logfile)
        area     = mesh_area(offp)
        header, data = read_log(logfile)
        if header is None: continue

        t = data[:, header.index("Time")]
        E = data[:, header.index("AdhesionEnergy")]
        wrap = -E/(strength*area)
        t_norm = t/t.max()

        ax.plot(t_norm, wrap, '-', color=color, label=rf"$u_{{mod}}={u_mod}$")

    ax.set_ylabel("Wrapping fraction")
    ax.set_xlim(0,1)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', direction='in', which='major', length=6)
    ax.tick_params(axis='both', direction='in', which='minor', length=4)
    if show_legend:
        ax.legend(frameon=False, loc='best')

# --------------------
# Main
# --------------------
def main():
    # adjust these to your actual run directories:
    y_path = "/mnt/c/Users/didarula/Desktop/production_run/rod_flat_sperical"
    z_path = "/mnt/c/Users/didarula/Desktop/production_run/rod_spherical"

    fig, axs = plt.subplots(
        2, 3,
        sharex=True,
        constrained_layout=True
    )

    # Column 0: Bending
    plot_bending(axs[0,0], y_path, show_legend=False)
    plot_bending(axs[1,0], z_path, show_legend=True)

    # Column 1: Reduced Volume
    plot_reduced(axs[0,1], y_path, show_legend=False)
    plot_reduced(axs[1,1], z_path, show_legend=True)

    # Column 2: Wrapping Fraction
    plot_wrapping(axs[0,2], y_path, show_legend=False)
    plot_wrapping(axs[1,2], z_path, show_legend=True)

    # Only bottom row has x-labels
    for ax in axs[0, :]:
        ax.tick_params(labelbottom=False)
    for ax in axs[1, :]:
        ax.set_xlabel("Normalized time")

    
    plt.subplots_adjust(hspace=0.0)
    plt.tight_layout(pad=0.0)
    plt.savefig("combined_2x3_normalized_time.png", dpi=600, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()
