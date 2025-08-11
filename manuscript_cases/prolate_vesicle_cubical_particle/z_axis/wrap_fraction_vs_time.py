import os
import numpy as np
import igl  # Make sure you have the libigl Python binding installed (pip install igl)
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
import matplotlib.font_manager as fm
# --------------------
# Set style parameters
# --------------------
plt.rcParams['font.family'] = 'Arial'
mpl.rcParams['figure.dpi'] = 600
plt.rcParams['font.size'] = 14

# Choose a qualitative colormap for distinct curves
cmap = plt.get_cmap('tab10')

# --------------------
# Functions to read and process data
# --------------------
def read_off(filename):
    """Reads a triangulated OFF file and returns vertices and faces."""
    with open(filename, 'r') as file:
        lines = file.readlines()
    lines = [l.strip() for l in lines if l.strip() and not l.strip().startswith('#')]
    if lines[0] != 'OFF':
        raise ValueError("Not a valid OFF file: missing 'OFF' header.'")
    parts = lines[1].split()
    n_vertices, n_faces = int(parts[0]), int(parts[1])
    vertices = np.array([list(map(float, lines[i].split()))
                         for i in range(2, 2 + n_vertices)])
    faces = []
    for i in range(2 + n_vertices, 2 + n_vertices + n_faces):
        parts = list(map(int, lines[i].split()))
        if parts[0] != 3:
            raise ValueError(f"Non-triangular face in {filename}.")
        faces.append(parts[1:4])
    return vertices, np.array(faces)

def compute_mesh_area_libigl(off_file):
    """Computes total surface area of a triangulated mesh via libigl."""
    V, F = read_off(off_file)
    result = igl.doublearea(V, F)
    double_areas = result[0] if isinstance(result, tuple) else result
    return np.sum(double_areas) / 2.0

def read_simulation_parameters(filename):
    """
    Reads simulation parameters from the log file.
    Expects a line with 'Particle adhesion strength:' and extracts its value.
    """
    strength = None
    with open(filename, 'r') as f:
        for line in f:
            if "Particle adhesion strength:" in line:
                try:
                    strength = float(line.split(":")[1].strip())
                except Exception as e:
                    raise ValueError(f"Cannot parse adhesion strength: {line}") from e
                break
    if strength is None:
        raise ValueError("Particle adhesion strength not found.")
    return strength

def read_log_file(filename):
    """
    Reads the simulation logfile.
    Expects a header line beginning with 'Iteration' followed by rows of numerical data.
    """
    header = None
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("Iteration"):
                header = line.split()
                continue
            if header:
                try:
                    data.append(list(map(float, line.split())))
                except Exception:
                    pass
    if header is None:
        return None, None
    return header, np.array(data)

# --------------------
# Main plotting function
# --------------------
def main_multi():
    base_path = "/mnt/c/Users/didarula/Desktop/production_run/hauser_cube_prolate/y_axis"
    u_dirs = ["u=2.0", "u=4.0", "u=6.0", "u=8.0", "u=10.0"]
    
    fig, ax = plt.subplots(figsize=(10, 7), dpi=600)
    markers = ['o', 's', '^', 'D', 'v', 'p', 'h']
    
    for i, u_folder in enumerate(u_dirs):
        folder_path = os.path.join(base_path, u_folder)
        logfile    = os.path.join(folder_path, "logfile.txt")
        off_file   = os.path.join(folder_path, "shifted_cube_prolate_y_axis.off")
        
        if not os.path.exists(logfile) or not os.path.exists(off_file):
            print(f"Missing files in {u_folder}, skipping.")
            continue
        
        # Read parameters & geometry
        u_strength = read_simulation_parameters(logfile)
        area       = compute_mesh_area_libigl(off_file)
        header, data = read_log_file(logfile)
        if header is None:
            print(f"No header in {u_folder}, skipping.")
            continue
        
        # Extract columns
        try:
            t_idx = header.index("Time")
            e_idx = header.index("AdhesionEnergy")
        except ValueError:
            print(f"Required columns not found in {u_folder}.")
            continue
        
        times = data[:, t_idx]
        energies = data[:, e_idx]
        wrap_frac = -energies / (u_strength * area)
        
        # —— Normalize time to [0,1] —— 
        norm_times = times / times.max()
        
        color  = cmap(i % cmap.N)
        marker = markers[i % len(markers)]
        
        ax.plot(
            norm_times,
            wrap_frac,
            '-.',
            color=color,
            # marker=marker,
            # markerfacecolor='none',
            # markersize=6,
            # markevery=100000,
            label=u_folder
        )
    
    ax.set_xlabel("Normalized time ($t_{normalized}$)",fontsize=14)
    ax.set_ylabel("Wrapping fraction ($\chi$)",fontsize=14)
    ax.legend(frameon=False, fontsize=16)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', direction='in', which='major', length=6)
    ax.tick_params(axis='both', direction='in', which='minor', length=4)
    
    plt.tight_layout()
    plt.savefig("wrapping_fraction_vs_norm_time.jpg", format="jpg", dpi=600)
    plt.close()
    print("Figure saved as wrapping_fraction_vs_norm_time.jpg")

if __name__ == '__main__':
    main_multi()
