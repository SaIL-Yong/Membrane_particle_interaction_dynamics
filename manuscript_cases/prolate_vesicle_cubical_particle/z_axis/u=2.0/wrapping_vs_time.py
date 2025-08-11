import numpy as np
import igl  # Ensure you have the libigl Python binding installed (pip install igl)
import matplotlib.pyplot as plt
import os

def compute_mesh_area_libigl(off_file):
    """
    Computes the total surface area of a triangulated mesh using libigl.
    
    This function reads the mesh using igl.readOFF(), computes per-face double areas
    with igl.doublearea(), and returns the total mesh area.
    """
    # Load vertices (V) and faces (F) using libigl's built-in function.
    V, F = igl.readOFF(off_file)
    # Compute per-face double areas (each value equals 2 * area of the face)
    total_area = np.sum(igl.doublearea(V, F)) / 2.0
    return total_area

def read_log_file(filename):
    """
    Reads the simulation log file.
    Expects a header line starting with 'Iteration' followed by rows of numerical data.
    """
    header = None
    data = []
    with open(filename, 'r') as f:
        found_header = False
        for line in f:
            line = line.strip()
            if line.startswith("Iteration"):
                header = line.split()
                found_header = True
                continue
            if not found_header:
                continue
            try:
                row = list(map(float, line.split()))
                data.append(row)
            except Exception:
                continue
    data = np.array(data)
    return header, data

def main():
    # File paths (update if needed)
    logfile = "/mnt/c/Users/didarula/Desktop/production_run/hauser_cube_prolate/z_axis/u=2.0/logfile.txt"
    off_file = "/mnt/c/Users/didarula/Desktop/production_run/hauser_cube_prolate/z_axis/u=2.0/shifted_cube_prolate_z_axis.off"

    # Verify file existence
    if not os.path.exists(off_file):
        print("OFF file not found:", off_file)
        return
    if not os.path.exists(logfile):
        print("Log file not found:", logfile)
        return

    # Simulation parameter from log (particle adhesion strength)
    particle_adhesion_strength = 0.168827

    # Compute the particle area using libigl's built-in functions
    try:
        particle_area = compute_mesh_area_libigl(off_file)
    except Exception as e:
        print("Error computing mesh area with libigl:", e)
        return
    print(f"Particle area computed from '{off_file}': {particle_area}")

    # Read simulation log data
    header, data = read_log_file(logfile)
    if header is None:
        print("Header not found in the log file. Please check that the file contains a proper header.")
        return

    # Identify column indices for 'Time' and 'AdhesionEnergy'
    try:
        time_idx = header.index("Time")
        adhesion_energy_idx = header.index("AdhesionEnergy")
    except ValueError:
        print("Could not find 'Time' and/or 'AdhesionEnergy' in the log file header.")
        return

    # Extract time and adhesion energy from the data
    times = data[:, time_idx]
    adhesion_energies = data[:, adhesion_energy_idx]

    # Compute the wrapping fraction using:
    # wrapping_fraction = -(AdhesionEnergy) / (Adhesion Strength * particle_area)
    wrapping_fraction = -adhesion_energies / (particle_adhesion_strength * particle_area)

    # Plot time vs. wrapping fraction
    plt.figure(figsize=(8, 6))
    plt.plot(times, wrapping_fraction, marker='o', linestyle='-', color='blue')
    plt.xlabel("Time")
    plt.ylabel("Wrapping Fraction")
    plt.title("Time vs. Wrapping Fraction")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
