import numpy as np
import matplotlib.pyplot as plt

def plot_com_distance(filename='com_distance.txt'):
    """
    Reads a text file with two columns (time, distance) and plots distance vs. time.
    """
    # 1) Load the data, ignoring lines that start with '%'
    data = np.loadtxt(filename, comments='%')
    
    # 2) Separate columns
    time = data[:, 0]
    dist = data[:, 1]
    
    # 3) Create the plot
    plt.figure(figsize=(8,6))
    plt.plot(time, dist, '-o', markersize=5, linewidth=2)
    plt.xlabel('Time')
    plt.ylabel('Distance between COMs')
    plt.title('Vesicle-Particle COM Distance vs. Time')
    plt.grid(True)
    
    # 4) Show the plot (or save with plt.savefig('distance_vs_time.png'))
    plt.show()

if __name__ == '__main__':
    plot_com_distance()
