import numpy as np

# Function to extract the initial principal axes matrix from the file content
def extract_initial_principal_axes(file_lines):
    principal_axes = []
    for idx, line in enumerate(file_lines):
        if "Principal axes:" in line:
            # Extract the next three lines containing the principal axes
            matrix_lines = file_lines[idx + 1: idx + 4]
            # Convert the lines to a 3x3 matrix and ensure column-based interpretation
            principal_axes = np.array([list(map(float, row.split())) for row in matrix_lines]).T
            break
    return principal_axes

# Function to parse rotation matrices from the file content
def extract_rotation_matrices(file_lines):
    rotation_matrices = []
    reading_matrix = False
    matrix_lines = []
    
    for line in file_lines:
        if "Rotation Matrix:" in line:
            reading_matrix = True
            matrix_lines = []
        elif reading_matrix:
            matrix_lines.append(line.strip())
            if len(matrix_lines) == 3:
                # Convert to a numpy array, assuming the axes are in columns
                matrix = np.array([list(map(float, row.split())) for row in matrix_lines]).T
                rotation_matrices.append(matrix)
                reading_matrix = False
    
    return rotation_matrices

# Read the log file content
file_path = '/mnt/c/Users/didarula/Desktop/production_run/hauser_cube_prolate/z_axis/u=8.0/u=8.0.log'
with open(file_path, 'r') as file:
    file_content = file.readlines()

# Extract the initial principal axes matrix from the log file
initial_principal_axes = extract_initial_principal_axes(file_content)

# Extract rotation matrices from the file content
rotation_matrices = extract_rotation_matrices(file_content)
#rotation_matrices = rotation_matrices[1:]

# Transform the principal axes using rotation matrices
transformed_principal_axes_list = [initial_principal_axes]  # Include initial axes
for rotation_matrix in rotation_matrices:
    transformed_principal_axes = rotation_matrix# @ initial_principal_axes
    transformed_principal_axes_list.append(transformed_principal_axes)

# Construct the data file format with timesteps and transformed principal axes
timesteps = np.arange(0, len(transformed_principal_axes_list) * 1, 1)  # Assuming a timestep increment of 100
output_data = []

for i, transformed_axes in enumerate(transformed_principal_axes_list):
    flattened_axes = transformed_axes.flatten()
    output_data.append([timesteps[i], *flattened_axes])

output_data = np.array(output_data)

# Save the transformed principal axes for each timestep to a text file
formatted_output_file_path = '/mnt/c/Users/didarula/Desktop/production_run/hauser_cube_prolate/z_axis/u=8.0/transformed_principal_axes.txt'
header = "# TimeStep v_ax1 v_ay1 v_az1 v_ax2 v_ay2 v_az2 v_ax3 v_ay3 v_az3"
np.savetxt(formatted_output_file_path, output_data, header=header, comments='', fmt='%.5f')
print(f"Angular changes have been successfully saved to: {formatted_output_file_path}")

