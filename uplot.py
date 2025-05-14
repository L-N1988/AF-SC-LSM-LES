import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Base directory
base_dir = "./postProcessing/velocityPlaneX0"

# Get list of time directories (only numeric ones)
time_dirs = sorted(
        [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)) and d.replace('.', '', 1).isdigit()],
        key=lambda x: float(x)
        )

# Iterate through each time directory
for time_dir in time_dirs:
    raw_path = os.path.join(base_dir, time_dir, "U_x0Plane.raw")

    if not os.path.exists(raw_path):
        print(f"Skipping {time_dir}: File not found - {raw_path}")
        continue

    with open(raw_path, 'r') as file:
        # Skip the first metadata line
        file.readline()

        # Read the header line and split into column names
        headers = file.readline().strip().split()

        # Read and parse numerical data
        data = []
        for line in file:
            line = line.strip()
            if line:
                data.append([float(value) for value in line.split()])

    print(f"Processed time = {time_dir}")
    # Convert to NumPy array
    data = np.array(data)

    # Extract columns
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    Ux = data[:, 3]
    Uy = data[:, 4]
    Uz = data[:, 5]

    # Assume y and z form a regular grid
    unique_y = np.unique(y)
    unique_z = np.unique(z)

    # Reshape U fields to match the grid
    Ux_grid = Ux.reshape(len(unique_y), len(unique_z))
    Uy_grid = Uy.reshape(len(unique_y), len(unique_z))
    Uz_grid = Uz.reshape(len(unique_y), len(unique_z))

    output_folder = 'centerPlaneU'
    # Plotting function
    def plot_and_save(field, title, filename):
        plt.figure(figsize=(8, 6))
        plt.imshow(field, extent=[unique_z.min(), unique_z.max(), unique_y.max(), unique_y.min()],
                aspect='auto', cmap='viridis')
        plt.colorbar(label=title)
        plt.xlabel('Z')
        plt.ylabel('Y')
        plt.title(f'Velocity Component: {title}')
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, filename))
        plt.close()

    # Generate and save plots
    plot_and_save(Ux_grid, 'U_x', 'Ux{time_dir}.png')
    plot_and_save(Uy_grid, 'U_y', 'Uy{time_dir}.png')
    plot_and_save(Uz_grid, 'U_z', 'Uz{time_dir}.png')
    print(f"Plots saved in folder '{output_folder}'")
    print(f"Error reading {time_dir}: {e}")

