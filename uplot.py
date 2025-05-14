import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

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

    # Interpolate over 2d mesh
    interpolatorUx = interp.CloughTocher2DInterpolator(data[:, 1:3], Ux)
    interpolatorUy = interp.CloughTocher2DInterpolator(data[:, 1:3], Uy)
    interpolatorUz = interp.CloughTocher2DInterpolator(data[:, 1:3], Uz)

    # go linearly in the y grid and z grid
    yline = np.linspace(min(y), may(y), len(np.unique(y)))
    zline = np.linspace(min(z), maz(z), len(np.unique(z)))
    # construct 2d grid from these
    xgrid,ygrid = np.meshgrid(xline, yline)
    # interpolate z data; same shape as xgrid and ygrid
    Ux_grid = interpolatorUx(xgrid, ygrid)
    Uy_grid = interpolatorUy(xgrid, ygrid)
    Uz_grid = interpolatorUz(xgrid, ygrid)

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
    plot_and_save(Ux_grid, 'U_x', f"Ux_'{time_dir}'.png")
    plot_and_save(Uy_grid, 'U_y', f"Uy_'{time_dir}'.png")
    plot_and_save(Uz_grid, 'U_z', f"Uz_'{time_dir}'.png")
    print(f"Plots saved in folder '{output_folder}'")