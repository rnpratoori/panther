"""Create an initial-condition Exodus file from the last step of another file."""

import meshio
import numpy as np
import sys

# --- USER INPUTS --- #
# This is the only section you need to modify.

# List of circle definitions: [(center_x, center_y, radius), ...]
CIRCLES = [
    (0.25, 0.25, 0.1),
    (0.75, 0.75, 0.1),
    (0.25, 0.75, 0.15),
]

# The name of your multi-timestep exodus file.
INPUT_FILE = 'output/2phase_0.2.e'

# The name of the new exodus file for the initial condition.
OUTPUT_FILE = 'initial_condition.e'

# --- END OF INPUTS --- #


def create_ic_from_last_step(input_path, output_path, circle_params):
    """
    Reads the last time step from a multi-timestep Exodus file, adds a new 
    variable 'eta' based on circular domains, modifies 'c' and 'c2' to be 
    zero within those domains, and saves it as a new single-state Exodus file.
    """
    print(f"➡️ Reading mesh and LAST time step data from: {input_path}")
    try:
        # meshio.read on a transient Exodus file automatically loads the data 
        # from the final time step into mesh.point_data.
        mesh = meshio.read(input_path)
    except FileNotFoundError:
        print(f"❌ ERROR: Input file not found at '{input_path}'")
        sys.exit(1)

    # Extract node coordinates
    points = mesh.points
    x_coords = points[:, 0]
    y_coords = points[:, 1]
    
    # --- Step 1: Create and define the 'eta' variable ---
    print("   - Adding 'eta' variable and initializing all nodes to -1.")
    eta = np.full(x_coords.shape, -1.0)

    # Create a boolean mask to identify nodes within any of the specified circles.
    in_any_circle_mask = np.full(x_coords.shape, False, dtype=bool)
    for cx, cy, r in circle_params:
        print(f"   - Processing circle at ({cx}, {cy}) with radius {r}.")
        distance_sq = (x_coords - cx)**2 + (y_coords - cy)**2
        in_any_circle_mask |= (distance_sq < r**2)
    
    # Set eta = 1 for all nodes inside any circle.
    eta[in_any_circle_mask] = 1.0

    # --- Step 2: Modify 'c' and 'c2' based on 'eta' ---
    if 'c' not in mesh.point_data or 'c2' not in mesh.point_data:
        print("❌ ERROR: Input file's last time step must contain 'c' and 'c2' as nodal variables.")
        print("   Available point_data:", list(mesh.point_data.keys()))
        sys.exit(1)

    c_var = mesh.point_data['c']
    c2_var = mesh.point_data['c2']
    print("   - Modifying 'c' and 'c2' variables...")
    
    # Where the mask is True (eta > 0), set c and c2 to zero.
    c_var[in_any_circle_mask] = 0.0
    c2_var[in_any_circle_mask] = 0.0

    # --- Step 3: Write the new single-state Exodus file ---
    # Add the new/modified variables back to the mesh object.
    mesh.point_data['eta'] = eta
    mesh.point_data['c'] = c_var
    mesh.point_data['c2'] = c2_var

    print(f"➡️ Writing single-state initial condition file to: {output_path}")
    
    # Explicitly writing in "exodus" format ensures Paraview compatibility.
    meshio.write(output_path, mesh, file_format="exodus")
    
    print("\n✅ Done! Your initial condition file based on the last time step is ready.")


if __name__ == '__main__':
    create_ic_from_last_step(INPUT_FILE, OUTPUT_FILE, CIRCLES)
