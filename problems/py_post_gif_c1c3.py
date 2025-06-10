import os
import netCDF4
import numpy as np
import pyvista as pv
import vtk
from tqdm import tqdm

def create_lookup_table(color='red'):
    """
    Create a custom vtkLookupTable that maps a scalar value in [0,1]
    to an RGBA color. The LUT is defined such that:
      - value 0 → light color with low opacity (transparent)
      - value 1 → dark color with full opacity (opaque)
    """
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(256)
    lut.SetRange(0, 1)
    for i in range(256):
        x = i / 255.0  # normalized scalar
        if color == 'red':
            r = 0.5 + 0.5 * x  # from 0.5 (light) to 1.0 (dark red)
            g = 0.0
            b = 0.0
        elif color == 'green':
            r = 0.0
            g = 0.5 + 0.5 * x  # from 0.5 (light) to 1.0 (dark green)
            b = 0.0
        elif color == 'blue':
            r = 0.0
            g = 0.0
            b = 0.5 + 0.5 * x  # from 0.5 (light) to 1.0 (dark blue)
        else:
            r = g = b = x
        a = x  # opacity: 0 for 0, 1 for 1
        lut.SetTableValue(i, r, g, b, a)
    lut.Build()
    return lut

def create_custom_colormap(color='red'):
    """
    Create a custom colormap that maps values from [0,1] to a color gradient
    with increasing opacity.
    """
    if color == 'red':
        base_color = [1, 0, 0]  # Red
    elif color == 'green':
        base_color = [0, 1, 0]  # Green
    elif color == 'blue':
        base_color = [0, 0, 1]  # Blue
    else:
        base_color = [1, 1, 1]  # White

    def custom_cmap(x):
        rgba = np.zeros((len(x), 4))
        rgba[:, :3] = np.array(base_color)  # RGB values
        rgba[:, 3] = x  # Alpha values
        return rgba

    return custom_cmap

def main(exodus_filename, output_filename):
    # Open the Exodus file
    ds = netCDF4.Dataset(exodus_filename)
    
    # Get the number of timesteps from the "time_whole" variable
    times = ds.variables['time_whole'][:]
    nt = len(times)
    
    # Read coordinates (assumed 2D; Z is set to zero)
    X = np.ma.getdata(ds.variables['coordx'][:])
    Y = np.ma.getdata(ds.variables['coordy'][:])
    Z = np.zeros_like(X)
    points = np.vstack([X, Y, Z]).T
    
    # Read all element blocks (connectivity) dynamically
    connect_vars = [v for v in ds.variables if v.startswith('connect')]
    block_meshes = []
    
    for connect_var in connect_vars:
        elem_node = np.ma.getdata(ds.variables[connect_var][:]) - 1
        mesh = pv.UnstructuredGrid({vtk.VTK_QUAD: elem_node}, points)
        block_meshes.append(mesh)
    
    # Combine all blocks into one mesh
    base_mesh = block_meshes[0]
    for mesh in block_meshes[1:]:
        base_mesh = base_mesh.merge(mesh, merge_points=True)
    
    # Assume that the Exodus file has nodal variables for c1 and c2 stored as:
    # 'vals_nod_var1' for c1 and 'vals_nod_var2' for c2.
    # Adjust these indices if needed.
    c1_all = ds.variables['vals_nod_var1'][:]  # shape: (nt, n_nodes)
    c3_all = ds.variables['vals_nod_var2'][:]  # c3 is now read from var2
    c2_all = 1 - c1_all - c3_all               # c2 is now computed
    
    ds.close()
    
    # Create three copies of the base mesh, one for each variable
    mesh_c1 = base_mesh.copy(deep=True)
    mesh_c3 = base_mesh.copy(deep=True)  # mesh_c3 now for c2
    mesh_c2 = base_mesh.copy(deep=True)  # mesh_c2 now for c3
    
    # Initialize scalar data for each mesh
    mesh_c1.point_data['c1'] = c1_all[0]  # Set initial data
    mesh_c3.point_data['c3'] = c3_all[0]  # mesh_c3 now for c2
    mesh_c2.point_data['c2'] = c2_all[0]  # mesh_c2 now for c3
    
    # Create custom colormaps
    cmap_red = create_custom_colormap('red')
    cmap_green = create_custom_colormap('green')
    cmap_blue = create_custom_colormap('blue')
    
    # Create a PyVista plotter in off-screen mode
    plotter = pv.Plotter(off_screen=True)
    plotter.open_movie(output_filename, framerate=60)

    # Add the three meshes with custom colormaps
    plotter.add_mesh(mesh_c1, scalars='c1', clim=(0, 1), 
                     show_scalar_bar=False, opacity=0.7,
                     cmap=cmap_red, render=False)
    plotter.add_mesh(mesh_c3, scalars='c3', clim=(0, 1), 
                     show_scalar_bar=False, opacity=0.7,
                     cmap=cmap_green, render=False)
    plotter.add_mesh(mesh_c2, scalars='c2', clim=(0, 1), 
                     show_scalar_bar=False, opacity=0.7,
                     cmap=cmap_blue, render=False)

    # Set a fixed camera view (2D view)
    plotter.view_xy()

    # Loop over timesteps and update scalar data for each layer
    for i in tqdm(range(nt), desc="Generating frames"):
        mesh_c1.point_data['c1'] = c1_all[i]
        mesh_c3.point_data['c3'] = c3_all[i]
        mesh_c2.point_data['c2'] = c2_all[i]
        plotter.write_frame()

    plotter.close()
    
if __name__ == "__main__":
    # Process all .e files in the results/output_dump_3p directory
    input_dir = "output"
    for exodus_file in sorted(os.listdir(input_dir)):
        if exodus_file.endswith('dis_spline.e'):
            exodus_path = os.path.join(input_dir, exodus_file)
            # Create output filename by replacing .e with .gif
            output_filename = os.path.join(input_dir, exodus_file.replace('.e', '_v3.mp4'))
            print(f"Processing {exodus_file}...")
            main(exodus_path, output_filename)
            print(f"Created {output_filename}")