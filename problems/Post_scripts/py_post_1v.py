"""Render a single nodal variable from Exodus output to an off-screen movie."""

import os
import netCDF4
import numpy as np
import pyvista as pv
import vtk
from tqdm import tqdm

def main(exodus_filename, output_filename, variable_name):
    ds = netCDF4.Dataset(exodus_filename)

    times = ds.variables['time_whole'][:]
    nt = len(times)

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

    # Read c1 and c2 data (assumes var1 is c1, var2 is c2)
    c1_all = ds.variables['vals_nod_var1'][:]
    c2_all = ds.variables['vals_nod_var2'][:]
    c3_all = 1 - c1_all - c2_all
    ds.close()

    mesh_c3 = base_mesh.copy(deep=True)
    mesh_c3.point_data[variable_name] = c3_all[0]

    plotter = pv.Plotter(off_screen=True)
    plotter.open_gif(output_filename)

    vmin = float(np.min(c3_all[-1]))
    vmax = float(np.max(c3_all[-1]))

    plotter.add_mesh(mesh_c3, scalars=variable_name, clim=(vmin, vmax),
                     show_scalar_bar=True, cmap='rainbow', show_edges=True,
                     edge_color='black', line_width=0.5)
    plotter.show_grid()
    plotter.view_xy()

    for i in tqdm(range(len(times)), desc=f"Generating frames for {variable_name}"):
        mesh_c3.point_data[variable_name] = c3_all[i]
        plotter.write_frame()

    plotter.close()
    return True
        
if __name__ == "__main__":
    input_dir = "/home/rnp/MOOSE/projects/panther/problems/output"
    variable = "c3"  # Start with a known variable name
    
    for exodus_file in sorted(os.listdir(input_dir)):
        if exodus_file.endswith('_bcic_t2.e'):
            exodus_path = os.path.join(input_dir, exodus_file)
            output_filename = os.path.join(input_dir, 
                                         exodus_file.replace('.e', f'_{variable}.gif'))
            
            print(f"\nProcessing {exodus_file}...")
            if main(exodus_path, output_filename, variable):
                print(f"Successfully created {output_filename}")
            else:
                print(f"Failed to process {exodus_file}")
