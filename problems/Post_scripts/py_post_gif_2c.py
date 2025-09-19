import os
import netCDF4
import numpy as np
import pyvista as pv
import vtk
from tqdm import tqdm

def main(exodus_filename, output_filename):
    # Open the Exodus file
    ds = netCDF4.Dataset(exodus_filename)
    times = ds.variables['time_whole'][:]
    nt = len(times)

    X = np.ma.getdata(ds.variables['coordx'][:])
    Y = np.ma.getdata(ds.variables['coordy'][:])
    Z = np.zeros_like(X)
    points = np.vstack([X, Y, Z]).T

    connect_vars = [v for v in ds.variables if v.startswith('connect')]
    block_meshes = []
    for connect_var in connect_vars:
        elem_node = np.ma.getdata(ds.variables[connect_var][:]) - 1
        mesh = pv.UnstructuredGrid({vtk.VTK_QUAD: elem_node}, points)
        block_meshes.append(mesh)

    base_mesh = block_meshes[0]
    for mesh in block_meshes[1:]:
        base_mesh = base_mesh.merge(mesh, merge_points=True)

    c_all = ds.variables['vals_nod_var1'][:]   # shape: (nt, n_nodes)
    c2_all = ds.variables['vals_nod_var2'][:]  # shape: (nt, n_nodes)
    ds.close()

    mesh = base_mesh.copy(deep=True)

    # PyVista plotter (off-screen movie mode)
    plotter = pv.Plotter(off_screen=True)
    plotter.open_movie(output_filename, framerate=48)
    plotter.view_xy()

    camera_set = False
    for i in tqdm(range(nt), desc="Generating frames"):
        c = c_all[i]
        c2 = c2_all[i]

        # Map c to red, c2 to green, blue=0
        rgb_colors = np.clip(np.stack([c, c2, np.zeros_like(c)], axis=1), 0, 1)
        mesh.point_data['rgb'] = (rgb_colors * 255).astype(np.uint8)

        plotter.add_mesh(mesh, scalars='rgb', rgb=True, show_scalar_bar=False)
        if not camera_set:
            plotter.view_xy()
            plotter.camera.Zoom(1.0)
            camera_set = True

        plotter.write_frame()
        plotter.clear_actors()

    plotter.close()

if __name__ == "__main__":
    input_dir = "output"
    for exodus_file in sorted(os.listdir(input_dir)):
        if exodus_file.endswith('2phase_taylor.e'):
            exodus_path = os.path.join(input_dir, exodus_file)
            output_filename = os.path.join(input_dir, exodus_file.replace('.e', '_2c_rg.mp4'))
            print(f"Processing {exodus_file}...")
            main(exodus_path, output_filename)
            print(f"Created {output_filename}")