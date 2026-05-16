"""Render ternary Exodus fields to an RGB composition movie."""

import os
import netCDF4
import numpy as np
import pyvista as pv
import vtk
from tqdm import tqdm
import matplotlib.colors as mcolors

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

    c1_all = ds.variables['vals_nod_var1'][:]
    c2_all = ds.variables['vals_nod_var2'][:]
    c3_all = 1 - c1_all - c2_all
    ds.close()

    mesh = base_mesh.copy(deep=True)

    # PyVista plotter (off-screen movie mode)
    plotter = pv.Plotter(off_screen=True)
    plotter.open_movie(output_filename, framerate=48)
    plotter.view_xy()

    camera_set = False
    for i in tqdm(range(nt), desc="Generating frames"):
        c1 = c1_all[i]
        c2 = c2_all[i]
        c3 = c3_all[i]

        # Map c1 to magenta, c2 to yellow, c3 to cyan
        # rgb_colors = np.clip(
        #     np.stack([
        #         c1 * 1 + c2 * 1 + c3 * 0,  # Red: c1 + c2
        #         c1 * 0 + c2 * 1 + c3 * 1,  # Green: c2 + c3
        #         c1 * 1 + c2 * 0 + c3 * 1   # Blue: c1 + c3
        #     ], axis=1), 0, 1
        # )
        rgb_colors = np.clip(np.stack([c1,c2,c3], axis=1), 0, 1)
        rgb_colors = np.clip(rgb_colors, 0, 1)
        hsv = mcolors.rgb_to_hsv(rgb_colors)
        hsv[:, 1] = 1  # set saturation to max
        rgb_colors = mcolors.hsv_to_rgb(hsv)
        mesh.point_data['rgb'] = (rgb_colors * 255).astype(np.uint8)

        plotter.add_mesh(mesh, scalars='rgb', rgb=True, show_scalar_bar=False)
        if not camera_set:
            plotter.view_xy()
            plotter.camera.Zoom(1.0)  # Optional: adjust zoom if needed
            camera_set = True

        plotter.write_frame()
        plotter.clear_actors()  # Only remove mesh, keep camera

    plotter.close()

def main_block0(exodus_filename, output_filename):
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

    # Get block 0 node indices
    elem_node_block0 = np.ma.getdata(ds.variables['connect1'][:]) - 1
    block0_nodes = np.unique(elem_node_block0.flatten())

    c1_all = ds.variables['vals_nod_var1'][:]
    c2_all = ds.variables['vals_nod_var2'][:]
    c3_all = 1 - c1_all - c2_all
    ds.close()

    mesh = base_mesh.copy(deep=True)

    plotter = pv.Plotter(off_screen=True)
    plotter.open_movie(output_filename, framerate=48)
    plotter.view_xy()

    camera_set = False
    n_nodes = mesh.n_points
    for i in tqdm(range(nt), desc="Generating frames (block0)"):
        c1 = np.zeros(n_nodes)
        c2 = np.zeros(n_nodes)
        c3 = np.zeros(n_nodes)
        c1[block0_nodes] = c1_all[i][block0_nodes]
        c2[block0_nodes] = c2_all[i][block0_nodes]
        c3[block0_nodes] = c3_all[i][block0_nodes]

        rgb_colors = np.clip(np.stack([c1, c2, c3], axis=1), 0, 1)
        hsv = mcolors.rgb_to_hsv(rgb_colors)
        hsv[:, 1] = 1
        rgb_colors = mcolors.hsv_to_rgb(hsv)
        mesh.point_data['rgb'] = (rgb_colors * 255).astype(np.uint8)

        plotter.add_mesh(mesh, scalars='rgb', rgb=True, show_scalar_bar=False)
        if not camera_set:
            plotter.view_xy()
            plotter.camera.Zoom(1.0)
            camera_set = True

        plotter.write_frame()
        plotter.clear_actors()

    plotter.close()

def main_block0_only(exodus_filename, output_filename):
    ds = netCDF4.Dataset(exodus_filename)
    times = ds.variables['time_whole'][:]
    nt = len(times)

    X = np.ma.getdata(ds.variables['coordx'][:])
    Y = np.ma.getdata(ds.variables['coordy'][:])
    Z = np.zeros_like(X)
    points = np.vstack([X, Y, Z]).T

    # Only use connectivity for block 0 (usually 'connect1')
    elem_node_block0 = np.ma.getdata(ds.variables['connect1'][:]) - 1
    mesh_block0 = pv.UnstructuredGrid({pv.CellType.QUAD: elem_node_block0}, points)

    c1_all = ds.variables['vals_nod_var1'][:]
    c2_all = ds.variables['vals_nod_var2'][:]
    c3_all = 1 - c1_all - c2_all
    ds.close()

    plotter = pv.Plotter(off_screen=True)
    plotter.open_movie(output_filename, framerate=48)
    plotter.view_xy()

    camera_set = False
    for i in tqdm(range(nt), desc="Generating frames (block0 only)"):
        c1 = c1_all[i]
        c2 = c2_all[i]
        c3 = c3_all[i]

        # Only assign data for block 0 nodes
        mesh_block0.point_data['rgb'] = (np.clip(np.stack([c1, c2, c3], axis=1), 0, 1) * 255).astype(np.uint8)

        plotter.add_mesh(mesh_block0, scalars='rgb', rgb=True, show_scalar_bar=False)
        if not camera_set:
            plotter.view_xy()
            plotter.camera.Zoom(1.0)
            camera_set = True

        plotter.write_frame()
        plotter.clear_actors()

    plotter.close()

# In your __main__ block, add:
if __name__ == "__main__":
    # Process all .e files in the results/output_dump_3p directory
    input_dir = "output"
    for exodus_file in sorted(os.listdir(input_dir)):
        if exodus_file.endswith('_test200_diff.e'):
            exodus_path = os.path.join(input_dir, exodus_file)
            # Create output filename by replacing .e with .gif
            output_filename = os.path.join(input_dir, exodus_file.replace('.e', '_v3_Mexp.mp4'))
            print(f"Processing {exodus_file} (all blocks)...")
            main(exodus_path, output_filename)
            print(f"Created {output_filename}")

            # # Block 0 only video
            # output_filename_block0 = os.path.join(input_dir, exodus_file.replace('.e', '_v3_Mexp_block0.mp4'))
            # print(f"Processing {exodus_file} (block 0 only)...")
            # main_block0(exodus_path, output_filename_block0)
            # print(f"Created {output_filename_block0}")

            # Block 0 only video (alternative method)
            output_filename_block0_only = os.path.join(input_dir, exodus_file.replace('.e', '_v3_Mexp_block0_only.mp4'))
            print(f"Processing {exodus_file} (block 0 only, alternative method)...")
            main_block0_only(exodus_path, output_filename_block0_only)
            print(f"Created {output_filename_block0_only}")
