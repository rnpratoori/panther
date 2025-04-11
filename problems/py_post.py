import os
from os import listdir
from os.path import isfile, join
from pathlib import Path
from tqdm import tqdm
import netCDF4
import numpy as np
import pyvista as pv
import vtk
from pyvirtualdisplay import Display

# # Start virtual display using pyvirtualdisplay
# display = Display(visible=0, size=(800, 600))
# display.start()

# Ensure off-screen rendering for PyVista
os.environ["PYVISTA_VTK_OFFSCREEN"] = "true"
# pv.start_xvfb()

# Function to convert Exodus file to PyVista UnstructuredGrid for a given timestep
def exodus2PyVista(filename, nstep=-1):  # Default is the last timestep (nstep=-1)
    # Open the Exodus file using netCDF4
    model = netCDF4.Dataset(filename)

    # Check if 'pvf' exists as a node variable
    node_var_names = getNames(model, 'name_nod_var')
    if 'pvf' not in node_var_names:
        print("Node variable 'pvf' does not exist in the Exodus file.")
        model.close()
        return None
    
    # Read the coordinates
    X_all = np.ma.getdata(model.variables['coordx'][:])
    Y_all = np.ma.getdata(model.variables['coordy'][:])
    Z_all = np.zeros_like(X_all)  # Z is set to zero for a 2D grid

    # Combine coordinates into a points array
    points = np.vstack([X_all, Y_all, Z_all]).T

    # Get the element connectivity (node mapping)
    elem_node = np.ma.getdata(model.variables['connect1'][:]) - 1

    # Create a PyVista UnstructuredGrid with the mesh data
    grid = pv.UnstructuredGrid({vtk.VTK_QUAD: elem_node}, points)

    # Set to the last timestep if nstep=-1
    if nstep == -1:
        nstep = len(model.variables['time_whole']) - 1

    # Add pvf nodal data to the PyVista mesh
    pvf_data = model.variables['vals_nod_var{}'.format(node_var_names.index('pvf') + 1)][:][nstep]
    grid['pvf'] = pvf_data

    # Close the Exodus model
    model.close()

    return grid

# Function to get names of the nodal/element variables from the Exodus file
def getNames(model, key='name_nod_var'):
    # Get names of variables from the Exodus file catalog
    name_var = []
    for vname in np.ma.getdata(model.variables[key][:]).astype('U8'):
        name_var.append(''.join(vname))
    return name_var

# Function to save the last timestep as an image
def save_image(exodus_file, image_filename):
    # Convert Exodus file to PyVista grid for the last timestep
    grid = exodus2PyVista(exodus_file)

    # Create a PyVista plotter
    plotter = pv.Plotter(off_screen=True)

    # Add the grid to the plotter
    plotter.add_mesh(grid, scalars='pvf', clim=(0,1), cmap='Greys', scalar_bar_args={'title':'PVF', 'vertical':'True'})

    # Set up the camera and axes for 2D
    plotter.view_xy()
    # plotter.show_bounds(grid='front')
    plotter.show_axes()

    # Save the plot as an image
    plotter.screenshot(image_filename)

    # Close the plotter to free-up memory
    plotter.close()

# Example usage
curpath=Path(__file__).parent.absolute()
dirname=os.path.basename(curpath)
output_path=Path(join(curpath,'output'))
save_path=Path("/mnt/c/Users/rnp/Documents/Research/")

subdirs = [d for d in output_path.iterdir() if d.is_dir()]
for subdir in tqdm(subdirs):
    filelist=[f for f in listdir(subdir) if isfile(join(subdir, f)) and f.endswith('_dbc.e')]
    exodus_file = join(subdir, filelist[0])  # Replace with your Exodus file path
    image_filename = join(save_path, filelist[0].replace('.e', '') + '.png')  # Image file to save the output
    save_image(exodus_file, image_filename)

# display.stop()
