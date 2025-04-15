import os
import shutil
import netCDF4
import numpy as np

# Append a new nodal variable to an Exodus file
def overwrite_nodal_variable(exo_file, var_name_to_overwrite, new_var_name, values):
    with netCDF4.Dataset(exo_file, 'r+') as ds:
        # Get number of nodes and timesteps
        nnodes = ds.dimensions['num_nodes'].size
        nsteps = ds.dimensions['time_step'].size

        if len(values) != nnodes:
            raise ValueError("Length of values must match number of nodes")

        # Read nodal variable names
        node_var_names = []
        for name in ds.variables['name_nod_var'][:]:
            name_str = ''.join(x.decode() if isinstance(x, bytes) else ' ' for x in name)
            node_var_names.append(name_str.strip())

        if var_name_to_overwrite not in node_var_names:
            raise ValueError(f"Variable '{var_name_to_overwrite}' not found")

        var_index = node_var_names.index(var_name_to_overwrite)

        # Update the variable name
        name_var = ds.variables['name_nod_var']
        name_length = name_var.shape[1]
        padded_name = np.zeros(name_length, dtype='S1')

        for i, char in enumerate(new_var_name[:name_length]):
            padded_name[i] = char.encode('ascii')

        name_var[var_index, :] = padded_name

        # Overwrite values at the last timestep
        var = ds.variables[f'vals_nod_var{var_index + 1}']
        var[nsteps - 1, :] = values

def get_min_max(exodus_file, variable='c'):
    model = netCDF4.Dataset(exodus_file)
    # Handle masked array properly
    node_var_names = []
    for name in model.variables['name_nod_var'][:]:
        # Convert masked array to string, handling masked values
        name_str = ''.join(x.decode() if isinstance(x, bytes) else ' ' for x in name)
        node_var_names.append(name_str.strip())
    
    var_index = node_var_names.index(variable) + 1
    last_step = model.variables['time_whole'].shape[0] - 1
    var_values = model.variables[f'vals_nod_var{var_index}'][:][last_step]

    model.close()
    return np.min(var_values), np.max(var_values)

def min_max_scale(values, target_min, target_max):
    """Scale values to target range"""
    original_min = np.min(values)
    original_max = np.max(values)
    scaled = (values - original_min) / (original_max - original_min)
    return scaled * (target_max - target_min) + target_min

def get_values(exodus_file, variable='c'):
    """Get values of a variable at the last timestep"""
    model = netCDF4.Dataset(exodus_file)
    
    # Handle masked array properly
    node_var_names = []
    for name in model.variables['name_nod_var'][:]:
        # Convert masked array to string, handling masked values
        name_str = ''.join(x.decode() if isinstance(x, bytes) else ' ' for x in name)
        node_var_names.append(name_str.strip())
    
    var_index = node_var_names.index(variable) + 1
    last_step = model.variables['time_whole'].shape[0] - 1
    var_values = model.variables[f'vals_nod_var{var_index}'][last_step, :].copy()
    
    model.close()
    return var_values

# === MAIN EXECUTION ===
if __name__ == '__main__':
    source_file = "output/2phase.e"
    copy_file = "output/2phase_copy.e"
    reference_file = "output/3phase_0.3_0.3.e"
    
    # Create a copy of the source file
    if os.path.exists(copy_file):
        os.remove(copy_file)
    shutil.copy2(source_file, copy_file)
    
    # First variable: c1_rescale
    target_min, target_max = get_min_max(reference_file, 'c1')
    print(f"c1 target range: [{target_min:.6f}, {target_max:.6f}]")
    
    source_values = get_values(source_file, 'c')
    source_min, source_max = np.min(source_values), np.max(source_values)
    print(f"c source range: [{source_min:.6f}, {source_max:.6f}]")
    
    scaled_values = min_max_scale(source_values, target_min, target_max)
    print(f"c1_rescale range: [{np.min(scaled_values):.6f}, {np.max(scaled_values):.6f}]")
    
    overwrite_nodal_variable(copy_file, 'c', 'c1_rescale', scaled_values)
    
    # Second variable: c2_rescale
    complement_values = 1 - source_values
    complement_min, complement_max = np.min(complement_values), np.max(complement_values)
    print(f"\n1-c range: [{complement_min:.6f}, {complement_max:.6f}]")
    
    # Scale complement_values (1-c) to match target range
    target_min_c2, target_max_c2 = get_min_max(reference_file, 'c2')
    scaled_complement = min_max_scale(complement_values, target_min_c2, target_max_c2)
    print(f"c2_rescale range: [{np.min(scaled_complement):.6f}, {np.max(scaled_complement):.6f}]")
    
    overwrite_nodal_variable(copy_file, 'c2', 'c2_rescale', scaled_complement)
    print(f"\nCreated {os.path.basename(copy_file)} with modified variables 'c1_rescale' and 'c2_rescale'")
