import os
import netCDF4
import numpy as np

# Ensure off-screen rendering for PyVista
os.environ["PYVISTA_VTK_OFFSCREEN"] = "true"

# Get names of variables in the Exodus file
def getNames(model, key='name_nod_var'):
    name_var = []
    for vname in np.ma.getdata(model.variables[key][:]).astype('U8'):
        name_var.append(''.join(vname))
    return name_var

# Get min and max value of a variable at last timestep
def get_min_max(exodus_file, variable='c'):
    model = netCDF4.Dataset(exodus_file)
    node_var_names = getNames(model, 'name_nod_var')

    if variable not in node_var_names:
        print(f"Variable '{variable}' not found in {exodus_file}")
        model.close()
        return None, None

    var_index = node_var_names.index(variable) + 1
    last_step = model.variables['time_whole'].shape[0] - 1
    var_values = model.variables[f'vals_nod_var{var_index}'][:][last_step]

    model.close()
    return np.min(var_values), np.max(var_values)

# === MAIN EXECUTION ===
if __name__ == '__main__':
    exodus_file = "output/3phase_0.3_0.3.e"
    variable = 'c2'

    # Print min and max values
    min_val, max_val = get_min_max(exodus_file, variable=variable)
    print(f"{os.path.basename(exodus_file)} --> Min: {min_val:.6f}, Max: {max_val:.6f}")
