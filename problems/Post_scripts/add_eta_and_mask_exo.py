#!/usr/bin/env python3
"""
add_eta_and_mask_exo.py

Reads an Exodus (netCDF) file, creates an `eta` nodal variable that is +1
inside specified circles and -1 elsewhere, zeroes nodal vars 'c' and 'c2'
where eta > threshold, and writes a new Exodus file with the new/modified
variables.

Usage:
  python add_eta_and_mask_exo.py \
    --in input.e --out output_with_eta.e \
    --circles "0.25,0.25,0.1 0.75,0.6,0.08" \
    --eta-name eta --threshold 0.0

Dependencies:
  pip install netCDF4 numpy
"""

import argparse
import os
import netCDF4
import numpy as np
import shutil
import sys

def decode_name_row(row):
    # row is an array of bytes like [b'a', b'b', b'c', ...]
    try:
        return b"".join(row).decode("ascii").strip()
    except Exception:
        # fallback
        return "".join([ch.decode("ascii") if isinstance(ch, bytes) else str(ch) for ch in row]).strip()

def read_node_var_names(ds):
    name_var = ds.variables['name_nod_var'][:]
    names = [decode_name_row(name_var[i, :]) for i in range(name_var.shape[0])]
    return names

def read_coords(ds):
    # Typical Exodus uses 'coordx', 'coordy', 'coordz'
    if 'coordx' in ds.variables:
        x = np.array(ds.variables['coordx'][:]).reshape(-1)
        y = np.array(ds.variables['coordy'][:]).reshape(-1) if 'coordy' in ds.variables else np.zeros_like(x)
        if 'coordz' in ds.variables:
            z = np.array(ds.variables['coordz'][:]).reshape(-1)
        else:
            z = np.zeros_like(x)
        pts = np.vstack([x, y, z]).T
        return pts
    # fallback: try 'coord' or error
    for vn in ds.variables:
        if vn.lower().startswith('coord') and ds.variables[vn].ndim == 1:
            # not robust; prefer explicit names above
            pass
    raise RuntimeError("Could not find coordinate variables 'coordx','coordy' in Exodus file.")

def build_eta_from_circles(points, circles, inside_val=1.0, outside_val=-1.0):
    xy = points[:, :2]
    eta = np.full(points.shape[0], outside_val, dtype=float)
    if not circles:
        return eta
    for (cx, cy, r) in circles:
        dx = xy[:,0] - cx
        dy = xy[:,1] - cy
        mask = (dx*dx + dy*dy) <= (r*r)
        eta[mask] = inside_val
    return eta

def parse_circles(s):
    circles = []
    s = s.strip()
    if s == "":
        return circles
    for token in s.split():
        parts = token.split(",")
        if len(parts) != 3:
            raise ValueError(f"Bad circle token '{token}'. Expect 'x,y,r'")
        circles.append(tuple(float(p) for p in parts))
    return circles

def create_copy_with_extra_nodal_var(srcfile, dstfile, new_var_name, new_values, modified_c=None, modified_c2=None):
    """
    Create a new NetCDF+Exodus file copying srcfile but with:
      - an extra nodal variable named new_var_name whose values at the last timestep
        are taken from new_values (others set to 0 for earlier times)
      - optionally replace last-timestep nodal arrays for 'c' and 'c2' with
        modified_c/modified_c2 (if provided)
    """
    with netCDF4.Dataset(srcfile, 'r') as src:
        # Read old name list and sizes
        old_name_var = src.variables['name_nod_var']
        old_nnodvar = old_name_var.shape[0]
        name_strlen = old_name_var.shape[1]
        # find an existing vals_nod_var to figure out time/node dims and dtype
        vals_keys = [k for k in src.variables.keys() if k.startswith('vals_nod_var')]
        if len(vals_keys) == 0:
            raise RuntimeError("No 'vals_nod_var#' variables found in source Exodus.")
        sample_vals = src.variables[vals_keys[0]]
        vals_dims = sample_vals.dimensions  # e.g. ('time_step','num_nodes')
        time_dim_name = vals_dims[0]
        node_dim_name = vals_dims[1]
        ntime = src.dimensions[time_dim_name].size
        nnodes = src.dimensions[node_dim_name].size

        # create dst file and copy dimensions (but increase num_nod_var by 1)
        with netCDF4.Dataset(dstfile, 'w') as dst:
            # copy dimensions except num_nod_var treated specially
            for dname, dim in src.dimensions.items():
                if dname == old_name_var.dimensions[0]:  # typically 'num_nod_var'
                    newsize = dim.size + 1 if dim.size is not None else None
                    dst.createDimension(dname, newsize)
                else:
                    dst.createDimension(dname, None if dim.isunlimited() else len(dim))
            # copy global attributes
            for attr in src.ncattrs():
                dst.setncattr(attr, src.getncattr(attr))

            # copy variables except name_nod_var (we will create new one) and exclude nothing else
            for vname, var in src.variables.items():
                if vname == 'name_nod_var':
                    continue
                # create variable in dst with same dtype and dims
                dst_var = dst.createVariable(vname, var.datatype, var.dimensions)
                # copy variable attributes
                for a in var.ncattrs():
                    try:
                        dst_var.setncattr(a, var.getncattr(a))
                    except Exception:
                        # ignore weird attributes that can't be copied
                        pass
                # copy data
                dst_var[:] = var[:]

            # build new name_nod_var array
            new_nnodvar = old_nnodvar + 1
            # Create empty array of S1 bytes
            new_name_arr = np.empty((new_nnodvar, name_strlen), dtype='S1')
            # copy old bytes
            old_raw = src.variables['name_nod_var'][:]
            new_name_arr[:old_nnodvar, :] = old_raw[:]
            # pad new var name
            padded = np.array(list(new_var_name.ljust(name_strlen)[:name_strlen].encode('ascii')), dtype='S1')
            new_name_arr[old_nnodvar, :] = padded
            # create and write name_nod_var in dst
            name_dims = src.variables['name_nod_var'].dimensions
            name_var = dst.createVariable('name_nod_var', src.variables['name_nod_var'].datatype, name_dims)
            name_var[:] = new_name_arr

            # create the new vals_nod_var{new_index}
            new_index = new_nnodvar  # 1-based index (vals_nod_var1 ... vals_nod_varN)
            new_vals_name = f'vals_nod_var{new_index}'
            # dtype and dims same as sample_vals
            new_vals_var = dst.createVariable(new_vals_name, sample_vals.datatype, sample_vals.dimensions)
            # copy any sample var attributes if desired (units etc) - none strictly needed
            for a in sample_vals.ncattrs():
                try:
                    new_vals_var.setncattr(a, sample_vals.getncattr(a))
                except Exception:
                    pass
            # initialize to zero
            new_vals_var[:] = np.zeros(new_vals_var.shape, dtype=new_vals_var.datatype)

            # fill the last timestep for the new var
            last_step = ntime - 1
            if new_values.shape[0] != nnodes:
                raise ValueError("new_values length does not match number of nodes in mesh")
            new_vals_var[last_step, :] = new_values

            # If user requested overwriting c and/or c2 at last timestep, update those copied variables
            # Find indices for 'c' and 'c2' in the old name list
            old_names = [decode_name_row(old_raw[i, :]) for i in range(old_nnodvar)]
            def overwrite_if_present(varname, new_arr):
                if varname in old_names:
                    vidx = old_names.index(varname) + 1
                    vname = f'vals_nod_var{vidx}'
                    if new_arr is None:
                        return
                    if dst.variables[vname].shape[0] != ntime or dst.variables[vname].shape[1] != nnodes:
                        raise RuntimeError(f"Unexpected shape for {vname}")
                    dst.variables[vname][last_step, :] = new_arr
                else:
                    # If 'c' not present, warn but continue
                    print(f"Warning: variable '{varname}' not found in source file; skipping overwrite", file=sys.stderr)

            overwrite_if_present('c', modified_c)
            overwrite_if_present('c2', modified_c2)

            # done
    return

def main():
    p = argparse.ArgumentParser(description="Add eta nodal var to Exodus and zero c/c2 where eta>threshold")
    p.add_argument("--in", dest="infile", required=True)
    p.add_argument("--out", dest="outfile", required=True)
    p.add_argument("--circles", required=True,
                   help="quoted string: 'x1,y1,r1 x2,y2,r2 ...' (units same as mesh coords)")
    p.add_argument("--eta-name", default="eta")
    p.add_argument("--threshold", type=float, default=0.0,
                   help="Zero c/c2 where eta > threshold (default 0.0)")
    args = p.parse_args()

    src = args.infile
    dst = args.outfile
    circles = parse_circles(args.circles)

    # read coords and existing c/c2
    with netCDF4.Dataset(src, 'r') as ds:
        pts = read_coords(ds)
        nnodes = pts.shape[0]

        node_var_names = read_node_var_names(ds)
        # find c,c2 (if absent raise/print)
        if 'c' not in node_var_names:
            print("Error: 'c' not found in source nodal variables", file=sys.stderr)
            # proceed: maybe user wants to create eta only; allow continuing but set c array to None
            c_vals = None
        else:
            idx = node_var_names.index('c') + 1
            varname = f'vals_nod_var{idx}'
            last_step = ds.variables[ds.variables.keys().__iter__().__next__()].shape[0] - 1  # placeholder, we'll compute properly below
            # find time dim via that var's dimensions
            sample_vals = ds.variables[f'vals_nod_var{1}']
            time_dim = sample_vals.dimensions[0]
            last_step = ds.dimensions[time_dim].size - 1
            c_vals = np.array(ds.variables[varname][last_step, :]).copy()

        if 'c2' not in node_var_names:
            c2_vals = None
        else:
            idx2 = node_var_names.index('c2') + 1
            varname2 = f'vals_nod_var{idx2}'
            sample_vals = ds.variables[f'vals_nod_var{1}']
            time_dim = sample_vals.dimensions[0]
            last_step = ds.dimensions[time_dim].size - 1
            c2_vals = np.array(ds.variables[varname2][last_step, :]).copy()

    # build eta
    eta = build_eta_from_circles(pts, circles, inside_val=1.0, outside_val=-1.0)
    # mask c and c2 where eta > threshold
    mask = eta > args.threshold
    if c_vals is not None:
        c_mod = c_vals.copy()
        c_mod[mask] = 0.0
    else:
        c_mod = None
    if c2_vals is not None:
        c2_mod = c2_vals.copy()
        c2_mod[mask] = 0.0
    else:
        c2_mod = None

    # create new file with extra nodal var and modified c/c2 at last timestep
    if os.path.exists(dst):
        os.remove(dst)
    create_copy_with_extra_nodal_var(src, dst, args.eta_name, eta, modified_c=c_mod, modified_c2=c2_mod)
    print(f"Wrote new Exodus file: {dst}")
    print(f"Added nodal variable '{args.eta_name}' and zeroed 'c'/'c2' where eta > {args.threshold}")

if __name__ == "__main__":
    main()
