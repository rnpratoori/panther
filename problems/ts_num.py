"""Print timestep counts and time ranges for Exodus files in the current folder."""

from netCDF4 import Dataset
import glob

for f in sorted(glob.glob("*.e*")):
    try:
        ds = Dataset(f)
        times = ds.variables['time_whole'][:]
        print(f"{f}: {len(times)} timesteps, time range = [{times[0]:.4e}, {times[-1]:.4e}]")
        ds.close()
    except Exception as e:
        print(f"Skipping {f} ({e})")
