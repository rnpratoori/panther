import os
import netCDF4
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter

def get_block_nodes(ds, block_idx):
    connect_var = f'connect{block_idx+1}'
    if connect_var in ds.variables:
        elem_node = np.ma.getdata(ds.variables[connect_var][:]) - 1
        return np.unique(elem_node.flatten())
    else:
        return None

def compute_averages(ds, c1_all, c2_all, cs_all, node_indices, times):
    nt = len(times)
    averages = []
    for i in range(nt):
        avg_c1 = np.mean(c1_all[i][node_indices])
        avg_c2 = np.mean(c2_all[i][node_indices])
        avg_cs = np.mean(cs_all[i][node_indices])
        averages.append([times[i], avg_c1, avg_c2, avg_cs])
    return pd.DataFrame(averages, columns=['time', 'avg_c1', 'avg_c2', 'avg_cs'])

def compute_and_save(exodus_filename, excel_filename):
    ds = netCDF4.Dataset(exodus_filename)
    times = ds.variables['time_whole'][:]
    c1_all = ds.variables['vals_nod_var1'][:]
    c2_all = ds.variables['vals_nod_var2'][:]
    cs_all = 1 - c1_all - c2_all

    X = np.ma.getdata(ds.variables['coordx'][:])
    Y = np.ma.getdata(ds.variables['coordy'][:])

    # Block 0
    block0_nodes = get_block_nodes(ds, 0)
    # Block 1
    # block1_nodes = get_block_nodes(ds, 1)
    # # All nodes
    # all_nodes = np.arange(c1_all.shape[1])

    ds.close()

    dfs = {}
    if block0_nodes is not None:
        dfs['block0'] = compute_averages(ds, c1_all, c2_all, cs_all, block0_nodes, times)
    # if block1_nodes is not None:
    #     dfs['block1'] = compute_averages(ds, c1_all, c2_all, cs_all, block1_nodes, times)
    # dfs['all_blocks'] = compute_averages(ds, c1_all, c2_all, cs_all, all_nodes, times)

    with pd.ExcelWriter(excel_filename, engine='openpyxl') as writer:
        for sheet, df in dfs.items():
            df.to_excel(writer, sheet_name=sheet, index=False)
    print(f"Saved averages to {excel_filename}")

    # Save plots for each case
    for sheet, df in dfs.items():
        # plot_averages(df, f"Averages vs Time ({sheet})", excel_filename.replace('.xlsx', f'_{sheet}.png'))
        animate_averages(df, f"Average volume fractions", excel_filename.replace('.xlsx', f'_{sheet}.mp4'))

    # Row-wise averages and animations
    unique_y, row_indices = get_row_indices(Y)
    row_cases = {
        'block0': block0_nodes,
        # 'block1': block1_nodes,
        # 'all_blocks': all_nodes
    }
    for case, node_mask in row_cases.items():
        if node_mask is None or len(node_mask) == 0:
            continue
        row_dfs = compute_row_averages(c1_all, c2_all, cs_all, row_indices, Y, node_mask=node_mask)
        # Plot for first and last timestep as example
        # plot_row_averages(row_dfs[0], f"Row averages at t=0 ({case})", excel_filename.replace('.xlsx', f'_{case}_row0.png'))
        # plot_row_averages(row_dfs[-1], f"Row averages at final t ({case})", excel_filename.replace('.xlsx', f'_{case}_rowfinal.png'))
        # Animate
        animate_row_averages(row_dfs, f"Depth profile", excel_filename.replace('.xlsx', f'_{case}_rowanim.mp4'))

def plot_averages(df, title, filename):
    plt.figure()
    plt.plot(df['time'], df['avg_c1'], label=r'$\bar{c}_1$', color='red')
    plt.plot(df['time'], df['avg_c2'], label=r'$\bar{c}_2$', color='green')
    plt.plot(df['time'], df['avg_cs'], label=r'$\bar{c}_s$', color='blue')
    plt.xlabel('Time')
    plt.ylabel('Average value')
    plt.title(title)
    plt.legend()
    plt.xlim(0, 1E-4)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def animate_averages(df, title, filename, fps=48, as_gif=False):
    fig, ax = plt.subplots()
    ax.set_xlim(0, 1E-4)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Time')
    ax.set_ylabel('Average value')
    ax.set_title(title)

    line1, = ax.plot([], [], 'r-', label=r'$\bar{c}_1$')
    line2, = ax.plot([], [], 'g-', label=r'$\bar{c}_2$')
    line3, = ax.plot([], [], 'b-', label=r'$\bar{c}_s$')
    ax.legend()

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        return line1, line2, line3

    def update(frame):
        x = df['time'][:frame+1]
        y1 = df['avg_c1'][:frame+1]
        y2 = df['avg_c2'][:frame+1]
        y3 = df['avg_cs'][:frame+1]
        line1.set_data(x, y1)
        line2.set_data(x, y2)
        line3.set_data(x, y3)
        return line1, line2, line3

    anim = FuncAnimation(fig, update, frames=len(df), init_func=init, blit=True, interval=1000/fps)
    if as_gif:
        anim.save(filename, writer=PillowWriter(fps=fps))
    else:
        anim.save(filename, writer=FFMpegWriter(fps=fps))
    plt.close(fig)
    print(f"Saved animation to {filename}")

def get_row_indices(Y, tol=1e-8):
    """Group node indices by unique Y coordinate (row). Returns a list of arrays, one per row."""
    unique_y = np.unique(Y)
    row_indices = []
    for y in unique_y:
        indices = np.where(np.abs(Y - y) < tol)[0]
        row_indices.append(indices)
    return unique_y, row_indices

def compute_row_averages(c1_all, c2_all, cs_all, row_indices, Y, node_mask=None):
    """
    Returns: list of DataFrames, one per timestep, with columns: y, avg_c1, avg_c2, avg_cs
    If node_mask is given, only use those nodes (e.g. for block 0/1).
    """
    nt = c1_all.shape[0]
    results = []
    for i in range(nt):
        row_avgs = []
        for indices in row_indices:
            if node_mask is not None:
                indices = np.intersect1d(indices, node_mask)
            if len(indices) == 0:
                row_avgs.append([np.nan, np.nan, np.nan, np.nan])
            else:
                y_val = np.mean(Y[indices])
                avg_c1 = np.mean(c1_all[i][indices])
                avg_c2 = np.mean(c2_all[i][indices])
                avg_cs = np.mean(cs_all[i][indices])
                row_avgs.append([y_val, avg_c1, avg_c2, avg_cs])
        df = pd.DataFrame(row_avgs, columns=['y', 'avg_c1', 'avg_c2', 'avg_cs'])
        results.append(df)
    return results

def plot_row_averages(df, title, filename):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    # Left subplot: c1 and c2
    ax1.plot(df['avg_c1'], df['y'], label=r'$\bar{c}_1$', color='red')
    ax1.plot(df['avg_c2'], df['y'], label=r'$\bar{c}_2$', color='green')
    ax1.set_xlabel('Average volume fraction (c1, c2)')
    ax1.set_ylabel('Y (height)')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.legend()
    ax1.set_title(r'$\bar{c}_1$ & $\bar{c}_2$')
    # Right subplot: cs
    ax2.plot(df['avg_cs'], df['y'], label=r'$\bar{c}_s$', color='blue')
    ax2.set_xlabel('Average volume fraction (cs)')
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.legend()
    ax2.set_title(r'$\bar{c}_s$')
    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def animate_row_averages(row_dfs, title, filename, fps=48, as_gif=False):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_xlabel('Average volume fraction (c1, c2)')
    ax1.set_ylabel('Y (height)')
    ax1.set_title(r'$\bar{c}_1$ & $\bar{c}_2$')
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.set_xlabel(r'Average volume fraction ($\bar{c}_s$)')
    ax2.set_title(r'$\bar{c}_s$')
    line1, = ax1.plot([], [], 'r-', label=r'$\bar{c}_1$')
    line2, = ax1.plot([], [], 'g-', label=r'$\bar{c}_2$')
    line3, = ax2.plot([], [], 'b-', label=r'$\bar{c}_s$')
    ax1.legend()
    ax2.legend()
    fig.suptitle(title)

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        return line1, line2, line3

    def update(frame):
        df = row_dfs[frame]
        y = df['y']
        line1.set_data(df['avg_c1'], y)
        line2.set_data(df['avg_c2'], y)
        line3.set_data(df['avg_cs'], y)
        return line1, line2, line3

    anim = FuncAnimation(fig, update, frames=len(row_dfs), init_func=init, blit=True, interval=1000/fps)
    if as_gif:
        anim.save(filename, writer=PillowWriter(fps=fps))
    else:
        anim.save(filename, writer=FFMpegWriter(fps=fps))
    plt.close(fig)
    print(f"Saved row animation to {filename}")

# Example usage after your Excel/plot saving loop:
# animate_averages(df, "Averages vs Time (block0)", "block0_anim.mp4")

if __name__ == "__main__":
    input_dir = "output"
    for exodus_file in sorted(os.listdir(input_dir)):
        if exodus_file.endswith('_test200_diff.e'):
            exodus_path = os.path.join(input_dir, exodus_file)
            excel_filename = os.path.join(input_dir, exodus_file.replace('.e', '_averages.xlsx'))
            print(f"Processing {exodus_file}...")
            compute_and_save(exodus_path, excel_filename)