import pandas as pd
import matplotlib.pyplot as plt
import glob
import numpy as np

# Define labels for the plots
labels = [
    r'$1\times10^{-6}s$',
    r'$1\times10^{-5}s$',
    r'$2\times10^{-5}s$',
    r'$3\times10^{-5}s$',
    r'$4\times10^{-5}s$',
    r'$5\times10^{-5}s$',
    r'$6\times10^{-5}s$',
    r'$7\times10^{-5}s$',
    r'$8\times10^{-5}s$',
    r'$9\times10^{-5}s$',
    r'$1\times10^{-4}s$'
]

# Define colors for the plots
colors = [
    'tab:blue',
    'tab:orange',
    'tab:green',
    'tab:red',
    'tab:purple',
    # Adding more in case they are needed
    'tab:brown',
    'tab:pink',
    'tab:gray',
    'tab:olive',
    'tab:cyan',
    'gold'
]

# --- NEW: Define linestyles and markers for B&W printing ---
linestyles = ['-', '--', '-.', ':', '-'] # Solid, dashed, dash-dot, dotted, solid
# markers = ['o', 's', 'v', '^', 'D'] # Circle, square, triangle_down, triangle_up, diamond

# --- Plot 1: Combined Stress vs. Stretch ---
# plt.figure(figsize=(8, 6)) # Create a figure for the stress-stretch curves

# Process mech files
csv_files1 = sorted(glob.glob('/home/rnp/MOOSE/projects/panther/problems/output/mech_void/mech*.csv'))
areas1 = []

for idx, csv_file in enumerate(csv_files1):
    # This check prevents errors if there are more csv files than defined styles
    # if idx >= len(linestyles):
    #     print(f"Warning: More files than linestyles defined. Skipping {csv_file}")
    #     break

    df = pd.read_csv(csv_file)
    stretch = 1 + df.iloc[:, 1] / 2.02
    stress = df.iloc[:, 2]
    c_stress = df.iloc[:, 2] / 1e8
    
    # --- MODIFIED: Added linestyle and marker arguments ---
    # Plot the stress vs. stretch for the current file
    # 'markevery' is used to place markers at intervals to avoid clutter
    # plt.plot(stretch, c_stress, 
    #          color=colors[idx], 
    #          linestyle=linestyles[idx], 
    #         #  marker=markers[idx],
    #          markevery=25, # Place a marker every 25 data points
    #          label=labels[idx])
    
    # Calculate the area under the curve for the second plot
    area = np.trapz(stress, stretch)
    areas1.append(area)

# # Configure and save the combined stress-stretch plot
# plt.xlabel('Stretch',fontsize=18)
# plt.ylabel(r'Cauchy Stress ($\times 10^8$)',fontsize=18)
# plt.title('Cauchy Stress vs. Stretch Curves',fontsize=20)
# plt.legend(fontsize=14)
# plt.xlim(1,1.5)
# plt.ylim(0,1.5)
# plt.yticks([0, 0.5, 1, 1.5],fontsize=16)
# plt.xticks(fontsize=16)
# plt.grid(True, linestyle='--', alpha=0.6) # A grid also helps with readability
# plt.savefig('stress_vs_stretch_combined.png', dpi=500)
# plt.show()


# --- Plot 2: Area Under Curve vs. Time Step (Original Plot) ---
plt.figure(figsize=(8, 6))
plt.plot(labels[:len(areas1)], areas1, marker='o', color='black', linestyle='-', label='case1')
plt.ylabel('Strength',fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.title('Strength vs. Time',fontsize=20)
# Rotate x-axis labels for better readability
plt.xticks([labels[0], labels[2], labels[5], labels[8], labels[10]], rotation=45, ha='right', fontsize=16)
plt.yticks(fontsize=16) 
# plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('mech_area_vs_timestep_compare.png', dpi=500)
plt.show()
