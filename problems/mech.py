import pandas as pd
import matplotlib.pyplot as plt
import glob
import numpy as np

labels = [
    r'$1\times10^{-8}s$',
    r'$1\times10^{-7}s$',
    r'$1\times10^{-6}s$',
    r'$1\times10^{-5}s$',
    r'$1\times10^{-4}s$'
]

colors = ['red', 'green', 'blue', 'orange', 'purple']

# Process mech_ files
csv_files1 = sorted(glob.glob('mech_*.csv'))
areas1 = []
for idx, csv_file in enumerate(csv_files1):
    df = pd.read_csv(csv_file)
    stretch = 1 + df.iloc[:, 1] / 1.02
    stress = df.iloc[:, 2]
    area = np.trapz(stress, stretch)
    areas1.append(area)

# Process mech2_ files
csv_files2 = sorted(glob.glob('mech2_*.csv'))
areas2 = []
for idx, csv_file in enumerate(csv_files2):
    df = pd.read_csv(csv_file)
    stretch = 1 + df.iloc[:, 1] / 1.02
    stress = df.iloc[:, 2]
    area = np.trapz(stress, stretch)
    areas2.append(area)

# Plot both area curves in one plot
plt.figure(figsize=(7,5))
plt.plot(labels[:len(areas1)], areas1, marker='o', color='black', linestyle='-', label='case1')
plt.plot(labels[:len(areas2)], areas2, marker='s', color='blue', linestyle='--', label='case2')
plt.ylabel('Area under Stress-Stretch Curve')
plt.title('Strength vs Time')
plt.legend()
plt.tight_layout()
plt.savefig('mech_area_vs_timestep_compare.png', dpi=500)
plt.show()