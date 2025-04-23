import pandas as pd
import matplotlib.pyplot as plt
import os

# Load the CSV file
csv_file = 'output/3p_dis_t4.csv'
data = pd.read_csv(csv_file)

# Create output filename by replacing .csv with .jpg
output_file = os.path.splitext(csv_file)[0] + '.jpg'

# Skip the second row by creating new dataframe without it
data = pd.concat([data.iloc[:1], data.iloc[2:]], ignore_index=True)

# Filter data between t=0 and t=1
data = data[(data[data.columns[0]] >= 0) & (data[data.columns[0]] <= 1)]

# Get column names from first row
columns = data.columns

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(data[columns[0]], data[columns[2]], label=columns[2])
plt.plot(data[columns[0]], data[columns[3]], label=columns[3])
# Calculate difference between columns 2 and 1 for third plot
local_energy = data[columns[3]] - data[columns[2]]
plt.plot(data[columns[0]], local_energy, label='local_energy')

# Add labels and formatting
plt.xlabel(columns[0])
plt.ylabel('Energy')
plt.title('Energy Evolution')
plt.legend()

# Set x-axis limits explicitly
plt.xlim(0, 1)

# Save the plot
plt.tight_layout()
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()