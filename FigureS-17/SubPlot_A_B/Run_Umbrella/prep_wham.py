import numpy as np

# ./wham 0.8 6.2 250 1e-6 300 0 metadata.txt free_energy_distance.txt 1000 12345
# Constants and file name
kappa = 2  # kcal / nm
T = 300  # K
file = 'separation.txt'
N = 1000  # Define this as the number of points per window

# Load the data
data = np.loadtxt(file)

# Initialize lists for data separation
x = []
c = []

# Separate the data into windows and store the values
for index in range(0, len(data), N):
    window = data[index:index + N]
    if len(window) > 0:
        x.append(window[0,1])  # Window minimum
        c.append(window[:, 2])  # Separation distance

# Write the time series files and metadata file for WHAM
metadata_lines = []

for i, (x_val, c_val) in enumerate(zip(x, c)):
    # Define the time series file name for this window
    timeseries_filename = f'timeseries_window_{i + 1}.txt'
    np.savetxt(timeseries_filename, np.column_stack((np.arange(len(c_val)), c_val)), fmt='%d %.6f')
    
    # Prepare the line for the metadata file
    metadata_lines.append(f"{timeseries_filename} {x_val} {kappa}")

# Write the metadata file
with open('metadata.txt', 'w') as meta_file:
    meta_file.write('\n'.join(metadata_lines))

print("Time series files and metadata file have been created.")
