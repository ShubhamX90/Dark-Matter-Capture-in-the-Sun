import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Read the CSV file
import os
this_dir, this_filename = os.path.split(__file__)

planet_path = os.path.join(this_dir, "escVel_vs_radius.csv")

planet_file = pd.read_csv(planet_path)

interpolate_func = interp1d(planet_file['Radius'], planet_file['Escape_Velocity'], kind='linear')

# Normalize the data

integral = np.trapz(interpolate_func(planet_file['Radius']), planet_file['Radius'])
normalized_dist = interpolate_func(planet_file['Radius']) / integral

print("Normalization factor : ", integral)

# Plot the escape velocity distribution
plt.figure(figsize=(10, 6))
plt.plot(planet_file['Radius'], normalized_dist)
plt.title('Esc_Vel Distribution in the Sun')
plt.xlabel('Radius')
plt.ylabel('Normalized Escape Velocity Distribution')
plt.grid(True)
plt.show()

# Print the normalized data
print("Normalized Data:")
print(planet_file)