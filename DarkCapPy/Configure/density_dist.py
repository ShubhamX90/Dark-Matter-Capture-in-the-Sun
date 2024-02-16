import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Read the CSV file
import os
this_dir, this_filename = os.path.split(__file__)

planet_path = os.path.join(this_dir, "struct_b16_agss09.csv")

planet_file = pd.read_csv(planet_path)

interpolate_func = interp1d(planet_file['Radius'], planet_file['Rho'], kind='cubic')

# Normalize the data

integral = np.trapz(interpolate_func(planet_file['Radius']), planet_file['Radius'])
normalized_density = interpolate_func(planet_file['Radius']) / integral

print("Normalization factor : ", integral)

# Plot the density distribution
plt.figure(figsize=(10, 6))
plt.plot(planet_file['Radius'], normalized_density)
plt.title('Density Distribution in the Sun')
plt.xlabel('Radius')
plt.ylabel('Normalized Density')
plt.grid(True)
plt.show()

# Print the normalized data
print("Normalized Data:")
print(planet_file)