import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Read the CSV file
import os
this_dir, this_filename = os.path.split(__file__)

planet_path = os.path.join(this_dir, "struct_b16_agss09.csv")

planet_file = pd.read_csv(planet_path, delim_whitespace=True, header = 8)

interpolate_func = interp1d(planet_file['Radius'], planet_file['Rho'], kind='cubic')

# Normalize the data
integrand = interpolate_func(planet_file['Radius']) * (planet_file['Radius'] ** 2)
integral = np.trapz(integrand, planet_file['Radius'])
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
print("Velocity Range\tNormalized Distribution")
for Radius, Rho in zip(planet_file['Radius'], normalized_density):
    print(f"{Radius}\t{Rho}")

# Create a new DataFrame for the normalized data
normalized_data = pd.DataFrame({'Radius': planet_file['Radius'], 'Normalized_Density': normalized_density})

# Save the DataFrame to a new CSV file
output_file_path = os.path.join(this_dir, "Normalised_density_dist.csv")
normalized_data.to_csv(output_file_path, index=False)

print(f"Normalized data saved to: {output_file_path}")