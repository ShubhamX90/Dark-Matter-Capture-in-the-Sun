import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#########################
# Importing the CSV file
#########################

import os
this_dir, this_filename = os.path.split(__file__)

Vel_Dist_Path = os.path.join(this_dir, "SunDMVelDist.csv")

vel_dist_file = pd.read_csv(Vel_Dist_Path)

#########################
# Interpolating the data
#########################

interpolate_func = interp1d(vel_dist_file['Velocity_Range'], vel_dist_file['VelocityDist_Planet_Frame'], kind='cubic')

############################
# Finding the normalization
############################

integrand = interpolate_func(vel_dist_file['Velocity_Range']) * (vel_dist_file['Velocity_Range'] ** 2)
integral = np.trapz(integrand, vel_dist_file['Velocity_Range'])
normalized_distribution = interpolate_func(vel_dist_file['Velocity_Range']) / integral

print("Normalization factor:", integral)


###############################
# The plot after normalization
###############################

plt.figure(figsize=(12, 8))
plt.plot(vel_dist_file['Velocity_Range'], normalized_distribution, label='Normalized Distribution')
plt.xlabel('Velocity Range')
plt.ylabel('Normalized Distribution')
plt.title('Normalized Velocity Distribution')
plt.legend()
plt.grid(True,  which='both')
plt.gca().ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.show()

###################################
# Printing the normalized vel dist
###################################

print("Normalized distribution:")
print("Velocity Range\tNormalized Distribution")
for Velocity_Range, VelocityDist_Planet_Frame in zip(vel_dist_file['Velocity_Range'], normalized_distribution):
    print(f"{Velocity_Range}\t{VelocityDist_Planet_Frame}")

# Create a new DataFrame for the normalized data
normalized_data = pd.DataFrame({'Velocity_Range': vel_dist_file['Velocity_Range'],
                                'Normalized_Distribution': normalized_distribution})

# Save the DataFrame to a new CSV file
output_file_path = os.path.join(this_dir, "Normalized_Velocity_Distribution.csv")
normalized_data.to_csv(output_file_path, index=False)

print(f"Normalized velocity distribution data saved to: {output_file_path}")