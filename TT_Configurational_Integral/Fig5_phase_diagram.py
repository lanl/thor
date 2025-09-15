#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Path to the saved transition data
data_file = './data_and_figures/Sn-phase_transition_data.txt'

# Load the data (skip header)
pressure, md_temp, ttci_temp = np.loadtxt(data_file, unpack=True, comments='#')

# Create the plot
plt.figure(figsize=(5, 4))
plt.plot(pressure, md_temp, '-r', label='MD')
plt.plot(pressure, ttci_temp, '--ob', label='TTCI (r=1)')

# Add annotations
plt.xlim([0, 2])
plt.ylim([0, 600])
plt.xlabel('Pressure (GPa)')
plt.ylabel('Temperature (K)')
plt.text(0.75, 200, 'DC', fontsize=14)
plt.text(1.50, 500, 'Beta', fontsize=14)
plt.legend(loc='lower right')
plt.grid(True)

# Save and show plot
output_path = './data_and_figures/phase_transition_plot.pdf'
plt.savefig(output_path, format='pdf', bbox_inches='tight')
plt.show()

print(f"Plot saved to: {output_path}")