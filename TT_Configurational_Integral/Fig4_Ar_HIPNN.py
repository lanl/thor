import numpy as np
import matplotlib.pyplot as plt

# Load data
txt_path = './data_and_figures/Ar_FCC_HIPNN_energy_pressure.txt'
data = np.loadtxt(txt_path, comments='#')

# Extract columns
densities = np.unique(data[:, 0])
T = data[:, 1]
E_MD = data[:, 2]
E_r1 = data[:, 3]
E_r2 = data[:, 4]
P_MD = data[:, 5]
P_r1 = data[:, 6]
P_r2 = data[:, 7]

colors = ['k', 'r', 'b']
save_path = './data_and_figures/Fig_Ar_FCC_HIPNN_energy_pressure_from_txt.pdf'

plt.figure(figsize=(8, 5), dpi=300)

# --- Energy Plot ---
plt.subplot(1, 2, 1)
for i, d in enumerate(densities):
    idx = data[:, 0] == d
    plt.plot(T[idx], E_r1[idx], 's'+colors[i], markerfacecolor='none', markersize=6, label=f'TT-r1-{d:.2f}')
    plt.plot(T[idx], E_r2[idx], '*'+colors[i], markerfacecolor='none', markersize=6, label=f'TT-r2-{d:.2f}')
    plt.plot(T[idx], E_MD[idx], '--'+colors[i], label=f'MD-{d:.2f}')
plt.xlabel('Temperature (K)')
plt.ylabel('Energy (eV)')
plt.grid(True)
plt.legend(fontsize=8)

# --- Pressure Plot ---
plt.subplot(1, 2, 2)
for i, d in enumerate(densities):
    idx = data[:, 0] == d
    plt.plot(T[idx], P_r1[idx], 's'+colors[i], markerfacecolor='none', markersize=6, label=f'TT-r1-{d:.2f}')
    plt.plot(T[idx], P_r2[idx], '*'+colors[i], markerfacecolor='none', markersize=6, label=f'TT-r2-{d:.2f}')
    plt.plot(T[idx], P_MD[idx], '--'+colors[i], label=f'MD-{d:.2f}')
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (GPa)')
plt.grid(True)
plt.legend(fontsize=8)

# Layout
plt.suptitle('Ar - FCC - HIPNN')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(save_path, format='pdf', bbox_inches='tight')
plt.show()

print(f"Plot saved to: {save_path}")