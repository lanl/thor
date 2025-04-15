import numpy as np
import matplotlib.pyplot as plt

# Load data
txt_path = './data_and_figures/Cu_FCC_smatb_energy_pressure.txt'
data = np.loadtxt(txt_path, comments='#')

densities = np.unique(data[:, 0])
colors = ['k', 'r', 'b']
save_path = './data_and_figures/Fig_Cu_FCC_smatb_energy_pressure.pdf'

plt.figure(figsize=(8, 5), dpi=300)

# Energy plot
plt.subplot(1, 2, 1)
for i, d in enumerate(densities):
    idx = np.where(data[:, 0] == d)
    T = data[idx, 1].squeeze()
    E_MD = data[idx, 2].squeeze()
    E_r1 = data[idx, 3].squeeze()
    E_r2 = data[idx, 4].squeeze()
    
    plt.plot(T, E_r1, 's'+colors[i], markerfacecolor='none', label=f'TT-r1-{d:.2f}')
    plt.plot(T, E_r2, '*'+colors[i], markerfacecolor='none', label=f'TT-r2-{d:.2f}')
    plt.plot(T, E_MD, '--'+colors[i], label=f'MD-{d:.2f}')

plt.xlabel('Temperature (K)')
plt.ylabel('Energy (eV)')
plt.legend(fontsize=8)
plt.grid()

# Pressure plot
plt.subplot(1, 2, 2)
for i, d in enumerate(densities):
    idx = np.where(data[:, 0] == d)
    T = data[idx, 1].squeeze()
    P_MD = data[idx, 5].squeeze()
    P_r1 = data[idx, 6].squeeze()
    P_r2 = data[idx, 7].squeeze()

    plt.plot(T, P_r1, 's'+colors[i], markerfacecolor='none', label=f'TT-r1-{d:.2f}')
    plt.plot(T, P_r2, '*'+colors[i], markerfacecolor='none', label=f'TT-r2-{d:.2f}')
    plt.plot(T, P_MD, '--'+colors[i], label=f'MD-{d:.2f}')

plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (GPa)')
plt.grid()

plt.suptitle('Cu - FCC - smatb')
plt.tight_layout()
plt.savefig(save_path, format='pdf', bbox_inches='tight')
plt.show()
print(f"Plot saved to: {save_path}")