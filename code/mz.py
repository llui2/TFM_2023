import matplotlib.pyplot as plt
import numpy as np

plt.rc('font', family='Times')
plt.rc('mathtext', fontset='cm')

# Read data from file
data = np.loadtxt('results/observables/MZ.dat')

# Set up arrays
TEMP = data[:, 0]
H = data[:, 1]
p = data[:, 2]
MZ = data[:, 3]

# Check unique values of H
H_values = np.unique(H)

for i in range(len(H_values)):

       plot_TEMP = []
       plot_MZ = []

       for j in range(len(H)):
              if H[j] == H_values[i]:
                     plot_TEMP.append(TEMP[j])
                     plot_MZ.append(MZ[j])

       # Set up figure and subplots
       fig, ax = plt.subplots(1, 1, figsize=(5, 5))

       ax.tick_params(direction='in', top=True, right=True)

       # Plot
       plt.title(f"$\\Gamma={H_values[i]}$")
       ax.plot(plot_TEMP, plot_MZ, linestyle='-', marker='o', color='red')
       ax.set_xlabel('$T$', fontfamily='Times')
       ax.set_ylabel('$\\langle |\\hat{\\mathrm{m}}_z| \\rangle$', fontfamily='Times')
       ax.set_xlim(0, 6.00)
       ax.set_ylim(0, 1)

       plt.margins(0, 0)

       # Show plot
       plt.savefig(f'results/fig{i+1}.pdf')
