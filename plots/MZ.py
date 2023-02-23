import matplotlib.pyplot as plt
import numpy as np

# Read data from file
data = np.loadtxt('results/observables/MZ.dat')

plt.rc('font', family='Times')
plt.rc('mathtext', fontset='cm')

# Extract columns from data
TEMP = data[:, 0]
H = data[:, 1]
p = data[:, 2]
MZ = data[:, 3]

# Set up figure and subplots
fig, ax = plt.subplots(1, 1, figsize=(5, 5))

ax.tick_params(direction='in', top=True, right=True)

# Plot
plt.title(f"$\\Gamma={H[0]}$")
ax.plot(TEMP, MZ, linestyle='-', marker='o', color='red')
ax.set_xlabel('$T$', fontfamily='Times')
ax.set_ylabel('$\\langle |\\hat{\\mathrm{m}}_z| \\rangle$', fontfamily='Times')
ax.set_xlim(0, 6.00)
ax.set_ylim(0, 1)

plt.margins(0, 0)

# Show plot
plt.savefig('plots/fig1.pdf')
