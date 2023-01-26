import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

# Load the data from the XY plot into lists for the X and Y values
time = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
rmsd = [0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]

# Use the gaussian_kde function from scipy.stats to estimate the PDF of the data
density = gaussian_kde(rmsd)

# Generate a range of values for the X axis
xs = np.linspace(min(rmsd), max(rmsd), 200)

# Plot the PDF
plt.plot(xs, density(xs))

# Add labels to the X and Y axes
plt.xlabel('RMSD (Å)')
plt.ylabel('Probability density')

# Show the plot
plt.show()
