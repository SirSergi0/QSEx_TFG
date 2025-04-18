import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Define the function
def f(x, y):
    return 4*x**2+4*y**2+4*x*y-12*x-12*y+9
# Create a grid of x and y values
x = np.linspace(-10, 10, 400)
y = np.linspace(-10, 10, 400)
X, Y = np.meshgrid(x, y)

# Evaluate the function on the grid
Z = f(X, Y)

# Plot the implicit curve f(x, y) = 0
plt.figure(figsize = (10, 7))
plt.contourf(X, Y, Z, levels=[-1000, 0], colors="aqua")
plt.contour(X, Y, Z, levels=[0], colors="turquoise")
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.xlabel('$\lambda_1$')
plt.ylabel('$\lambda_2$')
plt.xlim(0, 3)
plt.ylim(0, 3)
filled_region = mpatches.Patch(color='aqua', label=r'$\Delta$')
plt.legend(handles=[filled_region], loc='upper right')
plt.show()
plt.savefig("PerfectExclusionZoneZ3Z.pdf")
