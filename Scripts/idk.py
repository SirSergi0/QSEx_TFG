import numpy as np
import matplotlib.pyplot as plt

# Example data: list of [x, y] pairs
data = [
    [1.0, 2.0],
    [1.0, 2.0],
    [2.0, 3.0],
    [1.5, 2.5],
    [2.0, 3.0],
    [1.0, 2.0],
    [2.0, 3.0],
    [1.5, 2.5],
]

# Convert to NumPy array
data = np.array(data)
x = data[:, 0]
y = data[:, 1]

# Create 2D histogram
x_bins = 10
y_bins = 10
heatmap, xedges, yedges = np.histogram2d(x, y, bins=[x_bins, y_bins])

# Plot
plt.figure(figsize=(8, 6))
plt.imshow(
    heatmap.T,
    origin='lower',
    aspect='auto',
    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
    cmap='hot'
)
plt.colorbar(label='Frequency')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('2D Heatmap of Frequency')
plt.show()
