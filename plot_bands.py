import matplotlib.pyplot as plt
import numpy as np
import sys

# Check if the correct number of arguments is provided
if len(sys.argv) != 2:
    print("Usage: python plot_bands.py <eigenvalues.dat>")
    sys.exit(1)

# Get the filename from the command-line argument
filename = sys.argv[1]

try:
    # Load the data, skipping the header line
    data = np.loadtxt(filename, comments='#')

    # Extract kx values (first column)
    kx = data[:, 0]

    # Extract eigenvalues (remaining columns)
    eigenvalues = data[:, 1:]

    # Plot each eigenvalue band
    for i in range(eigenvalues.shape[1]):
        plt.plot(kx, eigenvalues[:, i])

    plt.xlabel('kx')
    plt.ylabel('Eigenvalue')
    plt.title('Band Structure')
    plt.grid(True)
    plt.show()

except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {e}")
    sys.exit(1)