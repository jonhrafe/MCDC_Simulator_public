import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

def plot_sphere(ax, center, radius, resolution=30):
    """
    Add a sphere to a 3D plot.
    :param ax: Matplotlib 3D axis.
    :param center: Center of the sphere (x, y, z).
    :param radius: Radius of the sphere.
    :param resolution: Resolution for the sphere mesh.
    """
    u = np.linspace(0, 2 * np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)
    x = radius * np.outer(np.cos(u), np.sin(v)) + center[0]
    y = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]
    ax.plot_surface(x, y, z, color='b', alpha=0.3, edgecolor='gray')

def visualize_spheres(file_path):
    """
    Visualize spheres from a file.
    :param file_path: Path to the file containing sphere data.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Skip the first line (scale) and parse the spheres
    spheres = []
    for line in lines[1:]:
        x, y, z, r = map(float, line.split())
        spheres.append((x, y, z, r))

    # Set up the 3D plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio

    # Plot each sphere
    for sphere in spheres:
        plot_sphere(ax, center=sphere[:3], radius=sphere[3])

    # Set axis limits based on the spheres
    all_coords = np.array([(x, y, z) for x, y, z, _ in spheres])
    min_coord = all_coords.min(axis=0)
    max_coord = all_coords.max(axis=0)
    padding = (max_coord - min_coord).max() * 0.1
    ax.set_xlim(min_coord[0] - padding, max_coord[0] + padding)
    ax.set_ylim(min_coord[1] - padding, max_coord[1] + padding)
    ax.set_zlim(min_coord[2] - padding, max_coord[2] + padding)

    # Add labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Show the plot
    plt.show()

# Example usage
file_path = "/home/diffusion/Dropbox/Documents/EPFL/Main_Studies/Simulators/MCDC_Simulator_public/dev/dev_outputs/benchmark_test3_spheres_2025-01-20_21-27-03_gamma_distributed_sphere_list.txt"
visualize_spheres(file_path)
