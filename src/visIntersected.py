import trimesh
import numpy as np
import resource


def set_memory_limit(max_mem_gb):
    """
    Set the maximum address space (virtual memory) limit for the process.

    Parameters:
    - max_mem_gb (float): Maximum memory in gigabytes.
    """
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    max_mem_bytes = max_mem_gb * 1024 ** 3  # Convert GB to Bytes
    resource.setrlimit(resource.RLIMIT_AS, (max_mem_bytes, hard))

def visualize_intersected_trimesh(ply_file, mask_file):
    """
    1) Load the original .ply mesh via trimesh
    2) Load a mask file (0/1 per face)
    3) Filter and display only the faces where mask == 1
    """
    # 1) Load the full mesh
    mesh = trimesh.load(ply_file)
    if not isinstance(mesh, trimesh.Trimesh):
        print("Error: The file doesn't seem to be a single mesh or was not read properly.")
        return

    # 2) Load the mask, one integer per line
    mask = np.loadtxt(mask_file, dtype=int)
    if len(mask) != len(mesh.faces):
        print(f"Error: Mask length ({len(mask)}) != #faces in mesh ({len(mesh.faces)}).")
        return

    # 3) Filter the faces
    face_indices = np.where(mask == 1)[0]
    if len(face_indices) == 0:
        print("No faces marked as '1'. Nothing to display.")
        return

    # Extract the sub-mesh of intersected faces
    submesh = mesh.submesh([face_indices], append=True)
    # Alternatively:
    # submesh = trimesh.Trimesh(vertices=mesh.vertices,
    #                           faces=mesh.faces[mask == 1],
    #                           process=False)

    # 4) Show the submesh in an interactive viewer
    #submesh.show()
    submesh.export('intersected_model_CPU_output_mesh.ply')


if __name__ == "__main__":
    # Example usage:
    ply_file = "/home/jonathan/Downloads/output_mesh.ply"     # path to your original PLY filew //hp_outer.ply"; // output_mesh.ply

    set_memory_limit(20)
    mask_file = "intersected_faces_CPU_output_mesh.txt"     # path to your mask file (0/1 per triangle)
    visualize_intersected_trimesh(ply_file, mask_file)
    #submesh.export('intersected_model.ply')
