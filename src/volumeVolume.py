import numpy as np
import nibabel as nib
import argparse
import os

def create_nifti_from_volume_fraction(input_file, n, output_prefix):
    """
    Convert volume fraction data to two NIfTI files:
    1. density.nii.gz: Total particles
    2. volume_fractions.nii.gz: Relative intra and extra volume fractions
    
    Args:
        input_file (str): Path to the input volume fraction text file
        n (int): Number of subdivisions in each dimension
        output_prefix (str): Prefix for the output files
    """
    # Read the input file
    data = np.loadtxt(input_file)
    
    # Reshape the data into nxnxnx3
    # First reshape to (n^3, 3) to separate the three values per voxel
    data = data.reshape(-1, 3)
    
    # Calculate relative volume fractions
    total = data[:, 0]
    # Avoid division by zero by adding a small epsilon
    epsilon = 1e-10
    data[:, 1] = data[:, 1] / (total + epsilon)  # intra/total
    data[:, 2] = data[:, 2] / (total + epsilon)  # extra/total
    
    # Then reshape to nxnxnx3
    data = data.reshape(n, n, n, 3)
    
    # Create and save density file (first channel)
    density_img = nib.Nifti1Image(data[..., 0], np.eye(4))
    density_file = output_prefix + "_density.nii.gz"
    nib.save(density_img, density_file)
    print(f"Saved density file to: {density_file}")
    
    # Create and save volume fractions file (second and third channels)
    vol_frac_img = nib.Nifti1Image(data[..., 1:], np.eye(4))
    vol_frac_file = output_prefix + "_volume_fractions.nii.gz"
    nib.save(vol_frac_img, vol_frac_file)
    print(f"Saved volume fractions file to: {vol_frac_file}")

def main():
    # parser = argparse.ArgumentParser(description='Convert volume fraction data to NIfTI format')
    # parser.add_argument('input_file', help='Path to the input volume fraction text file')
    # parser.add_argument('n', type=int, help='Number of subdivisions in each dimension')
    # parser.add_argument('output_prefix', help='Prefix for the output files')
    
    # args = parser.parse_args()

    input_file = "/home/diffusion/Dropbox/Documents/EPFL/Main_Studies/Simulators/MCDC_Simulator_public/dev/initialization_test/output/init_test_rep_33_volFractions.txt"
    n = 20
    output_prefix = "/home/diffusion/Dropbox/Documents/EPFL/Main_Studies/Simulators/MCDC_Simulator_public/dev/initialization_test/output/init_test_rep_33"
    
    create_nifti_from_volume_fraction(input_file, n, output_prefix)
    
    # create_nifti_from_volume_fraction(args.input_file, args.n, args.output_prefix)

if __name__ == "__main__":
    main()
