#!/usr/bin/env python
"""Convert nifti to minc format."""
import os
import numpy as np
import nibabel as nib
import argparse
from pyminc.volumes.factory import volumeFromDescription

def convert_nifti_to_minc(nifti_file, minc_file):
    # Load NIfTI file
    nii = nib.load(nifti_file)
    data = nii.get_fdata()

    # Extract affine transformation
    affine = nii.affine
    sizes = data.shape  # MINC sizes (same as data shape)
    
    # Compute starts (origin) and steps (voxel sizes)
    starts = affine[:3, 3]  # The translation components (origin)
    steps = np.sqrt((affine[:3, :3] ** 2).sum(axis=0))  # The voxel sizes

    # Create MINC file with correct metadata
    vol = volumeFromDescription(minc_file, ("xspace", "zspace", "yspace"), sizes, starts, steps)

    # Assign data to MINC volume and write the file
    vol.data = data
    vol.writeFile()

    print(f"Converted: {nifti_file} → {minc_file}")



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')

    args = parser.parse_args()

    # class args:
    #     input = '/Volumes/projects/MINDLAB2021_MR-APOE-Microcirculation/scratch/masksdce_all_jun24/0006/20211007_112959/MR/T1UNI/DCESPACE/hippo_clean.nii'
    #     output = '/Volumes/projects/MINDLAB2021_MR-APOE-Microcirculation/scratch/masksdce_all_jun24/0006/20211007_112959/MR/T1UNI/DCESPACE/hippo_clean.mnc'

    convert_nifti_to_minc(args.input, args.output)
