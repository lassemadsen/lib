#!/usr/bin/env python
"""Convert minc to nifti format."""
import os
from nibabel import load, save, Nifti1Image
import argparse

def convert_minc_to_nifti(input_file, output_file):
    minc = load(input_file)

    if output_file is None:
        basename = minc.get_filename().split(os.extsep, 1)[0]
        output_file = basename + '.nii'

    out = Nifti1Image(minc.get_fdata(), minc.affine, minc.header)

    save(out, output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')

    args = parser.parse_args()

    # class args:
    #     input = '/Volumes/projects/MINDLAB2021_MR-APOE-Microcirculation/scratch/image_data/0040/bl/0040_20230322_PET_OEF.mnc'
    #     output = '/Volumes/projects/MINDLAB2021_MR-APOE-Microcirculation/scratch/image_data/0040/bl/0040_20230322_PET_OEF.nii'
    
    convert_minc_to_nifti(args.input, args.output)
