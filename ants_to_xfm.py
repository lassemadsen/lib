import numpy as np
import ants
import os

def ants_to_xfm(transform, output_file, clobber=False):
    """
    Conversion from ants transform to .xfm

    Parameters
    ----------
    transform : ants.ANTsTransform or str
        Can be ANTsTransform type or str (location of transform .mat file)
    output_file: str
        Location of output .xfm file.
    """
    if os.path.isfile(output_file) and not clobber:
        print(f'{output_file} exists. Use -clobber to overwrite.')


    if isinstance(transform, str):
        transform = ants.read_transform(transform)
    elif isinstance(transform, ants.ANTsTransform):
        pass
    else:
        print('Error. Transform should be str or ANTsTransform')
        return

    if transform.transform_type != "AffineTransform" or transform.dimension != 3:
        raise ValueError("Only 3D AffineTransform is supported.")
    
    # Extract raw 3x3 matrix and translation from parameters
    params = np.array(transform.parameters)
    matrix = params[:9].reshape(3, 3)
    translation = params[9:]
    
    # Try to get the actual center used by ANTs
    fixed_params = np.array(transform.fixed_parameters)
    if len(fixed_params) >= 3:
        center_phys = fixed_params[:3]
        
    # Build the full 4x4 transformation matrix
    # Transform: y = R*(x - c) + t + c = R*x + (t + c - R*c)
    full_affine = np.eye(4)
    full_affine[:3, :3] = matrix
    full_affine[:3, 3] = translation + center_phys - matrix @ center_phys

    # FOR .XFM FORMAT: Take the inverse!
    inverse_affine = np.linalg.inv(full_affine)
    
    # Write XFM file
    with open(output_file, 'w') as f:
        f.write("MNI Transform File\n")
        f.write("%Generated from ANTs transform\n")
        f.write("\n")
        f.write("Transform_Type = Linear;\n")
        f.write("Linear_Transform =\n")
        for row in inverse_affine[:3, :]:
            f.write(" " + " ".join(f"{val:.14f}" for val in row) + "\n")
        f.write(";\n")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('ants_transform', help='ANTS .mat transform file')
    parser.add_argument('output_xfm', help='File to save the converted .xfm file.')
    parser.add_argument('-clobber', action='store_true')
    args = parser.parse_args()

    ants_to_xfm(args.ants_transform, args.output_xfm, argparse.clobber)