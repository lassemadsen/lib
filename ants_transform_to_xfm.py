import numpy as np
import ants

def ants_transform_to_xfm(transform: ants.ANTsTransform, output_file):
    """
    Conversion from ants transform to .xfm
    """
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