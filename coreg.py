#!/usr/bin/env python
import bash_helper
import os
import argparse
import tempfile
from mask_qc import mask_qc

# TODO See if outfile can be set with _0_GenericAffine ending
# TODO Add option to perform non-linear and rigid coreg
# TODO Check minc? 

def coreg(moving_image, target_image, outdir, outname, resampled_file=None, qc_file=None, clobber=False):

    outfile_short = f'{outdir}/{outname}'
    outfile = f'{outdir}/{outname}_0_GenericAffine.xfm'

    if os.path.isfile(outfile) and not clobber:
        print('Outfile exits.. Used -clobber to overwrite')
        return

    # if moving_image.split('.')[-1] == 'nii':
    #     bash_helper.run_shell(f'nii2mnc {moving_image} -clobber -float')
    #     moving_image = moving_image.split('nii')[0] + 'mnc'
    # elif moving_image.split('.')[-1] == 'mnc':
    #     pass
    # else:
    #     print('Error. "moving_image" has to be .nii or .mnc')
    #     return

    # if target_image.split('.')[-1] == 'nii':
    #     bash_helper.run_shell(f'nii2mnc {target_image} -clobber -float')
    #     target_image = target_image.split('nii')[0] + 'mnc'
    # elif target_image.split('.')[-1] == 'mnc':
    #     pass
    # else:
    #     print('Error. "target_image" has to be .nii or .mnc')
    #     return


    print(f'Performing coregistration from {moving_image} to {target_image}...')
    
    bash_helper.run_shell('antsRegistration --dimensionality 3 --float 0 '
                        f'--output [{outfile_short}_] '
                        '--interpolation Linear '
                        '--winsorize-image-intensities [0.005,0.995] '
                        '--use-histogram-matching 0 '
                        '--transform Rigid[0.1] '
                        f'--metric MI[{moving_image}, {target_image},1,32,Regular,0.25] '
                        '--convergence [500x250x100,1e-6,10] '
                        '--shrink-factors 4x2x1 '
                        '--smoothing-sigmas 2x1x0vox '
                        '--minc')
    
    if resampled_file is not None:
        bash_helper.run_shell(f'itk_resample --like {target_image} --transform {outfile_short}_0_GenericAffine.xfm {moving_image} {resampled_file} --clobber')

    if qc_file is not None:
        if resampled_file is None:
            with tempfile.TemporaryDirectory() as tmp_dir:
                bash_helper.run_shell(f'itk_resample --like {target_image} --transform {outfile_short}_0_GenericAffine.xfm {moving_image} {tmp_dir}/tmp_img.nii --clobber')
            resampled_file = f'{tmp_dir}/tmp_img.nii'

        mask_qc(target_image, resampled_file, qc_file, clobber=clobber, mask_qrange=[0.01, 0.99])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('moving_image', help='File location of moving image.')
    parser.add_argument('target_image', help='File location of target image.')
    parser.add_argument('outdir', help='Directory of output transformation file.')
    parser.add_argument('outname', help='Name of transfomration file (prefix, i.e. without file type).')
    parser.add_argument('-qc_file', default=None, help='Optinal. Name of QC file.')
    parser.add_argument('-clobber', action='store_true', help='If -clobber, exsiting output files will be overwritten')

    args = parser.parse_args()

    coreg(args.moving_image, args.target_image, args.outdir, args.outname, args.qc_file, args.clobber)
