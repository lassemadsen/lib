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
    xfm_file = f'{outfile_short}_0_GenericAffine.xfm'
    
    if resampled_file is not None:
        bash_helper.run_shell(f'itk_resample --like {target_image} --transform {xfm_file} {moving_image} {resampled_file} --clobber')

    if qc_file is not None:
        if resampled_file is None:
            with tempfile.TemporaryDirectory() as tmp_dir:
                bash_helper.run_shell(f'itk_resample --like {target_image} --transform {xfm_file} {moving_image} {tmp_dir}/tmp_img.nii --clobber')
                mask_qc(target_image, f'{tmp_dir}/tmp_img.nii', qc_file, clobber=clobber, mask_qrange=[0.01, 0.99])
        else:
            mask_qc(target_image, resampled_file, qc_file, clobber=clobber, mask_qrange=[0.01, 0.99])

    return xfm_file


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

    # target = '/projects/MINDLAB2013_18-MR-Amnestic-MCI/scratch/datakurtosis_jun25/0004/20131113_091044/MR/KURTOSIS_DKITOOLS139_MD/NATSPACE/0001.nii'
    # moving = '/projects/MINDLAB2013_18-MR-Amnestic-MCI/scratch/datakurtosis_jun25/0004/20131113_091044/MR/T1/NATSPACE/0001.nii'

    # coreg(moving, target, '/public/lama/lib', 'test', qc_file='/public/lama/lib/test.jpg', clobber=True)