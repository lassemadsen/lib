#!/usr/bin/env python
import bash_helper
import os
import argparse
import tempfile
import ants

# TODO See if outfile can be set with _0_GenericAffine ending
# TODO Add option to perform non-linear and rigid coreg

def coreg(from_file, to_file, outdir, outname, qc_file=None, clobber=False):

    if from_file.split('.')[-1] == 'nii':
        bash_helper.run_shell(f'nii2mnc {from_file} -clobber -float')
        from_file = from_file.split('nii')[0] + 'mnc'
    elif from_file.split('.')[-1] == 'mnc':
        pass
    else:
        print('Error. "from_file" has to be .nii or .mnc')
        return

    if to_file.split('.')[-1] == 'nii':
        bash_helper.run_shell(f'nii2mnc {to_file} -clobber -float')
        to_file = to_file.split('nii')[0] + 'mnc'
    elif to_file.split('.')[-1] == 'mnc':
        pass
    else:
        print('Error. "to_file" has to be .nii or .mnc')
        return

    outfile_short = f'{outdir}/{outname}'
    outfile = f'{outdir}/{outname}_0_GenericAffine.xfm'

    print(f'Performing coregistration from {from_file} to {to_file}...')
    if not os.path.isfile(outfile) or clobber:
        bash_helper.run_shell('antsRegistration --dimensionality 3 --float 0 '
                            f'--output [{outfile_short}_] '
                            '--interpolation Linear '
                            '--winsorize-image-intensities [0.005,0.995] '
                            '--use-histogram-matching 0 '
                            '--transform Rigid[0.1] '
                            f'--metric MI[{to_file},{from_file},1,32,Regular,0.25] '
                            '--convergence [500x250x100,1e-6,10] '
                            '--shrink-factors 4x2x1 '
                            '--smoothing-sigmas 2x1x0vox '
                            '--minc')
        
    if qc_file is not None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            bash_helper.run_shell(f'itk_resample --like {to_file} --transform {outfile_short}_0_GenericAffine.xfm {from_file} {tmp_dir}/tmp_img.nii --clobber')
            resampled = ants.image_read(f'{tmp_dir}/tmp_img.nii')
            to_file = ants.image_read(to_file)

            resampled.plot(overlay = to_file, overlay_alpha=0.7, filename=qc_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('from_file', help='Location of file to coregister from.')
    parser.add_argument('to_file', help='Location of file to coregister to.')
    parser.add_argument('outdir', help='Directory of output transformation file.')
    parser.add_argument('outname', help='Name of transfomration file (prefix, i.e. without file type).')
    parser.add_argument('-qc_file', default=None, help='Optinal. Name of QC file.')
    parser.add_argument('-clobber', action='store_true', help='If -clobber, exsiting output files will be overwritten')

    args = parser.parse_args()

    coreg(args.from_file, args.to_file, args.outdir, args.outname, args.qc_file, args.clobber)
