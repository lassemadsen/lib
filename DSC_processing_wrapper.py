#!/usr/bin/env python3
from pathlib import Path
import nibabel as nib
import ants
from DSC_processing import DSC_process
from stormdb.access import Query
from stormdb_functions import convert_pwi
import argparse

if str(Path('__file__').resolve()).startswith('/Volumes'):
    path_prefix = '/Volumes'
else:
    path_prefix = ''

def main(args):
    project = args.project
    branch = args.branch
    structural_branch = args.structural_branch
    subject = args.subject
    timepoint = args.timepoint
    pwi_type = args.pwi_type

    q = Query(project)
    filtered_series = q.filter_series(types=pwi_type, subjects=subject, return_files=True)
    filtered_series = [f for f in filtered_series if f['study'] == timepoint][0]

    project_dir = f'{path_prefix}/projects/{project}/scratch'

    sub_id = filtered_series['subject'].zfill(4)
    tp = filtered_series['study']
    pwi_type = filtered_series['type']
    img_file = f'{project_dir}/data{branch}/{sub_id}/{tp}/MR/{pwi_type}/NATSPACE/0001.nii'

    print(f'Converting {pwi_type} image...', end='', flush=True)
    info = convert_pwi(filtered_series, img_file)

    print(f'Running DSC {pwi_type} processing for {sub_id} - {timepoint}...', flush=True)
    proc = DSC_process(sub_id, tp, pwi_type, img_file, info, project_dir, branch)
    
    print(f'Performing slice-time correction...', end='', flush=True)
    proc.slice_time_correction()
    print(' \u2713')
    print(f'Masking image...', end='', flush=True)
    proc.mask_image()
    print(' \u2713')
    print(f'Detecting baseline...', end='', flush=True)
    proc.baseline_detection()
    print(' \u2713')
    print(f'Truncating signal...', end='', flush=True)
    proc.truncate_signal()
    print(' \u2713')
    print(f'Performing motion correction...', end='', flush=True)
    proc.motion_correction()
    print(' \u2713')
    print(f'Calculating concentration...', end='', flush=True)
    proc.calc_concentration()
    print(' \u2713')

    t1_type = 'T1UNI'
    t1_file = f'{project_dir}/data{structural_branch}/{sub_id}/{tp}/MR/{t1_type}/NATSPACE/0001.nii'
    t1_mask_dir = f'{project_dir}/masks{structural_branch}/{sub_id}/{tp}/MR/{t1_type}/NATSPACE'

    dsc_mean = ants.image_read(f'{proc.data_dir}/{pwi_type}MEAN/NATSPACE/0001.nii')
    t1 = ants.image_read(t1_file)

    transform = ants.registration(dsc_mean, t1, type_of_transform = 'SyN')
    # ants.write_transform(transform['fwdtransforms'], f'{outdir_info}/t1_{pwi_type}MEAN_transform')

    t1_resampled = ants.apply_transforms(fixed=dsc_mean, moving=t1, transformlist=transform['fwdtransforms'])
    Path(f'{proc.data_dir}/{t1_type}/{pwi_type}SPACE').mkdir(parents=True, exist_ok=True)
    ants.image_write(t1_resampled, f'{proc.data_dir}/{t1_type}/{pwi_type}SPACE/0001.nii')

    mask_dir_pwi_space = f'{proc.mask_dir}/{t1_type}/{pwi_type}SPACE'
    Path(mask_dir_pwi_space).mkdir(parents=True, exist_ok=True)

    for mask_file in ['struc_gm.nii', 'struc_wm.nii', 'AIFsearchMask.nii']:
        mask = ants.image_read(f'{t1_mask_dir}/{mask_file}')

        interpolator = 'nearestNeighbor'

        mask_resampled = ants.apply_transforms(fixed=dsc_mean, moving=mask, transformlist=transform['fwdtransforms'], interpolator=interpolator)

        ants.image_write(mask_resampled, f'{mask_dir_pwi_space}/{mask_file}')

    # Create smooth mask
    gm_mask_file = f'{mask_dir_pwi_space}/struc_gm.nii'
    wm_mask_file = f'{mask_dir_pwi_space}/struc_wm.nii'
    gm_mask = nib.load(gm_mask_file)
    wm_mask = nib.load(wm_mask_file)
    smooth_mask = (gm_mask.get_fdata()) + (wm_mask.get_fdata())

    smooth_mask_file = f'{mask_dir_pwi_space}/smooth_mask.nii'
    smooth_mask_to_save = nib.Nifti1Image(smooth_mask, gm_mask.affine, gm_mask.header)
    nib.save(smooth_mask_to_save, smooth_mask_file)

    aif_search_mask_file = f'{mask_dir_pwi_space}/AIFsearchMask.nii'

    proc.smooth_data(smooth_mask_file)
    proc.aif_selection(aif_search_mask_file, gm_mask_file)
    proc.calc_perfusion()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('project', help='Name of project')
    parser.add_argument('branch', help='Name of branch')
    parser.add_argument('structural_branch', help='Name of structural branch (location of T1 image and T1 masks)')
    parser.add_argument('subject', help='Subject ID (from stormdb)')
    parser.add_argument('pwi_type', help='Type of DSC image (can be PWI or SEPWI)')
    parser.add_argument('timepoint', help='Timepoint of DSC image')

    args = parser.parse_args()

    main(args)