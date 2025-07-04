#!/usr/bin/env python
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import nibabel as nib
import ants
from nibabel.orientations import aff2axcodes, io_orientation, axcodes2ornt, ornt_transform, apply_orientation
import numpy as np
# import matplotlib.cm as cm
import argparse
import os
from PIL import Image
import tempfile

def mask_qc(img_file : str, mask_file : str, outfile : str, img_cmap : str = 'gray', mask_cmap : str = 'Spectral', 
                 title : str = '', mask_qrange = None, mask_alpha=0.5, clobber : bool = False):
    
    if os.path.isfile(outfile) and not clobber:
        print(f'{outfile} exists. Use clobber to overwrite.')
        return

    img = ants.image_read(img_file)
    mask = ants.image_read(mask_file)

    if mask_qrange is None:
        vminol = vmaxol = None
    else:
        vminol, vmaxol = mask.quantile(mask_qrange[0]), mask.quantile(mask_qrange[1])
    mask_cmap = transparent_cmap(mask_cmap)

    image_paths = []
    with tempfile.TemporaryDirectory() as tmp_dir:
        for axis in range(len(img.shape)):
            if axis == 0:
                t = title
            else:
                if title == '': # If title is empty, the other titles can also be empty. Else ' ', to allign figures when combining. 
                    t = ''
                else:
                    t = ' '
            filename_axis = f'{tmp_dir}/axis_{axis}.jpg'
            img.plot(
                overlay=mask,
                filename=filename_axis,
                cmap=img_cmap,
                overlay_cmap=mask_cmap,
                overlay_alpha=mask_alpha,
                axis=axis,
                vminol = vminol,
                vmaxol = vmaxol,
                title = t
            )
            image_paths.append(filename_axis)

        # Open and concatenate images horizontally
        images = [Image.open(path) for path in image_paths]
        widths, heights = zip(*(im.size for im in images))

        total_width = sum(widths)
        max_height = max(heights)

        concatenated = Image.new('RGB', (total_width, max_height))

        x_offset = 0
        for im in images:
            concatenated.paste(im, (x_offset, 0))
            x_offset += im.width

        # Save final combined image
        concatenated.save(outfile)

def mask_qc_old(img_file : str, mask_file : str, outfile : str, img_cmap : str = 'gray', mask_cmap : str = 'Spectral', 
         img_range : tuple = None, mask_range : tuple = None, file_id : str = '', mask_id : str = '', 
         info : str = '', mask_alpha=0.4, mask_is_image=False, clobber : bool = False):
    
    if not clobber and os.path.exists(outfile):
        print(f'{outfile} already exists. Use -clobber to overwrite.')
        return

    img = nib.load(img_file)
    img_data = img.get_fdata()
    img_affine = img.affine

    mask = nib.load(mask_file)
    mask_data = mask.get_fdata()
    mask_affine = mask.affine

    # Check identical affine by comparing dimension and determinant:
    det1 = np.linalg.det(img_affine[:3, :3])
    det2 = np.linalg.det(mask_affine[:3, :3])

    if not np.isclose(det1,det2):
        print(f'Error. Image ({img_file}) and mask ({mask_file}) do not have same dimension..')
        return

    axcodes = aff2axcodes(img_affine)

    dim_names = []
    target_orientation = []
    for a in axcodes:
        if a in ['L', 'R']:
            dim_names.extend(['x'])
            target_orientation.extend(['R'])
        elif a in ['A', 'P']:
            dim_names.extend(['y'])
            target_orientation.extend(['A'])
        elif a in ['S', 'I']:
            dim_names.extend(['z'])
            target_orientation.extend(['S'])
    
    target_ornt = axcodes2ornt(target_orientation)
    current_ornt = io_orientation(img_affine)
    transform = ornt_transform(current_ornt, target_ornt)
    img_data = apply_orientation(img_data, transform)
    mask_data = apply_orientation(mask_data, transform)
        
    if img_range is None:
        img_range = np.nanpercentile(img_data,[1,99])
    else:
        assert len(img_range) == 2, f"Image range should be [min, max]. Got: {img_range}"

    if mask_range is None:
        if mask_is_image:
            mask_range = np.nanpercentile(mask_data,[1,99])
        else:
            mask_range = [np.nanmin(mask_data), np.nanmax(mask_data)] # Dont display 0 values 
    else:
        assert len(mask_range) == 2, f"Mask range should be [min, max]. Got: {mask_range}"
    mask_cmap = transparent_cmap(mask_cmap)

    views = {dim_names[0]: [round(i) for i in np.linspace(img_data.shape[0]/4,img_data.shape[0]-(img_data.shape[0]/4),9)],
             dim_names[1]: [round(i) for i in np.linspace(img_data.shape[1]/4,img_data.shape[1]-(img_data.shape[1]/4),9)],
             dim_names[2]: [round(i) for i in np.linspace(img_data.shape[2]/4,img_data.shape[2]-(img_data.shape[2]/4),9)]}

    if img_file.split('.')[-1] == 'mnc':
        img_data_show = {dim_names[0]: [img_data[v,:,:] for v in views[dim_names[0]]],
                         dim_names[1]: [img_data[:,v,:] for v in views[dim_names[1]]],
                         dim_names[2]: [img_data[:,:,v] for v in views[dim_names[2]]]}
    else:
        img_data_show = {dim_names[0]: [np.rot90(img_data[v,:,:]) for v in views[dim_names[0]]],
                         dim_names[1]: [np.rot90(img_data[:,v,:]) for v in views[dim_names[1]]],
                         dim_names[2]: [np.rot90(img_data[:,:,v]) for v in views[dim_names[2]]]}
    
    if mask_file.split('.')[-1] == 'mnc':
        mask_data_show = {dim_names[0]: [mask_data[v,:,:] for v in views[dim_names[0]]],
                          dim_names[1]: [mask_data[:,v,:] for v in views[dim_names[1]]],
                          dim_names[2]: [mask_data[:,:,v] for v in views[dim_names[2]]]}
    else:
        mask_data_show = {dim_names[0]: [np.rot90(mask_data[v,:,:]) for v in views[dim_names[0]]],
                          dim_names[1]: [np.rot90(mask_data[:,v,:]) for v in views[dim_names[1]]],
                          dim_names[2]: [np.rot90(mask_data[:,:,v]) for v in views[dim_names[2]]]}

    voxel_sizes = {dim_names[0]: img.header.get_zooms()[:3][0],
                   dim_names[1]: img.header.get_zooms()[:3][1],
                   dim_names[2]: img.header.get_zooms()[:3][2]}
       
    fig = plt.figure(figsize=(17,10))

    ax = []

    for row in range(6):
        for col in range(9):
            if (row == 0) & (col == 0):
                a = plt.subplot2grid(shape=(6,9), loc=(row, col), colspan=3, rowspan=3)
            elif (row == 0) & (col == 3):
                a = plt.subplot2grid(shape=(6,9), loc=(row, col), colspan=3, rowspan=3)
            elif (row == 0) & (col == 6):
                a = plt.subplot2grid(shape=(6,9), loc=(row, col), colspan=3, rowspan=3)
            elif row > 2:
                a = plt.subplot2grid(shape=(6,9), loc=(row, col))
            else:
                continue

            ax.extend([a])

    fig.set_facecolor('k')

    x_views = list(range(len(views['x'])))
    y_views = list(range(len(views['y'])))
    z_views = list(range(len(views['z'])))

    extent = {'x': [0, img_data_show['x'][0].shape[0]*voxel_sizes['z'], 0, img_data_show['x'][0].shape[1]*voxel_sizes['y']],
              'y': [0, img_data_show['y'][0].shape[0]*voxel_sizes['z'], 0, img_data_show['y'][0].shape[1]*voxel_sizes['x']],
              'z': [0, img_data_show['z'][0].shape[0]*voxel_sizes['y'], 0, img_data_show['z'][0].shape[1]*voxel_sizes['x']]}
    
    aspect = {'x': extent['x'][1]/extent['x'][3],
              'y': extent['y'][1]/extent['y'][3],
              'z': extent['z'][1]/extent['z'][3]}
    
    origin = {'image': 'lower' if img_file.split('.')[-1] == 'mnc' else 'upper',
              'mask': 'lower' if mask_file.split('.')[-1] == 'mnc' else 'upper',} 

    for a in ax:
        large_view = 5 # Location of the large (top row) view from x,y,z_views lists. 

        # Top row - Large view
        if a.get_subplotspec().get_geometry()[2] == 0:
            a.imshow(img_data_show['x'][large_view], cmap=img_cmap, vmin=img_range[0], vmax=img_range[1], interpolation='nearest', extent=extent['x'], aspect=aspect['x'], origin=origin['image'])
            a.imshow(mask_data_show['x'][large_view], alpha=mask_alpha, cmap=mask_cmap, vmin=mask_range[0], vmax=mask_range[1], interpolation='nearest', extent=extent['x'], aspect=aspect['x'], origin=origin['mask'])
        elif a.get_subplotspec().get_geometry()[2] == 3:
            a.imshow(img_data_show['y'][large_view], cmap=img_cmap, vmin=img_range[0], vmax=img_range[1], interpolation='nearest', extent=extent['y'], aspect=aspect['y'], origin=origin['image'])
            a.imshow(mask_data_show['y'][large_view], alpha=mask_alpha, cmap=mask_cmap, vmin=mask_range[0], vmax=mask_range[1], interpolation='nearest', extent=extent['y'], aspect=aspect['y'], origin=origin['mask'])
        elif a.get_subplotspec().get_geometry()[2] == 6:
            a.imshow(img_data_show['z'][large_view], cmap=img_cmap, vmin=img_range[0], vmax=img_range[1], interpolation='nearest', extent=extent['z'], aspect=aspect['x'], origin=origin['image'])
            a.imshow(mask_data_show['z'][large_view], alpha=mask_alpha, cmap=mask_cmap, vmin=mask_range[0], vmax=mask_range[1], interpolation='nearest', extent=extent['z'], aspect=aspect['z'], origin=origin['mask'])

        # Small views
        elif a.get_subplotspec().get_geometry()[2] in [27, 28, 29, 36, 37, 38, 45, 46, 47]:
            x_view = x_views.pop(0)
            a.imshow(img_data_show['x'][x_view], cmap=img_cmap, vmin=img_range[0], vmax=img_range[1], interpolation='nearest', extent=extent['x'], aspect=aspect['x'], origin=origin['image'])
            a.imshow(mask_data_show['x'][x_view], alpha=mask_alpha, cmap=mask_cmap, vmin=mask_range[0], vmax=mask_range[1], interpolation='nearest', extent=extent['x'], aspect=aspect['x'], origin=origin['mask'])
        elif a.get_subplotspec().get_geometry()[2] in [30, 31, 32, 39, 40, 41, 48, 49, 50]:
            y_view = y_views.pop(0)
            a.imshow(img_data_show['y'][y_view], cmap=img_cmap, vmin=img_range[0], vmax=img_range[1], interpolation='nearest', extent=extent['y'], aspect=aspect['y'], origin=origin['image'])
            a.imshow(mask_data_show['y'][y_view], alpha=mask_alpha, cmap=mask_cmap, vmin=mask_range[0], vmax=mask_range[1], interpolation='nearest', extent=extent['y'], aspect=aspect['y'], origin=origin['mask'])
        elif a.get_subplotspec().get_geometry()[2] in [33, 34, 35, 42, 43, 44, 51, 52, 53]:
            z_view = z_views.pop(0)
            a.imshow(img_data_show['z'][z_view], cmap=img_cmap, vmin=img_range[0], vmax=img_range[1], interpolation='nearest', extent=extent['z'], aspect=aspect['z'], origin=origin['image'])
            a.imshow(mask_data_show['z'][z_view], alpha=mask_alpha, cmap=mask_cmap, vmin=mask_range[0], vmax=mask_range[1], interpolation='nearest', extent=extent['z'], aspect=aspect['z'], origin=origin['mask'])

        a.axis('off')

    plt.tight_layout()
    # plt.subplots_adjust(wspace=-0.1, hspace=-0.1)
    fig.text(.02, .98, f'{file_id} {img_file.split("/")[-1]}', horizontalalignment='left', verticalalignment='center', color='white')
    fig.text(.02, .96, f'{mask_id} {mask_file.split("/")[-1]}', horizontalalignment='left', verticalalignment='center', color='white')
    fig.text(.02, .94, f'{info}', horizontalalignment='left', verticalalignment='center', color='white')

    plt.savefig(outfile, dpi=300)
    plt.close()
    

    # img = ants.image_read(img_file)
    # mask = ants.image_read(mask_file)

    # img_data = img.numpy()
    # mask_data = mask.numpy()

    # spacing = img.spacing  # (x, y, z)

    # # Orientation check (basic)
    # if not np.allclose(img.spacing, mask.spacing, rtol=1e-2):
    #     print("Error. Image and mask have different spacing.")
    #     return


def transparent_cmap(cmap_str):

    cmap = matplotlib.colormaps.get_cmap(cmap_str)  #cm.get_cmap(cmap_str)
    # Modify the colormap's alpha values to make lowest values transparent
    cmap_colors = cmap(np.arange(cmap.N))
    cmap_colors[:, 3] = np.where(np.arange(cmap.N) < 1, 0, cmap_colors[:, 3])
    cmap = cmap.from_list(cmap.name + '_transparent', cmap_colors, cmap.N)

    return cmap


if __name__ == "__main__": 
    # mask_qc('/Volumes/projects/MINDLAB2021_MR_ENIGMA/scratch/dataDWI_seg/0004/20210427_124908/dwi_trace.nii', '/Volumes/projects/MINDLAB2021_MR_ENIGMA/scratch/masksDWI_seg/0004/24h.nii', '/Users/au483096/Desktop/test_DWI.jpg', 
    #         clobber=True, mask_alpha=1, mask_cmap='Reds')
    # mask_qc('/Volumes/projects/MINDLAB2013_18-MR-Amnestic-MCI/scratch/datakurtosis_jun25/0004/20131113_091044/MR/T1/DKISPACE/0001.nii', 
    #         '/Volumes/projects/MINDLAB2013_18-MR-Amnestic-MCI/scratch/datakurtosis_jun25/0004/20131113_091044/MR/KURTOSIS_DKITOOLS139_MD/NATSPACE/0001.nii', 
    #         '/Users/au483096/Desktop/test_DKI.jpg', clobber=True, mask_qrange = [0.001, 0.999])
    # mask_qc('/Volumes/projects/MINDLAB2021_MR_ENIGMA/scratch/dataSEP_MR_feb25/0004/20210427_124908/MR/SEPWIMEAN/NATSPACE/0001.nii', '/Volumes/projects/MINDLAB2021_MR_ENIGMA/scratch/dataSEP_MR_feb25/0004/20210427_124908/MR/T1UNI/SEPWIMEAN/T1_resampled.nii', '/Users/au483096/Desktop/test_SEPWI.jpg', clobber=True)

    parser = argparse.ArgumentParser(description='Script to plot QC of minc image and mask')
    parser.add_argument('img_file', help='Location of image file to show.')
    parser.add_argument('mask_file', help='Location of mask file to overlay on image.')
    parser.add_argument('outfile', help='Location of output. Should be .jpg or .png.')
    parser.add_argument('-img_cmap', default='gray', help='Matlplotlib cmap to plot image. Default is gray.')
    parser.add_argument('-mask_cmap', default='jet', help='Matlplotlib cmap to plot mask. Default is jet.')
    parser.add_argument('-title', default='', help='Add additional information to be shown on the figure.')
    parser.add_argument('-mask_alpha', default=0.4, help='Alpha (transparancy) of mask. Between 0 and 1.')
    parser.add_argument('-clobber', action='store_true', help='If -clobber, existing file will be overwritten.')

    args = parser.parse_args()    

    mask_qc(args.img_file, args.mask_file, args.outfile, img_cmap=args.img_cmap, mask_cmap=args.mask_cmap, 
            title=args.title, mask_alpha=args.mask_alpha, clobber=args.clobber)

