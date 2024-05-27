#!/usr/bin/env python
import matplotlib
# matplotlib.use('agg')
import matplotlib.pyplot as plt
import pyminc.volumes.factory as minc
import numpy as np
import matplotlib.cm as cm
import argparse
import os

def main(args):
    img_file = args.img_file
    mask_file = args.mask_file
    outfile = args.outfile
    img_cmap = args.img_cmap
    mask_cmap = args.mask_cmap
    img_range = args.img_range
    mask_range = args.mask_range
    file_id = args.file_id
    mask_id = args.mask_id
    info = args.info
    clobber = args.clobber

    minc_qc(img_file, mask_file, outfile, img_cmap=img_cmap, mask_cmap=mask_cmap, img_range=img_range, mask_range=mask_range, file_id=file_id, mask_id=mask_id, info=info, clobber=clobber)

def minc_qc(img_file, mask_file, outfile, img_cmap='gray', mask_cmap='Spectral', img_range=None, mask_range=None, file_id='', mask_id='', info='', clobber=False):

    if not clobber and os.path.exists(outfile):
        print(f'{outfile} already exists. Use -clobber to overwrite.')
        return

    img_minc = minc.volumeFromFile(img_file)
    img_data = np.array(img_minc[:])

    mask_minc = minc.volumeFromFile(mask_file)
    mask_data = np.array(mask_minc[:])
    assert img_data.shape == mask_data.shape, 'Error: Dimensions of images do not match!'

    if img_range is None:
        img_range = [img_data.min(), img_data.max()]
    else:
        assert len(img_range) == 2, f"Image range should be [min, max]. Got: {img_range}"

    if mask_range is None:
        mask_range = [mask_data.min(), mask_data.max()]
    else:
        assert len(mask_range) == 2, f"Mask range should be [min, max]. Got: {mask_range}"
        mask_cmap = transparent_cmap(mask_cmap)

    dims = {img_minc.dimnames[0]: img_data.shape[0],
            img_minc.dimnames[1]: img_data.shape[1],
            img_minc.dimnames[2]: img_data.shape[2]}
    
    steps = {img_minc.dimnames[0]: img_minc.separations[0],
             img_minc.dimnames[1]: img_minc.separations[1],
             img_minc.dimnames[2]: img_minc.separations[2]}
    
    views = {'xspace': [round(i) for i in np.linspace(dims['xspace']/4,dims['xspace']-(dims['xspace']/4),9)],
             'yspace': [round(i) for i in np.linspace(dims['yspace']/4,dims['yspace']-(dims['yspace']/4),9)],
             'zspace': [round(i) for i in np.linspace(dims['zspace']/4,dims['zspace']-(dims['zspace']/4),9)]}

    img_data_show = {img_minc.dimnames[0]: [img_data[v,:,:] for v in views[img_minc.dimnames[0]]],
                    img_minc.dimnames[1]: [img_data[:,v,:] for v in views[img_minc.dimnames[1]]],
                    img_minc.dimnames[2]: [img_data[:,:,v] for v in views[img_minc.dimnames[2]]]}
    
    mask_data_show = {mask_minc.dimnames[0]: [mask_data[v,:,:] for v in views[mask_minc.dimnames[0]]],
                      mask_minc.dimnames[1]: [mask_data[:,v,:] for v in views[mask_minc.dimnames[1]]],
                      mask_minc.dimnames[2]: [mask_data[:,:,v] for v in views[mask_minc.dimnames[2]]]}
    
    # Flip data  
    if steps['xspace'] > 0 and steps['yspace'] > 0 and steps['zspace'] > 0:
        img_data_show['zspace'] = [np.flip(v,0) for v in img_data_show['zspace']]
        mask_data_show['zspace'] = [np.flip(v,0) for v in mask_data_show['zspace']]

        img_data_show['yspace'] = [np.flip(v,0) for v in img_data_show['yspace']]
        mask_data_show['yspace'] = [np.flip(v,0) for v in mask_data_show['yspace']]

        img_data_show['xspace'] = [np.flip(v,0) for v in img_data_show['xspace']]
        mask_data_show['xspace'] = [np.flip(v,0) for v in mask_data_show['xspace']]
    elif steps['xspace'] > 0 and steps['yspace'] < 0 and steps['zspace'] > 0:
        img_data_show['zspace'] = [np.transpose(v) for v in img_data_show['zspace']]
        # img_data_show['zspace'] = [np.flip(v) for v in img_data_show['zspace']]
        mask_data_show['zspace'] = [np.transpose(v) for v in mask_data_show['zspace']]
        # mask_data_show['zspace'] = [np.flip(v) for v in mask_data_show['zspace']]

        img_data_show['yspace'] = [np.transpose(v) for v in img_data_show['yspace']]
        img_data_show['yspace'] = [np.flip(v,0) for v in img_data_show['yspace']]
        mask_data_show['yspace'] = [np.transpose(v) for v in mask_data_show['yspace']]
        mask_data_show['yspace'] = [np.flip(v,0) for v in mask_data_show['yspace']]
        
        img_data_show['xspace'] = [np.flip(v) for v in img_data_show['xspace']]
        mask_data_show['xspace'] = [np.flip(v) for v in mask_data_show['xspace']]
    else: 
        print('Warning! Direction of minc files not implemented. Images may be flipped.. ')

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

    alpha = 0.5

    x_views = list(range(len(views['xspace'])))
    y_views = list(range(len(views['yspace'])))
    z_views = list(range(len(views['zspace'])))

    for a in ax:
        large_view = 5 # Location of the large (top row) view from x,y,z_views lists. 

        # Top row - Large view
        if a.get_subplotspec().get_geometry()[2] == 0:
            a.imshow(img_data_show['xspace'][large_view], cmap=img_cmap, vmin=img_range[0], vmax=img_range[1], interpolation='nearest')
            a.imshow(mask_data_show['xspace'][large_view], alpha=alpha, cmap=mask_cmap, vmin=mask_range[0], vmax=mask_range[1], interpolation='nearest')
        elif a.get_subplotspec().get_geometry()[2] == 3:
            a.imshow(img_data_show['yspace'][large_view], cmap=img_cmap, vmin=img_range[0], vmax=img_range[1], interpolation='nearest')
            a.imshow(mask_data_show['yspace'][large_view], alpha=alpha, cmap=mask_cmap, vmin=mask_range[0], vmax=mask_range[1], interpolation='nearest')
        elif a.get_subplotspec().get_geometry()[2] == 6:
            a.imshow(img_data_show['zspace'][large_view], cmap=img_cmap, vmin=img_range[0], vmax=img_range[1], interpolation='nearest')
            a.imshow(mask_data_show['zspace'][large_view], alpha=alpha, cmap=mask_cmap, vmin=mask_range[0], vmax=mask_range[1], interpolation='nearest')

        # Small views
        elif a.get_subplotspec().get_geometry()[2] in [27, 28, 29, 36, 37, 38, 45, 46, 47]:
            x_view = x_views.pop(0)
            a.imshow(img_data_show['xspace'][x_view], cmap=img_cmap, vmin=img_range[0], vmax=img_range[1], interpolation='nearest')
            a.imshow(mask_data_show['xspace'][x_view], alpha=alpha, cmap=mask_cmap, vmin=mask_range[0], vmax=mask_range[1], interpolation='nearest')
        elif a.get_subplotspec().get_geometry()[2] in [30, 31, 32, 39, 40, 41, 48, 49, 50]:
            y_view = y_views.pop(0)
            a.imshow(img_data_show['yspace'][y_view], cmap=img_cmap, vmin=img_range[0], vmax=img_range[1], interpolation='nearest')
            a.imshow(mask_data_show['yspace'][y_view], alpha=alpha, cmap=mask_cmap, vmin=mask_range[0], vmax=mask_range[1], interpolation='nearest')
        elif a.get_subplotspec().get_geometry()[2] in [33, 34, 35, 42, 43, 44, 51, 52, 53]:
            z_view = z_views.pop(0)
            a.imshow(img_data_show['zspace'][z_view], cmap=img_cmap, vmin=img_range[0], vmax=img_range[1], interpolation='nearest')
            a.imshow(mask_data_show['zspace'][z_view], alpha=alpha, cmap=mask_cmap, vmin=mask_range[0], vmax=mask_range[1], interpolation='nearest')

        a.axis('off')

    plt.tight_layout()
    plt.subplots_adjust(wspace=-0.1, hspace=-0.1)
    fig.text(.02, .98, f'{file_id} {img_file.split("/")[-1]}', horizontalalignment='left', verticalalignment='center', color='white')
    fig.text(.02, .96, f'{mask_id} {mask_file.split("/")[-1]}', horizontalalignment='left', verticalalignment='center', color='white')
    fig.text(.02, .94, f'{info}', horizontalalignment='left', verticalalignment='center', color='white')

    plt.savefig(outfile)
    

def transparent_cmap(cmap_str):

    cmap = cm.get_cmap(cmap_str)
    # Modify the colormap's alpha values to make lowest values transparent
    cmap_colors = cmap(np.arange(cmap.N))
    cmap_colors[:, 3] = np.where(np.arange(cmap.N) < 1, 0, cmap_colors[:, 3])
    cmap = cmap.from_list(cmap.name + '_transparent', cmap_colors, cmap.N)

    return cmap


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to plot QC of minc image and mask')
    parser.add_argument('img_file', help='Location of image file to show.')
    parser.add_argument('mask_file', help='Location of mask file to overlay on image.')
    parser.add_argument('outfile', help='Location of output. Should be .jpg or .png.')
    parser.add_argument('-img_cmap', default='gray', help='Matlplotlib cmap to plot image. Default is gray.')
    parser.add_argument('-mask_cmap', default='jet', help='Matlplotlib cmap to plot mask. Default is jet.')
    parser.add_argument('-img_range', nargs=2, default=None, help='Min Max range of image. Autoscale is none provided.', metavar=('min', 'max'))
    parser.add_argument('-mask_range', nargs=2, default=None, help='Min Max range of mask. Autoscale is none provided.', metavar=('min', 'max'))
    parser.add_argument('-file_id', default='', help='Add additional information about the image, e.g. sub_id, to be shown with the filename.')
    parser.add_argument('-mask_id', default='', help='Add additional information about the mask, e.g. sub_id, to be shown with the filename.')
    parser.add_argument('-info', default='', help='Add additional information to be shown on the figure.')
    parser.add_argument('-clobber', action='store_true', help='If -clobber, existing file will be overwritten.')

    args = parser.parse_args()

    # # DEBUGGNING
    # class args:
    #     # img_file = '/Users/au483096/Desktop/stx2_0004_20211007_082130_t1.mnc' # ZYX
    #     # mask_file = '/Users/au483096/Desktop/gm.mnc' # ZYX
    #     # outfile = '/Users/au483096/Desktop/test_flip_zyx.jpg'
    #     img_file = '/Users/au483096/Desktop/0001.mnc' # XZY
    #     mask_file = '/Users/au483096/Desktop/perfusion_roi.mnc' # XZY
    #     outfile = '/Users/au483096/Desktop/test_flip_xzy.jpg'
    #     img_cmap = 'gray'
    #     mask_cmap = 'jet'
    #     img_range = None
    #     mask_range = None
    #     file_id = ''
    #     mask_id = ''
    #     clobber = True

    main(args)
