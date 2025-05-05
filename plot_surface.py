#!/usr/bin/env python
import argparse
import numpy as np
import os

def main(args):
    # Check clobber
    if os.path.isfile(args.output) and not args.clobber:
        print(f'{args.output} exists! Use clobber to overwrite.')
        return

    from surface_plot import plot_surface
    views = args.views

    # Load data
    data_left = np.genfromtxt(args.data_left, dtype=float, missing_values="")
    data_right = np.genfromtxt(args.data_right, dtype=float, missing_values="")

    # Skip header row # TODO check if there is a header row
    data_left = data_left[1:]
    data_right = data_right[1:]

    data = {'left': data_left,  
            'right': data_right}
    # check_data_header(data)

    surface = {'left': args.surface_left,
               'right': args.surface_right}

    # Define mask
    if args.mask is not None:
        mask_left = np.loadtxt(args.mask[0], skiprows=1)
        mask_right = np.loadtxt(args.mask[1], skiprows=1)
        mask = {'left': mask_left,
                'right': mask_right}
        # check_data_header(mask)
    else:
        if args.mask_value is not None:
            mask = {'left': data_left != args.mask_value,
                    'right': data_right != args.mask_value}
        else:
            mask = {'left': np.ones_like(data_left),
                    'right': np.ones_like(data_right)}
    
    mask['left'][np.isnan(data_left)] = 0
    mask['right'][np.isnan(data_right)] = 0

    mask['left'] =  mask['left'].astype(bool)
    mask['right'] =  mask['right'].astype(bool)

    plot_surface.plot_surface(data, args.output, surface=surface, vlim=args.vlim, mask=mask, cbar_loc=args.cbar_loc, 
                              cbar_title=args.cbar_title, title=args.title, cmap=args.cmap, views=views, clobber=args.clobber)

# def check_data_header(data):
#     for hemisphere in ['left', 'right']:
#         if len(data[hemisphere]) == N_VERTEX[hemisphere]+1:
#             data[hemisphere] = data[hemisphere][1:] # Remove header line
#         elif len(data[hemisphere]) != N_VERTEX[hemisphere]:
#             raise ImportError(f'Data shape does not fit number of vertices. Should be left: {N_VERTEX["left"]}, right: {N_VERTEX["right"]}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Plot data')
    parser.add_argument('surface_left', help='Location of left surface.')
    parser.add_argument('surface_right', help='Location of right surface.')
    parser.add_argument('data_left', help='Location of data to plot on left hemisphere.')
    parser.add_argument('data_right', help='Location of data to plot on right hemisphere.')
    parser.add_argument('output', help='Location of output (.jpg or .png)')
    parser.add_argument('-vlim', default=None, nargs=2, metavar=('v_min', 'v_max'), type=float, help='Value limits on the plot.')
    parser.add_argument('-mask', default=None, nargs=2, metavar=('left_mask', 'right_mask'), help='Location of mask for left and right hemisphere')
    parser.add_argument('-mask_value', default=None, type=float, help='Mask values on surface. Can be float or nan. If mask is set, this argument is ignored')
    parser.add_argument('-cbar_loc', default='left', choices=['left', 'bottom', 'none'], help='Location of colorbar. If None, no colorbar is displayed')
    parser.add_argument('-cbar_title', default=None, type=str, help='Title of the colorbar')
    parser.add_argument('-title', default=None, type=str, help='Title of the plot')
    parser.add_argument('-cmap', default='turbo', help='Colormap of the plot. Can be any cmap from matplotlib')
    parser.add_argument('-views', default='standard', choices=['compact', 'standard', 'complete'], help='Option to select which views to show the surfaces.')
    parser.add_argument('-clobber', action='store_true', help='Overwrite existing files.')

    args = parser.parse_args()

    # class args():
    #     surface_left = '/Users/au483096/data/surface/mni_icbm152_t1_tal_nlin_sym_09c_left_smooth.obj'
    #     surface_right = '/Users/au483096/data/surface/mni_icbm152_t1_tal_nlin_sym_09c_right_smooth.obj'
    #     data_left = '/Volumes/projects/MINDLAB2013_18-MR-Amnestic-MCI/scratch/surface_data/0049/bl/0049_20150626_082005_mid_left_SEPWI_RTH_std_blur20.dat'
    #     data_right = '/Volumes/projects/MINDLAB2013_18-MR-Amnestic-MCI/scratch/surface_data/0049/bl/0049_20150626_082005_mid_right_SEPWI_RTH_std_blur20.dat'
    #     output = '/Users/au483096/Desktop/test.jpg'
    #     vlim = None
    #     mask = None
    #     mask_value = None
    #     cbar_loc = None
    #     cbar_title = None
    #     title = None
    #     cmap = 'turbo'
    #     views = 'standard'
    #     clobber = False

    # Convert none string to None
    if args.cbar_loc == 'none':
        args.cbar_loc = None

    main(args)


    
