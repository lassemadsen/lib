#!/usr/bin/env python
import argparse
from surface_plot import plot_stats
import numpy as np

from surface_plot.config import N_VERTEX

def main(args):
    # Load data
    pval_left = np.loadtxt(args.pval_left)
    pval_right = np.loadtxt(args.pval_right)

    pval_data = {'left': pval_left,
                 'right': pval_right}
    check_data_header(pval_data)

    # --- T-values ---
    if args.tval is not None:
        tval_left = np.loadtxt(args.tval[0])
        tval_right = np.loadtxt(args.tval[1])

        tval_data = {'left': tval_left,
                    'right': tval_right}
        check_data_header(tval_data)
    else:
        tval_data = None

    # --- Second threshold mask ---
    if args.second_threshold_mask is not None:
            # Load data
        second_threshold_mask_left = np.loadtxt(args.second_threshold_mask[0])
        second_threshold_mask_right = np.loadtxt(args.second_threshold_mask[1])

        second_threshold_mask = {'left': second_threshold_mask_left,
                                 'right': second_threshold_mask_right}
        check_data_header(second_threshold_mask)
    else:
        second_threshold_mask = None

    # --- Mask ---
    if args.mask is not None:
        mask = {'left': args.mask[0],
                'right': args.mask[1]}
    else:
        if args.mask_value is not None:
            if np.isnan(args.mask_value):
                        mask = {'left': ~np.isnan(pval_left),
                                'right': ~np.isnan(pval_right)}
            else:
                mask = {'left': pval_left!=args.mask_value,
                        'right': pval_right!=args.mask_value}
        else:
            mask = None
    
    if args.narrow_edge:
        expand_edge = False
    else:
        expand_edge = True

    plot_stats.plot_pval(pval_data, args.output, tval=tval_data, p_threshold=args.p_threshold, mask=mask,
                         cbar_loc=args.cbar_loc, titles=args.titles, second_threshold_mask=second_threshold_mask,
                         expand_edge=expand_edge, clobber=args.clobber)

def check_data_header(data):
    for hemisphere in ['left', 'right']:
        if len(data[hemisphere]) == N_VERTEX[hemisphere]+1:
            data[hemisphere] = data[hemisphere][1:] # Remove header line
        elif len(data[hemisphere]) != N_VERTEX[hemisphere]:
            raise ImportError(f'Data shape does not fit number of vertices. Should be left: {N_VERTEX["left"]}, right: {N_VERTEX["right"]}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Plot data')
    parser.add_argument('pval_left', help='Location of pval data to plot on left hemisphere.')
    parser.add_argument('pval_right', help='Location of pval data to plot on right hemisphere.')
    parser.add_argument('output', help='Location of output')
    parser.add_argument('-p_threshold', default=0.01, type=float, help='Statistical threshold of the plot.')
    parser.add_argument('-tval', default=None, nargs=2, metavar=('left_tval', 'right_tval'), help='Location of tval data. Used to define "positive" or "negative" p-values. If not given, p-values are plottet as they are.')
    parser.add_argument('-df', default=None, type=int, help='Degrees of freedom of the data. Only used if -p_threshold is set and -pval is not')
    parser.add_argument('-mask', default=None, nargs=2, metavar=('left_mask', 'right_mask'), help='Location of mask for left and right hemisphere')
    parser.add_argument('-mask_value', default=None, type=float, help='Mask values on surface. Can be float or nan. If mask is set, this argument is ignored')
    parser.add_argument('-cbar_loc', default='left', choices=['left', 'bottom', 'none'], help='Location of colorbar. If None, no colorbar is displayed')
    parser.add_argument('-titles', nargs=2, default=None, help='Title of the plots.')
    parser.add_argument('-second_threshold_mask', default=None, nargs=2, metavar=('left_mask', 'right_mask'), help='Non-zero values of mask are surrounded by a white line one plot')
    parser.add_argument('-narrow_edge', action='store_true', help='If set, the second_threshold_mask is only surrounded by one vertex instead of two as default')
    parser.add_argument('-clobber', action='store_true', help='Overwrite existing files.')

    args = parser.parse_args()

    # Convert none string to None
    if args.cbar_loc == 'none':
        args.cbar_loc = None

    main(args)

