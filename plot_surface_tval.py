#!/usr/bin/env python
import argparse
from surface_plot import plot_stats
import numpy as np

from surface_plot.config import N_VERTEX

def main(args):
    # Load data
    tval_left = np.loadtxt(args.tval_left)
    tval_right = np.loadtxt(args.tval_right)

    tval_data = {'left': tval_left,
                 'right': tval_right}
    check_data_header(tval_data)

    # --- P-values ---
    if args.pval is not None:
        pval_left = np.loadtxt(args.pval[0])
        pval_right = np.loadtxt(args.pval[1])

        pval_data = {'left': pval_left,
                    'right': pval_right}
        check_data_header(pval_data)
    else:
        pval_data = None

    # --- Second threshold mask --- 
    if args.second_threshold_mask is not None:
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
                        mask = {'left': ~np.isnan(tval_left),
                                'right': ~np.isnan(tval_right)}
            else:
                mask = {'left': tval_left!=args.mask_value,
                        'right': tval_right!=args.mask_value}
        else:
            mask = None
    
    if args.one_tailed:
        two_tailed = False
    else:
        two_tailed = True
    
    if args.narrow_edge:
        expand_edge = False
    else:
        expand_edge = True

    plot_stats.plot_tval(tval_data, args.output, t_lim=args.t_lim, t_threshold=args.t_threshold, mask=mask,
                         p_threshold=args.p_threshold, pval=pval_data, df=args.df, two_tailed=two_tailed,
                         title=args.title, cbar_loc=args.cbar_loc, second_threshold_mask=second_threshold_mask,
                         expand_edge=expand_edge, clobber=args.clobber)

def check_data_header(data):
    for hemisphere in ['left', 'right']:
        if len(data[hemisphere]) == N_VERTEX[hemisphere]+1:
            data[hemisphere] = data[hemisphere][1:] # Remove header line
        elif len(data[hemisphere]) != N_VERTEX[hemisphere]:
            raise ImportError(f'Data shape does not fit number of vertices. Should be left: {N_VERTEX["left"]}, right: {N_VERTEX["right"]}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Plot data')
    parser.add_argument('tval_left', help='Location of tval data to plot on left hemisphere.')
    parser.add_argument('tval_right', help='Location of tval data to plot on right hemisphere.')
    parser.add_argument('output', help='Location of output')
    parser.add_argument('-t_lim', default=None, nargs=2, metavar=('t_min', 't_max'), type=float, help='Value limits for the t-values.')
    parser.add_argument('-t_threshold', default=2, type=float, help='Statistical threshold of t-values')
    parser.add_argument('-p_threshold', default=None, type=float, help='Statistical threshold of t-values given in p-values. If this is set, t_threshold is ignored. Requires either pval or df')
    parser.add_argument('-pval', default=None, nargs=2, metavar=('left_pval', 'right_pval'), help='Location of pval data to plot on left and right hemisphere. Only used if p_threshold is set.')
    parser.add_argument('-df', default=None, type=int, help='Degrees of freedom of the data. Only used if -p_threshold is set and -pval is not')
    parser.add_argument('-one_tailed', action='store_true', help='One-sided t-values')
    parser.add_argument('-mask', default=None, nargs=2, metavar=('left_mask', 'right_mask'), help='Location of mask for left and right hemisphere')
    parser.add_argument('-mask_value', default=None, type=float, help='Mask values on surface. Can be float or nan. If mask is set, this argument is ignored')
    parser.add_argument('-cbar_loc', default='left', choices=['left', 'bottom', 'none'], help='Location of colorbar. If None, no colorbar is displayed')
    parser.add_argument('-title', default=None, help='Title of the plot')
    parser.add_argument('-second_threshold_mask', default=None, nargs=2, metavar=('left_mask', 'right_mask'), help='Non-zero values of mask are surrounded by a white line one plot')
    parser.add_argument('-narrow_edge', action='store_true', help='If set, the second_threshold_mask is only surrounded by one vertex instead of two as default')
    parser.add_argument('-clobber', action='store_true', help='Overwrite existing files.')

    args = parser.parse_args()

    if args.p_threshold is not None and (args.pval is None and args.df is None):
        parser.error('The -p_threshold argument requires either the -pval or -df arguments')
    elif args.p_threshold is not None and (args.pval is not None and args.df is not None):
        print('Warning: Both pval and df is given. Only pval files will be used, df is ignored.')

    # Convert none string to None
    if args.cbar_loc == 'none':
        args.cbar_loc = None

    main(args)

