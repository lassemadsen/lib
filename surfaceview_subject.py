#!/usr/local/bin/python3
""" Script to show surface of specific subject provided with data from given data file in the MCI study (MNI space)
"""
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)

import argparse
from visbrain.objects import BrainObj
from visbrain.gui import Brain
import pandas as pd
import numpy as np

def main(args):
    """
    Visualize parameter data MNI surface using
    the Brain gui from Visbrain. 

    Parameters in args
    ------------------

    Notes
    -----
    Works with numpy version 1.16.2
    Seems not to work with numpy version 1.16.4 for unknown reasons
    """

    if args.parameter in ['pib_suvr', 'pk_bp', 'ftp_suvr']:
        data_path = '/Users/au483096/Documents/PhD/BANIMCI/scratch/data/surface_data_PET'
    else:
        data_path = f'/Users/au483096/Documents/PhD/BANIMCI/scratch/data/surface_data_{args.perf_pipeline}'

    data = pd.read_csv(f'{data_path}/{args.parameter}_{args.tp}_{args.hemisphere}.csv', usecols=[args.subject])
    data = data.to_numpy().ravel()

    b_obj = BrainObj(f'mni_icbm152_t1_tal_nlin_sym_09c_both_smooth', translucent=False, hemisphere=args.hemisphere)

    b_obj.add_activation(data=data, hemisphere=args.hemisphere, cmap='jet')

    vb = Brain(brain_obj=b_obj, bgcolor='white')
    vb.cbar_control('brain', txtcolor='black', bgcolor='white')
    
    if args.clim is not None:
        vb.cbar_control('brain', clim=args.clim)
    else:    
        vb.cbar_autoscale('brain')

    vb.cbar_select('brain')

    vb.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-subject', help='Subject to show')
    parser.add_argument('-parameter', help='Parameter to show')
    parser.add_argument('-hemisphere', help='Hemisphere to show')
    parser.add_argument('-tp', help='Timepoint to show')
    parser.add_argument('-clim', nargs=2, default=None, help='Colorbar limits, if None: autoscale', type=float)
    parser.add_argument('-perf_pipeline', default='lsmDec20_longi', help='Name of MR perfusion data')
    args = parser.parse_args()
    main(args)