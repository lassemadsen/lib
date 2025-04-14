#!/usr/bin/env python
import warnings

warnings.simplefilter(action='ignore', category=UserWarning)

import argparse
import numpy as np
import pandas as pd
from visbrain.objects import BrainObj
from visbrain.gui import Brain
import sys

from os.path import expanduser
home = expanduser("~")

surface_dir = f'{home}/data/surface'

SURFACE_MNI = {'left': f'{surface_dir}/mni_icbm152_t1_tal_nlin_sym_09c_left_smooth.gii',
               'right': f'{surface_dir}/mni_icbm152_t1_tal_nlin_sym_09c_right_smooth.gii',
               'both': f'{surface_dir}/mni_icbm152_t1_tal_nlin_sym_09c_both_smooth.gii'}

N_VERTEX_MNI = {'left': 81349,
                'right': 81233}

SURFACE_FS = {'left': f'{surface_dir}/lh.pial.gii',
               'right': f'{surface_dir}/rh.pial.gii',
               'both': f'{surface_dir}/both.pial.gii'}

N_VERTEX_FS = {'left': 163842,
                'right': 163842}

def get_surface(nv):
    """Return the surface template, that fits with the number of vertices in the left and right hemisphere

    Parameters
    ----------
    nv : int
        Number of vertices
    """

    header = False

    # MNI surface
    if nv == N_VERTEX_MNI['left']:
        surface = SURFACE_MNI
        hemisphere = 'left'
    elif nv == N_VERTEX_MNI['left']+1:
        surface = SURFACE_MNI
        hemisphere = 'left'
        header = True
    elif nv == N_VERTEX_MNI['right']:
        surface = SURFACE_MNI
        hemisphere = 'right'
    elif nv == N_VERTEX_MNI['right']+1:
        surface = SURFACE_MNI
        hemisphere = 'right'
        header = True

    # Freesurfer
    elif nv == N_VERTEX_FS['left']:
        surface = SURFACE_FS
        hemisphere = None
    elif nv == N_VERTEX_FS['left']+1:
        surface = SURFACE_FS
        hemisphere = None
        header = True
    else:
        sys.exit('Error: Number of vertices do not match any surface...')

    return surface, hemisphere, header

def main(args):
    """
    Visualize mni_icbm152_t1_tal_nlin_sym_09c_both_smooth using
    the Brain gui from Visbrain. 
    If data is provided, this will be added as activation on the surface template. 

    Parameters in args
    ------------------
    -data : str
        The data to add as activation on the surface
        Only one data file can be added. Either left or right 
        hemisphere.
    -surface : str
        Basename of surface file located in visbrain_data/templates
        If not set, the defualt mni surface will be loaded.

    Notes
    -----
    Works with numpy version 1.16.2
    Seems not to work with numpy version 1.16.4 for unknown reasons
    """

    # ===== If data is provided, this is loaded ===== 
    if args.surface is None:
        if args.data is None:
            surface, _, _ = get_surface(81349) # If no input args, show MNI surface
            b_obj = BrainObj(surface['both'], translucent=False, hemisphere='both')
        else:
            try:
                data = pd.read_csv(args.data)
                data = data.replace(np.nan, -1)
                data = np.array(data.iloc[:,0])
            except :
                print('Error: Data file could not be loaded')
                exit()

            n_vert = len(data)

            surface, hemisphere, header = get_surface(n_vert)

            if hemisphere is None:
                if args.hemis is None:
                    sys.exit('Error. When plotting freesurfer data, -hemisphere needs to be given')
                else:
                    hemisphere = args.hemis

            if header:
                data = data[1:] # Remove header

            b_obj = BrainObj(surface[hemisphere], translucent=False)
            b_obj.add_activation(data=data, cmap='jet')

    else:
        # TODO update this to automatically convert surface 
        b_obj = BrainObj(args.surface, translucent=False, hemisphere='both')

        if args.data is not None:
            try:
                data = np.loadtxt(args.data)
            except:
                print('Error: Data file could not be loaded')
                exit()

            if len(b_obj.vertices) != len(data):
                # Check for n_vertices header (first line) and remove if present            
                if len(data) == len(b_obj.vertices) + 1:    
                    data = data[1:]
                else:
                    print('Dimensions of data not matching dimensions of template...')
                    print('Data should have dimensions: (%s,)' % len(b_obj.vertices))
                    print('Dimensions given: ' + str(data.shape))
                    exit()

            b_obj.add_activation(data=data, hemisphere='both', cmap='jet')
    
    # ===== Create Brain gui ===== 
    vb = Brain(brain_obj=b_obj, bgcolor='white')
    if args.data is not None:
        vb.cbar_control('brain', txtcolor='black', bgcolor='white')
        vb.cbar_autoscale('brain')
        vb.cbar_select('brain')

    vb.show()
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-data', help='Datafile to overlay on surface map')
    parser.add_argument('-hemis', help='When plotting freesurfer, hemisphere needs to be given, as it cannot be implied from the data')
    parser.add_argument('-surface', help='Basename of *.npz file containing surface file. Should be in the visbrain/template folder. Default is left/right MNI surface in visbrain folder')
    args = parser.parse_args()
    main(args)
