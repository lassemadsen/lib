from visbrain.objects import BrainObj
from visbrain.gui import Brain

# from os.path import expanduser
# surface_dir = f'{expanduser("~")}/data/surface'

# SURFACE_MNI = {'left': f'{surface_dir}/mni_icbm152_t1_tal_nlin_sym_09_left_smooth.gii',
#                'right': f'{surface_dir}/mni_icbm152_t1_tal_nlin_sym_09c_right_smooth.gii',
#                'both': f'{surface_dir}/mni_icbm152_t1_tal_nlin_sym_09c_both_smooth.gii'}

b_obj = BrainObj('/Users/au483096/data/surface/mni_icbm152_t1_tal_nlin_sym_09c_both_smooth.obj')

vb = Brain(brain_obj=b_obj)
vb.show()
