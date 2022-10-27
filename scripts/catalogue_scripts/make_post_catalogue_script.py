import sys
import os
sys.path.append('/home/pearsonw/dev/XID_plus/') #CHANGE ME!
import xidplus
from xidplus import extra_functions as extra

if len(sys.argv) == 1:
    start_tile = 0
    end_tile = 0
elif len(sys.argv) < 3 or len(sys.argv) > 3:
    print('sys.argv is the wrong size')
    print('call as:')
    print('python make_post_catalogue_script.py <start_tile> <end_tile>')
    exit()
else:
    start_tile = int(sys.argv[1])
    end_tile = int(sys.argv[2])


output_folder = '/home/pearsonw/dev/XID_plus/output_numpyro/'
tile_file_name = 'MIPS_24_gaussian_'
Master_filename = 'Master_prior_gaussian.pkl'
instrument = 'MIPS'     #SPIRE or PACS or MIPS

if instrument == 'SPIRE':
    extra.make_master_SPIRE_post_catalogue_HEALpix(output_folder, Master_filename, tile_file_name, start_tile=start_tile, end_tile=end_tile)
elif instrument == 'PACS':
    extra.make_master_PACS_post_catalogue_HEALpix(output_folder, Master_filename, tile_file_name, start_tile=start_tile, end_tile=end_tile)
elif instrument == 'MIPS':
    extra.make_master_MIPS_post_catalogue_HEALpix(output_folder, Master_filename, tile_file_name, start_tile=start_tile, end_tile=end_tile)

