import sys
import os
sys.path.append('/home/pearsonw/dev/XID_plus/') #CHANGE ME!
import xidplus
from xidplus import extra_functions as extra

output_folder = '/home/pearsonw/dev/XID_plus/output/'
tile_file_name = 'MIPS_24_vanilla_'
Master_filename = 'Master_prior_vanilla.pkl'
instrument = 'MIPS'     #SPIRE or PACS or MIPS

start_tile = 0
end_tile = 1000

if instrument == 'SPIRE':
    extra.make_master_SPIRE_post_catalogue_HEALpix(output_folder, Master_filename, tile_file_name, start_tile=start_tile, end_tile=end_tile)
elif instrument == 'PACS':
    extra.make_master_PACS_post_catalogue_HEALpix(output_folder, Master_filename, tile_file_name, start_tile=start_tile, end_tile=end_tile)
elif instrument == 'MIPS':
    extra.make_master_MIPS_post_catalogue_HEALpix(output_folder, Master_filename, tile_file_name, start_tile=start_tile, end_tile=end_tile)

