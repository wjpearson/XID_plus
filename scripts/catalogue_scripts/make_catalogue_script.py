import sys
import os
sys.path.append('/home/pearsonw/dev/XID_plus/') #CHANGE ME!
import xidplus
from xidplus import extra_functions as extra

output_folder = '/home/pearsonw/dev/XID_plus/output_numpyro/'
tile_file_name = 'MIPS_24_gaussian_'
Master_filename = 'Master_prior_gaussian.pkl'
instrument = 'MIPS' #SPIRE or PACS or MIPS

if instrument == 'SPIRE':
    extra.make_master_SPIRE_catalogue_HEALpix(output_folder, Master_filename, tile_file_name)
elif instrument == 'PACS':
    extra.make_master_PACS_catalogue_HEALpix(output_folder, Master_filename, tile_file_name)
elif instrument == 'MIPS':
    extra.make_master_MIPS_catalogue_HEALpix(output_folder, Master_filename, tile_file_name)

