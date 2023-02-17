import sys
import os
sys.path.append('/home/pearsonw/dev/XID_plus/') #CHANGE ME!
import xidplus
from xidplus import extra_functions as extra
from astropy.io import fits
from astropy import wcs

output_folder = '/mnt/disk1/sky_maps/XID_plus_output/'
tile_file_name = 'SPIRE_vanilla_'
Master_filename = 'Master_SPIRE_prior_vanilla.pkl'
band = 'psw'    #mips 24 - mips24
                #pacs 100 - green
                #pacs 160 - red
                #spire 250 - psw
                #spire 350 - pmw
                #spire 500 - plw

imfolder='/mnt/disk1/sky_maps/SIDES/SPIRE/'
origional_map = 'pySIDES_from_original_SPIRE250_smoothed_Jy_beam.fits'

hdulist = fits.open(imfolder+origional_map)
im = hdulist[0].data #Do not scale the map here...
imhdu = hdulist[0].header
w_pri = wcs.WCS(hdulist[0].header)
hdulist.close()

extra.make_residual_map_HEALpix(output_folder, Master_filename, tile_file_name, band, im, imhdu, w_pri, scale=1e3) #...scale the map with scale= here
#Optional parameters:
#   scale      - value used to scale the image when XID+ was run, default=1.
#   start_tile - tile to start from, default=0
#   noise      - add noise to the replicated image, default=False
#   sample     - posterior sample to build map for, default=None (50th percentile values will be used)
