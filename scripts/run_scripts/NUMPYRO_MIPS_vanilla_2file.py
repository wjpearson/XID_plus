import numpy as np
from astropy.io import fits
from astropy import wcs
import pickle
import dill
import sys
import os
sys.path.append('/home/pearsonw/dev/XID_plus/') #CHANGE ME!
import xidplus

'''Original Script written by Peter Hurley
    Adapted by William J. Pearson'''

#Folder containing maps
imfolder='/home/pearsonw/dev/XID_plus/sky_maps/MIPS/' #'./input/'
#field
field='COSMOS'
#SMAP version
SMAPv='6.0'


# In[7]:

pswfits=imfolder+'pySIDES_from_original_MIPS24_smoothed_Jy_beam.fits' #MIPS 24 map
pswnois=imfolder+'mips_24_GO3_pySIDES.fits' #MIPS 24 noise map


#----output folder-----------------
output_folder='/home/pearsonw/dev/XID_plus/output_numpyro/'

#Folder containing prior input catalogue
folder='/home/pearsonw/dev/XID_plus/input/'
#prior catalogue
prior_cat='pySIDES_from_original.fits'

hdulist = fits.open(folder+prior_cat)
fcat=hdulist[1].data
hdulist.close()
inra=fcat['RA']
indec=fcat['Dec']
f_src=fcat['qflag']
df_src=f_src
nrealcat=fcat.size
bkg24=-5


ids = fcat['id']
PSW_mu = fcat['SMIPS24']*1e3 # * hasIRAC
PSW_sigma = (fcat['S24']*1e3/10) * 2

# Open images and noise maps and use WCS module in astropy to get header information

# In[9]:

#------24-------------
hdulist = fits.open(pswfits)
im24phdu=hdulist[0].header
im24hdu=hdulist[0].header

im24=hdulist[0].data *1.0E3
w_24 = wcs.WCS(hdulist[0].header)
pixsize24=3600.0*w_24.wcs.cdelt[1] #pixel size (in arcseconds)
hdulist.close()

hdulist = fits.open(pswnois)
nim24=hdulist[0].data *1.0E3
hdulist.close()


#only use all aources with predicted PSW > 0.7
sgood=f_src == False
sgood1=PSW_mu > 0.005
sgood *= sgood1

inra=inra[sgood]
indec=indec[sgood]
n_src=sgood.sum()



# Point response information, at the moment its 2D Gaussian,

#pixsize array (size of pixels in arcseconds)
pixsize=np.array([pixsize24])
#point response function for the three bands
prfsize=np.array([5.9])
#use Gaussian2DKernel to create prf (requires stddev rather than fwhm hence pfwhm/2.355)
from astropy.convolution import Gaussian2DKernel


#Set prior classes
#---prior24---------
prior24=xidplus.prior(im24,nim24,im24phdu,im24hdu) #Initialise with map, uncertianty map, wcs info and primary header
prior24.prior_cat(inra,indec,prior_cat,ID=ids[sgood]) #Set input catalogue
prior24.prior_bkg(bkg24,5) #Set prior on background

print('fitting '+ str(prior24.nsrc)+' sources \n')


##---------fit using Gaussian beam-----------------------
prf24=Gaussian2DKernel(prfsize[0]/2.355,x_size=101,y_size=101)
prf24.normalize(mode='peak')

pind24=np.arange(0,101,1)*1.0/pixsize[0] #get 24 scale in terms of pixel scale of map

prior24.set_prf(prf24.array,pind24,pind24)

#from moc, get healpix pixels at a given order
from xidplus import moc_routines
order=11
tiles=moc_routines.get_HEALPix_pixels(order,prior24.sra,prior24.sdec,unique=True)

try:
    if sys.argv[1] == 'Master':
        print('----- There are '+str(len(tiles))+' tiles required for input catalogue')
        outfile=output_folder+'Master_prior_vanilla.pkl'
        with open(outfile, 'wb') as f:
            pickle.dump({'psw':prior24,'tiles':tiles,'order':order},f)
        raise SystemExit()


except:
    pass

try:
    #This is setup to be run on a HPC running SLURM scheduler
    taskid = np.int(os.environ['SLURM_ARRAY_TASK_ID'])
    #This is the setup for what ever they use on their HPC in Sussex
    #task_first=np.int(os.environ['SGE_TASK_FIRST'])
    #task_last=np.int(os.environ['SGE_TASK_LAST'])

except KeyError:
    #If not running on a HPC, you can choose the tile here
    print("Error: could not read SLURM_ARRAY_TASK_ID from environment")
    taskid = int(input("Please enter task id: "))
    print("you entered", taskid)


moc=moc_routines.get_fitting_region(order,tiles[taskid-1])
prior24.set_tile(moc)

print('fitting '+ str(prior24.nsrc)+' sources \n')
print('there are '+ str(prior24.snpix)+' pixels')

prior24.get_pointing_matrix()

print('set prior flux scale')
prior24.upper_lim_map()
prior24.lower_lim_flux(0.0)

#RUN XID+
from xidplus.numpyro_fit import MIPS
fit=MIPS.MIPS_24(prior24,num_samples=750,num_warmup=750)
posterior=xidplus.posterior_numpyro(fit,[prior24])
#Save output
outfile=output_folder+'MIPS_24_vanilla_'+str(tiles[taskid-1])+'_'+str(order)+'.pkl'
with open(outfile, 'wb') as f:
   pickle.dump({'psw':prior24,'posterior':posterior},f)

