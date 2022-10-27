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
imfolder='/home/pearsonw/dev/XID_plus/sky_maps/PACS/' #'./input/'
#field
field='COSMOS'
#SMAP version
SMAPv='6.0'


# In[7]:

pswfits=imfolder+'COSMOS-Nest_image_100_SMAP_v6.0.fits' #PACS 100 map
pmwfits=imfolder+'COSMOS-Nest_image_160_SMAP_v6.0.fits' #PACS 160 map


#----output folder-----------------
output_folder='/home/pearsonw/dev/XID_plus/output/'

#Folder containing prior input catalogue
folder='/home/pearsonw/dev/XID_plus/input/'
#prior catalogue
prior_cat='COSMOS2015_CIGALE_MC_GA_DR_x_PACS.fits'

hdulist = fits.open(folder+prior_cat)
fcat=hdulist[1].data
hdulist.close()
inra=fcat['RA']
indec=fcat['Dec']
f_src=fcat['PSW']
df_src=f_src
nrealcat=fcat.size
bkg100=-5
bkg160=-5


ids = fcat['id']

# Open images and noise maps and use WCS module in astropy to get header information

# In[9]:

#-----100-------------
hdulist = fits.open(pswfits)
im100phdu=hdulist[0].header
im100hdu=hdulist[1].header

im100=hdulist[1].data*1.0E3
nim100=hdulist[2].data*1.0E3
w_100 = wcs.WCS(hdulist[1].header)
pixsize100=3600.0*w_100.wcs.cd[1,1] #pixel size (in arcseconds)
hdulist.close()
#-----160-------------
hdulist = fits.open(pmwfits)
im160phdu=hdulist[0].header
im160hdu=hdulist[1].header

im160=hdulist[1].data*1.0E3
nim160=hdulist[2].data*1.0E3
w_160 = wcs.WCS(hdulist[1].header)
pixsize160=3600.0*w_160.wcs.cd[1,1] #pixel size (in arcseconds)
hdulist.close()


#only use all sources with PSW > 0.7
sgood=f_src > 0.7

inra=inra[sgood]
indec=indec[sgood]
n_src=sgood.sum()



# Point response information, at the moment its 2D Gaussian,

#pixsize array (size of pixels in arcseconds)
pixsize=np.array([pixsize100,pixsize160])
#point response function for the three bands
prfsize=np.array([9.76,13.51])
#use Gaussian2DKernel to create prf (requires stddev rather than fwhm hence pfwhm/2.355)
from astropy.convolution import Gaussian2DKernel


#Set prior classes
#---prior100--------
prior100=xidplus.prior(im100,nim100,im100phdu,im100hdu) #Initialise with map, uncertianty map, wcs info and primary header
prior100.prior_cat(inra,indec,prior_cat,ID=ids[sgood]) #Set input catalogue
prior100.prior_bkg(bkg100,5) #Set prior on background
#---prior160--------
prior160=xidplus.prior(im160,nim160,im160phdu,im160hdu)
prior160.prior_cat(inra,indec,prior_cat,ID=ids[sgood])
prior160.prior_bkg(bkg160,5)


print('fitting '+ str(prior160.nsrc)+' sources \n')


##---------fit using Gaussian beam-----------------------
prf100=Gaussian2DKernel(prfsize[0]/2.355,x_size=101,y_size=101)
prf100.normalize(mode='peak')
prf160=Gaussian2DKernel(prfsize[1]/2.355,x_size=101,y_size=101)
prf160.normalize(mode='peak')

pind100=np.arange(0,101,1)*1.0/pixsize[0] #get 100 scale in terms of pixel scale of map
pind160=np.arange(0,101,1)*1.0/pixsize[1] #get 160 scale in terms of pixel scale of map

prior100.set_prf(prf100.array,pind100,pind100)
prior160.set_prf(prf160.array,pind160,pind160)

#from moc, get healpix pixels at a given order
from xidplus import moc_routines
order=11
tiles=moc_routines.get_HEALPix_pixels(order,prior100.sra,prior100.sdec,unique=True)

try:
    if sys.argv[1] == 'Master':
        print('----- There are '+str(len(tiles))+' tiles required for input catalogue')
        outfile=output_folder+'Master_prior_vanilla.pkl'
        with open(outfile, 'wb') as f:
            pickle.dump({'psw':prior100,'pmw':prior160,'tiles':tiles,'order':order},f)
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
prior100.set_tile(moc)
prior160.set_tile(moc)

print('fitting '+ str(prior100.nsrc)+' sources \n')
print('there are '+ str(prior100.snpix)+' pixels')

prior100.get_pointing_matrix()
prior160.get_pointing_matrix()

print('set prior flux scale')
prior100.upper_lim_map()
prior100.lower_lim_flux(0.0)
prior160.upper_lim_map()
prior160.lower_lim_flux(0.0)


#Run XID+
from xidplus.stan_fit import PACS
fit=PACS.all_bands(prior100,prior160,iter=1500)
posterior=xidplus.posterior_stan(fit,[prior100,prior160])
#Save results
outfile=output_folder+'PACS_vanilla_'+str(tiles[taskid-1])+'_'+str(order)+'.pkl'
with open(outfile, 'wb') as f:
   pickle.dump({'psw':prior100,'pmw':prior160,'posterior':posterior},f)

