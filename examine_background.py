from astropy.io import fits
from astropy.table import Table, Column, vstack 
from astropy.io import ascii
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors
import glob
import os 

def examine_background(canonical_filename): 
	
        canonical = ascii.read(canonical_filename+'.canonical') 
        dirlist = canonical['col2'][np.where(canonical['col1'] == 1)]

	exposures = Table(names=('Rootname','Grating', 'Mean Flux', 'Mean Net', 'Mean Back'), dtype=('S20','S10','f10','f10','f10')) 

        number = 1 

 	for dirname in dirlist: 
                if (os.path.isdir(dirname)): 
			os.chdir(dirname) 

			filelist = glob.glob(os.path.join('.', '*x1d.fits'))

			#print "There are ", np.size(filelist), " exposures for target ", dirname 

			for filename in filelist: 
        		    #### Grab the first file and create html page for the "sample" 
		            hdulist = fits.open(filelist[0])
		            hdr0 = hdulist[0].header
		            hdr1 = hdulist[1].header

		            targname = hdr0['TARGNAME']

			    data = hdulist[1].data 

			    wave_a = data['WAVELENGTH'][0] 
			    flux_a = data['FLUX'][0] 
			    back_a = data['BACKGROUND'][0] 
			    net_a  = data['NET'][0] 

		 	    print filelist[0], hdr0['OPT_ELEM'], 'mean flux, net, bck :', np.mean(flux_a[7200:8200]), np.mean(net_a[7200:8200]), np.mean(back_a[7200:8200]) 
			    exposures.add_row([filelist[0], hdr0['OPT_ELEM'], np.mean(flux_a[7200:8200]), np.mean(net_a[7200:8200]), np.mean(back_a[7200:8200])] ) 

			os.chdir('..')

	exposures.write('exposures_background.fits', format='fits', overwrite=True) 


