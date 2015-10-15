from astropy.io import fits
from astropy.io import ascii 
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors
import glob
import os 
import sys  

import quick_look 
import line_vs_continuum 
import add 

mpl.rcParams['agg.path.chunksize'] = 10000 
 

def drive_quicklook(canonical_filename, clobber): 
	
	canonical = ascii.read(canonical_filename+'.canonical') 
	print canonical['col1'] 
	print canonical['col2'] 
 	dirlist = canonical['col2'][np.where(canonical['col1'] == 1)]
	print "dirlist in driver" 
	print dirlist 

	print 'clobber ', clobber 

 	for dirname in dirlist: 
		print "Driving target:  ", dirname
                if (os.path.isdir(dirname)): 
			os.chdir(dirname) 

		# PUT ANYTHING THAT WILL BE DONE IN TARGET DIRECTORY HERE 

			filelist = os.listdir('.') 
			print dirname, ":  ", filelist 
 
  			if(clobber !=1 and os.path.exists(dirname+'_quicklook.html')): 
  				print dirname+ ":  Quicklook already exists, skipping  "
			else: 
				print dirname+":  Creating Quick Look for ", dirname  
			quick_look.get_quick_look() 
  
		# NOW DONE WITH EVERYTHING THAT WILL BE DONE IN TARGET DIRECTORY 
 
			os.chdir('..')  		# go back to "datapile" 








