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
import sys  
import subprocess

import quick_look 
import line_vs_continuum 
#import harvest_headers 
 

def sample_download(canonical_filename): 
	
    canonical = ascii.read(canonical_filename+'.canonical') 
    print canonical['col1']
    print canonical['col2']
    dirlist = canonical['col2'][np.where(canonical['col1'] == 1)]
    print "dirlist in driver"
    print dirlist 
    os.system('tar -cvf '+canonical_filename+'.tar '+canonical_filename+'.*') 

    for dirname in dirlist: 
        print "Driving target:  ", dirname
        if (os.path.isdir(dirname)): 
            print 'tar -rvf '+canonical_filename+'.tar '+dirname  
            os.system('tar -rvf '+canonical_filename+'.tar '+dirname) 

    os.system('gzip '+canonical_filename+'.tar') 
