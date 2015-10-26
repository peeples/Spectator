from astropy.io import ascii
import os 
import numpy as np 

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
