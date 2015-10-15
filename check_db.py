from astropy.io import fits
from astropy.table import Table, Column, vstack 
from astropy.io import ascii
import numpy as np
import glob
import os 
import sys  
import copy 

def check_db(): 

    sample = Table.read('/grp/hst/HST_spectro/datapile/all_targets_sample.fits')
    websample = Table.read('/grp/hst/HST_spectro/datapile/all_targets_websample.fits') 

    canonical_files = ['STAR_EARLY', 'STAR_LATE', 'STAR_DWARF_NOVAE', 'STAR_TTAURI', 'STAR_POST_AGB', 'STAR_WD', 'STAR_LMXB',\
	'solar_system', 'QSOALS', 'GALAXY_STARBURST','GALAXY_SPIRAL','GALAXY_STARFORMING', 'GALAXY_DWARF_COMPACT','GALAXY_EMISSION_LINE', 'GALAXY_IRREGULAR', 'SUPERNOVAE'] 

    canonical_sum = ascii.read(canonical_files[0]+'.canonical') 

    for cfile in canonical_files[1:]: 
	print cfile+'.canonical' 
        c = ascii.read(cfile+'.canonical') 
        canonical_sum = vstack([canonical_sum, c])  

    all_canonical = ascii.read('/grp/hst/HST_spectro/datapile/all_targets.canonical') 
    copy_canonical = copy.deepcopy(all_canonical) 

    for ii in np.arange(np.size(all_canonical)): 
        if all_canonical['Target Name'][ii] in canonical_sum['col2']: 
	    print 'yes', all_canonical['Target Name'][ii] 
            jj = np.where(copy_canonical['Target Name'] == all_canonical['Target Name'][ii])[0] 
            #print jj.dtype 
            #print np.size(jj[0]) 
	    #help(jj[0]) 
            copy_canonical.remove_row(jj[0]) 
	else: 
	    print 'no', all_canonical[ii] 

    copy_canonical.write('unused_targets.html', format='jsviewer') 
    copy_canonical.write('UNUSED.temp',format='ascii.fixed_width', delimiter=',') 

    os.system('cat UNUSED.temp | grep -v arget |  cut -b 4-121 > UNUSED.canonical') 
    os.system('rm UNUSED.temp')

    return copy_canonical 
