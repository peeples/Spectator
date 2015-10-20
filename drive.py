#! /usr/bin/env python

from astropy.io import ascii
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import os
import multiprocessing as mp

import quick_look
import add
import sys


def drive_targets(targets):
    clobber = targets[1]
    canonical = ascii.read(targets[0]+'.canonical')
    print canonical['col1']
    print canonical['col2']
    dirlist = canonical['col2'][np.where(canonical['col1'] == 1)]
    print "dirlist in driver"
    print dirlist
    print 'clobber ', clobber
    dirc = []
    for d in dirlist:
        dirc.append((d,clobber))

    #### need to rewrite this for loop to its own function to parallelize .....
    mp_drive_dirlist(dirc)

#-----------------------------------------------------------------------------------------------------

def drive_dirlist(dirc):
    dirname = dirc[0]
    clobber = dirc[1]
    print "Driving target:  ", dirname, "with clobber = ",clobber
    if os.path.isdir(dirname):
        os.chdir(dirname)
        # PUT ANYTHING THAT WILL BE DONE IN TARGET DIRECTORY HERE
        filelist = os.listdir('.')
            #if(clobber != 1 and os.path.exists(dirname+'_FUV_M_coadd.dat') or os.path.exists(dirname+'_FUV_L_coadd.dat')):
            #    print dirname+ ":  Coadd already exists, skipping  "
            #else:
            #    print dirname+ ":  Running coadd_x1d"
            #    os.system('/Applications/exelis/idl82/bin/idl -e "@../../code/coadd_script.pro"')

        if clobber != 1 and (os.path.exists(dirname+'_coadd_G130M_final_all.fits') or os.path.exists(dirname+'_coadd_G140L_final_all.fits')):
            print dirname+ ":  JT Coadd already exists, skipping  "
        else:
            print dirname+ ":  Running add"
            coadd = add.main('FUVM')
            coadd = add.main('FUVL')

        if clobber != 1 and os.path.exists(dirname+'_quicklook.html'):
            print dirname+ ":  Quicklook already exists, skipping  "
        else:
            print dirname+":  Creating Quick Look for ", dirname
            quick_look.get_quick_look()

#  			if(os.path.exists(filelist[0]+'_s_to_n.png')):
#  				print dirname+ ":  S/N estimates exist, skipping  "
#			else:
#				print "Assessing Line vs. Continuum for   ", dirname
#				line_vs_continuum.determine_line()


            # NOW DONE WITH EVERYTHING THAT WILL BE DONE IN TARGET DIRECTORY

        os.chdir('..')  		# go back to "datapile"


#-----------------------------------------------------------------------------------------------------

def mp_drive_targets(targets):
    pool = mp.Pool(processes=1)
    pool.map(drive_targets, targets)


#-----------------------------------------------------------------------------------------------------

def mp_drive_dirlist(dirc):
    pool = mp.Pool(processes=6)
    pool.map(drive_dirlist, dirc)


#-----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    targets = ( "qso_absorber" , 1 ) 
    clobber = 1
    
    drive_targets(targets)
    sys.exit("""
    
    ~~~~~~~*~*~*~*~
    ~~~~~~~*~*~*~*~  all done!!!! spectra are fun!
    ~~~~~~~*~*~*~*~
    """)
