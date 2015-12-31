#! /usr/bin/env python

from astropy.io import ascii 
import numpy as np
import multiprocessing as mp
import argparse
import os 
import sys  

import quick_look 

 

def parse_args():
    '''
    Parse command line arguments.  Returns args object.
    '''
    parser = argparse.ArgumentParser(description="makes quicklooks for everything in 'targets.list' file")
    parser.add_argument('targets', metavar='list', type=str, action='store',
                        help="""targets.list is the file to be read in;
                              first column = flag (0,1) if target is to be used,
                              second column = target/directory name""")

    parser.add_argument('--clobber', dest='clobber', action='store_true')
    parser.add_argument('--no-clobber', dest='clobber', action='store_false', help="default is no clobber")
    parser.set_defaults(clobber=False)

    args = parser.parse_args()
    return args

#-----------------------------------------------------------------------------------------------------

def drive_quick_look(targets): 
    clobber = targets[1]
    canonical = ascii.read(targets[0]+'.list') 
    print canonical
    dirlist = canonical['targname'][np.where(canonical['flag'] == 1)]
    print "dirlist in driver" 
    print dirlist 
    print 'clobber ', clobber 
    dirc = []
    for d in dirlist:
        dirc.append((d, clobber))

    mp_drive_ql_dirlist(dirc)

#-----------------------------------------------------------------------------------------------------

def drive_ql_dirlist(dirc):
    dirname = dirc[0]
    clobber = dirc[1]
    print "Driving target:  ", dirname
    if (os.path.isdir(dirname)): 
        os.chdir(dirname) 
        # PUT ANYTHING THAT WILL BE DONE IN TARGET DIRECTORY HERE 
        filelist = os.listdir('.') 
        print dirname, ":  ", filelist 
        if(clobber != 1 and os.path.exists(dirname+'_quicklook.html')): 
            print dirname+ ":  Quicklook already exists, skipping  "
        else: 
            print dirname+":  Creating Quick Look for ", dirname  
        quick_look.get_quick_look() 
        # NOW DONE WITH EVERYTHING THAT WILL BE DONE IN TARGET DIRECTORY 
    os.chdir('..')          # go back to "datapile" 


#-----------------------------------------------------------------------------------------------------

def mp_drive_ql_dirlist(dirc):
    pool = mp.Pool(processes=6)
    pool.map(drive_ql_dirlist, dirc)


#-----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()
    targets = (args.targets, args.clobber)
    
    drive_quick_look(targets)
    sys.exit("""
    
    ~~~~~~~*~*~*~*~
    ~~~~~~~*~*~*~*~  all done!!!! spectra are fun!
    ~~~~~~~*~*~*~*~
    """)


