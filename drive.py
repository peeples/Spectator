
from astropy.io import ascii
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import os

import quick_look
#import line_vs_continuum
import add

def drive_targets(canonical_filename, clobber):

    canonical = ascii.read(canonical_filename+'.canonical')
    print canonical['col1']
    print canonical['col2']
    dirlist = canonical['col2'][np.where(canonical['col1'] == 1)]
    print "dirlist in driver"
    print dirlist
    print 'clobber ', clobber

    for dirname in dirlist:
        print "Driving target:  ", dirname
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

if __name__ == "__main__":
    drive_targets()
    sys.exit("""
    ~~~~~~~*~*~*~*~
    ~~~~~~~*~*~*~*~  all done!!!! spectra are fun!
    ~~~~~~~*~*~*~*~""")
