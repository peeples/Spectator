#! /usr/bin/env python

from astropy.io import fits
from astropy.table import Table, vstack 
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import argparse
import glob
import os 
import sys  

def parse_args():
    '''
    Parse command line arguments.  Returns args object.
    '''
    parser = argparse.ArgumentParser(description="scrapes headers for everything in 'filename.list' file")
    parser.add_argument('filename', metavar='canonical', type=str, action='store',
                        help='filename.list is the file to be read in')

    args = parser.parse_args()
    return args

#-----------------------------------------------------------------------------------------------------

def scrape_headers(targets): 
    canonical_filename = targets
    canonical = ascii.read(canonical_filename+'.list')

    print canonical
    dirlist = canonical['targname'][np.where(canonical['flag'] == 1)]
    print "dirlist in driver"
    print dirlist 

    sample_fitstable = Table(names=('Number','Target Name', 'RA','DEC','N_exposures','Target Category', \
        'Target Description','S/N (130)'), dtype=('i4','S200','f4','f4','I8','S20','S20','f1')) 

    sample_webtable = Table(names=('Number','Target Name', 'RA','DEC','N_exposures','Target Category', \
        'Target Description','Simbad', 'MAST', 'S/N (130)', 'FUV M','FUV L', 'FUV M Download (by LP)'), \
        dtype=('i4','S200','f4','f4','S240','S20','S20','S240','S240','f1','S250','S250','S350')) 

    targets = Table(names=('Flag','Target Name', 'Target Category', 'Target Description'), dtype=('i4','S200','S25', 'S350')) 

    # this will contain a global database of all exposures 
    exposures = Table(names=('Rootname','Target Name', 'RA','DEC','PropID','PI Name','Detector','Segment',\
        'LP','Grating', 'Cenwave','FPPOS','Exptime','Nevents','Extended','Date','Target Description'),   
        dtype=('S20','S35','f4','f4','I5','S20','S4','S5','S2','S10','S10','I2','f10','f8','S4','S12','S200'))


    #### set up the master "header table" 
    hdulist=fits.open('generic_x1d.fits')
    h0=hdulist[0].header
    for key in ["HISTORY", "COMMENT", ""]:
        del h0[key]
    names = sorted(h0.keys())
    rows = [h0[k] for k in names]
    header_table0 = Table(rows=[rows], names=names)

    h1=hdulist[1].header
    for key in ["HISTORY", "COMMENT", ""]:
        del h1[key]
    names = sorted(h1.keys())
    rows = [h1[k] for k in names]
    header_table1 = Table(rows=[rows], names=names)

    number = 1 

    total_number_of_headers = 0 

    for dirname in dirlist: 
        print "Driving target:  ", dirname
        if (os.path.isdir(dirname)): 
            os.chdir(dirname) 

            filelist = glob.glob(os.path.join('.', '*x1d.fits'))

            print "There are ", np.size(filelist), " exposures for target ", dirname 

            #### Grab the first file and create html page for the "sample_webtable" 
            hdulist = fits.open(filelist[0])
            hdr0 = hdulist[0].header
            hdr1 = hdulist[1].header

            targname = hdr0['TARGNAME']
            targdesc = hdr0['TARDESCR']
            ra = hdr0['RA_TARG']
            dec = hdr0['DEC_TARG']
            print ra, dec, targname 

            median_sn = -9.99 
        
            targname_urlstring = '<a href="'+targname+'/'+targname+'_quicklook.html">'+targname+'</a>'

            simbad_string = '<a href="http://simbad.u-strasbg.fr/simbad/sim-coo?CooDefinedFrames=none&CooEpoch=2000&Coord='+str(ra)+'d'+str(dec)+'d&submit=submit%20query&Radius.unit=arcsec&CooEqui=2000&CooFrame=FK5&Radius=4"> SIMBAD </a>'  

            mast_string = '<a href="https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery='+str(ra)+','+str(dec)+'"> MAST  </a>'  
            print mast_string 

            n_exp_string = '<a href="'+targname+'/all_exposures.html">'+str(np.size(filelist))+'</a>' 

            fuv_m_quicklook_urlstring = '...' 
            fuv_l_quicklook_urlstring = '...' 
            if (os.path.exists(hdr0['targname']+'_coadd_final_all.png')): 
                fuv_m_quicklook_urlstring = '<a href="'+targname+'/'+targname+'_coadd_G130M_final_all.fits"><img height="40" src="'+targname+'/'+targname+'_coadd_final_all.png"></a>'

            if (os.path.exists(hdr0['targname']+'_coadd_G140L_final_all.png')): 
                fuv_l_quicklook_urlstring = '<a href="'+targname+'/'+targname+'_coadd_G140L_final_all.fits"><img height="40" src="'+targname+'/'+targname+'_coadd_G140L_final_all.png"></a>'
                this_coadd = Table.read(targname+'_coadd_G140L_final_all.fits') 
                i_good = np.where(this_coadd['FLUX'] > 0) 
                median_sn = np.median(this_coadd['SN'][i_good]) 
                print 'Median SN for ', targname, ' = ', median_sn 

            if (os.path.exists(hdr0['targname']+'_coadd_G130M_final_all.fits')): 
                download_string = '<a href="'+targname+'/'+targname+'_coadd_G130M_final_all.fits">ALL</a> |'

                this_coadd = Table.read(targname+'_coadd_G130M_final_all.fits') 
                i_good = np.where(this_coadd['FLUX'] > 0) 
                median_sn = np.median(this_coadd['SN'][i_good]) 
                print 'Median SN for ', targname, ' = ', median_sn 
                
                if (os.path.exists(hdr0['targname']+'_coadd_G130M_final_lp1.fits')): 
                    download_string = download_string  + \
                                '  '+'<a href="'+targname+'/'+targname+'_coadd_G130M_final_lp1.fits">LP1</a> | '
                else: 
                    download_string = download_string+'  ' + \
                                '. . . .  | ' 
                if (os.path.exists(hdr0['targname']+'_coadd_G130M_final_lp2.fits')): 
                    download_string = download_string + \
                                '  '+'<a href="'+targname+'/'+targname+'_coadd_G130M_final_lp2.fits">LP2</a> | '
                else: 
                    download_string = download_string + \
                                '  '+'. . . .  | ' 
                if (os.path.exists(hdr0['targname']+'_coadd_G130M_final_lp3.fits')): 
                    download_string = download_string + \
                                '  '+'<a href="'+targname+'/'+targname+'_coadd_G130M_final_lp3.fits">LP3</a>   '
                else: 
                    download_string = download_string+'  '+'. . . .  ' 
            else: 
                download_string = '. . . | . . . | . . . | . . . ' 

            print """
            """,targname,"""DOWNLOAD_STRING
            """, download_string,"""
            """

            sample_webtable.add_row([number,targname_urlstring,ra, dec,n_exp_string,str.split(targdesc,';')[0],targdesc,simbad_string, mast_string, median_sn, fuv_m_quicklook_urlstring, fuv_l_quicklook_urlstring, download_string]) 

            sample_fitstable.add_row([number,targname, ra, dec, np.size(filelist),str.split(targdesc,';')[0],targdesc,median_sn]) 

            targets.add_row([1,targname,str.split(targdesc,';')[0], targdesc]) 

            dataset_list = glob.glob(os.path.join('.', '*x1d.fits'))
            print "Making Exposure Catalog: " , filelist
            
            # exposure_cat contains database of all exposures for this target 
            exposure_cat = Table(\
                names=('Flag', 'Rootname','Target Name', 'RA','DEC','PropID',\
                       'PI Name','Detector','Segment','LP','Grating', 'Cenwave','FPPOS',\
                       'Exptime','Nevents','Extended','Date','Target Description'),   
                dtype=('I3', 'S20','S35','f4','f4','I5',\
                       'S20','S4','S5','S2','S10','S10','I2',\
                       'f10','f8','S4','S12','S200'))

            for filename in filelist:
                hdulist = fits.open(filename) 
                hdr0 = hdulist[0].header
                hdr1 = hdulist[1].header
                print "Obtaining headers for :", filename 
                if (np.shape(hdulist[1].data)[0] < 1):
                    print "no data:",filename
                else: 
                     exposure_cat.add_row([1, hdr0['ROOTNAME'], hdr0['TARGNAME'], hdr0['RA_TARG'], hdr0['DEC_TARG'], \
                        hdr0['PROPOSID'], hdr0['PR_INV_L'], hdr0['DETECTOR'], hdr0['SEGMENT'], hdr0['LIFE_ADJ'],  \
                        hdr0['OPT_ELEM'], hdr0['CENWAVE'], hdr0['FPPOS'], hdr1['EXPTIME'], hdr1['NEVENTS'], \
                        hdr0['EXTENDED'], hdr1['DATE-OBS'], hdr0['TARDESCR']] )  

                if (False): 
                    for key in ["HISTORY", "COMMENT",""]: 
                        del hdr0[key] 
                    names = sorted(hdr0.keys()) 
                    rows = [hdr0[k] for k in names]
                    header_table0 = vstack([header_table0, Table(rows=[rows], names=names)])
    
                    for key in ["HISTORY", "COMMENT",""]: 
                        del hdr1[key] 
                    names = sorted(hdr1.keys()) 
                    rows = [hdr1[k] for k in names]
                    header_table1 = vstack([header_table1, Table(rows=[rows], names=names)])
                    total_number_of_headers = total_number_of_headers + 1 
                    print 'TOTAL NUMBER OF HEADERS : ', total_number_of_headers 

            ascii.write(exposure_cat, 'all_exposures.txt')  # write out the exposures for this target by itself 
            exposure_cat.write('all_exposures.html', format='jsviewer') # write out the exposures for this target by itself 
            print "ALL EXPOSURES"   
            print exposure_cat 

            number = number + 1 

            os.chdir('..')          # go back to "datapile" 

    print exposures 

    print 'TOTAL NUMBER OF HEADERS : ', total_number_of_headers 

    sample_fitstable.write(canonical_filename+'_sample.fits', format='fits', overwrite=True) 
    sample_webtable.write(canonical_filename+'_websample.fits', format='fits', overwrite=True) 
    sample_webtable.write('sample_webtable.temp', format='jsviewer') 
    sample_webtable.write(canonical_filename+'_sample_webtable.txt' ,format='ascii') 
    os.system('sed "s/&lt;/</g" sample_webtable.temp | sed "s/&gt;/>/g" > '+canonical_filename+'_sample.html') 
    os.system('rm sample_webtable.temp')

    targets.write(canonical_filename+'.info',format='ascii.fixed_width', delimiter=',') 

    exposures.write(canonical_filename+'_exposures.fits', format='fits', overwrite=True) 
    exposures.write(canonical_filename+'_exposures.html',format='jsviewer') 

    header_table0.write(canonical_filename+'_headers0.fits', format='fits', overwrite=True) 
    header_table0.write(canonical_filename+'_headers0.html', format='jsviewer')
    
    header_table1.write(canonical_filename+'_headers1.fits', format='fits', overwrite=True) 
    header_table1.write(canonical_filename+'_headers1.html', format='jsviewer')




#-----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()
    targets = args.filename
    
    scrape_headers(targets)
    sys.exit("""
    
    ~~~~~~~*~*~*~*~
    ~~~~~~~*~*~*~*~  all done!!!! spectra are fun!
    ~~~~~~~*~*~*~*~
    """)
