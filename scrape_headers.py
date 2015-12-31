#! /usr/bin/env python

import fitsio
from astropy.table import Table, vstack
from astropy.io import ascii
import numpy as np
import argparse
import glob
import os 
import sys  

def parse_args():
    '''
    Parse command line arguments.  Returns args object.
    '''
    parser = argparse.ArgumentParser(description="scrapes headers for everything in 'filename.list' file")
    parser.add_argument('filename', metavar='filename', type=str, action='store',
                        help='filename.list is the file to be read in')
    parser.add_argument('--altnames', metavar='altnames', type=str, default='blank', action='store',
                        help='altnames is the file containing the NED names list')
    parser.add_argument('--redshifts', metavar='redshifts', type=bool, default=False, action='store',
                        help='do or do not include redshifts, default if not')

    args = parser.parse_args()
    return args

#-----------------------------------------------------------------------------------------------------

def scrape_headers(targets,altnames,redshifts): 
    canonical_filename = targets
    canonical = ascii.read(canonical_filename+'.list')

    have_alt_names = False 
    if altnames != 'blank': 
       alt_namefile = altnames
       altname_table = fitsio.read(alt_namefile, ext=1, header=True)[0] 
       print 'Successfully opened '+alt_namefile+' with the alt names in it.' 
       have_alt_names = True 

    print canonical
    dirlist = canonical['targname'][np.where(canonical['flag'] == 1)]

    sample_fitstable = Table(names=('Number','Target Name', 'RA','DEC','Nexp','Target Category', \
        'Target Description','Alt Name','median S/N'), dtype=('i4','S200','f4','f4','I8','S20','S20','S30','f4')) 

    sample_webtable = Table(names=('Number','Target Name', 'RA','DEC','Nexp','Target Category', \
        'Target Description','AltName','AltClass', 'Redshift','MAST', 'Median S/N', 'FUV M', 'FUV L', 'Download', 'L Download'), \
        dtype=('i4','S200','f4','f4','S240','S20','S20','S600','S30', 'S5','S240','f4','S350','S350','S350','S350')) 

    targets = Table(names=('Flag','Target Name', 'Target Category', 'Target Description'), dtype=('i4','S200','S25', 'S350')) 

    ## exposures = Table(names=('Rootname','Target Name', 'RA','DEC','PropID','PI Name','Detector','Segment',\
        # 'LP','Grating', 'Cenwave','FPPOS','Exptime','Nevents','Extended','Date','Target Description'),   
       # dtype=('S20','S35','f4','f4','I5','S20','S4','S5','S2','S10','S10','I2','f10','f8','S4','S12','S200'))
    exposures = []
       
    #### set up the master "header table"
    ###### ---->>>>>> why does this have to be done each time ????? <<<<<---------
    # generic = "generic_x1d.fits"
    # hdr0 = fitsio.read_header(generic, 0) 
    # hdr1 = fitsio.read_header(generic, 1)
    # for key in ["HISTORY", "COMMENT", ""]:
    #     del hdr0[key]
    # names = sorted(hdr0.keys())
    # rows = [hdr0[k] for k in names]
    # header_table0 = Table(rows=[rows], names=names)

    # for key in ["HISTORY", "COMMENT", ""]:
    #     del hdr1[key]
    # names = sorted(hdr1.keys())
    # rows = [hdr1[k] for k in names]
    # header_table1 = Table(rows=[rows], names=names)

    redshift_string = ' . . . ' 
    altname_string = ' skldjf' 

    counter = 1 
    total_number_of_headers = 0 

    for dirname in dirlist: 
        print "Driving target:  ", dirname
        if (os.path.isdir(dirname)): 
            os.chdir(dirname) 

            filelist = glob.glob(os.path.join('.', '*x1d.fits.gz'))

            nfiles = np.size(filelist)
            print "There are ", nfiles, " exposures for target ", dirname 

            if nfiles > 0:			### if there are no files, then this target was aliased or something else happened and you won't be using it. 
                webtable_row, fitstable_row, targetstable_row = get_webtable_info(filelist[0], nfiles, counter)

                if have_alt_names: 
                    i_alt = np.where(altname_table[:]['Target Name'] == dirname)
                    print 'Redshift = ', altname_table[int(i_alt[0])][12] 
                    altname_string = altname_table[int(i_alt[0])][16] 
                    if 'No match' in altname_table[int(i_alt[0])][11]: altname_string = ' . . . ' 
                    altname_class = altname_table[int(i_alt[0])][14] 
                    redshift_string = str(altname_table[int(i_alt[0])][12]) 

                    webtable_row[7] = altname_string 
                    webtable_row[8] = altname_class 
                    webtable_row[9] = redshift_string 

                sample_webtable.add_row(webtable_row) 
                sample_fitstable.add_row(fitstable_row)
                targets.add_row(targetstable_row) 

                dataset_list = glob.glob(os.path.join('.', '*x1d.fits.gz'))
                print "SCRAPE_HEADERS: Making Exposure Catalog: " , filelist
  
                exposure_cat = make_exposure_catalog(filelist)
                print len(exposures)
                print "this is exposure cat:"
                print exposure_cat
                if len(exposures) == 0:
                    exposures = exposure_cat
                else:
                    exposure_tmp = exposures
                    exposures = vstack([exposure_tmp, exposure_cat])
                    # exposures.add_row(exposure_cat)
                print "THIS IS THE EXPOSURES WHAT IT LOOKS LIKE NOW:"
                print exposures
                print "DID ANYTHING PRINT"

                counter = counter + 1 

            os.chdir('..')          # go back to "datapile" 

    if not have_alt_names:  
      if ('AltName' in sample_webtable.keys()): del sample_webtable['AltName'] 
      if ('AltClass' in sample_webtable.keys()): del sample_webtable['AltClass'] 

    if not redshifts:  
      del sample_webtable['Redshift'] 

    sample_fitstable.write(canonical_filename+'_sample.fits', format='fits', overwrite=True) 

    del sample_webtable['Target Category','L Download']				# just omit this for now 
    sample_webtable.write(canonical_filename+'_websample.fits', format='fits', overwrite=True) 
    #sample_webtable.write(canonical_filename+'_sample_webtable.txt' ,format='ascii') 
    sample_webtable.write('sample_webtable.temp', format='jsviewer') 
    os.system('sed "s/&lt;/</g" sample_webtable.temp | sed "s/&gt;/>/g" > '+canonical_filename+'_sample.html') 
    os.system('rm sample_webtable.temp')

    targets.write(canonical_filename+'.info',format='ascii.fixed_width', delimiter=',') 

    exposures.write(canonical_filename+'_exposures.fits', format='fits', overwrite=True) 
    exposures.write(canonical_filename+'_exposures.html',format='jsviewer') 
 
    print 
    print "END OF SCRAPE_HEADERS" 
    print 

    # only do this if you have created the master header table in the commented out bits above 
    # header_table0.write(canonical_filename+'_headers0.fits', format='fits', overwrite=True) 
    # header_table0.write(canonical_filename+'_headers0.html', format='jsviewer')
    
    # header_table1.write(canonical_filename+'_headers1.fits', format='fits', overwrite=True) 
    # header_table1.write(canonical_filename+'_headers1.html', format='jsviewer')



#-----------------------------------------------------------------------------------------------------

def get_webtable_info(filename, nfiles, counter):
    hdr0 = fitsio.read_header(filename, 0) 
    hdr1 = fitsio.read_header(filename, 1)

    targname = hdr0['TARGNAME'].strip() 
    targdesc = hdr0['TARDESCR'].strip() 
    ra = hdr0['RA_TARG']
    dec = hdr0['DEC_TARG']

    median_sn = 0.0     
    
    targname_urlstring = '<a href="../datapile/'+targname+'/'+targname+'_quicklook.html">'+targname+'</a>'

    altname_string = ' . . . ' 
    altname_class = ' . . . ' 
    redshift_string = ' . . . ' 

    mast_string = '<a href="https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery='+str(ra)+','+str(dec)+'"> MAST  </a>'  

    n_exp_string = '<a href="../datapile/'+targname+'/all_exposures.html">'+str(nfiles)+'</a>' 

    fuv_m_quicklook_urlstring = ' . . . ' 
    fuv_l_quicklook_urlstring = ' . . . ' 

    download_string = ' . . . ' 

    G130M_coadd_exists = False 
    G160M_coadd_exists = False 
    G140L_coadd_exists = False 
    if os.path.exists(hdr0['targname'].strip()+'_coadd_G130M_final_all.fits.gz'): G130M_coadd_exists = True 
    if os.path.exists(hdr0['targname'].strip()+'_coadd_G160M_final_all.fits.gz'): G160M_coadd_exists = True 
    if os.path.exists(hdr0['targname'].strip()+'_coadd_G140L_final_all.fits.gz'): G140L_coadd_exists = True 

    if G130M_coadd_exists or G160M_coadd_exists: 
        fuv_m_quicklook_urlstring = '<a href="../datapile/'+targname+'/'+targname+'_coadd_final_all.png"><img height="40" src="../datapile/'+targname+'/'+targname+'_coadd_final_all.png"></a>'

    if (os.path.exists(hdr0['targname'].strip()+'_coadd_G140L_final_all.png')): 
        fuv_l_quicklook_urlstring = '<a href="../datapile/'+targname+'/'+targname+'_coadd_G140L_final_all.png"><img height="40" src="../datapile/'+targname+'/'+targname+'_coadd_G140L_final_all.png"></a>'
        l_download_string = '<a href="../datapile/'+targname+'/'+targname+'_coadd_G140L_final_all.fits.gz">ALL</a> |'

        this_coadd = Table.read(targname+'_coadd_G140L_final_all.fits.gz') 
        i_good = np.where((this_coadd['FLUX'] > 0) & (this_coadd['WAVE'] > 1100) & (this_coadd['WAVE'] < 1900))  ## take regions that are not lousy S/N 
        median_sn_140 = np.median(this_coadd['SN'][i_good]) 

        if (os.path.exists(hdr0['targname'].strip()+'_coadd_G140L_final_lp1.fits.gz')):
            print 'I found a coadd for G140L LP=1'
            l_download_string = l_download_string  + \
            '  '+'<a href="../datapile/'+targname+'/'+targname+'_coadd_G140L_final_lp1.fits.gz">LP1</a> | '
        else:
            l_download_string = l_download_string+'  ' + \
            '. . . .  | '
        if (os.path.exists(hdr0['targname'].strip()+'_coadd_G140L_final_lp2.fits.gz')):
            print 'I found a coadd for G140L LP=2'
            l_download_string = l_download_string + \
            '  '+'<a href="../datapile/'+targname+'/'+targname+'_coadd_G140L_final_lp2.fits.gz">LP2</a> | '
        else:
            l_download_string = l_download_string + \
            '  '+'. . . .  | '
        if (os.path.exists(hdr0['targname'].strip()+'_coadd_G140L_final_lp3.fits.gz')):
            print 'I found a coadd for G140L LP=3'
            l_download_string = l_download_string + \
            '  '+'<a href="../datapile/'+targname+'/'+targname+'_coadd_G140L_final_lp3.fits.gz">LP3</a>   '
        else:
            l_download_string = l_download_string+'  '+'. . . .  '
    else:
        l_download_string = '. . . | . . . | . . . | . . . '


#    if (os.path.exists(hdr0['targname'].strip()+'_coadd_G130M_final_all.fits')): 
#        download_string = '<a href="../datapile/'+targname+'/'+targname+'_coadd_G130M_final_all.fits">ALL</a> |'
#
#        this_coadd = Table.read(targname+'_coadd_G130M_final_all.fits') 
#        i_good = np.where(this_coadd['FLUX'] > 0) 
#        median_sn_130 = np.median(this_coadd['SN'][i_good]) 
#        print 'Median G130M SN for ', targname, ' = ', str(median_sn)[0:6] 
#    
#        if (os.path.exists(hdr0['targname'].strip()+'_coadd_G130M_final_lp1.fits')): 
#            print 'I found a coadd for G130M LP=1' 
#            download_string = download_string  + \
#            '  '+'<a href="../datapile/'+targname+'/'+targname+'_coadd_G130M_final_lp1.fits">LP1</a> | '
#        else: 
#            download_string = download_string+'  ' + \
#            '. . . .  | ' 
#        if (os.path.exists(hdr0['targname'].strip()+'_coadd_G130M_final_lp2.fits')): 
#            print 'I found a coadd for G130M LP=2' 
#            download_string = download_string + \
#            '  '+'<a href="../datapile/'+targname+'/'+targname+'_coadd_G130M_final_lp2.fits">LP2</a> | '
#        else: 
#            download_string = download_string + \
#            '  '+'. . . .  | ' 
#        if (os.path.exists(hdr0['targname'].strip()+'_coadd_G130M_final_lp3.fits')): 
#            print 'I found a coadd for G130M LP=3' 
#            download_string = download_string + \
#            '  '+'<a href="../datapile/'+targname+'/'+targname+'_coadd_G130M_final_lp3.fits">LP3</a>   '
#        else: 
#            download_string = download_string+'  '+'. . . .  ' 
#    else: 
#        download_string = '. . . | . . . | . . . | . . . ' 

### why is SN a string here ??? ####

    if G140L_coadd_exists:
        this_coadd = Table.read(targname+'_coadd_G140L_final_all.fits.gz') 
        i_good = np.where((this_coadd['FLUX'] > 0) & (this_coadd['WAVE'] > 1100) & (this_coadd['WAVE'] < 1900))  ## take regions that are not lousy S/N 
        median_sn = np.median(this_coadd['SN'][i_good]) 
        print 'Median G130M SN for ', targname, ' = ', str(median_sn)[0:6] 
        #median_sn = str(median_sn)[0:6] 

    if G160M_coadd_exists:
        this_coadd = Table.read(targname+'_coadd_G160M_final_all.fits.gz') 
        i_good = np.where(this_coadd['FLUX'] > 0) 
        median_sn = np.median(this_coadd['SN'][i_good]) 
        print 'Median G160M SN for ', targname, ' = ', str(median_sn)[0:6] 
        #median_sn = str(median_sn)[0:6] 

    if G130M_coadd_exists:
        this_coadd = Table.read(targname+'_coadd_G130M_final_all.fits.gz') 
        i_good = np.where(this_coadd['FLUX'] > 0) 
        median_sn = np.median(this_coadd['SN'][i_good]) 
        print 'Median G130M SN for ', targname, ' = ', str(median_sn)[0:6] 
        #median_sn = str(median_sn)[0:6] 



    if (os.path.exists(hdr0['targname'].strip()+'.tar.gz')): 
        download_string = '<a href="../datapile/'+targname+'/'+targname+'.tar.gz">ALL</a>'

    webtable_row = [counter, targname_urlstring, ra, dec, n_exp_string, str.split(targdesc,';')[0],
        targdesc, altname_string, altname_class, redshift_string, mast_string, median_sn, 
        fuv_m_quicklook_urlstring, fuv_l_quicklook_urlstring, download_string, l_download_string]

    print 'WEBTABLE_ROW', webtable_row 
    fitstable_row = [counter, targname, ra, dec, nfiles, str.split(targdesc,';')[0], targdesc, altname_string, median_sn]
    print fitstable_row 
    
    targetstable_row = [1,targname,str.split(targdesc,';')[0], targdesc]
    return webtable_row, fitstable_row, targetstable_row


#-----------------------------------------------------------------------------------------------------

def make_exposure_catalog(filelist):
    # exposure_cat contains database of all exposures for this target
    exposure_cat = Table(\
    names=('Flag', 'Rootname', 'Target Name', 'RA', 'DEC', 'PropID',\
                'PI Name', 'Detector', 'Segment', 'LP', 'Grating', 'Cenwave', 'FPPOS',\
                'Exptime', 'Nevents', 'Mean Flux', 'Median Flux', 'Extended', 'Date', 'Target Description'),   
            dtype=('I3', 'S20', 'S35', 'f4', 'f4', 'I5',\
                   'S20', 'S4', 'S5', 'S2', 'S10', 'I4', 'I2',\
                    'f10', 'f8', 'f8', 'f8', 'S4', 'S12', 'S200'))

    if len(filelist) == 0:
        return []

                    
    ## does alias.txt exist?
    alias_file = "alias.txt"
    if (os.path.exists(alias_file)):
        af = open(alias_file, 'r')
        targname = af.readline()
        af.close()
    else:
        print "---->>>> make_expsoure_catalog can't find ",alias_file,"!!!!! making one instead.....  <<<<----------"
        hdr0 = fitsio.read_header(filelist[0], 0) 
        targname = hdr0['TARGNAME'].strip()
        af = open(alias_file, 'w')
        af.write(targname)
        af.close()
                    
    LYA_MIN = 1206 ## should this depend on M vs L grating?
    LYA_MAX = 1226

    for filename in filelist:
        hdr0 = fitsio.read_header(filename, 0) 
        data, hdr1 = fitsio.read(filename, ext=1, header=True) 
        print "Obtaining headers for :", filename 
        if (np.shape(data)[0] < 1):
            print "no data:",filename
        else:
            hdr0['TARGNAME'] = hdr0['TARGNAME'].strip() 
            indices = np.where((data["DQ_WGT"] > 0) & (data["DQ"] == 0) & ((data["WAVELENGTH"] > LYA_MAX) | (data["WAVELENGTH"] < LYA_MIN)))
            if hdr0['DETECTOR'].strip() == "FUV":
                flag = 1
            else:
                flag = 0
            
            exposure_cat.add_row([flag, hdr0['ROOTNAME'].strip(), hdr0['TARGNAME'].strip(), hdr0['RA_TARG'], hdr0['DEC_TARG'], \
                hdr0['PROPOSID'], hdr0['PR_INV_L'].strip(), hdr0['DETECTOR'].strip(), hdr0['SEGMENT'].strip(), hdr0['LIFE_ADJ'],  \
                hdr0['OPT_ELEM'].strip(), hdr0['CENWAVE'], hdr0['FPPOS'], hdr1['EXPTIME'], hdr1['NEVENTS'], \
                np.mean(data["FLUX"][indices]), np.median(data['FLUX'][indices]), \
                hdr0['EXTENDED'].strip(), hdr1['DATE-OBS'].strip(), hdr0['TARDESCR'].strip()] )  

            ## want a way to consolidate every header keyword for every exposure into single table / file,
            ## but this is really slow. method other than vstack?
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
    return exposure_cat



#-----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()
    targets = args.filename
    
    scrape_headers(targets,args.altnames,args.redshifts)
    sys.exit("""
    
    ~~~~~~~*~*~*~*~
    ~~~~~~~*~*~*~*~  all done!!!! spectra are fun!
    ~~~~~~~*~*~*~*~
    """)
