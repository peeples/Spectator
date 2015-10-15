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

def create_sample(canonical_filename): 

	sample = Table.read('/grp/hst/HST_spectro/datapile/all_targets_sample.fits')
	websample = Table.read('/grp/hst/HST_spectro/datapile/all_targets_websample.fits') 
	exposures = Table.read('/grp/hst/HST_spectro/datapile/all_targets_exposures.fits')
	
        canonical = ascii.read(canonical_filename+'.canonical') 
        targetlist = canonical['col2'][np.where(canonical['col1'] == 1)]

	sub_sample = sample[0:1] # copy over first row just to have the format 
	sub_websample = websample[0:1] # copy over first row just to have the format 

	number = 1 

 	for targetname in targetlist: 
 		print "Driving target:  ", targetname 

		i_target = np.where(sample['Target Name'] == targetname)[0] 
		if (np.size(i_target) > 0): 
		    print 'This target', targetname, ' is the entry ', i_target, np.size(i_target) 

		    median_sn = -9.99 

		    ra = sample['RA'][i_target][0] 
		    dec = sample['DEC'][i_target][0] 
		    targdesc = sample['Target Description'][i_target][0]
		    targ_category = str.split(targdesc,';')[0] 
		
	    	    targname_urlstring = '<a href="'+targetname+'/'+targetname+'_quicklook.html">'+targetname+'</a>'
		    print targname_urlstring 

		    simbad_urlstring = '<a href="http://simbad.u-strasbg.fr/simbad/sim-coo?CooDefinedFrames=none&CooEpoch=2000&Coord='+str(ra)+\
				'd'+str(dec)+'d&submit=submit%20query&Radius.unit=arcsec&CooEqui=2000&CooFrame=FK5&Radius=4"> SIMBAD </a>'  
		    print simbad_urlstring 

		    mast_urlstring = '<a href="https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html?searchQuery='+targetname+'"> MAST  </a>'  
		    print mast_urlstring 

		    fuv_m_quicklook_urlstring = '...' 
		    fuv_l_quicklook_urlstring = '...' 
 		    if (os.path.exists(targetname+'/'+targetname+'_coadd_final_all.png')): 
		        	fuv_m_quicklook_urlstring = '<a href="'+targetname+'/'+targetname+'_coadd_G130M_final_all.fits"><img height="40" src="'+targetname+'/'+targetname+'_coadd_final_all.png"></a>'
				print 'fuv_m_quicklook_string', fuv_m_quicklook_urlstring 

		    if (os.path.exists(targetname+'/'+targetname+'_coadd_G140L_final_all.png')): 
			fuv_l_quicklook_urlstring = '<a href="'+targetname+'/'+targetname+'_coadd_G140L_final_all.fits"><img height="40" src="'+targetname+'/'+targetname+'_coadd_G140L_final_all.png"></a>'
			this_coadd = Table.read(targetname+'/'+targetname+'_coadd_G140L_final_all.fits') 
			i_good = np.where(this_coadd['FLUX'] > 0) 
			median_sn = np.median(this_coadd['SN'][i_good]) 
			print 'Median SN for ', targetname, ' = ', median_sn 

 		    if (os.path.exists(targetname+'/'+targetname+'_coadd_G130M_final_all.fits')): 
			fuv_download_string = '<a href="'+targetname+'/'+targetname+'_coadd_G130M_final_all.fits">ALL</a> |'

			this_coadd = Table.read(targetname+'/'+targetname+'_coadd_G130M_final_all.fits') 
			i_good = np.where(this_coadd['FLUX'] > 0) 
			median_sn = np.median(this_coadd['SN'][i_good]) 
			print 'Median SN for ', targetname, ' = ', median_sn 
				
 			if (os.path.exists(targetname+'/'+targetname+'_coadd_G130M_final_lp1.fits')): 
				fuv_download_string = fuv_download_string+'  '+'<a href="'+targetname+'/'+targetname+'_coadd_G130M_final_lp1.fits">LP1</a> | '
			else: 
				fuv_download_string = fuv_download_string+'  '+'. . . .  | ' 
 			if (os.path.exists(targetname+'/'+targetname+'_coadd_G130M_final_lp2.fits')): 
				fuv_download_string = fuv_download_string+'  '+'<a href="'+targetname+'/'+targetname+'_coadd_G130M_final_lp2.fits">LP2</a> | '
			else: 
				fuv_download_string = fuv_download_string+'  '+'. . . .  | ' 
 			if (os.path.exists(targetname+'/'+targetname+'_coadd_G130M_final_lp3.fits')): 
				fuv_download_string = fuv_download_string+'  '+'<a href="'+targetname+'/'+targetname+'_coadd_G130M_final_lp3.fits">LP3</a>   '
			else: 
				fuv_download_string = fuv_download_string+'  '+'. . . .  ' 
		    else: 
			fuv_download_string = '. . . | . . . | . . . | . . . ' 



                    sample['S/N (130)'][i_target] = median_sn 
		    sub_sample = vstack([sub_sample, sample[i_target]]) 

                    websample['FUV M'][i_target] = fuv_m_quicklook_urlstring 
                    websample['FUV L'][i_target] = fuv_l_quicklook_urlstring 
                    websample['FUV M Download (by LP)'][i_target] = fuv_download_string 
                    websample['S/N (130)'][i_target] = median_sn 
		    sub_websample = vstack([sub_websample, websample[i_target]]) 

		    number = number + 1 
 
	sub_sample = sub_sample[1:]  
	sub_websample = sub_websample[1:]  

 	sub_sample.write(canonical_filename+'_sample.fits', format='fits', overwrite=True) 
 	sub_websample.write('sample.temp',format='jsviewer') 
 	sub_websample.write(canonical_filename+'_sample.txt',format='ascii') 
 	os.system('sed "s/&lt;/</g" sample.temp | sed "s/&gt;/>/g" > '+canonical_filename+'_sample.html') 
 	os.system('rm sample.temp')
        print sub_sample.keys() 
        print sub_websample.keys() 
 
#	targets.write(canonical_filename+'.list',format='ascii.fixed_width', delimiter=',') 
#
#	exposures.write(canonical_filename+'_exposures.fits', format='fits', overwrite=True) 
#	exposures.write(canonical_filename+'_exposures.html',format='jsviewer') 
#
#	header_table0.write(canonical_filename+'_headers0.fits', format='fits', overwrite=True) 
#	header_table0.write(canonical_filename+'_headers0.html', format='jsviewer')
#	
#	header_table1.write(canonical_filename+'_headers1.fits', format='fits', overwrite=True) 
#	header_table1.write(canonical_filename+'_headers1.html', format='jsviewer')
	


