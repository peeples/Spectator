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
 

def h(): 

#	h = fits.open('lbiy83i0q_x1d.fits')
#
#	h0 = h[0].header 
#
#	unique_keys = sorted(set(h0))
#        r = [[h0[k] for k in unique_keys]] 
#	#r = np.reshape(r,132)
#	u = unique_keys 
#	# u = np.reshape(u,  132)  
#	
#	print 'U = ', np.size(u), np.shape(u)
#	print 'R = ', np.size(r), np.shape(r)
#	
#	t = Table(rows=[r], names=u)

	hdulist=fits.open('lbiy83i0q_x1d.fits')
	h=hdulist[0].header
	for key in ["HISTORY", "COMMENT", ""]:
		del h[key] 
	names = sorted(h.keys())
	rows = [h[k] for k in names]
	print np.shape(rows), np.shape(names) 

	t = Table(rows=[rows], names=names)



	hdulist=fits.open('lbiy83i3q_x1d.fits') 
	h=hdulist[0].header
	for key in ["HISTORY", "COMMENT", ""]:
		del h[key] 
	names = sorted(h.keys())
	rows = [h[k] for k in names]
	print np.shape(rows), np.shape(names) 
	t.add_row(rows) 



	print t['FPPOS'] 
