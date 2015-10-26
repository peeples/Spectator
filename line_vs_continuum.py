from astropy.io import fits
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import glob
import os 
 

def line_vs_continuum(filename): 

	hdulist = fits.open(filename) 

	if (np.shape(hdulist[1].data)[0] < 1):
   		print filename, ":   no data:"
		return False 
        hdr = hdulist[0].header
	data = hdulist[1].data
        w = data[0]["wavelength"]
        f = data[0]["flux"]
        e = data[0]["error"]
   
        s_to_n = f / e 

	#### here's a new benign comment 

	fig = plt.figure(figsize=(6,6),dpi=300)
        ax = fig.add_subplot(111) 
        plt.plot(w, s_to_n) 
        indices = np.where(s_to_n > 1) 
        ratio = "{:4.2f}".format((1. * np.shape(indices)[1] / np.shape(w)[0]))
        plt.text(np.max(w)-50, np.median(s_to_n)*1.5, ratio, fontsize=16) 
        print filename, ":  line vs. Continuum:  median flux / error = ",  np.median(s_to_n), " and Fraction flux / error > 1:    ", ratio
        plt.savefig(hdr['rootname']+"_s_to_n.png")
	plt.close(fig)

def determine_line(): 
	
	dataset_list = glob.glob(os.path.join('.', '*x1d.fits'))

	for filename in dataset_list:
 		line_vs_continuum(filename) 

