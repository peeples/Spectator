from astropy.io import fits
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


import quick_look 
import line_vs_continuum 
import scrape_headers as sh 

mpl.rcParams['agg.path.chunksize'] = 10000 

def sample_script(): 

	#sh.scrape_headers("ASTEROID") 
	#sh.scrape_headers("COMET") 
	#sh.scrape_headers("SATELLITE") 
	#sh.scrape_headers("PLANET") 

	#sh.scrape_headers("CALIBRATION") 

	#sh.scrape_headers("GALAXY") 
	#sh.scrape_headers("CLUSTER-OF-GALAXIES") 

	#sh.scrape_headers("STAR") 
	#sh.scrape_headers("EXT-STAR") 
	#sh.scrape_headers("STELLAR-CLUSTER") 

	#sh.scrape_headers("ISM") 
	#sh.scrape_headers("EXT-MEDIUM") 
	
	#sh.scrape_headers("UNIDENTIFIED") 

	sh.scrape_headers("SOLAR_SYSTEM") 

	sh.scrape_headers("GALAXY_STARBURST") 
	sh.scrape_headers("GALAXY_SPIRAL") 
	sh.scrape_headers("GALAXY_STARFORMING") 
	sh.scrape_headers("GALAXY_DWARF_COMPACT") 
	sh.scrape_headers("GALAXY_EMISSION_LINE") 
	sh.scrape_headers("GALAXY_IRREGULAR") 

	sh.scrape_headers("STAR_EARLY") 
	sh.scrape_headers("STAR_LATE") 
	sh.scrape_headers("STAR_WD") 
	sh.scrape_headers("STAR_TTAURI") 
	sh.scrape_headers("STAR_DWARF_NOVAE") 
	sh.scrape_headers("STAR_POST_AGB") 
	sh.scrape_headers("STAR_LMXB") 
	
	sh.scrape_headers("QSOALS") 

	sh.scrape_headers("SUPERNOVAE") 

	sh.scrape_headers("UNALLOCATED") 
