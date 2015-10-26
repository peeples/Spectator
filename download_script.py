
from astropy.io import fits
from astropy.table import Table, Column, vstack 
from astropy.io import ascii
import numpy as np
import glob
import os 
import sys  
import copy 
import sample_download 

def download_script(): 

    canonical_files = ['STAR_EARLY', 'STAR_LATE', 'STAR_DWARF_NOVAE', 'STAR_TTAURI', 'STAR_POST_AGB', 'STAR_WD', 'STAR_LMXB',\
	'solar_system', 'QSOALS', 'GALAXY_STARBURST','GALAXY_SPIRAL','GALAXY_STARFORMING', 'GALAXY_DWARF_COMPACT','GALAXY_EMISSION_LINE', 'GALAXY_IRREGULAR', 'SUPERNOVAE'] 

    for cfile in canonical_files: 
        sample_download.sample_download(cfile) 
