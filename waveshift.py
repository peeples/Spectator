from astropy.io import fits as f
from astropy.table import Table, Column
from astropy.io import ascii
from datetime import datetime 
import astropy.units as u
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import copy 

def main(): 

this wil produce diagnostic plots of wavelength shifts for some kinds of targets 

