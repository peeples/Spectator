# Spectator

Quick look at spectra and demographics for Hubble/COS data.

## Installation

After cloning the repository, install via

```bash
python setup.py install
```

Alternatively, you may install directly from GitHub using `pip`

```bash
pip install git+https://github.com/peeples/Spectator.git
```

Required packages:
    - python 2
    - astropy: http://www.astropy.org/
    - fitsio: https://pypi.python.org/pypi/fitsio/ (soon to be optional if you 
    have astropy)
    - matplotlib

## Usage

Once installed, the package is callable directly from the command line

```bash
$ spectator targets
```

where targets is a file to be read in with first column denoting if target is 
to be used (boolean 0/1 flag) and the second column describing the 
target/directory name.

March 7, 2016

To run:
% python spectator.py targets
where targets.list is a file to be read in with
first column = flag (0,1) if target is to be used,
second column = target/directory name"  

This will run both scrape_headers (which generates tables of information 
describing the sample) and quick_look, which generates demographics plots and 
plots of spectra

## Known issues

Currently this only works with data taken with the Cosmic Origins Spectrograph 
(COS) aboard the Hubble Space Telescope, and many of the peculiarities, both in 
plotting and scraping the headers for information, are particular to COS. That 
said, it does work for _all_ COS data, both in the FUV and NUV, publicly 
available as of February 2016.
