# Spectator

Quick look at spectra and demographics for Hubble/COS data.

## Installation

After cloning the repository, install via

```bash
$ python setup.py install
```

Alternatively, you may install directly from GitHub using `pip`

```bash
$ pip install git+https://github.com/peeples/Spectator.git
```

Required packages:
- python 2
- astropy: http://www.astropy.org/
- fitsio: https://pypi.python.org/pypi/fitsio/ (soon to be optional if you have astropy)
- matplotlib

## Usage

### Command line

Once installed, the package is callable directly from the command line

```bash
$ spectator targets
```

where targets is a file to be read in with first column denoting if target is 
to be used (boolean 0/1 flag) and the second column describing the 
target/directory name.

### In code

Spectator can be imported into your code as follows

```python
import spectator
```

## Known issues

Currently this only works with data taken with the Cosmic Origins Spectrograph 
(COS) aboard the Hubble Space Telescope, and many of the peculiarities, both in 
plotting and scraping the headers for information, are particular to COS. That 
said, it does work for _all_ COS data, both in the FUV and NUV, publicly 
available as of February 2016.
