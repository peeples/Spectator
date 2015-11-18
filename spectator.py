#! /usr/bin/env python

'''
AUTHOR: Molly Peeples
DATE: 10/29/2015
NAME: spectator.py
DESCRIPTION: gets demographics and makes quick look plots and tables for Hubble/COS data (optimized for FUV)
'''

import argparse
import sys

import drive_quick_look 
import scrape_headers

def parse_args():
    '''
    Parse command line arguments.  Returns args object.
    '''
    parser = argparse.ArgumentParser(description = \
                                     "scrapes headers and makes quicklooks for everything in 'targets.list' file")
    parser.add_argument('targets', metavar='list', type=str, action='store',
                        help="""targets.list is the file to be read in;
                              first column = flag (0,1) if target is to be used,
                              second column = target/directory name""")

    parser.add_argument('--clobber', dest='clobber', action='store_true')
    parser.add_argument('--no-clobber', dest='clobber', action='store_false', help="default is no clobbering")
    parser.set_defaults(clobber=False)

    args = parser.parse_args()
    return args



#-----------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()
    targets = (args.targets, args.clobber)    

    scrape_headers.scrape_headers(args.targets)
    drive_quick_look.drive_quick_look(targets)

    sys.exit("""
    
    ~~~~~~~*~*~*~*~~~~~~~~                                 ~*~*~*~*~~~~~~~~~~~~~
    ~~~~~~~*~*~*~*~~~~~~    all done!!!! spectra are fun!    ~~~~~~~~~~~~~~~~~~~
    ~~~~~~~*~*~*~*~~~~~~~~                                 ~*~*~*~*~~~~~~~~~~~~~
    """)
