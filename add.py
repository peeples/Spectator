from astropy.io import fits as f
from astropy.table import Table, Column
from astropy.io import ascii
from datetime import datetime
import numpy as np
import copy
import sys
import os

def main(screen):

    #######################################################################################################
    #               STEP 1: obtain all exposures in the "all_exposures.txt" file                          #
    #######################################################################################################

    tt, h0, h1  = get_dict_of_exposures()             # obtain the headers and data for ALL exposures in the target directory

    targname = h0[h0.keys()[0]]['TARGNAME']
    print 'ADD.MAIN: We will now perform a coadd Target ', targname

    #######################################################################################################
    #               STEP 2: screen out the exposures that we don't want to have                           #
    #######################################################################################################

    lifetime = 'ALL'
    tt, h0, h1 = screen_dict_of_exposures(tt, h0, h1, screen=screen, lifetime=lifetime)
    number_of_remaining_exposures = np.size(h0.keys())
    print 'ADD.MAIN: Exposure cuts are: ', screen, ' and ', lifetime
    print 'ADD.MAIN: Number of exposures surviving the cuts: ', number_of_remaining_exposures

    #######################################################################################################
    #               STEP 3: create the output dictionary and then populate it with a things we have       #
    #######################################################################################################

    out_dict = {}             ##### start with an empty dictionary to contain the output stuff

    out_dict['exptable'] = tt  ##### these keys contain the *input* header and data file contents
    out_dict['h0table'] = h0
    out_dict['h1table'] = h1
    out_dict['targname'] = targname

    #######################################################################################################
    #               STEP 4: Add DQ_WGT entries to exposure tables AND do DQ screening                     #
    #######################################################################################################

    for nn in tt.keys():                                    #### Add a few columns to the exposure tables to include per-pixel exptime and DQ weights
        tt[nn]['EXP_PIX'] = tt[nn]['FLUX'] * 0.0 + tt[nn][0]['EXPTIME']            ####    every pixel in the output array gets the same exposure time, for now.
        tt[nn]['EXP_PIX'].unit = 's'
        tt[nn]['DQ_WEIGHT'] = tt[nn]['FLUX'] * 0.0 + 1.0                #### DQ_WGT is 1 by default
        tt[nn]['DQ_WEIGHT'].unit = '   '
        tt[nn]['FLUXFACTOR'] = np.nan_to_num(tt[nn]['FLUX'] / tt[nn]['NET'])      #### use nan_to_num to eat the nonsensical zeros in gross counts (denominator)

        number_of_segments = np.size(tt[nn]['SEGMENT'])
        print 'DQ number of segments', number_of_segments

        
        for nseg in np.arange(number_of_segments):
            test = (tt[nn]['DQ'][nseg] & 2**0) | (tt[nn]['DQ'][nseg] & 2**3)  | (tt[nn]['DQ'][nseg] & 2**6) | (tt[nn]['DQ'][nseg] & 2**7) | \
                (tt[nn]['DQ'][nseg] & 2**8) | (tt[nn]['DQ'][nseg] & 2**9) | (tt[nn]['DQ'][nseg] & 2**11) | (tt[nn]['DQ'][nseg] & 2**13)
            i_omit = np.where(test > 0)
            tt[nn]['FLUX'][nseg][i_omit] = 0.
            tt[nn]['ERROR'][nseg][i_omit] = 0.
            tt[nn]['GROSS'][nseg][i_omit] = 0.
            tt[nn]['GCOUNTS'][nseg][i_omit] = 0.
            tt[nn]['NET'][nseg][i_omit] = 0.
            tt[nn]['BACKGROUND'][nseg][i_omit] = 0.
            tt[nn]['EXP_PIX'][nseg][i_omit] = 0.
            tt[nn]['DQ_WEIGHT'][nseg][i_omit] = 0.
            #print 'We have omitted ', np.size(np.where(test > 0)), ' pixels from segment '+tt[nn]['SEGMENT'][nseg]+' of exposure ', h0[nn]['ROOTNAME'], ' for data quality reasons.'

            test = (tt[nn]['DQ'][nseg] & 2**2) | (tt[nn]['DQ'][nseg] & 2**4) | (tt[nn]['DQ'][nseg] & 2**5) | (tt[nn]['DQ'][nseg] & 2**10) | (tt[nn]['DQ'][nseg] & 2**12)
            i_deweight = np.where(test > 0)
            tt[nn]['EXP_PIX'][nseg][i_omit] = tt[nn]['EXP_PIX'][nseg][i_omit] / 2.             #### set exposure time down by 2x to deweight
            tt[nn]['DQ_WEIGHT'][nseg][i_omit] = 0.5
            #print 'We have deweighted ', np.size(i_deweight), ' pixels from segment '+tt[nn]['SEGMENT'][nseg]+' of exposure ', h0[nn]['ROOTNAME'], ' for data quality reasons.'

    #######################################################################################################
    #               STEP 5: Create the main dictionary entries for final spectra                          #
    #######################################################################################################

    count_exposures(out_dict)

    if (out_dict['number_of_G130M_exposures'] > 0):
        final_wave_130 = get_wavelength_grid('G130M')
        out_dict['G130M_final_all'] = create_coadd_format(final_wave_130)

    if (out_dict['number_of_G160M_exposures'] > 0):
        final_wave_160 = get_wavelength_grid('G160M')
        out_dict['G160M_final_all'] = create_coadd_format(final_wave_160)

    if (out_dict['number_of_G140L_exposures'] > 0):
        final_wave_140 = get_wavelength_grid('G140L')
        out_dict['G140L_final_all'] = create_coadd_format(final_wave_140)

    #######################################################################################################
    #               STEP 6: obtain the new pixel coordinates of the pixels in each x1d epxosure          #
    #######################################################################################################

    for nn in tt.keys():                                    ##### Add a few column to the exposure table for the "newpix"
        tt[nn]['NEWPIX'] = tt[nn]['FLUX'] * 0.0 +  np.arange(16384)
        tt[nn]['NEWPIX'].unit = ' '

    a = get_pix_to_wave(out_dict)

    #######################################################################################################
    #               STEPS 7: go ahead and do the final coadd                                              #
    #######################################################################################################

    print "Passing ", np.size(tt.keys()), ' exposures to cos_counts_coadd at LP = ALL '
    number_of_remaining_exposures = np.size(h0.keys())
    total_number_of_exposures = number_of_remaining_exposures
    print 'Prior to calling cos_counts_coadd, the exposure counts are: '
    print '       G130M: ', out_dict['number_of_G130M_exposures']
    print '       G160M: ', out_dict['number_of_G160M_exposures']
    print '       G140L: ', out_dict['number_of_G140L_exposures']

    if (out_dict['number_of_G130M_exposures'] > 0):
        out_dict = cos_counts_coadd(out_dict, 'G130M_final_all')
    if (out_dict['number_of_G160M_exposures'] > 0):
        out_dict = cos_counts_coadd(out_dict, 'G160M_final_all')
    if (out_dict['number_of_G140L_exposures'] > 0):
        out_dict = cos_counts_coadd(out_dict, 'G140L_final_all')

    lp1_dict = copy.deepcopy(out_dict)                             #### now we're going to do LP1 with cos_counts_coadd
    count_exposures(lp1_dict)
    print "Screening ", np.size(lp1_dict['exptable'].keys()), ' exposures to get LP = 1'
    tt_lp1, h0_lp1, h1_lp1 = screen_dict_of_exposures(lp1_dict['exptable'], lp1_dict['h0table'], lp1_dict['h1table'], screen=screen, lifetime=1)
    count_exposures(lp1_dict)
    print "There are ", number_of_remaining_exposures, ' exposures at LP = 1'
    if (lp1_dict['number_of_G130M_exposures'] > 0):
        print "We will now coadd ", lp1_dict['number_of_G130M_exposures'], " at LP = 1"
        lp1_dict['G130M_final_lp1'] = create_coadd_format(final_wave_130)
        lp1_dict = cos_counts_coadd(lp1_dict, 'G130M_final_lp1')
        if (lp1_dict['number_of_G130M_exposures'] > 0): out_dict['G130M_final_lp1'] = lp1_dict['G130M_final_all']
    if (lp1_dict['number_of_G160M_exposures'] > 0):
        print "We will now coadd ", lp1_dict['number_of_G160M_exposures'], " at LP = 1"
        lp1_dict['G160M_final_LP1'] = create_coadd_format(final_wave_160)
        lp1_dict = cos_counts_coadd(lp1_dict, 'G160M_final_LP1')
        if (lp1_dict['number_of_G160M_exposures'] > 0): out_dict['G160M_final_LP1'] = lp1_dict['G160M_final_all']
    if (lp1_dict['number_of_G140L_exposures'] > 0):
        print "We will now coadd ", lp1_dict['number_of_G140L_exposures'], " at LP = 1"
        lp1_dict['G140L_final_LP1'] = create_coadd_format(final_wave_140)
        lp1_dict = cos_counts_coadd(lp1_dict, 'G140L_final_LP1')
        if (lp1_dict['number_of_G140L_exposures'] > 0): out_dict['G140L_final_LP1'] = lp1_dict['G140L_final_all']


    lp2_dict = copy.deepcopy(out_dict)                             #### now we're going to do LP2 with cos_counts_coadd
    count_exposures(lp2_dict)
    print "Screening ", np.size(lp2_dict['exptable'].keys()), ' exposures to get LP = 2'
    tt_lp2, h0_lp2, h1_lp2 = screen_dict_of_exposures(lp2_dict['exptable'], lp2_dict['h0table'], lp2_dict['h1table'], screen=screen, lifetime=2)
    count_exposures(lp2_dict)
    number_of_remaining_exposures = np.size(h0_lp2.keys())
    print "There are ", number_of_remaining_exposures, ' exposures at LP = 2'
    if (lp2_dict['number_of_G130M_exposures'] > 0):
        print "We will now coadd ", lp2_dict['number_of_G130M_exposures'], " at LP = 2"
        lp2_dict['G130M_final_lp2'] = create_coadd_format(final_wave_130)
        lp2_dict = cos_counts_coadd(lp2_dict, 'G130M_final_lp2')
        if (lp2_dict['number_of_G130M_exposures'] > 0): out_dict['G130M_final_lp2'] = lp2_dict['G130M_final_all']
    if (lp2_dict['number_of_G160M_exposures'] > 0):
        print "We will now coadd ", lp2_dict['number_of_G160M_exposures'], " at LP = 2"
        lp2_dict['G160M_final_lp2'] = create_coadd_format(final_wave_160)
        lp2_dict = cos_counts_coadd(lp2_dict, 'G160M_final_lp2')
        if (lp2_dict['number_of_G160M_exposures'] > 0): out_dict['G160M_final_lp2'] = lp2_dict['G160M_final_all']
    if (lp2_dict['number_of_G140L_exposures'] > 0):
        print "We will now coadd ", lp2_dict['number_of_G140L_exposures'], " at LP = 2"
        lp2_dict['G140L_final_lp2'] = create_coadd_format(final_wave_140)
        lp2_dict = cos_counts_coadd(lp2_dict, 'G140L_final_lp2')
        if (lp2_dict['number_of_G140L_exposures'] > 0): out_dict['G140L_final_lp2'] = lp2_dict['G140L_final_all']

    lp3_dict = copy.deepcopy(out_dict)                             #### now we're going to do LP3 with cos_counts_coadd
    count_exposures(lp3_dict)
    print "Screening ", np.size(lp3_dict['exptable'].keys()), ' exposures to get LP = 3'
    tt_lp3, h0_lp3, h1_lp3 = screen_dict_of_exposures(lp3_dict['exptable'], lp3_dict['h0table'], lp3_dict['h1table'], screen=screen, lifetime=3)
    count_exposures(lp3_dict)
    number_of_remaining_exposures = np.size(h0_lp3.keys())
    print "There are ", number_of_remaining_exposures, ' exposures at LP = 3'
    if (lp3_dict['number_of_G130M_exposures'] > 0):
        print "We will now coadd ", lp3_dict['number_of_G130M_exposures'], " at LP = 3"
        lp3_dict['G130M_final_lp3'] = create_coadd_format(final_wave_130)
        lp3_dict = cos_counts_coadd(lp3_dict, 'G130M_final_lp3')
        if (lp3_dict['number_of_G130M_exposures'] > 0): out_dict['G130M_final_lp3'] = lp3_dict['G130M_final_all']
    if (lp3_dict['number_of_G160M_exposures'] > 0):
        print "We will now coadd ", lp3_dict['number_of_G160M_exposures'], " at LP = 3"
        lp3_dict['G160M_final_lp3'] = create_coadd_format(final_wave_160)
        lp3_dict = cos_counts_coadd(lp3_dict, 'G160M_final_lp3')
        if (lp3_dict['number_of_G160M_exposures'] > 0): out_dict['G160M_final_lp3'] = lp3_dict['G160M_final_all']
    if (lp3_dict['number_of_G140L_exposures'] > 0):
        print "We will now coadd ", lp3_dict['number_of_G140L_exposures'], " at LP = 3"
        lp3_dict['G140L_final_lp3'] = create_coadd_format(final_wave_140)
        lp3_dict = cos_counts_coadd(lp3_dict, 'G140L_final_lp3')
        if (lp3_dict['number_of_G140L_exposures'] > 0): out_dict['G140L_final_lp3'] = lp3_dict['G140L_final_all']


    #######################################################################################################
    #               STEP 8: write out the output file that the quick_look etc is going to use             #
    #######################################################################################################

    if (total_number_of_exposures > 0):
        write_output(out_dict)

    #######################################################################################################
    #               STEPS WHATEVER+1:                                                                     #
    #######################################################################################################

    return out_dict


######## END OF THE MAIN PROGRAM HERE

def count_exposures(out_dict):                                    ##### counts the number of exposures with certain properties and
                                                ##### places these in the dictionary
    counter130 = 0
    counter160 = 0
    counter140 = 0

    for i in out_dict['h0table'].keys():
        print 'COUNT_EXPOSURES: OPT_ELEM: ', out_dict['h0table'][i]['ROOTNAME'], out_dict['h0table'][i]['OPT_ELEM']
        if(out_dict['h0table'][i]['OPT_ELEM'] == 'G130M'):
            counter130 = counter130 + 1
        if(out_dict['h0table'][i]['OPT_ELEM'] == 'G160M'):
            counter160 = counter160 + 1
        if(out_dict['h0table'][i]['OPT_ELEM'] == 'G140L'):
            counter140 = counter140 + 1

    out_dict['number_of_G130M_exposures'] = counter130
    out_dict['number_of_G160M_exposures'] = counter160
    out_dict['number_of_G140L_exposures'] = counter140

    print 'COUNT_EXPOSURES: There are ', counter130, out_dict['number_of_G130M_exposures']
    print 'COUNT_EXPOSURES: There are ', counter160, out_dict['number_of_G160M_exposures']
    print 'COUNT_EXPOSURES: There are ', counter140, out_dict['number_of_G140L_exposures']


def create_coadd_format(wavegrid):                                ##### this creates a dictionary that holds the canonical outputs of the coadd
                                                ##### this has to be a dictionary because it contains a bunch of tables
                                                ##### with different levels of the coadd hierarchy in them
    flux = wavegrid * 0.0                                     ##### with the input wavelength grid
    error = wavegrid * 0.0

    print 'CREATE_OUTPUT: creating master output table'
    t = Table([wavegrid, flux, error], names=('WAVE','FLUX','ERROR'))                ##### these are the basics

    t['SIGMA_UP'] = t['WAVE'] * 0.0
    t['SIGMA_DOWN'] = t['WAVE'] * 0.0
    t['GROSSCOUNTS'] = t['WAVE'] * 0.0
    t['NETCOUNTS'] = t['WAVE'] * 0.0                                 ##### will contain the S/N per pixel
    t['BACK'] = t['WAVE'] * 0.0
    t['EXP_PIX'] = t['WAVE'] * 0.0
    t['DQ'] = t['WAVE'] * 0.0
    t['SN'] = t['WAVE'] * 0.0                                     ##### this is the factor multipled by gross count rate to obtain flux
    t['FLUXFACTOR'] = t['WAVE'] * 0.0
    t['WGTFLUX'] = t['WAVE'] * 0.0
    t['NEWPIX'] =  t['WAVE'] * 0.0
    t['CHI2'] = t['WAVE'] * 0.0                                    ##### will contain a measure of statistical fluctuation between coadd and inputs
    ##### add any other desired contents of the main output final data table here. . . .

    return t                                                ##### return the output dictionary to the calling routine


def get_pix_to_wave(out_dict):

    exptable = out_dict['exptable']                                        #### we'll just use this variable inside this routine

    number_of_exposures = np.size(exptable.keys())

    oldpix = np.empty([ 2, 16384, number_of_exposures ])
    newpix = np.empty([ 2, 16384, number_of_exposures ])

    if (out_dict.__contains__('G130M_final_all')):
        final_vector_130 = out_dict['G130M_final_all']['WAVE']                        #### this is the "Final" G130 wavelength vector to interpolate to
    if (out_dict.__contains__('G160M_final_all')):
        final_vector_160 = out_dict['G160M_final_all']['WAVE']                        #### this is the "Final" G160 wavelength vector to interpolate to
    if (out_dict.__contains__('G140L_final_all')):
        final_vector_140 = out_dict['G140L_final_all']['WAVE']                        #### this is the "Final" G160 wavelength vector to interpolate to

    exp_keys = exptable.keys()

    for nn in np.arange(number_of_exposures):
        exptable[exp_keys[nn]]['NEWPIX'] = exptable[exp_keys[nn]]['FLUX'] * 0.0 +  np.arange(16384)            #### add NEWPIX vectors to the table for each exposure, no units
        exptable[exp_keys[nn]]['NEWPIX'].unit = ' '

        print 'exptable', nn, exptable[exp_keys[nn]]

        for i in np.arange(2):                                    #### we may not actually need oldpix
            oldpix[i, :, nn] = np.arange(16384)

        number_of_segments = np.size(exptable[exp_keys[nn]]['SEGMENT'])

        if (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G130M'):
            print out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'], exptable[exp_keys[nn]]['SEGMENT']
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_130, np.arange(50000)))         #### segment A for 130
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]
        elif (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G160M'):
            print out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'], exptable[exp_keys[nn]]['SEGMENT']
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_160, np.arange(50000)))         #### segment A for 130
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]
        elif (out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'] == 'G140L'):
            print out_dict['h0table'][exp_keys[nn]]['OPT_ELEM'], exptable[exp_keys[nn]]['SEGMENT']
            for nseg in np.arange(number_of_segments):
                newpix[nseg,:,nn] = np.round(np.interp(exptable[exp_keys[nn]]['WAVELENGTH'][nseg], final_vector_140, np.arange(50000)))         #### segment A for 130
                exptable[exp_keys[nn]]['NEWPIX'][nseg] = newpix[nseg,:,nn]
        else:
            print 'OH CRAP GET_PIX_TO_WAVE DID NOT FIND THE GRATING YOU WANT'

    print "STILL NEED TO IMPLEMENT DQ FLAGGING FOR HITCHES, HERE, HERE, AND HERE"

def get_wavelength_grid(grating):

    wavefile = "xaa18189l_disp.fits"
    print "looking for wavelength grid file",wavefile,"......should really be passing in a reference/ dir........"
    if os.path.isfile("/grp/hst/HST_spectro/datapile/"+wavefile):
        a = Table.read('/grp/hst/HST_spectro/datapile/'+wavefile)
    elif os.path.isfile("../../reference/"+wavefile):
        a = Table.read('../../reference/'+wavefile)
    elif os.path.isfile("../../msp-code/reference/"+wavefile):
        a = Table.read('../../msp-code/reference/'+wavefile)
    #a = Table.read('/Users/tumlinson/Dropbox/MESA/HSTSpecLegacy/COSLegacyDatapile/code/reference/xaa18189l_disp.fits')
    

    if (grating == 'G130M' or grating == 'G160M'):
        a = a[np.where(a['OPT_ELEM'] == grating)]
        a = a[np.where(a['CENWAVE'] > 1250)]
        a = a[np.where(a['APERTURE'] == 'PSA')]

        n_settings = np.size(a['COEFF'][:,1])

        w = np.float64(np.empty([16384, n_settings]))      #### create array for default wavelength solutions

        for i in np.arange(n_settings):
            w[:,i] = (np.arange(16384) - a['D'][i] + a['D_TV03'][i]) * np.float64(a['COEFF'][i,1]) + a['COEFF'][i,0]

        mean_dispersion = np.mean(a['COEFF'][:,1])     ##### this will be the "dispersion" of the final wavelength vector
        nominal_wavelength = np.mean(a['COEFF'][:,1] * (16384/2.) + a['COEFF'][:,0])            ##### this will be the "zeropoint" of the final wavelength vector
        print 'Will use mean dispersion: ', mean_dispersion, '  and nominal wavelengths : ', nominal_wavelength
        final_vector = np.float64(np.arange(50000)-50000/2) * mean_dispersion + nominal_wavelength    ##### create the final wavelength vector
    else:
        a = a[np.where(a['OPT_ELEM'] == grating)]
        a = a[np.where(a['APERTURE'] == 'PSA')]
        a = a[np.where(a['SEGMENT'] == 'FUVA')]

        n_settings = np.size(a['COEFF'][:,1])
        print 'there are N ', n_settings

        w = np.float64(np.empty([16384, n_settings]))   #### create array for default wavelength solutions

        for i in np.arange(n_settings):
            w[:,i] = (np.arange(16384) - a['D'][i] + a['D_TV03'][i]) * np.float64(a['COEFF'][i,1]) + a['COEFF'][i,0]
            print 'tell ', a['COEFF'][i,0], a['COEFF'][i,1], a['COEFF'][i,2]

        mean_dispersion = np.mean(a['COEFF'][:,1])    ##### this will be the "dispersion" of the final wavelength vector
        nominal_wavelength = np.mean(a['COEFF'][:,1] * (16384/2.) + a['COEFF'][:,0])            ##### this will be the "zeropoint" of the final wavelength vector
        print 'Will use mean dispersion: ', mean_dispersion, '  and nominal wavelengths : ', nominal_wavelength
        final_vector = np.float64(np.arange(50000)-50000/2) * mean_dispersion + nominal_wavelength    ##### create the final wavelength vector

    return final_vector


def cos_counts_coadd(out_dict, key_to_combine):

    if(out_dict.__contains__(key_to_combine)):
        final_table = out_dict[key_to_combine]

    exptable = out_dict['exptable']       #### NOTE: exptable[i][0] refers to segment A, exptable[i][1] refers to segment B
    exp_keys = exptable.keys()         #### such that exptable[i][0]['WAVELENGTH'] is the segment A wavelength vector and so on
    h0 = out_dict['h0table']
    number_of_exposures = np.size(exptable.keys())
    print "COS_COUNTS_COADD: There are ", number_of_exposures, ' exposures to coadd for key_to_combine', key_to_combine

    temp_chi2_array = np.empty([np.size(final_table['WAVE']), number_of_exposures])*0.0            ##### an empty array that we'll place the interpolaated exposures into

    exp_counter = 0

    for i in np.arange(number_of_exposures):

        number_of_segments = np.size(exptable[exp_keys[i]]['SEGMENT'])

        if h0[exp_keys[i]]['OPT_ELEM'] in key_to_combine:

            print 'Stacking up : ', h0[exp_keys[i]]['ROOTNAME'], h0[exp_keys[i]]['OPT_ELEM'], key_to_combine

            for nseg in np.arange(number_of_segments):                #### loop over segments
                newpix = out_dict['exptable'][exp_keys[i]]['NEWPIX'][nseg,:]

                final_table['NEWPIX'][newpix.astype(int)] = newpix.astype(int)
                final_table['FLUX'][newpix.astype(int)] = final_table['FLUX'][newpix.astype(int)] \
                         + exptable[exp_keys[i]][nseg]['FLUX'] * exptable[exp_keys[i]][nseg]['EXP_PIX']
                final_table['EXP_PIX'][newpix.astype(int)] = \
                            final_table['EXP_PIX'][newpix.astype(int)] + exptable[exp_keys[i]][nseg]['EXP_PIX']
                final_table['ERROR'][newpix.astype(int)] = final_table['ERROR'][newpix.astype(int)] + \
                            (exptable[exp_keys[i]][nseg]['EXP_PIX'] * exptable[exp_keys[i]][nseg]['ERROR'])**2
                final_table['GROSSCOUNTS'][newpix.astype(int)] = final_table['GROSSCOUNTS'][newpix.astype(int)] + \
                            exptable[exp_keys[i]][nseg]['GROSS'] * exptable[exp_keys[i]][nseg]['EXP_PIX']
                final_table['NETCOUNTS'][newpix.astype(int)] = final_table['NETCOUNTS'][newpix.astype(int)] + \
                            exptable[exp_keys[i]][nseg]['NET'] * exptable[exp_keys[i]][nseg]['EXP_PIX']
                final_table['DQ'][newpix.astype(int)] = exptable[exp_keys[i]][nseg]['DQ'] #### this is not correct
                                                      #### this needs to get bitwise complicated
                final_table['FLUXFACTOR'][newpix.astype(int)] = final_table['FLUXFACTOR'][newpix.astype(int)] + \
                            exptable[exp_keys[i]][nseg]['FLUXFACTOR']*exptable[exp_keys[i]][nseg]['EXP_PIX']

                temp_chi2_array[newpix.astype(int),i] = exptable[exp_keys[i]][nseg]['NET'] # this is in COUNT RATES

            exp_counter = exp_counter + 1                        #### increment the number of exposures we've added
        else:
            print "YOU'RE IN BIG TROUBLE SINCE I DO NOT KNOW HOW TO COADD FILES WITH THAT OPT_ELEM", h0[exp_keys[i]]['ROOTNAME'], h0[exp_keys[i]]['OPT_ELEM']

    i_good_pixels = np.where(final_table['EXP_PIX'] > 0.1)

    #THERE MIGHT BE SOME PROBLEM WITH THE G140L walvenehgh solution leaving little gaps !!

    if (exp_counter > 0):
        final_table['FLUX'][i_good_pixels] = final_table['FLUX'][i_good_pixels] / final_table['EXP_PIX'][i_good_pixels]     #### exposure-weighted mean
        final_table['GROSSCOUNTS'][i_good_pixels] = final_table['GROSSCOUNTS'][i_good_pixels] / final_table['EXP_PIX'][i_good_pixels]
        final_table['NETCOUNTS'][i_good_pixels] = final_table['NETCOUNTS'][i_good_pixels] / final_table['EXP_PIX'][i_good_pixels]
        final_table['ERROR'][i_good_pixels] = np.sqrt(final_table['ERROR'][i_good_pixels]) / final_table['EXP_PIX'][i_good_pixels]
        final_table['SN'][i_good_pixels] = final_table['FLUX'][i_good_pixels] / final_table['ERROR'][i_good_pixels]        #### S/N per "pixel" or "bin"
        final_table['FLUXFACTOR'][i_good_pixels] = final_table['FLUXFACTOR'][i_good_pixels] / final_table['EXP_PIX'][i_good_pixels]
        final_table['WGTFLUX'][i_good_pixels] = final_table['NETCOUNTS'][i_good_pixels] * final_table['FLUXFACTOR'][i_good_pixels]


    for i in np.arange(number_of_exposures):
        temp_chi2_array[:,i] = temp_chi2_array[:,i] - final_table['NETCOUNTS']
        temp_chi2_array[:,i] = temp_chi2_array[:,i] * temp_chi2_array[:,i] / final_table['NETCOUNTS']          #### should divide by error on net counts instead?
        temp_chi2 = temp_chi2_array[:,0] * 0.0
        temp_chi2 =  np.sum((temp_chi2_array)**1, axis=1)
    out_dict['temp_chi2_array'] = temp_chi2_array
    out_dict['temp_chi2'] = temp_chi2

    final_table['CHI2'] = temp_chi2                                            ### this should be around 1, but isn't - why?


    if(out_dict.__contains__(key_to_combine)):    out_dict[key_to_combine] = final_table

    return out_dict


def get_exposure_list(exposure_list):

    catalog = ascii.read(exposure_list)
    if "Flag" in catalog.colnames:
        catalog.remove_rows(np.where(catalog['Flag'] == 0))    #### delete exposures flagged with 0 in the all_exposures.txt file
    return catalog



def get_dict_of_exposures():

    dataset_list =  get_exposure_list("all_exposures.txt")['Rootname']
    h0 = {}                            # this will hold the headers [fits extention 0]
    h1 = {}                            # this will hold the headers [fits extention 1]
    tt = {}                            # this will hold the data

    for i in range(np.size(dataset_list)):

        #get the header
        hdu = f.open(dataset_list[i]+'_x1d.fits')
        h0[i] = hdu[0].header
        h1[i] = hdu[1].header

        t_new = Table.read(dataset_list[i]+'_x1d.fits')
        #print 'GET_DICT_OF_EXPOSURES: We have read in ', dataset_list[i]+'_x1d.fits  for targname  ', h0[0]['TARGNAME']
        tt[i] = t_new

    print 'GET_DICT_OF_EXPOSURES: ', dataset_list
    print 'GET_DICT_OF_EXPOSURES: We have read in i = ', i+1, '  exposures '

    return tt, h0, h1


def screen_dict_of_exposures(tt, h0, h1, screen='FUVM', lifetime='ALL'):         # screen out the dictionary entries we don't want
                                            # for now this just cuts out all exposures that are NUV, L, or "blue mode"
                                            # this will need to be modified once the code needs to deal with L gratings and blue modes.
    #print 'SCREEN_DICT_OF_EXPOSURES: The requested screen is: ', screen
    #print 'SCREEN_DICT_OF_EXPOSURES: The requested LP is: ', lifetime

    if (screen == 'FUVM' and lifetime == 'ALL'):                # for FUVM and all LPs
         for i in h0.keys():
            if (h0[i]['DETECTOR'] == 'FUV' and (h0[i]['OPT_ELEM'] == 'G130M' or h0[i]['OPT_ELEM'] == 'G160M') and h0[i]['CENWAVE'] > 1250):
                print 'SCREEN_DICT_OF_EXPOSURES1: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
            else:
                print 'SCREEN_DICT_OF_EXPOSURES1: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'FUVM'):                     # for FUVM and a specific LP
         for i in h0.keys():
            if (h0[i]['DETECTOR'] == 'FUV' and (h0[i]['OPT_ELEM'] == 'G130M' or h0[i]['OPT_ELEM'] == 'G160M') and h0[i]['CENWAVE'] > 1250 and h0[i]['LIFE_ADJ'] == lifetime):
                print 'SCREEN_DICT_OF_EXPOSURES2: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
            else:
                print 'SCREEN_DICT_OF_EXPOSURES2: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'FUVL' and lifetime == 'ALL'):                     # for FUVM and a specific LP
         for i in h0.keys():
            if (h0[i]['DETECTOR'] == 'FUV' and h0[i]['OPT_ELEM'] == 'G140L'):
                print 'SCREEN_DICT_OF_EXPOSURES2: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
            else:
                print 'SCREEN_DICT_OF_EXPOSURES2: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
                del tt[i]
                del h0[i]
                del h1[i]
    elif (screen == 'FUVL'):                        # for FUVM and a specific LP
         for i in h0.keys():
            if (h0[i]['DETECTOR'] == 'FUV' and h0[i]['OPT_ELEM'] == 'G140L' and h0[i]['LIFE_ADJ'] == lifetime):
                print 'SCREEN_DICT_OF_EXPOSURES3: Keeping dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
            else:
                print 'SCREEN_DICT_OF_EXPOSURES3: Deleting dictionary entry for ', h0[i]['ROOTNAME'], h0[i]['DETECTOR'], h0[i]['OPT_ELEM'], h0[i]['CENWAVE'], h0[i]['LIFE_ADJ']
                del tt[i]
                del h0[i]
                del h1[i]
    else:
        print "SORRY BUDDY, BUT I CAN'T DO NON-FUV SCREENS"

    return tt, h0, h1




def write_output(out_dict):

    ##### note this currently does not write out the fits header, prob requires astropy.io.fits
    if (out_dict.__contains__('G130M_final_all')):
        targname = out_dict['h0table'][out_dict['h0table'].keys()[0]]['TARGNAME']
        #out_dict['G130M_final_all'].write(targname+'_coadd_G130M_final_all_old.fits', format='fits', overwrite=True)
        write_fits_output(out_dict,'G130M_final_all')

    if (out_dict.__contains__('G160M_final_all')):
        targname = out_dict['h0table'][out_dict['h0table'].keys()[0]]['TARGNAME']
        #out_dict['G160M_final_all'].write(targname+'_coadd_G160M_final_all.fits', format='fits', overwrite=True)
        write_fits_output(out_dict,'G160M_final_all')

    if (out_dict.__contains__('G140L_final_all')):
        targname = out_dict['h0table'][out_dict['h0table'].keys()[0]]['TARGNAME']
        #out_dict['G140L_final_all'].write(targname+'_coadd_G140L_final_all.fits', format='fits', overwrite=True)
        write_fits_output(out_dict,'G140L_final_all')

    if (out_dict.__contains__('G130M_final_lp1')):
        #out_dict['G130M_final_lp1'].write(targname+'_coadd_G130M_final_lp1.fits', format='fits', overwrite=True)
        write_fits_output(out_dict,'G130M_final_lp1')

    if (out_dict.__contains__('G160M_final_lp1')):
        #out_dict['G160M_final_lp1'].write(targname+'_coadd_G160M_final_lp1.fits', format='fits', overwrite=True)
        write_fits_output(out_dict,'G160M_final_lp1')

    if (out_dict.__contains__('G130M_final_lp2')):
        #out_dict['G130M_final_lp2'].write(targname+'_coadd_G130M_final_lp2.fits', format='fits', overwrite=True)
        write_fits_output(out_dict,'G130M_final_lp2')

    if (out_dict.__contains__('G160M_final_lp2')):
        #out_dict['G160M_final_lp2'].write(targname+'_coadd_G160M_final_lp2.fits', format='fits', overwrite=True)
        write_fits_output(out_dict,'G160M_final_lp2')

    if (out_dict.__contains__('G130M_final_lp3')):
        #out_dict['G130M_final_lp3'].write(targname+'_coadd_G130M_final_lp3.fits', format='fits', overwrite=True)
        write_fits_output(out_dict,'G130M_final_lp3')

    if (out_dict.__contains__('G160M_final_lp3')):
        #out_dict['G160M_final_lp3'].write(targname+'_coadd_G160M_final_lp3.fits', format='fits', overwrite=True)
        write_fits_output(out_dict,'G160M_final_lp3')



def write_fits_output(out_dict, key_to_write_out):
    # can eventually expand this to include multiple extendsion for each LP instead of keeping those in separate files

    table_to_write_out = out_dict[key_to_write_out]
    prihdr = f.Header()
    prihdr['CREATOR'] = 'JT'                     # OK, here's how we can add keywords to the primary header
    prihdr['TIMESTAMP'] = str(datetime.now())
    prihdr['LEGACY'] = 'HST-50MP'
    prihdu = f.PrimaryHDU(header=prihdr)

    table_hdu = f.BinTableHDU.from_columns(np.array(table_to_write_out.filled()))
    table_hdu.header['LP'] = 'ALL'                    # this is how to add header keywords to the table extension!
    if 'lp1' in key_to_write_out: table_hdu.header['LP'] = 'LP1'
    if 'lp2' in key_to_write_out: table_hdu.header['LP'] = 'LP2'
    if 'lp3' in key_to_write_out: table_hdu.header['LP'] = 'LP3'

    thdulist = f.HDUList([prihdu, table_hdu])

    thdulist.writeto(out_dict['targname']+'_coadd_'+key_to_write_out+'.fits',clobber=True)
