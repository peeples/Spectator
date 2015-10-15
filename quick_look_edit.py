from astropy.io import fits as f
from astropy.table import Table, Column
from astropy.io import ascii
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import collections
import scipy.io as scio
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors
import scipy.interpolate as sp
import glob
import os
import webbrowser

def get_quick_look(): 

    mpl.rcParams['font.family'] = 'stixgeneral'

    window, wc = get_default_windows()

    lpmoves = [56131.0, 57062.0]  ## COS FUV lifetime position moves, 2012-07-23 and 2015-02-09

    #### Find the x1d's in the directory
    DATA_DIR = '.'
    dataset_list = glob.glob(os.path.join(DATA_DIR, '*x1d.fits'))
    pathname = "file://"+os.getcwd()+"/"
    pathname = '' 

    #### Grab the first file and create html page
    hdulist = f.open(dataset_list[0])
    hdr = hdulist[0].header
    targname = hdr['targname']
    outfilename = targname + "_quicklook.html"
    exists = os.path.isfile(outfilename)
    if exists:
        command = "rm -f "+outfilename
        os.system(command)
    outfile = open(outfilename,"w")

    coord = SkyCoord(ra=hdr['RA_TARG']*u.degree, dec=hdr['DEC_TARG']*u.degree)
    info = """<html>
    <head><h1 style="font-size:350%">"""+targname+"""</h1></head>
    <body><p style="font-size=200%">"""+str(hdr['TARDESCR'])+"""<br>
    at &alpha; = """+str(hdr['RA_TARG'])+""", &delta; = """+str(hdr['DEC_TARG'])+""" ("""+coord.to_string('hmsdms')+""")</font></p>
    <hr />
    """
    outfile.write(info)

    #### Now loop through all the files to figure out the demographics
    t, r = get_demographics(dataset_list)
    LAMBDA_MIN = np.min(r['minwave']-2)
    LAMBDA_MAX = np.max(r['maxwave']+2)
    MAX_FLUX = 1.05*np.max(r['maxflux'])

    t_pid = t.group_by(['PID'])
    info="""<p style="font-size=200%">COS quick look of<br>
    Programs: """
    for k in t_pid.groups.keys:
        info = info+"""<a href="http://www.stsci.edu/cgi-bin/get-proposal-info?id="""+str(k[0])+"""&observatory=HST">"""+str(k[0])+"""</a> """
    info = info + """</p>"""
    outfile.write(info)

    width = 0.7
    colors = get_default_pid_colors()
    figname = "exptime_histogram.png"
    plot_exptime_histogram(t,figname,colors=colors,t_pid=t_pid,width=width)
    addfig = r"""<br><img src='"""+pathname+figname+r"""' style="width:600pix">"""
    outfile.write(addfig)

    figname = "exptime_fppos_histogram.png"
    plot_exptime_fppos_histogram(t,figname,t_pid=t_pid,width=width)
    addfig = r"""<img src='"""+pathname+figname+r"""' style="width:600pix">"""
    outfile.write(addfig)

    figname = "lifetime_position_histogram.png"
    plot_lifetime_position_histogram(t,figname,t_pid=t_pid,width=width)
    addfig = r"""<img src='"""+pathname+figname+r"""' style="width:600pix">"""
    outfile.write(addfig)
 


     #### add coadd legend
    info = """<p style="font-size=200%">Legend: <b><font color="black">flux in black</font></b>, <b><font color="grey">errors in grey</font></b>, both smoothed over 7 pixels (~1 resel). S/N&equiv;median(flux/error), per unsmoothed pixel, in shaded window.</p>"""
    outfile.write(info)


    #### coadd????? ######

    coadd_exists = False
    if(os.path.exists(hdr['targname']+'_coadd_G130M_final_jt.fits')):
        # if the "JT" (add.main) 
        coadd_exists = True
        output_name = hdr['targname']+'_coadd_G130M_final_jt.png'
        coadd = Table.read(hdr['targname']+'_coadd_G130M_final_jt.fits') 
        print  hdr['targname']+':  quick_look opened  ' + hdr['targname']+'_coadd_G130M_final_jt.fits' 
        labeltext = """coadd of """+str(hdr['targname'])+""" COS/FUV M"""
        plot_spectrum(output_name, coadd['wave'], coadd['flux'], LAMBDA_MIN, LAMBDA_MAX, 0, MAX_FLUX, window=window, wc=wc, labeltext=labeltext, smooth=7)
        ### plt.plot(coadd['WAVELENGTH'],coadd['FLUX'],'r.') 

    if(os.path.exists(hdr['targname']+'_FUV_M_coadd.dat')):
        # if the FUV M coadd exists, use it 
        coadd_exists = True
        output_name = hdr['targname']+"_FUV_M_coadd.png"
        coadd = scio.readsav(hdr['targname']+'_FUV_M_coadd.dat')
        print hdr['targname']+':  quick_look opened  ' + hdr['targname']+'_FUV_M_coadd.dat'
        labeltext = """coadd of """+str(hdr['targname'])+""" COS/FUV M"""    
        plot_spectrum(output_name, coadd['wave'], coadd['flux'], LAMBDA_MIN, LAMBDA_MAX, 0, MAX_FLUX, window=window, wc=wc, labeltext=labeltext, error=coadd['err'], smooth=7)
        addfig = r"""<br><img src='"""+pathname+output_name+r"""' style="width:1200pix">"""
        outfile.write(addfig)
    
    if(os.path.exists(hdr['targname']+'_FUV_L_coadd.dat')):
        # if the FUV L coadd exists, use it
        coadd_exists = True
        output_name = hdr['targname']+"_FUV_L_coadd.png"
        coadd = scio.readsav(hdr['targname']+'_FUV_L_coadd.dat')
        print  'quick_look opened ' + hdr['targname']+'_FUV_L_coadd.dat'
        d = Time(hdulist[1].header['expstart'], format='mjd')
        labeltext = """coadd of """+str(hdr['targname'])+""" COS/FUV L"""
        plot_spectrum(output_name, coadd['wave'], coadd['flux'], LAMBDA_MIN, LAMBDA_MAX, 0, MAX_FLUX, window=window, wc=wc, labeltext=labeltext, error=coadd['err'], smooth=7)
        addfig = r"""<br><img src='"""+pathname+output_name+r"""' style="width:1200pix">"""
        outfile.write(addfig)
    
    if not coadd_exists:
        print "could not find a coadd, so sad"
        
    #### add legend
    info = """<p style="font-size=200%">Individual exposures<br>Legend: <b><font color="black">flux in black</font></b>, <b><font color="grey">errors in grey</font></b>, both smoothed over 7 pixels (~1 resel). S/N&equiv;median(flux/error), per unsmoothed pixel, in shaded window.</p>"""
    outfile.write(info)


    #### Loop through files and make plots!
    time_flux = []
    for filename in dataset_list:
        ### plot_spectrum(filename)
        hdulist = f.open(filename)
        hdr = hdulist[0].header
        data = hdulist[1].data
        if (hdulist[1].data == None):
            print "no data:",filename
            continue
        if (np.shape(data['flux'])[0] == 0):
            print np.shape(data['flux']),filename
            continue
        output_name = hdr['targname']+"_"+str(hdr['proposid'])+"_"+str(hdr['cenwave'])+"_"+str(hdr['rootname'])+".png"
        print "FILENAME", output_name
        labeltext = """PID """+str(hdr['proposid'])+""" visit"""+str(hdr['linenum'].split('.')[0])+"""
        """+str(hdr['primesi'])+"""/"""+str(hdr['opt_elem']) + """/"""+str(hdr['cenwave']) + """/FPPOS"""+str(hdr['fppos'])+"""
        """+ hdulist[1].header['date-obs']+""", Exptime = """ +str(int(data['exptime'][0]))+"""s"""
        plot_spectrum(output_name, data['wavelength'], data['flux'], LAMBDA_MIN, LAMBDA_MAX, 0, MAX_FLUX, window=window, wc=wc, labeltext=labeltext, error=data['error'], wgt=data['dq_wgt'], smooth=7)

        addfig = r"""<br><img src='"""+pathname+output_name+r"""' style="width:1200pix">"""
        print "writing out file ", addfig 
        outfile.write(addfig)

    #     ## calculate flux in windows because not done in plot_spectrum YET
    #     for w in range(np.shape(window)[0]):
    #         indices = np.where((data[i]["dq_wgt"] > 0) & ((data[i]["wavelength"] > window[w][0]) & (data[i]["wavelength"] < window[w][1])))
    #         if(np.size(indices) > 0):
    #             time = np.average([hdulist[1].header['expstart'], hdulist[1].header['expend']])
    #             # flux = np.average(data[i]['flux'][indices],weights=1/data[i]['error'][indices])
    #             flux = np.median(data[i]['flux'][indices])
    #             error = np.sqrt(np.average(np.average((data[i]['flux'][indices]-flux)**2,weights=data[i]['error'][indices])))/np.sqrt(len(indices))
    #             print time,flux,error,w
    #             time_flux.append([time,flux,error,w])

    # print "flux vs. time isn't yet plotted, but here it is:"
    # tf = Table(rows=time_flux, names=('time','flux','error','window'))
    # tf_w = tf.group_by('window')
    # print tf_w


    outfile.write("</body></html>")
        
#-----------------------------------------------------------------------------------------------------


def movingaverage(values,weights,window):
    indices = np.where(weights > 0)
    mvavg = np.convolve(values[indices], np.ones(window)/window,'same')                    
    return mvavg

#-----------------------------------------------------------------------------------------------------


def plot_spectrum(output_name, wavelength, flux, wavemin, wavemax, fluxmin, fluxmax, **kwargs):
    error = kwargs.get("error",None)
    labeltext = kwargs.get("labeltext","")
    window = kwargs.get("window",get_default_windows()[0])
    wc = kwargs.get("wc",get_default_windows()[1])
    smooth = kwargs.get("smooth",1)
    wgt = kwargs.get("wgt",np.ones(np.shape(flux)))
 
    fig = plt.figure(figsize=(18,6),dpi=300)
    ax = fig.add_subplot(111)
    if np.asarray(np.shape(np.shape(flux))) > 1:
        for i in range(np.shape(flux)[0]):
            indices = np.where(wgt[i] > 0)
            f = movingaverage(flux[i],wgt[i],smooth)
            wave = wavelength[i][indices]
            ax.step(wave,f,lw=1,color='black')
    else:
        indices = np.where(wgt > 0)
        f = movingaverage(flux,wgt,smooth)
        wave = wavelength[indices]
        ax.step(wave,f,lw=1,color='black')
    if "error" in kwargs:
        if np.asarray(np.shape(np.shape(flux))) > 1:
            for i in range(np.shape(flux)[0]):
                indices = np.where(wgt[i] > 0)
                e = movingaverage(error[i],wgt[i],smooth)
                wave = wavelength[i][indices]
                ax.step(wave,e,lw=1,color='grey',alpha=0.6)
        else:
            e = movingaverage(error,wgt,smooth)
            ax.step(wave,e,lw=1,color='grey',alpha=0.6)
        for w in range(np.shape(window)[0]):
            indices = np.where((wgt > 0) & (wavelength > window[w][0]) & (wavelength < window[w][1]))
            if(np.shape(indices)[1] > 0):
                medflux = np.median(flux)
                err = np.sqrt(np.average(np.average((medflux-flux)**2,weights=error)))/np.sqrt(np.size(indices))
                sn = np.median(flux[indices]/error[indices])
                if (sn > 0):
                    plt.text(window[w][1],0.75*fluxmax,"S/N="+"{:.1f}".format(sn),fontsize=10)
            plt.axvspan(window[w][0],window[w][1],facecolor=wc[w],alpha=0.5)

    # do we want to overplot any other spectra?
    if "overwave" in kwargs and "overflux" in kwargs:
        print "overplotting something"
        overwavelength = kwargs.get("overwave",[])
        overflux = kwargs.get("overflux",[])
        overwgt = kwargs.get("overwgt",np.ones(np.shape(overflux)))
        overcolor = kwargs.get("overcolor","blue")
        if np.asarray(np.shape(np.shape(overflux))) > 1:
            for i in range(np.shape(overflux)[0]):
                indices = np.where(overwgt[i] > 0)
                f = movingaverage(overflux[i],overwgt[i],smooth)
                overwave = overwavelength[i][indices]
                print "plot here",i
                ax.step(overwave,f,lw=1,color=overcolor,zorder=10)
        else:
            print "or plot here"
            indices = np.where(overwgt > 0)
            f = movingaverage(overflux,overwgt,smooth)
            overwave = overwavelength[indices]
            ax.step(overwave,f,lw=1,color=overcolor,zorder=11,alpha=0.2)

        
    plt.xlim(wavemin, wavemax)
    plt.ylim(fluxmin, fluxmax)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel(r'flux [erg / s / cm$^2/$ \AA]',fontsize=20)
    plt.xlabel(r'wavelength [\AA]',fontsize=20)
    plt.text(0.75, 0.85, labeltext,fontsize=16, transform=ax.transAxes)
    plt.tight_layout()
    if os.path.isfile(output_name):
        command = "rm -f " + output_name
        os.system(command)
    plt.savefig(output_name)
    plt.close(fig)

    
#-----------------------------------------------------------------------------------------------------
    

def get_default_windows():
    window = np.array([[1145,1155],
                       [1245,1255],
                       [1345,1355],
                       [1445,1455],
                       [1545,1555],
                       [1645,1655],
                       [1745,1755],
                       [1845,1855]])
    wc = ['indigo','violet','blue','green','gold','orange','red','darkred']
    
    return window, wc


def get_default_pid_colors():
    return ['purple','yellow','red','cyan','orange','green','blue','magenta','tan','grey','pink','maroon']


#-----------------------------------------------------------------------------------------------------

def get_demographics(dataset_list):
    print "getting the demographics"
    LYA_MIN = 1206
    LYA_MAX = 1226
    demographics = []
    ranges = []
    for filename in dataset_list:
        hdulist = f.open(filename)
        hdr = hdulist[0].header
        if hdr['DETECTOR'] != 'FUV':
            print hdr['DETECTOR']
#            continue
        data = hdulist[1].data
        if (hdulist[1].data == None):
            print "no data:",filename
            continue
        if (np.shape(data['flux'])[0] == 0):
            print np.shape(data['flux']),filename
            continue
        demographics.append((hdr['proposid'],hdr['linenum'].split('.')[0],hdr['linenum'].split('.')[1],hdr['cenwave'],hdr['fppos'],hdr['life_adj'],hdr['pr_inv_l'],hdr['rootname'],data['exptime'][0]))
        indices = np.where((data["dq_wgt"] > 0) & (data["dq"] == 0) & ((data["wavelength"] > LYA_MAX) | (data["wavelength"] < LYA_MIN)))
        minwave = 1000 
        maxwave = 2000 
        maxflux = 1e-14 
        if(np.shape(indices)[1] > 0):
            if(hdr['cenwave'] == 1280 or hdr['cenwave'] == 1105):
                minwave = 900.
                maxwave = 2160.
                maxflux = 5.0 * np.mean(data["flux"][indices])   # changed from max of flux by JT 072015 
            if((hdr['OPT_ELEM'] == 'G130M') or (hdr['OPT_ELEM'] == 'G160M')): 
                minwave = 1100
                maxwave = 1900 
                maxflux = 3.0 * np.median(data["flux"][indices])   # changed from max of flux by JT 072015 
            ranges.append((minwave, maxwave, maxflux))
        print "Incorporating minwave, maxwave, maxflux from:    ", filename, minwave, maxwave, maxflux 
        hdulist.close()
    t = Table(rows=demographics, names=('PID','Visit','expnum','cenwave','FPPOS','lp','PI_NAME','rootname','exptime'))
    r = Table(rows=ranges, names=('minwave','maxwave','maxflux'), dtype=('f8','f8','f8'))
    print "--------------------- THESE ARE THE DEMOGRAPHICS-------------------" 
    print t
    print r
    return t, r


#-----------------------------------------------------------------------------------------------------

def plot_exptime_histogram(t, figname, **kwargs):
    print "plotting exposure time by cenwave histogram"
    colors = kwargs.get("colors",get_default_pid_colors())
    t_pid = kwargs.get("t_pid",t.group_by(['PID']))
    t_cen = kwargs.get("t_cen", t.group_by(['cenwave']))
    bins = kwargs.get("bins",t_cen.groups.keys['cenwave'])
    width = kwargs.get("width",0.7)

    fig = plt.figure(figsize=(6,6),dpi=300)
    ax = fig.add_subplot(111)
    ind = np.arange(len(bins))
    bottom = np.zeros(len(ind))
    ## initialize cenwaves_dict
    cenwaves_dict = {}
    for k in range(len(t_pid.groups.keys)):
        for c in t_cen.groups.keys:
            cenwaves_dict[str(c[0])] = 0
        v = t_pid.groups[k].group_by(['cenwave'])
        for group in v.groups:
            # print t_pid.groups.keys[k][0], group['cenwave'][0], np.sum(group['exptime'])
            cenwaves_dict[str(group['cenwave'][0])] += np.sum(group['exptime'])
        # print k, ind, cenwaves_dict.values(), colors[k], bottom
        p = plt.bar(ind, cenwaves_dict.values(), width, color=colors[k],bottom=bottom,label=str(v['PID'][0]))
        bottom += cenwaves_dict.values()    
    plt.xlim(0-width/2.,np.max(ind)+1.5*width)
    plt.xticks(ind+width/2., bins)
    plt.xlabel("cenwave", fontsize=16)
    plt.ylabel("exposure time [s]", fontsize=16)
    plt.title("distribution of exposure times by cenwave",fontsize=16)
    lg = ax.legend(loc='upper left',labelspacing=0.08,fontsize=16)
    lg.draw_frame(False)
    plt.tight_layout()
    if os.path.isfile(figname):
        command = "rm -f "+figname
        os.system(command)
    plt.savefig(figname)
    plt.close(fig)

    

#-----------------------------------------------------------------------------------------------------

def plot_exptime_fppos_histogram(t, figname, **kwargs): 
    print "plotting exposure time by cenwave+fppos histogram"
    colors = kwargs.get("colors",get_default_pid_colors())
    t_pid = kwargs.get("t_pid",t.group_by(['PID']))
    t_cen = kwargs.get("t_cen", t.group_by(['cenwave']))
    bins = kwargs.get("bins",t_cen.groups.keys['cenwave'])
    width = kwargs.get("width",0.7)

    if len(bins) > 8:
        fig = plt.figure(figsize=(8,6),dpi=300)
    else:
        fig = plt.figure(figsize=(6,6),dpi=300)
    ax = fig.add_subplot(111)
    ind = np.arange(len(bins))
    bottom = {}
    for i in range(4):
        bottom[str(i)] = np.zeros(len(bins))
    ## initialize cenwaves_dict
    cenwaves_dict = {}
    for k in range(len(t_pid.groups.keys)):
        fp = t_pid.groups[k].group_by(['FPPOS'])
        for fg in fp.groups:
            for c in t_cen.groups.keys:
                cenwaves_dict[str(c[0])] = 0
            v = fg.group_by(['cenwave'])
            for group in v.groups:
                # print t_pid.groups.keys[k][0], group['cenwave'][0], np.sum(group['exptime'])
                cenwaves_dict[str(group['cenwave'][0])] += np.sum(group['exptime'])
                p = plt.bar(ind+(fg['FPPOS'][0]-1)*width/4., np.array(cenwaves_dict.values()), width/4, color=colors[k],bottom=bottom[str(fg['FPPOS'][0]-1)],label=str(v['PID'][0]))
            bottom[str(fg['FPPOS'][0]-1)] += cenwaves_dict.values()
    xticks = []
    xlbls = []
    xticks = np.concatenate((xticks,ind+width/2.))
    xlbls = np.concatenate((xlbls,np.round(bins)))
    va = np.zeros(len(ind))-0.03
    for i in range(len(ind)):
        for j in range(4):
            # print j, i+(2*j+1)*width/8.
            xticks = np.concatenate((xticks,[i+(2*j+1)*width/8.]))
            xlbls = np.concatenate((xlbls,[j+1]))
            va = np.concatenate((va,[0]))
    xticks = np.ravel(xticks)
    xlbls = np.ravel(xlbls)
    xlbls = xlbls.astype(int)
    va = np.ravel(va)
    # print xticks, xlbls
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlbls)
    for tick, y in zip( ax.get_xticklabels( ), va ):
        tick.set_y( y )
    plt.xlim(0-width/2.,np.max(ind)+1.5*width)
    plt.xlabel("cenwave", fontsize=16)
    plt.ylabel("exposure time [s]", fontsize=16)
    plt.title("distribution of exposure times by cenwave and FP-POS",fontsize=16)
    plt.tight_layout()
    if os.path.isfile(figname):
        command = "rm -f "+figname
        os.system(command)
    plt.savefig(figname)
    plt.close(fig)

    
#-----------------------------------------------------------------------------------------------------

def plot_lifetime_position_histogram(t, figname, **kwargs): 
    print "plotting exposure time by cenwave+fppos histogram"
    colors = kwargs.get("colors",get_default_pid_colors())
    t_pid = kwargs.get("t_pid",t.group_by(['PID']))
    width = kwargs.get("width",0.7)

      
    fig = plt.figure(figsize=(4,6),dpi=300)
    ax = fig.add_subplot(111)
    t_lp = t.group_by(['lp'])
    if -1 in t_lp['lp']:
        bins = [-1,1,2,3]
    else:
        bins = np.arange(3) + 1
    ind = np.arange(len(bins))
    bottom = np.zeros(len(bins))
    lp_dict = {}
    for k in range(len(t_pid.groups.keys)):
        for c in bins:
            lp_dict[str(c)] = 0
        v = t_pid.groups[k].group_by(['lp'])
        for group in v.groups:
            print t_pid.groups.keys[k][0], group['cenwave'][0], np.sum(group['exptime'])
            lp_dict[str(group['lp'][0])] += np.sum(group['exptime'])
        lp_sort = collections.OrderedDict(sorted(lp_dict.items()))
        print k, ind, lp_sort.values(), colors[k], bottom
        p = plt.bar(ind, lp_sort.values(), width, color=colors[k],bottom=bottom,label=str(v['PID'][0]))
        bottom += lp_sort.values()
    plt.xlim(0-width/2.,np.max(ind)+1.5*width)
    plt.xticks(ind+width/2., bins)
    plt.xlabel("lifetime position", fontsize=16)
    plt.ylabel("exposure time [s]", fontsize=16)
    plottitle = """distribution of exposure times
    by lifetime position"""
    plt.title(plottitle,fontsize=16)
    plt.tight_layout()
    if os.path.isfile(figname):
        command = "rm -f "+figname
        os.system(command)
    plt.savefig(figname)
    plt.close(fig)
