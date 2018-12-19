## updating square wave test

import numpy as np
import pdb
import os
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
#from scipy.stats import signaltonoise as snr
from astropy.io import fits
import pandas as pd
import time
from array import *
import matplotlib.gridspec as gridspec
from astropy.io import fits
import argparse
import scipy
import timeconstantfunc as tcf

## if peaks are not detected correctly, redefine tolerance parameter

clrlst = ['b','r','g','k','deeppink','cyan','magenta']
xlbl = 'Frame(n)'
ylbl = '$\Delta$ ADU'

parser = argparse.ArgumentParser(description='Frequency Response Reduction & Ramp Extraction Plots')
parser.add_argument('-fits','--fitsfile',type=str,metavar='',required=True,help='Fits File you would like to reduce and plot.')
parser.add_argument('-tol','--peak_tolerance',type=int,metavar='',required=False,help='The amount of tolerance between frame n and frame n+npeaksteps to find peak',default=12)
parser.add_argument('-nsteps','--nsteps',metavar='',required=False,type=int,help='How many steps to loo ahead when trying to find a peak',default=50)
parser.add_argument('-save','--savefigures',metavar='',required=False,type=bool,help='True if you would like to save figures, False is you dont. Defaults to False.',default=False)
parser.add_argument('-write','--writeout',metavar='',required=False,type=bool,help='True if you want full statistics textfile on code and data, False if you dont. Defaults to False.',default=False)
args = parser.parse_args()

## Handling optional arguments:
print("\n Parsing arguments...")
dfile = args.fitsfile
tolerance = args.peak_tolerance
npeaksteps = args.nsteps
savefigs = args.savefigures
if savefigs:
    savefilename = dfile[0:7]
writefile = args.writeout
if writefile:
    txtfile = '{}.txt'.format(dfile[0:7])
print("\t...arguments parsed.")

def read_reduce(dfile):
    print(' Reducing data ...')
    ffile = fits.open(dfile)
    fdata = ffile[0].data
    d0,d1,d2 = fdata.shape[0],fdata.shape[1],fdata.shape[2]
    if d1 != 256 or d2 != 64:
        print('Warning: Data misshaped for subframe extraction. Please input data with shape n-by-256-by-64')
        pdb.set_trace()
    temp_subframe = fdata[:,63:189,32:] ## these are NOT magic numbers, these are hard coded for specific subframe extraxtion
    reduced_subframe = []
    ## Running CDS(1,1) reduction:
    for ii in range(0,d0-1,1):
        temp = fdata[ii+1,:,:]-fdata[ii,:,:]
        temp_med = np.mean(temp)
        reduced_subframe.append(temp_med)
    dataframe = pd.DataFrame(reduced_subframe,columns=['adu'])
    print('\t... data reduced.')
    return dataframe

def compress_df(dframe):
    print(' Compressing DataFrame ...')
    ds = dframe.values[:,0]
    ds = ds.tolist()
    print('\t... DataFrame compressed.')
    return ds
    
def detect_peaks(data,npeakstep=npeaksteps):
    print(' Detecting peaks ...')
    onidx = []
    for ii in range(16,len(data)-npeakstep,1):
        dadu = np.abs(data[ii])-np.abs(data[ii+npeakstep])
        if np.abs(dadu) >= tolerance and dadu > 0.:
            ## this also identifies random peaks, so we need to get rid of that
            ## using median to get rid of outliers
            median_val = np.median(data[ii-5:ii+5])
            if median_val * 1.5 <= dadu:
                dadu = 0
            else:
                onidx.append(ii)
    print('\t... peaks detected.')
    return onidx

def check_peakidx(data,peaks):
    print(' Extracting outliar peaks...')
    peak_idx = []
    idx_bool = np.diff(peaks) >= 350
    idx = np.where(idx_bool==True)[0].tolist()
    #plt.plot(data,c='r')
    #for ii in range(len(peaks)):
    #    plt.axvline(x=peaks[ii],linewidth=0.5,c='b')
    #plt.show()
    
    for ii in range(len(idx)):
        peak_idx.append(peaks[idx[ii]])
    print('\t... outlier peaks extracted.')
    adjusted_peaks = refine_peaks(data,peak_idx)

    return adjusted_peaks

def refine_peaks(data,peaks,ledoffset=True):
    print(" Refining peaks ...")
    adjusted_peaks = []
    for ii in range(len(peaks)):
        idx = peaks[ii]
        nadjust = 2
        dadu_down = False
        while dadu_down == False:
            temp_data = data[idx-nadjust:idx]
            #pdb.set_trace()
            if temp_data[0] >= temp_data[1]:
                dadu_down = False
                nadjust += 1
            else: ## meaning the previous index has a lower delta ADU value, meaning the delta ADU isn't constantly decreasing starting at this index ... 
                dadu_down = True
                print("\tPeak {} adjusted from {} to {}".format(ii,idx,idx-nadjust))
                adjusted_peaks.append(idx-nadjust)
            if nadjust == 100:
                print("Stuck in loop - restart program with more defined parameters for tolerance and npeaksteps")
    print("\t...all peaks refined.")
    return adjusted_peaks

def get_ramp_idx(onidx):
    print(' Extracting ramp indicies...')
    rampidx = []
    rampoffset = []
    dramps = np.diff(onidx).tolist()
    #rampoffset = int(np.median(dramps)/2)
    for ii in range(len(dramps)):
        rampoffset.append(dramps[ii]/2)
        rampidx.append(onidx[ii]+rampoffset[ii])
    print('\t... ramp indicies extracted.')
    return rampidx,rampoffset

def plot_data(data,offidx,ridx):
    print(' Plotting data...')
    plt.plot(data,color='k')
    for ii in range(len(offidx)):
        plt.axvline(x=offidx[ii],color='r',linestyle='--')
    for jj in range(len(ridx)):
        plt.axvline(x=ridx[jj],color='g',linestyle='-.')
    ## finding local max/mins for image formatting
    plt.ylim(top = np.max(data[int(len(data)*0.1):int(len(data)*0.9)]),bottom=np.min(data[int(len(data)*0.1):int(len(data)*0.9)]))
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.title(dfile)
    plt.savefig("{}{}_full.pdf".format(dfile[0:7],dfile[-7:-5]))
    plt.show()

    ## for ymax, median first half, for ymin median second hald 
    ## plotting all ramps on one subplot
    ## if len(ridx) == 4, code only extracts three so implement conditional for exception handling
    if len(ridx) > 4:
        nsubplots = 4
    else:
        nsubplots = len(ridx)-1

    for kk in range(nsubplots):
        ramp = data[int(ridx[kk]):int(ridx[kk+1])]
        plt.plot(ramp,c=clrlst[kk],linewidth=0.5,label='R{}'.format(kk))
    ## plotting last ramp bc python loops up to but not including
    lastramp = data[ridx[-1]:ridx[-1]+len(data[int(ridx[-2]):int(ridx[-1])])]
    plt.plot(lastramp,c=clrlst[kk+1],linewidth=0.5,label="R{}".format(kk+1))
    ytop = np.median(ramp[:len(ramp)/2])+5
    ybottom = np.median(ramp[len(ramp)/2:])-5
    plt.ylim(top=ytop,bottom=ybottom)
    plt.legend(prop={'size':6})
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.savefig("{}{}_ramps.pdf".format(dfile[0:7],dfile[-7:-5]))
    plt.show()

    ## plotting ramps as indicvidual subplots
    plt.clf()
    plt.ylabel(ylbl)
    gs = gridspec.GridSpec(nrows=nsubplots+1,ncols=1)
    for ii in range(nsubplots):
        ax = plt.subplot(gs[ii,0])
        ramp = data[int(ridx[ii]):int(ridx[ii+1])]
        ytop = np.median(ramp[:len(ramp)/2])+5
        ybottom = np.median(ramp[len(ramp)/2:])-5
        plt.plot(ramp,c=clrlst[ii],linewidth=0.5,label='R{}'.format(ii))
        plt.ylim(top=ytop,bottom=ybottom)
        plt.legend(prop={'size': 6})
        if ii != len(ridx)-2:
            ax.xaxis.set_visible(False)
    ## adding last subplot
    ax = plt.subplot(gs[ii+1,0])
    plt.plot(lastramp,linewidth=0.5,c=clrlst[ii+1],label="R{}".format(ii+1))
    plt.legend(prop={'size': 6})
    plt.subplots_adjust(wspace=0,hspace=0,top=0.94,bottom=0.1,left=0.1,right=0.975)
    plt.xlabel(xlbl)
    plt.savefig("{}{}_subplots.pdf".format(dfile[0:7],dfile[-7:-5]))
    #plt.ylabel(ylbl)
    plt.show()

    ## plotting co-added ramps to dampen noise
    minnramp = 2000
    for ii in range(nsubplots): ## finding shortest ramp
        ramp = data[int(ridx[ii]):int(ridx[ii+1])]
        if len(ramp) < minnramp:
            minnramp = len(ramp)
    addramp = np.zeros(minnramp)
    
    for ii in range(nsubplots):
        for jj in range(minnramp):
            addramp[jj] += data[int(ridx[ii]):int(ridx[ii+1])][jj]
    for ii in range(len(addramp)):
        addramp[ii] += lastramp[ii]
    plt.clf()
    #pdb.set_trace()
    ytop,ybottom = np.median(addramp[:len(addramp)/2])+10,np.median(addramp[len(addramp)/2:])-5
    plt.plot(addramp,c='k',label="R0-R{} Co-Added".format(nsubplots))
    plt.ylim(top=ytop,bottom=ybottom)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.legend(prop={'size':6})
    plt.savefig("{}{}_coadded.pdf".format(dfile[0:7],dfile[-7:-5]))
    plt.show()

    for ii in range(len(addramp)):
        addramp[ii] = addramp[ii]/(float(nsubplots)+1) ## to account for the added ramp
    ytop = np.median(addramp[:len(addramp)/2])+5
    ybottom = np.median(addramp[len(addramp)/2:])-5
    plt.plot(addramp,c='k',label="Co-Added Mean")
    plt.legend(prop={'size':6})
    plt.title(dfile)
    plt.ylim(top=ytop,bottom=ybottom)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.savefig("{}{}_coadded_mean.pdf".format(dfile[0:7],dfile[-7:-5]))
    plt.show()
    pdb.set_trace()


if __name__=='__main__':
    ## read fits file, reduce data as CDS(1,1) and return pandas dataframe
    dataframe = read_reduce(dfile)
    reduced_data = compress_df(dataframe)
    allpeakidx = detect_peaks(reduced_data)
    offsetidx = check_peakidx(reduced_data,allpeakidx)
    ramp_idx,ramp_offset = get_ramp_idx(offsetidx)

    timeconstant_onoff = tcf.calc_timeconstant(reduced_data,ramp_idx,offsetidx,ramp_offset,flipdata=True)
    timeconstant_offon = tcf.calc_timeconstant(reduced_data,ramp_idx,offsetidx,ramp_offset)
    #timeconstant_onoff = tcf.calc_timeconstant(reduced_data,ramp_idx,offsetidx,ramp_offset,flipdata=True)

    #get_timeconstant(reduced_data,ramp_idx,offsetidx)
    pdb.set_trace()
    plot_data(reduced_data,offsetidx,ramp_idx)

    ## compute avg_ledon for all three ramps [129-384]
    ## compute avg_ledoff for all three ramps [513-768]
    ## compute delta_onoff = avg_ledon - avg_ledoff for each ramp
    ## subtract avg_ledoff from all data in that ramp (do for all three ramps)
    ## normalize the data: divide each value by delta_onoff (for all three ramps)
    ## plot linear and log scales
    ## on log scale, calculate slope of function to get time constant
    #pdb.set_trace()

    ## pylatex for documentation
    ## SNR proportional root N 
    ## second reduction pipeline w/ columns, then MCMC QAL-EQW-ABS 3-D density code

    ## Plot LED-off-to-LED-on values of co-added ramps on natural log plot, use None-type for dADU <= int(zero)
    ## extract indices and differentiate and plot diff; diff should be constant and it is the value of the time constant

    ## update version to sq*_02.py