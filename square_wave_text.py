## updating square wave test

import numpy as np
import pdb
import os
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from astropy.io import fits
import pandas as pd
import time
import datetime
from array import *
import matplotlib.gridspec as gridspec
from astropy.io import fits
import argparse
import scipy
import timeconstantfunc as tcf


## if peaks are not detected correctly, redefine tolerance parameter

## Global constants for plotting labels/colors
clrlst = ['b','r','g','k','deeppink','cyan','magenta']
xlbl = 'Frame(n)'
ylbl = '$\Delta$ ADU'

parser = argparse.ArgumentParser(description='Frequency Response Reduction & Ramp Extraction Plots')
parser.add_argument('-fits','--fitsfile',type=str,metavar='',required=True,help='Fits File you would like to reduce and plot.')
parser.add_argument('-tol','--peak_tolerance',type=int,metavar='',required=False,help='The amount of tolerance between frame n and frame n+npeaksteps to find peak',default=12)
parser.add_argument('-nsteps','--nsteps',metavar='',required=False,type=int,help='How many steps to loo ahead when trying to find a peak',default=50)
parser.add_argument('-save','--savefigures',metavar='',required=False,type=bool,help='True if you would like to save figures, False is you do not. Defaults to True.',default=True)
parser.add_argument('-showfigs','--showfigures',metavar='',required=False,type=bool,help='True if you want to see figures as script compiled. Else other wise. Defaults to False to save time.',default=False)
parser.add_argument('-write','--writeout',metavar='',required=False,type=bool,help='True if you want full statistics textfile on code and data, False if you dont. Defaults to True.',default=True)
args = parser.parse_args()


## Handling optional arguments:
t0 = time.time()
print("\n Parsing arguments...")
dfile = args.fitsfile
tolerance = args.peak_tolerance
npeaksteps = args.nsteps
showfigs = args.showfigures
savefigs = args.savefigures

if savefigs:
    savefilename = dfile[0:7]
writefile = args.writeout
if writefile:
    start_time = time.time()
    time_sum = 0.0
    txtfile = 'info_{}.txt'.format(dfile[0:6])
    infofile = open(txtfile,"w+")
    now = datetime.datetime.now()
    header = "\n ########################\n ~~~ Square Wave Test ~~~\n ########################\n\n"
    header += "Date: \t {}-{}-{}\nTime:\t{}:{}:{}\n".format(now.year, now.month, now.day, now.hour, now.minute, now.second)
    header += "\nData File: '{}'\n".format(dfile)
    header += "\nParameters:\n\tTolerance:\t{}\t\tPeakSteps:\t{}\n".format(tolerance,npeaksteps)
    header += " \tShowFigures:\t{}\t\tSaveFigures:\t{}".format(showfigs,savefigs)
    infofile.write(header)
print("\t...arguments parsed.")

def write_file(sometext):
    infofile.write(sometext)
    return

def read_reduce(dfile):
    print(' Reducing data ...')
    if writefile:
        text = "\n\nReduction Info:"
        t0 = time.time()
    ## Opening the fits file and instantiating data to fdata variable
    ffile = fits.open(dfile)
    fdata = ffile[0].data

    ## Extracting data shape to verify it is correct. Must be 256-by-64 pixels per frame.
    d0,d1,d2 = fdata.shape[0],fdata.shape[1],fdata.shape[2]
    if d1 != 256 or d2 != 64:
        print('Warning: Data misshaped for subframe extraction. Please input data with shape n-by-256-by-64')
        pdb.set_trace() ## To prevent crashes if not correct shape
    
    temp_subframe = fdata[:,63:189,32:] ## These are NOT magic numbers, these are hard coded for specific subframe extraxtion
    reduced_subframe = []

    ## Running CDS(1,1) reduction:
    for ii in range(0,d0-1,1):
        temp = fdata[ii+1,:,:]-fdata[ii,:,:]
        temp_med = np.mean(temp)
        reduced_subframe.append(temp_med)
    
    ## Creating Pandas DataFrame to return
    dataframe = pd.DataFrame(reduced_subframe,columns=['adu'])

    if writefile:
        t1 = time.time()
        text += "\n\tReduction Time: {0:.5f} seconds".format(t1-t0)
        write_file(text)

    print('\t... data reduced.')
    return dataframe

def compress_df(dframe):
    ## Compressing DataFrame to a list of values and returning the list.
    print(' Compressing DataFrame ...')
    ds = dframe.values[:,0]
    ds = ds.tolist()
    print('\t... DataFrame compressed.')
    return ds
    
def detect_peaks(data,npeakstep=npeaksteps):
    print(' Detecting peaks ...')
    if writefile:
        text = "\n\nPeak Detection Info:"
        t0 = time.time()

    ## Creating an empty list to append
    onidx = []

    ## Excluding the ends of the data as they are often skewed.
    for ii in range(16,len(data)-npeakstep,1):
        ## Delta ADU as dadu
        dadu = np.abs(data[ii])-np.abs(data[ii+npeakstep])

        ## Checking if the change in ADU values is above the specified tolerance to be considered a peak (parameter, tol).
        ## The extra conditional, dadu > 0., is just a scanity check.

        ## This also identifies random peaks within the data, so there must be a process to get rid of these random peaks; I do this by checking the surrounding ADU values and seeing if it is an outliar
        if np.abs(dadu) >= tolerance and dadu > 0.:
            median_val = np.median(data[ii-5:ii+5])

            if median_val * 1.5 <= dadu:
                dadu = 0 ## This is a space filler, it does nothing.
            else:
                ## Meaning it passed all the checks and the index can be appended to a list.
                onidx.append(ii)

    print('\t... peaks detected.')
    if writefile:
        t1 = time.time()
        text += "\n\tPeak Detection Time: {0:.5f} seconds".format(t1-t0)
        write_file(text)
    return onidx

def check_peakidx(data,peaks):
    print(' Extracting outliar peaks...')
    if writefile:
        text = "\n\nCheck Peak Info:"
        t0 = time.time()
    ## Creating an empty list to append
    peak_idx = []

    ## This function looks 350 frames ahead to see if there are any peaks in front of it.
    ## If there are, it takes the last one as the adjusted peak. This peak is almost always the right index for LED-onset/LED-offset
    idx_bool = np.diff(peaks) >= 350
    idx = np.where(idx_bool==True)[0].tolist()
    
    for ii in range(len(idx)):
        ## Just appending the list here
        peak_idx.append(peaks[idx[ii]])

    print('\t... outlier peaks extracted.')

    ## Sending appended peaks to refine_peaks(data,peak_idx).
    ## The last index found above might be plus-or-minus 15 indicies of the exact LED_onset / LED-offset.
    ##  refine_peaks() finds exactly where LED-onset / LED-offset happens.
    if writefile:
        t1 = time.time()
        text += "\n\tCheck Peak Time: {0:.5f} seconds".format(t1-t0)
        write_file(text)

    adjusted_peaks = refine_peaks(data,peak_idx)

    return adjusted_peaks

def refine_peaks(data,peaks,ledoffset=True):
    print(" Refining peaks ...")
    if writefile:
        t0 = time.time()
        text = "\n\nRefining Peak Info:"
    
    ## Creating another list to append to
    adjusted_peaks = []

    ## Looping to find exactly where onset / offset happens.
    for ii in range(len(peaks)):
        idx = peaks[ii]
        nadjust = 2
        dadu_down = False
        while dadu_down == False:
            ## Looking nadjust indicies behind index
            temp_data = data[idx-nadjust:idx]

            ##  If executed below, there are greater Delta ADU values behind it, meaning the peak is a previous index. Change the step size and try again.
            if temp_data[0] >= temp_data[1]:
                dadu_down = False
                nadjust += 1
            
            ## If executing below, this means the previous index has a LOWER Delta ADU values, which means the Delta ADU is not constantly decreasing STARTING at this index, meaing there is a peak at this index
            else:
                dadu_down = True
                print("\tPeak {} adjusted from {} to {}".format(ii,idx,idx-nadjust))
                if writefile:
                    text += "\n\tPeak {} adjusted from index {} to {}".format(ii,idx,idx-nadjust)
                ## Appending the list
                adjusted_peaks.append(idx-nadjust)
            ## Scanity check to prevent infinite loop.
            if nadjust == 100:
                print("Stuck in loop - restart program with more defined parameters for tolerance and npeaksteps")
    print("\t...all peaks refined.")

    if writefile:
        t1 = time.time()
        text += "\n\n\tRefined Peaks Time: {0:.5f} seconds".format(t1-t0)
        write_file(text)
    return adjusted_peaks

def get_ramp_idx(onidx):
    print(' Extracting ramp indicies...')
    if writefile:
        t0 = time.time()
        text = "\n\nExtracting Ramp Index Info:"
    ## Creating list to append to
    rampidx = []
    rampoffset = []

    ## Checking the difference between onset indicies and creating a list
    dramps = np.diff(onidx).tolist()

    ## Then loop through calculating an offset. I define it as half the number of indicies from onset index n to n+1 (defined by dramps)
    for ii in range(len(dramps)):
        ## Appends an offset list
        rampoffset.append(int(dramps[ii]/2))

        ## Appends a starting ramp index
        rampidx.append(int(onidx[ii]+rampoffset[ii]))
    print('\t... ramp indicies extracted.')
    if writefile:
        t1 = time.time()
        text += "\n\tExtracted Ramp Index Time: {0:.5f} seconds".format(t1-t0)
        write_file(text)
    return rampidx,rampoffset

def plot_data(data,offidx,ridx):
    if writefile:
        t0 = time.time()
        text = "\n\nPlotting Info:"
    if showfigs or savefigs:
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
    if savefigs:
        text += "\n\tFigure Names:"
        print("\t... figure saved as '{}{}_full.pdf'".format(dfile[0:7],dfile[-7:-5]))
        if writefile:
            text += "\n\t\t - {}{}_full.pdf".format(dfile[0:7],dfile[-7:-5])
        plt.savefig("{}{}_full.pdf".format(dfile[0:7],dfile[-7:-5]))
    if showfigs == True:
        plt.show()

    ## for ymax, median first half, for ymin median second hald 
    ## plotting all ramps on one subplot
    ## if len(ridx) == 4, code only extracts three so implement conditional for exception handling
    plt.clf()
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
    ytop = np.median(ramp[:int(len(ramp)/2)])+5
    ybottom = np.median(ramp[int(len(ramp)/2):])-5
    plt.ylim(top=ytop,bottom=ybottom)
    plt.legend(prop={'size':6})
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    if savefigs:
        print("\t\t... figure saved as '{}{}_ramps.pdf'".format(dfile[0:7],dfile[-7:-5]))
        if writefile:
            text += "\n\t\t - {}{}_ramps.pdf".format(dfile[0:7],dfile[-7:-5])
        plt.savefig("{}{}_ramps.pdf".format(dfile[0:7],dfile[-7:-5]))
    if showfigs:
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
    ytop,ybottom = np.median(addramp[:int(len(addramp)/2)])+10,np.median(addramp[int(len(addramp)/2):])-5
    plt.plot(addramp,c='k',label="R0-R{} Co-Added".format(nsubplots))
    plt.ylim(top=ytop,bottom=ybottom)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.legend(prop={'size':6})
    if savefigs:
        print("\t\t... figure saved as '{}{}_coadded.pdf'".format(dfile[0:7],dfile[-7:-5]))
        if writefile:
            text += "\n\t\t - {}{}_coadded.pdf".format(dfile[0:7],dfile[-7:-5])
        plt.savefig("{}{}_coadded.pdf".format(dfile[0:7],dfile[-7:-5]))
    if showfigs:
        plt.show()

    for ii in range(len(addramp)):
        addramp[ii] = addramp[ii]/(float(nsubplots)+1) ## to account for the added ramp
    ytop = np.median(addramp[:int(len(addramp)/2)])+5
    ybottom = np.median(addramp[int(len(addramp)/2):])-5
    plt.plot(addramp,c='k',label="Co-Added Mean")
    plt.legend(prop={'size':6})
    plt.title(dfile)
    plt.ylim(top=ytop,bottom=ybottom)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    if savefigs:
        print("\t\t... figure saved as '{}{}_coadded_mean.pdf'".format(dfile[0:7],dfile[-7:-5]))
        if writefile:
            text += "\n\t\t - {}{}_coadded_mean.pdf".format(dfile[0:7],dfile[-7:-5])
        plt.savefig("{}{}_coadded_mean.pdf".format(dfile[0:7],dfile[-7:-5]))
    if showfigs:
        plt.show()

    if writefile:
        t1 = time.time()
        text += "\n\tPlotting Time: {0:.5f} seconds".format(t1-t0)
        write_file(text)


if __name__=='__main__':
    dataframe = read_reduce(dfile)
    reduced_data = compress_df(dataframe)
    allpeakidx = detect_peaks(reduced_data)
    offsetidx = check_peakidx(reduced_data,allpeakidx)
    ramp_idx,ramp_offset = get_ramp_idx(offsetidx)
    
    if writefile:
        timeconstant_onoff,onoff_text = tcf.calc_timeconstant(reduced_data,ramp_idx,offsetidx,ramp_offset,showfigs,savefigs,writefile,flipdata=True)
        timeconstant_offon,offon_text = tcf.calc_timeconstant(reduced_data,ramp_idx,offsetidx,ramp_offset,showfigs,savefigs,writefile)
        write_file(onoff_text)
        write_file(offon_text)

    else:
        timeconstant_onoff = tcf.calc_timeconstant(reduced_data,ramp_idx,offsetidx,ramp_offset,showfigs,savefigs,writefile,flipdata=True)
        timeconstant_offon = tcf.calc_timeconstant(reduced_data,ramp_idx,offsetidx,ramp_offset,showfigs,savefigs,writefile)

    plot_data(reduced_data,offsetidx,ramp_idx)

    if writefile:
        end_time = time.time()
        text = "\n\nRun Time: {0:.4f} seconds".format(end_time-start_time)
        text += "\n\n ########################\n\n ~~~ End Of Document ~~~\n\n ########################"
        write_file(text)
        infofile.close()

## Need to incorperate alpha and beta constants to writefile and everything is perfect!