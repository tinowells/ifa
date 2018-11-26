## updating square wave test

import numpy as np
import pdb
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from astropy.io import fits
import pandas as pd
import time
from array import *
import matplotlib.gridspec as gridspec
from astropy.io import fits
import argparse

clrlst = ['b','r','g','k','deeppink','cyan','magenta']

parser = argparse.ArgumentParser(description='Frequency Response Reduction & Ramp Extraction Plots')
parser.add_argument('-fits','--fitsfile',type=str,metavar='',required=True,help='Fits File you would like to reduce and plot.')
parser.add_argument('-tol','--peak_tolerance',type=int,metavar='',required=False,help='The amount of tolerance between frame n and frame n+npeaksteps to find peak')
parser.add_argument('-nsteps','--nsteps',metavar='',required=False,type=int,help='How many steps to loo ahead when trying to find a peak')
parser.add_argument('-save','--savefigures',metavar='',required=False,type=bool,help='True if you would like to save figures, False is you dont. Defaults to False.')
parser.add_argument('-write','--writeout',metavar='',required=False,type=bool,help='True if you want full statistics textfile on code and data, False if you dont. Defaults to False.')
args = parser.parse_args()

## Handling optional arguments:
dfile = args.fitsfile
if args.peak_tolerance is None:
    tolerance = 5
else:
    tolerance = args.peak_tolerance
if args.nsteps is None:
    npeaksteps = 50
else:
    npeaksteps = args.nsteps
if args.savefigures is None:
    savefigs = False
else:
    savefigs = True
    savefilename = dfile[0:7] ## just the date and _
if args.writeout is None:
    writefile = False
else:
    writefile = True
    txtfile = '{}.txt'.format(dfile[0:7])

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
        temp_med = np.median(temp)
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
    for ii in range(0,len(data)-npeakstep,1):
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
    for ii in range(len(idx)):
        peak_idx.append(peaks[idx[ii]])
    print('\t... outlier peaks extracted.')
    return peak_idx

def get_ramp_idx(onidx):
    print(' Extracting ramp indicies...')
    rampidx = []
    dramps = np.diff(onidx)[0].tolist()
    rampoffset = int(np.median(dramps)/2)-50
    for ii in range(len(onidx)):
        rampidx.append(onidx[ii]+rampoffset)
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
    plt.savefig('{}.pdf'.format(dfile.replace('.csv','')))
    plt.show()
    
    ## for ymax, median first half, for ymin median second hald 
    ## plotting all ramps on one subplot
    for kk in range(len(ridx)-1):
        ramp = data[int(ridx[kk]):int(ridx[kk+1])]
        plt.plot(ramp,c=clrlst[kk],linewidth=0.5,label='R{}'.format(kk))
    ytop = np.median(ramp[:len(ramp)/2])+5
    ybottom = np.median(ramp[len(ramp)/2:])-5
    plt.ylim(top=ytop,bottom=ybottom)
    plt.legend(prop={'size':6})
    plt.show()

    ## plotting ramps as indicvidual subplots
    plt.clf()
    gs = gridspec.GridSpec(nrows=len(ridx)-1,ncols=1)
    for ii in range(len(ridx)-1):
        ax = plt.subplot(gs[ii,0])
        ramp = data[int(ridx[ii]):int(ridx[ii+1])]
        ytop = np.median(ramp[:len(ramp)/2])+5
        ybottom = np.median(ramp[len(ramp)/2:])-5
        plt.plot(ramp,c=clrlst[ii],label='R{}'.format(ii))
        plt.ylim(top=ytop,bottom=ybottom)
        plt.legend(prop={'size': 6})
        if ii != len(ridx)-2:
            ax.xaxis.set_visible(False)
    plt.subplots_adjust(wspace=0,hspace=0,top=0.94,bottom=0.1,left=0.1,right=0.975)
    plt.show()

if __name__=='__main__':
    ## read fits file, reduce data as CDS(1,1) and return pandas dataframe
    dataframe = read_reduce(dfile)
    reduced_data = compress_df(dataframe)
    allpeakidx = detect_peaks(reduced_data)
    offsetidx = check_peakidx(reduced_data,allpeakidx)
    ramp_idx,ramp_offset = get_ramp_idx(offsetidx)
    plot_data(reduced_data,offsetidx,ramp_idx)