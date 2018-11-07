## Square Wave Test 

import numpy as np
import pdb
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from astropy.io import fits
import pandas as pd
import time

dfile = '180928_LEDonoff_8000f_6.6Vbias_cube.fits'
tolerance = 4
npeaksteps = 40
n_order = 79

def get_data(dfile):
    ffile = fits.open(dfile)
    gdata = ffile[0].data
    return gdata

def polyfit_peaks(y):
    #y = y.values.flatten() # Turn into 1D Array if pd.DataFrame
    x = np.arange(0,len(y),1) # Indice Array
    z = np.polyfit(x,y,n_order) # Fit polynomial of order n_order
    f = np.poly1d(z) # Function representing the fit polynomial
    y_pred = f(x) # y-values of that polynomial
    peaks = find_peaks(y_pred) # Indicies of peaks
    x_peaks = [x[p] for p in peaks[0]] # x_peaks is list of indices of peaks of polynomial fit to freq_data
    #pdb.set_trace()
    return x_peaks,y_pred

def reduce_data(rdata):
    ## rdata for raw data sent from get_data()

    ## As long as a 1-d array or list is returned, we can change this reduction to whatever we want
    nframes,ncols,npix = rdata.shape[0],rdata.shape[1],rdata.shape[2]
    mdata = [np.median(rdata[ii,:,:]) for ii in range(nframes)]
    s1_data = [mdata[ii+1]-mdata[ii] for ii in range(nframes-1)]
    s2_data = [np.median(s1_data[ii-16:ii+16]) for ii in range(16,nframes-16)] ## DONT USE MEAN - 16ms bias in mid-peak levels
    ## median each column of pixels
    return s2_data

def detect_peaks(pdata):
    tol = tolerance
    delta_adu = np.zeros(len(pdata))
    peak_idx = []
    for ii in range(250,len(pdata)-250,1):
        delta_adu[ii] = np.abs(pdata[ii])-np.abs(pdata[ii+npeaksteps]) ## This statement will have to change if data from 0-(-40)
        if delta_adu[ii] >= tol and pdata[ii] > pdata[ii+npeaksteps]: ## tol is negative
            maxadu = np.max(pdata[ii:ii+npeaksteps])
            max_idxs = np.where(maxadu==pdata[ii:ii+npeaksteps])
            peak_idx.append(max_idxs[0][-1]+ii)  
    return peak_idx

def compress_data(udata):
    norm_data = (udata-np.min(udata))/(np.max(udata)-np.min(udata))
    std_data = (udata-np.mean(udata)/np.std(udata,ddof=1))
    return norm_data,std_data

def get_mean_ramp_vals(rampdata):
    ## Finding minimum ramp length for averaging later
    min_nramp = len(rampdata[:][0])
    for ii in range(1,len(rampdata)):
        if len(rampdata[:][ii]) <= min_nramp:
            min_nramp = len(rampdata[:][ii])
    min_nramp -= 1 ## Python indexes with zero, must accout for in looping
    mr_vals = [] 
    for ii in range(min_nramp):
        mean = 0.
        for jj in range(len(rampdata)):
            mean += ramp_data[jj][ii]
            if jj == len(rampdata)-1:
                mr_vals.append(mean/(len(rampdata)))
    return mr_vals,min_nramp

def extract_neighbor_peaks(apeakidx):
    gpeaks = []
    gpeaks.append(apeakidx[0])
    for ii in range(len(apeakidx)-1):
        if apeakidx[ii]-apeakidx[ii+1] < -50: ## Extracting neighboring peaks in all_peakidx list
            gpeaks.append(apeakidx[ii+1])
    return gpeaks

if __name__=='__main__':
    t0 = time.time()
    rawdata = get_data(dfile)
    reduced_data = reduce_data(rawdata)
    #t0 = time.time()
    #plt.plot(reduced_data)
    t1=time.time()
    print('time = {}'.format(t1-t0))
    #plt.show()
    #pdb.set_trace()
    all_peakidx = detect_peaks(reduced_data)

    goodpeaks = extract_neighbor_peaks(all_peakidx)
    ramp_offset = (goodpeaks[1]-goodpeaks[0])/4
    ramp_idx = goodpeaks-ramp_offset

    
    polypeaks,poly = polyfit_peaks(np.array(reduced_data))
    # 1: Extracting ramp data
    ramp_data = []
    poly_ramps = []
    core_std_data = []
    core_norm_data = []

    for ii in range(len(ramp_idx)-1): ## This does not get last ramp so we will fix later
        ramp_data.append(reduced_data[ramp_idx[ii]:ramp_idx[ii+1]])
        poly_ramps.append(poly[polypeaks[ii]:polypeaks[ii+1]])

    ## Returning the mean values and the lowest number of frames per ramp so the plots allign correctly
    mean_ramp_vals,n_frames = get_mean_ramp_vals(ramp_data)
    mrval_poly,nframes_poly = get_mean_ramp_vals(poly_ramps)
    all_mean_peaks = detect_peaks(mean_ramp_vals)
    mean_peak = extract_neighbor_peaks(all_mean_peaks)
    t1=time.time()
    print('{} sec'.format(t1-t0))
        
    #pdb.set_trace()
    for ii in range(len(ramp_data)):
        plt.plot(ramp_data[ii][0:n_frames], linewidth=0.75,label='Ramp {}'.format(ii+1))
        plt.plot(poly_ramps[ii][0:n_frames],linewidth=0.75,linestyle='--',label='Poly {}'.format(ii+1))

    plt.plot(mean_ramp_vals,c='k',label='Mean Ramp Value')
    #plt.plot(mrval_poly,c='r',linestyle=':',label='PolyMeanVal')
    plt.axvline(x=mean_peak[0],linestyle=':',c='k',label='Onset') ## mean_peak should only be one element, but returns as list
    #plt.axvline(x=polypeaks[0],label='PolyOnset')
    plt.subplots_adjust(left=0.05,bottom=0.08,right=.975,top=0.975)
    plt.xlabel('t (n$^{th}$ Frame)')
    plt.ylabel('$\Delta$ ADU')
    plt.legend(prop={'size':6})
    plt.savefig('all_ramps.pdf')
    plt.show()

    ## Plotting: full frequency response w/ algorithm detected peaks and full ramps 
    plt.clf()
    plt.plot(reduced_data,linewidth=0.75,c='k',label='Frequency Data')
    for ii in range(len(goodpeaks)):
        if ii == 0:
            plt.axvline(x=goodpeaks[ii],c='r',linewidth=0.75,label='Peak')
            plt.axvline(x=ramp_idx[ii],c='g',linewidth=0.75,label='Ramps')
            plt.axvline(x=polypeaks[ii],c='b',linewidth=0.75,label='PolyRamps')
        else:
            plt.axvline(x=goodpeaks[ii],c='r',linewidth=0.75)
            plt.axvline(x=ramp_idx[ii],c='g',linewidth=0.75)
            plt.axvline(x=polypeaks[ii],c='b',linewidth=0.75)

    plt.subplots_adjust(left=0.05,bottom=0.08,right=.95,top=.95)
    plt.plot(poly,c='k',linestyle='-.',label='PolyFit',linewidth=0.75)
    plt.xlabel('t (n$^{th}$ Frame)')
    plt.ylabel('$\Delta$ ADU')
    plt.legend(prop={'size': 6})
    plt.show()
    

    ## make detect peaks func() return only good peaks (include extract peaks and extract_neighboring_peaks func() in detect_peaks())