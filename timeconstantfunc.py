## calc timeconstant

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
import pdb
from scipy.optimize import curve_fit
import time
import datetime

clrlst = ['b','r','g','k','deeppink','cyan','magenta']
xlbl = 'Frame(n)'
ylbl = '$\Delta$ ADU'

def func(xx,aa,bb):
    return aa*np.exp(xx*bb)

def calc_timeconstant(data,ridx,offidx,rampoffset,showfigs,savefigs,writefile,flipdata=False):
    ## Creating lists to append, and hard coding specific indicies for ramp normalization
    if writefile:
        t0 = time.time()
        text = "\n\nTime Constant Info:"
    avg_ledon = []
    onrange = [129,384]
    avg_ledoff = []
    offrange = [513,768]
    delta_onoff = []
    reduced_ramps = []

    ## Looping to normalize each ramp individually
    ## Must be len(ridx)-1 because it will fail when extracting temp list
    for ii in range(len(ridx)-1):
        ## Creating a temporary list from the reduced data to normalize
        temp = data[ridx[ii]:ridx[ii+1]]

        ## Extracting the Delta ADU values for LED-on / LED-off to compute averages
        onlist = temp[onrange[0]:onrange[1]]    
        offlist = temp[offrange[0]:offrange[1]]
        avg_ledon.append(np.mean(onlist))                 ## computes avg_ledon for all three ramps [129-384]
        avg_ledoff.append(np.mean(offlist))               ## computes avg_ledoff for all three ramps [513-768]
        delta_onoff.append(avg_ledon[ii]-avg_ledoff[ii])  ## computes delta_onoff for all three ramps
        for jj in range(len(temp)):
            temp[jj] -= avg_ledoff[ii]                    ## subtracts avg_ledoff from entire ramp
            temp[jj] /= delta_onoff[ii]

        ## Appending the list
        reduced_ramps.append(temp)

    ## Because python loops up-to-but-not-including the last ramp, this block of code is necessary to extract the last, not included, ramp
    #####################################################
    nlastramp = len(data[ridx[ii]:ridx[ii+1]])          #
    lasttemp = data[ridx[ii+1]:ridx[ii+1]+nlastramp]    #
    last_onlist = lasttemp[onrange[0]:onrange[1]]       #
    last_offlist = lasttemp[offrange[0]:offrange[1]]    #
    avg_ledon.append(np.mean(last_onlist))              #
    avg_ledoff.append(np.mean(last_offlist))            #
    delta_onoff.append(avg_ledon[-1]-avg_ledoff[-1])    #
    for ii in range(len(lasttemp)):                     #
        lasttemp[ii] -= avg_ledoff[-1]                  #
        lasttemp[ii] /= delta_onoff[-1]                 #
    reduced_ramps.append(lasttemp)                      #
    #####################################################

    ## creating list to append to
    rampfit = []
    flipped = []
    alphalist = []
    betalist = []

    print(' Calculating ...')
    for ii in range(len(reduced_ramps)):
        if flipdata:
            ## If executed, then the software is calculating time constant for LED-onset.
            ## It must flip to polarity of the data to fit the exponential.
            print('\t... Onset Time constant for ramp #{0} (i.e., R{1}) ... '.format(ii+1,ii))
            temp = reduced_ramps[ii][:rampoffset[ii]]
            for jj in range(len(temp)):
                temp[jj] -= 1.
                temp[jj] *= -1.
            flipped.append(temp)

        elif flipdata==False:
            ## If executed, then the software is calculating time constant for LED-offset.
            ## Regular exponential fitting does not require polarity flipping.
            temp = reduced_ramps[ii][rampoffset[ii]:]
            print('\t... Offset Time constant for ramp #{0} (i.e., R{1}) ... '.format(ii+1,ii))
        ## Scanity check
        else: 
            print("Parameter 'flipdata' not defined. Shown data may not represent what is expected.")
        
        ## Instantiating a frame array, from 0, 1, 2, ... , len(temp) as index 
        nframes = np.arange(0,len(temp))

        ## Calling the func() function to fit the exponential
        f = curve_fit(func,nframes,temp,p0=[-10e-5,10e-5])[0]

        ## Extracting alpha and beta constants defined in func() function
        a,b = f[0],f[1]

        ## Appending a list of the fitted exponential for overplotting for this specific ramp
        rampfit.append(func(nframes,a,b))

        ## Appending individiaul parameters for labels / extractions
        alphalist.append(a)
        betalist.append(b)

    ## Plotting
    if showfigs or savefigs:
        print("\tPlotting...")
    ## Plotting ramp subplots
    plt.clf()
    gs = gridspec.GridSpec(nrows=len(reduced_ramps),ncols=1)
    for ii in range(len(reduced_ramps)): 
        ax = plt.subplot(gs[ii,0])
        plt.plot(reduced_ramps[ii],c=clrlst[ii],linewidth=0.5,label='R{} Reduced'.format(ii))
        plt.axvline(x=rampoffset[ii],c=clrlst[ii],linestyle="-.",linewidth=0.5)
        plt.legend(prop={'size': 6})
    plt.xlabel(xlbl)
    if savefigs:
        print("\t\t... figure saved as 'rampsubplots.pdf'")
        if writefile:
            text += "\n\t\t - rampsubplots.pdf"
        plt.savefig('rampssubplots.pdf')
    if showfigs:
        plt.show()

    ## plotting linear scale in subplots
    plt.clf()
    gs = gridspec.GridSpec(nrows=len(reduced_ramps),ncols=1)
    if flipdata:
        for ii in range(len(reduced_ramps)):
            ax = plt.subplot(gs[ii,0])
            plt.plot(flipped[ii],c=clrlst[ii],linewidth=0.5,label='$R{}^*$ Reduced (flipped)'.format(ii))
            plt.axvline(x=0,c=clrlst[ii],linestyle='-.',linewidth=0.5,label='$R{}^*$ LED-Onset'.format(ii))
            plt.plot(rampfit[ii],c='k',linestyle='--',linewidth=0.5,label='R${0}^*$ fit ($\\alpha,\\beta$)=({1:.4f},{2:.4f})'.format(ii,alphalist[ii],betalist[ii]))
            plt.legend(prop={'size':6})
        plt.xlabel(xlbl)
        if savefigs:
            print("\t\t... figure saved as 'offon_tc.pdf'")
            if writefile:
                text += "\n\t\t - offon_tc.pdf"
            plt.savefig('offon_tc.pdf')
        if showfigs:
            plt.show()

        ## plotting flipped data
        plt.clf()
        for ii in range(len(reduced_ramps)):
            fit = rampfit[ii]
            for jj in range(len(fit)):
                fit[jj] -= 1.
                fit[jj] *= -1.
            ax = plt.subplot(gs[ii,0])
            plt.plot(reduced_ramps[ii][:rampoffset[ii]],c=clrlst[ii],linewidth=0.5,label='$R{}^*$ Reduced'.format(ii))
            plt.axvline(x=len(fit),c=clrlst[ii],linestyle='-.',linewidth=0.5,label='$R{}^*$ LED-Onset'.format(ii))
            plt.plot(fit,c='k',linestyle='--',linewidth=0.7,label='R${0}^*$ fit ($\\alpha,\\beta$)=({1:.4f},{2:.4f})'.format(ii,alphalist[ii],betalist[ii]))
            plt.legend(prop={'size':6})
        plt.xlabel(xlbl)
        if savefigs:
            print("\t\t... figure saved as offon_tc_flipped.pdf")
            if writefile:
                text += "\n\t\t - offon_tc_flipped.pdf"
            plt.savefig('offon_tc_flipped.pdf')
        if showfigs:
            plt.show()
        

    else:
        plt.clf()
        for ii in range(len(reduced_ramps)):
            ax = plt.subplot(gs[ii,0])
            plt.plot(reduced_ramps[ii][rampoffset[ii]:],c=clrlst[ii],linewidth=0.5,label='R{} Reduced'.format(ii))
            plt.axvline(x=0,c=clrlst[ii],linestyle='-.',linewidth=0.5,label='R{} LED-Offset'.format(ii))
            plt.plot(rampfit[ii],c='k',linestyle='--',linewidth=0.5,label='R{0} fit ($\\alpha,\\beta)$=({1:.4f},{2:.4f})'.format(ii,alphalist[ii],betalist[ii]))
            plt.legend(prop={'size':6})
        plt.xlabel(xlbl)
        if savefigs:
            print("\t\t... figure saved as 'onoff_tc.pdf'")
            if writefile:
                text += "\n\t\t - onoff_tc.pdf"
            plt.savefig('onoff_tc.pdf')
        if showfigs:
            plt.show()

    if writefile:
        return betalist,text
    else:
        return betalist
