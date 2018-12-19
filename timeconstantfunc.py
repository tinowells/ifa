## calc timeconstant

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
import pdb
from scipy.optimize import curve_fit

clrlst = ['b','r','g','k','deeppink','cyan','magenta']
xlbl = 'Frame(n)'
ylbl = '$\Delta$ ADU'

def func(xx,aa,bb):
    return aa*np.exp(xx*bb)

def calc_timeconstant(data,ridx,offidx,rampoffset,flipdata=False):
    avg_ledon = []
    onrange = [129,384] ## hard coded 
    avg_ledoff = []
    offrange = [513,768]
    delta_onoff = []
    reduced_ramps = []
    for ii in range(len(ridx)-1):
        temp = data[ridx[ii]:ridx[ii+1]]
        onlist = temp[onrange[0]:onrange[1]]    
        offlist = temp[offrange[0]:offrange[1]]
        avg_ledon.append(np.mean(onlist))                 ## computes avg_ledon for all three ramps [129-384]
        avg_ledoff.append(np.mean(offlist))               ## computes avg_ledoff for all three ramps [513-768]
        delta_onoff.append(avg_ledon[ii]-avg_ledoff[ii])  ## computes delta_onoff for all three ramps
        for jj in range(len(temp)):
            temp[jj] -= avg_ledoff[ii]                    ## subtracts avg_ledoff from entire ramp
            temp[jj] /= delta_onoff[ii]
        reduced_ramps.append(temp)

    ## adding last ramp
    nlastramp = len(data[ridx[ii]:ridx[ii+1]])
    lasttemp = data[ridx[ii+1]:ridx[ii+1]+nlastramp]
    last_onlist = lasttemp[onrange[0]:onrange[1]]
    last_offlist = lasttemp[offrange[0]:offrange[1]]
    avg_ledon.append(np.mean(last_onlist))
    avg_ledoff.append(np.mean(last_offlist))
    delta_onoff.append(avg_ledon[-1]-avg_ledoff[-1])
    for ii in range(len(lasttemp)):
        lasttemp[ii] -= avg_ledoff[-1]
        lasttemp[ii] /= delta_onoff[-1]
    reduced_ramps.append(lasttemp)

    ## plotting linear scale in subplots
    plt.clf()
    gs = gridspec.GridSpec(nrows=len(reduced_ramps),ncols=1)
    for ii in range(len(reduced_ramps)): 
        ax = plt.subplot(gs[ii,0])
        plt.plot(reduced_ramps[ii],c=clrlst[ii],linewidth=0.5,label='R{} Reduced'.format(ii))
        plt.axvline(x=rampoffset[ii],c=clrlst[ii],linestyle="-.",linewidth=0.5)
        plt.legend(prop={'size': 6})
    plt.xlabel(xlbl)
    plt.savefig('rampssubplots.pdf')
    plt.show()

    ## creating list to append to
    rampfit = []
    flipped = []
    alphalist = []
    betalist = []


    print(' Calculating ...')
    for ii in range(len(reduced_ramps)):
        if flipdata:
            print('\t... Onset Time constant for ramp #{0} (i.e., R{1}) ... '.format(ii+1,ii))
            temp = reduced_ramps[ii][:rampoffset[ii]]
            for ii in range(len(temp)):
                temp[ii] -= 1.
                temp[ii] *= -1.
            flipped.append(temp)
        elif flipdata==False:
            temp = reduced_ramps[ii][rampoffset[ii]:]
            print('\t... Offset Time constant for ramp #{0} (i.e., R{1}) ... '.format(ii+1,ii))
        else: 
            print("Parameter 'flipdata' not defined. Shown data may not represent what is expected.")
        
        nframes = np.arange(0,len(temp))
        f = curve_fit(func,nframes,temp,p0=[-10e-5,10e-5])[0]
        a,b = f[0],f[1]
        rampfit.append(func(nframes,a,b))
        alphalist.append(a)
        betalist.append(b)

    #pdb.set_trace()
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
        plt.savefig('offon_tc.pdf')
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
        plt.savefig('offon_tc_flipped.pdf')
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
        plt.savefig('onoff_tc.pdf')
        plt.show()

    return betalist
