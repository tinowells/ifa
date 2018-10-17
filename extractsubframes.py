## Extract subframes

import numpy as np
import pdb
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import fits


dfile = '181003_LEDonoff_8000f_6.6Vbias_cube.fits'
dir_path = 'ExtractSubframes/'
stat_dir_path = ['median/','mean/','both/']
stat_name = ['Median','Mean','Median (blue) & Mean (red)']
pltclrlst = ['b','r']

xlbl = 'n$^{th}$ Frame'
ylbl = '$\\Delta$ ADU'
nsubplots = 4

calc_type = 0 # zero for median, 1 for mean, 2 for both

def read_data(fitsfile):
    hdulist = fits.open(fitsfile)
    frame_data = hdulist[0].data
    ## Extracting first 1000 frames to run faster, will delete when compiled
    #framedata = frame_data[0:1000,:,:]
    return frame_data

def extractdata(rdata):
    nframes,npix,ncols = rdata.shape[0],rdata.shape[1],rdata.shape[2]
    edata = rdata[0:nframes,int(npix/4):int(3*npix/4),int(ncols/2):]
    return edata

def plot_columns(data):
    frameidx = np.arange(0,data.shape[0],1)
    for ii in range(data.shape[2]): ##0-31
        colvals = []
        if calc_type == 2:
            plotmedmean(data,frameidx,ii)
        elif calc_type == 0:
            for jj in frameidx: ##0-999
                colvals.append(np.median(data[jj,ii,:]))
        elif calc_type == 1:
            for jj in frameidx:
                colvals.append(np.mean(data[jj,ii,:]))

        if calc_type != 2:
            plt.clf()
            gs = gridspec.GridSpec(nrows=nsubplots,ncols=1)
            for kk in range(0,nsubplots,1):
                ax = plt.subplot(gs[kk,0])
                plt.plot(colvals[int(kk*len(colvals)/nsubplots):int((kk+1)*len(colvals)/nsubplots)],linewidth=0.25,c=pltclrlst[calc_type])
                if kk == 0:
                    plt.title('{};    ncol = {}'.format(dfile,ii))
                plt.xlabel(xlbl)
                if kk != nsubplots-1:
                    ax.xaxis.set_visible(False)
            #ax.text(-160.,75.,ylbl,ha='center',va='center',rotation='vertical')
            plt.subplots_adjust(wspace=0,hspace=0,top=0.94,bottom=0.1,left=0.1,right=0.975)
            plt.savefig(dir_path+stat_dir_path[calc_type]+'ncol{}.pdf'.format(ii))
            print('Plotting ii = {}\t\tcalc_type == {}'.format(ii,stat_name[calc_type]))
            #pdb.set_trace()


    done = '\tDone! Check directory for figures.'
    return done

def plotmedmean(bdata,bidx,icounter):
    medvals = []
    meanvals = []
    for jj in bidx:
        medvals.append(np.median(bdata[jj,icounter,:]))
        meanvals.append(np.mean(bdata[jj,icounter,:]))

    plt.clf()
    gs = gridspec.GridSpec(nrows=nsubplots,ncols=1)
    for kk in range(0,nsubplots,1):
        ax = plt.subplot(gs[kk,0])
        plt.plot(medvals[int(kk*len(medvals)/nsubplots):int((kk+1)*(len(medvals)/nsubplots))],linewidth=0.25,c=pltclrlst[0],label=stat_name[0])
        plt.plot(meanvals[int(kk*len(meanvals)/nsubplots):int((kk+1)*(len(meanvals)/nsubplots))],linewidth=0.25,c=pltclrlst[1],label=stat_name[1])
        if kk == 0:
            plt.title('{};    ncol = {}'.format(dfile,icounter))      
        plt.xlabel(xlbl)
        if kk != nsubplots-1:
            ax.xaxis.set_visible(False)
    plt.subplots_adjust(wspace=0,hspace=0,top=0.94,bottom=0.1,left=0.1,right=0.975)
    plt.legend(prop={'size': 6})
    plt.savefig(dir_path+stat_dir_path[calc_type]+'ncol{}.pdf'.format(icounter),dpi=1000)
    print('Plotting icounter = {}\t\tcalc_type == {}'.format(icounter,stat_name[calc_type]))            
    return 


if __name__=='__main__':
    rawdata = read_data(dfile)
    extracted_data = extractdata(rawdata)
    tino = plot_columns(extracted_data)
    print(tino)


    ## read in data
    ## extract subframes
    ## create figure for each column (put median/mean in **kwargs)
    ## - save figure in directory
    ## - Write to textfile list of figure files names
    ##