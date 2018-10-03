## template #1

## Packages:
import numpy as np
import pdb
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
 ###
 ### Remove array slice to get to all 64 pixels 
 ###
## Constants:
nframe_avg = 32
plusminus = int(nframe_avg/2) ## for "plus or minus" in looping
dfile='181003_LEDonoff_8000f_6.6Vbias_cube.fits'
nsubplots = 4

## Code:
hdulist = fits.open(dfile)
## Semi-assuming were working with a data from of 1 (technically zero; check with hdulist.info() should only have one row)
tdata = hdulist[0].data

## Extracting dimensions:
d1,d2,d3 = tdata.shape[0],tdata.shape[1],tdata.shape[2] ## d1,d2,d3==499000,1,32 or something along those lines
xidx = np.arange(0,d1,1) ## Setting x-axis to integer
frame_med_norm = np.empty(d1)
frame_row_med = np.empty((d1,d2))

## Set_trace() here to check new data shape
pdb.set_trace()

## each row of 64, median to get a 256 set of medians for each frame
## subtract 265-n+1 from 265-n frame, then average


for ii in range(0,d1-1,1):
    ## Extracting array dimensions for easier calculations and cleaner code
    frame_med_norm[ii] = np.median(tdata[ii,:,:])-np.median(tdata[ii+1,:,:])
    for jj in range(0,d2-1,1):
        frame_row_med[ii,jj] = np.mean(tdata[ii,jj,:])

    ## plot raw values
    if ii%100==0:
        print("\tii:\t{} / {}".format(ii,d1))

tino = np.array([frame_row_med[ii]-frame_row_med[ii+1] for ii in range(0,d1-1,1)])
tinoavg = np.array([np.median(tino[zz-plusminus:zz+plusminus-1]) for zz in range(plusminus,d1-plusminus+1,1)])
for ii in range(plusminus): ## must insert first plusminus indicies as zeros or data does not overlap with reduced
    tinoavg = np.insert(tinoavg,ii,0.)

plt.clf()
gs = gridspec.GridSpec(nrows=nsubplots,ncols=1)
frame_reduced = np.zeros(d1)
## Starting loop at p_m; throw away edges of plot
for frame_set in range(plusminus,d1-plusminus,1):
    frame_reduced[frame_set] = np.mean(frame_med_norm[frame_set-plusminus:frame_set+plusminus-1])
    if frame_set%100==0:
        print("\tframe_set:\t{} / {}".format(frame_set,d1-plusminus))
plt.clf()
for ii in range(0,nsubplots,1):
    ax = plt.subplot(gs[ii,0])
    plt.plot(frame_med_norm[int(ii*d1/nsubplots):int((ii+1)*(d1/nsubplots))],linewidth=0.2,c='r',label='Normalized')
    plt.plot(frame_reduced[int(ii*d1/nsubplots):int((ii+1)*(d1/nsubplots))],linewidth=0.5,c='k',label='Reduced')
    #plt.plot(tinoavg[int(ii*d1/nsubplots):int((ii+1)*(d1/nsubplots))],linewidth=0.5,c='b',label='New')

    plt.xlim(left=0,right=d1/nsubplots)
    if ii != nsubplots-1:
        ax.xaxis.set_visible(False)
plt.subplots_adjust(wspace=0,hspace=0,top=1.,bottom=0.026,left=0.026,right=0.989)

plt.legend()
#plt.tight_layout()
plt.savefig('{}_normrednew.pdf'.format(dfile[:-5])) ## Cutting off *.fits extension
plt.show()
##pdb.set_trace() #Good habit to have set_trace()