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
dfile='180928_LEDonoff_8000f_6.6Vbias_cube.fits'
nsubplots = 4
## Code:
hdulist = fits.open(dfile)
## Semi-assuming were working with a data from of 1 (technically zero; check with hdulist.info() should only have one row)
tdata = hdulist[0].data

## Extracting dimensions:
d1,d2,d3 = tdata.shape[0],tdata.shape[1],tdata.shape[2] ## d1,d2,d3==499000,1,32 or something along those lines
xidx = np.arange(0,d1,1) ## Setting x-axis to integer
frame_med_norm = np.empty(d1)

for ii in range(0,d1-1,1):
    ## Extracting array dimensions for easier calculations and cleaner code
    frame_med_norm[ii] = np.median(tdata[ii,:,:])-np.median(tdata[ii+1,:,:])
    ## plot raw values
    if ii%100==0:
        print("\tii:\t{} / {}".format(ii,d1))

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
    plt.plot(frame_med_norm[int(ii*d1/nsubplots):int((ii+1)*(d1/nsubplots))],linewidth=0.3,c='r',label='Normalized')
    plt.plot(frame_reduced[int(ii*d1/nsubplots):int((ii+1)*(d1/nsubplots))],linewidth=0.5,c='k',label='Reduced')
    plt.xlim(left=0,right=d1/nsubplots)
    if ii != nsubplots-1:
        ax.xaxis.set_visible(False)
plt.subplots_adjust(hspace=0)
#plt.ylim(top=100.,bottom=-125.)
plt.legend()
plt.tight_layout()
plt.savefig('{}_normred.pdf'.format(dfile[:-5]),dpi=10000) ## Cutting off *.fits extension
plt.show()
pdb.set_trace() #Good habit to have set_trace()