## template #1

## Packages:
import numpy as np
import pdb
import time
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
print("Normalizing Data...")
for ii in range(0,d1-1,1):
    print('\t\tNormalizing Iteration: {} / {}'.format(ii,d1-1))
    ## Extracting array dimensions for easier calculations and cleaner code
    frame_med_norm[ii] = np.median(tdata[ii,:,:])-np.median(tdata[ii+1,:,:])
print("\t...Data normalized.")

frame_reduced = np.zeros(d1)
frameavg = np.zeros((d1,d2,d3)) ## Extracting/duplicating frameset and size
print("Reducing Data ... ")
t0 = time.time()
tt = np.zeros(d1-plusminus)
for nframe in range(plusminus,d1-plusminus,1):
    if nframe%10==0:
        t1 = time.time()
        tt[nframe] = t1-t0
        print("\tIteration: {} / {} \t\t ETR: {} seconds".format(nframe,d1-plusminus,int((tt[nframe])*(d1-nframe-plusminus)/10)))
        t0 = time.time()
    frame_reduced[nframe] = np.mean(frame_med_norm[nframe-plusminus:nframe+plusminus])
    for jj in range(d2):
        for kk in range(d3):
            frameavg[nframe,jj,kk] = np.mean(tdata[nframe-plusminus:nframe+plusminus,jj,kk])

red2 = np.zeros(d1)
for nframe in range(0,d1-1,1):
    red2[nframe] = np.median(frameavg[nframe])-np.median(frameavg[nframe+1])
print("\t...Data reduced.")

print('... Plotting ...')
plt.clf()
gs = gridspec.GridSpec(nrows=nsubplots,ncols=1)

for ii in range(0,nsubplots,1):
    ax = plt.subplot(gs[ii,0])
    plt.plot(frame_med_norm[int(ii*d1/nsubplots):int((ii+1)*(d1/nsubplots))],linewidth=0.2,c='r',label='Normalized')
    plt.plot(frame_reduced[int(ii*d1/nsubplots):int((ii+1)*(d1/nsubplots))],linewidth=0.5,c='k',label='Reduced')
    plt.plot(red2[int(ii*d1/nsubplots):int((ii+1)*(d1/nsubplots))],linewidth=0.5,c='b',label='Red2')

    plt.xlim(left=0,right=d1/nsubplots)
    if ii != nsubplots-1:
        ax.xaxis.set_visible(False)
plt.subplots_adjust(wspace=0,hspace=0,top=1.,bottom=0.026,left=0.026,right=0.989)

plt.legend()
plt.savefig('{}_normrednew.pdf'.format(dfile[:-5])) ## Cutting off *.fits extension
plt.show()
