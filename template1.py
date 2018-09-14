## template #1

## Packages:
import numpy as np
import pdb
from astropy.io import fits
import matplotlib.pyplot as plt

## Constants:



## Code:

dfile='testfile.fits'
tfile = fits.open(dfile)

## Semi-assuming were working with a data from of 1 (technically zero; check with tdata.info() should only have one row)
tdata = tfile[0].data

## Extracting dimensions:
d1,d2,d3 = tdata.shape[0],tdata.shape[1],tdata.shape[2] ## d1,d2,d3==499000,1,32 or something along those lines
xidx = np.arange(0,d1,1) ## Setting x-axis to integer
ymed = np.empty(d1)
yavg = np.empty(d1)
for ii in range(d1):
    ymed[ii] = np.median(tdata[ii,:])
    yavg[ii] = np.mean(tdata[ii,:])
    if ii%1000==0:
        print("\tii:\t{}".format(ii))

error = np.abs(ymed-yavg)
pdb.set_trace()
plt.plot(xidx,ymed,linewidth=0.1,c='b',label='Median')
plt.plot(xidx,yavg,linewidth=0.1,c='g',label='Mean')

plterror = input("\n\tPlot Error? (y or n)")
if plterror=='y':
    plt.plot(xidx,error,linewidth=0.1,c='r',label='$\delta$')
#plt.plot(xidx,np.median(yaxis),linewidth=0.5,c='r',label='Median')
plt.xlabel("Index (increment of n==1)")
#plt.ylabel("np.median(data[ii,:])")
plt.legend()
plt.show()







pdb.set_trace() #Good habit to have set_trace()
