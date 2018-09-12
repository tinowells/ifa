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

## Printing dimesnions:
see_data = input("Would you like to print data? (y or n): ")

if see_data=='y':
    for ii in range(d1):
        print("ii:\t{0}\n\n{1}".format(ii,tdata[ii]))
else:
    print(' ')

nave = np.empty(d3)
for ii in range(d3):
    randint = np.random.randint(0,d1)
    print("Random int selected: {}".format(randint))
    plt.plot(tdata[randint,0,:],linewidth=0.5)
    nave[ii] = np.average(tdata[:,0,ii])
    

plt.plot(nave,".-",c='r',label="$\overline{nData}$",linewidth=0.5)
pdb.set_trace()
plt.legend()
plt.show()

