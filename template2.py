## template2

import pdb
import pandas as pd
import numpy as np
from astropy.io import fits
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

def read_file(filename = '180925_LED_pulse_3.0V_1000f_cube.fits'):
    fits_file = fits.open(filename)
    holder = fits_file[0].data
    cat1 = holder[:,0,:]
    cat = pd.DataFrame(cat1)
    return cat

def plot(meds_norm, sum_arr, nrows = 4, ncols = 1):
    gs = gridspec.GridSpec(nrows, ncols)
    split_reduced = np.array_split(sum_arr, nrows)
    split_norm = np.array_split(meds_norm,nrows)

    for ii in range(0, nrows):
        ax = plt.subplot(gs[ii,0])
        plt.plot(split_norm[ii], linewidth=0.3, c='r')
        plt.plot(split_reduced[ii], linewidth=0.7, c='k')

    plt.subplots_adjust(hspace=0)
    plt.show()


if __name__ == '__main__':

    df = read_file()
    meds = df.median(axis = 1).values
    meds_norm = np.array([meds[ii]-meds[ii+10] for ii in range(0,len(meds)-10)])
    sum_arr = np.array([meds_norm[ii-16:ii+15].mean() for ii in range(16, len(meds_norm)-15)])
    plot(meds_norm, sum_arr)
    
