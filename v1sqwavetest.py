## Updates square wave test

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from scipy.signal import find_peaks
import pdb
import time

rampclr = ['r','b','g','deeppink','cyan','magenta','k']

## This assumes data already reduced and saved to *.csv file.
dfile = 'freq_data.csv'
n_order = 79

def read_data(data_file):
    ddata = pd.read_csv(data_file)
    ddata = ddata.drop(['Unnamed: 0'], axis = 1)
    return ddata

def polyfit_peaks(y):
    y = y.values.flatten() # Turn into 1D Array
    x = np.arange(0,len(y),1) # Indice Array
    z = np.polyfit(x,y,n_order) # Fit polynomial of order n_order
    f = np.poly1d(z) # Function representing the fit polynomial
    y_pred = f(x) # y-values of that polynomial
    peaks = find_peaks(y_pred) # Indicies of peaks
    x_peaks = [x[p] for p in peaks[0]] # x_peaks is list of indices of peaks of polynomial fit to freq_data
    pdb.set_trace()
    return x_peaks,y_pred


if __name__=='__main__':
    data = read_data(dfile)
    xidx = np.arange(0,len(data),1)
    polypeaks,poly = polyfit_peaks(data)
    ## Plotting
    plt.plot(data)
    plt.plot(xidx,poly)
    for ii in range(len(polypeaks)):
        #plt.plot(data[polypeaks[ii]:polypeaks[ii+1]],c=rampclr[ii])
        plt.axvline(x=polypeaks[ii],linestyle='--',c=rampclr[ii])
    plt.show()
    
    #ramp_data = []
    #for ii in range(len(polypeaks)-1):
    #    ramp_data.append(data[polypeaks[ii]:polypeaks[ii+1]])

    pdb.set_trace()