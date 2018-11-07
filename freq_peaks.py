import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from scipy.signal import find_peaks
import pdb
import time

if __name__ == '__main__':
    t0 = time.time()
    n_order = 79
    y = pd.read_csv('freq_data.csv')

    # Correct weird formatting shit
    y = y.drop(['Unnamed: 0'], axis = 1)

    # Turn into 1D Array
    y = y.values.flatten()

    # Indice Array
    x = np.arange(0,len(y),1)

    # Fit polynomial of order n_order
    z = np.polyfit(x,y,n_order)

    # Function representing the fit polynomial
    f = np.poly1d(z)

    # y-values of that polynomial
    y_pred = f(x)

    # Indicies of peaks
    peaks = find_peaks(y_pred)

    # x_peaks is list of indices of peaks of polynomial fit to freq_data
    x_peaks = [x[p] for p in peaks[0]]
    
    print(x_peaks)
    t1 = time.time()
    print('{0} s to find peaks'.format(t1 - t0))
    
    pdb.set_trace()
    

