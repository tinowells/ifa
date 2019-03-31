## testing to extract and create subarrays, create figures and save files
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import scipy.special as sspecial
import warnings
warnings.filterwarnings("error")
import pdb
import os
import shutil ## for moving files
import glob
import argparse

parser = argparse.ArgumentParser(description='Reduction Pipeline Parameters:')
parser.add_argument('-ddir', '--datadirectory', metavar='', required=True, type=str)
parser.add_argument('-sbin', '--sizeofbins', metavar='', required=False, type=int, default=2)
parser.add_argument('-imupper', '--imshowupperlimit', metavar='', required=False, type=int, default=15)
parser.add_argument('-imlower', '--imshowlowerlimit', metavar='', required=False, type=int, default=-35)
parser.add_argument('-histupper', '--upperhistlimit', metavar='', required=False, type=int, default=5)
parser.add_argument('-histlow', '--lowerhistlimit', metavar='', required=False, type=int, default=-40)
parser.add_argument('-c_map', '--colormap', metavar='', required=False, type=str, default='jet')
args = parser.parse_args()
datadirectory = args.datadirectory
binsize = args.sizeofbins
imshowlim0 = args.imshowupperlimit
imshowlim1 = args.imshowlowerlimit
imshowlim = [imshowlim0, imshowlim1]
histlim1 = args.upperhistlimit
histlim0 = args.lowerhistlimit
histlim = [histlim0, histlim1]
c_map = args.colormap


## Dimensions for subarray extractions
xdim = [[67,93],[99,125],[131,157],[163,189],[195,221],[227,253]] ## 320 axis
ydim = [65,192] ## 256 axis
nsubarr = len(xdim)
biaslist = ['2.5', '4.5', '6.5', '7.5', '8.5', '9.5', '10.5', '11.5']

def extract_subarrays(fullarray):
    hdu = fits.open(fullarray)
    fullarr_data = hdu[0].data
    for ii in range(nsubarr):
        subarr = fullarr_data[ydim[0]:ydim[1],xdim[ii][0]:xdim[ii][1]]
        fits.writeto('{}/subarray_{}_{}_{}.fits'.format(os.getcwd(),file[0:-5],xdim[ii][0],xdim[ii][1]),subarr)

def sort_filegroup(group):
    holder = group
    holder.insert(1,holder[-1])
    holder.insert(1,holder[-2])
    holder.pop(-1)
    holder.pop(-1)
    return holder

def get_fits_data(ffile):
    hdu = fits.open(ffile)
    fimg = hdu[0].data
    fdata = fimg.flatten()
    hdu.close()
    return fdata, fimg

def get_bins(datalst):
    range = datalst.max()-datalst.min()
    number_bins = range/binsize
    return int(number_bins), range

def shift_to_center(ends):
    centers = ends - np.diff(ends)[0]/2.
    return centers[1:] ## Because of how were adjust

def plot_fullarray(fimage,numbins,xfit,yfit,centers,counts,mu,sig,savefilename):
    gs = gridspec.GridSpec(nrows=2, ncols=1)
    figure = plt.gcf() # get current figure
    figure.set_size_inches(13, 7)
    ax = plt.subplot(gs[0,0])
    plt.imshow(fimage, cmap=c_map, vmin=imshowlim[0], vmax=imshowlim[1])
    plt.colorbar(shrink=0.9)
    plt.ylim(top=65, bottom=192)
    plt.xlabel('Pixel Index', size=8)
    plt.ylabel('Pixel Index', size=8)

    ax = plt.subplot(gs[1,0])
    plt.hist(fimage.flatten(), bins=numbins, color='grey')
    plt.plot(centers, counts, color='r', label='Counts')
    plt.plot(xfit, yfit, ls=':', color='k', label='Fit')
    plt.axhline(y=np.max(yfit), color='orange', linewidth=0.75, label='Peak: {}'.format(int(np.max(yfit))))
    plt.axvline(x=mu, color='b', label='$\mu = {0:.2f}$'.format(mu))
    plt.axvline(x=mu+sig, color='m', label='$\sigma = {0:.2f}$'.format(sig))
    plt.axvline(x=mu-sig, color='m')
    plt.axvline(x=mu+2*sig, color='g')
    plt.axvline(x=mu-2*sig, color='g')
    plt.xlim(left=histlim[0], right=histlim[1])
    plt.legend(prop={'size':14},loc='center left')
    plt.subplots_adjust(left=0.05,bottom=0.07,right=0.99,top=0.93)
    plt.savefig('{}'.format(savefilename))

def plot_subarrays(simagelst, numbinslst, xfitlst, yfitlst, centerslst, countslst, mulst, siglst, subsavefilename):
    plt.clf()
    gs = gridspec.GridSpec(nrows=2, ncols=len(simagelst))
    figure = plt.gcf() # get current figure
    figure.set_size_inches(13, 7)
    for jj in range(len(simagelst)):
        simage = simagelst[jj]
        sbins = numbinslst[jj]
        sxfit, syfit = xfitlst[jj], yfitlst[jj]
        scenters, scounts = centerslst[jj], countslst[jj]
        smu, ssig = mulst[jj], siglst[jj]

        ax = plt.subplot(gs[0,jj])
        plt.imshow(simage, cmap=c_map, vmin=imshowlim[0], vmax=imshowlim[1])
        plt.colorbar(shrink=0.9)
        plt.xlabel('Pixel Index', size=8)
        plt.ylabel('Pixel Index', size=8)

        ax = plt.subplot(gs[1,jj])
        plt.hist(simage.flatten(), bins=sbins, color='grey')
        plt.plot(scenters, scounts, color='r')
        plt.plot(sxfit, syfit, ls=':',color='k')
        plt.axhline(y=np.max(syfit), color='orange', label='Peak: {}'.format(int(np.max(syfit))))
        plt.axvline(x=smu, color='b', label='$\mu = {0:.2f}$'.format(smu))
        plt.axvline(x=smu+ssig, color='m', label='$\sigma = {0:.2f}$'.format(ssig))
        plt.axvline(x=smu-ssig, color='m')
        plt.axvline(x=smu+2*ssig, color='g')
        plt.axvline(x=smu-2*ssig, color='g')
        plt.xlim(left=histlim[0],right=histlim[1])
        plt.legend(prop={'size':8},loc='upper left', bbox_to_anchor=(-0.25, 1.1),framealpha=1.)
        if ii == 0:
            plt.ylabel('Pixel Counts',size=8)
            plt.xlabel('Pixel Values',size=8)
    plt.subplots_adjust(left=0.05,bottom=0.07,right=0.97,top=0.96,wspace=0.25,hspace=0.1)
    plt.savefig('{}'.format(subsavefilename))

def get_file_name(file):
    ## bias, up/down/top
    settled = ['up', 'down', 'top']
    bbias, ssettle = '_', '_'
    ## Finding the bias voltage
    for bias in biaslist:
        if str.__contains__(file, bias):
            bbias = bias
    ## Finding if it is settled up, down or top
    for settle in settled:
        if str.__contains__(file, settle):
            ssettle = settle
    ## Finding if it is a subarray or full array
    if str.__contains__(file, 'subarray'):
        arraytype = 'subarray'
    else:
        arraytype = 'fullarray'
    ## Meaning information was not found
    if bbias == '_' or ssettle == '_':
        print('\n\tFILENAME NOT FOUND: returning regular file name')
        return file[:-5]
    else:
        fname = '{0}V-Bias_{1}_{2}.pdf'.format(bbias, ssettle, arraytype)
        return fname

def fit_gaussian(x,alpha,mu,sigma,offset):
    return alpha*np.exp(-(x-mu)**2/(2*sigma**2)) + offset

def main(filelst):
    ## Fitting Gaussian to full array First
    fullarray, fullimage = get_fits_data(filelst[0])
    nbins, range = get_bins(fullarray)
    full_np_hist = np.histogram(fullarray, bins=nbins)
    full_counts, full_ends = full_np_hist[0].tolist(), full_np_hist[1].tolist()
    full_centers = shift_to_center(full_ends)
    try:
        full_params, full_pconv = curve_fit(fit_gaussian, full_centers, full_counts, p0=[1, np.mean(full_counts), np.std(full_counts,ddof=1), 0.], maxfev=1000000000)
    except Exception:
        print('\tRe-fitting full array ...')
        diff = np.diff(full_counts)
        checkidxs = np.where(diff >= 150)[0].tolist()
        goodfit = False
        counter = 0
        while goodfit == False and counter < len(checkidxs):
            checkidx = checkidxs[counter]
            try:
                full_params, full_pconv = curve_fit(fit_gaussian, full_centers[checkidx:], full_counts[checkidx:], p0=[1, np.mean(full_counts[checkidx:]), np.std(full_counts[checkidx:],ddof=1), 0.], maxfev=1000000000)
                print(' ... extracting better fit.')
                goodfit = True
            except Exception:
                counter += 1
                goodfit = False
    full_xfit = np.linspace(full_centers[0], full_centers[-1], num=int(range*1000))
    full_yfit = fit_gaussian(full_xfit, *full_params)
    savefile = get_file_name(filelst[0])
    ## NOW PLOT FULL ARRAY INFO
    plot_fullarray(fullimage, nbins, full_xfit, full_yfit, full_centers, full_counts, full_params[1], full_params[2], savefile)

    ## Begin fitting on Subarrays
    sub_imagelst = []
    sub_centerslst, sub_countslst = [], []
    sub_xfitlst, sub_yfitlst = [], []
    sub_mulst, sub_siglst = [], []
    sub_binslst = []

    for subidx in np.arange(1,len(filelst)):
        if subidx == 1: ## Meaning the first subarray
            x0, x1, x2, x3 = 1, full_params[1], full_params[2], full_params[3]
        else:
            x0, x1, x2, x3 = sub_params[0], sub_params[1], sub_params[2], sub_params[3]
        subarr, subimg = get_fits_data(filelst[subidx])
        nbins, range = get_bins(subarr)
        sub_np_hist = np.histogram(subarr, bins=nbins)
        sub_counts, sub_ends = sub_np_hist[0].tolist(), sub_np_hist[1].tolist()
        sub_centers = shift_to_center(sub_ends)

        try:
            sub_params, sub_pconv = curve_fit(fit_gaussian, sub_centers, sub_counts, p0=[x0, x1, x2, x3], maxfev=1000000000)
        except Exception:
            print("\tRefitting subarray ....")
            diff = np.diff(sub_counts)
            checkidxs = np.where(diff >= 150)[0].tolist()
            goodfit = False
            counter = 0
            while goodfit == False and counter < len(checkidxs):
                checkidx = checkidxs[counter]
                try:
                    sub_params, sub_pconv = curve_fit(fit_gaussian, sub_centers[checkidx:], sub_counts[checkidx:], p0=[1, np.mean(sub_counts[checkidx:]), np.std(sub_counts[checkidx:],ddof=1), 0.], maxfev=1000000000)
                    print(' ... extracting better fit.')
                    goodfit = True
                except Exception:
                    counter += 1
                    goodfit = False
            if counter == len(checkidxs):
                print('SUBPARMETERS COULD NOT BE MET: RETURNING ONES')
                sub_params = [1, 1, 1, 1]
            #pdb.set_trace()

        sub_xfit = np.linspace(sub_centers[0], sub_centers[-1], num=int(range*1000))
        sub_yfit = fit_gaussian(sub_xfit, *sub_params)

        sub_imagelst.append(subimg)
        sub_centerslst.append(sub_centers)
        sub_countslst.append(sub_counts)
        sub_xfitlst.append(sub_xfit)
        sub_yfitlst.append(sub_yfit)
        sub_mulst.append(sub_params[1])
        sub_siglst.append(sub_params[2])
        sub_binslst.append(nbins)
        if subidx == len(filelst)-1:
            subsavefile = get_file_name(filelst[subidx])
            plot_subarrays(sub_imagelst, sub_binslst, sub_xfitlst, sub_yfitlst, sub_centerslst, sub_countslst, sub_mulst, sub_siglst, subsavefile)
    return


if __name__=='__main__':
    os.chdir(datadirectory)
    #os.mkdir('fullarrays')
    #os.mkdir('subarrays')
    #os.mkdir('figures')
    fullarrayfiles = glob.glob('*.fits')
    fullarrayfiles.sort()
    for file in fullarrayfiles:
        extract_subarrays(file)
    subarrayfiles = glob.glob('subarray*fits')
    subarrayfiles.sort()
        #shutil.move(file,'fullarrays/')
    ## All files are now moved, need to get full/subarrays to fit guassians
    for ii in range(len(fullarrayfiles)):
        full_arr = fullarrayfiles[ii]
        file_group = []
        file_group.append(full_arr)
        for sub_arr in subarrayfiles:
            if sub_arr.__contains__(full_arr[:-5]):
                file_group.append(sub_arr)
        filegroup = sort_filegroup(file_group)
        main(filegroup)

    os.mkdir('fullarrays')
    os.mkdir('subarrays')
    os.mkdir('figures')
    figurefiles = glob.glob('*.pdf')

    for fullfile in fullarrayfiles:
        shutil.move(fullfile, 'fullarrays/')
    for subfile in subarrayfiles:
        shutil.move(subfile, 'subarrays/')
    for figfile in figurefiles:
        shutil.move(figfile, 'figures/')

    print('\n\n\nDONE\n\n\n')
