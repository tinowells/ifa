README.txt for sqwavetest.py

Square Wave Test:
Characterization of infrared sensore precedures for SAPHIRA avalanche photo
diode array; Python software for the acquisition, reduction, & analysis of 
SAPHIRA test data.

Written by Tino Wells <tinow@hawaii.edu>
Fall-2018, Spring-2019

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

sqwavetest.py:
 - This software was written for python 3.6.7. It should be noted it may not
    work for alternate versions of python, although the author wrote this 
    software with that in mind to avoid failures, and was designed to be scripted.

 - To run program, make sure all necessary files are in working directory and 
    the following:

        someone@somemachine: ~/some_path& python3 sqwavetest.py -fits some_fits_file.fits

    *** User must profide the *.fits file name or the program will not run.
    *** There are other arguements, but these are optional. See sqwavetest.py for 
            which each parameter is or run 'python3 sqwavetest.py -h' to prompt
            the help menu.


## Keys to sqwavetest:
 (1) Parse arguements
 (2) Reduce data
 (3) Detect peaks (changes in intensity)
 (4) Filter peaks to only onset/offset indicies
 (5) Extract ramp indicies, starting from LED-onset of ramp n, to LED-onset of ramp n+1
 (6) Calculate timeconstants for LED-onset / LED-offset
 (7) Plot & Save data / figures

(1) Parsing Arguments:
     - This is where the user tells the script what to do. Again, see sqwavetest.py or
        run 'python3 sqwavetest.py -h' to prompt help menu.

(2) Reduce Data: read_reduce()
     - The reduction pipeline for this script follows a CDS(1,1) reduction, meaning 
        it is hard coded to take frame[n+1] and subtract frame[n].
     - In read_reduce(dfile) in the script, we read in, CDS(1,1) reduce, and return a
        Pandas DataFrame to compress_df(dataframe).


(3) Detect Peaks: check_peakidx()
     - First excludes the ends of the data. These can be very noise, annoying to deal 
        with in the long haul, and they're not important anyways, so they're ignored.
     - Calculate the Delta ADU values from frame[n] to frame[n+npeakstep] (npeakstep 
        defaults to 50 frames). If the change in ADU values is larger than the 
        specified tolerance (tol), then it has detected a peak. 
     - The indicie of detection is appended to a list and the list is passed to 
        check_peakidx(data,peaks) to filter out noise peaks.
     - Some noise still presists in the selected peaks. The script then passes the 
        list of peaks to refine_peaks(data,peak_idx) to find the exact location where
        LED-onset / LED-offset takes place.

(4) Filter Peaks: refine_peaks()
     - This function looks around the index of interest and decides if the given index
        a peak or not. It does this by looking at previous indicies to see if there 
        are higher Delta ADU values BEHIND the index of interest.
     - If there are, it increments until this is no longer true.
     - Once this is no longer true, or if it was not true in the first place, it 
        breaks the loop and appends the index to a list of exact indicies.

(5) Extract Ramp Indicies: get_ramp_idx()
     - Getting the difference between onset indicies and creates a list.
     - Then appends a list of ramp offsets by defining half the number of indicies 
        onset index n to n+1 (defined by dramps w/ np.diff())
     - Then appends a list of indicies where each ramp starts, as defined by what index
        the LED turns on.

(6) Calculate Time Constants: timeconstantfunc.py
     - Extracts data from reduced data w/ ramp indicies
     - Appends a ramp list
        *** Adds the last ramp because python loops up-to-but-not-including the last iteration
     - Determines if data needs to be flipped (for fitting LED-onset time constant) or not 
        (for LED-offset time constant).
            *** If data needs to be flipped, it is for LED-onset. The polarity is flipped
                    to fit the negative exponential function defined by func().
     - Fits as an exponential decay and appends fitted data and parameters to a list.

(7) Plot & Save Data / Figures:
     - Basically just plots everything.
        1 - Ramp subplots with onset-index plotted as dashed vertical lines.
            *** THIS IS FLIPPED DATA
        2 - Overplotting fitted exponential decay with alpha and beta parameters
                in legend. Plots start at LED-onset index.
        3 - Ramp subplots of reduced onset data. 
            *** The exponential was fitted with flipped polarity. This figure is 
                    showing the original polarity, and overplotting the fitted
                    function.
        4 - Same as figure (1).
        5 - Overplotting fitted exponential decay with alpha and beta parameters
                in legend. Plot starts at LED-offset index. 
        6 - The entire data file in one figure. Extracted LED-onset and LED-offset
                shown as dashed green and red lines, respectively.
        7 - Overplotting all extracted ramps.
        8 - Ramps co-added to reduce noise, in a single figure.
        9 - Co-added mean value of ramps in single figure.




