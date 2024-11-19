import numpy as np
import csv
import os
import msvcrt
import fcntl
from scipy.signal import find_peaks

class PeakDetector:

    def __init__(self, message, static_thr, higher_ratio): #def __init__(self, message, static_thr = False, higher_ratio = False):
        # Definition of the input variables for the functions 
        self.height = 25        # This parameter will be overwritten during the iterations
        self.prominence = 15    # This parameter will be overwritten during the iterations
        self.width = 2
        self.ratio = 2.0
        self.window = 60
        self.doublepeak_elim = 1e-3 # (1 microssecond)
        self.wlen = 180        
        self.message = message
        self.static_thr = static_thr
        self.higher_ratio = higher_ratio
    
    def __enter__(self):
        # Code executed when entering the 'with' statement
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Code executed when exiting the 'with' statement
        self.cleanup()

    def cleanup(self):
        # Clean up any variables or resources here
        self.ratio_check = None
        self.peak_detector = None

    #Function for eliminating peaks that don't meet the ratio criteria
    def ratio_check(self, data):
        maxindex = 0
        if np.average(data) != 0:
            if np.abs(np.max(data))/np.abs(np.average(data)) > self.ratio: # Evaluate if the maximum value divided by the average is greater than 1.7
                maxindex = 1  # window 
        return maxindex
    
    def calculate_threshold(self, arr):
        max_hf = np.max(arr)
        threshold = (20/512)*max_hf + 10               #line equation
        return int(threshold)
    
    # write peaks to file
    def write_csv(self,  file_path, output_file, station, ts, x_i, x_f, npeaks, hfts, arr, OS, mac :str = ""):

        # write to CSV in write mode
        output_path = os.path.join(file_path, output_file)
        if OS == 0:
            with open(output_path, 'w', newline='') as csvfile:
                msvcrt.locking(csvfile.fileno(), msvcrt.LK_LOCK, 1)  # Locks the file     
                writer = csv.writer(csvfile) # Creates a writer object using the csv module
                writer.writerow([ts,station,x_i,x_f,mac]) # Writes a header row in the file
                for i in range(len(npeaks)): # For each of the peaks already filtered
                    writer.writerow([hfts[npeaks[i]], arr[npeaks[i]]]) # Writes the time and value of the peak
                # Release the lock
                msvcrt.locking(csvfile.fileno(), msvcrt.LK_UNLCK, 1)  # Unlocks the file
                
                return 0 # Success   
        elif OS == 1:
            with open(output_path, 'w', newline='') as csvfile:
                fcntl.flock(csvfile.fileno(), fcntl.LOCK_EX)                  
                writer = csv.writer(csvfile) # Creates a writer object using the csv module
                writer.writerow([ts,station,x_i,x_f,mac]) # Writes a header row in the file
                for i in range(len(npeaks)): # For each of the peaks already filtered
                    writer.writerow([hfts[npeaks[i]], arr[npeaks[i]]]) # Writes the time and value of the peak
                # Release the lock
                fcntl.flock(csvfile.fileno(), fcntl.LOCK_UN)
                return 0 # Success       

    #Function for peak detection
    def peak_detector(self, file_path, output_file, hfts, hf, station, ts, OS, mac: str = ""):
        try:
            x_i = hfts[0]  # Takes the initial time value
            x_f = hfts[-1] # Takes the final time value

            # Converts it to array format
            arr = np.array(hf)  
            xarr = np.array(hfts)

            if self.static_thr:
                self.height = int(self.static_thr)
            else:
                self.height = self.calculate_threshold(arr)  # Calculate the dynamic threshold for each file
            self.prominence = self.height*0.6

            if self.higher_ratio:
                self.ratio = float(self.higher_ratio)

            # Initialization of variables 
            npeaks = []  
            peaks, properties = find_peaks(arr, height=self.height, wlen=self.wlen, prominence=self.prominence, width=self.width)
            peaksN, propertiesN = find_peaks(-arr, height=self.height+18, wlen=self.wlen, prominence=self.prominence, width=self.width) # add 20 units for positive trigger

            initial_i = 0
            
            for i in range(len(peaks)):                                              
                    new = arr[peaks[i] - (self.window - 1):peaks[i] + 1] # Receives a window with 40 samples with the detected peak at sample 39
                    peakcheck = self.ratio_check(new) # Calls the ratio_check function
                    if peakcheck > 0:
                        if xarr[peaks[i]] - initial_i > self.doublepeak_elim: # Eliminates peaks that are actually just bit fluctuations
                            npeaks.append(peaks[i]) # Adds the detected peak to the array
                            initial_i = xarr[peaks[i]]  # Updates the value at each iteration
            initial_i = 0
            # For each of the peaks found by the 'find_peaks' function.
            for i in range(len(peaksN)): 
                    new = -arr[peaksN[i] - (self.window - 1):peaksN[i] + 1] # Receives a window with 40 samples with the detected peak at sample 39
                    peakcheck = self.ratio_check(new)  # Calls the ratio_check function
                    if peakcheck > 0:
                        if xarr[peaksN[i]] - initial_i > self.doublepeak_elim: # Eliminates peaks that are actually just bit fluctuations
                            npeaks.append(peaksN[i]) # Adds the detected peak to the array
                            initial_i = xarr[peaksN[i]] # Updates the value at each iteration

            if len(npeaks) > 500:
                self.height = 20
                self.prominence = self.height * 0.6
                npeaks = []

                while True:
                    peaks, properties = find_peaks(arr, height=self.height, wlen=self.wlen, prominence=self.prominence, width=self.width)
                    peaksN, propertiesN = find_peaks(-arr, height=self.height+18, wlen=self.wlen, prominence=self.prominence, width=self.width) # add 20 units for positive trigger

                    for i in range(len(peaks)):
                        new = arr[peaks[i] - (self.window - 1):peaks[i] + 1]
                        peakcheck = self.ratio_check(new)
                        if peakcheck > 0 and xarr[peaks[i]] - initial_i > self.doublepeak_elim:
                            npeaks.append(peaks[i])
                            initial_i = xarr[peaks[i]]

                    initial_i = 0

                    for i in range(len(peaksN)):
                        new = -arr[peaksN[i] - (self.window - 1):peaksN[i] + 1]
                        peakcheck = self.ratio_check(new)
                        if peakcheck > 0 and xarr[peaksN[i]] - initial_i > self.doublepeak_elim:
                            npeaks.append(peaksN[i])
                            initial_i = xarr[peaksN[i]]

                    # Stopping condition
                    if len(npeaks) <= 500:
                        break
                    else:
                        self.height += 10
                        self.prominence = self.height * 0.6

            # Opens the CSV file in write mode
            self.write_csv(file_path, output_file, station, ts, x_i, x_f, npeaks, hfts, arr, OS, mac)
            return len(npeaks)
                
        except Exception as e: # Exception handling.
            bins_print(self.message + f' Exception while performing peak detection on file {file_path}: {e}')
            return len(npeaks)