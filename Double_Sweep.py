# -*- coding: utf-8 -*-
#This code makes double parameter sweep

import matplotlib
matplotlib.use('Agg') # This MUST come before any other imports!
from funcs import *
import glob, os
from netpyne.analysis.tools import loadData
from itertools import compress
#import matplotlib.gridspec as gridspec
import numpy as np
from scipy import signal
#import matplotlib.pyplot as plt
from funcs import *
# Import simulation and plot code to create & visualize data
#from neurodsp.plts import plot_timefrequency
import re

# Import function for Morlet Wavelets
from neurodsp.timefrequency.wavelets import compute_wavelet_transform
import sys
seed_num = sys.argv[1]
print('Seed_Num = ',seed_num)

#########################################################################################################################################

# Function to plot the ACF and mark the peak (Updated to find FIRST local maximum)
def plot_autocorrelation(signal_avg, dt_s, min_lag_samples):
    
    # 1. Compute ACF on the averaged signal
    demeaned_signal = signal_avg - np.mean(signal_avg)
    autocorr = np.correlate(demeaned_signal, demeaned_signal, mode='full')
    
    # Normalize the ACF (to 1.0 at Lag 0)
    autocorr = autocorr / autocorr.max()
    
    # Setup Lag Time Vector
    N = len(demeaned_signal)
    center_index = N - 1
    
    # Create a full lag vector from -(N-1) to (N-1) samples
    full_lags_samples = np.arange(-(N - 1), N) 
    # Convert to time in ms
    full_lags_time_ms = full_lags_samples * dt_s * 1e3
    
    # Positive lags only for finding the peak
    positive_lags_acorr = autocorr[center_index + 1:] 

    # Define the search region starting after the minimum lag threshold
    search_lags = positive_lags_acorr[min_lag_samples:]
    
    # --- FIND LOCAL MAXIMA ---
    # Find all local maxima in the search region. 
    # height=0.05 ensures the peak correlation is meaningfully positive.
    # distance=20 prevents picking up noise between major peaks (0.2 ms separation).
    peak_indices_relative_to_search, _ = signal.find_peaks(
        search_lags, 
        height=0.05, 
        distance=20 
    )
    # -------------------------
    
    # If a local peak is found:
    if len(peak_indices_relative_to_search) > 0:
        
        # --- MODIFIED LOGIC: Find the index corresponding to the peak with MAX correlation ---
        # 1. Get the correlation values (heights) for all detected peaks
        peak_correlations = search_lags[peak_indices_relative_to_search]
        
        # 2. Find the index within the 'peak_indices_relative_to_search' list that corresponds to the MAX correlation value
        max_peak_value_index = np.argmax(peak_correlations)
        
        # 3. Select the actual index (lag position) of the strongest peak
        peak_index_relative = peak_indices_relative_to_search[max_peak_value_index]
        # --- END MODIFIED LOGIC ---
        
        peak_index_absolute = peak_index_relative + min_lag_samples
        
        # Determine Period and Frequency
        period_time_ms = full_lags_time_ms[center_index + 1 + peak_index_absolute]
        frequency_Hz = 1000.0 / period_time_ms
        
        return frequency_Hz
        
    return None

#########################################################################################################################################

def get_conductances(fileName):
    """
    Extracts the E2I (i) and I2E (j) indices from the simulation filename 
    and calculates the final scaled conductance values (both in nS).

    Filename Pattern: '...Conductances_i_j_data.pkl'
    i (E2I Index) and j (I2E Index) both run from 0 to 10.

    Args:
        fileName (str): The filename string (e.g., 'Conductances_1_8_data.pkl' or 'Last_Batch_Conductances_1_8_data.pkl').
        weight_list (np.array): The base array for E2I.
        weight_i2e_array (np.array): The final scaled array for I2E.

    Returns:
        tuple: (E2I_conductance_nS, I2E_conductance_nS)
    """
    #weight_list = np.arange(0.00005,0.00105,0.00005)
    #weights = weight_list
    #scale_factors = np.arange(0.5, 10.0, 0.5)
    #base_weight_value = 1.65
    #weight_i2e_array = base_weight_value * scale_factors
    # 1. Extract the i and j indices using regex, anchored to 'Conductances_'
    match = re.search(r'GE_GI_Sims_\d+_(\d+)_(\d+)$', fileName)

    if match:
        # Group 1 is the E2I Index (i) -> 0 to 10
        e2i_index_i = int(match.group(1))
        print(e2i_index_i)
        
        # Group 2 is the I2E Index (j) -> 0 to 10
        i2e_index_j = int(match.group(2)) 
        print(i2e_index_j)

        """# --- 2. Calculate E2I Conductance (nS) ---
        if 0 <= e2i_index_i < len(weight_list):
            base_e2i = weight_list[e2i_index_i]
            # Scaling: Base value * 120,000 (to nS)
            E2I_conductance_nS = base_e2i * 100000
        else:
            print(f"Error: E2I index {e2i_index_i} out of bounds.")
            E2I_conductance_nS = np.nan

        # --- 3. Calculate I2E Conductance (nS) ---
        if 0 <= i2e_index_j < len(weight_i2e_array):
            base_i2e = weight_i2e_array[i2e_index_j]
            # Scaling: Base value (which is in uS) * 20 (to nS)
            I2E_conductance_nS = base_i2e * 20
        else:
            print(f"Error: I2E index {i2e_index_j} out of bounds.")
            I2E_conductance_nS = np.nan"""

        return e2i_index_i, i2e_index_j

    else:
        print(f"Error: Filename pattern not matched for {fileName}.")
        return np.nan, np.nan

#########################################################################################################################################
#Some Housekeeping parameters
dt = 1e-2 # in ms
dt_s = dt * 1e-3 # Time step in seconds (1e-5 s)
fTheta = 8
thetaPeriod = 1000/fTheta

#Firstly load the files

cellTrace={}

CapCurrentLeft = 70; CapCurrentRight = 80;
b, a = signal.butter(3, [50, 300], fs=1/(dt*1e-3), btype='band')
waveletPlotPoints = 0

FileID = "0711_GEGI_Batch_Data/" 
os.makedirs(FileID, exist_ok = True)

dt_rec = dt; # Down sample the rates (in ms)

#freqss = np.arange(50.,300.1,0.2);

chunkLength = int(np.round(thetaPeriod/dt_rec));

#print('Chunklength value and shape is ',chunkLength, np.shape(chunkLength))

CycToPlot = 13
offset = 40

# Generate a test sinusoidal signal
fs = 1/(dt*1e-3) # Sampling frequency in Hz
t = np.linspace(0, 125, chunkLength) # Time vector from 0 to 1 second

freqs = np.arange(50., 301, 1)
Power_Array = []
i = 0

# Define the minimum lag in samples to ignore high-frequency noise (350 samples > 3.33 ms, the period of 300 Hz)
MIN_LAG_SAMPLES = 350 
#{seed_num}
for file in glob.glob(f"../0709_GE_GI_Sims/*Sims_{seed_num}_*_*_data.pkl"):
    print('Loading the file from specifie folder')
    fileInfo = loadData(file)
    fileName = fileInfo['simConfig']['filename'][16:]
    print(fileName)
    sim_time = fileInfo['simConfig']['duration']
    cellTrace['0'] = fileInfo['simData']['SEClamp']['cell_130']#20
    tSim = fileInfo['simData']['t']
    E__I,I__E = get_conductances(fileName)
    print(E__I,I__E)
    
    cellTr0 = 1000*np.array(cellTrace['0'])

    y0 = cellTr0#(VC-cellTr0)/rs*Factor
    y0 -= y0[-1] #offsetting based on zero
    try:
        index = int(np.argwhere(y0<minV)[0])
    except:
        index = CapCurrentLeft
    print('Printing legths of cellTr, tSim, and y0 : ',len(cellTr0),len(tSim),len(y0))
    rates = signal.filtfilt(b, a, y0) #Here is where we are using the band pass filter and using the signal to get gamma waveform
    rates -= rates[-1]

    print('Printing length of y0 and rates ',len(y0),len(rates))
    print('Printing length of y0 and rates after offset cycles',len(y0[(offset * chunkLength):]),len(rates[offset * chunkLength:]))

    # n_cycles is the center frequency of the wavelet - w0
    # It's also the length of the filter, as the number of cycles of the oscillation with specified frequency
    # the analytic morlet wavelet in matlab uses w0 = 6

    # Compute wavelet transform using compute Morlet wavelet transform algorithm
    mwt_avg = []
    signal_chunks = []
    for j in range(offset+2,offset+CycToPlot+2):
        signal_1 = rates[(j)*chunkLength: (j+1)*chunkLength]
        signal_chunks.append(signal_1) # Store chunks for averaging later
        mwt = compute_wavelet_transform(signal_1, fs=fs, freqs=freqs, n_cycles = 6)
        mwt_avg.append(mwt)

    factor = 2
    mwt_avg = np.array(mwt_avg)
    #print('Shape of mwt_avg is', np.shape(mwt_avg))

    mwt_avg = np.mean(mwt_avg,axis = 0)
    z = np.abs(factor*mwt_avg)**2
    
    maxPowerAvg = np.max(np.max(z))

    # Find the frequency of peak power in each band
    max_power_index = np.argmax(np.max(z, axis=1))

    freq_at_peak_power = freqs[max_power_index]
    print('Wavelet Analysis (Average Power):')
    print(f'Maximum power of this signal is {maxPowerAvg:.2f}')
    print(f'Frequency at maximum power is {freq_at_peak_power:.2f} Hz')

    # ------------------ ACF Analysis & Plotting --------------------------
    
    # Time domain average of the 8 cycles
    min_len = min(len(s) for s in signal_chunks)
    signals_array = np.array([s[:min_len] for s in signal_chunks])
    signal_avg_8cycles = np.mean(signals_array, axis=0)
    
    # Calculate and plot the ACF on the averaged signal
    acf_frequency = plot_autocorrelation(signal_avg_8cycles, dt_s, MIN_LAG_SAMPLES)

    print('\nAutocorrelation Analysis:')
    if acf_frequency:
        print(f'ACF Frequency (from averaged 8 cycles - First Local Max): {acf_frequency:.2f} Hz')
    else:
        print('ACF analysis failed to detect a peak.')
    
    seed = int(seed_num)
    Power_Array.append([9+seed,E__I,I__E,maxPowerAvg,freq_at_peak_power,acf_frequency])
    print(Power_Array)

power_array = np.array(Power_Array)
#print(np.shape(power_array))
#power_array = power_array[power_array[:, 0].argsort()]
sort_keys = (power_array[:, 2],power_array[:, 1], power_array[:, 0]) 
power_array = power_array[np.lexsort(sort_keys)]
#power_array = np.round(power_array,2)
print(power_array)
np.savetxt(FileID+'0711_GEGI_Data_Seed_9.txt', power_array, fmt="%.2f", delimiter=",")
np.savetxt(FileID+'0711_GEGI_Data_Seed_9.csv',power_array, fmt="%.2f", delimiter=",")
