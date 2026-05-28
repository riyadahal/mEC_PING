# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from funcs import *
import glob, os
from netpyne.analysis.tools import loadData
import matplotlib.gridspec as gridspec
import numpy as np
from scipy import signal
# Import function for Morlet Wavelets
from neurodsp.timefrequency.wavelets import compute_wavelet_transform

# Function to plot the ACF and mark the peak (Updated to find FIRST local maximum)
def plot_autocorrelation(signal_avg, dt_s, min_lag_samples, fileName2, FileID):
    
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
    
    # --- NEW LOGIC: Find the FIRST local maximum ---
    # Find all local maxima in the search region. 
    # height=0.05 ensures the peak correlation is meaningfully positive.
    # distance=20 prevents picking up noise between major peaks (0.2 ms separation).
    peak_indices_relative_to_search, _ = signal.find_peaks(
        search_lags, 
        height=0.05, 
        distance=20 
    )
    # --- END NEW LOGIC ---
    
    # If a local peak is found:
    if len(peak_indices_relative_to_search) > 0:
        # Select the index of the FIRST local maximum
        peak_index_relative = peak_indices_relative_to_search[0]
        peak_index_absolute = peak_index_relative + min_lag_samples
        
        # Determine Period and Frequency
        period_time_ms = full_lags_time_ms[center_index + 1 + peak_index_absolute]
        frequency_Hz = 1000.0 / period_time_ms
        
        # Plotting
        plt.figure(figsize=(8, 5))
        plt.plot(full_lags_time_ms, autocorr, label='Autocorrelation Function', color='gray') 
        
        """# Mark the minimum lag threshold (only shown on positive side)
        min_lag_ms = full_lags_time_ms[center_index + 1 + min_lag_samples]
        plt.axvline(x=min_lag_ms, color='r', linestyle='--', 
                    label=f'Min Lag Threshold ({min_lag_ms:.2f} ms)')"""
        
        # Mark the detected peak
        peak_correlation_value = positive_lags_acorr[peak_index_absolute]
        plt.plot(period_time_ms, peak_correlation_value, 'ro', 
                 label=f'First Peak: {period_time_ms:.2f} ms ({frequency_Hz:.2f} Hz)')
        
        max_lag_ms = full_lags_time_ms[-1]
        plt.xlim(-max_lag_ms, max_lag_ms) 
        plt.ylim(-1.05, 1.05)
        plt.title(f'Autocorrelation of Averaged Filtered IPSC (First Local Max: {frequency_Hz:.2f} Hz)')
        plt.xlabel('Lag Time (ms)')
        plt.ylabel('Correlation Coefficient (Normalized)')
        #plt.legend(loc='upper right')
        #plt.grid(True, linestyle=':', alpha=0.6)
        
        plt.savefig(FileID + fileName + '_ACF.png')
        plt.savefig(FileID + fileName + '_ACF.eps')
        plt.close()
        
        return frequency_Hz
        
    return None

#########################################################################################################################################
# Some Housekeeping parameters

# --- Time Step Definitions ---
dt_ms = 1e-2        # Time step in milliseconds (0.01 ms)
dt_s = dt_ms * 1e-3 # Time step in seconds (1e-5 s)
# -----------------------------

fTheta = 8
thetaPeriod = 1000/fTheta # Period in ms

VC = 0 #Clamped voltage
rs=1e-4
Factor = 1e3

CapCurrentLeft = 70; CapCurrentRight = 80;
b, a = signal.butter(3, [50, 300], fs=1/dt_s, btype='band') 
waveletPlotPoints = 0

FileID = "1215_Compare/" 
if not os.path.exists(FileID):
    os.makedirs(FileID)

dt_rec = dt_ms; # Down sample the rates (in ms)
chunkLength = int(np.round(thetaPeriod/dt_rec));
print('chunkLength is ', chunkLength)
CycToPlot = 13
offset = 40
fs = 1.0 / dt_s # Sampling frequency in Hz (100,000 Hz)

freqs = np.arange(50., 301, 1)

# Define the minimum lag in samples to ignore high-frequency noise (350 samples > 3.33 ms, the period of 300 Hz)
MIN_LAG_SAMPLES = 350 

for file in glob.glob("../tau2runIPSC/*.pkl"):
    print('--------------------------------------------------')
    print('Loading the file from specified folder')
    fileInfo = loadData(file)
    
    # --- FILENAME FIX ---
    #base_name = os.path.basename(file)
    #current_fileName = os.path.splitext(base_name)[0]
    #fileName2 = current_fileName # Use the full, current filename for saving
    
    # Keep old variable for logging the original name (for comparison)
    fileName = fileInfo['simConfig']['filename'][-5:]
    # --------------------
    
    print(f"Analyzing files (Original Name): {fileName}")
    #print(f"Analyzing files (HD Name): {fileName2}")
    
    sim_time = fileInfo['simConfig']['duration']
    cellTrace = {}
    cellTrace['0'] = fileInfo['simData']['SEClamp']['cell_1']
    
    print(sim_time)
    cellTr0 = 1000*np.array(cellTrace['0'])

    y0 = cellTr0#(VC-cellTr0)/rs*Factor
    y0 -= y0[-1] #offsetting based on zero
    try:
        index = int(np.argwhere(y0<minV)[0])
    except:
        index = CapCurrentLeft
        
    rates = signal.filtfilt(b, a, y0) #Band pass filter (gamma waveform)
    rates -= rates[-1]
    print('Printing length of y0 and rates ',len(y0),len(rates))
    print('Printing length of y0 and rates after offset cycles',len(y0[(offset * chunkLength):]),len(rates[offset * chunkLength:]))
    # Initialize lists to store results
    mwt_avg_list = []
    signal_chunks = []
    
    # --- ACF and Wavelet Calculation (over 8 cycles: j=2 to j=9) ---
    for j in range(offset+2, CycToPlot+offset + 2):
        signal_1 = rates[j * chunkLength: (j+1) * chunkLength]
        signal_chunks.append(signal_1) # Store chunks for averaging later

        # Wavelet Calculation
        mwt = compute_wavelet_transform(signal_1, fs=fs, freqs=freqs, n_cycles = 6)
        mwt_avg_list.append(mwt)

    # ------------------ Wavelet Analysis Output ----------------------
    factor = 2
    mwt_avg = np.array(mwt_avg_list)
    mwt_avg = np.mean(mwt_avg,axis = 0) 
    z = np.abs(factor*mwt_avg)**2
    
    maxPowerAvg = np.max(np.max(z))
    
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
    acf_frequency = plot_autocorrelation(signal_avg_8cycles, dt_s, MIN_LAG_SAMPLES, fileName, FileID)

    print('\nAutocorrelation Analysis:')
    if acf_frequency:
        print(f'ACF Frequency (from averaged 8 cycles - First Local Max): {acf_frequency:.2f} Hz')
    else:
        print('ACF analysis failed to detect a peak.')
    # -----------------------------------------------------------------

    # --- Plotting Code Below (Contourf and Currents) ---
    plt.figure()
    #maxPowerAvg = 3500
    plotPow(z,Ncycs=1,dt=dt_rec,fs=freqs,levels=np.linspace(0.,maxPowerAvg,150,endpoint=True))
    plt.savefig(FileID+fileName+'_cell_230_Pow.eps')
    plt.savefig(FileID+fileName+'_cell_230_Pow.png')
    plt.clf();
    
    fig = plt.figure(figsize=(10, 8))
    outer = gridspec.GridSpec(1, 1, wspace=0.3, hspace=0.1)
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0], wspace=0.1, hspace=0.1, height_ratios=[1,1])

    ax0 = plt.Subplot(fig, inner[0])
    tSim = np.linspace(0,int(sim_time)+1,int(sim_time)*100)
    print('Length of the tSim ',len(tSim))
    ax0.plot(tSim[offset*12500:], y0[offset*12500:], 'r-',alpha = 1, label='Unfiltered')
    ax0.set_xlim(sim_time-2.*125.+30,sim_time-2.*125.+60);
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.spines['left'].set_visible(True)
    ax0.set_ylabel("Current (pA)")
    #ax0.set_ylim([-100,1200])
    ax0.set_xticks([])
    ax0.set_xticklabels([])
    ax0.legend(loc='upper right')
    fig.add_subplot(ax0)

    """ax1 = plt.Subplot(fig, inner[1])
    ax1.plot(tSim, rates, 'b-',label='Filtered 50-300 hz')
    ax1.set_xlim(sim_time-2.*125.,sim_time);
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['left'].set_visible(True)
    ax1.set_ylabel("Current (pA)")
    ax1.set_ylim([-350,350])
    ax1.set_xticks([])
    ax1.set_xticklabels([])
    ax1.legend(loc='upper right')
    fig.add_subplot(ax1)"""

    gsin = 5
    ax2 = plt.Subplot(fig, inner[1])
    t1 = np.linspace(sim_time-4.*125.,sim_time,int(4*125/dt_ms)+1)
    ax2.set_xlim(sim_time-2.*125.,sim_time)
    ax2.plot(t1, gsin*np.sin(2*np.pi*fTheta*t1*1e-3-np.pi/2),color = 'k')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.set_xticks([])
    ax2.set_yticks([])
    fig.add_subplot(ax2)

    fig.tight_layout()
    fig.savefig(FileID+fileName+'_cell_230_Currents.png', bbox_inches='tight')
    fig.savefig(FileID+fileName+'_cell_230_Currents.eps', bbox_inches='tight')
    plt.close(fig)

