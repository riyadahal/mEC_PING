# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import glob, os
import numpy as np
from scipy import signal
from netpyne.analysis.tools import loadData

# =============================================================================
# 1. Configuration & Parameters
# =============================================================================

folder_EPSC = "../tau2runEPSC/" 
folder_IPSC = "../tau2runIPSC/" 
cell_id_EPSC = 'cell_1'  
cell_id_IPSC = 'cell_230' 

dt_ms = 0.01             
dt_s = dt_ms * 1e-3      
fs = 1.0 / dt_s          
sim_duration_ms = 6875.0 
fTheta = 8.0             
theta_period_ms = 1000.0 / fTheta # 125 ms

# Filter 50-300 Hz
b, a = signal.butter(3, [50, 300], fs=fs, btype='band')

output_dir = "1215_Cross_Corr/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# =============================================================================
# 2. Helper Functions
# =============================================================================

def load_first_pkl(folder_path):
    files = glob.glob(os.path.join(folder_path, "*.pkl"))
    if not files:
        raise FileNotFoundError(f"No .pkl file found in {folder_path}")
    filename = files[0]
    print(f"Loading: {filename}")
    data = loadData(filename)
    return data, os.path.basename(filename)

def get_filtered_trace(file_data, cell_id, b, a):
    """
    Extracts SEClamp current for a specific cell, handles data types safely,
    and bandpass filters it.
    """
    # 1. Navigate to the SEClamp data
    try:
        # We assume simData is a dict-like object.
        # Check if 'simData' exists
        if 'simData' not in file_data:
            raise KeyError("key 'simData' not found in file.")
        
        sim_data = file_data['simData']
        
        if 'SEClamp' not in sim_data:
            raise KeyError("key 'SEClamp' not found in simData.")
            
        se_clamp_data = sim_data['SEClamp']
        
        # 2. Check and Extract the specific cell data
        if cell_id not in se_clamp_data:
            # Debugging: Print available keys if the specific one is missing
            available_keys = list(se_clamp_data.keys())
            # Show first 5 keys to avoid cluttering console
            print(f"Key {cell_id} missing. Available keys (first 5): {available_keys[:5]}")
            raise KeyError(f"Cell ID {cell_id} not found in SEClamp data.")
        
        # 3. FIX: Convert to Python List first
        # NetPyNE sometimes returns a custom object that confuses np.array()
        # converting to list() strips away the wrapper.
        raw_obj = se_clamp_data[cell_id]
        raw_list = list(raw_obj)
        
        # Now convert to numpy array
        raw_trace = np.array(raw_list)
        
    except Exception as e:
        print(f"CRITICAL ERROR extracting {cell_id}: {e}")
        raise

    # 4. Processing
    # Scale to pA (assuming NetPyNE output is nA, scaling by 1000 -> pA)
    raw_trace_pA = 1000 * raw_trace 
    
    # Remove DC offset 
    raw_trace_pA -= np.mean(raw_trace_pA)
    
    # Apply Bandpass Filter
    filtered_trace = signal.filtfilt(b, a, raw_trace_pA)
    
    return filtered_trace

# =============================================================================
# 3. Main Processing Logic
# =============================================================================

print("--- Starting Cross-Correlation Analysis ---")

# A. Load Data
data_EPSC, name_EPSC = load_first_pkl(folder_EPSC)
data_IPSC, name_IPSC = load_first_pkl(folder_IPSC)

# B. Get Filtered Traces (50-300 Hz)
trace_EPSC = -1*get_filtered_trace(data_EPSC, cell_id_EPSC, b, a)
trace_IPSC = get_filtered_trace(data_IPSC, cell_id_IPSC, b, a)

# C. Slicing Logic (Last n Theta Cycles)
time_window_ms = 12 * theta_period_ms
start_time_ms = sim_duration_ms - time_window_ms

start_idx = int(start_time_ms / dt_ms)
end_idx = int(sim_duration_ms / dt_ms)

# Slice the filtered arrays
segment_EPSC = trace_EPSC[start_idx:end_idx]
segment_IPSC = trace_IPSC[start_idx:end_idx]

print(f"Analyzing time window: {start_time_ms} ms to {sim_duration_ms} ms")

# D. Cycle-by-Cycle Cross-Correlation
samples_per_cycle = int(theta_period_ms / dt_ms)
num_cycles = 12
corrs_list = []

print("Computing Cross-Correlation for individual cycles...")

for i in range(num_cycles):
    # 1. Extract the chunk for the i-th cycle
    idx_a = i * samples_per_cycle
    idx_b = (i + 1) * samples_per_cycle
    
    chunk_epsc = segment_EPSC[idx_a:idx_b]
    chunk_ipsc = segment_IPSC[idx_a:idx_b]
    
    # 2. Compute Cross-Correlation
    # Reference = EPSC, Target = IPSC
    # Positive lag means IPSC is shifted to the right (Later) relative to EPSC
    cc = signal.correlate(chunk_ipsc, chunk_epsc, mode='full')
    
    # 3. Normalization
    norm_factor = np.sqrt(np.sum(chunk_ipsc**2) * np.sum(chunk_epsc**2))
    cc_normalized = cc / norm_factor
    
    corrs_list.append(cc_normalized)

# E. Average
avg_corr = np.mean(np.array(corrs_list), axis=0)

# F. Lags
lags_samples = signal.correlation_lags(len(chunk_ipsc), len(chunk_epsc), mode='full')
lags_ms = lags_samples * dt_ms

# =============================================================================
# 4. Plotting (Restricted to +/- 50 ms)
# =============================================================================

window_limit_ms = 50
mask = (lags_ms >= -window_limit_ms) & (lags_ms <= window_limit_ms)

lags_window = lags_ms[mask]
corr_window = avg_corr[mask]

# Find Peak
peak_idx = np.argmax(corr_window)
peak_lag = lags_window[peak_idx]
peak_val = corr_window[peak_idx]

print(f"Peak Correlation: {peak_val:.4f} at Lag: {peak_lag:.2f} ms")

plt.figure(figsize=(8, 6))
plt.plot(lags_window, corr_window, 'k-', linewidth=2, label='Avg Cross-Corr (IPSC vs EPSC)')
plt.axvline(0, color='gray', linestyle='--', alpha=0.7)
plt.axhline(0, color='gray', linestyle='-', linewidth=0.5)

plt.plot(peak_lag, peak_val, 'ro', label=f'Peak: {peak_lag:.2f} ms')

plt.xlim(-window_limit_ms, window_limit_ms)
plt.title(f"Cross-Correlogram (EPSC -> IPSC)\nAvg of Last 4 Theta Cycles")
plt.xlabel("Lag (ms)\n(Positive Lag = IPSC follows EPSC)")
plt.ylabel("Correlation Coefficient")
plt.legend()
plt.grid(True, linestyle=':', alpha=0.6)

save_base = output_dir + "CrossCorr_EPSC_IPSC"
plt.savefig(save_base + ".png")
plt.savefig(save_base + ".eps")
print(f"Figures saved to {output_dir}")

plt.show()