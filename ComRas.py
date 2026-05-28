# This script loads NetPyne simulation files and generates a single
# combined raster plot for a fixed subset of 30 PV cells and 60 SC cells,
# focusing on the last two 125ms cycles (250ms total).

import numpy as np
import matplotlib.pyplot as plt
import glob, os
from netpyne.analysis.tools import loadData

# Raster plot specific parameters
theta_cycles_duration_ms = 125 # ms (Each cycle is 125ms)
marker_size_val = 5 # Marker size set to 5 for all plots

# Base directory for saving plots
base_output_dir = '1215_Raster/'

# Ensure the base output directory exists
if not os.path.exists(base_output_dir):
    os.makedirs(base_output_dir)

# --- FIXED RANDOM SUBSET SELECTION ---
# We use a fixed seed to ensure the same 30 PV and 60 SC cells are selected
# for every simulation file processed.
np.random.seed(42)

# PV cell IDs range from 1 to 100
pv_all_ids = np.arange(0, 99)
# SC cell IDs range from 101 to 500
sc_potential_ids = np.arange(100, 499) 

# Select fixed subsets (randomly chosen, but consistent due to fixed seed)
fixed_pv_subset = np.sort(np.random.choice(pv_all_ids, size=20, replace=False))
fixed_sc_subset = np.sort(np.random.choice(sc_potential_ids, size=80, replace=False))

# Combined array for mapping (PVs first, then SCs)
combined_target_ids = np.concatenate([fixed_pv_subset, fixed_sc_subset])
num_pv = len(fixed_pv_subset) # 30
num_sc = len(fixed_sc_subset) # 60
total_neurons = num_pv + num_sc # 90
# -------------------------------------


# --- MAIN FUNCTION TO PLOT COMBINED RASTER ---
def plot_combined_raster(pv_ids, sc_ids, all_spike_x, all_spike_ids, plot_title, base_filename, file_id_prefix):
    """
    Generates and saves a raster plot combining two distinct groups (PV, SC).
    PV cells are plotted first (Y-axis 1-30, blue), followed by SC cells (Y-axis 31-90, red).
    """
    
    # 1. Define the overall mapping from actual NetPyne ID to Serial Plot ID (1-based)
    combined_ids = np.concatenate([pv_ids, sc_ids])
    master_serial_map = {cell_id: i + 1 for i, cell_id in enumerate(combined_ids)}
    
    # 2. Filter spikes that belong to the combined subset (PV or SC)
    combined_mask = np.isin(all_spike_ids, combined_ids)
    
    group_x_values = all_spike_x[combined_mask]
    group_cell_ids = all_spike_ids[combined_mask]
    
    if len(group_x_values) == 0:
        print(f"  No spikes for Combined Raster in the last two cycles. Skipping plot.")
        return

    # 3. Get plotting Y-positions and colors
    plot_y_positions = np.array([master_serial_map[cid] for cid in group_cell_ids])
    
    is_pv_spike = np.isin(group_cell_ids, pv_ids)
    is_sc_spike = np.isin(group_cell_ids, sc_ids)
    
    # 4. Plotting
    fig = plt.figure(figsize=(7, 7)) # Increased height for better visibility
    ax = plt.gca()
    
    # Plot PV spikes (Blue)
    plt.plot(group_x_values[is_pv_spike], plot_y_positions[is_pv_spike], '|',
             color='blue', markersize=marker_size_val, markeredgewidth=1.5, label='PV Cells')

    # Plot SC spikes (Red)
    plt.plot(group_x_values[is_sc_spike], plot_y_positions[is_sc_spike], '|',
             color='red', markersize=marker_size_val, markeredgewidth=1.5, label='SC Cells')

    # Add a dividing line between PV and SC groups
    ax.axhline(y=len(pv_ids) + 0.5, color='gray', linestyle='--', linewidth=1)
    
    # 5. Y-axis setup
    # Define labels (only mark the boundary between PV and SC groups)
    y_tick_positions = np.arange(1, total_neurons + 1)
    
    y_tick_labels = [str(pos) if pos % 10 == 0 or pos == 1 or pos == total_neurons else '' for pos in y_tick_positions]
    
    # Add labels for the group boundaries for clarity
    y_tick_labels[0] = '1 (PV)' # First PV
    y_tick_labels[num_pv - 1] = str(num_pv) # Last PV
    y_tick_labels[num_pv] = str(num_pv + 1) + ' (SC)' # First SC
    y_tick_labels[total_neurons - 1] = str(total_neurons) # Last SC

    plt.yticks(y_tick_positions, y_tick_labels, fontsize=8)
    plt.ylim([0.5, total_neurons + 0.5])
    
    # Remove the top and right spines ("box")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set X-axis to relative time (0-250ms for 4 cycles)
    plt.xlabel('Time (ms)', fontsize=12)
    plt.xlim([0, theta_cycles_duration_ms * 4]) 

    plt.ylabel(f'Neuron Serial Number ({num_pv} PV, {num_sc} SC)', fontsize=6)
    #plt.title(plot_title, fontsize=14)
    #plt.legend(loc='upper right', frameon=False)
    plt.tight_layout()

    # Construct the full file path for saving
    full_path_png = os.path.join(file_id_prefix, base_filename + '.png')
    full_path_eps = os.path.join(file_id_prefix, base_filename + '.eps')

    plt.savefig(full_path_png, dpi=300)
    plt.savefig(full_path_eps, format='eps')
    
    plt.close(fig)


# --- Main loop to process files ---
file_pattern = "../tau2runIPSC/*.pkl"
simulation_files = glob.glob(file_pattern)

if not simulation_files:
    print(f"No simulation files found matching pattern: {file_pattern}")
    print("Please check the path and filename pattern.")
else:
    print(f"Found {len(simulation_files)} simulation files.")

for file in simulation_files:
    print(f"\nProcessing file: {file}")
    fileInfo = loadData(file)
    
    # Extract fileName2 as the last 3 characters of the simConfig filename
    full_sim_filename = fileInfo['simConfig']['filename']
    # Safely get the last 3 characters
    fileName2 = full_sim_filename[-4:] if len(full_sim_filename) >= 3 else full_sim_filename
    print(f"  Extracted identifier: {fileName2}")

    sim_duration = fileInfo['simConfig']['duration'] # Simulation duration in ms

    # Get all spike times and cell IDs from the simulation
    spkt_all = np.array(fileInfo['simData']['spkt'])
    spkid_all = np.array(fileInfo['simData']['spkid'])

    # --- Identify and extract data for the LAST two theta cycles (Time 0-250ms relative) ---
    total_num_cycles_in_sim = int(np.floor(sim_duration / theta_cycles_duration_ms))
    num_cycles_to_plot = 4 

    if total_num_cycles_in_sim < num_cycles_to_plot:
        print(f"  Skipping {file}: Simulation duration only has {total_num_cycles_in_sim} full cycles. Need {num_cycles_to_plot}.")
        continue

    # 1. Define the index and time boundaries for the last two cycles
    start_cycle_idx = total_num_cycles_in_sim - num_cycles_to_plot # 0-based index of the start cycle
    
    start_time_ms = start_cycle_idx * theta_cycles_duration_ms
    end_time_ms = total_num_cycles_in_sim * theta_cycles_duration_ms
    
    print(f"  Plotting for the LAST {num_cycles_to_plot} cycles ({theta_cycles_duration_ms * num_cycles_to_plot:.0f}ms): {start_time_ms:.0f}ms - {end_time_ms:.0f}ms)")

    # 2. Filter all spikes for events within this specific last cycle
    spikes_in_mask = (spkt_all >= start_time_ms) & (spkt_all < end_time_ms)
    spkts_last_cycle = spkt_all[spikes_in_mask]
    spkids_last_cycle = spkid_all[spikes_in_mask]

    # 3. Convert spike times to be relative to the start of the two-cycle window (X-axis values 0 to 250ms)
    relative_spikes_ms = spkts_last_cycle - start_time_ms
    
    # The X-axis data is now relative_spikes_ms

    # --- Generate and save the COMBINED RASTER plot ---
    plot_combined_raster(fixed_pv_subset, fixed_sc_subset, relative_spikes_ms, spkids_last_cycle,
                         f'Combined Raster', 
                         f'Combined_PV{num_pv}_SC{num_sc}_Raster_{fileName2}', base_output_dir)

print("\nAll combined raster plot generation complete across all specified files.")

