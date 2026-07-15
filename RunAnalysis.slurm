#!/bin/bash
#SBATCH --job-name=GEGI_Analysis    # Job name
#SBATCH --array=0-9   #Make it 0-0 for single sweeps
#SBATCH --nodes=1                        # Run all processes on a single node.
#SBATCH --cpus-per-task=1                # Number of CPU cores per task (multithreaded tasks)
#SBATCH --mem=32G
#SBATCH  --output=/mnt/beegfs/home/avedur/NetPyNe_Sims/mEC_model2026/Analysis/output/Simuslurm-%A_%a.out

# Forcefully remove any display routing
unset DISPLAY
export QT_QPA_PLATFORM=offscreen
export XDG_RUNTIME_DIR=/tmp

python3 -u Double_Sweep.py 0

#Comment the above line and uncomment the below line for running single sweep
#python3 -u Single_Sweep.py
