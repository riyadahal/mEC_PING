"""
Sairam

batch.py 

Batch simulation for mEC model using NetPyNE

Contributors: @Ananth
"""
import os
from netpyne.batch import Batch
from netpyne import specs
import numpy as np

# ----------------------------------------------------------------------------------------------
# Parameter for Changing Weights of Synaptic Connections
# ----------------------------------------------------------------------------------------------

def Batch_EI():
    params = specs.ODict()

    # Define a list of seed dictionaries
    num_simulations = 1
    seed_list = []
    for i in range(num_simulations):
        seed_list.append({'conn': 4321 + i, 'stim': 1234 + i, 'loc': 4321 + i, 'Inet': 7894 + i , 'seed_index': i}) # Example of varying seeds

    params['seeds'] = seed_list
    #params['g_sinExc'] = [27.*1e-3,28.*1e-3,29.*1e-3,30.*1e-3]
    weight_list = np.linspace(0.0002,0.0050,25)
    weights = weight_list
    params['Weight_E2I'] = weights
    groupedParams = []
    initCfg = {}
    #b = Batch(params=params, initCfg=initCfg, groupedParams=groupedParams)
    b = Batch(params=params,
              netParamsFile='src/netParams.py',
              cfgFile='src/cfg.py')
    b.method = 'grid'
    return b
#'conn': 4321, 'stim': 1234, 'loc': 4321, 'cell': 4321, 'brian2': 7894

def EI_Strength():
    params = specs.ODict()
    num_simulations = 10
    seed_list = []
    for i in range(num_simulations):
        seed_list.append({'conn': 4321 + i, 'stim': 1234 + i, 'loc': 4321 + i, 'cell': 7894 + i, 'brian2': 7894 + i}) # Example of varying seeds
    params['seeds'] = seed_list
    weight_list = np.arange(0.00005,0.00105,0.00005) #np.linspace(0.0002,0.0050,25)
    weights = weight_list
    params['Weight_E2I'] = weights
    #params['g_sin'] = np.linspace(0.*1e-3,9.*1e-3,10)
    #params['g_sin_Exc'] = np.linspace(0.*1e-3,9.*1e-3,10)

    groupedParams = []
    initCfg = {}
    b = Batch(params=params,
              netParamsFile='src/netParams.py',
              cfgFile='src/cfg.py')
    b.method = 'grid'
    return b

def GE_GI():
    params = specs.ODict()
    num_simulations = 9
    seed_list = []
    for i in range(num_simulations):
    #i = 15
    	seed_list.append({'conn': 4321 + i, 'stim': 1234 + i, 'loc': 4321 + i, 'cell': 7894 + i, 'brian2': 7894 + i}) # Example of varying seeds

    params['seeds'] = seed_list
    
    weight_list = np.arange(0.00005,0.00105,0.00005)
    weights = weight_list
    params['Weight_E2I'] = weights
    #scale_factors = np.arange(0.0, 5.5, 0.5)
    params['Weight_I2E'] = ['0.5*lognormal(1.65,2.17)*1e-3','1.0*lognormal(1.65,2.17)*1e-3','1.5*lognormal(1.65,2.17)*1e-3','2.0*lognormal(1.65,2.17)*1e-3','2.5*lognormal(1.65,2.17)*1e-3','3.0*lognormal(1.65,2.17)*1e-3','3.5*lognormal(1.65,2.17)*1e-3','4.0*lognormal(1.65,2.17)*1e-3','4.5*lognormal(1.65,2.17)*1e-3','5.0*lognormal(1.65,2.17)*1e-3','5.5*lognormal(1.65,2.17)*1e-3','6.0*lognormal(1.65,2.17)*1e-3','6.5*lognormal(1.65,2.17)*1e-3','7.0*lognormal(1.65,2.17)*1e-3','7.5*lognormal(1.65,2.17)*1e-3','8.0*lognormal(1.65,2.17)*1e-3','8.5*lognormal(1.65,2.17)*1e-3','9.0*lognormal(1.65,2.17)*1e-3','9.5*lognormal(1.65,2.17)*1e-3','10.0*lognormal(1.65,2.17)*1e-3']
    groupedParams = []
    initCfg = {}
    b = Batch(params=params,
              netParamsFile='src/netParams.py',
              cfgFile='src/cfg.py')
    b.method = 'grid'
    return b

# ----------------------------------------------------------------------------------------------
# Run configurations
# ----------------------------------------------------------------------------------------------
def setRunCfg(b, type='mpi_bulletin', nodes=1, coresPerNode=8):
    if type=='mpi_bulletin':
        b.runCfg = {'type': 'mpi_bulletin', 
            'script': 'src/init.py', 
            'skip': True}

    elif type=='mpi_direct':
        b.runCfg = {'type': 'mpi_direct',
            'cores': 4,
            'script': 'init_cell.py',
            'mpiCommand': 'mpirun',
            'skip': True}

    elif type=='tigerfish':
        b.runCfg = {'type': 'hpc_slurm', 
            'nodes': nodes,
            'coresPerNode': coresPerNode,   
            'script': 'src/init.py', 
            'skip': True,
            'mpiCommand': 'mpiexec', 
            'vmem': str(coresPerNode * 2) + 'G',  # Dynamically requests 2GB per core 
            'walltime': "00:30:00",               # 30 mins (gives a safe buffer for your 6 min sim)
            'custom': '''
 source /mnt/beegfs/home/avedur/miniconda3/etc/profile.d/conda.sh
 conda activate Netpyne
 export PYTHONPATH="."
 export NEURON_MODULE_OPTIONS="-nogui"
 ''',                                              # Injects our fixes into every batch script!
            'skipCustom': '_raster.png'}

# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------
"""
if __name__ == '__main__': 

    b = Vary_Weight_Seeds()
    b.batchLabel = 'Varying_Seeds'
    path = '~/output/'
    os.makedirs(path, exist_ok = True) 
    b.saveFolder = path+b.batchLabel
    b.method = 'grid'  # evol
    setRunCfg(b, 'tigerfish', nodes=1, coresPerNode=1)  # cores = nodes * 8 
    b.run() # run batch """

if __name__ == '__main__': 

    b = GE_GI()
    b.batchLabel = '0710_GE_GI_Sims'
    folder = '0710_GE_GI_Sims'
    os.makedirs(folder, exist_ok = True) 
    b.saveFolder = folder
    b.method = 'grid'  # evol
    setRunCfg(b, 'tigerfish', nodes=1, coresPerNode=20)  # cores = nodes * 8 
    b.run() # run batch 
