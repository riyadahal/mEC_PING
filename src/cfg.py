from netpyne import specs
import numpy as np

# Simulation options
cfg = specs.SimConfig()       # object of class SimConfig to store simulation configuration

###############################################################################
## Simulation parameters
###############################################################################
cfg.ThetaCycles = 15        # Number of theta cycles to simulate
cfg.Theta2Plot = 4           # Number of theta cycles to plot
cfg.duration = cfg.ThetaCycles*125+5000.  # Duration of the simulation, in ms + 5s baseline needed for SC cells
cfg.dt = 1e-2                # Internal integration timestep to use
cfg.hParams = {'v_init': -65}  
cfg.saveFolder = 'output/today'  # Folder to save output
cfg.simLabel = 'PING'  # Simulation label, used in output file names
cfg.validateNetParams = False
cfg.verbose = False           # Shows detailed messages
cfg.progressBar = 0       # Shows progress bar
cfg.recordStep = cfg.dt        # Step size in ms to save data (e.g. V traces, LFP, etc)
cfg.seeds = {'conn': 4321, 'stim': 1234, 'loc': 4321, 'cell': 4321, 'brian2': 7894, 'opto': 42} # Random seeds for reproducibility. brian2 seed is for the PV network.
cfg.saveDataInclude = ['simData', 'simConfig', 'net']  # Which data to save in the output file
cfg.printRunTime = 0.1 # Print run time every 0.1 seconds
cfg.recordTime = False  
cfg.createNEURONObj = True
cfg.createPyStruct = True
cfg.backupCfgFile = ['src/cfg.py', cfg.saveFolder]
cfg.gatherOnlySimData = False
cfg.saveCellSecs = False
cfg.saveCellConns = True
cfg.saveJson = False
cfg.savePickle = True
cfg.recordStims = True

#------------------------------------------------------------------------------
# Current inputs 
#------------------------------------------------------------------------------
cfg.addIClamp = False

cfg.IClamp1 = {'pop': 'FS', 'sec': 'soma', 'loc': 0.5, 'start': 50, 'dur': 100, 'amp': 0.50}
cfg.IClamp2 = {'pop': 'SC', 'sec': 'soma', 'loc': 0.5, 'start': 50, 'dur': 100, 'amp': 0.50}

###############################################################################
## SimParams
############################################################################### 
# PV+ cells properties 
cfg.NPV=100 # Number of Inhibitory neurons
cfg.HOMOGENEOUS_PV = False 
cfg.NumModelsPV = 1 if (cfg.HOMOGENEOUS_PV and not cfg.GAP) else cfg.NPV 
cfg.FactorTau, cfg.FactorKv3, cfg.FactorKv7 = 1, 1, 1   # To modify the activation curves for ion channels and the membrane time constant

# Stellate cells properties 
cfg.Mittal = True # uses the Mittal et al. model for Stellate cells
cfg.NSC= 4*cfg.NPV # Number of Stellate cells
cfg.HOMOGENEOUS_SC = False 
cfg.NumModelsSC = 1 if cfg.HOMOGENEOUS_SC else 454 # Load all the valid SC models
cfg.SCidx = 0 # Which model to load if using homogeneous population

if cfg.Mittal==False: 
    cfg.NumModelsSC=1
    cfg.HOMOGENEOUS_SC = True 

# Optogenetic drive                                                                                                                                                                                                                                                                                                                            
cfg.OPTODRIVE= True #True                                                  
cfg.g_sin = 0.*1e-3 # Optogenetic conductance for the inhibitory population                 

cfg.g_sinExc = 0.015 #nS, Optogenetic conductance for the excitatory population.
cfg.fsin=8  # Optogenetic sinusoidal stimulation, in Hz
cfg.delayStim = 5000 
cfg.durationStim = cfg.duration - cfg.delayStim #100

# Heterogeneous optogenetic drive for the PV+ cells
cfg.HETERDRIVE = True #False #True
np.random.seed(cfg.seeds['opto'])  # Fixed seed
cfg.drives = np.random.normal(loc=cfg.g_sin, scale=0.1*cfg.g_sin, size=cfg.NPV)

# Noisy stimulation
cfg.NOISE=False # If true, adds a random fluctuating current simulating in-vivo like noise. TODO: Adjust parameters for mEC region 
cfg.MeanENoise, cfg.MeanINoise, cfg.StdENoise, cfg.StdINoise = 0, 0*0.0001, 0.000001, 0.000001                                                                                                                                                                  

# Voltage clamp parameters
cfg.Clamp= True #True                                                                  
cfg.Vclamp = 0 #0 for IPSC, -70 for EPSC # Voltage at which clamp, if cfg.Clamp==True   
cfg.vcDelay = 5000  #ms        

# Synapses parameters                                                                                                                                                                                                                                 
cfg.gmsScale = 1 # Scaling for the synaptic conductances                                                                                                                                         
cfg.ggsScale = 1 # Scaling for gap junction conductances                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
cfg.GAP=True                                                                                                                                                                            
cfg.WODepression=False # Without short-term depression                                                                                                                                                                                                                                                                                                      
cfg.SYNAPSES = 'Hyper' # ['Hyper','Shunt','Uniform'] 
cfg.Esyn_inh = {'Hyper': -75., 'Shunt': -55., 'Uniform': 'uniform(-70,-55)'}
cfg.GapJunctProb, cfg.ChemycalConnProb = 0.2, 0.25 # Probability connections from FS PV+ to FS PV+                                                           
cfg.ConnProbIE, cfg.ConnProbEI = 0.2, 0.25 # FS PV+ to SC and SC to FS PV+ probability connections
cfg.Weight_E2I = 0.0008
cfg.WeightI2E = '(lognormal(1.65,2.17)*1e-3/0.3)*1.5'
cfg.tau_fallExc = 2
###############################################################################
## Recording and analysis
###############################################################################
# The index of cells are: First 100 are PV, the rest are SC 
cfg.ClampCells = [1,20,40,50,90,105,110,120,130,230,250,270,290,300,350]
recordedCells2 = [1,3,8,20,40,49,50,60,72,80,90,105,110,112,118,120,130,200,201,202,230,250,270,290,300,350,370,371,380,381,385,386,387,388,389,390]#[1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,104,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,240,235,240,245,250,255,260,265,270,275,280,285,290,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,350,385,386,387,388,389,390,400]
for_raster = ['FS', 'SC'] # Cells to include in the raster plot
offset = 5000
cells_interest = sorted({
    1,2,9,17,20,37,40,53,66,70,71,87,88,92,101,105,110,125,132,140,141,142,143,144,145,146,147,148,149,150,151,152,
    163,174,185,196,200,201,211,222,243,244,248,252,260,270,300,350,375,380,388
})

timeRange = [6250,cfg.duration]
cfg.recordCells = recordedCells2
cfg.recordTraces = {'V_soma':{'sec':'soma','loc':0.5,'var':'v'},
                    'SEClamp': {'sec': 'soma', 'loc': 0.5, 'stim':'Vclamp->Cells', 'var': 'i'}}  # Dict with traces to record
cfg.analysis['plotRaster'] = {'include': for_raster,'saveFig': True, 'figName': 'raster.png', 'timeRange': timeRange}       # Plot a raster
cfg.analysis['plotSpikeHist'] = {'include': ['FS', 'SC'], 'saveFig': True, 'timeRange': timeRange, 'binSize': 1, 'measure': 'rate'}                  # Plot a Spike Histogram
cfg.analysis['plotTraces'] = {'include': recordedCells2, 'saveFig': True, 'timeRange': timeRange}  # Plot recorded traces for this list of cells
