import os
import sys

repo_root = os.path.abspath(os.getcwd())
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

import Inet.CreateNetworkParameters as Inet
import src.defs as defs
import numpy as np
from netpyne import specs
import os
try:
    from __main__ import cfg
except:
    from src.cfg import cfg

cwd = os.getcwd()  # Get current working directory

# Default network parameters
netParams = specs.NetParams()  # object of class NetParams to store the network parameters
netParams.defaultDelay = 0 #This is because it creates gap junctions with defaultDelay = 1 ms if not 
netParams.defaultThreshold = -30.0

###############################################################################
## Create and load parameters for the PV network
###############################################################################
# gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s, SharedParams, syns, delays, gms, synsgj, ggs = Inet.NetworkParams(NumNeurons=cfg.NPV, FactorTau=cfg.FactorTau, FactorKv3=cfg.FactorKv3, FactorKv7=cfg.FactorKv7, GapJunctProb=cfg.GapJunctProb, ChemycalConnProb=cfg.ChemycalConnProb, delaymin=.6, delaymax=1., meangms=0., sigmagms=1.,homogeneous=cfg.HOMOGENEOUS_PV, randomseed = cfg.seeds['brian2'])

gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s, SharedParams, synsgj, ggs = Inet.ViaModelParams(NumNeurons=cfg.NPV, FactorTau=cfg.FactorTau, FactorKv3=cfg.FactorKv3, FactorKv7=cfg.FactorKv7, GapJunctProb=cfg.GapJunctProb,homogeneous=cfg.HOMOGENEOUS_PV, randomseed = cfg.seeds['brian2'])
syns, delays, gms = Inet.NetConnectivity(ChemycalConnProb=cfg.ChemycalConnProb, delaymin=.6, delaymax=1., meangms=0., sigmagms=1.)

###############################################################################
## Cell types
###############################################################################
cellParamsSC = defs.SCell_HH(cfg)
cellParamsPV = defs.PVCell(cfg, gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s, SharedParams)
cellParamsSC_Mittal = defs.SC_Mittal(cwd, cfg)

netParams.cellParams = cellParamsPV | cellParamsSC_Mittal if cfg.Mittal else cellParamsPV | cellParamsSC # Combine all dictionaries

###############################################################################
# NETWORK PARAMETERS
###############################################################################
# Population parameters
netParams.popParams['FS'] = {'cellType': 'FS', 'numCells': cfg.NPV, 'diversity': True} # add dict with params for this pop
if cfg.HOMOGENEOUS_SC:
    netParams.popParams['SC'] = {'cellType': 'SC', 'numCells': cfg.NSC} 
else:
    netParams.popParams['SC'] = {'cellType': 'SC', 'numCells': cfg.NSC, 'diversity': True} # add dict with params for this pop
#netParams.popParams['PYR'] = {'cellType': 'PYR', 'numCells': 1} 

###############################################################################
## VoltageClamp
###############################################################################
if cfg.Clamp==True:
	netParams.stimSourceParams['Vclamp'] = {'type': 'SEClamp', 'dur1': 1e9, 'amp1': cfg.Vclamp, 'rs': 1e-5}
    # Stimulation mapping parameters
	netParams.stimTargetParams['Vclamp->Cells'] = {
        'source': 'Vclamp',
        'sec': 'soma',
        'loc': 0.5,
        'conds': {'cellList': cfg.ClampCells}}
#------------------------------------------------------------------------------
# Current inputs (IClamp)
#------------------------------------------------------------------------------
if cfg.addIClamp:
    for key in [k for k in dir(cfg) if k.startswith('IClamp')]:
        params = getattr(cfg, key, None)
        [pop,sec,loc,start,dur,amp] = [params[s] for s in ['pop','sec','loc','start','dur','amp']]

        #cfg.analysis['plotTraces']['include'].append((pop,0))  # record that pop

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'IClamp', 'delay': start, 'dur': dur, 'amp': amp}
        
        # connect stim source to target
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source': key, 
            'conds': {'pop': pop},
            'sec': sec, 
            'loc': loc}

#------------------------------------------------------------------------------
# Current inputs (IClamp)
#------------------------------------------------------------------------------
if cfg.addIClamp:
    for key in [k for k in dir(cfg) if k.startswith('IClamp')]:
        params = getattr(cfg, key, None)
        [pop,sec,loc,start,dur,amp] = [params[s] for s in ['pop','sec','loc','start','dur','amp']]

        #cfg.analysis['plotTraces']['include'].append((pop,0))  # record that pop

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'IClamp', 'delay': start, 'dur': dur, 'amp': amp}
        
        # connect stim source to target
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source': key, 
            'conds': {'pop': pop},
            'sec': sec, 
            'loc': loc}


###############################################################################
## Synaptic mechs
###############################################################################

## Inhibitory synapses FS-> FS
tau_rise=0.3
tau_fall=2.
c_fall = 1/tau_fall
c_rise = 1/tau_rise
f = 1/( np.exp(-c_fall*np.log(c_rise/c_fall)/(c_rise-c_fall)) - np.exp(-c_rise*np.log(c_rise/c_fall)/(c_rise-c_fall)) )
U_SE=0.3; tau_d=100.
# Let's add synapses  
if cfg.WODepression==True:
    netParams.synMechParams['inhFSFS'] = {'mod': 'synact', 'tau_rise': tau_rise, 'tau_fall': tau_fall, 'f' : f, 'Es': cfg.Esyn_inh[cfg.SYNAPSES]}        
else:
    gms = [i/U_SE for i in gms]
    netParams.synMechParams['inhFSFS'] = {'mod': 'synactdep', 'tau_rise': tau_rise, 'tau_fall': tau_fall, 'f' : f, 'xs' : 1, 'U_SE' : U_SE, 'tau_d' : tau_d, 'Es': cfg.Esyn_inh[cfg.SYNAPSES]}
gms = [cfg.gmsScale*i for i in gms]
ggs = [cfg.ggsScale*i for i in ggs]
# Connectivity parameters
netParams.connParams['FS->FS_chem'] = {
        'preConds': {'pop': 'FS'},         # presynaptic conditions
        'postConds': {'pop': 'FS'},        # postsynaptic conditions
        'sec':'soma',
        'connList': syns,
        'weight': gms,                      # weight of each connection. I'm subbing weight for conductance, gms was in nanosiemens and needs to be converted to uS
        'synMech': 'inhFSFS',                   # target inh synapse
         'delay': delays}                    # delay
if cfg.GAP==True:
    netParams.synMechParams['gap'] = {'mod': 'ElectSyn', 'g': 1}
    # Connectivity parameters
    netParams.connParams['FS->FS_gap'] = {
            'preConds': {'pop': 'FS'},         # presynaptic conditions
            'postConds': {'pop': 'FS'},        # postsynaptic conditions
            'connList': synsgj,
            'gapJunction': True, #Netpyne auto makes these junctions bidirectional so we don't want to read the direction in twice
            'sec':'soma',
            'weight': ggs,                      # weight of each connection. I'm subbing weight for conductance, gms was in nanosiemens and needs to be converted to uS
            'synMech': 'gap',                   # target inh synapse
            'delay': 0}
###############################################################################
## Inhibitory synapses FS-> SC

tau_riseExc=0.4
tau_fallExc= cfg.tau_fallExc
c_fall = 1./tau_fallExc; c_rise = 1./tau_riseExc
norm_synExc = 1./( np.exp(-c_fall*np.log(c_rise/c_fall)/(c_rise-c_fall)) - np.exp(-c_rise*np.log(c_rise/c_fall)/(c_rise-c_fall)) )
netParams.synMechParams['inhFSSC'] = {'mod': 'synactdep', 'tau_rise': tau_riseExc, 'tau_fall': tau_fallExc, 'f' : norm_synExc, 'xs' : 1, 'U_SE' : U_SE, 'tau_d' : tau_d, 'Es': -65.}
# Connectivity parameters
netParams.connParams['FS->SC'] = {
        'preConds': {'pop': 'FS'},         # presynaptic conditions
        'postConds': {'pop': 'SC'},        # postsynaptic conditions
        'sec':'soma',
        'probability': cfg.ConnProbIE,
        #I'm subbing weight for conductance, gms was in nanosiemens and needs to be converted to uS
        'weight': cfg.WeightI2E, #'(lognormal(1.65,2.17)*1e-3/0.3)*1.5',                      # weight of each connection. Take into account that numpy and NEURON arguments are different (numpy args are mean and std for subjacent normal distribution, not the lognorm as in NEURON)
        'synMech': 'inhFSSC',                   # target inh synapse
        'delay': '0.6+(1-0.6)*uniform(0,1)'}                    # delay

###############################################################################
## Excitatory synapses SC -> FS

netParams.synMechParams['AMPA'] = {'mod': 'ExpSyn', 'tau': 1, 'e': 0}  # excitatory synaptic mechanism
# Synaptic mechanism parameters
netParams.synMechParams['NMDA'] = {'mod': 'Exp2Syn', 'tau1': 0.1, 'tau2': 5.0, 'e': 0}  # NMDA synaptic mechanism (NOT IMPLEMENTED YET)
# Connectivity parameters
netParams.connParams['SC->FS'] = {
        'preConds': {'pop': 'SC'},         # presynaptic conditions
        'postConds': {'pop': 'FS'},        # postsynaptic conditions
        'sec':'soma',
        'probability': cfg.ConnProbEI,
        #am subbing weight for conductance, gms was in nanosiemens and needs to be converted to uS
        'weight': cfg.Weight_E2I,#'0.0039*lognormal(0.01,1e-1)',       #'0.8*lognormal(0.01,1e-1)'               # weight of each connection#gms*0.001
        'synMech': 'AMPA',                   # target exc synapse
         'delay': '0.6+(1-0.6)*uniform(0,1)'} #'0.6+(1-0.6)*uniform(0,1)'                   # delay

