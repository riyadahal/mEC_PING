"""
defs.py

Auxiliar file with some cell definitions and functions used in the mEC model.
"""
from netpyne import specs
import numpy as np
import gc
import os
import csv
from netpyne.specs.dicts import Dict
import copy

def PVCell(cfg, gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s, SharedParams):
    netParamsAux = specs.NetParams()
    #To get a length from net capacitance most literature assumes neurons have 1uf/cm^2. This capacitance is in nF
    SAOrig=CapsOrig*1e-3
    SAGapJunct=CapsMod*1e-3
    #neuron does length in um so lets make it um^2
    SAumOrig=SAOrig*1e8
    SAumGapJunct=SAGapJunct*1e8
    #we're going to assume a diameter of 20um to get a radius of 10um
    #SA of a cyclinder(not including ends)=2*pi*r*L
    lengthsOrig=SAumOrig/(20*np.pi)
    lengthsGapJunct=SAumGapJunct/(20*np.pi)
    # Active conductances
    # Heterogeneous peak conductances for active currents. This gave me nanosiemens. I need siemens/cm^2 for it's conductance, so I need to multiple by e-9 to get nanosiemens then divide by the surface area in cm^2
    gLsSA = gLs*1e-9/SAOrig
    ConductWithGapJunctSA = ConductWithGapJunct*1e-9/SAGapJunct
    gNasSA = gNas*1e-9/SAOrig
    gNasWithGapJunctSA = gNas*1e-9/SAGapJunct
    gKv3sSA = gKv3s*1e-9/SAOrig
    gKv3sWithGapJunctSA = gKv3s*1e-9/SAGapJunct
    gKv7sSA = gKv7s*1e-9/SAOrig
    gKv7sWithGapJunctSA = gKv7s*1e-9/SAGapJunct
    # Fast Spiking PV+ Basket Cell from mEC
    for k in range(cfg.NumModelsPV):
        cellRule = {'conds': {'cellType': 'FS'}, 'diversityFraction': 1./cfg.NumModelsPV , 'secs': {}}  # cell rule dict        
        if cfg.OPTODRIVE==False:
            cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}} 
        else:
            cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}, 'pointps': {}}
            drive = cfg.g_sin if cfg.HETERDRIVE==False else cfg.drives[k]
            cellRule['secs']['soma']['pointps']['optodrive'] = {'mod' : 'optodrive', 'gsin': drive, 'Ese': 0, 'f': cfg.fsin, 'del': cfg.delayStim, 'dur': cfg.durationStim}
        if cfg.NOISE==True:
            cellRule['secs']['soma']['pointps']['InVivoNoise'] = {'mod' : 'Gfluct','g_e0': cfg.MeanENoise, 'g_i0': cfg.MeanINoise, 'std_e': cfg.StdENoise, 'std_i': cfg.StdINoise}
        if cfg.GAP==True:
            cellRule['secs']['soma']['geom'] = {'diam': 20, 'L': lengthsGapJunct[k], 'cm': 1}   # soma geometry
            cellRule['secs']['soma']['mechs']['pas'] = {'g': ConductWithGapJunctSA[k], 'e': ReversPotWithGapJunct[k]}
            cellRule['secs']['soma']['mechs']['naG'] = {'gbar': gNasWithGapJunctSA[k], 'thm1': thm1s[k], 'thh2': thh2s[k] }
            cellRule['secs']['soma']['mechs']['kv7'] = {'gbar': gKv7sWithGapJunctSA[k], 'tha1': tha1s[k], 'ka1': SharedParams[13], 'ka2':SharedParams[14] }
            cellRule['secs']['soma']['mechs']['kv3'] = {'gbar': gKv3sWithGapJunctSA[k], 'thn1': thn1s[k], 'kn1':SharedParams[10], 'kn2':SharedParams[11] }
        else:
            cellRule['secs']['soma']['geom'] = {'diam': 20, 'L': lengthsOrig[k], 'cm': 1}   # soma geometry
            cellRule['secs']['soma']['mechs']['pas'] = {'g': gLsSA[k], 'e': ELs[k]}
            cellRule['secs']['soma']['mechs']['naG'] = {'gbar': gNasSA[k], 'thm1': thm1s[k], 'thh2': thh2s[k] }
            cellRule['secs']['soma']['mechs']['kv7'] = {'gbar': gKv7sSA[k], 'tha1': tha1s[k], 'ka1': SharedParams[13],'ka2':SharedParams[14]}
            cellRule['secs']['soma']['mechs']['kv3'] = {'gbar': gKv3sSA[k], 'thn1': thn1s[k], 'kn1':SharedParams[10], 'kn2':SharedParams[11] }    
        cellRule['secs']['soma']['vinit'] = np.random.uniform(-75,-65) # set initial membrane potential
        netParamsAux.cellParams['FS_'+str(k)] = cellRule
    
    gc.collect()
    return netParamsAux.cellParams

def SCell_HH(cfg):
    netParamsAux = specs.NetParams()
    # Stellate Cell (Modified Hodgkin Huxley)
    SCcell = {'secs': {}, 'conds': {'cellType': 'SC'}}
    if cfg.OPTODRIVE==False:
        SCcell['secs']['soma'] = {'geom': {}, 'mechs': {}} 
    else:
        SCcell['secs']['soma'] = {'geom': {}, 'mechs': {}, 'pointps': {}}
        SCcell['secs']['soma']['pointps']['optodrive'] = {'mod' : 'optodrive', 'gsin': cfg.g_sinExc, 'Ese': 0, 'f': cfg.fsin, 'del': cfg.delayStim, 'dur': cfg.durationStim}
    if cfg.NOISE==True:
        SCcell['secs']['soma']['pointps']['InVivoNoise'] = {'mod' : 'Gfluct', 'g_e0': cfg.MeanENoise, 'g_i0': cfg.MeanINoise, 'std_e': cfg.StdENoise, 'std_i': cfg.StdINoise}
    SCcell['secs']['soma']['geom'] = {'diam': 18.2, 'L': 18.2, 'Ra': 150, 'cm':1}                           # soma geometry
    SCcell['secs']['soma']['mechs']['hh'] = {'gnabar': '0.12*normal(1,3e-2)', 'gkbar': '0.036*normal(1,3e-2)', 'gl': '0.0000357*normal(1.2,1e-1)', 'el': -68}  
    SCcell['secs']['soma']['vinit'] = np.random.uniform(-65,-58) # set initial membrane potential

    netParamsAux.cellParams['SC'] =  SCcell

    gc.collect()
    return netParamsAux.cellParams

def createMechsDict(cwd):
    paramsFile = csv.DictReader(open(cwd+'/SCcells/combinedparameters.csv'))
    paramsList = [row for row in paramsFile]

    mechList = []
    for k in range(len(paramsList)):
        mechsDict = {}
        mechsDict['NaT'] = {'gmax': float(paramsList[k]['gnafbar'])*1e-09,
                            'mvhalf': float(paramsList[k]['nafmvhalf']),
                            'mslope': float(paramsList[k]['nafmslope']),
                            'mfactor': float(paramsList[k]['nafmfactor']),
                            'hvhalf': float(paramsList[k]['nafhvhalf']),
                            'hslope': float(paramsList[k]['nafhslope']),
                            'hfactor': float(paramsList[k]['nafhfactor'])}
        
        mechsDict['KDR'] = {'gmax': float(paramsList[k]['gkdrbar'])*1e-09,
                            'nvhalf': float(paramsList[k]['kdrnvhalf']),
                            'nslope': float(paramsList[k]['kdrnslope']),
                            'nfactor': float(paramsList[k]['kdrnfactor'])}
        
        mechsDict['ih'] = {'gslowbar': float(paramsList[k]['ghcnsbar'])*1e-09,  
                            'gfastbar': float(paramsList[k]['hcnratio'])*(float(paramsList[k]['ghcnsbar'])*1e-09),
                            'mifo': float(paramsList[k]['hcnmifo']),
                            'miso': float(paramsList[k]['hcnmiso']),
                            'mifd': float(paramsList[k]['hcnmifd']),
                            'misd': float(paramsList[k]['hcnmisd']),
                            'ffactor': float(paramsList[k]['hcnffactor']),
                            'sfactor': float(paramsList[k]['hcnsfactor'])}
        
        mechsDict['NaP'] = {'gbar': float(paramsList[k]['gnapbar'])*1e-09, 
                            'mvhalf': float(paramsList[k]['napmvhalf']),
                            'mslope': float(paramsList[k]['napmslope']),
                            'mfactor': float(paramsList[k]['napmfactor']),
                            'hvhalf': float(paramsList[k]['naphvhalf']),
                            'hslope': float(paramsList[k]['naphslope']),
                            'hfactor': float(paramsList[k]['naphfactor'])}
        
        mechsDict['KA'] = {'gmax': float(paramsList[k]['gkabar'])*1e-09,    
                            'mvhalf': float(paramsList[k]['kamvhalf']),
                            'mslope': float(paramsList[k]['kamslope']),
                            'mfactor': float(paramsList[k]['kamfactor']),
                            'hvhalf': float(paramsList[k]['kahvhalf']),
                            'hslope': float(paramsList[k]['kahslope']),
                            'hfactor': float(paramsList[k]['kahfactor'])}
        
        mechsDict['HVA'] = {'gmax': float(paramsList[k]['ghvabar'])*1e-09,
                            'mvhalf': float(paramsList[k]['hvamvhalf']),
                            'mslope': float(paramsList[k]['hvamslope']),
                            'mfactor': float(paramsList[k]['hvamfactor']),
                            'hvhalf': float(paramsList[k]['hvahvhalf']),
                            'hslope': float(paramsList[k]['hvahslope']),
                            'hfactor': float(paramsList[k]['hvahfactor'])}
        
        mechsDict['LVA'] = {'glvabar': float(paramsList[k]['glvabar'])*1e-09,
                            'mvhalf': float(paramsList[k]['lvamvhalf']),
                            'mslope': float(paramsList[k]['lvamslope']),
                            'mfactor': float(paramsList[k]['lvamfactor']),
                            'hvhalf': float(paramsList[k]['lvahvhalf']),
                            'hslope': float(paramsList[k]['lvahslope']),
                            'hfactor': float(paramsList[k]['lvahfactor'])}
        
        mechsDict['km'] = {'gbar': float(paramsList[k]['gkmbar'])*1e-09,
                            'vhalfl': float(paramsList[k]['kmvhalfl']),
                            'kl': float(paramsList[k]['kmkl']),
                            'mfactor': float(paramsList[k]['kmmfactor'])}

        mechsDict['skkin'] = {'gbar': float(paramsList[k]['gskbar'])*1e-09}

        mechsDict['pas'] = {'g': 1/float(paramsList[k]['Rm'])}
        mechsDict['geom'] = {'cm': float(paramsList[k]['CM'])}
        mechsDict['cad'] =  {'tauca': float(paramsList[k]['tauca'])}
        mechList.append(mechsDict)

    return mechList

def SC_Mittal(cwd, cfg):
    netParamsAux = specs.NetParams()
    cellRule = netParamsAux.importCellParams(label='SC_0', conds={'cellType': 'SC'},
    fileName=cwd+'/SCcells/Stochastic_Bifurcations_Mittal_Narayanan/Figure3-6_SCmodel/Noise_Osc_july.hoc',
    cellName='SC', cellInstance = True)
    mechList = createMechsDict(cwd)

    netParamsAux2 = specs.NetParams()

    for k in range(cfg.NumModelsSC):
        cellRuleAux = copy.deepcopy(cellRule.todict())
        if cfg.HOMOGENEOUS_SC: k = cfg.SCidx
        paramsAux = mechList[k]

        if cfg.OPTODRIVE==True:
            cellRuleAux['secs']['soma']['pointps'] = {}
            cellRuleAux['secs']['soma']['pointps']['optodrive'] = {'mod' : 'optodrive', 'gsin': cfg.g_sinExc, 'Ese': 0, 'f': cfg.fsin, 'del': cfg.delayStim, 'dur': cfg.durationStim}

        if cfg.NOISE==True:
            cellRuleAux['secs']['soma']['pointps']['InVivoNoise'] = {'mod' : 'Gfluct', 'g_e0': cfg.MeanENoise, 'g_i0': cfg.MeanINoise, 'std_e': cfg.StdENoise, 'std_i': cfg.StdINoise}
    

        cellRuleAux['diversityFraction'] = 1./cfg.NumModelsSC

        for mechName in cellRuleAux['secs']['soma']['mechs'].keys():
            if mechName in paramsAux.keys():
                cellRuleAux['secs']['soma']['mechs'][mechName].update(paramsAux[mechName])
            else:
                print(f"Mechanism {mechName} not found in params for SC model {k}.") 

        cellRuleAux['secs']['soma']['geom'].update(paramsAux['geom'])
        netParamsAux2.cellParams['SC_'+str(k)] = copy.deepcopy(cellRuleAux)

    gc.collect()
    return netParamsAux2.cellParams

def GapJunctSpatialConnectivity(sim, netParams, gLmin=2., ggapHigh=1.3, full_dist=False, sigma=.4):
  from scipy.stats import truncnorm
  NumNeurons = len(sim.net.pops['FS'].cellGids)
  Ngaps = int(NumNeurons*sim.cfg.GapJunctProb)
  ConductOrig = netParams.ConductOrig
  ReversPotOrig = netParams.ReversPotOrig
  ConductWithGapJunct = np.copy(ConductOrig)
  ReversPotWithGapJunctAux = ConductOrig*ReversPotOrig
  ReversPotWithGapJunct = np.zeros(np.shape(ReversPotOrig))
  AreaOrig = np.array([sim.net.cells[i].secs.soma.geom.diam*np.pi*sim.net.cells[i].secs.soma.geom.L*1e-2 for i in sim.net.pops['FS'].cellGids])
  synsgj = []
  ggs = []
  for presynGJ in range(NumNeurons):
    NumGapJuncts = 0
    NumIter = 0
    while NumGapJuncts <= int(Ngaps/2.) and NumIter<2*NumNeurons: # We start adding ggaps with Ngaps/2 partners (plus bidirectional, approx Ngaps partners) to each pre. We put a maximum in the iteration number, just in case the amount of possible gap junctions is less than Ngaps/2
      postsynGJ = np.random.randint(0, NumNeurons)
      distancePrePost = np.sqrt((sim.net.cells[presynGJ].tags['x'] - sim.net.cells[postsynGJ].tags['x'])**2 + (sim.net.cells[presynGJ].tags['y'] - sim.net.cells[postsynGJ].tags['y'])**2 + (sim.net.cells[presynGJ].tags['z'] - sim.net.cells[postsynGJ].tags['z'])**2) 
      if distancePrePost>sim.cfg.GapJunctMaxDist: continue # To avoid gap junctions between cells that are too far away
      if presynGJ==postsynGJ: continue # To avoid self-gap junctions
      if (presynGJ,postsynGJ) not in synsgj:
        if full_dist==True: ggap = [ sigma*truncnorm.rvs(a=.0,b=1.5/sigma) if np.random.rand()<.75 else ggapHigh ][0]
        if full_dist==False: ggap = sigma*truncnorm.rvs(a=.0,b=1.5/sigma)
        if ConductWithGapJunct[presynGJ]-ggap>=gLmin and ConductWithGapJunct[postsynGJ]-ggap>=gLmin:
          #Netpyne auto makes these junctions bidirectional so we don't want to read the direction in twice
          synsgj.append((presynGJ,postsynGJ)) # Bi-directional gap junctions
          #synsgj.append((postsynGJ,presynGJ)) # Bi-directional gap junctions
          ggs.append(ggap); #ggs.append(ggap) # Bi-directional gap junctions
          ConductWithGapJunct[presynGJ] -= ggap
          ConductWithGapJunct[postsynGJ] -= ggap # We reduce Leakage conductance
          ReversPotWithGapJunctAux[presynGJ] -= ggap*ReversPotOrig[postsynGJ]
          ReversPotWithGapJunctAux[postsynGJ] -= ggap*ReversPotOrig[presynGJ]
          NumGapJuncts+=1
      NumIter+=1

    ReversPotWithGapJunct[presynGJ] = ReversPotWithGapJunctAux[presynGJ]/ConductWithGapJunct[presynGJ]
  # Update values
  for i in sim.net.pops['FS'].cellGids:
    sim.net.cells[i].secs.soma.mechs.pas.g = ConductWithGapJunct[i]/AreaOrig[i]*1e-3
    sim.net.cells[i].secs.soma.mechs.pas.e = ReversPotWithGapJunct[i]
  ggs = [i*1e-3 for i in ggs] 
  if sim.cfg.GAP==True:
    # Add gap junction parameters to existing netParams
    sim.net.params.synMechParams['gap'] = {'mod': 'ElectSyn', 'g': 1}
    sim.net.params.connParams['FS->FS_gap'] = {
                'preConds': {'pop': 'FS'}, 
                'postConds': {'pop': 'FS'}, 
                'connList': synsgj,
                'gapJunction': True,
                'sec': 'soma',
                'weight': ggs,
                'synMech': 'gap',
                'delay': 0
    }

  return None
