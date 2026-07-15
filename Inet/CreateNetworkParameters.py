import numpy as np
from scipy.stats import truncnorm
import os

# See end of file for Usage
#Brian2 dependence has been removed from this set of simulations so the seeds are different here
def gen_syns_samepop(N,p):
    """
    Brian2-compatible replacement for:
    Synapses.connect(p=p, condition="i!=j")

    Brian2 rewrites this call as the generator expression:
        j = "_k for _k in sample(N_post, p=p) if i != _k"
    i.e. for every presynaptic index, it samples every postsynaptic index with
    probability p and only then filters out self-connections.
    """
    presynaptic = []
    postsynaptic = []

    for i in range(N):
        if p == 1:
            js = np.arange(N, dtype=np.int32)
        else:
            js = np.flatnonzero(np.random.random(N) < p).astype(np.int32)
        js = js[js != i]
        if js.size == 0:
            continue
        presynaptic.append(np.full(js.size, i, dtype=np.int32))
        postsynaptic.append(js)

    if not presynaptic:
        empty = np.array([], dtype=np.int32)
        return empty, empty

    return np.concatenate(presynaptic), np.concatenate(postsynaptic)

def gen_unif_ds(Nsyns,dmin=.6,dmax=1.): return dmin+(dmax-dmin)*np.random.random(size=(1,Nsyns))

def gen_lognormal_gms(Nsyns,mean=0.,sigma=1.): return np.random.lognormal(mean=mean,sigma=sigma,size=(1,Nsyns))

def _load_network_model_parameters():
  current_dir = os.getcwd()
  module_dir = os.path.dirname(os.path.abspath(__file__))
  candidate_dirs = [current_dir, os.path.join(current_dir, "Inet"), module_dir]

  def load_from_candidates(filename):
    for directory in candidate_dirs:
      path = os.path.join(directory, filename)
      if os.path.exists(path):
        return np.loadtxt(path)
    raise FileNotFoundError(f"Could not find '{filename}' in any expected location.")

  gLs = load_from_candidates("gLsLast.dat")
  ELs = load_from_candidates("ELsLast.dat")
  taums = load_from_candidates("taumsLast.dat")
  active_params = np.transpose(load_from_candidates("ActiveParamsLast.dat"))
  return gLs, ELs, taums, active_params

def _match_population_size(values, NumNeurons):
  values = np.asarray(values)
  NumModels = len(values)
  if NumNeurons > NumModels:
    q, mod = divmod(NumNeurons, NumModels)
    return np.concatenate((np.tile(values, q), values[:mod]))
  return values[:NumNeurons]

def ConductAndReversPotsWithGapJunct(ConductOrig,ReversPotOrig,gLmin=2.,ggapHigh=1.3,full_dist=False,Ngaps=20,sigma=.4):
  NumNeurons = len(ConductOrig)
  ConductWithGapJunct = np.copy(ConductOrig)
  ReversPotWithGapJunctAux = ConductOrig*ReversPotOrig
  ReversPotWithGapJunct = np.zeros(np.shape(ReversPotOrig))
  synsgj = []
  ggs = []
  for presynGJ in range(NumNeurons):
    NumGapJuncts = 0
    NumIter = 0
    while NumGapJuncts <= int(Ngaps/2.) and NumIter<2*NumNeurons: # We start adding ggaps with Ngaps/2 partners (plus bidirectional, approx Ngaps partners) to each pre. We put a maximum in the iteration number, just in case the amount of possible gap junctions is less than Ngaps/2
      postsynGJ = np.random.randint(0,NumNeurons)
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

  return ConductWithGapJunct, ReversPotWithGapJunct, ggs, synsgj

def NetworkParams(NumNeurons=100,
FactorTau=1,
FactorKv3=1,
FactorKv7=1,
GapJunctProb=0.2,
ChemycalConnProb=0.36,
delaymin=.6,
delaymax=1.,
meangms=0.,
sigmagms=1.,
homogeneous=False,
NeuronModel=43,
randomseed = 7894 # FIX IT FOR REPRODUCIBILITY
):
  np.random.seed(randomseed)
  gLs, ELs, taums, active_params = _load_network_model_parameters()
  taums = FactorTau*taums
  gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s = active_params

  if homogeneous:
    gLs = np.array([gLs[NeuronModel]])
    ELs = np.array([ELs[NeuronModel]])
    taums = np.array([taums[NeuronModel]])
    gNas = np.array([gNas[NeuronModel]])
    gKv3s = np.array([gKv3s[NeuronModel]])
    gKv7s = np.array([gKv7s[NeuronModel]])
    thm1s = np.array([thm1s[NeuronModel]])
    thh2s = np.array([thh2s[NeuronModel]])
    thn1s = np.array([thn1s[NeuronModel]])
    tha1s = np.array([tha1s[NeuronModel]])

  gLs = _match_population_size(gLs, NumNeurons)
  ELs = _match_population_size(ELs, NumNeurons)
  taums = _match_population_size(taums, NumNeurons)
  gNas = _match_population_size(gNas, NumNeurons)
  gKv3s = _match_population_size(gKv3s, NumNeurons)
  gKv7s = _match_population_size(gKv7s, NumNeurons)
  thm1s = _match_population_size(thm1s, NumNeurons)
  thh2s = _match_population_size(thh2s, NumNeurons)
  thn1s = _match_population_size(thn1s, NumNeurons)
  tha1s = _match_population_size(tha1s, NumNeurons)

  ENa=50; sigm1=4.; km2=.1; sigm2=13.; kh1=.012; kh2=.2; sigh1=-20.; sigh2=3.5; #Sodium current parameters
  EK=-90; sign1=12.; kn1 = FactorKv3*1.; kn2=FactorKv3*.001; sign2=-8.5; #Kv3 current parameters
  ka1 = FactorKv7*1; ka2=FactorKv7*.02; siga1=12.; siga2=-80.; #Kv7 current parameters (the rest of parameters are as in Kv3)
  SharedParams = [ENa, sigm1, km2, sigm2, kh1, kh2, sigh1, sigh2, EK, sign1, kn1, kn2, sign2, ka1, ka2, siga1, siga2]

  if GapJunctProb < 1e-3: 
    Ngaps=0
    ConductWithGapJunct, ReversPotWithGapJunct, ggs, synsgj = np.copy(gLs), np.copy(ELs), [], []
  else:
    Ngaps = int(NumNeurons*GapJunctProb)
    ConductWithGapJunct, ReversPotWithGapJunct, ggs, synsgj = ConductAndReversPotsWithGapJunct(gLs, ELs,Ngaps=Ngaps)
  
  syns = gen_syns_samepop(NumNeurons,ChemycalConnProb)
  Nsyns = len(syns[0])  
  delays = gen_unif_ds(Nsyns,dmin=delaymin,dmax=delaymax)[0]
  gms = gen_lognormal_gms(Nsyns,mean=meangms,sigma=sigmagms)[0]

  syns = [(i,j) for (i,j) in zip(syns[0],syns[1])]

  CapsOrig =  taums*gLs/1e3
  CapsMod =  taums*ConductWithGapJunct/1e3
  #print(CapsOrig); quit()
  gms = [i*1e-3 for i in gms]
  delays = [i for i in delays]
  ggs = [i*1e-3 for i in ggs] 

  return gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s, SharedParams, syns, delays, gms, synsgj, ggs

def ViaModelParams(NumNeurons=100,
FactorTau=1,
FactorKv3=1,
FactorKv7=1,
GapJunctProb=0.2,
homogeneous=False,
NeuronModel=43,
randomseed = 7894 # FIX IT FOR REPRODUCIBILITY
):
  np.random.seed(randomseed)
  gLs, ELs, taums, active_params = _load_network_model_parameters()
  taums = FactorTau*taums
  gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s = active_params

  if homogeneous:
    gLs = np.array([gLs[NeuronModel]])
    ELs = np.array([ELs[NeuronModel]])
    taums = np.array([taums[NeuronModel]])
    gNas = np.array([gNas[NeuronModel]])
    gKv3s = np.array([gKv3s[NeuronModel]])
    gKv7s = np.array([gKv7s[NeuronModel]])
    thm1s = np.array([thm1s[NeuronModel]])
    thh2s = np.array([thh2s[NeuronModel]])
    thn1s = np.array([thn1s[NeuronModel]])
    tha1s = np.array([tha1s[NeuronModel]])

  gLs = _match_population_size(gLs, NumNeurons)
  ELs = _match_population_size(ELs, NumNeurons)
  taums = _match_population_size(taums, NumNeurons)
  gNas = _match_population_size(gNas, NumNeurons)
  gKv3s = _match_population_size(gKv3s, NumNeurons)
  gKv7s = _match_population_size(gKv7s, NumNeurons)
  thm1s = _match_population_size(thm1s, NumNeurons)
  thh2s = _match_population_size(thh2s, NumNeurons)
  thn1s = _match_population_size(thn1s, NumNeurons)
  tha1s = _match_population_size(tha1s, NumNeurons)

  ENa=50; sigm1=4.; km2=.1; sigm2=13.; kh1=.012; kh2=.2; sigh1=-20.; sigh2=3.5; #Sodium current parameters
  EK=-90; sign1=12.; kn1 = FactorKv3*1.; kn2=FactorKv3*.001; sign2=-8.5; #Kv3 current parameters
  ka1 = FactorKv7*1; ka2=FactorKv7*.02; siga1=12.; siga2=-80.; #Kv7 current parameters (the rest of parameters are as in Kv3)
  SharedParams = [ENa, sigm1, km2, sigm2, kh1, kh2, sigh1, sigh2, EK, sign1, kn1, kn2, sign2, ka1, ka2, siga1, siga2]

  if GapJunctProb < 1e-3: 
    Ngaps=0
    ConductWithGapJunct, ReversPotWithGapJunct, ggs, synsgj = np.copy(gLs), np.copy(ELs), [], []
  else:
    Ngaps = int(NumNeurons*GapJunctProb)
    ConductWithGapJunct, ReversPotWithGapJunct, ggs, synsgj = ConductAndReversPotsWithGapJunct(gLs, ELs,Ngaps=Ngaps)
  ggs = [i*1e-3 for i in ggs] 

  CapsOrig =  taums*gLs/1e3
  CapsMod =  taums*ConductWithGapJunct/1e3

  return gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s, SharedParams, synsgj, ggs

def NetConnectivity(NumNeurons=100,
ChemycalConnProb=0.36,
delaymin=.6,
delaymax=1.,
meangms=0.,
sigmagms=1.,
randomseed = 7894 # FIX IT FOR REPRODUCIBILITY
):
  # Keep the historical behavior: this function advances the current NumPy RNG
  # stream instead of reseeding, so it stays equivalent to NetworkParams when
  # called after ViaModelParams.
  syns = gen_syns_samepop(NumNeurons,ChemycalConnProb)
  Nsyns = len(syns[0])  
  delays = gen_unif_ds(Nsyns,dmin=delaymin,dmax=delaymax)[0]
  gms = gen_lognormal_gms(Nsyns,mean=meangms,sigma=sigmagms)[0]

  syns = [(i,j) for (i,j) in zip(syns[0],syns[1])]

  gms = [i*1e-3 for i in gms]
  delays = [i for i in delays]

  return syns, delays, gms

### USAGE
#NumNetw = 2
#for k in range(NumNetw):
#  gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, syns, synsgj, SharedParams = NetParams(NumNeurons=100, FactorTau=1, FactorKv3=1, FactorKv7=1, GapJunctProb=0.1, ChemycalConnProb=0.36, delaymin=.6, delaymax=1., meangms=0., sigmagms=1.)

#  print(gLs[:2], ELs[:2], CapsOrig[:2], ConductWithGapJunct[:2], ReversPotWithGapJunct[:2], CapsMod[:2])
####
# gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s, SharedParams, syns, delays, gms, synsgj, ggs = NetworkParams()
# OR EQUIVALENTLY
# gLs2, ELs2, CapsOrig2, ConductWithGapJunct2, ReversPotWithGapJunct2, CapsMod2, gNas2, gKv3s2, gKv7s2, thm1s2, thh2s2, thn1s2, tha1s2, SharedParams2, synsgj2, ggs2 = ViaModelParams()
# syns2, delays2, gms2 = NetConnectivity()

# assert set(gLs) == set(gLs2), f"Assertion failed: gLs != gLs2 as sets\n  set(gLs): {set(gLs)}\n  set(gLs2): {set(gLs2)}"
# assert set(ELs) == set(ELs2), f"Assertion failed: gLs != gLs2 as sets\n  set(gLs): {set(gLs)}\n  set(gLs2): {set(gLs2)}"
# assert set(CapsOrig) == set(CapsOrig2), f"Assertion failed: gLs != gLs2 as sets\n  set(gLs): {set(gLs)}\n  set(gLs2): {set(gLs2)}"
# assert set(ConductWithGapJunct) == set(ConductWithGapJunct2), f"Assertion failed: gLs != gLs2 as sets\n  set(gLs): {set(gLs)}\n  set(gLs2): {set(gLs2)}"
# assert set(ReversPotWithGapJunct) == set(ReversPotWithGapJunct2), f"Assertion failed: gLs != gLs2 as sets\n  set(gLs): {set(gLs)}\n  set(gLs2): {set(gLs2)}"
# assert set(CapsMod) == set(CapsMod2), f"Assertion failed: gLs != gLs2 as sets\n  set(gLs): {set(gLs)}\n  set(gLs2): {set(gLs2)}"
# assert set(syns) == set(syns2), f"Assertion failed: gLs != gLs2 as sets\n  set(gLs): {set(syns)}\n  set(gLs2): {set(syns2)}"
# assert set(synsgj) == set(synsgj2), f"Assertion failed: gLs != gLs2 as sets\n  set(gLs): {set(gLs)}\n  set(gLs2): {set(gLs2)}"
# assert set(ggs) == set(ggs2), f"Assertion failed: gLs != gLs2 as sets\n  set(gLs): {set(gLs)}\n  set(gLs2): {set(gLs2)}"



