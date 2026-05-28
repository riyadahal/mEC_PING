import numpy as np
#from funcsAux import *
import matplotlib.pyplot as plt
import glob, os
from netpyne.analysis.tools import loadData
from itertools import compress
import matplotlib.gridspec as gridspec
gsin = 5
f = 8
dt_rec = 1e-2 # Down sample the rates
#minFreq=30.
#maxFreq=350.
#maxPower= 500 # 100. # Sets the maximum for the scalogram. Varies according to the power of the signal. Change it for better visualization
#colorbar=True
#numBins = 2*125+1
###################################
#fs = np.arange(minFreq,maxFreq,3.)
#levelsWav = np.linspace(0., maxPower, 20, endpoint=True)

i = 0
spktFSAux = {}
spktSCAux = {}
FScellTrace={}
SCcellTrace={}

fileID = '1215_Traces/'
if not os.path.exists(fileID):
    os.makedirs(fileID)


for file in glob.glob("../tau2runIPSC/*.pkl"):
    fileInfo = loadData(file)
    fileName = fileInfo['simConfig']['filename'][-5:]
    print(fileName)
    sim_time = fileInfo['simConfig']['duration']
    FScellTrace['0'] = fileInfo['simData']['V_soma']['cell_8']
    FScellTrace['1'] = fileInfo['simData']['V_soma']['cell_201']
    FScellTrace['2'] = fileInfo['simData']['V_soma']['cell_112']
    #maskFS = np.array(fileInfo['simData']['spkid'])< 100
    #spktFSAux['0'] = np.array( list( compress(fileInfo['simData']['spkt'], maskFS) ) )
    i += 1

#spktFS0 = spktFSAux['0']
cellTrFS0 = FScellTrace['0']
cellTrFS1 = FScellTrace['1']
cellTrFS2 = FScellTrace['2']

tSim = np.linspace(0,int(sim_time)+1,int(sim_time)*100)

fig = plt.figure(figsize=(10, 8))
outer = gridspec.GridSpec(1, 1, wspace=0.3, hspace=0.1)

inner = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=outer[0], wspace=0.1, hspace=0.1, height_ratios=[1,4,4,4])


ax0 = plt.Subplot(fig, inner[0])
t = np.linspace(sim_time-4.*125.,sim_time,2*125+1)
ax0.plot(t, gsin*np.sin(2*np.pi*f*t*1e-3-np.pi/2),color = 'k')
ax0.set_xlim(sim_time-4.*125.,sim_time)
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.set_xticks([])
ax0.set_yticks([])
fig.add_subplot(ax0)

ax1 = plt.Subplot(fig, inner[1], label='FS - Stimulate only FS')
y = cellTrFS0
ax1.plot(tSim, y, 'b-', label='FS - Stimulate only FS')
ax1.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax1)

ax2 = plt.Subplot(fig, inner[2], label='FS - Stimulate FS and SC')
y = cellTrFS1
ax2.plot(tSim, y, 'b-', label='FS - Stimulate FS and SC')
ax2.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax2)

ax3 = plt.Subplot(fig, inner[3], label='SC - Stimulate FS and SC')
y = cellTrFS2
ax3.plot(tSim, y, 'b-', label='SC - Stimulate FS and SC')
ax3.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax3)

fig.tight_layout()
fig.savefig(fileID+'Traces.png', bbox_inches='tight',dpi = 300)
fig.savefig(fileID+'Traces.eps', bbox_inches='tight', dpi = 300)
