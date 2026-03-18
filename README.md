# mEC Model
 
 This code was written using NetPyNE (Networks with Python and NEURON). For more information about the tool, see [NetPyNE docs](http://doc.netpyne.org/)

## Download the Repository

To clone the dev branch of the repository, open a terminal in the directory where you'd like to store the project and run:

```bash
git clone https://github.com/riyadahal/mEC_PING.git
````

After that, make sure to move to the repo folder using `cd mEC_PING`.

## Installation

Install the required Python packages using:

```bash
pip install -r requirements.txt
```

## Simulations Setup

Add the root directory to `PYTHONPATH`, and compile the mechanisms:

```bash
export PYTHONPATH=PYTHONPATH:$PWD
nrnivmodl mod
```

### Setting up the parameters

In `spatialModel/cfg.py` you can modify all the parameters, variables and conditions for the simulation. For example, you can change the values of the flags to have sinusoidal optogenetic drive, choose the stellate cell model, add gap junction connectivity between the PV cells, change the recorded cells, simulation time, etc.

## ▶Running Simulations

### Run a Single Simulation

All simulation parameters are defined in the `spatialModel/cfg.py` file. To run a single simulation:

```bash
python -u spatialModel/init.py
```
### Updates
Some new parameters added and changes to the optogenetic stimulus are:
1. Optogenetic drive now allows two new parameters (in the same fashion as `IClamp`): `cfg.delayStim` and 
`cfg.durationStim` which allow to give a delay to the stimulation and to finish it early too
2. Also, when `f=0` in the optogenetic drive, the current is a constant of value `g_sin` or `g_sinExc`
3. Now the current injected in the voltage-clamped neurons is being recorded and plotted. Current is in nA
4. The flag `cfg.PlotWavelet` allows to plot the injected current for this cells, to filter the signal and to plot its wavelet transform. It makes simulations end slower, so use it if you need it.

### Run Batch Simulations (Experimental)

To run multiple simulations on the terminal using the batch submitter:

```bash
python -u spatialModel/batch.py
```
The function `runNetworks(NumNetworks=X)` will run X number of networks with different connectivity but the same parameters as in `spatialModel/cfg.py`. It is possible to create another function with another set of 
```bash
params['property'] = [list of values]
```
where that property is any of the ones defined in `spatialModel/cfg.py`. That will create automatically all possible combination of values between all the params (be careful since the number of combinations grows exponentially). See for example the function `runDifferentReversalPot()`

> ⚠️ Note: Batch simulations have not been fully tested yet.

## Model Description

*(Add your model description here.)*

## To-Do List

- [ ] Check that the temperature for the simulation is correct (`cfg.hParams = {'celsius': 23, 'v_init': -80}` in `src/cfg.py`). That shouldn't change anything unless any of the ion channels have a temperature-dependent behavior
- [ ] Add instructions on how to run NEURON with `mpi` to speed up sims. Ask Chris since it depends on each cluster configuration.
- [ ] Correct the weight distributions, change the probability of connection adjusting the connection rules `cfg.E2IProb, cfg.I2EProb, cfg.I2IProbchem` to reproduce experimental findings
- [ ] Update density and connectivity parameters based on Bjerke et al. 2020 and Fernandez et al. 2022
- [x] Completed Item
