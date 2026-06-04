# microCircuit

A biophysical microcircuit model of thalamocortical (TC) and interneuron (IN) cell populations, built with [Blue Brain Project](https://www.epfl.ch/research/domains/bluebrain/) tools and simulated in [NEURON](https://neuron.yale.edu/).

## Overview

The network consists of:
- **TC neurons** (excitatory) — thalamocortical relay cells with detailed biophysical morphology
- **IN neurons** (inhibitory) — interneurons providing feedback inhibition to TC cells
- **Virtual input population** — 25 Poisson spike-train sources driving both TC and IN populations

Connectivity is distance-dependent (Euclidean separation threshold). Synaptic dynamics use `exp2syn` (dual-exponential) with GABA (IN→TC) and AMPA (input→TC, input→IN) kinetics.

## Ion Channel Mechanisms

Custom NMODL mechanisms compiled for NEURON (in `raw_model/`):

| Mechanism | Description |
|-----------|-------------|
| `TC_HH` | Hodgkin-Huxley Na/K channels |
| `TC_iT_Des98` / `TC_ITGHK_Des98` | T-type Ca²⁺ current (Destexhe 1998) |
| `TC_Ih_Bud97` | Hyperpolarization-activated (Ih) current (Budde 1997) |
| `TC_Nap_Et2` | Persistent Na⁺ current |
| `TC_iA` | A-type K⁺ current |
| `TC_iL` | L-type Ca²⁺ current |
| `TC_cadecay` | Ca²⁺ decay dynamics |
| `SK_E2` | SK-type Ca²⁺-activated K⁺ channel |

## Simulation Parameters

| Parameter | Value |
|-----------|-------|
| Simulator | NEURON (via BMTK BioNet) |
| Duration (`tstop`) | 3000 ms |
| Time step (`dt`) | 0.1 ms |
| Temperature | 34°C |
| Resting potential | −80 mV |
| Spike threshold | −15 mV |

## Repository Structure

```
microCircuit/
├── Interneuron_model.ipynb     # Main notebook: build network, generate inputs, run simulation
├── config.json                 # BMTK simulation config (SONATA format)
├── network/                    # Generated network files (SONATA H5 + CSV)
│   ├── recurrent_network/      # TC ↔ IN recurrent connections
│   └── source_input/           # Virtual Poisson input population
├── components/                 # Cell models, morphologies, synaptic params
├── raw_model/                  # Legacy/raw NEURON files and compiled mechanisms
│   ├── *.mod                   # NMODL mechanism source files
│   ├── memodels/               # Morpho-electrical models and .asc morphology files
│   └── history/                # Prior notebook versions
└── output/                     # Simulation outputs (spikes.h5, cell_vars.h5)
```

## Requirements

- Python 3.7+
- [BMTK](https://github.com/AllenInstitute/bmtk) (`pip install bmtk`)
- [NEURON](https://neuron.yale.edu/neuron/download) (with Python interface)
- numpy, pandas, matplotlib, h5py

## Usage

Open and run `Interneuron_model.ipynb` in order:

1. **Build the network** — generates `network/recurrent_network/` and `network/source_input/` SONATA files
2. **Generate Poisson inputs** — creates `network/source_input/poission_input_spk_train.h5`
3. **Run simulation** — executes BMTK BioNet with `config.json`; outputs written to `output/`
4. **Analyze results** — spike raster and membrane voltage plots

> **Note:** NEURON mechanisms in `raw_model/` must be compiled with `nrnivmodl` before running the simulation.
