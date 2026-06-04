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
| Simulator | NEURON |
| Duration (`tstop`) | 3000 ms |
| Time step (`dt`) | 0.1 ms |
| Temperature | 34°C |
| Resting potential | −80 mV |
| Spike threshold | −15 mV |

## Repository Structure

```
microCircuit/
├── config.json                         # Simulation config (SONATA format)
├── network/                            # Generated network files (SONATA H5 + CSV)
│   ├── recurrent_network/              # TC ↔ IN recurrent connections
│   └── source_input/                   # Poisson input population
├── components/                         # Cell models, morphologies, synaptic params
└── raw_model/                          # NEURON model files and compiled mechanisms
    ├── *.mod                           # NMODL mechanism source files
    └── memodels/
        ├── mechanisms/
        │   ├── cell.py                 # Cell model definitions
        │   ├── subnetwork.py           # Subnetwork construction
        │   ├── poisson.py              # Poisson spike-train generation
        │   ├── diffrun.py              # Run differential IN network arrangement
        │   └── samerun.py             # Run same IN network arrangement
        └── morphologies_*/             # Neurolucida .asc morphology files
```

## Requirements

- Python 3.7+
- [NEURON](https://neuron.yale.edu/neuron/download) (with Python interface)
- numpy, matplotlib

## Usage

1. **Compile mechanisms** — run `nrnivmodl` in `raw_model/memodels/mechanisms/` to compile the `.mod` files
2. **Generate Poisson inputs** — run `poisson.py` to create spike-train inputs
3. **Run simulation** — run `diffrun.py` (differential IN arrangement) or `samerun.py` (same IN arrangement); outputs written to `spike_output/`
