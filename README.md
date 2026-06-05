Developed May to September 2020

# microCircuit

> Originally developed May–September 2020. Uploaded to GitHub in 2026.

Biophysical model of a thalamic microcircuit consisting of thalamocortical relay (TC) and local interneuron (IN) cell populations, driven by whisker stimulus-derived Poisson spike trains. Built with [Blue Brain Project](https://www.epfl.ch/research/domains/bluebrain/) tools and simulated in [NEURON](https://neuron.yale.edu/).

The core experiment tests two network arrangements — **differential** (inputs routed to different INs) vs. **same-opponent** (inputs to the same IN) — and includes a set of ablation variants to isolate the contribution of specific synaptic motifs.

---

## Network Architecture

Each subnetwork contains:
- **1 IN** (inhibitory interneuron) — bAC morpho-electrical type, dendro-dendritic output synapses
- **4 TC cells** (excitatory thalamocortical relay) — positioned radially around the IN
- **External Br inputs** — two whisker-tuned Poisson spike train sources (Branch 1, Branch 2)

Two subnetworks are coupled via IN–IN synapses (axonal and dendritic), creating either:
- `diffrun` — Br1 → subnetwork 1 IN, Br2 → subnetwork 2 IN (differential arrangement)
- `samerun` — Br1 + Br2 both → same IN (same-opponent arrangement)

### Ablation Variants

Seven circuit variants test the role of specific synaptic motifs:

| Script | Removed motif |
|--------|---------------|
| `diffsamerun_no_triadic.py` | Triadic (dendro-dendritic) synapses |
| `diffsamerun_no_axonal.py` | Axonal synapses |
| `diffsamerun_no_distriadic.py` | Distal triadic synapses |
| `diffsamerun_no_disaxonal.py` | Distal axonal synapses |
| `diffsamerun_no_triadic_no_axonal.py` | Both triadic and axonal |
| `diffsamerun_no_distriadic_no_disaxonal.py` | Both distal triadic and distal axonal |
| `diffsamerun_no_inhib.py` | All inhibition |

---

## Ion Channel Mechanisms

Custom NMODL mechanisms in `memodels/mechanisms/` and `memodels/dAD_ltb/mechanisms/`:

| Mechanism | Description |
|-----------|-------------|
| `TC_HH` | Hodgkin-Huxley Na/K channels |
| `TC_ITGHK_Des98` | T-type Ca²⁺ current (Destexhe 1998) |
| `TC_Ih_Bud97` | Hyperpolarization-activated Ih current (Budde 1997) |
| `TC_Nap_Et2` | Persistent Na⁺ current |
| `TC_iA` | A-type K⁺ current |
| `TC_iL` | L-type Ca²⁺ current |
| `TC_cadecay` | Ca²⁺ decay dynamics |
| `SK_E2` | Ca²⁺-activated K⁺ channel |

---

## Simulation Parameters

| Parameter | Value |
|-----------|-------|
| Simulator | NEURON |
| Duration | 1000 ms |
| Temperature | 34°C |
| Resting potential | −78 mV |
| Spike threshold | −29 mV (cell) / −34 mV (release) |

---

## Repository Structure

```
scripts/                          # Python simulation scripts
├── cell_200619.py                # IN and TC cell model classes
├── subnetwork_200619.py          # Subnetwork construction (TC + IN + synapses)
├── poisson_whisker.py            # Whisker stimulus Poisson spike train generator
├── convert_inputs.py             # Input format conversion
├── runall_200719.py              # Main simulation runner (diff + same arrangements)
├── diffsamerun_*.py              # Ablation variant scripts (7 total)
├── subnetwork_200619_*.py        # Subnetwork variants for each ablation
├── *_stp.py                      # Short-term plasticity variants
├── *_heiberg.py                  # Heiberg model variants
└── stgen.py                      # Vendored NeuroTools spike train generator

memodels/
├── mechanisms/                   # NMODL ion channel definitions (*.mod)
│   └── history/                  # Archived earlier script versions
├── dAD_ltb/                      # BluePyOpt single-cell optimization pipeline
│   ├── opt_model.py              # Optimization runner (genetic algorithm)
│   ├── finals.py                 # Post-optimization release finalization
│   ├── pick_features.py          # Feature selection and weighting
│   ├── analyse.py                # Post-run analysis
│   ├── setup/                    # Evaluator, protocols, templates
│   ├── config/                   # Features, params, protocols, recipes (JSON)
│   └── mechanisms/               # Ion channel .mod files for optimization
├── morphologies_IN_bAC/          # Interneuron Neurolucida morphologies (.asc)
├── morphologies_TC_dAD_ltb/      # TC relay cell morphologies (.asc)
├── morphologies_TC_dNAD_ltb/
└── bAC_IN_legacy.hoc             # NEURON HOC cell template

assets/                           # Supplementary files
├── README.pdf                    # PDF version of documentation
└── interneurondendrites*.txt     # Interneuron dendrite reference data

spike_input/                      # Poisson spike train inputs (gitignored, placeholder only)
spike_output/                     # Simulation outputs (gitignored, placeholder only)
figures/                          # Output figures (gitignored, placeholder only)
```

---

## Requirements

- Python 3.7+
- [NEURON](https://neuron.yale.edu/neuron/download) (with Python interface)
- [NeuroTools](https://github.com/NeuralEnsemble/NeuroTools) (for spike train generation)
- [BluePyOpt](https://github.com/BlueBrain/BluePyOpt) (for `dAD_ltb` optimization only)
- numpy, matplotlib

---

## Usage

### 1. Compile mechanisms

```bash
cd memodels/mechanisms
nrnivmodl
```

### 2. Generate whisker stimulus inputs

```bash
cd scripts
python poisson_whisker.py
```

Generates Poisson spike trains for 32 stimulus diameters and saves them to `spike_input/`.

### 3. Run simulation

```bash
cd scripts

# Main differential vs. same-opponent experiment
python runall_200719.py

# Or run an ablation variant, e.g.:
python diffsamerun_no_triadic.py
```

Outputs saved to `spike_output/`.

### 4. Single-cell optimization (optional)

```bash
cd memodels/dAD_ltb
python opt_model.py --etype dAD_ltb
```
