# radcopoly-agent

Tools for radical copolymerization informatics and simulation.

## Overview

`radcopoly-agent` combines literature-derived copolymerization parameters from CoPolDB with stochastic copolymerization simulations.

Current workflow:

```text
Monomer names or SMILES
        ↓
CoPolDB lookup
        ↓
r1, r2 extraction
        ↓
Mayo-Lewis stochastic simulation
        ↓
Mn, Mw, dispersity
sequence statistics
block-length statistics
plots and CSV export
```

## Current capabilities

* CoPolDB lookup by monomer name
* RDKit canonical SMILES matching
* Local caching of CoPolDB pages
* Molecular-weight extraction
* Mayo-Lewis terminal-model propagation
* Optional stochastic termination
* Mn, Mw, dispersity calculations
* DP and molecular-weight histograms
* CSV export of simulation results

## Installation

```bash
conda create -n radcopoly python=3.11 -y
conda activate radcopoly
conda install -c conda-forge rdkit -y
pip install -e .
```

## Running the workflow

```bash
python -m radcopoly_agent.query_and_simulate
```

## Scientific assumptions

Current model assumptions:

* Binary copolymerization
* Constant feed composition
* Mayo-Lewis terminal model
* Optional stochastic termination probability
* No monomer depletion
* No chain transfer
* No initiator kinetics

See `docs/science_notes.md` for details.

## Roadmap

### Near term

* Composition drift
* Monomer depletion
* Improved caching and local datasets
* Better plotting and analysis tools

### Long term

* Full Gillespie SSA implementation
* ORCA integration
* MLIP workflows
* HPC execution on Andes and Frontier
* Autonomous parameter-estimation workflows

```
```
