# radcopoly-agent

Tools for radical copolymerization informatics and simulation.

## Overview

`radcopoly-agent` combines literature-derived copolymerization parameters from CoPolDB with stochastic copolymerization simulations.
Reactivity-ratio data are obtained from CoPolDB when available.

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

## Data Source

This project uses data obtained from **CoPolDB**, a database of radical copolymerization reactions and monomer reactivity ratios.

If you use `radcopoly-agent` in research, please cite both this software and the CoPolDB publication:

> T. Yamamoto, T. Nakao, H. Uehara, T. Kubo, and K. Ute, *CoPolDB: a database for radical copolymerization of vinyl monomers*, Polymer Chemistry, 2024.
>
> DOI: https://doi.org/10.1039/D3PY01372C

### CoPolDB Website

https://www.copoldb.jp/

### Disclaimer

`radcopoly-agent` is not affiliated with the CoPolDB project. CoPolDB remains the authoritative source of the copolymerization data used by this software.

## Citation

If you use this software in a publication, please cite:

### Software

Carrillo, J. M. *radcopoly-agent* (2026)

GitHub:
https://github.com/janmikel/radcopoly-agent

### Database

Yamamoto, T.; Nakao, T.; Uehara, H.; Kubo, T.; Ute, K.
*CoPolDB: a database for radical copolymerization of vinyl monomers*.
Polymer Chemistry (2024).

DOI: https://doi.org/10.1039/D3PY01372C


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
