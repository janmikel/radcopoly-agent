# Science notes

This document explains the scientific assumptions behind `radcopoly-agent`.

## Current model

The current simulator is a stochastic Mayo-Lewis copolymerization model. It uses terminal-model propagation probabilities derived from monomer reactivity ratios `r1` and `r2` and a constant feed fraction `f1`.

The model is useful for sequence-level simulation, copolymer composition estimates, block-length statistics, and approximate molecular-weight distributions when stochastic termination is enabled.

## Propagation model

For a chain ending in monomer 1, the probability of adding monomer 1 is approximated as:

```text
P(add M1 | chain ends in M1) = r1 f1 / (r1 f1 + f2)
```

For a chain ending in monomer 2:

```text
P(add M1 | chain ends in M2) = f1 / (f1 + r2 f2)
```

where `f2 = 1 - f1`.

## Stochastic termination

When `p_terminate` is provided, each chain has a constant probability of stopping after each propagation event. This creates a geometric chain-length distribution and gives dispersities near the ideal free-radical limit under simple assumptions.

This is not yet a full Gillespie stochastic simulation algorithm. It does not explicitly model initiator decomposition, radical concentration, bimolecular termination, chain transfer, or monomer depletion.

## Current limitations

- Binary copolymerization only.
- Constant feed composition.
- No monomer depletion or composition drift.
- No chain transfer.
- No explicit initiator kinetics.
- No explicit reaction clock.
- No full Gillespie SSA yet.

## Planned scientific upgrades

1. Composition drift and monomer depletion.
2. Full Gillespie SSA with explicit propagation and termination propensities.
3. Chain transfer events.
4. Quantum-chemistry fallback for missing reactivity ratios.
5. Integration with ORCA, MLIP, and HPC workflows.
