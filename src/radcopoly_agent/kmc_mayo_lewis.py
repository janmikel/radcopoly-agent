"""
kmc_mayo_lewis.py

Stochastic binary copolymerization simulator based on the Mayo-Lewis
terminal model.

This module generates polymer chains as strings of "1" and "2", where:
    "1" = monomer 1
    "2" = monomer 2

The simulator currently supports two chain-length modes:

1. Fixed-DP mode
   Each chain is grown to exactly ``target_dp`` repeat units.
   This is useful for sequence-statistics studies but gives dispersity
   close to 1 because all chains have nearly the same length.

2. Stochastic-termination mode
   Each chain has a constant probability ``p_terminate`` of stopping
   after each propagation event. This creates an approximately geometric
   chain-length distribution and gives dispersities near the ideal
   free-radical polymerization limit under simple assumptions.

Scientific assumptions
----------------------
- Binary copolymerization only.
- Constant feed composition.
- Terminal Mayo-Lewis propagation probabilities.
- Optional stochastic termination.
- No monomer depletion.
- No composition drift.
- No chain transfer.
- No initiator kinetics.
- No explicit reaction clock.
- Not a full Gillespie SSA model.

The goal of this module is to provide a simple, interpretable first-stage
KMC-like simulator for exploring how reactivity ratios affect sequence
statistics, block lengths, composition, and approximate molecular-weight
distributions.
"""

from dataclasses import dataclass
from collections import Counter
import random


@dataclass
class KMCResult:
    """Container for simulation outputs.

    Attributes
    ----------
    chains
        List of generated chain sequences. Each sequence is a string of
        "1" and "2" characters.
    mn
        Number-average molecular weight computed from chain masses.
    mw
        Weight-average molecular weight computed from chain masses.
    dispersity
        Molecular-weight dispersity, Mw / Mn.
    mean_dp
        Mean degree of polymerization across generated chains.
    fraction_m1
        Overall fraction of monomer 1 units across all chains.
    fraction_m2
        Overall fraction of monomer 2 units across all chains.
    avg_block_m1
        Average contiguous block length of monomer 1.
    avg_block_m2
        Average contiguous block length of monomer 2.
    """

    chains: list[str]
    mn: float
    mw: float
    dispersity: float
    mean_dp: float
    fraction_m1: float
    fraction_m2: float
    avg_block_m1: float
    avg_block_m2: float


def choose_next_monomer(last, f1, r1, r2):
    """Choose the next monomer using Mayo-Lewis terminal probabilities.

    Parameters
    ----------
    last
        Last unit in the growing chain. Use "1", "2", or None for an
        uninitialized chain.
    f1
        Feed mole fraction of monomer 1.
    r1
        Reactivity ratio of monomer 1, r1 = k11 / k12.
    r2
        Reactivity ratio of monomer 2, r2 = k22 / k21.

    Returns
    -------
    str
        "1" if monomer 1 is added, otherwise "2".
    """

    f2 = 1.0 - f1

    if last == "1":
        p_add_1 = (r1 * f1) / (r1 * f1 + f2)
    elif last == "2":
        p_add_1 = f1 / (f1 + r2 * f2)
    else:
        # For the first unit, there is no terminal radical identity yet.
        # We initialize from the feed composition.
        p_add_1 = f1

    return "1" if random.random() < p_add_1 else "2"


def simulate_chain(length, f1, r1, r2):
    """Generate one fixed-length copolymer chain.

    Parameters
    ----------
    length
        Target degree of polymerization.
    f1, r1, r2
        Feed fraction and reactivity ratios passed to
        :func:`choose_next_monomer`.

    Returns
    -------
    str
        Chain sequence encoded as a string of "1" and "2".
    """

    chain = []
    last = None

    for _ in range(length):
        monomer = choose_next_monomer(last, f1, r1, r2)
        chain.append(monomer)
        last = monomer

    return "".join(chain)


def simulate_chain_with_termination(max_dp, f1, r1, r2, p_terminate):
    """Generate one chain with stochastic termination.

    Parameters
    ----------
    max_dp
        Safety cap on the maximum chain length. This prevents runaway
        chains if ``p_terminate`` is very small.
    f1, r1, r2
        Feed fraction and reactivity ratios passed to
        :func:`choose_next_monomer`.
    p_terminate
        Probability that the chain terminates after each propagation step.

    Returns
    -------
    str
        Chain sequence encoded as a string of "1" and "2".
    """

    chain = []
    last = None

    while len(chain) < max_dp:
        monomer = choose_next_monomer(last, f1, r1, r2)
        chain.append(monomer)
        last = monomer

        if random.random() < p_terminate:
            break

    return "".join(chain)


def block_lengths(chain):
    """Return contiguous block lengths for monomer 1 and monomer 2.

    Parameters
    ----------
    chain
        Chain sequence encoded as a string of "1" and "2".

    Returns
    -------
    dict
        Dictionary with keys "1" and "2". Each value is a list of block
        lengths for that monomer type.
    """

    if not chain:
        return {"1": [], "2": []}

    blocks = {"1": [], "2": []}
    current = chain[0]
    count = 1

    for monomer in chain[1:]:
        if monomer == current:
            count += 1
        else:
            blocks[current].append(count)
            current = monomer
            count = 1

    blocks[current].append(count)
    return blocks


def chain_mass(chain, mw1, mw2):
    """Compute chain molecular weight from repeat-unit molecular weights.

    Notes
    -----
    This sums repeat-unit masses and does not yet add initiator fragments,
    end groups, or mass changes due to specific polymerization chemistry.
    """

    return sum(mw1 if unit == "1" else mw2 for unit in chain)


def simulate_copolymerization(
    n_chains=1000,
    target_dp=100,
    f1=0.5,
    r1=1.0,
    r2=1.0,
    mw1=None,
    mw2=None,
    seed=123,
    p_terminate=None,
    max_dp=10000,
):
    """Simulate a population of binary copolymer chains.

    Parameters
    ----------
    n_chains
        Number of polymer chains to generate.
    target_dp
        Fixed chain length used when ``p_terminate`` is None.
    f1
        Feed mole fraction of monomer 1.
    r1
        Reactivity ratio of monomer 1, r1 = k11 / k12.
    r2
        Reactivity ratio of monomer 2, r2 = k22 / k21.
    mw1
        Molecular weight of monomer/repeat unit 1.
    mw2
        Molecular weight of monomer/repeat unit 2.
    seed
        Random seed for reproducibility.
    p_terminate
        Optional probability of stopping after each propagation event.
        If None, fixed-DP chains are generated.
    max_dp
        Safety maximum DP used only in stochastic-termination mode.

    Returns
    -------
    KMCResult
        Dataclass containing chain sequences and summary statistics.
    """

    random.seed(seed)

    if mw1 is None or mw2 is None:
        raise ValueError("mw1 and mw2 must be provided")

    if p_terminate is None:
        chains = [
            simulate_chain(length=target_dp, f1=f1, r1=r1, r2=r2)
            for _ in range(n_chains)
        ]
    else:
        chains = [
            simulate_chain_with_termination(
                max_dp=max_dp,
                f1=f1,
                r1=r1,
                r2=r2,
                p_terminate=p_terminate,
            )
            for _ in range(n_chains)
        ]

    dps = [len(chain) for chain in chains]
    mean_dp = sum(dps) / len(dps)

    masses = [chain_mass(chain, mw1, mw2) for chain in chains]

    mn = sum(masses) / len(masses)
    mw = sum(m * m for m in masses) / sum(masses)
    dispersity = mw / mn

    counts = Counter("".join(chains))
    total_units = counts["1"] + counts["2"]

    fraction_m1 = counts["1"] / total_units
    fraction_m2 = counts["2"] / total_units

    all_blocks_1 = []
    all_blocks_2 = []

    for chain in chains:
        blocks = block_lengths(chain)
        all_blocks_1.extend(blocks["1"])
        all_blocks_2.extend(blocks["2"])

    avg_block_m1 = sum(all_blocks_1) / len(all_blocks_1) if all_blocks_1 else 0.0
    avg_block_m2 = sum(all_blocks_2) / len(all_blocks_2) if all_blocks_2 else 0.0

    return KMCResult(
        chains=chains,
        mn=mn,
        mw=mw,
        dispersity=dispersity,
        mean_dp=mean_dp,
        fraction_m1=fraction_m1,
        fraction_m2=fraction_m2,
        avg_block_m1=avg_block_m1,
        avg_block_m2=avg_block_m2,
    )


if __name__ == "__main__":
    # Example run using Styrene / Methoxymethyl methacrylate values.
    # These values are intended only as a smoke test for the module.
    result = simulate_copolymerization(
        n_chains=1000,
        target_dp=100,
        f1=0.5,
        r1=0.395,
        r2=0.586,
        mw1=104.152,
        mw2=130.143,
        seed=123,
        p_terminate=0.01,
        max_dp=10000,
    )

    print(f"Mn: {result.mn:.2f}")
    print(f"Mw: {result.mw:.2f}")
    print(f"Dispersity: {result.dispersity:.3f}")
    print(f"Mean DP: {result.mean_dp:.2f}")
    print(f"Fraction M1: {result.fraction_m1:.3f}")
    print(f"Fraction M2: {result.fraction_m2:.3f}")
    print(f"Average M1 block length: {result.avg_block_m1:.3f}")
    print(f"Average M2 block length: {result.avg_block_m2:.3f}")
