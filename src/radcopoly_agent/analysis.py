"""
analysis.py

Post-processing utilities for radcopoly-agent simulations.

This module converts simulated copolymer chains into tabular data and plots.
It is intentionally lightweight: the simulator produces chain sequences, and
this module handles reporting, CSV export, and basic diagnostic figures.

Current outputs
---------------
- Per-chain degree of polymerization
- Per-chain molecular weight
- Per-chain composition
- Summary CSV files
- Molecular-weight distribution histograms
- Degree-of-polymerization distribution histograms

Chain encoding
--------------
Chains are represented as strings of "1" and "2", where:
    "1" = monomer 1
    "2" = monomer 2

The molecular-weight calculations use monomer/repeat-unit molecular weights
provided by the CoPolDB workflow.
"""

from pathlib import Path
from collections import Counter
import csv
import matplotlib.pyplot as plt


def chain_mass(chain, mw1, mw2):
    """Compute the approximate molecular weight of one chain.

    Parameters
    ----------
    chain
        Chain sequence encoded as a string of "1" and "2".
    mw1
        Molecular weight of monomer/repeat unit 1.
    mw2
        Molecular weight of monomer/repeat unit 2.

    Returns
    -------
    float
        Sum of repeat-unit masses in the chain.

    Notes
    -----
    This does not currently add initiator fragments, end groups, or mass
    corrections for specific polymerization chemistry.
    """

    return sum(mw1 if unit == "1" else mw2 for unit in chain)


def chain_composition(chain):
    """Return monomer fractions for one chain.

    Parameters
    ----------
    chain
        Chain sequence encoded as a string of "1" and "2".

    Returns
    -------
    tuple[float, float]
        Fraction of monomer 1 and fraction of monomer 2.
    """

    counts = Counter(chain)
    total = len(chain)
    return counts["1"] / total, counts["2"] / total


def collect_chain_records(chains, mw1, mw2):
    """Convert chain sequences into per-chain dictionary records.

    Parameters
    ----------
    chains
        Iterable of chain sequences.
    mw1, mw2
        Molecular weights for monomer/repeat units 1 and 2.

    Returns
    -------
    list[dict]
        Records containing chain ID, DP, mass, composition, and sequence.
    """

    records = []

    for i, chain in enumerate(chains, start=1):
        mass = chain_mass(chain, mw1, mw2)
        f1, f2 = chain_composition(chain)

        records.append(
            {
                "chain_id": i,
                "dp": len(chain),
                "mass": mass,
                "fraction_m1": f1,
                "fraction_m2": f2,
                "sequence": chain,
            }
        )

    return records


def export_summary_csv(
    path,
    sim,
    copoldb_result,
    f1,
    n_chains,
    target_dp,
    max_dp,
    p_terminate,
):
    """Write a one-column simulation summary CSV.

    Parameters
    ----------
    path
        Output CSV path.
    sim
        KMCResult object returned by ``simulate_copolymerization``.
    copoldb_result
        ReactivityRatioResult containing monomer names, r1/r2, DOI, and
        molecular weights.
    f1
        Feed fraction of monomer 1.
    n_chains
        Number of simulated chains.
    target_dp
        Fixed DP used in fixed-DP mode.
    max_dp
        Safety maximum DP used in stochastic-termination mode.
    p_terminate
        Termination probability per propagation step, or None for fixed-DP
        mode.
    """

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    rows = [
        ["monomer1", copoldb_result.monomer1],
        ["monomer2", copoldb_result.monomer2],
        ["r1", copoldb_result.r1],
        ["r2", copoldb_result.r2],
        ["doi", copoldb_result.doi],
        ["monomer1_mw", copoldb_result.monomer1_mw],
        ["monomer2_mw", copoldb_result.monomer2_mw],
        ["feed_fraction_m1", f1],
        ["n_chains", n_chains],
        ["target_dp", target_dp],
        ["max_dp", max_dp],
        ["p_terminate", p_terminate],
        ["mn", sim.mn],
        ["mw", sim.mw],
        ["dispersity", sim.dispersity],
        ["mean_dp", sim.mean_dp],
        ["fraction_m1", sim.fraction_m1],
        ["fraction_m2", sim.fraction_m2],
        ["avg_block_m1", sim.avg_block_m1],
        ["avg_block_m2", sim.avg_block_m2],
    ]

    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def export_chains_csv(path, chains, mw1, mw2):
    """Write per-chain simulation data to CSV.

    The output contains one row per chain with:
    chain ID, DP, mass, monomer composition, and raw sequence.
    """

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    records = collect_chain_records(chains, mw1, mw2)

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "chain_id",
                "dp",
                "mass",
                "fraction_m1",
                "fraction_m2",
                "sequence",
            ],
        )
        writer.writeheader()
        writer.writerows(records)


def plot_molecular_weight_distribution(path, chains, mw1, mw2, bins=50):
    """Save a molecular-weight distribution histogram.

    Parameters
    ----------
    path
        Output image path.
    chains
        Simulated chain sequences.
    mw1, mw2
        Molecular weights for monomer/repeat units 1 and 2.
    bins
        Number of histogram bins.
    """

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    masses = [chain_mass(chain, mw1, mw2) for chain in chains]

    plt.figure()
    plt.hist(masses, bins=bins)
    plt.xlabel("Molecular weight")
    plt.ylabel("Count")
    plt.title("Molecular weight distribution")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def plot_dp_distribution(path, chains, bins=50):
    """Save a degree-of-polymerization distribution histogram."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    dps = [len(chain) for chain in chains]

    plt.figure()
    plt.hist(dps, bins=bins)
    plt.xlabel("Degree of polymerization")
    plt.ylabel("Count")
    plt.title("DP distribution")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()
