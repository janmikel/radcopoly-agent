from pathlib import Path
from collections import Counter
import csv
import matplotlib.pyplot as plt


def chain_mass(chain, mw1, mw2):
    return sum(mw1 if unit == "1" else mw2 for unit in chain)


def chain_composition(chain):
    counts = Counter(chain)
    total = len(chain)
    return counts["1"] / total, counts["2"] / total


def collect_chain_records(chains, mw1, mw2):
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


def export_summary_csv(path, sim, copoldb_result, f1, n_chains, target_dp, max_dp, p_terminate):
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