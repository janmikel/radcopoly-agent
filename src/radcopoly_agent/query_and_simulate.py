"""
query_and_simulate.py

Interactive end-to-end workflow for radcopoly-agent.

This script connects the main pieces of the project:

1. Query CoPolDB by monomer name or SMILES.
2. Select a reactivity-ratio result.
3. Retrieve monomer molecular weights.
4. Run the Mayo-Lewis stochastic copolymerization simulator.
5. Print simulation statistics.
6. Optionally export CSV files and plots.

This is intended as a command-line research workflow rather than a library API.
For reusable functions, import from:

- radcopoly_agent.copoldb_client
- radcopoly_agent.kmc_mayo_lewis
- radcopoly_agent.analysis

Example
-------
Run from the repository root:

    python -m radcopoly_agent.query_and_simulate
"""

from pathlib import Path

from radcopoly_agent.analysis import (
    export_chains_csv,
    export_summary_csv,
    plot_dp_distribution,
    plot_molecular_weight_distribution,
)
from radcopoly_agent.copoldb_client import CoPolDBClient
from radcopoly_agent.kmc_mayo_lewis import simulate_copolymerization


def choose_result(results):
    """Choose one CoPolDB result from a list of candidates.

    Parameters
    ----------
    results
        List of ReactivityRatioResult objects.

    Returns
    -------
    ReactivityRatioResult or None
        The selected result. Returns None if the result list is empty.

    Notes
    -----
    CoPolDB searches can return multiple literature measurements or related
    monomer pairs. The workflow asks the user to choose which result should be
    used for the simulation.
    """

    if not results:
        return None

    if len(results) == 1:
        return results[0]

    print("\nMultiple CoPolDB results found:")
    for i, result in enumerate(results, start=1):
        print(
            f"{i}. {result.monomer1} / {result.monomer2} "
            f"r1={result.r1}, r2={result.r2}, DOI={result.doi}"
        )

    choice = input("\nChoose result number: ").strip()
    index = int(choice) - 1
    return results[index]


def ensure_monomer_molecular_weights(client, result):
    """Populate molecular weights on a CoPolDB result if they are missing.

    Parameters
    ----------
    client
        CoPolDBClient instance.
    result
        ReactivityRatioResult selected for simulation.

    Returns
    -------
    ReactivityRatioResult
        The same result object, possibly updated with monomer molecular weights.
    """

    if result.monomer1_mw is None and result.monomer1_url:
        result.monomer1_mw = client.extract_molecular_weight_from_monomer_page(
            result.monomer1_url
        )

    if result.monomer2_mw is None and result.monomer2_url:
        result.monomer2_mw = client.extract_molecular_weight_from_monomer_page(
            result.monomer2_url
        )

    return result


def main():
    """Run the interactive CoPolDB lookup and simulation workflow."""

    client = CoPolDBClient()

    mode = input("Search by name or SMILES? [name/smiles]: ").strip().casefold()

    if mode == "smiles":
        smiles1 = input("Monomer 1 SMILES: ").strip()
        smiles2 = input("Monomer 2 SMILES: ").strip()
        results = client.search_by_smiles_exact(smiles1, smiles2)
    else:
        monomer1 = input("Monomer 1 name: ").strip()
        monomer2 = input("Monomer 2 name: ").strip()
        results = client.search_by_names(monomer1, monomer2)

    result = choose_result(results)

    if result is None:
        print("No CoPolDB result found. Cannot run KMC.")
        return

    result = ensure_monomer_molecular_weights(client, result)

    print("\nSelected CoPolDB result:")
    print(f"Monomer 1: {result.monomer1}")
    print(f"Monomer 2: {result.monomer2}")
    print(f"r1: {result.r1}")
    print(f"r2: {result.r2}")
    print(f"DOI: {result.doi}")
    print(f"Monomer 1 MW: {result.monomer1_mw}")
    print(f"Monomer 2 MW: {result.monomer2_mw}")

    if result.monomer1_mw is None or result.monomer2_mw is None:
        print("Could not determine monomer molecular weights. Cannot run KMC.")
        return

    f1 = float(input("\nFeed fraction f1 for monomer 1 [0.5]: ") or 0.5)
    n_chains = int(input("Number of chains [1000]: ") or 1000)
    target_dp = int(input("Target DP [100]: ") or 100)

    p_term_text = input(
        "Termination probability per propagation step [blank = fixed DP]: "
    ).strip()
    p_terminate = float(p_term_text) if p_term_text else None

    if p_terminate is not None:
        max_dp = int(input("Safety max DP [10000]: ") or 10000)
    else:
        max_dp = target_dp

    seed = int(input("Random seed [123]: ") or 123)

    sim = simulate_copolymerization(
        n_chains=n_chains,
        target_dp=target_dp,
        f1=f1,
        r1=result.r1,
        r2=result.r2,
        mw1=result.monomer1_mw,
        mw2=result.monomer2_mw,
        seed=seed,
        p_terminate=p_terminate,
        max_dp=max_dp,
    )

    print("\nKMC Mayo-Lewis simulation results:")
    print(f"Mn: {sim.mn:.2f}")
    print(f"Mw: {sim.mw:.2f}")
    print(f"Dispersity: {sim.dispersity:.3f}")
    print(f"Mean DP: {sim.mean_dp:.2f}")
    print(f"Fraction M1: {sim.fraction_m1:.3f}")
    print(f"Fraction M2: {sim.fraction_m2:.3f}")
    print(f"Average M1 block length: {sim.avg_block_m1:.3f}")
    print(f"Average M2 block length: {sim.avg_block_m2:.3f}")

    print("\nExample chains:")
    for chain in sim.chains[:5]:
        print(chain)

    save_outputs = input("\nSave CSV and plot outputs? [y/n]: ").strip().casefold()

    if save_outputs == "y":
        output_dir = Path("outputs")
        output_dir.mkdir(exist_ok=True)

        export_summary_csv(
            output_dir / "simulation_summary.csv",
            sim=sim,
            copoldb_result=result,
            f1=f1,
            n_chains=n_chains,
            target_dp=target_dp,
            max_dp=max_dp,
            p_terminate=p_terminate,
        )

        export_chains_csv(
            output_dir / "chains.csv",
            chains=sim.chains,
            mw1=result.monomer1_mw,
            mw2=result.monomer2_mw,
        )

        plot_molecular_weight_distribution(
            output_dir / "molecular_weight_distribution.png",
            chains=sim.chains,
            mw1=result.monomer1_mw,
            mw2=result.monomer2_mw,
        )

        plot_dp_distribution(
            output_dir / "dp_distribution.png",
            chains=sim.chains,
        )

        print("\nSaved outputs to:")
        print(output_dir.resolve())


if __name__ == "__main__":
    main()
