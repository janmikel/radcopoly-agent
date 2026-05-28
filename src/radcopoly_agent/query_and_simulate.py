from radcopoly_agent.copoldb_client import CoPolDBClient
from radcopoly_agent.kmc_mayo_lewis import simulate_copolymerization


def choose_result(results):
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


def main():
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

    print("\nSelected CoPolDB result:")
    print(f"Monomer 1: {result.monomer1}")
    print(f"Monomer 2: {result.monomer2}")
    print(f"r1: {result.r1}")
    print(f"r2: {result.r2}")
    print(f"DOI: {result.doi}")

    f1 = float(input("\nFeed fraction f1 for monomer 1 [0.5]: ") or 0.5)
    n_chains = int(input("Number of chains [1000]: ") or 1000)
    target_dp = int(input("Target DP [100]: ") or 100)
    seed = int(input("Random seed [123]: ") or 123)

    sim = simulate_copolymerization(
        n_chains=n_chains,
        target_dp=target_dp,
        f1=f1,
        r1=result.r1,
        r2=result.r2,
        seed=seed,
    )

    print("\nKMC Mayo-Lewis simulation results:")
    print(f"Mn: {sim.mn:.2f}")
    print(f"Mw: {sim.mw:.2f}")
    print(f"Dispersity: {sim.dispersity:.3f}")
    print(f"Fraction M1: {sim.fraction_m1:.3f}")
    print(f"Fraction M2: {sim.fraction_m2:.3f}")
    print(f"Average M1 block length: {sim.avg_block_m1:.3f}")
    print(f"Average M2 block length: {sim.avg_block_m2:.3f}")

    print("\nExample chains:")
    for chain in sim.chains[:5]:
        print(chain)


if __name__ == "__main__":
    main()