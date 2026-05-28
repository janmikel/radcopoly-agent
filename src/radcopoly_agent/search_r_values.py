from radcopoly_agent.copoldb_client import CoPolDBClient


def print_results(results):
    if not results:
        print("No CoPolDB results found.")
        return

    for result in results:
        print()
        print(f"Monomer 1: {result.monomer1}")
        print(f"Monomer 2: {result.monomer2}")
        print(f"r1: {result.r1}")
        print(f"r2: {result.r2}")
        print(f"DOI: {result.doi}")
        print(f"Source: {result.source_url}")

        if result.monomer1_smiles:
            print(f"Monomer 1 SMILES: {result.monomer1_smiles}")
        if result.monomer2_smiles:
            print(f"Monomer 2 SMILES: {result.monomer2_smiles}")


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

    print_results(results)


if __name__ == "__main__":
    main()