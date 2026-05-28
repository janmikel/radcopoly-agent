from radcopoly_agent.copoldb_client import CoPolDBClient


def main():
    client = CoPolDBClient()

    mode = input("Search by name or SMILES? [name/smiles]: ").strip().casefold()

    if mode == "smiles":
        smiles1 = input("Monomer 1 SMILES: ").strip()
        smiles2 = input("Monomer 2 SMILES: ").strip()
        results = client.search_candidates_by_smiles(smiles1, smiles2)
        print("WARNING: CoPolDB SMILES search is candidate/fuzzy search, not exact matching.")
        print("Use name search or RDKit canonical matching before trusting r1/r2.")

    else:
        monomer1 = input("Monomer 1 name: ").strip()
        monomer2 = input("Monomer 2 name: ").strip()
        results = client.search_by_names(monomer1, monomer2)

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


if __name__ == "__main__":
    main()