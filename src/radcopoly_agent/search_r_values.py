from radcopoly_agent.copoldb_client import CoPolDBClient


def main():
    client = CoPolDBClient()

    monomer1 = input("Monomer 1 name or SMILES: ")
    monomer2 = input("Monomer 2 name or SMILES: ")

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