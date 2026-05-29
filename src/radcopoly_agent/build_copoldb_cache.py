import csv
import time
from pathlib import Path

from radcopoly_agent.copoldb_client import CoPolDBClient


def main():
    input_csv = Path("data/cache_seed_pairs.csv")

    if not input_csv.exists():
        print("Create data/cache_seed_pairs.csv with columns: monomer1,monomer2")
        return

    client = CoPolDBClient(
        verify_ssl=False,
        use_cache=True,
        sleep_seconds=1.0,
    )

    with input_csv.open() as f:
        reader = csv.DictReader(f)

        for row in reader:
            monomer1 = row["monomer1"].strip()
            monomer2 = row["monomer2"].strip()

            print(f"\nQuerying: {monomer1} / {monomer2}")

            try:
                results = client.search_by_names(
                    monomer1,
                    monomer2,
                    exact=False,
                    populate_details=True,
                )

                print(f"Found {len(results)} candidate results")

                for result in results:
                    print(
                        f"  {result.monomer1} / {result.monomer2} "
                        f"r1={result.r1}, r2={result.r2}"
                    )

            except Exception as exc:
                print(f"ERROR: {exc}")

            time.sleep(1.0)


if __name__ == "__main__":
    main()