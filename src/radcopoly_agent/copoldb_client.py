from dataclasses import dataclass
from typing import Optional
import requests
from bs4 import BeautifulSoup


@dataclass
class ReactivityRatioResult:
    monomer1: str
    monomer2: str
    r1: Optional[float]
    r2: Optional[float]
    doi: Optional[str]
    source_url: str


class CoPolDBClient:
    base_url = "https://www.copoldb.jp"

    def _same_name(self, a: str, b: str) -> bool:
        return a.strip().casefold() == b.strip().casefold()

    def _parse_results(self, html: str, source_url: str):
        soup = BeautifulSoup(html, "html.parser")
        results = []

        for table in soup.find_all("table"):
            for row in table.find_all("tr"):
                cells = row.find_all("td")
                if not cells:
                    continue

                text = row.get_text(" ", strip=True)
                parts = text.split()

                numbers = []
                for part in parts:
                    try:
                        numbers.append(float(part))
                    except ValueError:
                        pass

                if len(numbers) < 2:
                    continue

                doi = None
                for link in row.find_all("a", href=True):
                    href = link["href"]
                    if "doi.org" in href:
                        doi = href

                monomer_links = [
                    link.get_text(" ", strip=True)
                    for link in row.find_all("a", href=True)
                    if "/monomer/detail" in link["href"]
                ]

                if len(monomer_links) < 2:
                    continue

                results.append(
                    ReactivityRatioResult(
                        monomer1=monomer_links[0],
                        monomer2=monomer_links[1],
                        r1=numbers[0],
                        r2=numbers[1],
                        doi=doi,
                        source_url=source_url,
                    )
                )

        return results

    def search_by_names(self, monomer1: str, monomer2: str, exact: bool = True):
        url = f"{self.base_url}/copolym/list"
        params = {
            "m1": monomer1,
            "m2": monomer2,
            "mode": 0,
        }

        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()

        results = self._parse_results(response.text, response.url)

        if exact:
            results = [
                result
                for result in results
                if self._same_name(result.monomer1, monomer1)
                and self._same_name(result.monomer2, monomer2)
            ]

        return results

    def search_candidates_by_smiles(self, smiles1: str, smiles2: str):
        url = f"{self.base_url}/copolym/list"
        params = {
            "sm1": smiles1,
            "sm2": smiles2,
            "mode": 0,
        }

        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()

        return self._parse_results(response.text, response.url)