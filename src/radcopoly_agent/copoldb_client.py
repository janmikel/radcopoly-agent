from dataclasses import dataclass
from typing import Optional
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup
from rdkit import Chem

import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

@dataclass
class ReactivityRatioResult:
    monomer1: str
    monomer2: str
    r1: Optional[float]
    r2: Optional[float]
    doi: Optional[str]
    source_url: str
    monomer1_url: Optional[str] = None
    monomer2_url: Optional[str] = None
    monomer1_smiles: Optional[str] = None
    monomer2_smiles: Optional[str] = None
    monomer1_mw: Optional[float] = None
    monomer2_mw: Optional[float] = None


class CoPolDBClient:
    base_url = "https://www.copoldb.jp"

    def __init__(self, verify_ssl: bool = False):
        self.verify_ssl = verify_ssl

    def canonicalize_smiles(self, smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return Chem.MolToSmiles(mol)

    def _same_name(self, a: str, b: str) -> bool:
        return a.strip().casefold() == b.strip().casefold()

    def extract_smiles_from_monomer_page(self, monomer_url: str) -> Optional[str]:
        response = requests.get(monomer_url, timeout=30, verify=self.verify_ssl)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")
        text = soup.get_text("\n", strip=True)
        lines = [line.strip() for line in text.splitlines() if line.strip()]

        for i, line in enumerate(lines):
            if line.casefold() == "smiles" and i + 1 < len(lines):
                return lines[i + 1]

        return None

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

                monomer_links = []
                for link in row.find_all("a", href=True):
                    href = link["href"]
                    if "/monomer/detail" in href:
                        monomer_links.append(
                            (
                                link.get_text(" ", strip=True),
                                urljoin(self.base_url, href),
                            )
                        )

                if len(monomer_links) < 2:
                    continue

                results.append(
                    ReactivityRatioResult(
                        monomer1=monomer_links[0][0],
                        monomer2=monomer_links[1][0],
                        r1=numbers[0],
                        r2=numbers[1],
                        doi=doi,
                        source_url=source_url,
                        monomer1_url=monomer_links[0][1],
                        monomer2_url=monomer_links[1][1],
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

        response = requests.get(url, params=params, timeout=30, verify=self.verify_ssl)
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

        response = requests.get(url, params=params, timeout=30, verify=self.verify_ssl)
        response.raise_for_status()

        return self._parse_results(response.text, response.url)

    def search_by_smiles_exact(self, smiles1: str, smiles2: str):
        target1 = self.canonicalize_smiles(smiles1)
        target2 = self.canonicalize_smiles(smiles2)

        candidates = self.search_candidates_by_smiles(smiles1, smiles2)
        exact_matches = []

        for result in candidates:
            if not result.monomer1_url or not result.monomer2_url:
                continue

            db_smiles1 = self.extract_smiles_from_monomer_page(result.monomer1_url)
            db_smiles2 = self.extract_smiles_from_monomer_page(result.monomer2_url)

            if not db_smiles1 or not db_smiles2:
                continue

            try:
                db_can1 = self.canonicalize_smiles(db_smiles1)
                db_can2 = self.canonicalize_smiles(db_smiles2)
            except ValueError:
                continue

            result.monomer1_smiles = db_smiles1
            result.monomer2_smiles = db_smiles2

            result.monomer1_mw = (
                self.extract_molecular_weight_from_monomer_page(
                    result.monomer1_url
                )
            )

            result.monomer2_mw = (
                self.extract_molecular_weight_from_monomer_page(
                    result.monomer2_url
                )
            )

            if db_can1 == target1 and db_can2 == target2:
                exact_matches.append(result)

        return exact_matches
    
    def extract_molecular_weight_from_monomer_page(
        self,
        monomer_url: str,
    ) -> Optional[float]:

        response = requests.get(monomer_url, timeout=10, verify=self.verify_ssl)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")

        text = soup.get_text("\n", strip=True)
        lines = [line.strip() for line in text.splitlines() if line.strip()]

        for i, line in enumerate(lines):
            if "molecular weight" in line.casefold():

                for j in range(i + 1, min(i + 5, len(lines))):
                    try:
                        return float(lines[j])
                    except ValueError:
                        continue

        return None