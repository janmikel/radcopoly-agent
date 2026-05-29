from dataclasses import dataclass
from typing import Optional
from urllib.parse import urljoin

import hashlib
import json
import re
import time
from pathlib import Path

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

    def __init__(
        self,
        verify_ssl: bool = False,
        cache_dir: str | Path = "data/cache/copoldb",
        use_cache: bool = True,
        sleep_seconds: float = 0.2,
    ):
        self.verify_ssl = verify_ssl
        self.cache_dir = Path(cache_dir)
        self.use_cache = use_cache
        self.sleep_seconds = sleep_seconds

        self.search_cache_dir = self.cache_dir / "search_pages"
        self.monomer_cache_dir = self.cache_dir / "monomer_pages"

        self.search_cache_dir.mkdir(parents=True, exist_ok=True)
        self.monomer_cache_dir.mkdir(parents=True, exist_ok=True)

    def canonicalize_smiles(self, smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return Chem.MolToSmiles(mol)

    def _same_name(self, a: str, b: str) -> bool:
        return a.strip().casefold() == b.strip().casefold()

    def _cache_key(self, url: str, params: Optional[dict] = None) -> str:
        payload = {
            "url": url,
            "params": params or {},
        }
        text = json.dumps(payload, sort_keys=True)
        return hashlib.sha256(text.encode("utf-8")).hexdigest()

    def _cache_path(self, cache_dir: Path, url: str, params: Optional[dict] = None) -> Path:
        return cache_dir / f"{self._cache_key(url, params)}.html"

    def _metadata_path(self, html_path: Path) -> Path:
        return html_path.with_suffix(".json")

    def _get_html(
        self,
        url: str,
        params: Optional[dict] = None,
        timeout: int = 30,
        cache_dir: Optional[Path] = None,
    ) -> tuple[str, str]:
        """
        Fetch HTML from CoPolDB with a small local cache.

        Returns:
            (html_text, final_url)
        """

        cache_dir = cache_dir or self.cache_dir
        html_path = self._cache_path(cache_dir, url, params)
        metadata_path = self._metadata_path(html_path)

        if self.use_cache and html_path.exists() and metadata_path.exists():
            html = html_path.read_text(encoding="utf-8")
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            return html, metadata.get("final_url", url)

        response = requests.get(
            url,
            params=params,
            timeout=timeout,
            verify=self.verify_ssl,
        )
        response.raise_for_status()

        html = response.text

        if self.use_cache:
            html_path.write_text(html, encoding="utf-8")
            metadata_path.write_text(
                json.dumps(
                    {
                        "url": url,
                        "params": params or {},
                        "final_url": response.url,
                        "fetched_at_unix": time.time(),
                    },
                    indent=2,
                    sort_keys=True,
                ),
                encoding="utf-8",
            )

        if self.sleep_seconds:
            time.sleep(self.sleep_seconds)

        return html, response.url

    def extract_smiles_from_monomer_page(self, monomer_url: str) -> Optional[str]:
        html, _ = self._get_html(
            monomer_url,
            timeout=10,
            cache_dir=self.monomer_cache_dir,
        )

        soup = BeautifulSoup(html, "html.parser")
        text = soup.get_text("\n", strip=True)
        lines = [line.strip() for line in text.splitlines() if line.strip()]

        for i, line in enumerate(lines):
            if line.casefold() == "smiles" and i + 1 < len(lines):
                return lines[i + 1]

        return None

    def extract_molecular_weight_from_monomer_page(
        self,
        monomer_url: str,
    ) -> Optional[float]:
        html, _ = self._get_html(
            monomer_url,
            timeout=10,
            cache_dir=self.monomer_cache_dir,
        )

        soup = BeautifulSoup(html, "html.parser")
        text = soup.get_text("\n", strip=True)
        lines = [line.strip() for line in text.splitlines() if line.strip()]

        for i, line in enumerate(lines):
            if "molecular weight" in line.casefold():
                for j in range(i + 1, min(i + 5, len(lines))):
                    match = re.search(r"[-+]?\d*\.?\d+", lines[j])
                    if match:
                        return float(match.group())

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

    def _populate_monomer_details(self, result: ReactivityRatioResult) -> ReactivityRatioResult:
        if result.monomer1_url:
            result.monomer1_smiles = self.extract_smiles_from_monomer_page(
                result.monomer1_url
            )
            result.monomer1_mw = self.extract_molecular_weight_from_monomer_page(
                result.monomer1_url
            )

        if result.monomer2_url:
            result.monomer2_smiles = self.extract_smiles_from_monomer_page(
                result.monomer2_url
            )
            result.monomer2_mw = self.extract_molecular_weight_from_monomer_page(
                result.monomer2_url
            )

        return result

    def search_by_names(
        self,
        monomer1: str,
        monomer2: str,
        exact: bool = True,
        populate_details: bool = False,
    ):
        url = f"{self.base_url}/copolym/list"
        params = {
            "m1": monomer1,
            "m2": monomer2,
            "mode": 0,
        }

        html, final_url = self._get_html(
            url,
            params=params,
            timeout=30,
            cache_dir=self.search_cache_dir,
        )

        results = self._parse_results(html, final_url)

        if exact:
            results = [
                result
                for result in results
                if self._same_name(result.monomer1, monomer1)
                and self._same_name(result.monomer2, monomer2)
            ]

        if populate_details:
            results = [self._populate_monomer_details(result) for result in results]

        return results

    def search_candidates_by_smiles(self, smiles1: str, smiles2: str):
        url = f"{self.base_url}/copolym/list"
        params = {
            "sm1": smiles1,
            "sm2": smiles2,
            "mode": 0,
        }

        html, final_url = self._get_html(
            url,
            params=params,
            timeout=30,
            cache_dir=self.search_cache_dir,
        )

        return self._parse_results(html, final_url)

    def search_by_smiles_exact(
        self,
        smiles1: str,
        smiles2: str,
        max_candidates: int = 100,
    ):
        target1 = self.canonicalize_smiles(smiles1)
        target2 = self.canonicalize_smiles(smiles2)

        candidates = self.search_candidates_by_smiles(smiles1, smiles2)
        exact_matches = []

        n_check = min(len(candidates), max_candidates)

        for i, result in enumerate(candidates[:max_candidates], start=1):
            print(
                f"Checking candidate {i}/{n_check}: "
                f"{result.monomer1} / {result.monomer2}"
            )

            if not result.monomer1_url or not result.monomer2_url:
                continue

            result = self._populate_monomer_details(result)

            if not result.monomer1_smiles or not result.monomer2_smiles:
                continue

            try:
                db_can1 = self.canonicalize_smiles(result.monomer1_smiles)
                db_can2 = self.canonicalize_smiles(result.monomer2_smiles)
            except ValueError:
                continue

            if db_can1 == target1 and db_can2 == target2:
                exact_matches.append(result)

        return exact_matches

    def clear_cache(self):
        for path in self.cache_dir.glob("**/*"):
            if path.is_file():
                path.unlink()
