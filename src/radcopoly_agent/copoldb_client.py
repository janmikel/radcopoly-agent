"""
copoldb_client.py

Client interface for CoPolDB (https://www.copoldb.jp).

This module is responsible for retrieving copolymerization parameters from
CoPolDB and converting the returned HTML pages into Python data structures
that can be consumed by the rest of radcopoly-agent.

Responsibilities
----------------
- Query CoPolDB by monomer name.
- Query CoPolDB by SMILES.
- Canonicalize SMILES using RDKit.
- Extract monomer reactivity ratios r1 and r2.
- Extract DOI references.
- Extract monomer SMILES and molecular weights from monomer detail pages.
- Cache downloaded CoPolDB pages locally.

Caching
-------
Downloaded HTML pages are stored under:

    data/cache/copoldb/

This reduces repeated requests to CoPolDB, makes repeated simulations faster,
and improves reproducibility when working with the same monomer pairs.

Scientific role
---------------
This module provides the literature-derived reactivity-ratio data used by the
simulation workflows in radcopoly-agent. CoPolDB remains the authoritative data
source.

If you use CoPolDB data in research, cite:

Yamamoto et al.
CoPolDB: a database for radical copolymerization of vinyl monomers.
Polymer Chemistry (2024).

DOI: https://doi.org/10.1039/D3PY01372C
"""

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

# Temporary workaround:
# CoPolDB has periodically presented SSL certificate issues. We disable
# warnings because this client may run with verify_ssl=False when reading
# public CoPolDB pages. Do not use verify_ssl=False for private credentials
# or sensitive data.
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


@dataclass
class ReactivityRatioResult:
    """Container for one CoPolDB copolymerization result.

    Attributes
    ----------
    monomer1, monomer2
        Human-readable monomer names returned by CoPolDB.
    r1, r2
        Monomer reactivity ratios.
    doi
        Literature DOI associated with the measurement, if available.
    source_url
        CoPolDB URL used to retrieve the result.
    monomer1_url, monomer2_url
        URLs for the monomer detail pages.
    monomer1_smiles, monomer2_smiles
        SMILES strings extracted from monomer detail pages.
    monomer1_mw, monomer2_mw
        Molecular weights extracted from monomer detail pages.
    """

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
    """Small HTML client for CoPolDB.

    Parameters
    ----------
    verify_ssl
        Whether requests should verify the CoPolDB HTTPS certificate. This is
        False by default because the site has produced certificate errors in
        testing. Set True once certificate validation is reliable.
    cache_dir
        Directory where cached search pages and monomer pages are stored.
    use_cache
        If True, previously downloaded HTML pages are reused.
    sleep_seconds
        Delay after network fetches. This keeps the client polite and avoids
        sending rapid repeated requests.

    Notes
    -----
    This client intentionally does not bulk-mirror CoPolDB. It caches only the
    pages requested by user workflows.
    """

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
        """Return an RDKit canonical SMILES string.

        CoPolDB's SMILES search can return fuzzy or partial matches. Canonical
        SMILES comparison lets the client filter candidate results down to
        exact molecular-structure matches.

        Raises
        ------
        ValueError
            If RDKit cannot parse the provided SMILES.
        """

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return Chem.MolToSmiles(mol)

    def _same_name(self, a: str, b: str) -> bool:
        """Case-insensitive exact name comparison."""

        return a.strip().casefold() == b.strip().casefold()

    def _cache_key(self, url: str, params: Optional[dict] = None) -> str:
        """Build a stable SHA-256 cache key from a URL and query parameters."""

        payload = {
            "url": url,
            "params": params or {},
        }
        text = json.dumps(payload, sort_keys=True)
        return hashlib.sha256(text.encode("utf-8")).hexdigest()

    def _cache_path(self, cache_dir: Path, url: str, params: Optional[dict] = None) -> Path:
        """Return the HTML cache path for a request."""

        return cache_dir / f"{self._cache_key(url, params)}.html"

    def _metadata_path(self, html_path: Path) -> Path:
        """Return the JSON metadata path associated with a cached HTML file."""

        return html_path.with_suffix(".json")

    def _get_html(
        self,
        url: str,
        params: Optional[dict] = None,
        timeout: int = 30,
        cache_dir: Optional[Path] = None,
    ) -> tuple[str, str]:
        """Fetch HTML from CoPolDB using a small local cache.

        Parameters
        ----------
        url
            URL to request.
        params
            Optional query parameters.
        timeout
            Request timeout in seconds.
        cache_dir
            Specific cache directory to use. Search pages and monomer pages
            use separate cache folders.

        Returns
        -------
        tuple[str, str]
            ``(html_text, final_url)``. The final URL includes resolved query
            parameters and redirects from the request.
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
        """Extract the SMILES string from a CoPolDB monomer detail page.

        The function uses a simple text scan over the rendered HTML. It expects
        a line labeled ``SMILES`` followed by the corresponding SMILES value.
        """

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
        """Extract molecular weight from a CoPolDB monomer detail page.

        Returns
        -------
        float or None
            Molecular weight if found, otherwise None.
        """

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
        """Parse CoPolDB search-result HTML into ReactivityRatioResult objects.

        This parser relies on the current CoPolDB table layout. It extracts:
        monomer names, monomer detail links, r1, r2, and DOI links.
        """

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
        """Attach SMILES and molecular weights to a search result."""

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
        """Search CoPolDB by monomer names.

        Parameters
        ----------
        monomer1, monomer2
            Monomer names to search.
        exact
            If True, require exact case-insensitive name matches in the
            returned CoPolDB results. This prevents fuzzy matches such as
            returning substituted methacrylates for methyl methacrylate.
        populate_details
            If True, also fetch monomer detail pages to populate SMILES and
            molecular weights.

        Returns
        -------
        list[ReactivityRatioResult]
            Matching CoPolDB results.
        """

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
        """Return CoPolDB candidate results for a pair of SMILES strings.

        CoPolDB's SMILES endpoint is treated as a candidate generator rather
        than an exact matcher. Use :meth:`search_by_smiles_exact` when exact
        canonical RDKit matching is required.
        """

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
        """Search by SMILES and filter candidates by exact canonical SMILES.

        CoPolDB's SMILES search can return related structures. This method:
        1. canonicalizes the user-provided SMILES with RDKit,
        2. queries CoPolDB for candidate results,
        3. fetches candidate monomer detail pages,
        4. canonicalizes database SMILES,
        5. returns only exact canonical matches.

        Parameters
        ----------
        smiles1, smiles2
            Input monomer SMILES strings.
        max_candidates
            Maximum number of candidate rows to inspect. This prevents slow
            scans when CoPolDB returns many fuzzy candidates.
        """

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
        """Delete cached CoPolDB HTML and metadata files."""

        for path in self.cache_dir.glob("**/*"):
            if path.is_file():
                path.unlink()
