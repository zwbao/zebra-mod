"""UniProt REST (protein identity and function) and AlphaFold DB (predicted structure).

UniProt: https://rest.uniprot.org/uniprotkb/{accession} (no key).
AlphaFold DB: https://alphafold.ebi.ac.uk/api/prediction/{accession} (no key).
"""

from __future__ import annotations

import re
import urllib.parse
from typing import Any, Dict, List, Optional

from zebra.core import Outcome
from zebra.http import get_json, source_record

UNIPROT = "https://rest.uniprot.org/uniprotkb"
ALPHAFOLD = "https://alphafold.ebi.ac.uk"
FIELDS = "accession,id,gene_primary,protein_name,length,cc_function,cc_disease"
CACHE_TTL = 30 * 86400
# rest.uniprot.org gzips its JSON twice when asked for gzip; ask for identity.
NO_GZIP = {"Accept-Encoding": "identity"}
_REFS = re.compile(r"\s*\((?:PubMed|By similarity|Ref\.)[^)]*\)")


def _short(text: str, limit: int = 600) -> str:
    text = re.sub(r"\s*\[\s*\]", "", _REFS.sub("", text or "")).strip()
    if len(text) <= limit:
        return text
    cut = text[:limit]
    end = cut.rfind(". ")
    return (cut[: end + 1] if end > 200 else cut.rstrip() + "…")


def parse_entry(d: Dict[str, Any]) -> Dict[str, Any]:
    pd = d.get("proteinDescription") or {}
    name = ((pd.get("recommendedName") or {}).get("fullName") or {}).get("value")
    if not name:
        sub = pd.get("submissionNames") or [{}]
        name = ((sub[0] or {}).get("fullName") or {}).get("value")
    genes = [((g.get("geneName") or {}).get("value")) for g in d.get("genes") or []]
    function = []
    diseases = []
    for c in d.get("comments") or []:
        if c.get("commentType") == "FUNCTION":
            function.extend(t.get("value", "") for t in c.get("texts") or [])
        elif c.get("commentType") == "DISEASE" and c.get("disease"):
            dis = c["disease"]
            xref = dis.get("diseaseCrossReference") or {}
            diseases.append({"name": dis.get("diseaseId"), "acronym": dis.get("acronym"),
                             "omim": f"OMIM:{xref['id']}" if xref.get("database") == "MIM" and xref.get("id") else None})
    acc = d.get("primaryAccession")
    return {
        "accession": acc, "entry_name": d.get("uniProtkbId"), "reviewed": "reviewed" in str(d.get("entryType", "")).lower()
        and "unreviewed" not in str(d.get("entryType", "")).lower(),
        "protein": name, "gene": next((g for g in genes if g), None),
        "length": (d.get("sequence") or {}).get("length"),
        "function": _short(" ".join(function)) or None,
        "diseases": diseases[:10],
        "url": f"https://www.uniprot.org/uniprotkb/{acc}/entry" if acc else None,
    }


def entry(accession: str) -> Outcome:
    acc = accession.strip()
    resp = get_json(f"{UNIPROT}/{urllib.parse.quote(acc)}", source="UniProt", params={"fields": FIELDS, "format": "json"},
                    cache_ttl=CACHE_TTL, timeout=60, headers=NO_GZIP)
    res = parse_entry(resp.json())
    return Outcome(res, sources=[source_record("UniProt", acc, resp, url=res.get("url"))])


def by_gene(symbol: str) -> Outcome:
    """Reviewed human entry for a gene symbol (when the accession is not known)."""
    q = f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true"
    resp = get_json(f"{UNIPROT}/search", source="UniProt", params={"query": q, "fields": FIELDS, "format": "json", "size": 3},
                    cache_ttl=CACHE_TTL, timeout=60, headers=NO_GZIP)
    results = resp.json().get("results") or []
    if not results:
        raise ValueError(f"no reviewed human UniProt entry for gene {symbol}")
    res = parse_entry(results[0])
    warnings = [f"UniProt: {len(results)} reviewed entries for {symbol}; used {res['accession']}"] if len(results) > 1 else []
    return Outcome(res, sources=[source_record("UniProt", res["accession"], resp, url=res.get("url"))], warnings=warnings)


def parse_alphafold(models: List[Dict[str, Any]], accession: str) -> Optional[Dict[str, Any]]:
    if not models:
        return None
    want = f"AF-{accession}-F1"
    m = next((x for x in models if (x.get("modelEntityId") or x.get("entryId")) == want), models[0])
    mid = m.get("modelEntityId") or m.get("entryId")
    return {
        "model": mid, "mean_plddt": m.get("globalMetricValue"), "version": m.get("latestVersion"),
        "residues": [m.get("sequenceStart") or m.get("uniprotStart"), m.get("sequenceEnd") or m.get("uniprotEnd")],
        "pdb_url": m.get("pdbUrl"), "cif_url": m.get("cifUrl"),
        "page": f"{ALPHAFOLD}/entry/{accession}",
        "note": None if mid == want else f"canonical model {want} not served; showing {mid}",
    }


def alphafold(accession: str) -> Outcome:
    resp = get_json(f"{ALPHAFOLD}/api/prediction/{urllib.parse.quote(accession)}", source="AlphaFold DB",
                    cache_ttl=CACHE_TTL, timeout=60, ok_statuses=(200, 404))
    if resp.status == 404:
        return Outcome(None, sources=[source_record("AlphaFold DB", accession, resp)],
                       warnings=[f"AlphaFold DB has no model for {accession}"])
    res = parse_alphafold(resp.json(), accession)
    return Outcome(res, sources=[source_record("AlphaFold DB", accession, resp, url=(res or {}).get("page"))])
