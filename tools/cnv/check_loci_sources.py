"""Check that every threshold and reading zebra.cnv quotes is verbatim in the source it cites.

    python3 tools/cnv/check_loci_sources.py <dir-with-saved-sources>

zebra/cnv.py carries the repeat categories (FMR1, HTT, DMPK, FXN, C9orf72), the
SMN1/SMN2 readings and the ACMG section-3 bands as code, each with the sentence it
was read from and the sha256 of the page it was read in (retrieved 2026-10-06).
This script re-checks both against saved copies of those pages:

  NBK1352.txt NBK1384.txt NBK1305.txt NBK1165.txt NBK1281.txt NBK268647.txt
      GeneReviews chapters, saved as the browser's page text (document.body.innerText)
      of https://www.ncbi.nlm.nih.gov/books/<NBK>/ — NCBI Bookshelf answers scripted
      clients with a reCAPTCHA page, so these are saved from a browser.
  PMC3499739.xml PMC11275604.xml
      Europe PMC full text: https://www.ebi.ac.uk/europepmc/webservices/rest/<PMCID>/fullTextXML
  pmc7313390.html
      https://pmc.ncbi.nlm.nih.gov/articles/PMC7313390/ (Riggs et al. 2020)

The GeneReviews texts are not redistributed with zebra (GeneReviews' terms allow
excerpts with credit, not modified copies). A sha256 mismatch means the page changed
since 2026-10-06: re-read the sections the quotes come from and update zebra/cnv.py.
Exit status 1 when any quote is missing.
"""

from __future__ import annotations

import hashlib
import html
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from zebra import cnv as C  # noqa: E402

FILES = {"GR_SMA": "NBK1352.txt", "GR_FMR1": "NBK1384.txt", "GR_HD": "NBK1305.txt", "GR_DM1": "NBK1165.txt",
         "GR_FRDA": "NBK1281.txt", "GR_C9": "NBK268647.txt", "EMQN_DM": "PMC3499739.xml",
         "SMN_EXONS": "PMC11275604.xml", "RIGGS_2020": "pmc7313390.html"}


def text_of(path: Path) -> str:
    raw = path.read_text("utf-8", errors="replace")
    if path.suffix in (".xml", ".html"):
        raw = html.unescape(re.sub(r"<[^>]+>", " ", raw))
    return raw


def squash(text: str) -> str:
    return " ".join(text.split())


def quotes_by_source() -> dict:
    out: dict = {key: [] for key in C.SOURCES}
    for key, src in C.SOURCES.items():
        out[key].extend(src.get("quotes") or [])
    out["GR_SMA"].extend(C.SMA_QUOTES.values())
    out["GR_SMA"].extend(C.SMN2_TABLE_QUOTES)
    for locus in C.REPEAT_LOCI.values():
        bucket = out[locus["source"]]
        bucket.append(locus["inheritance_quote"])
        bucket.extend(locus.get("note_quotes") or [])
        for cat in locus["categories"]:
            bucket.append(cat["quote"])
            bucket.extend(cat.get("risk_quotes") or [])
        for span in (locus.get("span_readings") or {}).values():
            bucket.extend(span["quotes"])
        bucket.extend(rn["quote"] for rn in locus.get("range_notes") or [])
        if locus.get("apparent_homozygosity"):
            bucket.append(locus["apparent_homozygosity"])
    return out


def main(argv: list) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    folder = Path(argv[1])
    missing = 0
    checked = 0
    for key, quotes in quotes_by_source().items():
        path = folder / FILES[key]
        if not path.exists():
            print(f"[skip] {key}: {path.name} not in {folder}")
            continue
        body = path.read_bytes()
        digest = hashlib.sha256(body).hexdigest()
        want = C.SOURCES[key].get("sha256_of_page_text") or C.SOURCES[key].get("sha256_of_body")
        status = "sha256 matches" if digest == want else f"sha256 DIFFERS ({digest[:12]}… vs {want[:12]}…)"
        text = text_of(path)
        flat = squash(text)
        bad = [q for q in quotes if q not in text and squash(q) not in flat]
        checked += len(quotes)
        missing += len(bad)
        print(f"[{'ok' if not bad else 'MISSING'}] {key} ({path.name}, {status}): {len(quotes) - len(bad)}/{len(quotes)} "
              "quotes found")
        for q in bad:
            print(f"    not found: {q!r}")
    print(f"{checked} quotes checked, {missing} missing")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
