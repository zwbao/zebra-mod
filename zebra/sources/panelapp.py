"""PanelApp (Genomics England) and PanelApp Australia: diagnostic panels that contain a gene.

API: /api/v1/genes/{symbol}/ on both hosts (no key). Confidence level:
3 = green (diagnostic grade), 2 = amber, 1 = red, 0 = no rating.
"""

from __future__ import annotations

import urllib.parse
from typing import Any, Dict, List

from zebra.core import Outcome
from zebra.http import get_json, source_record

HOSTS = {
    "GE": ("PanelApp (Genomics England)", "https://panelapp.genomicsengland.co.uk"),
    "AU": ("PanelApp Australia", "https://panelapp-aus.org"),
}
LEVELS = {3: "green", 2: "amber", 1: "red", 0: "no rating"}
CACHE_TTL = 7 * 86400


def parse_results(results: List[Dict[str, Any]], base: str) -> List[Dict[str, Any]]:
    out = []
    for r in results:
        if r.get("entity_type") not in (None, "gene"):
            continue
        panel = r.get("panel") or {}
        try:
            level = int(r.get("confidence_level"))
        except (TypeError, ValueError):
            level = 0
        pid = panel.get("id")
        out.append({
            "panel_id": pid, "panel": panel.get("name"), "version": panel.get("version"),
            "confidence": level, "rating": LEVELS.get(level, str(level)),
            "moi": r.get("mode_of_inheritance") or None,
            "phenotypes": [p for p in (r.get("phenotypes") or []) if p][:3],
            "disease_group": panel.get("disease_group") or None,
            "url": f"{base}/panels/{pid}/gene/{r.get('entity_name')}/" if pid else None,
        })
    out.sort(key=lambda x: (-x["confidence"], x["panel"] or ""))
    return out


def gene_panels(symbol: str, source: str = "GE", limit: int = 20) -> Outcome:
    name, base = HOSTS[source]
    url = f"{base}/api/v1/genes/{urllib.parse.quote(symbol.strip())}/"
    results: List[Dict[str, Any]] = []
    sources = []
    pages = 0
    while url and pages < 4:
        resp = get_json(url, source=name, cache_ttl=CACHE_TTL, timeout=60)
        data = resp.json()
        results.extend(data.get("results") or [])
        sources.append(source_record(name, symbol, resp))
        url = data.get("next")
        pages += 1
    panels = parse_results(results, base)
    summary = {"panels": len(panels), "green": sum(1 for p in panels if p["confidence"] == 3),
               "amber": sum(1 for p in panels if p["confidence"] == 2), "red": sum(1 for p in panels if p["confidence"] == 1)}
    return Outcome({"source": name, "summary": summary, "panels": panels[:limit],
                    "shown": min(limit, len(panels))}, sources=sources)
