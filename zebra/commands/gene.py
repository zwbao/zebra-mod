"""zebra gene: one gene for rare disease (HGNC, Monarch, ClinGen, gnomAD constraint, PanelApp, UniProt)."""

from __future__ import annotations

import argparse
from typing import Any, Dict, List

from zebra.core import Outcome, UsageError

MOI_SHORT = {
    "Autosomal dominant inheritance": "AD", "Autosomal recessive inheritance": "AR", "X-linked inheritance": "XL",
    "X-linked recessive inheritance": "XLR", "X-linked dominant inheritance": "XLD", "Mitochondrial inheritance": "MT",
    "Semidominant inheritance": "SD", "Y-linked inheritance": "YL",
}


def _moi_panelapp(moi: Any) -> str:
    m = str(moi or "")
    if not m:
        return "-"
    return m.split(",")[0].strip()


def render(r: Dict[str, Any]) -> str:
    h = r.get("hgnc") or {}
    lines: List[str] = []
    lines.append(f"{r['symbol']}  {h.get('hgnc_id') or '-'}  {h.get('name') or ''}  ({h.get('location') or '-'}; {h.get('locus_type') or '-'})")
    e = r.get("ensembl") or {}
    ids = []
    if e.get("id"):
        ids.append(f"{e['id']} chr{e.get('chrom')}:{e.get('start')}-{e.get('end')} ({'+' if e.get('strand') == 1 else '-'}) {e.get('assembly') or ''}".strip())
    if h.get("mane_select"):
        ids.append("MANE " + " / ".join(h["mane_select"]))
    if h.get("omim"):
        ids.append(" ".join(h["omim"]))
    if h.get("uniprot_ids"):
        ids.append("UniProt " + ",".join(h["uniprot_ids"]))
    if ids:
        lines.append("  ".join(ids))
    if h.get("previous_symbols") or h.get("aliases"):
        lines.append(f"previous: {', '.join(h.get('previous_symbols') or []) or '-'}; aliases: {', '.join(h.get('aliases') or []) or '-'}")
    ds = r.get("diseases")
    if ds is None:
        lines.append("diseases (Monarch): unavailable")
    else:
        lines.append(f"diseases (Monarch causal, {len(ds)}):")
        for d in ds[:15]:
            inh = MOI_SHORT.get(d.get("inheritance") or "", d.get("inheritance") or "-")
            xr = f" [{', '.join(d['xrefs'])}]" if d.get("xrefs") else ""
            lines.append(f"  {d['id']:<14} {d.get('name')}{xr}  {inh}  ({'/'.join(d.get('sources') or [])})")
    cg = r.get("clingen") or {}
    val = cg.get("validity")
    if val is None:
        lines.append("ClinGen validity: unavailable")
    elif not val:
        lines.append("ClinGen validity: no curation")
    else:
        lines.append("ClinGen validity:")
        for v in val[:10]:
            lines.append(f"  {v['classification']:<12} {v['disease']} ({v.get('mondo')}, {v.get('moi')}, {v.get('date')}, {v.get('gcep')})")
    dos = cg.get("dosage")
    if dos:
        hi, ts = dos.get("haploinsufficiency") or {}, dos.get("triplosensitivity") or {}
        lines.append(f"ClinGen dosage: HI {hi.get('score')} ({hi.get('description')}); TS {ts.get('score')} ({ts.get('description')}); {dos.get('last_evaluated')}")
    elif cg.get("dosage_note"):
        lines.append(f"ClinGen dosage: {cg['dosage_note']}")
    else:
        lines.append("ClinGen dosage: unavailable")
    c = r.get("constraint")
    if c:
        lines.append(f"constraint ({c.get('version')}): pLI {c.get('pLI')}, LOEUF {c.get('loeuf')} (o/e LoF {c.get('oe_lof')}, "
                     f"{c.get('obs_lof')}/{c.get('exp_lof')}), missense Z {c.get('mis_z')} (o/e {c.get('oe_mis')}), syn Z {c.get('syn_z')}"
                     + (f"; flags {','.join(c['flags'])}" if c.get("flags") else ""))
    else:
        lines.append("constraint: unavailable")
    pa = r.get("panelapp") or {}
    for key, label in (("genomics_england", "PanelApp GE"), ("australia", "PanelApp AU")):
        x = pa.get(key)
        if not x:
            lines.append(f"{label}: unavailable")
            continue
        s = x.get("summary") or {}
        top = "; ".join(f"{p['panel']} ({p['rating']}, {_moi_panelapp(p.get('moi'))})" for p in (x.get("panels") or [])[:4])
        lines.append(f"{label}: {s.get('panels', 0)} panels, {s.get('green', 0)} green, {s.get('amber', 0)} amber, {s.get('red', 0)} red"
                     + (f" — {top}" if top else ""))
    ex = r.get("expression")
    if ex:
        ts = ex.get("tissues") or []
        lines.append("expression (GTEx median TPM): " + (", ".join(
            f"{t.get('tissue')} {t.get('median_tpm')}" for t in ts[:6]) or "no tissue data")
            + (f"  [{ex.get('gencode_id')}]" if ex.get("gencode_id") else ""))
    pr = r.get("protein")
    if pr:
        lines.append(f"protein: {pr.get('accession')} {pr.get('protein')} ({pr.get('length')} aa)")
        if pr.get("function"):
            fn = pr["function"]
            lines.append("  function: " + (fn if len(fn) <= 300 else fn[:300].rsplit(" ", 1)[0] + " …"))
        af = pr.get("alphafold")
        if af:
            lines.append(f"  AlphaFold: {af.get('model')} mean pLDDT {af.get('mean_plddt')} {af.get('page')}"
                         + (f" ({af['note']})" if af.get("note") else ""))
    return "\n".join(lines)


def _gene(args: argparse.Namespace) -> Outcome:
    from zebra.sources import gene as gene_src

    sym = (args.symbol or "").strip()
    if not sym:
        raise UsageError("give an HGNC gene symbol, e.g. SCN1A")
    try:
        out = gene_src.card(sym)
    except ValueError as err:
        raise UsageError(str(err)) from None
    out.text = render(out.result)
    out.query = {"symbol": sym}
    return out


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("gene", help="one gene: HGNC, diseases + inheritance, ClinGen validity/dosage, gnomAD constraint, PanelApp, UniProt")
    p.add_argument("symbol", help="HGNC symbol, e.g. SCN1A")
    p.set_defaults(func=_gene)
