"""Gene card: HGNC, Monarch, ClinGen, PanelApp, UniProt/AlphaFold parsing (offline, real captured
responses), card assembly with sources replaced, and SCN1A / NGLY1 / DMD live (pytest -m live)."""

import json
import os

import pytest

from zebra.core import Outcome, UsageError
from zebra.http import Response
from zebra.sources import clingen, gene, gnomad, panelapp, uniprot

FIX = os.path.join(os.path.dirname(__file__), "fixtures")


def load(*parts):
    with open(os.path.join(FIX, *parts)) as fh:
        return json.load(fh)


def text(*parts):
    with open(os.path.join(FIX, *parts)) as fh:
        return fh.read()


# ---------------------------------------------------------------- HGNC

def _fake_hgnc(monkeypatch, mapping):
    def fake(url, source, **kw):
        for key, name in mapping.items():
            if url.endswith(key):
                return Response(url, 200, json.dumps(load("gene", name)), "2026-10-05T00:00:00+00:00", True)
        raise AssertionError(f"unexpected HGNC call {url}")
    monkeypatch.setattr(gene, "get_json", fake)


def test_parse_hgnc():
    rec = gene.parse_hgnc(load("gene", "hgnc_fetch_scn1a.json")["response"]["docs"][0])
    assert rec["hgnc_id"] == "HGNC:10585" and rec["ensembl_gene_id"] == "ENSG00000144285"
    assert "NM_001165963.4" in rec["mane_select"] and rec["uniprot_ids"] == ["P35498"]
    assert rec["omim"] == ["OMIM:182389"]


def test_hgnc_previous_symbol_followed(monkeypatch):
    _fake_hgnc(monkeypatch, {"fetch/symbol/C3orf72": "hgnc_fetch_c3orf72.json",
                             "search/prev_symbol/C3orf72": "hgnc_prev_c3orf72.json",
                             "fetch/symbol/FOXL2NB": "hgnc_fetch_foxl2nb.json"})
    out = gene.hgnc("C3orf72")
    assert out.result["symbol"] == "FOXL2NB"
    assert any("previous symbol" in w for w in out.warnings)


def test_hgnc_ambiguous_alias_is_not_guessed(monkeypatch):
    empty = "hgnc_fetch_c3orf72.json"  # a real zero-hit fetch response
    _fake_hgnc(monkeypatch, {"fetch/symbol/NAC1": empty, "search/prev_symbol/NAC1": empty,
                             "search/alias_symbol/NAC1": "hgnc_alias_nac1.json"})
    with pytest.raises(UsageError) as err:
        gene.hgnc("NAC1")
    assert "NACC1" in str(err.value) and "SCN1A" in str(err.value)


# ---------------------------------------------------------------- Monarch

def test_monarch_dedupes_on_disease():
    items = load("gene", "monarch_assoc_scn1a.json")["items"]
    rows = gene.parse_monarch_associations(items)
    ids = [r["id"] for r in rows]
    assert len(ids) == len(set(ids)) == 8 and len(items) == 9
    dee6a = next(r for r in rows if r["id"] == "MONDO:0100079")
    assert dee6a["xrefs"] == ["OMIM:607208"] and dee6a["sources"] == ["omim"]
    assert next(r for r in rows if r["id"] == "MONDO:0100135")["sources"] == ["clingen"]
    ent = load("gene", "monarch_entity_MONDO_0100079.json")
    assert ent["inheritance"]["name"] == "Autosomal dominant inheritance"


# ---------------------------------------------------------------- ClinGen

def test_clingen_validity_csv():
    rows = clingen.parse_validity_csv(text("clingen", "validity_excerpt.csv"), symbol="SCN1A", hgnc_id="HGNC:10585")
    assert rows[0]["disease"] == "Dravet syndrome" and rows[0]["mondo"] == "MONDO:0100135"
    assert rows[0]["classification"] == "Definitive" and rows[0]["moi"] == "AD"
    assert rows[0]["gcep"].startswith("Epilepsy") and rows[0]["url"].startswith("https://search.clinicalgenome.org/")
    assert {r["classification"] for r in rows} >= {"Definitive", "Moderate"}
    assert all(r["gene"] == "SCN1A" for r in rows)
    ngly1 = clingen.parse_validity_csv(text("clingen", "validity_excerpt.csv"), symbol="NGLY1")
    assert ngly1 and ngly1[0]["moi"] == "AR"
    assert clingen.parse_validity_csv(text("clingen", "validity_excerpt.csv"), symbol="NOSUCHGENE") == []


def test_clingen_dosage_tsv():
    t = text("clingen", "dosage_excerpt_grch38.tsv")
    s = clingen.parse_dosage_tsv(t, "SCN1A")
    assert s["haploinsufficiency"]["score"] == "3" and s["triplosensitivity"]["score"] == "0"
    assert s["haploinsufficiency"]["pmids"]
    assert clingen.parse_dosage_tsv(t, "NGLY1")["haploinsufficiency"]["score"] == "30"
    assert clingen.parse_dosage_tsv(t, "NOSUCHGENE") is None


def test_clingen_bad_file_raises():
    with pytest.raises(ValueError):
        clingen.parse_validity_csv("<html>error</html>", symbol="SCN1A")


# ---------------------------------------------------------------- PanelApp, UniProt, AlphaFold

def test_panelapp_parse():
    d = load("panelapp", "scn1a_ge.json")
    rows = panelapp.parse_results(d["results"], panelapp.HOSTS["GE"][1])
    assert [r["confidence"] for r in rows] == sorted((r["confidence"] for r in rows), reverse=True)
    assert rows[0]["rating"] == "green" and rows[-1]["rating"] == "red"
    assert rows[0]["url"].startswith("https://panelapp.genomicsengland.co.uk/panels/")
    au = panelapp.parse_results(load("panelapp", "ngly1_au.json")["results"], panelapp.HOSTS["AU"][1])
    assert all("BIALLELIC" in (r["moi"] or "") for r in au if r["confidence"] == 3)


def test_uniprot_parse():
    p = uniprot.parse_entry(load("uniprot", "P35498.json"))
    assert p["accession"] == "P35498" and p["gene"] == "SCN1A" and p["length"] == 2009
    assert p["function"] and "PubMed" not in p["function"] and len(p["function"]) <= 601
    assert any(d["omim"] == "OMIM:607208" for d in p["diseases"])


def test_alphafold_parse():
    a = uniprot.parse_alphafold(load("uniprot", "alphafold_P35498.json"), "P35498")
    assert a["model"] == "AF-P35498-F1" and a["note"] is None and a["pdb_url"].endswith(".pdb")
    d = uniprot.parse_alphafold(load("uniprot", "alphafold_P11532.json"), "P11532")
    assert d["model"] != "AF-P11532-F1" and "not served" in d["note"]


# ---------------------------------------------------------------- card assembly, offline

@pytest.fixture
def offline(monkeypatch):
    hg = gene.parse_hgnc(load("gene", "hgnc_fetch_scn1a.json")["response"]["docs"][0])
    monkeypatch.setattr(gene, "hgnc", lambda s: Outcome(hg, sources=[{"db": "HGNC"}]))
    monkeypatch.setattr(gene.ensembl, "lookup_symbol", lambda s: Outcome(
        {"id": "ENSG00000144285", "seq_region_name": "2", "start": 165984641, "end": 166182806, "strand": -1,
         "biotype": "protein_coding", "assembly_name": "GRCh38"}, sources=[{"db": "Ensembl lookup"}]))

    def monarch(hid):
        rows = gene.parse_monarch_associations(load("gene", "monarch_assoc_scn1a.json")["items"])
        ent = load("gene", "monarch_entity_MONDO_0100079.json")
        for r in rows:
            if r["id"] == ent["id"]:
                r["inheritance"] = ent["inheritance"]["name"]
        return Outcome({"total": 9, "diseases": rows}, sources=[{"db": "Monarch"}])

    monkeypatch.setattr(gene, "monarch_diseases", monarch)
    monkeypatch.setattr(gene.clingen, "validity", lambda symbol=None, hgnc_id=None: Outcome(
        clingen.parse_validity_csv(text("clingen", "validity_excerpt.csv"), symbol=symbol, hgnc_id=hgnc_id)))
    monkeypatch.setattr(gene.clingen, "dosage", lambda s, assembly="GRCh38": Outcome(
        clingen.parse_dosage_tsv(text("clingen", "dosage_excerpt_grch38.tsv"), s)))
    monkeypatch.setattr(gene.gnomad, "gene_constraint", lambda s, a="GRCh38", gene_id=None: Outcome(
        gnomad.parse_constraint(load("gnomad", "constraint_scn1a.json")["data"]["gene"], "GRCh38")))

    def pa(s, source="GE", limit=20):
        if source == "AU":
            from zebra.http import SourceError
            raise SourceError("PanelApp Australia", "https://panelapp-aus.org/api/v1/genes/SCN1A/", 503, "down")
        rows = panelapp.parse_results(load("panelapp", "scn1a_ge.json")["results"], panelapp.HOSTS["GE"][1])
        return Outcome({"source": "PanelApp (Genomics England)", "summary": {"panels": len(rows), "green": 3, "amber": 0, "red": 1},
                        "panels": rows, "shown": len(rows)})

    monkeypatch.setattr(gene.panelapp, "gene_panels", pa)
    monkeypatch.setattr(gene.uniprot, "entry", lambda acc: Outcome(uniprot.parse_entry(load("uniprot", "P35498.json"))))
    monkeypatch.setattr(gene.uniprot, "alphafold", lambda acc: Outcome(uniprot.parse_alphafold(load("uniprot", "alphafold_P35498.json"), acc)))


def test_card_offline_scn1a(offline):
    out = gene.card("SCN1A")
    r = out.result
    assert r["hgnc"]["hgnc_id"] == "HGNC:10585" and r["ensembl"]["id"] == "ENSG00000144285"
    dravet = next(v for v in r["clingen"]["validity"] if v["mondo"] == "MONDO:0100135")
    assert dravet["classification"] == "Definitive"
    # Monarch has no inheritance for MONDO:0100135; ClinGen's MOI fills it, labelled as such
    d = next(x for x in r["diseases"] if x["id"] == "MONDO:0100135")
    assert d["inheritance"] == "AD" and d["inheritance_source"] == "ClinGen MOI"
    assert next(x for x in r["diseases"] if x["id"] == "MONDO:0100079")["inheritance_source"].startswith("Monarch")
    assert r["constraint"]["loeuf"] < 0.2
    assert r["clingen"]["dosage"]["haploinsufficiency"]["score"] == "3"
    assert r["protein"]["alphafold"]["model"] == "AF-P35498-F1"
    # a failing source is a named warning, never a silent gap
    assert r["panelapp"]["australia"] is None
    assert any(w.startswith("PanelApp Australia unavailable") for w in out.warnings)


def test_render_gene(offline):
    from zebra.commands.gene import render

    t = render(gene.card("SCN1A").result)
    assert "Definitive   Dravet syndrome" in t and "LOEUF 0.107" in t and "PanelApp AU: unavailable" in t


# ---------------------------------------------------------------- live

@pytest.mark.live
def test_live_scn1a():
    r = gene.card("SCN1A").result
    assert r["hgnc"]["hgnc_id"] == "HGNC:10585"
    assert any(v["mondo"] == "MONDO:0100135" and v["classification"] == "Definitive" for v in r["clingen"]["validity"])
    assert r["constraint"]["loeuf"] < 0.3 and r["clingen"]["dosage"]["haploinsufficiency"]["score"] == "3"
    assert r["panelapp"]["genomics_england"]["summary"]["green"] >= 1


@pytest.mark.live
def test_live_ngly1_recessive():
    r = gene.card("NGLY1").result
    assert any(v["moi"] == "AR" for v in r["clingen"]["validity"])
    assert any((d.get("inheritance") or "").startswith(("Autosomal recessive", "AR")) for d in r["diseases"])
    assert r["protein"]["accession"] == "Q96IV0"


@pytest.mark.live
def test_live_dmd_xlinked():
    r = gene.card("DMD").result
    assert any(v["moi"] == "XL" for v in r["clingen"]["validity"])
    assert any(d["id"] == "MONDO:0010679" for d in r["diseases"])  # Duchenne
    assert r["constraint"]["pLI"] > 0.9


@pytest.mark.live
def test_live_unknown_symbol():
    with pytest.raises(UsageError):
        gene.card("NOTAGENE123")
