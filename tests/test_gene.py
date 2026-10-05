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
                             "fetch/prev_symbol/C3orf72": "hgnc_prev_c3orf72.json",
                             "fetch/symbol/FOXL2NB": "hgnc_fetch_foxl2nb.json"})
    out = gene.hgnc("C3orf72")
    assert out.result["symbol"] == "FOXL2NB"
    assert any("previous symbol" in w for w in out.warnings)


def test_hgnc_ambiguous_alias_is_not_guessed(monkeypatch):
    empty = "hgnc_fetch_c3orf72.json"  # a real zero-hit fetch response
    _fake_hgnc(monkeypatch, {"fetch/symbol/NAC1": empty, "fetch/prev_symbol/NAC1": empty,
                             "fetch/alias_symbol/NAC1": "hgnc_alias_nac1.json"})
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


# ---------------------------------------------------------------- F15: refuse a junk symbol

@pytest.mark.parametrize("bad,why", [
    ("../search/symbol/SCN1A", "character HGNC symbols do not use"),
    ("SCN1A extra", "space"),
    ("", "give a gene symbol"),
    ("   ", "give a gene symbol"),
    ("SCN1A/CFTR", "character HGNC symbols do not use"),
    ("a" * 40, "longer than any HGNC symbol"),
])
def test_F15_junk_symbols_are_refused_before_any_request(bad, why, monkeypatch):
    def boom(*a, **kw):
        raise AssertionError("no request may be made for a symbol that cannot be one")

    monkeypatch.setattr(gene, "get_json", boom)
    with pytest.raises(UsageError) as err:
        gene.card(bad)
    assert why in str(err.value)


def test_F15_check_symbol_accepts_real_hgnc_spellings():
    for ok in ("SCN1A", "C3orf72", "HLA-A", "MT-TL1", "RNU4ATAC", "NKX2-5", "IGH@", "ATP6V0A2"):
        assert gene.check_symbol(ok) == ok
    assert gene.check_symbol(" SCN1A ") == "SCN1A"


def test_F15_previous_symbol_resolution_uses_the_exact_field_fetch(monkeypatch):
    urls = []

    def fake(url, source, **kw):
        urls.append(url)
        name = {"fetch/symbol/C3orf72": "hgnc_fetch_c3orf72.json",
                "fetch/prev_symbol/C3orf72": "hgnc_prev_c3orf72.json",
                "fetch/symbol/FOXL2NB": "hgnc_fetch_foxl2nb.json"}
        for key, fx in name.items():
            if url.endswith(key):
                return Response(url, 200, json.dumps(load("gene", fx)), "2026-10-06T00:00:00+00:00", True)
        raise AssertionError(f"unexpected HGNC call {url}")

    monkeypatch.setattr(gene, "get_json", fake)
    out = gene.hgnc("C3orf72")
    assert out.result["symbol"] == "FOXL2NB"
    # never the Solr text search, which scored ESPL1/GSDMC/SH3GL1 for "SCN1A extra"
    assert not any("/search/" in u for u in urls)
    assert any("exact prev_symbol match" in w for w in out.warnings)


def test_F15_hgnc_404_is_a_usage_error_not_unavailable(monkeypatch):
    def fake(url, source, **kw):
        return Response(url, 404, "", "2026-10-06T00:00:00+00:00", False)

    monkeypatch.setattr(gene, "get_json", fake)
    with pytest.raises(UsageError) as err:
        gene.card("ZZZZZZ9")
    assert "HGNC has no record" in str(err.value)


def test_F15_an_unreachable_hgnc_says_the_symbol_is_unverified(monkeypatch):
    from zebra.http import SourceError

    def fake(url, source, **kw):
        raise SourceError("HGNC", url, 503, "maintenance")

    monkeypatch.setattr(gene, "get_json", fake)
    monkeypatch.setattr(gene, "_run", lambda tasks: {})
    out = gene.card("SCN1A")
    assert any("unverified symbol" in w and "may mean the symbol is wrong" in w for w in out.warnings)


# ---------------------------------------------------------------- E8: no key in any recorded URL

def test_E8_the_api_key_is_sent_but_never_recorded(monkeypatch):
    from zebra.sources import clinvar

    monkeypatch.setenv("NCBI_API_KEY", "SECRET_NCBI_KEY_123")
    wire = []

    def fake_request(url, source, **kw):
        wire.append({"url": url, "method": kw.get("method", "GET"), "body": kw.get("body"),
                     "headers": kw.get("headers")})
        if "esearch" in url:
            payload = {"esearchresult": {"count": "1", "idlist": ["7105"]}}
        else:
            payload = {"result": {"uids": []}}
        return Response(url, 200, json.dumps(payload), "2026-10-06T00:00:00+00:00", False)

    monkeypatch.setattr("zebra.http.request", fake_request)
    out = clinvar.lookup(rsid="rs113993960")
    # the key really was sent: in the POST body, not the query string
    assert wire and all(w["method"] == "POST" for w in wire)
    assert all("api_key=SECRET_NCBI_KEY_123" in (w["body"] or "") for w in wire)
    assert all(w["headers"]["Content-Type"] == "application/x-www-form-urlencoded" for w in wire)
    # and it is in no URL, no source row and no warning
    assert not any("SECRET_NCBI_KEY_123" in w["url"] for w in wire)
    assert "SECRET_NCBI_KEY_123" not in json.dumps(out.sources)
    assert "SECRET_NCBI_KEY_123" not in json.dumps(out.warnings)
    assert any("eutils.ncbi.nlm.nih.gov" in (s.get("url") or "") for s in out.sources)


def test_E8_without_a_key_the_call_stays_a_get(monkeypatch):
    from zebra.sources import clinvar

    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    seen = []

    def fake_request(url, source, **kw):
        seen.append(kw.get("method", "GET"))
        return Response(url, 200, json.dumps({"esearchresult": {"count": "0", "idlist": []}}),
                        "2026-10-06T00:00:00+00:00", False)

    monkeypatch.setattr("zebra.http.request", fake_request)
    clinvar.search("rs1[VRID]")
    assert seen == ["GET"]


def test_E8_a_rejected_key_is_redacted_from_the_error(monkeypatch):
    from zebra.http import SourceError
    from zebra.sources import clinvar

    monkeypatch.setenv("NCBI_API_KEY", "SECRET_NCBI_KEY_123")

    def fake_request(url, source, **kw):
        raise SourceError("ClinVar", url, 400,
                          '{"error":"API key invalid","api-key":"SECRET_NCBI_KEY_123"}')

    monkeypatch.setattr("zebra.http.request", fake_request)
    with pytest.raises(SourceError) as err:
        clinvar.search("rs1[VRID]")
    assert "SECRET_NCBI_KEY_123" not in str(err.value)
    assert "<NCBI_API_KEY redacted>" in err.value.message


def test_E8_public_url_strips_every_credential_parameter():
    from zebra.sources import public_url, record

    url = "https://e.ncbi/x.fcgi?db=clinvar&api_key=K&term=t&email=a%40b&token=T"
    out = public_url(url)
    assert "api_key" not in out and "email" not in out and "token" not in out
    assert "db=clinvar" in out and "term=t" in out
    assert public_url("https://e/x") == "https://e/x"
    assert public_url(None) is None
    assert record("db", "r", url=url)["url"] == out


def test_E8_no_source_url_in_any_owned_module_can_carry_a_key():
    """Nothing in zebra/sources or the owned commands may put a credential in `params`."""
    import pathlib
    import re as _re

    root = pathlib.Path(gene.__file__).parent.parent
    offenders = []
    for p in list((root / "sources").glob("*.py")) + list((root / "commands").glob("*.py")):
        text = p.read_text("utf-8")
        for m in _re.finditer(r'"(api_key|apikey|key|token|access_token|password|secret)"\s*\]?\s*=', text):
            line = text[:m.start()].count("\n") + 1
            if p.name == "clinvar.py":
                continue  # clinvar builds a POST body, checked above
            offenders.append(f"{p.name}:{line} {m.group(1)}")
    assert offenders == [], offenders


# ---------------------------------------------------------------- F11: one bad shape, one card section

def test_F11_source_layer_attempt_catches_attribute_error():
    from zebra.sources import attempt as src_attempt

    warnings = []

    def upstream_changed_shape():
        return ["not", "a", "dict"].get("field")

    assert src_attempt("Europe PMC", upstream_changed_shape, warnings) is None
    assert len(warnings) == 1 and "AttributeError" in warnings[0]


def test_F11_a_usage_error_still_escapes_so_a_bad_symbol_is_not_hidden():
    from zebra.sources import attempt as src_attempt

    warnings = []
    with pytest.raises(UsageError):
        src_attempt("HGNC", lambda: (_ for _ in ()).throw(UsageError("bad symbol")), warnings)


# ---------------------------------------------------------------- F10: never cache an error body

@pytest.mark.parametrize("body,status,why", [
    ("", 200, "empty body"),
    ("<html>maintenance</html>", 200, "HTML page"),
    ("not json at all", 200, "not JSON"),
    ('{"errors":[{"message":"upstream timeout"}]}', 200, "error body"),
    ('{"error":"nope"}', 200, "error body"),
])
def test_F10_validated_json_rejects_a_200_that_is_really_an_error(body, status, why):
    from zebra.http import SourceError
    from zebra.sources import validated_json

    resp = Response("https://u/x", status, body, "2026-10-06T00:00:00+00:00", False)
    with pytest.raises(SourceError) as err:
        validated_json(resp, "Upstream")
    assert why in err.value.message


def test_F10_validated_json_requires_the_key_the_parser_needs():
    from zebra.http import SourceError
    from zebra.sources import validated_json

    resp = Response("https://u/x", 200, '{"something": 1}', "2026-10-06T00:00:00+00:00", False)
    with pytest.raises(SourceError) as err:
        validated_json(resp, "ClinVar", require="esearchresult")
    assert "has no 'esearchresult'" in err.value.message


def test_F10_a_bad_body_from_the_cache_is_refetched_once():
    from zebra.sources import validated_json

    cached = Response("https://u/x", 200, '{"errors":[{"message":"timeout"}]}', "2026-10-05T00:00:00+00:00", True)
    fresh = Response("https://u/x", 200, '{"data": 1}', "2026-10-06T00:00:00+00:00", False)
    calls = []

    def refetch():
        calls.append(1)
        return fresh

    assert validated_json(cached, "Upstream", refetch=refetch) == {"data": 1}
    assert len(calls) == 1


def test_F10_a_bad_fresh_body_is_not_refetched():
    from zebra.http import SourceError
    from zebra.sources import validated_json

    fresh_bad = Response("https://u/x", 200, "", "2026-10-06T00:00:00+00:00", False)
    with pytest.raises(SourceError):
        validated_json(fresh_bad, "Upstream", refetch=lambda: pytest.fail("must not refetch a fresh body"))


def test_F10_validated_text_rejects_html_and_a_missing_header():
    from zebra.http import SourceError
    from zebra.sources import validated_text

    html = Response("https://u/x", 200, "<html>maintenance</html>", "2026-10-06T00:00:00+00:00", False)
    with pytest.raises(SourceError) as err:
        validated_text(html, "ClinGen validity", must_contain="GENE SYMBOL")
    assert "HTML page" in err.value.message
    no_header = Response("https://u/x", 200, "a,b,c\n1,2,3\n", "2026-10-06T00:00:00+00:00", False)
    with pytest.raises(SourceError) as err:
        validated_text(no_header, "ClinGen validity", must_contain="GENE SYMBOL")
    assert "not the expected table" in err.value.message
    good = Response("https://u/x", 200, "GENE SYMBOL,X\nSCN1A,1\n", "2026-10-06T00:00:00+00:00", False)
    assert validated_text(good, "ClinGen validity", must_contain="GENE SYMBOL").startswith("GENE SYMBOL")


def test_F10_a_clingen_maintenance_page_is_rejected_and_refetched(monkeypatch):
    from zebra.http import SourceError
    from zebra.sources import clingen as cg

    ttls = []

    def fake_request(url, source, **kw):
        ttls.append(kw.get("cache_ttl"))
        return Response(url, 200, "<html>maintenance</html>", "2026-10-06T00:00:00+00:00", len(ttls) == 1)

    monkeypatch.setattr(cg, "request", fake_request)
    with pytest.raises(SourceError) as err:
        cg.validity(symbol="SCN1A")
    assert len(ttls) == 2 and ttls[1] == 0
    assert "HTML page" in err.value.message  # not "header row 'GENE SYMBOL' not found"


def test_F10_an_empty_genereviews_map_is_rejected_not_parsed_to_nothing(monkeypatch):
    from zebra.http import SourceError
    from zebra.sources import genereviews as gr

    def fake_request(url, source, **kw):
        return Response(url, 200, "", "2026-10-06T00:00:00+00:00", False)

    monkeypatch.setattr(gr, "request", fake_request)
    with pytest.raises(SourceError) as err:
        gr._fetch_maps()
    assert "empty body" in err.value.message
