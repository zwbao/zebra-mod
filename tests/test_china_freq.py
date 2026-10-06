"""Chinese-population frequencies (zebra.sources.china_freq).

Offline tests replay exchanges recorded live on 2026-10-06 (tests/fixtures/china_freq/*.json:
every request the module made and the body each service returned) through a stand-in for
`zebra.http.request` that applies the same body checks the real one does. `pytest -m live`
runs two lookups against the real services.
"""

import json
import os

import pytest

from zebra import http
from zebra.core import Outcome, UsageError
from zebra.http import Response, SourceError
from zebra.sources import china_freq as cf

FIX = os.path.join(os.path.dirname(__file__), "fixtures", "china_freq")
ALL6 = ["NyuWa", "WBBC", "1000 Genomes CHB", "1000 Genomes CHS", "1000 Genomes CDX", "Taiwan Biobank WGS"]
CONTRACT_KEYS = {"name", "population", "af", "ac", "an", "url"}


def fixture(name):
    with open(os.path.join(FIX, name)) as fh:
        return json.load(fh) if name.endswith(".json") else fh.read()


def exchange(fx, needle, method=None):
    hits = [e for e in fx["exchanges"] if needle in e["url"] and (method is None or e["method"] == method)]
    assert len(hits) == 1, (needle, [e["url"] for e in fx["exchanges"]])
    return hits[0]


def fake_request(fx, override=None, calls=None):
    """`zebra.http.request` answering from recorded exchanges, with the real body checks.

    `override(method, url, body)` may return (status, text), an exception to raise, or None
    to fall through to the recording.
    """
    def fake(url, *, source, params=None, method="GET", body=None, headers=None, timeout=30.0, retries=3,
             cache_ttl=0, accept="application/json", ok_statuses=(200,), validate=None, refresh=False,
             not_found_ttl=None):
        assert params is None, "china_freq builds its URLs itself"
        if calls is not None:
            calls.append((method, url, body))
        got = override(method, url, body) if override else None
        if isinstance(got, Exception):
            raise got
        if got is None:
            hits = [e for e in (fx or {"exchanges": []})["exchanges"]
                    if e["method"] == method and e["url"] == url and e["body"] == body]
            if not hits:
                raise AssertionError(f"unrecorded request {method} {url} {body!r}")
            e = hits[0]
            if "error" in e:
                raise SourceError(source, url, e["error"]["status"], e["error"]["message"])
            got = (e["status"], e["text"])
        status, text = got
        if status not in ok_statuses:
            raise SourceError(source, url, status, text[:300])
        if status == 200:  # zebra.http.request's own checks, in its order
            bad = None
            if accept == "application/json" and text.lstrip()[:1] == "<":
                bad = "HTML page where JSON was expected"
            elif validate is not None:
                bad = validate(text)
            if bad is not None:
                raise SourceError(source, url, status, bad)
        return Response(url, status, text, "2026-10-06T10:00:00+00:00", True)
    return fake


def run(monkeypatch, name, override=None, calls=None, args=None):
    fx = fixture(name)
    monkeypatch.setattr(cf, "request", fake_request(fx, override, calls))
    return cf.lookup(*(args or fx["args"]))


def by_name(out):
    return {d["name"]: d for d in out.result["datasets"]}


def absent_names(out):
    return [d["name"] for d in out.result["not_found_in"]]


def not_checked(out):
    return {d["name"]: d["reason"] for d in out.result["not_checked"]}


# ---------------------------------------------------------------- found

def test_found_common_variant_rs671(monkeypatch):
    out = run(monkeypatch, "rs671_grch38.json")
    assert isinstance(out, Outcome)
    r = out.result
    assert r["variant"] == "12-111803962-G-A" and r["assembly"] == "GRCh38"
    assert r["checked"] == ALL6 and r["not_found_in"] == []
    d = by_name(out)
    assert set(d) == set(ALL6)
    for row in r["datasets"]:
        assert CONTRACT_KEYS <= set(row), row["name"]
        assert row["url"].startswith("http")
    assert (d["NyuWa"]["ac"], d["NyuWa"]["an"], d["NyuWa"]["homozygotes"]) == (1323, 5998, 151)
    assert d["NyuWa"]["af"] == pytest.approx(1323 / 5998)
    assert (d["WBBC"]["ac"], d["WBBC"]["an"], d["WBBC"]["samples"]) == (2138, 8960, 4480)
    assert d["WBBC"]["af"] == pytest.approx(0.238616071)
    assert {s["population"]: s["af"] for s in d["WBBC"]["subpopulations"]}["Lingnan"] == pytest.approx(0.357142857)
    assert (d["1000 Genomes CHB"]["ac"], d["1000 Genomes CHB"]["an"]) == (33, 206)
    assert (d["1000 Genomes CHS"]["ac"], d["1000 Genomes CHS"]["an"]) == (57, 210)
    assert (d["1000 Genomes CDX"]["ac"], d["1000 Genomes CDX"]["an"]) == (8, 186)
    assert "small sample" in d["1000 Genomes CHB"]["population"]
    twb = d["Taiwan Biobank WGS"]
    assert twb["af"] == pytest.approx(0.2825) and twb["ac"] is None and twb["samples"] == 1492
    # provenance for every request, through zebra.sources.record
    assert len(out.sources) == 6
    assert all(s["db"] and s["url"] and s["retrieved_at"] for s in out.sources)
    assert {s["db"] for s in out.sources} >= {"NyuWa NCVD", "WBBC", "Ensembl variation", "Taiwan Biobank (TaiwanView)"}
    # the resources that are not queried are named, never silent
    nc = not_checked(out)
    for name in ("ChinaMAP", "CMDB", "NyuWa Global 10K T2T", "PGG.Han 2.0", "PGG.SNV"):
        assert name in nc and nc[name]
        assert any(name in w for w in out.warnings)
    assert any("gnomAD" in n and "not Chinese-specific" in n for n in r["notes"])


def test_found_rare_variant_slc26a4(monkeypatch):
    """SLC26A4 c.919-2A>G: rare, in the Chinese cohorts, without 1000 Genomes genotypes."""
    out = run(monkeypatch, "slc26a4_grch38.json")
    d = by_name(out)
    assert (d["NyuWa"]["ac"], d["NyuWa"]["an"]) == (33, 5998)
    assert (d["WBBC"]["ac"], d["WBBC"]["an"]) == (36, 8960)
    assert d["Taiwan Biobank WGS"]["af"] == pytest.approx(0.008)
    for code in ("CHB", "CHS", "CDX"):
        assert f"1000 Genomes {code}" not in d  # no genotypes in Ensembl: not reported as af 0
        assert f"1000 Genomes {code}" in absent_names(out)


def test_indel_is_left_aligned_and_matched_everywhere(monkeypatch):
    """CFTR F508del given 3'-shifted (7-117559591-TCTT-T, Ensembl's spelling) is asked as 7-117559590-ATCT-A."""
    out = run(monkeypatch, "f508del_unaligned_grch38.json")
    r = out.result
    assert r["variant"] == "7-117559590-ATCT-A" and r["input"] == "7-117559591-TCTT-T"
    nyuwa = exchange(fixture("f508del_unaligned_grch38.json"), "NyuWa")
    assert nyuwa["url"].endswith("variant=7-117559590-ATCT-A")
    d = by_name(out)
    assert (d["NyuWa"]["ac"], d["NyuWa"]["an"]) == (1, 5998)
    assert (d["WBBC"]["ac"], d["WBBC"]["an"]) == (2, 8960)
    # 1000 Genomes genotyped the site (CHB: 206 x TCTT): a real count of 0, labelled as such
    chb = d["1000 Genomes CHB"]
    assert (chb["ac"], chb["an"], chb["af"]) == (0, 206, 0.0) and "not observed" in chb["note"]
    assert chb["rsid"] == "rs113993960"
    assert "Taiwan Biobank WGS" in absent_names(out)


# ---------------------------------------------------------------- not found

def test_not_found_is_not_reported_as_af_zero(monkeypatch):
    out = run(monkeypatch, "scn1a_absent_grch38.json")
    r = out.result
    assert r["datasets"] == []
    assert r["checked"] == ALL6
    assert absent_names(out) == ALL6
    for item in r["not_found_in"]:
        assert "not evidence the allele is absent" in item["note"] and "covers" in item["note"]
    nyuwa = r["not_found_in"][0]["note"]
    assert "2,999" in nyuwa
    # every dataset that answered "no record" still left provenance
    assert {s["db"] for s in out.sources} >= {"NyuWa NCVD", "WBBC", "Taiwan Biobank (TaiwanView)"}


def test_position_with_another_allele_is_not_found(monkeypatch):
    """WBBC holds G>A at rs671; G>T there is "no record", not af 0, and names what it does hold."""
    fx = fixture("rs671_grch38.json")
    nyuwa_absent = exchange(fixture("scn1a_absent_grch38.json"), "NyuWa")["text"]

    def override(method, url, body):
        if "NyuWa" in url:
            assert url.endswith("variant=12-111803962-G-T")
            return 200, nyuwa_absent
        return None
    out = run(monkeypatch, "rs671_grch38.json", override, args=("12", 111803962, "G", "T"))
    assert out.result["datasets"] == []
    notes = {d["name"]: d["note"] for d in out.result["not_found_in"]}
    assert "it holds G>A there" in notes["WBBC"]
    assert "no dbSNP record in Ensembl is this allele" in notes["1000 Genomes CHB"]
    # Taiwan Biobank is searchable by rsID only and no rsID is this allele: named, not silent
    assert "Taiwan Biobank WGS" in not_checked(out)
    assert any(w.startswith("Taiwan Biobank WGS not queried") for w in out.warnings)
    assert fx  # fixture loaded


# ---------------------------------------------------------------- failing resources

def test_failing_resource_is_a_warning_and_the_others_answer(monkeypatch):
    def override(method, url, body):
        if "wbbc" in url:
            return SourceError("WBBC", url, None, "network error: <urlopen error timed out>")
        return None
    out = run(monkeypatch, "rs671_grch38.json", override)
    assert any(w.startswith("WBBC unavailable: network error") for w in out.warnings)
    assert "WBBC" not in out.result["checked"] and "WBBC" not in by_name(out)
    assert "WBBC" in not_checked(out)
    assert set(by_name(out)) == set(ALL6) - {"WBBC"}


def test_failing_ensembl_rsid_lookup_drops_1000g_and_taiwan_with_reasons(monkeypatch):
    def override(method, url, body):
        if "/variation/human/" in url:
            return SourceError("Ensembl variation", url, 503, "Service Unavailable")
        return None
    out = run(monkeypatch, "rs671_grch38.json", override)
    nc = not_checked(out)
    for code in ("CHB", "CHS", "CDX"):
        assert f"1000 Genomes {code}" in nc
    assert any(w.startswith("1000 Genomes (Ensembl) unavailable") for w in out.warnings)
    assert set(by_name(out)) == {"NyuWa", "WBBC", "Taiwan Biobank WGS"}  # TaiwanView used the rsID from overlap


def test_reference_unavailable_still_answers_with_a_warning(monkeypatch):
    def override(method, url, body):
        if "/sequence/region/" in url:
            return SourceError("Ensembl sequence", url, None, "network error: timed out")
        return None
    out = run(monkeypatch, "rs671_grch38.json", override)
    assert any("not checked against GRCh38" in w for w in out.warnings)
    assert set(by_name(out)) == set(ALL6)


# ---------------------------------------------------------------- HTML with HTTP 200

def test_html_200_from_wbbc_is_rejected(monkeypatch):
    page = fixture("wbbc_search_page.html")

    def override(method, url, body):
        return (200, page) if "wbbc" in url else None
    out = run(monkeypatch, "rs671_grch38.json", override)
    assert any(w.startswith("WBBC unavailable: HTML page where JSON was expected") for w in out.warnings)
    assert "WBBC" not in by_name(out) and "WBBC" in not_checked(out)
    assert "NyuWa" in by_name(out)


def test_unexpected_html_page_from_nyuwa_is_rejected(monkeypatch):
    """NyuWa answers HTML by design; a page with neither a record nor its no-result message is refused."""
    page = fixture("nyuwa_home_page.html")
    assert cf._nyuwa_check(page) is not None

    def override(method, url, body):
        return (200, page) if "NyuWa" in url else None
    out = run(monkeypatch, "rs671_grch38.json", override)
    assert any(w.startswith("NyuWa (NCVD) unavailable: neither a NyuWa variant record") for w in out.warnings)
    assert "NyuWa" not in by_name(out) and "NyuWa" in not_checked(out)
    assert "WBBC" in by_name(out)


def test_wbbc_and_taiwanview_error_text_rejected():
    for check in (cf._json_list_or_empty("WBBC"), cf._json_list_or_empty("TaiwanView")):
        assert check("ERROR:0") and check("<html><body>502</body></html>") and check("[{bad json")
        assert check("") is None and check("[]") is None


def test_nyuwa_page_for_another_variant_is_refused():
    page = exchange(fixture("rs671_grch38.json"), "NyuWa")["text"]
    assert cf.parse_nyuwa(page, "12-111803962-G-A")["ac"] == 1323
    with pytest.raises(ValueError, match="not the 12-111803962-G-T"):
        cf.parse_nyuwa(page, "12-111803962-G-T")


def test_html_body_is_never_cached_by_the_real_client(monkeypatch, tmp_path):
    """Through the real zebra.http.request: the HTML page is refused and not written to the cache."""
    monkeypatch.delenv("ZEBRA_DEADLINE_MS", raising=False)
    monkeypatch.delenv("ZEBRA_NO_CACHE", raising=False)
    monkeypatch.setattr(http.time, "sleep", lambda s: None)
    bodies = [fixture("nyuwa_home_page.html")]

    class Resp:
        status = 200
        headers = {}

        def __init__(self, data):
            self.data = data

        def read(self):
            return self.data

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    class Opener:
        def open(self, req, timeout=None):
            return Resp(bodies[0].encode("utf-8"))

    monkeypatch.setattr(http, "_opener", lambda: Opener())
    v = ("12", 111803962, "G", "A")
    sources = []
    with pytest.raises(SourceError, match="neither a NyuWa variant record"):
        cf._nyuwa(v, None, sources)
    cache = http.cache_dir() / "http"
    assert not cache.exists() or not any(cache.rglob("*.json"))
    bodies[0] = exchange(fixture("rs671_grch38.json"), "NyuWa")["text"]
    got = cf._nyuwa(v, None, sources)
    assert got[0]["row"]["ac"] == 1323 and any(cache.rglob("*.json"))


# ---------------------------------------------------------------- input

@pytest.mark.parametrize("args", [
    ("23", 1, "A", "G"), ("chrZ", 1, "A", "G"), ("", 1, "A", "G"), (None, 1, "A", "G"), (True, 1, "A", "G"),
    ("12", 0, "A", "G"), ("12", -5, "A", "G"), ("12", "12a", "A", "G"), ("12", True, "A", "G"),
    ("12", 1.5, "A", "G"), ("12", 300_000_000, "A", "G"),
    ("12", 1, "N", "G"), ("12", 1, "-", "G"), ("12", 1, "", "G"), ("12", 1, "A", "A"), ("12", 1, "A", "a"),
    ("12", 1, "A", "<DEL>"), ("12", 1, "A", "T,C"), ("12", 1, None, "G"),
])
def test_bad_input_is_a_usage_error_before_any_request(args, monkeypatch):
    calls = []
    monkeypatch.setattr(cf, "request", fake_request(None, calls=calls))
    with pytest.raises(UsageError):
        cf.lookup(*args)
    assert calls == []


@pytest.mark.parametrize("assembly", ["hg19", "GRCh36", "", None])
def test_bad_assembly_is_a_usage_error(assembly, monkeypatch):
    monkeypatch.setattr(cf, "request", fake_request(None, calls=[]))
    with pytest.raises(UsageError, match="GRCh38 or GRCh37"):
        cf.lookup("12", 111803962, "G", "A", assembly)


def test_structural_variant_is_answered_with_named_gaps_not_an_error(monkeypatch):
    """A valid 1,001-bp deletion must not raise: the variant card calls lookup() and a UsageError would end it."""
    calls = []
    monkeypatch.setattr(cf, "request", fake_request(None, calls=calls))
    out = cf.lookup("12", 1000, "A" * 1001, "A")
    assert calls == [] and out.result["datasets"] == [] and out.result["checked"] == []
    nc = not_checked(out)
    assert all("structural variants" in nc[n] for n in ALL6)
    assert sum(w.endswith("not looked up") for w in out.warnings) == 6


def test_accepted_input_spellings():
    assert cf.parse_input("chr12", "111803962", "g", "a", "grch38") == (("12", 111803962, "G", "A"), "GRCh38")
    assert cf.parse_input("chrM", 3243, "A", "G", "GRCh37")[0][0] == "MT"
    assert cf.parse_input("M", 3243, "A", "G", "GRCh38")[0][0] == "MT"
    assert cf.parse_input("CHRX", 5, "A", "G", "GRCh38")[0][0] == "X"
    assert cf.parse_input(7, 5, "A", "G", "GRCh38")[0][0] == "7"


def test_position_past_chromosome_end_is_a_usage_error(monkeypatch):
    def override(method, url, body):
        return SourceError("Ensembl sequence", url, 400,
                           '{"error":"Cannot request a slice whose start (48000000) is greater than 46709983 for 21."}')
    monkeypatch.setattr(cf, "request", fake_request(None, override))
    with pytest.raises(UsageError, match="beyond the end of chromosome 21"):
        cf.lookup("21", 48000000, "A", "G")


# ---------------------------------------------------------------- assembly

def test_wrong_build_is_refused_and_a_single_base_match_is_not_called_a_build(monkeypatch):
    """rs671's GRCh38 coordinates sent as GRCh37: GRCh37 has C there, GRCh38 has G.

    A one-base REF matches the other build by chance one time in four, so the message
    says so instead of telling the caller to switch builds (adversarial review P1-10)."""
    with pytest.raises(UsageError) as err:
        run(monkeypatch, "rs671_wrong_build.json")
    msg = str(err.value)
    assert "does not match the GRCh37 reference at 12:111803962, which has C" in msg
    assert "GRCh38 has G at the same number" in msg and "by chance" in msg
    assert "use assembly" not in msg


def test_grch37_input_is_lifted_for_grch38_only_datasets(monkeypatch):
    out = run(monkeypatch, "rs671_grch37.json")
    r = out.result
    assert r["variant"] == "12-112241766-G-A" and r["assembly"] == "GRCh37" and r["grch38"] == "12-111803962-G-A"
    d = by_name(out)
    nyuwa = d["NyuWa"]
    assert nyuwa["queried_as"] == "12-111803962-G-A" and nyuwa["lifted_from"] == "GRCh37 12-112241766-G-A"
    assert "lifted from GRCh37" in nyuwa["note"] and nyuwa["ac"] == 1323
    # WBBC and 1000 Genomes answer natively on GRCh37
    assert d["WBBC"]["assembly"] == "GRCh37" and (d["WBBC"]["ac"], d["WBBC"]["an"]) == (2132, 8960)
    assert d["1000 Genomes CHB"]["url"].startswith("https://grch37.rest.ensembl.org/")
    assert d["Taiwan Biobank WGS"]["queried_as"] == "12-111803962-G-A"
    assert any(s["db"] == "Ensembl assembly map" for s in out.sources)


def test_grch37_liftover_failure_skips_only_the_grch38_datasets(monkeypatch):
    def override(method, url, body):
        return (200, '{"mappings": []}') if "/map/human/" in url else None
    out = run(monkeypatch, "rs671_grch37.json", override)
    nc = not_checked(out)
    assert "could not lift" in nc["NyuWa"] and "could not lift" in nc["Taiwan Biobank WGS"]
    assert any(w.startswith("NyuWa not queried: GRCh38 only") for w in out.warnings)
    assert set(by_name(out)) == {"WBBC", "1000 Genomes CHB", "1000 Genomes CHS", "1000 Genomes CDX"}
    assert "NyuWa" not in out.result["checked"]


def test_grch37_allele_that_differs_in_grch38_is_not_sent(monkeypatch):
    fx = fixture("rs671_grch37.json")
    seq38 = exchange(fx, "rest.ensembl.org/sequence/region/human/12:111803762")
    i = 111803962 - 111803762
    changed = seq38["text"][:i] + "T" + seq38["text"][i + 1:]

    def override(method, url, body):
        return (200, changed) if url == seq38["url"] else None
    out = run(monkeypatch, "rs671_grch37.json", override)
    nc = not_checked(out)
    assert "reference is T, not G" in nc["NyuWa"] and "differs between builds" in nc["NyuWa"]
    assert "NyuWa" not in by_name(out) and "WBBC" in by_name(out)


@pytest.mark.parametrize("chrom", ["X", "MT"])
def test_datasets_without_that_chromosome_are_named_not_asked(chrom, monkeypatch):
    absent_page = exchange(fixture("scn1a_absent_grch38.json"), "NyuWa")["text"]
    calls = []

    def override(method, url, body):
        if "/sequence/region/" in url:
            return 200, "A" * 401
        if "NyuWa" in url:
            return 200, absent_page
        if "/overlap/region/" in url:
            return 200, "[]"
        raise AssertionError(f"not expected for chr{chrom}: {url}")
    monkeypatch.setattr(cf, "request", fake_request(None, override, calls))
    out = cf.lookup(chrom, 1000, "A", "G")
    nc = not_checked(out)
    assert "autosomes only" in nc["WBBC"]
    assert any(w.startswith("WBBC not queried") for w in out.warnings)
    assert not any("wbbc" in u for _, u, _ in calls)
    if chrom == "MT":
        assert "no mitochondrial variants" in nc["NyuWa"]
    else:
        assert "NyuWa" in absent_names(out)


def _synthetic(monkeypatch, seq, calls=None):
    """Reference `seq` for every window, NyuWa/WBBC/overlap answering "no record"."""
    absent_page = exchange(fixture("scn1a_absent_grch38.json"), "NyuWa")["text"]

    def override(method, url, body):
        if "/sequence/region/" in url:
            return 200, seq
        if "NyuWa" in url:
            return 200, absent_page
        if "wbbc" in url:
            return 200, ""  # WBBC's own "no record" answer: an empty body
        if "/overlap/region/" in url:
            return 200, "[]"
        raise AssertionError(f"unexpected {url}")
    monkeypatch.setattr(cf, "request", fake_request(None, override, calls))


def test_unfinished_left_alignment_is_a_warning(monkeypatch):
    """A deletion inside a poly-A run longer than the fetched window cannot be fully left-aligned."""
    _synthetic(monkeypatch, "A" * 402)
    out = cf.lookup("3", 1000, "AA", "A")
    assert any("could not be left-aligned within 200 bp" in w for w in out.warnings)
    assert out.result["datasets"] == [] and "NyuWa" in absent_names(out)


def test_mnv_gets_a_note(monkeypatch):
    _synthetic(monkeypatch, "G" * 200 + "AC" + "G" * 200)
    out = cf.lookup("3", 1000, "AC", "GT")
    assert out.result["variant"] == "3-1000-AC-GT"
    assert any("multi-nucleotide variant" in n for n in out.result["notes"])


def test_wbbc_record_without_af_uses_its_genotype_counts(monkeypatch):
    row = {"Chr": "5", "Position": "100", "Ref": "C", "Alt": "T", "ID": None, "Genotype": "RR=4470|RA=10|AA=0",
           "WBBC_AF": None, "North_AF": None, "Central_AF": "0", "South_AF": None, "Lingnan_AF": None}
    monkeypatch.setattr(cf, "request", fake_request(None, lambda m, u, b: (200, json.dumps([row]))))
    got = cf._wbbc(("5", 100, "C", "T"), "GRCh38", [])
    assert got[0]["row"]["ac"] == 10 and got[0]["row"]["an"] == 8960
    row.update(Genotype="", WBBC_AF=None)
    monkeypatch.setattr(cf, "request", fake_request(None, lambda m, u, b: (200, json.dumps([row]))))
    with pytest.raises(ValueError, match="no usable genotype counts or AF"):
        cf._wbbc(("5", 100, "C", "T"), "GRCh38", [])


# ---------------------------------------------------------------- normalisation helpers

def test_left_normalize_and_right_end():
    seq = "GGATCTTTGG"  # 1-based: A at 3, T at 4, C at 5, T at 6..8
    base = lambda p: seq[p - 1]  # noqa: E731
    assert cf.left_normalize(4, "TCTT", "T", base) == (3, "ATCT", "A")
    assert cf.left_normalize(5, "CTT", "", base) == (3, "ATCT", "A")
    assert cf.right_end(3, "ATCT", "A", base) == 7  # the deletion can slide to TCT|T -> CTT at 5..7
    assert cf.right_end(3, "A", "G", base) == 3
    assert cf.trim(10, "ATTG", "AG") == (10, "ATT", "A")


# ---------------------------------------------------------------- live

@pytest.mark.live
def test_live_rs671_in_the_chinese_cohorts():
    out = cf.lookup("12", 111803962, "G", "A")
    d = by_name(out)
    assert (d["NyuWa"]["ac"], d["NyuWa"]["an"]) == (1323, 5998), out.warnings
    assert (d["WBBC"]["ac"], d["WBBC"]["an"]) == (2138, 8960), out.warnings
    if "1000 Genomes CHB" in d:  # Ensembl can be slow; a failure there must be a named warning
        assert (d["1000 Genomes CHB"]["ac"], d["1000 Genomes CHB"]["an"]) == (33, 206)
    else:
        assert "1000 Genomes CHB" in not_checked(out)
    assert out.result["not_found_in"] == []


@pytest.mark.live
def test_live_absent_variant_is_not_af_zero():
    out = cf.lookup("2", 166042334, "G", "A")
    names = absent_names(out)
    assert "NyuWa" in names and "WBBC" in names, out.warnings
    assert all(row["name"] not in ("NyuWa", "WBBC") for row in out.result["datasets"])
