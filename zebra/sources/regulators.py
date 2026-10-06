"""Regulatory status from the agencies' own records: FDA (openFDA) and EMA (its medicines export).

Why (review E-4): Open Targets lists cannabidiol and fenfluramine for Dravet
syndrome at maxClinicalStage APPROVAL but carries no agency report for that
pair, and zebra used to tell the reader "do not read this as an approved
therapy". Epidiolex (FDA 2018, EMA 2019) and Fintepla (FDA and EMA 2020) are
approved for Dravet syndrome. Before anything is called approved or not
approved, the agencies' own records are read:

* FDA labels — openFDA `drug/label.json` (no key; https://open.fda.gov/apis/drug/label/).
  A label is only counted when it carries an NDA, BLA or ANDA application number
  (unapproved marketed products have labels too) and its INDICATIONS AND USAGE
  section names the disease. openFDA answers a search with no hit as HTTP 404
  `{"error": {"code": "NOT_FOUND"}}`, which is a clean "no label", not a failure.
* FDA approval date — openFDA `drug/drugsfda.json`, the ORIG submission's AP date
  and the products' marketing status, for the application the label names.
* EMA — the agency's published JSON export of every centrally authorised,
  refused and withdrawn human medicine
  (https://www.ema.europa.eu/en/documents/report/medicines-output-medicines_json-report_en.json,
  ~7 MB, fields checked live 2026-10-06: `international_non_proprietary_name_common_name`,
  `active_substance`, `medicine_status` (Authorised / Refused / Withdrawn / ...),
  `therapeutic_indication`, `marketing_authorisation_date`,
  `refusal_of_marketing_authorisation_date`, `medicine_url`). A refusal is a
  record too: eteplirsen (Exondys) was refused by the EU in 2018.

Nothing found is reported as "no record found here", never as "not approved":
an FDA label that words the disease differently, or a medicine authorised only
nationally in an EU member state, is not in these records.
"""

from __future__ import annotations

import re
import threading
import unicodedata
import urllib.parse
from typing import Any, Dict, Iterable, List, Optional, Sequence

from zebra.core import Outcome
from zebra.http import SourceError, get_json
from zebra.sources import record as source_record
from zebra.sources import validated_json

FDA_LABEL = "https://api.fda.gov/drug/label.json"
FDA_DRUGSFDA = "https://api.fda.gov/drug/drugsfda.json"
EMA_MEDICINES = "https://www.ema.europa.eu/en/documents/report/medicines-output-medicines_json-report_en.json"
DAILYMED = "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid={}"
DRUGSFDA_PAGE = "https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=overview.process&ApplNo={}"
APPROVED_APPLICATION = re.compile(r"^(NDA|BLA|ANDA)\d+$")
# salt, ester and vehicle words that do not name the active moiety (FENFLURAMINE HYDROCHLORIDE -> fenfluramine)
SALT_WORDS = frozenset((
    "sodium", "disodium", "trisodium", "potassium", "calcium", "magnesium", "lithium", "zinc", "hydrochloride",
    "dihydrochloride", "hcl", "hydrobromide", "bromide", "chloride", "acetate", "citrate", "phosphate", "sulfate",
    "sulphate", "mesylate", "mesilate", "maleate", "tartrate", "bitartrate", "succinate", "fumarate", "besylate",
    "besilate", "tosylate", "lactate", "gluconate", "nitrate", "oxalate", "malate", "monohydrate", "dihydrate",
    "hydrate", "anhydrous", "hemihydrate", "free", "base", "acid",
))


# anions that are the drug when the cation is lithium / sodium (LITHIUM CARBONATE, SODIUM BENZOATE):
# stripping the cation would leave "carbonate", which matched sevelamer carbonate (adversarial review)
ANION_WORDS = frozenset(("carbonate", "benzoate", "chloride", "bromide", "iodide", "citrate", "phosphate", "sulfate",
                         "sulphate", "acetate", "gluconate", "lactate", "bicarbonate", "fluoride", "nitrate",
                         "hydroxide", "oxide"))
CATION_WORDS = frozenset(("sodium", "disodium", "trisodium", "potassium", "calcium", "magnesium", "lithium", "zinc"))


def moiety(drug: str) -> str:
    """The active-moiety words of a drug name: 'NUSINERSEN SODIUM' -> 'nusinersen', 'LITHIUM CARBONATE' kept whole."""
    words = [w for w in re.split(r"[^a-z0-9\-]+", (drug or "").lower()) if w]
    kept = [w for w in words if w not in SALT_WORDS]
    if not kept or (all(w in ANION_WORDS for w in kept if w not in CATION_WORDS)
                    and any(w in CATION_WORDS for w in words)):
        return " ".join(words)
    if all(w in ANION_WORDS for w in kept):
        return " ".join(words)
    return " ".join(kept)


def _terms(disease_terms: Iterable[str]) -> List[str]:
    out: List[str] = []
    for t in disease_terms:
        t = " ".join(str(t or "").split())
        if len(t) >= 5 and t.lower() not in (x.lower() for x in out):
            out.append(t)
    return out


# A disease counts as an indication only inside a sentence that says "indicated" (positive anchor)
# and does not exclude it. Bullet items stay with the sentence that introduces them: azithromycin's
# label lists "• patients with cystic fibrosis" under "should not be used in patients ... such as
# any of the following" (adversarial review P0-1).
_NEGATED = re.compile(r"\bnot\s+(?:been\s+)?(?:indicated|approved|recommended|intended|effective)\b|"
                      r"\bnot\s+been\s+(?:established|studied|demonstrated|evaluated|shown)\b|"
                      r"\blimitations?\s+of\s+use\b|\bnot\s+for\s+(?:use|the\s+treatment)\b|"
                      r"\bshould\s+not\s+be\s+used\b|\bdo\s+not\s+use\b|\bcontraindicated\b|"
                      r"\bnot\s+supported\b|\bbut\s+not\s+(?:for|in)\b|\bexcept\b", re.I)
_INDICATED = re.compile(r"\bindicated\b", re.I)


def _fold(text: str) -> str:
    t = unicodedata.normalize("NFKD", text or "")
    t = "".join(c for c in t if not unicodedata.combining(c)).lower()
    t = re.sub(r"['’]s\b", "", t)  # Pompe's disease ~ Pompe disease
    t = t.replace("ae", "e").replace("oe", "e")  # polycythaemia (EU labels) ~ polycythemia
    t = re.sub("[\u2010-\u2015\u2212\u00ad]", "-", t)  # label hyphens: "post‑polycythemia" (U+2011)
    return t


# "post-polycythemia vera myelofibrosis", "myelofibrosis secondary to polycythemia vera": the disease
# named as the origin of another condition is not the indication
_BEFORE_SECONDARY = re.compile(r"(?:\bpost[- ]?|\bsecondary\s+to\s+|\bfollowing\s+|\bafter\s+)$")
_AFTER_SECONDARY = re.compile(r"^(?:[- ](?:associated|related|induced)\b|\s+myelofibrosis\b)")


def names_disease(text: str, terms: Sequence[str], anchor: bool = True) -> Optional[str]:
    """The first term an affirmative sentence of the text names as a phrase (word-bounded).

    `anchor` (FDA labels) also requires the sentence to say "indicated". The EMA export's
    indication field is often a bare noun phrase ("Treatment of Duchenne muscular dystrophy.":
    697 of 2,249 human rows on 2026-10-06 lack the word), so EMA text is read without it; the
    negation and secondary-condition rules apply to both.
    """
    sentences = [x for x in re.split(r"(?<=[.;])\s*(?=[A-Z•(])|\n{2,}", text or "") if x.strip()]
    for t in terms:
        pat = re.compile(r"(?<![a-z0-9])" + re.escape(_fold(t)) + r"(?![a-z0-9])")
        for sent in sentences:
            if (anchor and not _INDICATED.search(sent)) or _NEGATED.search(sent):
                continue
            folded = _fold(sent)
            for m in pat.finditer(folded):
                if _BEFORE_SECONDARY.search(folded[:m.start()]) or _AFTER_SECONDARY.match(folded[m.end():]):
                    continue
                return t
    return None


def _excerpt(text: str, term: Optional[str], width: int = 300) -> str:
    flat = " ".join((text or "").split())
    if not term:
        return flat[:width]
    i = flat.lower().find(term.lower())
    start = max(0, i - width // 3) if i >= 0 else 0
    return ("…" if start else "") + flat[start:start + width] + ("…" if start + width < len(flat) else "")


MAX_TERMS = 10  # disease names put in the label query (name + synonyms; Pompe disease is a synonym)


def _q(text: str) -> str:
    """A phrase safe inside an openFDA quoted search term."""
    return re.sub(r'["\\]+', " ", text or "").strip()


def _label_search(q: str, limit: int):
    resp = get_json(FDA_LABEL, source="openFDA drug label", params={"search": q, "limit": limit},
                    cache_ttl=7 * 86400, not_found_ttl=86400, ok_statuses=(200, 404))
    if resp.status == 404 and "NOT_FOUND" not in (resp.text or ""):
        # only openFDA's own "No matches found!" is a clean miss; any other 404 is a broken request
        raise SourceError("openFDA drug label", resp.url, 404, (resp.text or "")[:200] or "HTTP 404")
    return resp


def _active_here(base: str, label: Dict[str, Any]) -> bool:
    """The moiety is an active ingredient of the labelled product, not an excipient or a diluent."""
    o = label.get("openfda") or {}
    pat = re.compile(r"(?<![a-z0-9])" + re.escape(base) + r"(?![a-z0-9])")
    if o.get("generic_name"):
        return any(pat.search(_fold(g)) for g in o["generic_name"])
    # no harmonised block: the product data read "<BRAND> <generic name> <ACTIVE> <INACTIVE INGREDIENTS…>",
    # with the generic name the one part not written in capitals ("EVRYSDI Risdiplam RISDIPLAM MANNITOL …")
    raw = " ".join(label.get("spl_product_data_elements") or []).split()
    generic = [w for w in raw if not w.isupper()]
    head = " ".join(generic) if generic else " ".join(raw[:3])
    return bool(pat.search(_fold(head)))


def fda_label(drug: str, disease_terms: Iterable[str], limit: int = 25) -> Outcome:
    """FDA labels for `drug` (generic or substance name) whose indications name one of `disease_terms`."""
    base = moiety(drug)
    terms = _terms(disease_terms)
    if not base:
        raise ValueError("no drug name")
    # spl_product_data_elements as well: some labels carry no harmonised `openfda` block
    # (EVRYSDI / risdiplam, checked 2026-10-06), so only the label text names the moiety
    names = f'(openfda.generic_name:"{base}" openfda.substance_name:"{base}" spl_product_data_elements:"{base}")'
    # the disease goes into the query: cannabidiol has dozens of labels (cosmetics, OTC) and
    # Epidiolex is not among the first ten of a name-only search
    ind = " ".join(f'indications_and_usage:"{_q(t)}"' for t in terms[:MAX_TERMS])
    q = f"{names} AND ({ind})" if ind else names
    resp = _label_search(q, limit)
    data = validated_json(resp, "openFDA drug label") if resp.status == 200 else {}
    if not data.get("results") and ind:
        # nothing names the disease: one name-only page says whether the drug has labels at all
        resp = _label_search(names, 25)
        data = validated_json(resp, "openFDA drug label") if resp.status == 200 else {}
    labels: List[Dict[str, Any]] = []
    for r in data.get("results") or []:
        o = r.get("openfda") or {}
        if not _active_here(base, r):
            continue  # the moiety is only an excipient or diluent of this product (sodium chloride, mannitol)
        apps = [a for a in o.get("application_number") or [] if APPROVED_APPLICATION.match(a)]
        ind = " ".join(r.get("indications_and_usage") or [])
        hit = names_disease(ind, terms)
        spl = " ".join(r.get("spl_product_data_elements") or [])
        labels.append({"brand": (o.get("brand_name") or [None])[0], "generic": (o.get("generic_name") or [None])[0],
                       "application": apps[0] if apps else None,
                       "applications_all": o.get("application_number") or [],
                       "spl_name": " ".join(spl.split()[:2]) or None,
                       "label_date": r.get("effective_time"), "set_id": r.get("set_id"),
                       "url": DAILYMED.format(r.get("set_id")) if r.get("set_id") else None,
                       "names_disease": hit, "indication": _excerpt(ind, hit)})
    # the originator's NDA/BLA before a generic's ANDA: its approval date is the drug's first approval
    labels.sort(key=lambda x: (not x["names_disease"], not x["application"],
                               not str(x["application"] or "").startswith(("NDA", "BLA")), -int(x["label_date"] or 0)))
    naming = [x for x in labels if x["names_disease"]]
    result = {"drug": drug, "searched_as": base, "terms": terms, "labels_found": len(labels),
              "labels_naming_disease": naming,
              "other_labels": [{k: x[k] for k in ("brand", "application", "label_date", "url")}
                               for x in labels if x not in naming][:5]}
    note = ("no FDA label for this name" if resp.status == 404 else
            f"{len(labels)} label(s), {len(naming)} naming the disease in their indications")
    return Outcome(result, sources=[source_record("openFDA drug label", f"{base} / {', '.join(terms[:3])}", resp,
                                                  note=note)])


def fda_by_ingredient(base: str) -> Outcome:
    """Drugs@FDA applications whose products contain this active ingredient, approved ones first."""
    b = moiety(base)
    if not b:
        raise ValueError("no ingredient")
    resp = get_json(FDA_DRUGSFDA, source="openFDA Drugs@FDA",
                    params={"search": f'products.active_ingredients.name:"{b}"', "limit": 20},
                    cache_ttl=7 * 86400, not_found_ttl=86400, ok_statuses=(200, 404))
    data = validated_json(resp, "openFDA Drugs@FDA") if resp.status == 200 else {}
    apps = []
    for r in data.get("results") or []:
        app = r.get("application_number") or ""
        if not APPROVED_APPLICATION.match(app):
            continue
        ap = sorted(s.get("submission_status_date") or "" for s in r.get("submissions") or []
                    if s.get("submission_type") == "ORIG" and s.get("submission_status") == "AP")
        d = ap[0] if ap and ap[0] else None
        apps.append({"application": app, "sponsor": r.get("sponsor_name"),
                     "approved_on": f"{d[:4]}-{d[4:6]}-{d[6:8]}" if d and len(d) == 8 else d,
                     "brands": sorted({p.get("brand_name") for p in r.get("products") or [] if p.get("brand_name")}),
                     "marketing_status": sorted({p.get("marketing_status") for p in r.get("products") or []
                                                 if p.get("marketing_status")}),
                     "url": DRUGSFDA_PAGE.format(re.sub(r"\D", "", app))})
    apps.sort(key=lambda a: (a["approved_on"] is None, not a["application"].startswith(("NDA", "BLA")),
                             a["approved_on"] or ""))
    return Outcome({"ingredient": b, "applications": apps},
                   sources=[source_record("openFDA Drugs@FDA", f"ingredient {b}", resp,
                                          note=f"{len(apps)} NDA/BLA/ANDA application(s)")])


def fda_application(application: str) -> Outcome:
    """Original approval date and marketing status of one FDA application (NDA/BLA/ANDA number)."""
    app = (application or "").strip().upper()
    if not APPROVED_APPLICATION.match(app):
        raise ValueError(f"not an FDA application number: {application!r}")
    resp = get_json(FDA_DRUGSFDA, source="openFDA Drugs@FDA", params={"search": f'application_number:"{app}"', "limit": 1},
                    cache_ttl=7 * 86400, not_found_ttl=86400, ok_statuses=(200, 404))
    data = validated_json(resp, "openFDA Drugs@FDA") if resp.status == 200 else {}
    rows = data.get("results") or []
    if not rows:
        return Outcome({"application": app, "found": False},
                       sources=[source_record("openFDA Drugs@FDA", app, resp, note="no record")])
    r = rows[0]
    orig = [s for s in r.get("submissions") or [] if s.get("submission_type") == "ORIG"]
    ap = sorted(s.get("submission_status_date") or "" for s in orig if s.get("submission_status") == "AP")
    d = ap[0] if ap and ap[0] else None
    result = {"application": app, "found": True, "sponsor": r.get("sponsor_name"),
              "approved_on": f"{d[:4]}-{d[4:6]}-{d[6:8]}" if d and len(d) == 8 else d,
              "marketing_status": sorted({p.get("marketing_status") for p in r.get("products") or [] if p.get("marketing_status")}),
              "brands": sorted({p.get("brand_name") for p in r.get("products") or [] if p.get("brand_name")}),
              "url": DRUGSFDA_PAGE.format(re.sub(r"\D", "", app))}
    return Outcome(result, sources=[source_record("openFDA Drugs@FDA", app, resp, url=result["url"],
                                                  note=f"ORIG approval {result['approved_on']}")])


_EMA_LOCK = threading.Lock()  # one download of the 7-MB export, however many threads ask at once


def _ema_rows() -> Outcome:
    with _EMA_LOCK:
        resp = get_json(EMA_MEDICINES, source="EMA medicines", cache_ttl=7 * 86400, timeout=150)
    data = validated_json(resp, "EMA medicines", require="data")
    meta = data.get("meta") or {}
    return Outcome(data.get("data") or [], sources=[source_record(
        "EMA medicines (EPAR register export)", f"{meta.get('total_records')} records", resp,
        note=f"dataset {meta.get('timestamp')}")])


def _components(field: Optional[str]) -> List[str]:
    """The active moieties named in an EMA INN / active-substance field (combinations split)."""
    parts = re.split(r"\s*(?:/|,|;|\+|\band\b)\s*", _fold(field or ""))
    return [moiety(p) for p in parts if p.strip()]


def ema(drug: str, disease_terms: Iterable[str] = ()) -> Outcome:
    """EMA records (any status) whose INN or active substance is `drug`'s active moiety."""
    base = moiety(drug)
    if not base:
        raise ValueError("no drug name")
    terms = _terms(disease_terms)
    got = _ema_rows()
    rows: List[Dict[str, Any]] = []
    for r in got.result:
        if (r.get("category") or "") != "Human":
            continue
        comps = _components(r.get("international_non_proprietary_name_common_name")) + \
            _components(r.get("active_substance"))
        if base not in comps:  # the same moiety, not a word inside another one (glycerol phenylbutyrate)
            continue
        ind = r.get("therapeutic_indication") or ""
        hit = names_disease(ind, terms, anchor=False)
        rows.append({
            "medicine": r.get("name_of_medicine"), "inn": r.get("international_non_proprietary_name_common_name"),
            "status": r.get("medicine_status"), "orphan": r.get("orphan_medicine"),
            "conditional": r.get("conditional_approval"),
            "authorised_on": _dmy(r.get("marketing_authorisation_date")),
            "refused_on": _dmy(r.get("refusal_of_marketing_authorisation_date")),
            "withdrawn_on": _dmy(r.get("withdrawal_expiry_revocation_lapse_of_marketing_authorisation_date"))
            or _dmy(r.get("withdrawal_of_application_date")),
            "suspended_on": _dmy(r.get("suspension_of_marketing_authorisation_date")),
            "names_disease": hit, "indication": _excerpt(ind, hit), "url": r.get("medicine_url")})
    rows.sort(key=lambda x: (not x["names_disease"], x["status"] != "Authorised", x["medicine"] or ""))
    return Outcome({"drug": drug, "searched_as": base, "terms": terms, "records": rows[:8]}, sources=got.sources)


def _dmy(text: Optional[str]) -> Optional[str]:
    m = re.match(r"^(\d{2})/(\d{2})/(\d{4})$", str(text or "").strip())
    return f"{m.group(3)}-{m.group(2)}-{m.group(1)}" if m else None


def check(drug: str, disease_terms: Iterable[str], with_fda_date: bool = True) -> Outcome:
    """FDA and EMA records for one drug and disease, with a one-line reading per agency.

    `result = {"drug", "FDA": {...}, "EMA": {...}, "approved_in": [...], "refused_or_withdrawn_in": [...],
    "headline"}`. Each agency block is `status` ∈ approved_for_disease | record_other_indication |
    refused | withdrawn | suspended | no_record | unavailable, with the evidence row.
    """
    from zebra.sources import attempt

    out = Outcome({"drug": drug})
    terms = _terms(disease_terms)
    fda: Dict[str, Any] = {"agency": "FDA (United States)", "status": "unavailable"}
    got = attempt("openFDA drug label", lambda: fda_label(drug, terms), out.warnings)
    if got is not None:
        lab = out.add(got)
        naming = lab["labels_naming_disease"]
        if naming and not naming[0]["application"]:
            # the label names the disease but openFDA links it to no application: ask Drugs@FDA
            # by ingredient (EVRYSDI's label has no `openfda` block; NDA213535 is in Drugs@FDA)
            best = naming[0]
            by_ing = attempt("openFDA Drugs@FDA (by ingredient)", lambda: fda_by_ingredient(lab["searched_as"]),
                             out.warnings)
            apps = [a for a in (out.add(by_ing)["applications"] if by_ing is not None else []) if a["approved_on"]]
            # only the application whose brand IS this label's product: an OTC/homeopathic CBD label
            # naming "anxiety" must not borrow Epidiolex's NDA (adversarial review P0-2)
            name = (best.get("spl_name") or "").upper()
            apps = [a for a in apps if any(b and name.startswith(b.upper()) for b in a["brands"])]
            if apps:
                a = apps[0]
                fda = {"agency": "FDA (United States)", "status": "approved_for_disease",
                       "reading": f"an FDA label for {lab['searched_as']} names {best['names_disease']} in its "
                                  f"indications, and Drugs@FDA holds approved application {a['application']} "
                                  f"({', '.join(a['brands'][:2])}) for that ingredient",
                       "brand": (a["brands"] or [None])[0], "application": a["application"],
                       "approved_on": a["approved_on"], "marketing_status": a["marketing_status"],
                       "label_date": best["label_date"], "indication": best["indication"], "url": best["url"],
                       "drugsfda_url": a["url"]}
            else:
                fda = {"agency": "FDA (United States)", "status": "label_without_application",
                       "reading": f"an FDA label names {best['names_disease']}, but no approved NDA/BLA/ANDA was "
                                  f"found for {lab['searched_as']} — approval not confirmed here",
                       "label_date": best["label_date"], "url": best["url"]}
            naming = []
        if naming:
            best = naming[0]
            fda = {"agency": "FDA (United States)", "status": "approved_for_disease",
                   "reading": f"the FDA-approved label of {best['brand'] or best['generic']} ({best['application']}"
                              + (", a generic application: its date is not the drug's first approval"
                                 if str(best["application"]).startswith("ANDA") else "")
                              + f") names {best['names_disease']} in its indications",
                   "brand": best["brand"], "application": best["application"], "label_date": best["label_date"],
                   "indication": best["indication"], "url": best["url"]}
            if with_fda_date and best["application"]:
                app = attempt("openFDA Drugs@FDA", lambda: fda_application(best["application"]), out.warnings)
                if app is not None:
                    a = out.add(app)
                    if a.get("found"):
                        fda.update({"approved_on": a.get("approved_on"), "marketing_status": a.get("marketing_status"),
                                    "drugsfda_url": a.get("url")})
                        if a.get("marketing_status") and all(s == "Discontinued" for s in a["marketing_status"]):
                            fda["status"] = "withdrawn"
                            fda["reading"] += "; every product of this application is marked Discontinued"
        elif fda.get("status") in ("approved_for_disease", "label_without_application"):
            pass
        elif lab["labels_found"]:
            fda = {"agency": "FDA (United States)", "status": "record_other_indication",
                   "reading": f"{lab['labels_found']} FDA label(s) for {lab['searched_as']}, none with an approved "
                              f"application whose indications name {', '.join(terms[:2]) or 'this disease'}",
                   "labels": lab["other_labels"][:3]}
        else:
            fda = {"agency": "FDA (United States)", "status": "no_record",
                   "reading": f"no FDA label found for {lab['searched_as']} (openFDA label search)"}
    em: Dict[str, Any] = {"agency": "EMA (European Union)", "status": "unavailable"}
    got = attempt("EMA medicines", lambda: ema(drug, terms), out.warnings)
    if got is not None:
        e = out.add(got)
        recs = e["records"]
        named = [r for r in recs if r["names_disease"]]
        pick = named[0] if named else (recs[0] if recs else None)
        if pick is None:
            em = {"agency": "EMA (European Union)", "status": "no_record",
                  "reading": f"no centrally authorised, refused or withdrawn EU medicine with INN {e['searched_as']}"}
        else:
            st = (pick["status"] or "").lower()
            status = ("approved_for_disease" if st == "authorised" and named else
                      "refused" if st == "refused" else
                      "suspended" if st == "suspended" else
                      "withdrawn" if st in ("withdrawn", "revoked", "lapsed", "expired", "application withdrawn",
                                            "withdrawn from rolling review") else
                      "opinion_pending" if st.startswith("opinion") else
                      "record_other_indication" if st == "authorised" else st.replace(" ", "_") or "unknown")
            # the date of the status shown: Translarna "Expired" is dated by its expiry, not its 2014 authorisation
            when = {"refused": pick.get("refused_on"), "suspended": pick.get("suspended_on"),
                    "withdrawn": pick.get("withdrawn_on")}.get(status) or pick.get("authorised_on") \
                or pick.get("refused_on") or pick.get("withdrawn_on")
            em = {"agency": "EMA (European Union)", "status": status,
                  "reading": f"EU: {pick['medicine']} ({pick['inn']}) — {pick['status']}"
                             + (f" {when}" if when else "")
                             + ("; its indication names " + pick["names_disease"] if pick["names_disease"] else
                                "; its indication does not name this disease"),
                  "medicine": pick["medicine"], "medicine_status": pick["status"], "date": when,
                  "indication": pick["indication"], "url": pick["url"],
                  "other_records": [{k: r[k] for k in ("medicine", "status", "url")} for r in recs if r is not pick][:3]}
    approved = [b["agency"] for b in (fda, em) if b["status"] == "approved_for_disease"]
    gone = [f"{b['agency']} ({b['status']})" for b in (fda, em) if b["status"] in ("refused", "withdrawn", "suspended")]
    down = [b["agency"] for b in (fda, em) if b["status"] == "unavailable"]
    parts = []
    if approved:
        parts.append("approved for this disease per the agency record: " + ", ".join(approved))
    if gone:
        parts.append("refused/withdrawn/suspended: " + ", ".join(gone))
    if any(b["status"] == "opinion_pending" for b in (fda, em)):
        parts.append("EMA: an opinion, not yet an authorisation")
    if down:
        parts.append("could not be checked (source unavailable): " + ", ".join(down))
    if not approved and not gone:
        parts.append("no FDA label or EMA record naming this disease was found"
                     + (" in the sources that answered" if down else "")
                     + " — approval status unknown here, not 'not approved'")
    out.result = {"drug": drug, "terms": terms, "FDA": fda, "EMA": em, "approved_in": approved,
                  "refused_or_withdrawn_in": gone, "headline": "; ".join(parts)}
    return out
