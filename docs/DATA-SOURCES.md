# Data sources, access and terms

zebra queries these resources live (no key unless stated) and records each response in the case's evidence ledger (database, record, URL, retrieval time). Endpoints were verified from a test machine on 2026-10-05; services change, so `zebra doctor` re-checks reachability. Each resource keeps its own licence and terms — the links are the authority; the notes are a summary for orientation, not legal advice.

| Area | Source | Used for | Access | Terms to know |
|---|---|---|---|---|
| Phenotype | **HPO** — JAX ontology API (`ontology.jax.org/api`) and release files (GitHub `obophenotype/human-phenotype-ontology`) | term lookup/verification, disease↔phenotype annotations, offline ranking (`zebra hpo fetch`) | none | HPO licence: free with attribution, no redistribution of a modified ontology — https://hpo.jax.org/license |
| Phenotype (Chinese) | **HPO official Chinese labels** (`obophenotype/hpo-translations`, `hp-zh.babelon.tsv`) | Chinese phenotype search | none (downloaded by `zebra hpo fetch`) | as HPO |
| Phenotype ranking | **Monarch Initiative** v3 API (`api-v3.monarchinitiative.org`) | semantic-similarity disease ranking, entities, gene↔disease associations | none | see monarchinitiative.org (aggregated sources keep their licences) |
| Phenotype ranking | **PubCaseFinder** (DBCLS, `pubcasefinder.dbcls.jp`) | phenotype → disease and gene ranking | none; slow (12–50 s) | provider terms |
| Disease | **Orphanet / Orphadata API** (`api.orphadata.com`) | names (incl. Chinese), definitions, prevalence, onset, inheritance, genes, cross-references | none | CC BY 4.0 |
| Disease ids | **EBI OLS4** (Mondo) | MONDO terms, xrefs, obsolete → replacement | none | Mondo: CC BY 4.0 |
| Disease | **GeneReviews** via NCBI E-utilities (PubMed/Bookshelf) | chapter lookup | none (≤ 3 req/s; `NCBI_API_KEY` raises it) | NCBI policies; GeneReviews © University of Washington |
| China | **National rare disease lists** (2018, 121 diseases; 2023, 86 diseases) from gov.cn | list membership | bundled data with provenance | official publications |
| Variant | **Ensembl REST** (`rest.ensembl.org`, `grch37.rest.ensembl.org`) — VEP with AlphaMissense, REVEL, CADD, SpliceAI; variant_recoder; sequence; lookup | consequence on MANE, predictors, normalisation, reference sequence | none (≈ 15 req/s; intermittent 5xx, retried) | Ensembl: open; plugin data keep their licences (AlphaMissense CC BY-NC-SA 4.0 for predictions; CADD free for non-commercial use; SpliceAI scores non-commercial) |
| Variant | **gnomAD** GraphQL (`gnomad.broadinstitute.org/api`) | allele frequencies by genetic ancestry group, filtering AF, constraint (pLI, LOEUF) | none; rate-limited | see gnomad.broadinstitute.org terms |
| Variant | **ClinVar** via NCBI E-utilities | germline classification, review status, conditions | none (≤ 3 req/s) | public domain (NCBI) |
| Gene | **ClinGen** (gene–disease validity, dosage sensitivity) | validity class, MOI, HI/TS scores | none | see clinicalgenome.org terms of use |
| CNV | **ClinGen region curation list** (GRCh38 and GRCh37 TSVs) | curated recurrent regions (22q11.2, 1p36, 7q11.23, 15q11-q13, 16p11.2 …) with HI/TS scores, checked before caching | none | clinicalgenome.org terms of use |
| Repeats, SMA | **GeneReviews** chapters NBK1352 (SMA), NBK1384 (FMR1), NBK1305 (HD), NBK1165 (DM1), NBK1281 (FRDA), NBK268647 (C9orf72), with revision dates; EMQN DM1 guidelines 2012 (PMC3499739); SMN exon numbering PMC11275604; Riggs 2020 CNV standard (PMC7313390) | size bands and the sentences behind each, quoted verbatim and checked by sha256 (`tools/cnv/check_loci_sources.py`) | text read once by a person-operated browser (NCBI answers scripts with a captcha); not redistributed | GeneReviews: no redistribution of modified text |
| Gene | **PanelApp** (Genomics England) and **PanelApp Australia** | diagnostic panel membership and confidence | none | provider terms |
| Gene / protein | **UniProt** REST, **AlphaFold DB** | protein function, length, structure model | none | CC BY 4.0 |
| Literature | **Europe PMC** REST | search, abstracts, open-access full text | none | metadata free; article licences vary |
| Literature | **PubTator3**, **LitVar2** (NCBI) | entity-normalised papers; papers mentioning a variant | none | NCBI policies |
| Trials | **ClinicalTrials.gov** API v2 | trials by condition, keyword, country, status | none | public |
| Therapy | **Open Targets Platform** GraphQL | drugs and clinical candidates, tractability | none | CC0 |
| Therapy | **EMA orphan designations** (JSON report) | EU orphan designations (the FDA OOPD database blocks automated access) | none | EMA terms |
| S2F | **SpliceAI / Pangolin** Broad lookup backends | splice deltas | none; the hosted service is a courtesy — cache, do not hammer | SpliceAI weights CC BY-NC 4.0; research use |
| S2F (optional) | **AlphaGenome** API via `s2f` CLI | expression/splicing/chromatin effects by tissue | `ALPHAGENOME_API_KEY` | non-commercial only; not for clinical decision-making; outputs not to train other models |
| S2F (optional) | **Evo 2** on NVIDIA NIM via `s2f` CLI | zero-shot sequence likelihood | `NVCF_RUN_KEY` | NVIDIA terms |
| S2F (optional) | **GPN-MSA** precomputed scores via `s2f` CLI | conservation/constraint (hg38 SNVs) | none; needs `tabix` | provider licence |
| Variant (China) | **NyuWa** (`bigdata.ibp.ac.cn`), **WBBC** (Westlake), **1000 Genomes** CHB/CHS/CDX (via Ensembl), **Taiwan Biobank** (`taiwanview.twbiobank.org.tw`) | allele frequencies in Chinese cohorts; a cohort without a record is "not found in", never AF 0 | none; NyuWa answers over plain HTTP only (the result says so) | provider terms |
| Variant | **MyVariant.info** (BioThings) | whole-exome frequency prefilter for `vcf triage --prefilter myvariant` (gnomAD v2.1.1 exomes / v3 genomes as served, ClinVar) | none; ≤ 1,000 ids per request | sends every quality-passing variant position — the mod asks first |
| Variant (mtDNA) | **MITOMAP** via the MITOMASTER web service | disease associations of mtDNA variants (curation status not returned) | none; the MITOMAP pages themselves refuse scripts | MITOMAP terms |
| Function | **MaveDB** API v1 (`api.mavedb.org`) | multiplexed assay scores, labelled only by the dataset's own calibration (Brnich 2019 → PS3/BS3 strength) | none | per-dataset licences |
| Expression | **GTEx Portal** API v2 (`gtexportal.org`) | median expression per tissue with UBERON/EFO ids (AlphaGenome tissue, RNA-test tissue) | none | GTEx open-access terms |
| S2F | **NCBI BLAST** URL API | genomic uniqueness of antisense target windows (`aso_screen --uniqueness`, opt-in) | none; one request per 10 s, 30–700 s per search | NCBI policies |
| Regulators | **openFDA** (drug labels), **Drugs@FDA**, **EMA medicines** export | FDA/EMA approval status, withdrawals and refusals, from the agencies' own records | none (openFDA: 240 requests/min without a key) | public |
| China | **NMPA / CDE** official documents (three rare-disease drug lists, CDE urgently-needed overseas drug lists, annual review reports 2018–2025, NMPA news) | `approved_in_china` / `named_not_approved` / `not_in_bundled_list` for 150 drugs | bundled data with URL, sha256 and row counts per source | not a complete approval register; no exact approval dates |
| China | **National Reimbursement Drug List 2025** (医保发〔2025〕33号) and the commercial innovative-drug list | reimbursement entry with its restriction text verbatim | bundled data with provenance | Western-medicine part only |
| China | **国家罕见病诊疗协作网** hospital lists (国卫办医函〔2019〕157号; 国卫办医政函〔2024〕64号) | hospitals by province, lead hospitals | bundled data with provenance | official publication |

Not used, and why:

- **OMIM** API needs a personal key and forbids redistribution; zebra links to OMIM ids returned by other sources. If you hold a key, set `OMIM_API_KEY` (reserved for a future release).
- **DECIPHER** needs a key.
- **FDA OOPD** blocks automated access; EMA's list is used instead.
- **ChiCTR** answers every scripted request with HTTP 405 (an Aliyun WAF page); WHO ICTRP's robots.txt disallows automated access. zebra gives the ChiCTR search address instead of claiming there are no Chinese trials.
- **ChinaMAP** needs a login and a captcha; **CMDB** needs a token.
- **UCSC BLAT** is behind a Cloudflare challenge; **SpliceVault** is a Shiny app with no API.
- **LIRICAL / Exomiser** have no web API (Java + ~37 GB data). zebra's offline ranker is a transparent baseline, not a replacement; run them locally if you can.

Privacy: zebra sends identifiers of biology (HPO ids, gene symbols, variants, disease names) to these services — never names, dates of birth, record numbers or raw genome files. The mod's privacy gate enforces this for every outgoing tool call.
