// The tools zebra-mod gives the model. Each one is a thin, typed door onto a
// `zebra` CLI command (bin/zebra, Python standard library only), so the CLI
// stays the single implementation: the same answers come back whether the
// model calls the tool, a skill runs the command in Bash, or a researcher
// scripts it.

export type ToolDef = {
  name: string
  description: string
  inputSchema: Record<string, unknown>
  argv: (input: Record<string, unknown>, casePath: string | null) => string[]
  deferred?: boolean
  touchesCase?: boolean
  timeoutMs?: number
}

const str = (v: unknown): string | undefined => (typeof v === 'string' && v.trim() ? v.trim() : undefined)
const num = (v: unknown): string | undefined => (typeof v === 'number' && Number.isFinite(v) ? String(v) : undefined)
const list = (v: unknown): string[] =>
  Array.isArray(v) ? v.filter((x): x is string => typeof x === 'string' && x.trim() !== '').map(x => x.trim()) : []
const flag = (name: string, v: string | undefined): string[] => (v === undefined ? [] : [name, v])

// What each statistics method accepts, so a params key cannot smuggle --help or --case.
const STATS_FLAGS: Record<string, string[]> = {
  segregation: ['ad-meioses', 'xlr-male-meioses', 'ar-affected-sibs', 'ar-unaffected-sibs', 'nonsegregations',
    'unaffected-carriers', 'full-penetrance'],
  maxaf: ['prevalence', 'allelic', 'genetic', 'penetrance', 'inheritance', 'an', 'faf95'],
  carrier: ['prevalence', 'allele-freqs'],
  recurrence: ['mode', 'penetrance', 'mosaic', 'prior', 'unaffected-sons', 'affected-sons'],
  fisher: ['a', 'b', 'c', 'd'],
  burden: ['case-carriers', 'case-n', 'control-carriers', 'control-n'],
  denovo: ['observed', 'trios', 'mu'],
  km: ['csv', 'time-col', 'event-col', 'group-col', 'event-coding'],
  nof1: ['effect', 'sd-diff', 'alpha', 'power', 'treatment', 'control'],
}

const HPO = { type: 'string', pattern: '^HP:\\d{7}$' }
const ASSEMBLY = { type: 'string', enum: ['GRCh38', 'GRCh37'], description: 'Genome build of genomic coordinates; ask when unknown.' }

export const TOOLS: ToolDef[] = [
  {
    name: 'case_status',
    description:
      'Read the active zebra case: phenotypes (HPO), variants, hypotheses, therapy leads, open questions and how many evidence rows back them. Call before relying on what the case holds.',
    inputSchema: { type: 'object', properties: {} },
    argv: (_i, casePath) => {
      if (!casePath) throw new Error('no active case (the person can run /zebra new <dir> or /zebra case <dir>)')
      return ['case', 'summary', casePath]
    },
  },
  {
    name: 'case_update',
    description:
      'Record findings in the active case. identifiers FIRST when records name the patient: the name (every spelling: 汉字, pinyin), date of birth, record/ID numbers — the privacy gate then keeps them out of every outgoing call; they are never echoed back. profile: role (family/patient/clinician/researcher), language, proband sex/age, consanguinity. phenotypes: HPO terms present or excluded (labels are verified against HPO; never invent an id — find it with hpo_search). variants, hypotheses (with ORPHA/OMIM/MONDO ids and ledger evidence ids for and against), ACMG readings of recorded variants (codes in, class computed by zebra), therapy leads, tests already done (CMA, panel, exome…, with result), family members (relation, affected, genotype — never names), timeline events, questions for the care team, and removals. Several at once.',
    inputSchema: {
      type: 'object',
      properties: {
        profile: {
          type: 'object',
          description: 'who the case is for and how to write to them; proband basics (no names or birth dates)',
          properties: {
            title: { type: 'string' },
            role: { type: 'string', enum: ['family', 'patient', 'clinician', 'researcher'] },
            language: { type: 'string', description: 'zh or en' },
            sex: { type: 'string', enum: ['female', 'male', 'unknown'] },
            age: { type: 'string', description: 'age or age band, e.g. "2y11m"' },
            ancestry: { type: 'string' },
            consanguinity: { type: 'boolean' },
          },
        },
        phenotypes: {
          type: 'array',
          items: {
            type: 'object',
            properties: {
              id: HPO,
              status: { type: 'string', enum: ['present', 'excluded'] },
              onset: { type: 'string' },
              source: { type: 'string', description: 'where it is documented, e.g. records/neuro-2023.pdf p2' },
              note: { type: 'string' },
            },
            required: ['id'],
          },
        },
        variants: {
          type: 'array',
          items: {
            type: 'object',
            properties: {
              gene: { type: 'string' },
              hgvs_c: { type: 'string', description: 'with transcript, e.g. NM_001165963.4:c.2134C>T' },
              hgvs_g: { type: 'string' },
              hgvs_p: { type: 'string' },
              vcf: { type: 'string', description: 'chrom-pos-ref-alt, e.g. 2-166001-C-T' },
              assembly: ASSEMBLY,
              zygosity: { type: 'string', enum: ['het', 'hom', 'hemi', 'mosaic', 'unknown'] },
              inheritance: { type: 'string', enum: ['de_novo', 'maternal', 'paternal', 'biparental', 'unknown'] },
              classification_lab: { type: 'string', description: 'as written on the lab report' },
              source: { type: 'string' },
              description: { type: 'string' },
              kind: {
                type: 'string',
                enum: ['small', 'cnv', 'exon_cnv', 'copy_number', 'repeat_expansion'],
                description: 'small (default) for an SNV/indel; cnv for a CMA/CNV-seq result; exon_cnv for an exon-level deletion or duplication; copy_number for SMN1-type dosage; repeat_expansion for a repeat expansion',
              },
              region: { type: 'string', description: 'cnv: chr15:23123715-28193120' },
              iscn: { type: 'string', description: 'cnv: the ISCN string as reported, e.g. arr[GRCh38] 22q11.21(18648855_21800471)x1' },
              cnv_type: { type: 'string', enum: ['loss', 'gain'] },
              copy_number: { type: 'number', description: 'copies reported (SMN1 exon 7 = 0, 1, 2 …)' },
              exons: { type: 'string', description: 'exon_cnv: 45-50, or a single exon number' },
              genes: { type: 'array', items: { type: 'string' }, description: 'genes the finding covers, as a tool returned them' },
              motif: { type: 'string', description: 'repeat: the repeat unit, e.g. CGG' },
              repeat_count: { type: 'number', description: 'repeat: the number of units reported' },
              method: { type: 'string', description: 'how it was measured: CMA, CNV-seq, MLPA, ddPCR, repeat-primed PCR …' },
            },
          },
        },
        hypotheses: {
          type: 'array',
          items: {
            type: 'object',
            properties: {
              disease: { type: 'string' },
              status: { type: 'string', enum: ['leading', 'considered', 'excluded', 'confirmed'] },
              ids: { type: 'array', items: { type: 'string' }, description: 'ORPHA:33069, OMIM:607208, MONDO:0100135 as tools returned them' },
              support: { type: 'array', items: { type: 'string' }, description: 'ledger ids, e.g. E4' },
              against: { type: 'array', items: { type: 'string' } },
              note: { type: 'string' },
            },
            required: ['disease'],
          },
        },
        leads: {
          type: 'array',
          items: {
            type: 'object',
            properties: {
              name: { type: 'string' },
              kind: { type: 'string', enum: ['approved', 'trial', 'repurposing', 'n-of-1', 'supportive', 'other'] },
              status: { type: 'string' },
              evidence: { type: 'array', items: { type: 'string' } },
              note: { type: 'string' },
            },
            required: ['name', 'kind'],
          },
        },
        acmg: {
          type: 'array',
          description: 'store a research-grade ACMG reading on a recorded variant: the codes you justified (classified by zebra, not by you)',
          items: {
            type: 'object',
            properties: {
              variant_id: { type: 'string', description: 'v1, v2 ... as case_status lists them' },
              codes: { type: 'array', items: { type: 'string' } },
              note: { type: 'string' },
            },
            required: ['variant_id', 'codes'],
          },
        },
        questions: { type: 'array', items: { type: 'string' } },
        identifiers: {
          type: 'array',
          items: { type: 'string' },
          description: 'protected identifiers to add (name in every spelling, date of birth, record/ID numbers); kept out of every outgoing call, never shown back',
        },
        tests: {
          type: 'array',
          description: 'tests already done and their reported result',
          items: {
            type: 'object',
            properties: {
              type: { type: 'string', description: 'CMA, karyotype, gene panel, trio exome, genome, MLPA, repeat sizing, metabolic screen, MRI, EEG…' },
              date: { type: 'string' },
              result: { type: 'string', description: 'as reported, e.g. "normal", "VUS SCN1A c.…", "arr[GRCh38] 22q11.21(…)x1"' },
              lab: { type: 'string' },
              method: { type: 'string' },
              source: { type: 'string', description: 'which record file it came from' },
              note: { type: 'string' },
            },
            required: ['type'],
          },
        },
        family: {
          type: 'array',
          description: 'relatives as the records describe them — relation, affected or not, genotype if tested; never names',
          items: {
            type: 'object',
            properties: {
              relation: { type: 'string', enum: ['mother', 'father', 'sibling', 'brother', 'sister', 'half-sibling', 'child', 'son', 'daughter', 'maternal grandmother', 'maternal grandfather', 'paternal grandmother', 'paternal grandfather', 'maternal aunt', 'maternal uncle', 'paternal aunt', 'paternal uncle', 'cousin', 'twin', 'other'] },
              sex: { type: 'string' },
              affected: { type: ['boolean', 'string'], description: 'true, false or "unknown"' },
              status: { type: 'string', description: 'alive, deceased, …' },
              genotype: { type: 'string', description: 'e.g. "het for v1", "not carrier of v1", "not tested"' },
              age: { type: 'string' },
              note: { type: 'string' },
              source: { type: 'string' },
            },
          },
        },
        timeline: {
          type: 'array',
          items: {
            type: 'object',
            properties: { date: { type: 'string' }, event: { type: 'string' }, source: { type: 'string' } },
            required: ['event'],
          },
        },
        remove: {
          type: 'array',
          items: {
            type: 'object',
            properties: {
              kind: { type: 'string', enum: ['phenotype', 'variant', 'hypothesis', 'lead', 'test', 'relative'] },
              id: { type: 'string' },
            },
            required: ['kind', 'id'],
          },
        },
      },
    },
    argv: (input, casePath) => {
      if (!casePath) throw new Error('no active case (the person can run /zebra new <dir> or /zebra case <dir>)')
      return ['case', 'apply', casePath, '--ops', JSON.stringify(input)]
    },
    touchesCase: true,
  },
  {
    name: 'hpo_search',
    description:
      'Find Human Phenotype Ontology terms for a clinical phrase. Chinese: official Chinese labels once the HPO release is fetched, and the lay phrases families use (走路晚, 抽风, 不会说话, 发热惊厥, 头围小, 听力下降) with or without it — also inside a short sentence ("孩子走路晚"); a negation (无明显抽搐) is not read as the feature and is warned about. English works either way. Returns verified ids with the label, the Chinese label where there is one, and what matched (label, synonym or lay phrase). Use it for every phenotype before recording one.',
    inputSchema: {
      type: 'object',
      properties: { text: { type: 'string' }, limit: { type: 'number', default: 8 } },
      required: ['text'],
    },
    argv: i => ['hpo', 'search', str(i.text) ?? '', ...flag('--limit', num(i.limit))],
  },
  {
    name: 'phenotype_rank',
    description:
      'Phenotype-driven differential diagnosis: rank diseases and genes for a set of HPO terms (present, and excluded). Sources: local (offline Resnik best-match over HPO annotations), monarch (Monarch semantic similarity), pubcasefinder (PubCaseFinder). Each source ranks on its own; agreement across them is the signal. Excluded terms are checked against each disease and listed in `excluded_hits` (a disease that usually has a feature the patient lacks) but do not lower the local score — on the phenopacket-store benchmark that ranked better (docs/BENCHMARK.md); read the hits and weigh them yourself, or set excluded_weight 1 for the 0.1.0 penalty. Ties at the top are reported. Measured accuracy (held-out cases whose own paper is not an HPO annotation source): correct disease in the local top 10 for 28%; all held-out cases 70% (an upper bound). from_case uses the active case\'s phenotypes.',
    inputSchema: {
      type: 'object',
      properties: {
        present: { type: 'array', items: HPO },
        excluded: { type: 'array', items: HPO },
        from_case: { type: 'boolean' },
        sources: { type: 'array', items: { type: 'string', enum: ['local', 'monarch', 'pubcasefinder'] } },
        top: { type: 'number', default: 15 },
        local_method: { type: 'string', enum: ['resnik', 'lr'], description: 'local scoring (default resnik)' },
        excluded_weight: { type: 'number', description: 'local penalty for excluded terms present in a disease: 0 (default, flagged only) to 10; 1 = zebra 0.1.0' },
      },
    },
    argv: i => {
      const out = ['phenotype', 'rank']
      const present = list(i.present)
      const excluded = list(i.excluded)
      if (present.length) out.push('--present', ...present)
      if (excluded.length) out.push('--exclude', ...excluded)
      if (i.from_case === true) out.push('--from-case')
      const sources = list(i.sources)
      if (sources.length) out.push('--sources', sources.join(','))
      return [...out, ...flag('--top', num(i.top)), ...flag('--local-method', str(i.local_method)),
        ...flag('--excluded-weight', num(i.excluded_weight))]
    },
    timeoutMs: 180_000,
  },
  {
    name: 'gene_card',
    description:
      'Everything that matters about one gene for rare disease: HGNC identity (aliases resolved), causal disease associations with inheritance (Monarch: OMIM and ClinGen), ClinGen gene–disease validity classes and dosage sensitivity (haploinsufficiency / triplosensitivity), gnomAD constraint (pLI, LOEUF, missense Z), PanelApp (England and Australia) panels with confidence, protein function (UniProt) and AlphaFold model. Use symbols as HGNC spells them.',
    inputSchema: { type: 'object', properties: { symbol: { type: 'string' } }, required: ['symbol'] },
    argv: i => ['gene', str(i.symbol) ?? ''],
  },
  {
    name: 'variant_card',
    description:
      'Annotate one variant: normalized forms (HGVS, GRCh38/37 coordinates), consequence on the MANE transcript (Ensembl VEP), population frequency (gnomAD, by genetic ancestry group, with filtering AF), ClinVar classification with review status, in-silico predictors (REVEL, AlphaMissense, CADD, SpliceAI where available) and literature mentions (LitVar). Input: HGVS with transcript (NM_...:c.), rsID, or chrom-pos-ref-alt with assembly.',
    inputSchema: {
      type: 'object',
      properties: {
        variant: { type: 'string', description: 'HGVS, chrom-pos-ref-alt, rsID, or an mtDNA change such as "m.3243A>G 35%"' },
        assembly: ASSEMBLY,
        gene: { type: 'string' },
        heteroplasmy: { type: 'number', description: 'mtDNA only: the reported heteroplasmy, in percent' },
      },
      required: ['variant'],
    },
    argv: i => [
      'variant', str(i.variant) ?? '', ...flag('--assembly', str(i.assembly)), ...flag('--gene', str(i.gene)),
      ...flag('--heteroplasmy', num(i.heteroplasmy)),
    ],
    timeoutMs: 180_000,
  },
  {
    name: 'disease_card',
    description:
      'One rare disease: identifiers across ORPHA, OMIM, MONDO, ICD; definition, prevalence, inheritance, age of onset, associated genes (Orphanet / Monarch), GeneReviews chapter, whether it is on China\'s national rare disease lists. Input: a name or an id (ORPHA:33069, OMIM:607208, MONDO:0100135).',
    inputSchema: { type: 'object', properties: { query: { type: 'string' } }, required: ['query'] },
    argv: i => ['disease', str(i.query) ?? ''],
    timeoutMs: 180_000,
  },
  {
    name: 'acmg',
    description:
      'ACMG/AMP classification arithmetic. classify: give the evidence codes you have justified (e.g. PVS1, PS2, PM2_Supporting, PP3_Strong, BS1) and get the points (Tavtigian 2020) and the 2015 combining-rule class, with rule warnings. suggest: for a variant, the codes that follow from numbers alone (frequency, REVEL, SpliceAI) with their ClinGen thresholds — judgement codes (PVS1, PS3, PM3, PP1, PS4...) stay yours to justify.',
    inputSchema: {
      type: 'object',
      properties: {
        mode: { type: 'string', enum: ['classify', 'suggest'] },
        codes: { type: 'array', items: { type: 'string' } },
        variant: { type: 'string' },
        assembly: ASSEMBLY,
        inheritance: { type: 'string', enum: ['AD', 'AR', 'XLD', 'XLR', 'unknown'] },
        prevalence: { type: 'number', description: 'suggest: disease prevalence (e.g. 0.0001) for the Whiffin maximum credible AF behind BS1' },
        allelic: { type: 'number', description: 'suggest: maximum allelic contribution of one variant (0-1]' },
        genetic: { type: 'number', description: 'suggest: maximum genetic contribution of this gene (0-1]' },
        penetrance: { type: 'number', description: 'suggest: penetrance (0-1]' },
        inheritance_mode: { type: 'string', enum: ['monoallelic', 'biallelic'], description: 'suggest: only for XLR, XLD or unknown inheritance; AD/AR choose it themselves' },
        pm2_max_af: { type: 'number', description: 'suggest: a gene-specific PM2 ceiling from a ClinGen VCEP' },
        heteroplasmy: { type: 'number', description: 'suggest, mtDNA: heteroplasmy in percent' },
        no_splice_lookup: { type: 'boolean', description: 'suggest: skip the SpliceAI/Pangolin lookup (±4,999 nt) and use VEP\'s precomputed SpliceAI' },
      },
      required: ['mode'],
    },
    argv: i => {
      if (i.mode === 'suggest') {
        return [
          'acmg', 'suggest', str(i.variant) ?? '', ...flag('--assembly', str(i.assembly)), ...flag('--inheritance', str(i.inheritance)),
          ...flag('--prevalence', num(i.prevalence)), ...flag('--allelic', num(i.allelic)), ...flag('--genetic', num(i.genetic)),
          ...flag('--penetrance', num(i.penetrance)), ...flag('--inheritance-mode', str(i.inheritance_mode)),
          ...flag('--pm2-max-af', num(i.pm2_max_af)), ...flag('--heteroplasmy', num(i.heteroplasmy)),
          ...(i.no_splice_lookup === true ? ['--no-splice-lookup'] : []),
        ]
      }
      return ['acmg', 'classify', ...list(i.codes)]
    },
    timeoutMs: 240_000,
  },
  {
    name: 's2f_predict',
    description:
      'Sequence-to-function predictions for one variant, each model reported separately with its claim ceiling (never averaged): spliceai and pangolin (Broad lookup service; GRCh38 and GRCh37) — the default; and, through the s2f CLI when installed with keys, gpn_msa (conservation, hg38 SNVs, ~1 min), alphagenome (expression/splicing/chromatin in one tissue: needs ontology, e.g. UBERON:0000955 brain; ~5 s) and evo2 (zero-shot likelihood; 3–10 min). A model that cannot run is not_run with the reason. Use for splice-region, deep intronic, UTR, promoter and other non-coding variants, and to test a mechanism.',
    inputSchema: {
      type: 'object',
      properties: {
        variant: { type: 'string', description: 'chrom-pos-ref-alt, transcript HGVS or rsID' },
        assembly: ASSEMBLY,
        models: { type: 'array', items: { type: 'string', enum: ['spliceai', 'pangolin', 'gpn_msa', 'alphagenome', 'evo2'] }, description: 'default spliceai + pangolin' },
        ontology: { type: 'string', description: 'tissue/cell CURIE for alphagenome (UBERON:/CL:), chosen for the disease' },
        distance: { type: 'number', description: 'SpliceAI window around the variant (default 500)' },
        timeout: { type: 'number', description: 'seconds per s2f model run (default 540)' },
      },
      required: ['variant'],
    },
    argv: i => {
      const models = list(i.models)
      return [
        's2f', 'predict', str(i.variant) ?? '',
        ...flag('--assembly', str(i.assembly)),
        '--models', (models.length ? models : ['spliceai', 'pangolin']).join(','),
        ...flag('--ontology', str(i.ontology)),
        ...flag('--distance', num(i.distance)),
        ...flag('--timeout', num(i.timeout) ?? '540'),
      ]
    },
    timeoutMs: 600_000,
    deferred: true,
  },
  {
    name: 'therapy_landscape',
    description:
      'Therapy landscape for a disease or gene (data only, no judgement): approved and investigational drugs and clinical candidates with stage and mechanism (Open Targets / ChEMBL), target tractability for a gene, and EU orphan designations (EMA; US designations are not checked). Input: a disease name or MONDO/EFO id, or a gene symbol (OMIM:/ORPHA: ids are not accepted — use the name). Mechanism fit and N-of-1 routes are for the zebra-therapy skill to judge.',
    inputSchema: { type: 'object', properties: { query: { type: 'string', description: 'disease name or MONDO/EFO id, or gene symbol' } }, required: ['query'] },
    argv: i => ['therapy', str(i.query) ?? ''],
    timeoutMs: 180_000,
    deferred: true,
  },
  {
    name: 'trials_search',
    description:
      'Clinical trials from ClinicalTrials.gov (v2 API): by condition and optionally a keyword (gene, drug, modality), country and status (default RECRUITING; ANY for all). With a country, only trials with a site there are returned, and those sites are listed. Returns NCT ids, phase, status, interventions, ages, eligibility criteria, contacts, sites and URL, and flags a trial whose status looks stale or inconsistent. ChiCTR (the Chinese registry) is not covered: for China, say so and point to it.',
    inputSchema: {
      type: 'object',
      properties: {
        condition: { type: 'string' },
        term: { type: 'string', description: 'extra keyword: gene, drug, modality' },
        country: { type: 'string', description: 'e.g. China, United States; only trials with a site there, and those sites are listed' },
        full_eligibility: { type: 'boolean', description: 'the whole eligibility text instead of its first 600 characters' },
        keep_unrelated: { type: 'boolean', description: 'keep trials whose conditions do not name this disease (normally moved to result.filtered)' },
        status: { type: 'string', enum: ['RECRUITING', 'NOT_YET_RECRUITING', 'ACTIVE_NOT_RECRUITING', 'COMPLETED', 'ANY'] },
        limit: { type: 'number', default: 20 },
      },
      required: ['condition'],
    },
    argv: i => [
      'trials', str(i.condition) ?? '',
      ...flag('--term', str(i.term)), ...flag('--country', str(i.country)), ...flag('--status', str(i.status)),
      ...(i.full_eligibility === true ? ['--full-eligibility'] : []), ...flag('--limit', num(i.limit)),
      ...(i.keep_unrelated === true ? ['--keep-unrelated'] : []),
    ],
    deferred: true,
  },
  {
    name: 'literature_search',
    description:
      'Search the literature (Europe PMC; with gene/variant, PubTator3 and LitVar for papers that mention them). Returns PMIDs, titles, years, journals and open-access links; read a paper before citing what it says.',
    inputSchema: {
      type: 'object',
      properties: {
        query: { type: 'string' },
        gene: { type: 'string' },
        variant: { type: 'string', description: 'rsID or HGVS protein/cDNA change' },
        sort: { type: 'string', enum: ['relevance', 'date', 'cited'] },
        abstract: { type: 'array', items: { type: 'string' }, description: 'PMIDs whose abstracts to return' },
        limit: { type: 'number', default: 15 },
      },
    },
    argv: i => [
      'lit', ...(str(i.query) ? [str(i.query) as string] : []),
      ...flag('--gene', str(i.gene)), ...flag('--variant', str(i.variant)), ...flag('--sort', str(i.sort)),
      ...(list(i.abstract).length ? ['--abstract', ...list(i.abstract)] : []), ...flag('--limit', num(i.limit)),
    ],
    deferred: true,
  },
  {
    name: 'rare_stats',
    description:
      'Classical statistics for rare-disease genetics (no network): segregation (cosegregation LR and PP1 strength), maxaf (Whiffin maximum credible allele frequency), carrier (Hardy–Weinberg carrier frequency from prevalence or allele frequencies), recurrence, fisher, burden (case/control collapsing test), denovo (de novo enrichment, Poisson), km (Kaplan–Meier from a CSV), nof1 (N-of-1 trial design or analysis). params are the CLI flags of `zebra stats <method>` without dashes.',
    inputSchema: {
      type: 'object',
      properties: {
        method: { type: 'string', enum: ['segregation', 'maxaf', 'carrier', 'recurrence', 'fisher', 'burden', 'denovo', 'km', 'nof1'] },
        params: { type: 'object', additionalProperties: true },
      },
      required: ['method'],
    },
    argv: i => {
      const method = str(i.method) ?? ''
      const allowed = STATS_FLAGS[method]
      if (!allowed) throw new Error(`unknown method ${JSON.stringify(method)}`)
      const out = ['stats', method]
      const params = (i.params && typeof i.params === 'object' ? i.params : {}) as Record<string, unknown>
      for (const [k, v] of Object.entries(params)) {
        const flag = k.replace(/_/g, '-')
        if (!allowed.includes(flag)) {
          throw new Error(`${method} does not take ${JSON.stringify(k)}; it takes ${allowed.join(', ')}`)
        }
        const key = `--${flag}`
        if (Array.isArray(v)) out.push(key, ...v.map(String))
        else if (typeof v === 'boolean') {
          if (v) out.push(key)
        } else if (v !== null && v !== undefined) out.push(key, String(v))
      }
      return out
    },
    deferred: true,
  },
  {
    name: 'edit_check',
    description:
      'Base-editing feasibility screen for an SNV: can an adenine or cytosine base editor revert the patient allele to reference, and which SpCas9 protospacers (PAM NGG, or NG for relaxed-PAM variants) put the base in the editing window, with bystander bases (annotate_bystanders checks their coding effect via VEP). Reference sequence from Ensembl. A screen for researchers, not a guide design; delivery to the tissue is the hard part.',
    inputSchema: {
      type: 'object',
      properties: {
        variant: { type: 'string', description: 'chrom-pos-ref-alt (SNV), transcript HGVS or rsID' },
        assembly: ASSEMBLY,
        pam: { type: 'string', enum: ['NGG', 'NG'], description: 'NG when no NGG protospacer exists' },
        window: { type: 'string', description: 'editing window in protospacer positions, default "4-8"' },
        annotate_bystanders: { type: 'boolean' },
      },
      required: ['variant'],
    },
    argv: i => [
      'edit', str(i.variant) ?? '',
      ...flag('--assembly', str(i.assembly)), ...flag('--pam', str(i.pam)), ...flag('--window', str(i.window)),
      ...(i.annotate_bystanders === true ? ['--annotate-bystanders'] : []),
    ],
    deferred: true,
  },
  {
    name: 'cnv_interpret',
    description:
      'Read the report forms that are not a single sequence variant. A CNV or microarray result (region or ISCN, with build, mosaicism such as x1~2 or 30%, de novo/maternal/paternal) → the genes it spans, every one checked against ClinGen gene dosage, and the ClinGen curated regions it overlaps (22q11.2, 1p36, Williams, PWS/AS, 16p11.2 …) with HI/TS scores, coverage and the ACMG/ClinGen 2020 section-2 row they point to; section 3 by gene count for losses and gains separately — inputs, never a classification (sections 4-5 need the case). An exon-level deletion or duplication (DMD exon 45-50 deletion, NM_004006.3:c.6439-?_7309+?del, intronic breakpoints) → exons, coding bases, whether the frame is kept, and — only when it is not — which additional exon skip restores it. SMN1/SMN2 copy number (clinical exon numbering 1, 2a, 2b, 3-8) read with the SMN2 count, framed as population data, not a prognosis. Repeat expansions (FMR1, HTT, DMPK, FXN, C9orf72) placed in their published size bands with sex-specific wording and what else is needed (methylation, AGG interruptions, the other allele) — a referenced reading, not a classification. Coordinates from Ensembl.',
    inputSchema: {
      type: 'object',
      properties: {
        result: {
          type: 'string',
          description: 'the finding as the report gives it: "chr15:23123715-28193120 loss" | "arr[GRCh38] 22q11.21(18648855_21800471)x1" | "DMD exon 45-50 deletion" | "NM_004006.3:c.6439-?_7309+?del" | "SMN1 exon 7 copy number 0" | "FMR1 CGG 230"',
        },
        assembly: ASSEMBLY,
        gene: { type: 'string', description: 'when the string does not name one' },
        copies: { type: 'number', description: 'copy number the report gives, when it is not in the string' },
        inheritance: { type: 'string', enum: ['de_novo', 'maternal', 'paternal', 'biparental', 'unknown'] },
        method: { type: 'string', description: 'CMA, CNV-seq, MLPA, ddPCR, repeat-primed PCR …' },
        record: { type: 'boolean', description: 'also record it in the active case (a write: asked about like any case change)' },
        sex: { type: 'string', enum: ['male', 'female'], description: 'decides what an X/Y copy number means (default: the case profile)' },
        smn2_copies: { type: 'number', description: 'SMN2 copy number, read together with SMN1' },
        related: { type: 'array', items: { type: 'string' }, description: 'other findings in the same report (the second allele of a repeat, a second CNV)' },
      },
      required: ['result'],
    },
    argv: i => [
      'cnv', str(i.result) ?? '',
      ...flag('--assembly', str(i.assembly)), ...flag('--gene', str(i.gene)), ...flag('--copies', num(i.copies)),
      ...flag('--inheritance', str(i.inheritance)), ...flag('--method', str(i.method)),
      ...flag('--sex', str(i.sex)), ...flag('--smn2-copies', num(i.smn2_copies)),
      ...(list(i.related).length ? ['--related', ...list(i.related)] : []),
      ...(i.record === true ? ['--record'] : []),
    ],
    touchesCase: true,
    timeoutMs: 180_000,
  },
  {
    name: 'china_rare',
    description:
      "China's national rare disease lists (第一批罕见病目录 2018, 第二批 2023): is a disease on them, its Chinese name and list number. status: on_list (the published entry, or a named member of a listed group); qualified (only a subtype or form is listed — the result says which); possible (closest entries to verify, e.g. an acronym or a shared stretch of text — not a match); not_found. With query \"hospitals\": the 国家罕见病诊疗协作网 hospitals (NHC 2024 list, 419), by province with the lead hospitals. Input: Chinese or English name.",
    inputSchema: {
      type: 'object',
      properties: {
        query: { type: 'string', description: 'a disease name, or "hospitals" for the collaboration-network hospitals' },
        province: { type: 'string', description: 'with "hospitals": e.g. 浙江 or 浙江省' },
      },
      required: ['query'],
    },
    argv: i => ['china', str(i.query) ?? '', ...flag('--province', str(i.province))],
    deferred: true,
  },
  {
    name: 'case_recheck',
    description:
      'Ask the active case\'s questions again and say what changed since the last recheck: ClinVar classification and stars of each recorded sequence variant, gnomAD frequency, ClinGen gene-disease validity of its genes, recruiting trials (new and no longer recruiting) and papers first published since then for each open hypothesis. The first run records the baseline. Each change says what to look at again (an ACMG reading, a trial\'s eligibility); nothing is concluded. Queries carry biology only and skip anything holding a protected identifier. plan: true lists the questions without sending them.',
    inputSchema: { type: 'object', properties: { plan: { type: 'boolean' } } },
    argv: (i, casePath) => {
      if (!casePath) throw new Error('no active case (the person can run /zebra new <dir> or /zebra case <dir>)')
      return ['case', 'recheck', casePath, ...(i.plan === true ? ['--plan'] : [])]
    },
    touchesCase: true,
    deferred: true,
    timeoutMs: 600_000,
  },
  {
    name: 'access',
    description:
      'Access to a treatment, for a drug, a disease or a gene, with the source of every line: FDA and EMA status from the agencies\' own records (openFDA labels, Drugs@FDA, EMA medicine data — withdrawn and refused shown as such), approval in China from the bundled official NMPA/CDE documents (approved_in_china, named_not_approved, or not_in_bundled_list — which is not "not approved"), China\'s 2025 national reimbursement list (NRDL) entry with its restriction text verbatim, trials with sites in China (ClinicalTrials.gov; ChiCTR cannot be queried by a script — it blocks them), and the collaboration-network hospitals for a province. A name it cannot resolve exactly comes back unresolved with candidates; it never answers for a near match.',
    inputSchema: {
      type: 'object',
      properties: {
        query: { type: 'string', description: 'drug (INN or Chinese name), disease (Chinese or English) or gene symbol' },
        as: { type: 'string', enum: ['auto', 'drug', 'disease', 'gene'] },
        province: { type: 'string', description: 'list the collaboration-network hospitals of this province, e.g. 浙江' },
        trials: { type: 'number', description: 'how many trials with a site in China (default a few)' },
        status: { type: 'string', description: 'trial status filter, e.g. RECRUITING or ANY' },
      },
      required: ['query'],
    },
    argv: i => [
      'access', str(i.query) ?? '', ...flag('--as', str(i.as)), ...flag('--province', str(i.province)),
      ...flag('--trials', num(i.trials)), ...flag('--status', str(i.status)),
    ],
    timeoutMs: 180_000,
  },
  {
    name: 'expression',
    description:
      'Median expression of a gene per tissue (GTEx v8/v10) with each tissue\'s ontology id (UBERON/EFO) — for choosing the AlphaGenome tissue and the tissue an RNA test can use: the result states whether blood, lymphoblastoid cells, fibroblasts, skin or muscle reach 1 TPM. Bulk adult medians, not a detection limit; an NMD-degraded transcript reads low.',
    inputSchema: {
      type: 'object',
      properties: {
        gene: { type: 'string' },
        top: { type: 'number' },
        dataset: { type: 'string', enum: ['gtex_v8', 'gtex_v10'] },
      },
      required: ['gene'],
    },
    argv: i => ['expression', str(i.gene) ?? '', ...flag('--top', num(i.top)), ...flag('--dataset', str(i.dataset))],
    deferred: true,
  },
  {
    name: 'aso_screen',
    description:
      'Splice-switching antisense feasibility screen for researchers: from SpliceAI for the variant, the aberrant event (pseudoexon / cryptic acceptor / cryptic donor, boundaries checked for AG/GT in the patient sequence; or exon_skip to restore a reading frame), then candidate target windows on the pre-mRNA with coordinates, target and antisense sequence, GC, hairpin and homopolymer flags, whether the variant lies inside, every ranking component shown, and published N-of-1 precedents retrieved from Europe PMC. When no aberrant splicing is predicted it says so and stops. Optional genomic uniqueness via NCBI BLAST (slow: 30-700 s). A screen, not a design: no chemistry, dose or delivery; RNA evidence of the aberrant splicing comes first.',
    inputSchema: {
      type: 'object',
      properties: {
        variant: { type: 'string', description: 'chrom-pos-ref-alt, transcript HGVS or rsID' },
        assembly: ASSEMBLY,
        event: { type: 'string', enum: ['auto', 'pseudoexon', 'cryptic_acceptor', 'cryptic_donor', 'exon_skip'] },
        lengths: { type: 'string', description: 'target lengths LO-HI between 12 and 40 (default 18-25)' },
        min_delta: { type: 'number', description: 'SpliceAI delta threshold (default 0.2)' },
        top: { type: 'number', description: '1-60 (default 12)' },
        uniqueness: { type: 'boolean', description: 'one NCBI BLAST search of the shortlist; adds 30-700 s' },
        distance: { type: 'number' },
      },
      required: ['variant'],
    },
    argv: i => [
      'aso', str(i.variant) ?? '', ...flag('--assembly', str(i.assembly)), ...flag('--event', str(i.event)),
      ...flag('--lengths', str(i.lengths)), ...flag('--min-delta', num(i.min_delta)), ...flag('--top', num(i.top)),
      ...flag('--distance', num(i.distance)), ...(i.uniqueness === true ? ['--uniqueness'] : []),
    ],
    deferred: true,
    timeoutMs: 600_000,
  },
  {
    name: 'report_export',
    description:
      'Turn a Markdown report from the case (family letter, visit-preparation sheet, clinician summary) into Word (.docx) and PDF — and HTML — with Chinese typography, for a family that does not use a terminal: the person running zebra-mod hands them the files. Written next to the report. Refuses a report that still contains the case\'s protected identifiers or a resident ID number. PDF is printed by a browser on this machine (Chrome, Edge, Chromium, Brave) or LibreOffice; without one the result says so and the other formats are still made.',
    inputSchema: {
      type: 'object',
      properties: {
        report: { type: 'string', description: 'the .md report: an absolute path, or relative to the case folder (reports/family-letter-2026-10-06.md)' },
        formats: { type: 'array', items: { type: 'string', enum: ['docx', 'pdf', 'html'] }, description: 'default docx and pdf' },
        title: { type: 'string', description: 'document title (default: the first heading)' },
      },
      required: ['report'],
    },
    argv: (i, casePath) => {
      const report = str(i.report) ?? ''
      const path = report && casePath && !/^(?:\/|~)/.test(report) ? `${casePath.replace(/\/$/, '')}/${report}` : report
      const formats = list(i.formats)
      return ['report', 'export', path, ...flag('--to', formats.length ? formats.join(',') : undefined), ...flag('--title', str(i.title))]
    },
    deferred: true,
    timeoutMs: 150_000,
  },
]

/** The keys the tool's own schema declares: the event also carries `tool`, `tool_use_id` and more. */
export function schemaArgs(def: ToolDef, input: Record<string, unknown>): Record<string, unknown> {
  const props = (def.inputSchema.properties ?? {}) as Record<string, unknown>
  const out: Record<string, unknown> = {}
  for (const key of Object.keys(props)) if (input[key] !== undefined) out[key] = input[key]
  return out
}

/** No value the model supplies legitimately starts with "-": the CLI would read it as a flag. */
function flagShaped(value: unknown, path = ''): string | undefined {
  if (typeof value === 'string') {
    return value.trimStart().startsWith('-') ? `${path || 'value'}: ${JSON.stringify(value)}` : undefined
  }
  if (Array.isArray(value)) {
    for (let i = 0; i < value.length; i++) {
      const hit = flagShaped(value[i], `${path}[${i}]`)
      if (hit) return hit
    }
    return undefined
  }
  if (value && typeof value === 'object') {
    for (const [k, v] of Object.entries(value as Record<string, unknown>)) {
      const hit = flagShaped(v, path ? `${path}.${k}` : k)
      if (hit) return hit
    }
  }
  return undefined
}

export async function toolArgv(def: ToolDef, input: Record<string, unknown>, casePath: string | null): Promise<string[]> {
  const clean = schemaArgs(def, input)
  // `case_update` carries the person's own prose (notes, questions), where a leading
  // dash is harmless: it is passed as one JSON argument, never as a command-line value.
  if (def.name !== 'case_update') {
    const shaped = flagShaped(clean)
    if (shaped) throw new Error(`${shaped} is not a value this tool takes (it reads as a command-line flag)`)
  }
  const argv = def.argv(clean, casePath)
  if (argv.some(a => a === '')) throw new Error('a required argument is empty')
  return argv
}
