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
      'Record findings in the active case. phenotypes: HPO terms present or excluded (labels are verified against HPO; never invent an id — find it with hpo_search). variants, hypotheses (with ORPHA/OMIM/MONDO ids and ledger evidence ids for and against), ACMG readings of recorded variants (codes in, class computed by zebra), therapy leads, questions for the care team, and removals. Several at once.',
    inputSchema: {
      type: 'object',
      properties: {
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
        remove: {
          type: 'array',
          items: {
            type: 'object',
            properties: { kind: { type: 'string', enum: ['phenotype', 'variant', 'hypothesis', 'lead'] }, id: { type: 'string' } },
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
      'Find Human Phenotype Ontology terms for a clinical phrase (English; translate Chinese descriptions first). Returns verified ids and labels. Use it for every phenotype before recording it.',
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
      'Phenotype-driven differential diagnosis: rank diseases and genes for a set of HPO terms (present, and excluded). Sources: local (offline Resnik best-match over HPO annotations), monarch (Monarch semantic similarity), pubcasefinder (PubCaseFinder). Each source ranks on its own; agreement across them is the signal. from_case uses the active case\'s phenotypes.',
    inputSchema: {
      type: 'object',
      properties: {
        present: { type: 'array', items: HPO },
        excluded: { type: 'array', items: HPO },
        from_case: { type: 'boolean' },
        sources: { type: 'array', items: { type: 'string', enum: ['local', 'monarch', 'pubcasefinder'] } },
        top: { type: 'number', default: 15 },
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
      return [...out, ...flag('--top', num(i.top))]
    },
    timeoutMs: 180_000,
  },
  {
    name: 'gene_card',
    description:
      'Everything that matters about one gene for rare disease: HGNC identity, associated diseases with inheritance (Monarch/OMIM/Orphanet), ClinGen gene–disease validity and dosage sensitivity, gnomAD constraint (pLI, LOEUF, missense Z), PanelApp panels, protein (UniProt). Use symbols as HGNC spells them.',
    inputSchema: { type: 'object', properties: { symbol: { type: 'string' } }, required: ['symbol'] },
    argv: i => ['gene', str(i.symbol) ?? ''],
  },
  {
    name: 'variant_card',
    description:
      'Annotate one variant: normalized forms (HGVS, GRCh38/37 coordinates), consequence on the MANE transcript (Ensembl VEP), population frequency (gnomAD, by genetic ancestry group, with filtering AF), ClinVar classification with review status, in-silico predictors (REVEL, AlphaMissense, CADD, SpliceAI where available) and literature mentions (LitVar). Input: HGVS with transcript (NM_...:c.), rsID, or chrom-pos-ref-alt with assembly.',
    inputSchema: {
      type: 'object',
      properties: { variant: { type: 'string' }, assembly: ASSEMBLY, gene: { type: 'string' } },
      required: ['variant'],
    },
    argv: i => ['variant', str(i.variant) ?? '', ...flag('--assembly', str(i.assembly)), ...flag('--gene', str(i.gene))],
    timeoutMs: 180_000,
  },
  {
    name: 'disease_card',
    description:
      'One rare disease: identifiers across ORPHA, OMIM, MONDO, ICD; definition, prevalence, inheritance, age of onset, associated genes (Orphanet / Monarch), GeneReviews chapter, whether it is on China\'s national rare disease lists. Input: a name or an id (ORPHA:33069, OMIM:607208, MONDO:0100135).',
    inputSchema: { type: 'object', properties: { query: { type: 'string' } }, required: ['query'] },
    argv: i => ['disease', str(i.query) ?? ''],
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
      },
      required: ['mode'],
    },
    argv: i => {
      if (i.mode === 'suggest') {
        return ['acmg', 'suggest', str(i.variant) ?? '', ...flag('--assembly', str(i.assembly)), ...flag('--inheritance', str(i.inheritance))]
      }
      return ['acmg', 'classify', ...list(i.codes)]
    },
    timeoutMs: 180_000,
  },
  {
    name: 's2f_predict',
    description:
      'Sequence-to-function predictions for a variant: splicing (SpliceAI and Pangolin via the Broad lookup service), and when keys or the s2f CLI are present, AlphaGenome (expression/splicing/chromatin tracks) and Evo 2 (zero-shot likelihood). Returns each model separately with its claim ceiling; use for splice-region, deep intronic, UTR, promoter and other non-coding variants, and to test a mechanism.',
    inputSchema: {
      type: 'object',
      properties: {
        variant: { type: 'string', description: 'chrom-pos-ref-alt or HGVS' },
        assembly: ASSEMBLY,
        models: { type: 'array', items: { type: 'string', enum: ['spliceai', 'pangolin', 'alphagenome', 'evo2'] } },
        distance: { type: 'number', description: 'SpliceAI window around the variant (default 500)' },
      },
      required: ['variant'],
    },
    argv: i => {
      const models = list(i.models)
      return [
        's2f', 'predict', str(i.variant) ?? '',
        ...flag('--assembly', str(i.assembly)),
        ...(models.length ? ['--models', models.join(',')] : []),
        ...flag('--distance', num(i.distance)),
      ]
    },
    timeoutMs: 300_000,
  },
  {
    name: 'therapy_landscape',
    description:
      'Genotype-to-therapy landscape for a disease or gene: approved and investigational drugs with mechanism and phase (Open Targets / ChEMBL), target tractability, and pointers for N-of-1 routes (antisense, gene therapy, base editing) to check. Results are leads with their evidence, not recommendations.',
    inputSchema: { type: 'object', properties: { query: { type: 'string', description: 'disease name/id or gene symbol' } }, required: ['query'] },
    argv: i => ['therapy', str(i.query) ?? ''],
    timeoutMs: 180_000,
  },
  {
    name: 'trials_search',
    description:
      'Clinical trials from ClinicalTrials.gov (v2 API): by condition and optionally gene/intervention keyword, country and recruitment status. Returns NCT ids, phase, status, locations count and contacts link.',
    inputSchema: {
      type: 'object',
      properties: {
        condition: { type: 'string' },
        term: { type: 'string', description: 'extra keyword: gene, drug, modality' },
        country: { type: 'string', description: 'e.g. China, United States' },
        status: { type: 'string', enum: ['RECRUITING', 'NOT_YET_RECRUITING', 'ACTIVE_NOT_RECRUITING', 'COMPLETED', 'ANY'] },
        limit: { type: 'number', default: 20 },
      },
      required: ['condition'],
    },
    argv: i => [
      'trials', str(i.condition) ?? '',
      ...flag('--term', str(i.term)), ...flag('--country', str(i.country)), ...flag('--status', str(i.status)), ...flag('--limit', num(i.limit)),
    ],
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
        limit: { type: 'number', default: 15 },
      },
    },
    argv: i => [
      'lit', ...(str(i.query) ? [str(i.query) as string] : []),
      ...flag('--gene', str(i.gene)), ...flag('--variant', str(i.variant)), ...flag('--limit', num(i.limit)),
    ],
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
      const out = ['stats', str(i.method) ?? '']
      const params = (i.params && typeof i.params === 'object' ? i.params : {}) as Record<string, unknown>
      for (const [k, v] of Object.entries(params)) {
        const key = `--${k.replace(/_/g, '-')}`
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
      'Is this variant a candidate for base editing? Checks whether an adenine or cytosine base editor could revert it (A•T↔G•C or C•G↔T•A), lists SpCas9-family protospacers that put the base in the editing window with bystander bases, from the Ensembl reference sequence. A feasibility screen for researchers, not a design.',
    inputSchema: {
      type: 'object',
      properties: { variant: { type: 'string', description: 'chrom-pos-ref-alt (SNV)' }, assembly: ASSEMBLY },
      required: ['variant'],
    },
    argv: i => ['edit', str(i.variant) ?? '', ...flag('--assembly', str(i.assembly))],
    deferred: true,
  },
  {
    name: 'china_rare',
    description:
      "China's national rare disease lists (第一批罕见病目录 2018, 第二批 2023): is a disease on them, its Chinese name and list number. Input: Chinese or English name.",
    inputSchema: { type: 'object', properties: { query: { type: 'string' } }, required: ['query'] },
    argv: i => ['china', str(i.query) ?? ''],
    deferred: true,
  },
]

export async function toolArgv(def: ToolDef, input: Record<string, unknown>, casePath: string | null): Promise<string[]> {
  const argv = def.argv(input, casePath)
  if (argv.some(a => a === '')) throw new Error('a required argument is empty')
  return argv
}
