export type BoardPhenotype = { id: string; label: string | null; status: string }

export type BoardVariant = {
  id: string
  gene: string | null
  label: string | null
  zygosity: string | null
  classification: string | null
}

export type BoardHypothesis = { id: string; disease: string; status: string; support: number; against: number }

export type BoardLead = { id: string; name: string; kind: string; status: string | null }

export type BoardTest = { id: string; type: string; date: string | null; result: string | null }

export type BoardRelative = { id: string; relation: string | null; affected: boolean | string | null; genotype: string | null }

export type Board = {
  path: string
  id: string
  title: string
  role: string
  language: string | null
  updated_at: string
  phenotypes: BoardPhenotype[]
  variants: BoardVariant[]
  hypotheses: BoardHypothesis[]
  therapy_leads: BoardLead[]
  questions: string[]
  tests?: BoardTest[]
  family?: BoardRelative[]
  evidence_count: number
  identifiers: number
}

/** A zebra tool call in flight, for the spinner, the band and the row. */
export type ZebraRun = { id: string; tool: string; queryZh: string; queryEn: string; sources: string[]; startedAt: number }

/** What zebra-mod did in the current turn: the line left at its end. */
export type ZebraTurn = {
  calls: number
  failed: number
  blocked: number
  sources: string[]
  cached: number
  evidence: number
  firstEvidence: number | null
  lastEvidence: number | null
}

export type Ready = { python: string | null; version: string | null; error: string | null }

declare module 'claude-code' {
  interface PluginState {
    'zebra-mod': {
      board: Board | null
      casePath: string | null
      guard: string[]
      ready: Ready | null
      trusted: string[]
      sessionIds: string[]
      running: ZebraRun[]
      turn: ZebraTurn
      blockedAt: number | null
      lang: string | null
      blockedWhy: string
    }
  }
}
