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
  evidence_count: number
  identifiers: number
}

export type Ready = { python: string | null; version: string | null; error: string | null }

declare module 'claude-code' {
  interface PluginState {
    'zebra-mod': {
      board: Board | null
      casePath: string | null
      guard: string[]
      ready: Ready | null
    }
  }
}
