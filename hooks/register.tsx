import { atom, read, update } from 'claude-code'
import type { EngineInterface, Register } from 'claude-code'

import type { Board, Ready } from '../types'
import { DOCTRINE, renderDoctrine } from './doctrine'
import { TOOLS, toolArgv, type ToolDef } from './tools'
import { guardInput, isOutboundShell, shellWords, uploadedPaths, uploadsGenome } from './privacy'

const PLUGIN = 'zebra-mod'
const PANE = 'zebra-board'
const TOOL_PREFIX = `mcp__${PLUGIN}__`
const MAX_RESULT_CHARS = 60_000
const MAX_ARTIFACT_SCAN = 200_000

// Tools that never reach the network: they read and write the local case only.
const LOCAL_TOOLS = new Set(['case_status', 'case_update'])
// Tools whose answer is a lookup with no side effect, safe to run without asking.
const READ_ONLY_TOOLS = new Set([
  'case_status', 'hpo_search', 'phenotype_rank', 'gene_card', 'variant_card', 'disease_card', 'acmg',
  's2f_predict', 'therapy_landscape', 'trials_search', 'literature_search', 'rare_stats', 'edit_check', 'china_rare',
  'cnv_interpret',
])

const board = atom({ plugin: 'zebra-mod', key: 'board' } as const, null as Board | null)
const casePath = atom({ plugin: 'zebra-mod', key: 'casePath' } as const, null as string | null)
const guard = atom({ plugin: 'zebra-mod', key: 'guard' } as const, [] as string[])
const ready = atom({ plugin: 'zebra-mod', key: 'ready' } as const, null as Ready | null)

type Envelope = {
  ok: boolean
  result?: unknown
  sources?: unknown[]
  warnings?: string[]
  ledger?: string[]
  error?: { type: string; message: string }
}

export const register: Register = (on, options) => {
  const python = String(options.python ?? 'python3')
  const doctrineMode = String(options.doctrine ?? 'always')
  const privacyOn = options.privacyGate !== false

  // ------------------------------------------------------------ session

  on('session.start', async ($, e, next) => {
    const started = await next(e)

    // The tools first: a slow interpreter must not delay the first prompt.
    for (const def of TOOLS) {
      await $.tool.register({ name: def.name, description: def.description, inputSchema: def.inputSchema })
    }
    await $.command.register({
      name: 'zebra',
      description: 'zebra-mod: rare-disease workspace — /zebra [board|case <dir>|new <dir> [title]|close|doctor|ledger]',
      argumentHint: '[board|case <dir>|new <dir> [title]|close|doctor|ledger]',
    })

    // bin/ on PATH, and the interpreter the tools use, so `zebra` in Bash is the same CLI
    const path = (await $.env.get('PATH')) ?? ''
    const bin = `${$.plugin.root}/bin`
    if (!path.split(':').includes(bin)) await $.env.set('PATH', `${bin}:${path}`)
    await $.env.set('ZEBRA_PYTHON', python)

    void (async () => {
      try {
        const v = await $.process.run([python, `${$.plugin.root}/bin/zebra`, '--version'], { timeoutMs: 20_000 })
        const version = v.exitCode === 0 ? v.stdout.trim() : null
        await update($, ready, () => ({ python, version, error: version ? null : (v.stderr || 'zebra CLI failed').slice(0, 300) }))
        if (!version) $.ui.toast(`zebra-mod: ${python} could not run the zebra CLI — run /zebra doctor`)
      } catch (err) {
        await update($, ready, () => ({ python: null, version: null, error: String(err).slice(0, 300) }))
        $.ui.toast(`zebra-mod: ${python} not found; set the plugin's python option (needs Python 3.9+)`)
      }

      // A case is adopted only when this directory IS one, or lies inside the one
      // remembered for it. A case is never inherited from another project or session.
      const held = await read($, casePath)
      let active: string | null = held
      if (active === null) {
        const own = await caseAt($, e.cwd)
        if (own) active = own
        else {
          const remembered = await rememberedFor($, e.cwd)
          if (remembered && (await caseAt($, remembered))) active = remembered
          else if (remembered) {
            $.ui.toast(`zebra-mod: the case remembered here (${remembered}) is gone; /zebra case <dir> to pick one`)
          }
        }
      }
      if (active) await setCase($, python, active, e.isInteractive)
      else $.ui.status(undefined)

      // keep the board in step with edits made outside the mod's tools (the CLI in Bash, an editor)
      let lastSeen = ''
      $.clock.every(4000, async () => {
        const now = await read($, casePath)
        if (!now) return
        try {
          const st = await $.fs.stat(`${now}/case.json`)
          const led = await $.fs.stat(`${now}/evidence/ledger.jsonl`).catch(() => undefined)
          const stamp = `${st.mtimeMs}:${led?.mtimeMs ?? 0}:${led?.size ?? 0}`
          if (stamp !== lastSeen) {
            lastSeen = stamp
            await refreshBoard($, python)
          }
        } catch {
          // the case folder moved or was deleted: leave the board as it was
        }
      })
    })()

    return started
  })

  // ------------------------------------------------------------ the doctrine

  on('prompt.compose', async ($, e, next) => {
    const composed = await next(e)
    if (doctrineMode === 'off' || e.traits.includes('bare')) return composed
    const active = await read($, casePath)
    if (doctrineMode === 'case' && !active) return composed
    const b = await read($, board)
    return {
      sections: [
        ...composed.sections,
        { id: 'zebra-mod:doctrine', text: renderDoctrine(DOCTRINE, active, b), scope: 'session' as const },
      ],
    }
  })

  // ------------------------------------------------------------ tools

  on('tool.call', async ($, e, next) => {
    if (!e.tool.startsWith(TOOL_PREFIX)) return next(e)
    const def = TOOLS.find(t => `${TOOL_PREFIX}${t.name}` === e.tool)
    if (!def) return next(e)
    const input = e as unknown as Record<string, unknown>
    let argv: string[]
    try {
      argv = await toolArgv(def, input, await read($, casePath))
    } catch (err) {
      return { result: `zebra ${def.name}: ${String(err instanceof Error ? err.message : err)}` }
    }
    const got = await runZebra($, python, argv, def.timeoutMs ?? 120_000)
    if (def.touchesCase) await refreshBoard($, python)
    return { result: formatEnvelope(def, got) }
  })

  on('tool.describe', async ($, e, next) => {
    if (!e.tool.startsWith(TOOL_PREFIX)) return next(e)
    const def = TOOLS.find(t => `${TOOL_PREFIX}${t.name}` === e.tool)
    if (!def) return next(e)
    const described = await next(e)
    return { ...described, isDeferred: def.deferred === true }
  })

  // ------------------------------------------------------------ permissions and privacy

  on('tool.check', async ($, e, next) => {
    const own = e.tool.startsWith(TOOL_PREFIX) ? TOOLS.find(t => `${TOOL_PREFIX}${t.name}` === e.tool) : undefined
    const isLocalOwn = own !== undefined && LOCAL_TOOLS.has(own.name)
    const command = e.tool === 'Bash' ? String((e.input as { command?: unknown })?.command ?? '') : ''
    const isRemoteAgent = e.tool === 'Agent' && (e.input as { isolation?: unknown })?.isolation === 'remote'
    const outbound =
      (own !== undefined && !isLocalOwn) ||
      e.tool === 'WebFetch' ||
      e.tool === 'WebSearch' ||
      e.tool === 'Artifact' ||
      e.tool === 'SendMessage' ||
      isRemoteAgent ||
      (e.tool.startsWith('mcp__') && own === undefined) ||
      (e.tool === 'Bash' && isOutboundShell(command))

    if (privacyOn && outbound) {
      const active = await read($, casePath)
      const held = await identifiersNow($, active)
      if (held.isClosed) {
        return {
          decision: 'deny' as const,
          reason: `zebra-mod privacy gate: ${active}/case.json cannot be read, so the protected identifiers are unknown and the gate is closed. Fix or re-create the case file, or /zebra close to work without one.`,
        }
      }
      const hit = guardInput(e.input, held.ids)
      if (hit) {
        return {
          decision: 'deny' as const,
          reason: `zebra-mod privacy gate: this call would send ${hit} off this machine. Query public databases with HPO ids, gene symbols, variants and disease names only.`,
        }
      }
      // a file from the case folder about to be sent somewhere
      const paths = e.tool === 'Bash' ? uploadedPaths(command) : artifactPaths(e.input)
      if (paths.length > 0 && active) {
        const inside = await firstInsideCase($, active, paths)
        if (inside) {
          return {
            decision: 'ask' as const,
            reason: `zebra-mod: this would send ${inside}, a file inside the case folder, off this machine. Its contents are not checked for names or record numbers — confirm only if you have read it and trust the destination.`,
          }
        }
      }
      if (e.tool === 'Artifact') {
        const leak = await artifactLeak($, e.input, held.ids)
        if (leak) {
          return {
            decision: 'deny' as const,
            reason: `zebra-mod privacy gate: the page about to be published contains ${leak}. Publishing puts it on the web.`,
          }
        }
      }
      if (e.tool === 'Bash' && uploadsGenome(command)) {
        return {
          decision: 'ask' as const,
          reason: 'zebra-mod: this command may move raw genome data (VCF/BAM/CRAM/FASTQ) to another machine. A genome identifies a person and their relatives — confirm the destination is one you trust.',
        }
      }
    }

    // The engine's own decision stands; this only spares the person a prompt for a
    // read-only lookup, and never overrides a deny, a rule or the mode's own answer.
    const below = await next(e)
    if (own && below.decision === 'ask' && below.rule === undefined && READ_ONLY_TOOLS.has(own.name)) {
      return { decision: 'allow' as const, reason: 'zebra-mod: read-only lookup in public databases and local case files' }
    }
    return below
  })

  // ------------------------------------------------------------ /zebra

  on('command.run', { command: 'zebra' }, async ($, e) => {
    const words = shellWords(e.args.trim())
    const verb = words[0] ?? ''
    const rest = words.slice(1)
    if (verb === '' || verb === 'help') {
      const r = await read($, ready)
      const active = await read($, casePath)
      if (active) await $.ui.open({ id: PANE, title: 'zebra · case board' })
      return { text: helpText(r, active) }
    }
    if (verb === 'board') {
      const active = await read($, casePath)
      if (!active) return { text: 'No active case. /zebra new <dir> [title] creates one; /zebra case <dir> opens one.' }
      await refreshBoard($, python)
      const opened = await $.ui.open({ id: PANE, title: 'zebra · case board' })
      return { text: opened.isPlaced ? 'Case board opened.' : 'Case board waits for a wider terminal.' }
    }
    if (verb === 'case') {
      const arg = rest.join(' ')
      if (!arg) return { text: `Active case: ${(await read($, casePath)) ?? 'none'}` }
      const dir = await absolute($, rest.length === 1 ? (rest[0] as string) : arg)
      if (!(await $.fs.exists(`${dir}/case.json`))) return { text: `No case.json in ${dir}. /zebra new ${arg} creates one.` }
      await setCase($, python, dir, e.origin.kind === 'composer')
      await $.ui.open({ id: PANE, title: 'zebra · case board' })
      return { text: `Active case: ${dir}` }
    }
    if (verb === 'new') {
      const dirArg = rest[0]
      const title = rest.slice(1).join(' ')
      if (!dirArg) return { text: 'Usage: /zebra new <dir> [title]   (quote a path or title that contains spaces)' }
      const dir = await absolute($, dirArg)
      const made = await runZebra($, python, ['case', 'init', dir, ...(title ? ['--title', title] : [])])
      if (!made.ok) return { text: `Could not create the case: ${made.error?.message ?? 'unknown error'}` }
      await setCase($, python, dir, e.origin.kind === 'composer')
      await $.ui.open({ id: PANE, title: 'zebra · case board' })
      return {
        text: `Case created at ${dir} (case.json, records/, evidence/, reports/). Put reports and lab results in records/.`,
        context: [`A new zebra case is active at ${dir}. Next useful step: /zebra-mod:zebra-intake to turn records into an HPO profile.`],
      }
    }
    if (verb === 'close') {
      await setCase($, python, null, e.origin.kind === 'composer')
      await $.ui.close({ id: PANE })
      return { text: 'No active case.' }
    }
    if (verb === 'doctor') {
      const got = await runZebra($, python, ['doctor'], 90_000)
      return { text: got.ok ? doctorText(got.result) : `zebra doctor failed: ${got.error?.message}` }
    }
    if (verb === 'ledger') {
      const active = await read($, casePath)
      if (!active) return { text: 'No active case.' }
      const got = await runZebra($, python, ['case', 'ledger', '--case', active])
      if (!got.ok) return { text: `Could not read the ledger: ${got.error?.message ?? 'unknown error'}` }
      const rows = (got.result as Array<Record<string, unknown>> | undefined) ?? []
      const tail = rows.slice(-25).map(r => `${r.eid}  ${r.db}  ${r.record ?? ''}  ${r.url ?? ''}`)
      return { text: rows.length ? `${rows.length} evidence rows (last 25):\n${tail.join('\n')}` : 'The evidence ledger is empty.' }
    }
    return { text: `Unknown: /zebra ${verb}. Try /zebra help.` }
  })

  // ------------------------------------------------------------ the case board

  on('ui.render', { component: 'Pane', requestId: PANE }, async ($, e) => {
    const { Box, Text, Button, Markdown } = $.ui.resolve(e)
    const b = await read($, board)
    if (!b) {
      return (
        <Box flexDirection="column">
          <Text dimColor>No active case. /zebra new &lt;dir&gt; [title]</Text>
        </Box>
      )
    }
    return (
      <Box flexDirection="column">
        <Text bold>🦓 {b.title}</Text>
        <Text dimColor>
          {b.role} · {b.evidence_count} evidence rows · {b.identifiers} protected identifiers
        </Text>
        <Markdown key="body" text={boardMarkdown(b)} />
        <Box flexDirection="row" gap={1}>
          <Button key="refresh" label="Refresh" hotkey="r" onPress={() => refreshBoard($, python)} />
          <Button key="close" label="Close" role="dismiss" onPress={() => $.ui.close({ id: PANE })} />
        </Box>
      </Box>
    )
  })
}

// ------------------------------------------------------------ $-using helpers (top level, as the engine requires)

async function runZebra($: EngineInterface, python: string, args: readonly string[], timeoutMs = 120_000): Promise<Envelope> {
  const active = await read($, casePath)
  const env: Record<string, string> = { ZEBRA_DEADLINE_MS: String(Math.max(5_000, timeoutMs - 10_000)) }
  if (active) env.ZEBRA_CASE = active
  try {
    const ran = await $.process.run([python, `${$.plugin.root}/bin/zebra`, '--json', ...args], { env, timeoutMs })
    if (ran.isStdoutTruncated) {
      return { ok: false, error: { type: 'OutputTooLarge', message: 'the CLI wrote more than 4 MiB; narrow the query (fewer items, a shorter range)' } }
    }
    const out = ran.stdout.trim()
    if (!out) {
      // a traceback ends with the exception, so keep the tail, not the head
      const err = (ran.stderr || `exit ${ran.exitCode}`).trim()
      return { ok: false, error: { type: 'NoOutput', message: err.slice(-2000) } }
    }
    return JSON.parse(out) as Envelope
  } catch (err) {
    return { ok: false, error: { type: 'ProcessError', message: String(err).slice(0, 2000) } }
  }
}

/** The case directory at `dir`, or null when it holds no zebra case. */
async function caseAt($: EngineInterface, dir: string): Promise<string | null> {
  try {
    if (!(await $.fs.exists(`${dir}/case.json`))) return null
    const raw = JSON.parse(await $.fs.read(`${dir}/case.json`)) as { schema?: string }
    return raw.schema === 'zebra.case/1' ? dir : null
  } catch {
    return null
  }
}

/** The case last used in this working directory, if any. Cases are never shared between projects. */
async function rememberedFor($: EngineInterface, cwd: string): Promise<string | null> {
  const saved = await $.store.get('casesByDir')
  if (saved && typeof saved === 'object') {
    const held = (saved as Record<string, unknown>)[cwd]
    if (typeof held === 'string') return held
  }
  return null
}

async function refreshBoard($: EngineInterface, python: string): Promise<void> {
  const active = await read($, casePath)
  if (!active) {
    await update($, board, () => null)
    await update($, guard, () => [])
    $.ui.status(undefined)
    return
  }
  const got = await runZebra($, python, ['case', 'summary', active], 30_000)
  if (got.ok && got.result) {
    const b = got.result as Board
    await update($, board, () => b)
    $.ui.status(statusLine(b))
  } else {
    $.ui.status(`zebra: case ${active} unreadable`)
  }
  const held = await identifiersNow($, active)
  if (!held.isClosed) await update($, guard, () => held.ids)
}

/**
 * The case's protected identifiers, read now rather than from the last poll, so a
 * `zebra case identifiers --add` in the same turn is already in force.
 * `isClosed` means the file could not be read: the caller must refuse outbound calls.
 */
async function identifiersNow($: EngineInterface, active: string | null): Promise<{ ids: string[]; isClosed: boolean }> {
  if (!active) return { ids: [], isClosed: false }
  try {
    const raw = JSON.parse(await $.fs.read(`${active}/case.json`)) as { privacy?: { identifiers?: unknown } }
    const listed = raw.privacy?.identifiers
    if (listed === undefined || listed === null) return { ids: [], isClosed: false }
    if (!Array.isArray(listed)) return { ids: [], isClosed: true }
    return { ids: listed.filter((s): s is string => typeof s === 'string' && s.trim().length >= 2), isClosed: false }
  } catch {
    const kept = await read($, guard)
    return { ids: kept, isClosed: true }
  }
}

async function setCase($: EngineInterface, python: string, path: string | null, isInteractive: boolean): Promise<void> {
  await update($, casePath, () => path)
  if (isInteractive) {
    const cwd = await $.session.cwd()
    const saved = await $.store.get('casesByDir')
    const map: Record<string, unknown> = saved && typeof saved === 'object' ? { ...(saved as Record<string, unknown>) } : {}
    if (path) map[cwd] = path
    else delete map[cwd]
    await $.store.set('casesByDir', map)
  }
  if (path) await $.env.set('ZEBRA_CASE', path)
  else await $.env.set('ZEBRA_CASE', undefined)
  await refreshBoard($, python)
}

/** The first of `paths` that lies inside the case folder, resolved through links. */
async function firstInsideCase($: EngineInterface, active: string, paths: readonly string[]): Promise<string | undefined> {
  const root = await $.fs.stat(active, { resolve: true }).catch(() => undefined)
  const realRoot = root?.realPath
  if (realRoot === undefined) return undefined
  for (const p of paths) {
    const stat = await $.fs.stat(p, { resolve: true }).catch(() => undefined)
    const real = stat?.realPath
    if (real !== undefined && (real === realRoot || real.startsWith(`${realRoot}/`))) return p
  }
  return undefined
}

function artifactPaths(input: unknown): string[] {
  const i = (input ?? {}) as { file_path?: unknown; file_paths?: unknown; files?: unknown }
  const out: string[] = []
  if (typeof i.file_path === 'string') out.push(i.file_path)
  if (Array.isArray(i.file_paths)) for (const p of i.file_paths) if (typeof p === 'string') out.push(p)
  if (i.files && typeof i.files === 'object') {
    for (const v of Object.values(i.files as Record<string, unknown>)) {
      if (typeof v === 'string') out.push(v)
      else if (v && typeof v === 'object' && typeof (v as { from?: unknown }).from === 'string') out.push((v as { from: string }).from)
    }
  }
  return out
}

/** What a page about to be published carries of the case's identifiers. */
async function artifactLeak($: EngineInterface, input: unknown, ids: readonly string[]): Promise<string | undefined> {
  if (ids.length === 0) return undefined
  for (const p of artifactPaths(input)) {
    const text = await $.fs.read(p).catch(() => undefined)
    if (text === undefined) continue
    const hit = guardInput(text.slice(0, MAX_ARTIFACT_SCAN), ids)
    if (hit) return `${hit} (in ${p})`
  }
  return undefined
}

// ------------------------------------------------------------ text

async function absolute($: EngineInterface, p: string): Promise<string> {
  const home = (await $.env.get('HOME')) ?? ''
  const expanded = p === '~' ? home : p.startsWith('~/') ? `${home}/${p.slice(2)}` : p
  const raw = expanded.startsWith('/') ? expanded : `${await $.session.cwd()}/${expanded}`
  const parts: string[] = []
  for (const part of raw.split('/')) {
    if (part === '' || part === '.') continue
    if (part === '..') parts.pop()
    else parts.push(part)
  }
  return `/${parts.join('/')}`
}

function statusLine(b: Board): string {
  const present = b.phenotypes.filter(p => p.status === 'present').length
  const lead = b.hypotheses.find(h => h.status === 'confirmed' || h.status === 'leading')
  return `🦓 ${b.title} · HPO ${present} · variants ${b.variants.length}${lead ? ` · ${lead.status}: ${lead.disease}` : ''} · E${b.evidence_count}`
}

/** Case text is the person's own; it must not be read as Markdown when the board draws it. */
function plain(text: string | null | undefined, fallback = ''): string {
  const t = (text ?? fallback).replace(/[\r\n]+/g, ' ')
  return t.replace(/[\\`*_[\]<>|#~]/g, m => `\\${m}`)
}

function boardMarkdown(b: Board): string {
  const lines: string[] = []
  const present = b.phenotypes.filter(p => p.status === 'present')
  const excluded = b.phenotypes.filter(p => p.status === 'excluded')
  lines.push(`**Phenotypes** (${present.length} present, ${excluded.length} excluded)`)
  for (const p of present.slice(0, 20)) lines.push(`- ${plain(p.label, p.id)} (${plain(p.id)})`)
  if (present.length > 20) lines.push(`- … ${present.length - 20} more`)
  for (const p of excluded.slice(0, 8)) lines.push(`- ${plain(p.label, p.id)} (${plain(p.id)}, absent)`)
  lines.push('', `**Variants** (${b.variants.length})`)
  for (const v of b.variants) {
    lines.push(`- ${plain(v.gene)} ${plain(v.label)} ${plain(v.zygosity)} — ${plain(v.classification, 'unclassified')}`)
  }
  lines.push('', `**Hypotheses** (${b.hypotheses.length})`)
  for (const h of b.hypotheses) lines.push(`- [${plain(h.status)}] ${plain(h.disease)} (+${h.support} / −${h.against})`)
  if (b.therapy_leads.length) {
    lines.push('', `**Therapy leads** (${b.therapy_leads.length})`)
    for (const t of b.therapy_leads) lines.push(`- [${plain(t.kind)}] ${plain(t.name)}${t.status ? ` — ${plain(t.status)}` : ''}`)
  }
  if (b.questions.length) {
    lines.push('', `**Questions for the care team** (${b.questions.length})`)
    for (const q of b.questions.slice(0, 10)) lines.push(`- ${plain(q)}`)
  }
  return lines.join('\n')
}

function helpText(r: Ready | null, active: string | null): string {
  return [
    '🦓 zebra-mod — rare-disease research workstation',
    `   CLI: ${r?.version ?? `not ready (${r?.error ?? 'python not checked'})`}`,
    `   active case: ${active ?? 'none'}`,
    '',
    '   /zebra new <dir> [title]   start a case folder (records stay on this machine)',
    '   /zebra case <dir>          switch to an existing case',
    '   /zebra board               open the case board',
    '   /zebra ledger              list the evidence the answers stand on',
    '   /zebra doctor              check Python, data files, API reachability and keys',
    '   /zebra close               no active case',
    '',
    '   Skills: /zebra-mod:zebra-start (start here), zebra-safety (urgent red flags, drug and',
    '   anaesthesia hazards), zebra-intake, zebra-diagnose, zebra-variant, zebra-reanalysis,',
    '   zebra-s2f, zebra-therapy, zebra-stats, zebra-literature, zebra-family, zebra-report',
  ].join('\n')
}

function doctorText(result: unknown): string {
  const r = (result ?? {}) as { checks?: Array<{ name: string; ok: boolean; detail?: string }> }
  const rows = (r.checks ?? []).map(c => `${c.ok ? '✓' : '✗'} ${c.name}${c.detail ? ` — ${c.detail}` : ''}`)
  return rows.length ? rows.join('\n') : JSON.stringify(result, null, 1)
}

/** The longest list inside `result`, so a too-large answer loses items rather than its provenance. */
function trimLongestList(result: unknown): { result: unknown; dropped: number; key: string } | undefined {
  if (!result || typeof result !== 'object') return undefined
  let best: { holder: Record<string, unknown>; key: string; list: unknown[] } | undefined
  const walk = (value: unknown, depth: number): void => {
    if (depth > 4 || !value || typeof value !== 'object') return
    for (const [key, v] of Object.entries(value as Record<string, unknown>)) {
      if (Array.isArray(v) && v.length > 1 && (best === undefined || v.length > best.list.length)) {
        best = { holder: value as Record<string, unknown>, key, list: v }
      }
      walk(v, depth + 1)
    }
  }
  walk(result, 0)
  if (!best) return undefined
  const keep = Math.max(1, Math.floor(best.list.length / 2))
  const dropped = best.list.length - keep
  best.holder[best.key] = best.list.slice(0, keep)
  return { result, dropped, key: best.key }
}

function formatEnvelope(def: ToolDef, got: Envelope): string {
  if (!got.ok) {
    return `zebra ${def.name} failed: ${got.error?.type ?? 'error'}: ${got.error?.message ?? 'no message'}`
  }
  const warnings = [...(got.warnings ?? [])]
  // warnings, ledger and sources first: they are what the answer must cite, so a
  // trim never costs them. The result is trimmed until the whole envelope fits.
  let body = { warnings, ledger: got.ledger ?? [], sources: got.sources ?? [], result: got.result }
  let text = JSON.stringify(body)
  for (let i = 0; i < 12 && text.length > MAX_RESULT_CHARS; i++) {
    const trimmed = trimLongestList(body.result)
    if (!trimmed) break
    warnings.push(`zebra-mod trimmed ${trimmed.dropped} of "${trimmed.key}" to fit the tool result; run the same query with the zebra CLI and --json for all of it`)
    body = { ...body, warnings, result: trimmed.result }
    text = JSON.stringify(body)
  }
  if (text.length > MAX_RESULT_CHARS) {
    body = {
      warnings: [...warnings, 'zebra-mod could not fit this result; only its provenance is shown. Run the same query with the zebra CLI and --json.'],
      ledger: body.ledger,
      sources: body.sources,
      result: null,
    }
    text = JSON.stringify(body)
  }
  return text
}
