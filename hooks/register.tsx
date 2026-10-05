import { atom, read, update } from 'claude-code'
import type { EngineInterface, Register } from 'claude-code'

import type { Board, Ready } from '../types'
import { DOCTRINE, renderDoctrine } from './doctrine'
import { TOOLS, toolArgv, type ToolDef } from './tools'
import { guardInput, uploadsGenome } from './privacy'

const PLUGIN = 'zebra-mod'
const PANE = 'zebra-board'
const TOOL_PREFIX = `mcp__${PLUGIN}__`
const MAX_RESULT_CHARS = 60_000

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

    // python and the CLI
    try {
      const v = await $.process.run([python, `${$.plugin.root}/bin/zebra`, '--version'], { timeoutMs: 20_000 })
      const version = v.exitCode === 0 ? v.stdout.trim() : null
      await update($, ready, () => ({ python, version, error: version ? null : (v.stderr || 'zebra CLI failed').slice(0, 300) }))
      if (!version) $.ui.toast(`zebra-mod: ${python} could not run the zebra CLI — run /zebra doctor`)
    } catch (err) {
      await update($, ready, () => ({ python: null, version: null, error: String(err).slice(0, 300) }))
      $.ui.toast(`zebra-mod: ${python} not found; set the plugin's python option (needs Python 3.9+)`)
    }

    // bin/ on PATH so skills can call `zebra` from Bash
    const path = (await $.env.get('PATH')) ?? ''
    const bin = `${$.plugin.root}/bin`
    if (!path.split(':').includes(bin)) await $.env.set('PATH', `${bin}:${path}`)

    // active case: the working directory if it is one, else the last one used
    let active: string | null = null
    if (await $.fs.exists(`${e.cwd}/case.json`)) {
      try {
        const raw = JSON.parse(await $.fs.read(`${e.cwd}/case.json`)) as { schema?: string }
        if (raw.schema === 'zebra.case/1') active = e.cwd
      } catch {
        active = null
      }
    }
    if (!active) {
      const saved = await $.store.get('casePath')
      if (typeof saved === 'string' && (await $.fs.exists(`${saved}/case.json`))) active = saved
    }
    await setCase($, python, active)

    for (const def of TOOLS) {
      await $.tool.register({ name: def.name, description: def.description, inputSchema: def.inputSchema })
    }
    await $.command.register({
      name: 'zebra',
      description: 'zebra-mod: rare-disease workspace — /zebra [board|case <dir>|new <dir> [title]|close|doctor|ledger]',
      argumentHint: '[board|case <dir>|new <dir> [title]|close|doctor|ledger]',
    })

    // keep the board in step with edits made outside the mod's tools (zebra CLI in Bash, an editor)
    let lastSeen = ''
    $.clock.every(4000, async () => {
      const active = await read($, casePath)
      if (!active) return
      try {
        const st = await $.fs.stat(`${active}/case.json`)
        const led = await $.fs.stat(`${active}/evidence/ledger.jsonl`).catch(() => undefined)
        const stamp = `${st.mtimeMs}:${led?.mtimeMs ?? 0}:${led?.size ?? 0}`
        if (stamp !== lastSeen) {
          lastSeen = stamp
          await refreshBoard($, python)
        }
      } catch {
        // case folder moved or deleted: leave the board as it was
      }
    })

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
    const ids = privacyOn ? await read($, guard) : []
    const isOwn = e.tool.startsWith(TOOL_PREFIX)
    const outbound =
      isOwn ||
      e.tool === 'WebFetch' ||
      e.tool === 'WebSearch' ||
      (e.tool.startsWith('mcp__') && !isOwn) ||
      (e.tool === 'Bash' && /\b(curl|wget|http|https|nc|ncat|scp|sftp|rsync|ssh|ftp|gh|aws|gsutil|rclone|python3?\s+-c)\b/.test(String((e.input as { command?: unknown })?.command ?? '')))

    if (privacyOn && outbound) {
      const hit = guardInput(e.input, ids)
      if (hit) {
        return {
          decision: 'deny' as const,
          reason: `zebra-mod privacy gate: this call would send ${hit} off this machine. Query with HPO ids, gene symbols and variants only; protected identifiers are listed in the case (zebra case identifiers).`,
        }
      }
      if (e.tool === 'Bash') {
        const cmd = String((e.input as { command?: unknown })?.command ?? '')
        if (uploadsGenome(cmd)) {
          return {
            decision: 'ask' as const,
            reason: 'zebra-mod: this command may move raw genome data (VCF/BAM/CRAM/FASTQ) to another machine. Genomes identify people and their relatives — confirm the destination is one you trust.',
          }
        }
      }
    }
    if (isOwn) {
      const def = TOOLS.find(t => `${TOOL_PREFIX}${t.name}` === e.tool)
      if (def) return { decision: 'allow' as const, reason: 'zebra-mod tool: public database lookups and local case files only' }
    }
    return next(e)
  })

  // ------------------------------------------------------------ /zebra

  on('command.run', { command: 'zebra' }, async ($, e) => {
    const [verb = '', ...rest] = e.args.trim().split(/\s+/).filter(Boolean)
    const arg = rest.join(' ')
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
      if (!arg) return { text: `Active case: ${(await read($, casePath)) ?? 'none'}` }
      const dir = await absolute($, arg)
      if (!(await $.fs.exists(`${dir}/case.json`))) return { text: `No case.json in ${dir}. /zebra new ${arg} creates one.` }
      await setCase($, python, dir)
      await $.ui.open({ id: PANE, title: 'zebra · case board' })
      return { text: `Active case: ${dir}` }
    }
    if (verb === 'new') {
      const [dirArg, ...titleParts] = rest
      if (!dirArg) return { text: 'Usage: /zebra new <dir> [title]' }
      const dir = await absolute($, dirArg)
      const made = await runZebra($, python, ['case', 'init', dir, ...(titleParts.length ? ['--title', titleParts.join(' ')] : [])])
      if (!made.ok) return { text: `Could not create the case: ${made.error?.message ?? 'unknown error'}` }
      await setCase($, python, dir)
      await $.ui.open({ id: PANE, title: 'zebra · case board' })
      return {
        text: `Case created at ${dir} (case.json, records/, evidence/, reports/). Put reports and lab results in records/.`,
        context: [`A new zebra case is active at ${dir}. Next useful step: /zebra-mod:zebra-intake to turn records into an HPO profile.`],
      }
    }
    if (verb === 'close') {
      await setCase($, python, null)
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
  const env: Record<string, string> = {}
  if (active) env.ZEBRA_CASE = active
  try {
    const ran = await $.process.run([python, `${$.plugin.root}/bin/zebra`, '--json', ...args], { env, timeoutMs })
    const out = ran.stdout.trim()
    if (!out) {
      return { ok: false, error: { type: 'NoOutput', message: (ran.stderr || `exit ${ran.exitCode}`).slice(0, 2000) } }
    }
    return JSON.parse(out) as Envelope
  } catch (err) {
    return { ok: false, error: { type: 'ProcessError', message: String(err).slice(0, 2000) } }
  }
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
  // the protected identifiers themselves, for the privacy gate
  try {
    const raw = JSON.parse(await $.fs.read(`${active}/case.json`)) as { privacy?: { identifiers?: string[] } }
    const ids = (raw.privacy?.identifiers ?? []).filter(s => typeof s === 'string' && s.trim().length >= 2)
    await update($, guard, () => ids)
  } catch {
    await update($, guard, () => [])
  }
}

async function setCase($: EngineInterface, python: string, path: string | null): Promise<void> {
  await update($, casePath, () => path)
  await $.store.set('casePath', path)
  if (path) await $.env.set('ZEBRA_CASE', path)
  else await $.env.set('ZEBRA_CASE', undefined)
  await refreshBoard($, python)
}

// ------------------------------------------------------------ text

async function absolute($: { session: { cwd: () => Promise<string> } }, p: string): Promise<string> {
  if (p.startsWith('/')) return p.replace(/\/+$/, '')
  const cwd = await $.session.cwd()
  return `${cwd}/${p}`.replace(/\/\.\//g, '/').replace(/\/+$/, '')
}

function statusLine(b: Board): string {
  const present = b.phenotypes.filter(p => p.status === 'present').length
  const lead = b.hypotheses.find(h => h.status === 'confirmed' || h.status === 'leading')
  return `🦓 ${b.title} · HPO ${present} · variants ${b.variants.length}${lead ? ` · ${lead.status}: ${lead.disease}` : ''} · E${b.evidence_count}`
}

function boardMarkdown(b: Board): string {
  const lines: string[] = []
  const present = b.phenotypes.filter(p => p.status === 'present')
  const excluded = b.phenotypes.filter(p => p.status === 'excluded')
  lines.push(`**Phenotypes** (${present.length} present, ${excluded.length} excluded)`)
  for (const p of present.slice(0, 20)) lines.push(`- ${p.label ?? p.id} \`${p.id}\``)
  if (present.length > 20) lines.push(`- … ${present.length - 20} more`)
  for (const p of excluded.slice(0, 8)) lines.push(`- ~~${p.label ?? p.id}~~ \`${p.id}\` (absent)`)
  lines.push('', `**Variants** (${b.variants.length})`)
  for (const v of b.variants) lines.push(`- ${v.gene ?? ''} ${v.label ?? ''} ${v.zygosity ?? ''} — *${v.classification ?? 'unclassified'}*`)
  lines.push('', `**Hypotheses** (${b.hypotheses.length})`)
  for (const h of b.hypotheses) lines.push(`- [${h.status}] ${h.disease} (+${h.support} / −${h.against})`)
  if (b.therapy_leads.length) {
    lines.push('', `**Therapy leads** (${b.therapy_leads.length})`)
    for (const t of b.therapy_leads) lines.push(`- [${t.kind}] ${t.name}${t.status ? ` — ${t.status}` : ''}`)
  }
  if (b.questions.length) {
    lines.push('', `**Questions for the care team** (${b.questions.length})`)
    for (const q of b.questions.slice(0, 10)) lines.push(`- ${q}`)
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
    '   Skills: /zebra-mod:zebra (start here), zebra-intake, zebra-diagnose, zebra-variant,',
    '   zebra-reanalysis, zebra-s2f, zebra-therapy, zebra-stats, zebra-literature, zebra-family, zebra-report',
  ].join('\n')
}

function doctorText(result: unknown): string {
  const r = (result ?? {}) as { checks?: Array<{ name: string; ok: boolean; detail?: string }> }
  const rows = (r.checks ?? []).map(c => `${c.ok ? '✓' : '✗'} ${c.name}${c.detail ? ` — ${c.detail}` : ''}`)
  return rows.length ? rows.join('\n') : JSON.stringify(result, null, 1)
}

function formatEnvelope(def: ToolDef, got: Envelope): string {
  if (!got.ok) {
    return `zebra ${def.name} failed: ${got.error?.type ?? 'error'}: ${got.error?.message ?? 'no message'}`
  }
  const body = {
    result: got.result,
    sources: got.sources ?? [],
    warnings: got.warnings ?? [],
    ledger: got.ledger ?? [],
  }
  let text = JSON.stringify(body)
  if (text.length > MAX_RESULT_CHARS) {
    text = `${text.slice(0, MAX_RESULT_CHARS)}… [truncated ${text.length - MAX_RESULT_CHARS} chars; run the same query with the zebra CLI and --json to see all]`
  }
  return text
}
