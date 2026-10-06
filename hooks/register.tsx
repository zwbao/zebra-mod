import { atom, read, update } from 'claude-code'
import type { EngineInterface, Register } from 'claude-code'

import type { Board, Ready } from '../types'
import { DOCTRINE, DOCTRINE_BRIEF, renderDoctrine } from './doctrine'
import { schemaArgs, TOOLS, toolArgv, type ToolDef } from './tools'
import { guardInput, isOutboundShell, outboundText, sendsVariantList, shellWords, uploadedPaths, uploadsGenome } from './privacy'

const PLUGIN = 'zebra-mod'
const PANE = 'zebra-board'
const TOOL_PREFIX = `mcp__${PLUGIN}__`
const MAX_RESULT_CHARS = 60_000
const MAX_ARTIFACT_SCAN = 200_000

// Tools that never reach the network: they read and write the local case only.
const LOCAL_TOOLS = new Set(['case_status', 'case_update', 'report_export'])
// Tools whose answer is a lookup with no side effect, safe to run without asking.
// tools that can carry text off the machine: what a crashed gate refuses instead of passing
const OUTBOUND_HINT = /^(?:Bash|WebFetch|WebSearch|Artifact|SendMessage|Agent|mcp__)/
const READ_ONLY_TOOLS = new Set([
  'case_status', 'hpo_search', 'phenotype_rank', 'gene_card', 'variant_card', 'disease_card', 'acmg',
  's2f_predict', 'therapy_landscape', 'trials_search', 'literature_search', 'rare_stats', 'edit_check', 'china_rare',
  'cnv_interpret', 'access', 'expression', 'aso_screen',
])

const board = atom({ plugin: 'zebra-mod', key: 'board' } as const, null as Board | null)
const casePath = atom({ plugin: 'zebra-mod', key: 'casePath' } as const, null as string | null)
const guard = atom({ plugin: 'zebra-mod', key: 'guard' } as const, [] as string[])
const ready = atom({ plugin: 'zebra-mod', key: 'ready' } as const, null as Ready | null)
// zebra tools the person allowed for the rest of this session in zebra's own approval dialog
const trusted = atom({ plugin: 'zebra-mod', key: 'trusted' } as const, [] as string[])
// identifiers registered while no case was open: protected for this session, kept only in memory,
// and added to the next case that is opened or created
const sessionIds = atom({ plugin: 'zebra-mod', key: 'sessionIds' } as const, [] as string[])

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
  const doctrineMode = String(options.doctrine ?? 'auto')
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
      const r = await checkCli($, python)
      if (!r.version) {
        $.ui.toast(r.python
          ? `zebra-mod: ${python} could not run the zebra CLI — run /zebra doctor`
          : `zebra-mod: ${python} not found; set the plugin's python option (needs Python 3.9+)`)
      }

      // A case is adopted only when this directory is one or lies inside one, or when it was
      // the case last used in this directory. A case is never inherited from another project.
      const held = await read($, casePath)
      let active: string | null = held
      if (active === null) {
        const own = (await caseAt($, e.cwd)) ?? (await enclosingCase($, e.cwd))
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
    // auto (the default): the full doctrine with a case open; otherwise a short section that
    // applies only to rare-disease questions, so other work in the same Claude Code is untouched
    const text = doctrineMode === 'auto' && !active ? DOCTRINE_BRIEF : renderDoctrine(DOCTRINE, active, await read($, board))
    return {
      sections: [...composed.sections, { id: 'zebra-mod:doctrine', text, scope: 'session' as const }],
    }
  })

  // ------------------------------------------------------------ tools

  on('tool.call', async ($, e, next) => {
    if (!e.tool.startsWith(TOOL_PREFIX)) return next(e)
    const def = TOOLS.find(t => `${TOOL_PREFIX}${t.name}` === e.tool)
    if (!def) return next(e)
    const input = e as unknown as Record<string, unknown>
    // A tool this plugin registers is answered here, and nothing beneath this hook runs
    // the engine's permission chain for it. So the decision is asked for explicitly:
    // the person's rules, the session's mode and the privacy gate apply to zebra's own
    // tools exactly as they do to any other tool.
    const verdict = await $.tool.check({ tool: e.tool, input: schemaArgs(def, input) })
    if (verdict.decision === 'deny') return { deny: verdict.reason ?? `${def.name} was refused` }
    if (verdict.decision === 'ask' && !(await approved($, def, verdict))) {
      return { deny: `${def.name} was not approved` }
    }
    if (def.name === 'case_update' && !(await read($, casePath))) {
      const held = await holdForSession($, input)
      if (held) return { result: held }
    }
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
      const now = await identifiersNow($, active)
      // a copy: the session list comes from engine state, which must not be mutated in place
      const held = { isClosed: now.isClosed, ids: [...now.ids, ...(e.tool === 'Bash' ? await namedCaseIdentifiers($, command, active) : [])] }
      if (held.isClosed) {
        return {
          decision: 'deny' as const,
          reason: `zebra-mod privacy gate: ${active}/case.json cannot be read, so the protected identifiers are unknown and the gate is closed. Fix or re-create the case file, or /zebra close to work without one.`,
        }
      }
      // Scan what actually leaves: a shell command's outbound segments without local paths
      // or sample names; a published page's text below; anything else as given.
      const scanned = e.tool === 'Bash' ? { command: outboundText(command) }
        : e.tool === 'Artifact' ? artifactFields(e.input)
        : own !== undefined ? schemaArgs(own, (e.input ?? {}) as Record<string, unknown>)
        : e.input
      const hit = guardInput(scanned, held.ids)
      if (hit && !hit.startsWith('protected identifier') && !active && own === undefined) {
        // No case is open, so this is not known to be patient work: an email or a phone number
        // in an ordinary call (a PR body, an API request) is asked about, never refused outright.
        return {
          decision: 'ask' as const,
          reason: `zebra-mod privacy gate: this call would send ${hit} off this machine. If it belongs to a patient, do not send it; if it is yours or public, confirm.`,
        }
      }
      if (hit) {
        return {
          decision: 'deny' as const,
          reason: `zebra-mod privacy gate: this call would send ${hit} off this machine. Query public databases with HPO ids, gene symbols, variants and disease names only.`,
        }
      }
      // a file from the case folder about to be sent somewhere
      const paths = e.tool === 'Bash' ? uploadedPaths(command) : artifactPaths(e.input)
      if (paths.length > 0 && active) {
        const inside = await firstInsideCase($, active, paths, (e.input as { root?: unknown })?.root)
        if (inside) {
          return {
            decision: 'ask' as const,
            reason: `zebra-mod: this would send ${inside}, a file inside the case folder, off this machine. Its contents are not checked for names or record numbers — confirm only if you have read it and trust the destination.`,
          }
        }
      }
      if (e.tool === 'Artifact') {
        const scan = await artifactLeak($, e.input, held.ids)
        if (scan.leak) {
          return {
            decision: 'deny' as const,
            reason: `zebra-mod privacy gate: the page about to be published contains ${scan.leak}. Publishing puts it on the web.`,
          }
        }
        if (scan.unread.length > 0) {
          return {
            decision: 'ask' as const,
            reason: `zebra-mod: ${scan.unread.join(', ')} could not be read to check for names or record numbers before publishing — confirm only if you have checked it yourself.`,
          }
        }
      }
      if (e.tool === 'Bash' && uploadsGenome(command)) {
        return {
          decision: 'ask' as const,
          reason: 'zebra-mod: this command may move raw genome data (VCF/BAM/CRAM/FASTQ) to another machine. A genome identifies a person and their relatives — confirm the destination is one you trust.',
        }
      }
      if (e.tool === 'Bash' && sendsVariantList(command)) {
        return {
          decision: 'ask' as const,
          reason: 'zebra-mod: --prefilter myvariant sends every quality-passing variant position in this VCF (often tens of thousands) to myvariant.info. Taken together they are the person\'s genome — confirm, or run triage without the prefilter (only the candidates that survive local filtering go out).',
        }
      }
    }

    // The engine's own decision stands; this only spares the person a prompt for a
    // read-only lookup, and never overrides a deny, a rule or the mode's own answer.
    const below = await next(e)
    const writes = own?.name === 'cnv_interpret' && (e.input as { record?: unknown })?.record === true
    if (own && !writes && below.decision === 'ask' && below.rule === undefined && READ_ONLY_TOOLS.has(own.name)) {
      return { decision: 'allow' as const, reason: 'zebra-mod: read-only lookup in public databases and local case files' }
    }
    return below
  }).catch(async ($, e, next) => {
    if (OUTBOUND_HINT.test(e.tool)) {
      return { decision: 'deny' as const, reason: `zebra-mod privacy gate could not check this call (${next.error.kind}); it was refused rather than let through unchecked` }
    }
    return next(e)
  })

  // ------------------------------------------------------------ /zebra

  on('command.run', { command: 'zebra' }, async ($, e) => {
    const words = shellWords(e.args.trim())
    const verb = words[0] ?? ''
    const rest = words.slice(1)
    if (verb === '' || verb === 'help') {
      if ((await read($, ready)) === null) await checkCli($, python)
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
      // --tail: the CLI keeps the newest rows and reports the true total (an output trim would keep the oldest)
      const got = await runZebra($, python, ['case', 'ledger', '--case', active, '--tail', '25'])
      if (!got.ok) return { text: `Could not read the ledger: ${got.error?.message ?? 'unknown error'}` }
      const r = (got.result ?? {}) as { total?: number; rows?: Array<Record<string, unknown>> }
      const rows = r.rows ?? []
      const tail = rows.map(row => `${row.eid}  ${row.db}  ${row.record ?? ''}  ${row.url ?? ''}`)
      return { text: rows.length ? `${r.total ?? rows.length} evidence rows (newest ${rows.length}):\n${tail.join('\n')}` : 'The evidence ledger is empty.' }
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

/** Run the CLI's version check now and record the answer. */
async function checkCli($: EngineInterface, python: string): Promise<Ready> {
  let r: Ready
  try {
    const v = await $.process.run([python, `${$.plugin.root}/bin/zebra`, '--version'], { timeoutMs: 20_000 })
    const version = v.exitCode === 0 ? v.stdout.trim() : null
    r = { python, version, error: version ? null : (v.stderr || 'zebra CLI failed').slice(0, 300) }
  } catch (err) {
    r = { python: null, version: null, error: String(err).slice(0, 300) }
  }
  await update($, ready, () => r)
  return r
}

/**
 * Ask the person whether a zebra tool may run; no one to ask means no. The engine's permission
 * dialog only opens for tools it runs itself, so this is an AskUserQuestion. An ask that comes
 * from the session's mode may be answered once for the session; one a settings rule asks for
 * (`rule` set) is asked every time, as the person configured.
 */
async function approved($: EngineInterface, def: ToolDef, verdict: { reason?: string; rule?: unknown }): Promise<boolean> {
  const sessionable = verdict.rule === undefined
  if (sessionable && (await read($, trusted)).includes(def.name)) return true
  const always = `Allow ${def.name} for this session`
  const options = sessionable ? ['Allow', always, 'Deny'] : ['Allow', 'Deny']
  try {
    const why = verdict.reason ? `${verdict.reason.replace(/[?？.。]\s*$/, '')}. ` : ''
    const answer = await $.ui.ask(`${why}Allow zebra-mod to run ${def.name}?`, options)
    if (sessionable && answer === always) {
      await update($, trusted, t => (t.includes(def.name) ? t : [...t, def.name]))
      return true
    }
    return answer === 'Allow'
  } catch {
    return false // dismissed, or a headless run with nobody to ask
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

/** The case folder this directory lies inside (a session started in <case>/records), up to four levels up. */
async function enclosingCase($: EngineInterface, cwd: string): Promise<string | null> {
  let dir = cwd.replace(/\/+$/, '')
  for (let i = 0; i < 4; i++) {
    const parent = dir.slice(0, dir.lastIndexOf('/'))
    if (!parent || parent === dir) return null
    const found = await caseAt($, parent)
    if (found) return found
    dir = parent
  }
  return null
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
  const session = await read($, sessionIds)
  if (!active) return { ids: session, isClosed: false }
  try {
    const raw = JSON.parse(await $.fs.read(`${active}/case.json`)) as { privacy?: { identifiers?: unknown } }
    const listed = raw.privacy?.identifiers
    if (listed === undefined || listed === null) return { ids: session, isClosed: false }
    if (!Array.isArray(listed)) return { ids: session, isClosed: true }
    return { ids: [...session, ...listed.filter((s): s is string => typeof s === 'string' && s.trim().length >= 2)], isClosed: false }
  } catch {
    const kept = await read($, guard)
    return { ids: [...session, ...kept], isClosed: true }
  }
}

/**
 * case_update with no case open. Identifiers are the one thing worth keeping without a case: a
 * parent pastes a clinic note before any folder exists, and the gate must know the child's name
 * from that moment. They are held in memory for this session (written nowhere) and added to the
 * next case opened or created. Anything else in the call needs a case and is reported as not
 * recorded. Returns the tool's answer, or undefined when the call carries no identifiers.
 */
async function holdForSession($: EngineInterface, input: Record<string, unknown>): Promise<string | undefined> {
  const raw = Array.isArray(input.identifiers) ? input.identifiers : []
  const values = raw.filter((v): v is string => typeof v === 'string' && v.trim().length >= 2).map(v => v.trim())
  if (values.length === 0) return undefined
  const held = await update($, sessionIds, list => [...new Set([...list, ...values])])
  const other = Object.keys(input).filter(k => k !== 'identifiers' && k !== 'tool' && k !== 'tool_use_id' && input[k] !== undefined)
  const fields = TOOLS.find(t => t.name === 'case_update')
  const declared = new Set(Object.keys(((fields?.inputSchema as { properties?: object })?.properties ?? {}) as object))
  const lost = other.filter(k => declared.has(k))
  return JSON.stringify({
    warnings: lost.length ? [`no case is open, so ${lost.join(', ')} ${lost.length === 1 ? 'was' : 'were'} not recorded: /zebra new <dir> starts a case`] : [],
    result: {
      identifiers: { protected_for_this_session: held.length },
      note: 'No case is open: these identifiers are protected for the rest of this session (kept in memory only, never written) and will be added to the next case opened or created.',
    },
  })
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
  if (path) {
    const pending = await read($, sessionIds)
    if (pending.length) {
      const got = await runZebra($, python, ['case', 'apply', path, '--ops', JSON.stringify({ identifiers: pending })], 30_000)
      if (!got.ok) $.ui.toast(`zebra-mod: the identifiers held for this session could not be added to ${path}: ${got.error?.message ?? 'unknown error'}`)
    }
  }
  await refreshBoard($, python)
}

/** The first of `paths` that lies inside the case folder, resolved through links. */
async function firstInsideCase($: EngineInterface, active: string, paths: readonly string[], root?: unknown): Promise<string | undefined> {
  const caseRoot = await $.fs.stat(active, { resolve: true }).catch(() => undefined)
  const realRoot = caseRoot?.realPath
  if (realRoot === undefined) return undefined
  const home = (await $.env.get('HOME')) ?? ''
  const base = typeof root === 'string' && root ? root.replace(/\/$/, '') : ''
  for (const p of paths) {
    let q = p.replace(/^\$HOME(?=\/|$)/, home).replace(/^~(?=\/|$)/, home)
    if (base && !q.startsWith('/')) q = `${base}/${q}`
    const stat = await $.fs.stat(q, { resolve: true }).catch(() => undefined)
    const real = stat?.realPath
    if (real !== undefined && (real === realRoot || real.startsWith(`${realRoot}/`))) return p
  }
  return undefined
}

function artifactPaths(input: unknown): string[] {
  const i = (input ?? {}) as { file_path?: unknown; file_paths?: unknown; files?: unknown; root?: unknown }
  const root = typeof i.root === 'string' && i.root ? i.root.replace(/\/$/, '') : ''
  const under = (p: string) => (root && !p.startsWith('/') ? `${root}/${p}` : p)
  const out: string[] = []
  if (typeof i.file_path === 'string') out.push(i.file_path)
  if (Array.isArray(i.file_paths)) for (const p of i.file_paths) if (typeof p === 'string') out.push(p)
  if (Array.isArray(i.files)) {
    for (const f of i.files) {
      if (typeof f === 'string') out.push(under(f))
      else if (f && typeof f === 'object' && typeof (f as { path?: unknown }).path === 'string') out.push(under((f as { path: string }).path))
    }
  } else if (i.files && typeof i.files === 'object') {
    for (const v of Object.values(i.files as Record<string, unknown>)) {
      if (typeof v === 'string') out.push(under(v))
      else if (v && typeof v === 'object' && typeof (v as { from?: unknown }).from === 'string') out.push(under((v as { from: string }).from))
    }
  }
  return out
}

/** The parts of an Artifact call that are published themselves (title, description), not paths. */
const ARTIFACT_PATH_KEYS = new Set(['file_path', 'file_paths', 'files', 'root', 'out_dir', 'url', 'from_url', 'type_url',
  'path', 'paths'])

function artifactFields(input: unknown): Record<string, unknown> {
  const i = (input ?? {}) as Record<string, unknown>
  const out: Record<string, unknown> = {}
  for (const [k, v] of Object.entries(i)) if (!ARTIFACT_PATH_KEYS.has(k)) out[k] = v
  return out
}

/**
 * The identifiers of every other case a shell command names with --case: a clinician with several
 * cases open in turn may query for one while another is active. Unreadable case files add nothing
 * here (the active case's own unreadable file still closes the gate).
 */
async function namedCaseIdentifiers($: EngineInterface, command: string, active: string | null): Promise<string[]> {
  const words = shellWords(command)
  const named = new Set<string>()
  words.forEach((w, i) => {
    if (w === '--case' && words[i + 1]) named.add(words[i + 1] as string)
    else if (w.startsWith('--case=')) named.add(w.slice('--case='.length))
  })
  if (named.size === 0) return []
  const home = (await $.env.get('HOME')) ?? ''
  const out: string[] = []
  for (const raw of named) {
    const dir = raw.replace(/^\$HOME(?=\/|$)/, home).replace(/^~(?=\/|$)/, home).replace(/\/+$/, '')
    if (!dir || dir === active?.replace(/\/+$/, '')) continue
    try {
      const data = JSON.parse(await $.fs.read(`${dir}/case.json`)) as { privacy?: { identifiers?: unknown } }
      const listed = data.privacy?.identifiers
      if (Array.isArray(listed)) out.push(...listed.filter((x): x is string => typeof x === 'string' && x.trim().length >= 2))
    } catch {
      // not a case, or not readable: nothing to add
    }
  }
  return out
}

/**
 * What a page about to be published carries: the case's identifiers and the ID/phone/email
 * patterns, checked in every file it publishes. A file that cannot be read (missing, over the
 * 4 MiB read limit) is reported as unread, so the caller asks instead of passing it unseen.
 */
async function artifactLeak($: EngineInterface, input: unknown, ids: readonly string[]): Promise<{ leak?: string; unread: string[] }> {
  const unread: string[] = []
  const root = (input as { root?: unknown })?.root
  const base = typeof root === 'string' && root ? root.replace(/\/$/, '') : ''
  for (const p of artifactPaths(input)) {
    const path = base && !p.startsWith('/') ? `${base}/${p}` : p
    const text = await $.fs.read(path).catch(() => undefined)
    if (text === undefined) {
      unread.push(p)
      continue
    }
    // paragraph by paragraph (tags, entities and soft line breaks folded first): numbers from all
    // over a page must never join into a record number, but a name split by a line break still counts
    const chunks = text.split(/\n\s*\n|<\/(?:p|div|li|tr|td|th|h[1-6]|section|table)>|<br\s*\/?>/i)
      .map(c => c.replace(/\s*\n\s*/g, ' '))
    const hit = guardInput(chunks, ids)
    if (hit) return { leak: `${hit} (in ${p})`, unread }
  }
  return { unread }
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
    '   /zebra new <dir> [title]   start a case folder (files are written only on this machine)',
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
  const r = (result ?? {}) as {
    checks?: Array<{ name: string; ok: boolean; optional?: boolean; detail?: string }>
    summary?: { ok?: number; failed?: number; optional_missing?: number }
  }
  // ○ marks what zebra works without (research keys, offline data): not a failure
  const rows = (r.checks ?? []).map(c => `${c.ok ? '✓' : c.optional ? '○' : '✗'} ${c.name}${c.detail ? ` — ${c.detail}` : ''}`)
  const s = r.summary
  if (s) rows.push(`${s.ok ?? 0} ok, ${s.failed ?? 0} failed, ${s.optional_missing ?? 0} optional not set up (○)`)
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
  // Pair each source with the evidence id the ledger gave it: one claim, one row to cite.
  const ledger = got.ledger ?? []
  const sources = (got.sources ?? []).map((src, i) =>
    ledger.length === (got.sources ?? []).length && src && typeof src === 'object'
      ? { eid: ledger[i], ...(src as Record<string, unknown>) }
      : src)
  // warnings, ledger and sources first: they are what the answer must cite, so a
  // trim never costs them. The result is trimmed until the whole envelope fits.
  let body = { warnings, ledger, sources, result: got.result }
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
