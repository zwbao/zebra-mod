// What zebra-mod shows in Claude Code, so a person can see when it is at work and what it did:
//   - its own rows for its tool calls (what is looked up, in which databases, how much evidence came back)
//     in place of a raw JSON envelope; a case_update row never shows the identifiers it registers
//   - the spinner says which databases are being queried
//   - a galloping zebra in the band above the prompt while a query runs, and a standing one with a
//     shield for a few seconds after the privacy gate stopped an outgoing call
//   - a footer label (🦓 and the open case) so it is always clear the mod is installed
//   - one line at the end of a turn that used it: calls, databases, evidence rows, calls refused
// Everything is additive and scoped to zebra-mod's own activity: unrelated work draws as before.
// These hooks are registered first, so they wrap the tool, gate and drawing hooks in register.tsx.

import { atom, read, update } from 'claude-code'
import type { EngineInterface, Register, ToolCallResult } from 'claude-code'

import type { Board, ZebraRun, ZebraTurn } from '../types'
import { FRAME_COUNT, RASTER_COLUMNS, RASTER_ROWS, zebraCells } from './ui-sprite'

type On = Parameters<Register>[0]
type Options = Parameters<Register>[1]
export type UiMode = 'full' | 'quiet' | 'off'

const PREFIX = 'mcp__zebra-mod__'
const LOCAL = 'local'
const FRAME_MS = 110

const EMPTY_TURN: ZebraTurn = { calls: 0, failed: 0, blocked: 0, sources: [], cached: 0, evidence: 0, firstEvidence: null, lastEvidence: null }

const board = atom({ plugin: 'zebra-mod', key: 'board' } as const, null as Board | null)
const running = atom({ plugin: 'zebra-mod', key: 'running' } as const, [] as ZebraRun[])
const turnStats = atom({ plugin: 'zebra-mod', key: 'turn' } as const, EMPTY_TURN)
const blockedAt = atom({ plugin: 'zebra-mod', key: 'blockedAt' } as const, null as number | null)
const lang = atom({ plugin: 'zebra-mod', key: 'lang' } as const, null as string | null)
// why the gate refused: an identifier in the call, or a case file it could not read (the gate closed)
const blockedWhy = atom({ plugin: 'zebra-mod', key: 'blockedWhy' } as const, 'identifier' as string)

// ------------------------------------------------------------ words

// [Chinese, English]
const LABELS: Record<string, [string, string]> = {
  case_status: ['病例概况', 'Case status'],
  case_update: ['更新病例', 'Case update'],
  case_recheck: ['病例复查', 'Case recheck'],
  hpo_search: ['HPO 术语检索', 'HPO term search'],
  phenotype_rank: ['表型排序', 'Phenotype ranking'],
  gene_card: ['基因卡', 'Gene card'],
  variant_card: ['变异卡', 'Variant card'],
  disease_card: ['疾病卡', 'Disease card'],
  acmg: ['ACMG 计分', 'ACMG points'],
  cnv_interpret: ['拷贝数变异解读', 'CNV reading'],
  s2f_predict: ['序列功能预测', 'Sequence-to-function'],
  therapy_landscape: ['治疗全景', 'Therapy landscape'],
  trials_search: ['临床试验检索', 'Trial search'],
  literature_search: ['文献检索', 'Literature search'],
  rare_stats: ['统计计算', 'Statistics'],
  edit_check: ['碱基编辑可行性', 'Base-editing screen'],
  china_rare: ['中国罕见病目录', 'China rare disease lists'],
  access: ['药物可及性', 'Drug access'],
  expression: ['组织表达', 'Tissue expression'],
  aso_screen: ['反义寡核苷酸筛查', 'Antisense screen'],
  report_export: ['导出报告', 'Report export'],
}

// What each tool is expected to query, shown while it runs; the row afterwards names what it did query.
const EXPECTED: Record<string, string[]> = {
  case_status: [LOCAL], case_update: [LOCAL], report_export: [LOCAL], rare_stats: [LOCAL],
  case_recheck: ['ClinVar', 'ClinGen', 'ClinicalTrials.gov', 'Europe PMC'],
  hpo_search: ['HPO'],
  phenotype_rank: ['HPO', 'Monarch', 'PubCaseFinder'],
  gene_card: ['HGNC', 'ClinGen', 'gnomAD', 'PanelApp', 'UniProt'],
  variant_card: ['Ensembl VEP', 'gnomAD', 'ClinVar', 'LitVar2'],
  disease_card: ['Orphanet', 'Monarch', 'GeneReviews', 'HPO'],
  cnv_interpret: ['Ensembl', 'ClinGen'],
  s2f_predict: ['SpliceAI', 'Pangolin'],
  therapy_landscape: ['Open Targets', 'EMA'],
  trials_search: ['ClinicalTrials.gov'],
  literature_search: ['Europe PMC', 'PubTator3', 'LitVar2'],
  edit_check: ['Ensembl'],
  china_rare: ['China national rare disease list'],
  access: ['Open Targets', 'EMA', 'China drug approvals', '国家医保药品目录', 'ClinicalTrials.gov'],
  expression: ['GTEx'],
  aso_screen: ['Ensembl', 'SpliceAI', 'Europe PMC'],
}

// Source names whose short form is not their first words.
const RENAMED: Array<[RegExp, string, string]> = [
  [/^China drug approvals/, 'NMPA/CDE', 'NMPA/CDE'],
  [/^China national rare disease list/, '国家罕见病目录', 'China rare disease lists'],
  [/^China rare disease list alias/, '', ''],
  [/^国家医保药品目录/, '国家医保目录', 'China NRDL'],
  [/^全国罕见病诊疗协作网/, '罕见病诊疗协作网', 'China rare disease network'],
  [/^GeneReviews/, 'GeneReviews', 'GeneReviews'],
  [/^local$/, '本机', 'this computer'],
]

/** "ClinVar esearch" → ClinVar, "HPO annotations (phenotype.hpoa)" → HPO: the database, not the endpoint. */
export function sourceName(db: string, zh: boolean): string {
  for (const [re, z, en] of RENAMED) if (re.test(db)) return zh ? z : en
  const head = (db.split(' (')[0] ?? db).trim()
  const kept: string[] = []
  for (const word of head.split(/\s+/)) {
    if (kept.length > 0 && /^[a-z]/.test(word)) break
    kept.push(word)
  }
  return kept.join(' ') || head
}

function sourceNames(dbs: readonly string[], zh: boolean): string[] {
  return [...new Set(dbs.map(d => sourceName(d, zh)).filter(Boolean))]
}

export function labelOf(tool: string, zh: boolean): string {
  const l = LABELS[tool]
  return l ? (zh ? l[0] : l[1]) : tool
}

function clip(text: string, max: number): string {
  const one = text.replace(/\s+/g, ' ').trim()
  return one.length > max ? `${one.slice(0, Math.max(1, max - 1))}…` : one
}

/** What a call is about, for its row. Never the identifiers a case_update registers. */
export function queryOf(tool: string, input: unknown, zh: boolean): string {
  const i = (input && typeof input === 'object' ? input : {}) as Record<string, unknown>
  const s = (k: string): string => (typeof i[k] === 'string' ? (i[k] as string).trim() : '')
  const n = (k: string): number => (Array.isArray(i[k]) ? (i[k] as unknown[]).length : 0)
  // the first plain argument, when the usual ones are absent (never the reserved keys or identifiers)
  const any = (): string => {
    for (const [k, v] of Object.entries(i)) {
      if (['tool', 'tool_use_id', 'consent', 'agentId', 'identifiers'].includes(k)) continue
      if (typeof v === 'string' && v.trim()) return v.trim()
      if (Array.isArray(v) && v.length && v.every(x => typeof x === 'string')) return (v as string[]).join(' ')
    }
    return ''
  }
  return pick() || (tool === 'case_update' || tool === 'case_status' || tool === 'case_recheck' ? '' : any())
  function pick(): string {
  switch (tool) {
    case 'case_update': {
      const parts: string[] = []
      const ids = n('identifiers')
      if (ids) parts.push(zh ? `登记受保护身份信息 ${ids} 项（内容不显示）` : `${ids} protected identifiers registered (not shown)`)
      const kinds: Array<[string, string, string]> = [
        ['phenotypes', '表型', 'phenotypes'], ['variants', '变异', 'variants'], ['hypotheses', '诊断假设', 'hypotheses'],
        ['acmg', 'ACMG 判读', 'ACMG readings'], ['leads', '治疗线索', 'therapy leads'], ['tests', '已做检查', 'tests'],
        ['family', '家庭成员', 'relatives'], ['timeline', '时间线', 'timeline events'], ['questions', '待问问题', 'questions'],
      ]
      for (const [k, z, en] of kinds) {
        const c = n(k)
        if (c) parts.push(zh ? `${z} ${c} 项` : `${c} ${en}`)
      }
      if (i.profile && typeof i.profile === 'object') parts.push(zh ? '基本信息' : 'profile')
      if (n('remove')) parts.push(zh ? `删除 ${n('remove')} 项` : `${n('remove')} removals`)
      return parts.join(' · ')
    }
    case 'case_status':
    case 'case_recheck':
      return ''
    case 'hpo_search':
      return s('text')
    case 'phenotype_rank': {
      if (i.from_case === true) return zh ? '使用病例中的表型' : 'phenotypes from the case'
      const p = n('present')
      const x = n('excluded')
      return zh ? `${p} 项表型${x ? ` · 排除 ${x} 项` : ''}` : `${p} present${x ? ` · ${x} excluded` : ''}`
    }
    case 'gene_card':
      return s('symbol')
    case 'acmg': {
      const codes = Array.isArray(i.codes) ? (i.codes as unknown[]).filter((c): c is string => typeof c === 'string') : []
      return codes.length ? codes.join(' + ') : s('variant')
    }
    case 'trials_search':
      return [s('condition'), s('term'), s('country')].filter(Boolean).join(' · ')
    case 'literature_search':
      return s('query') || [s('gene'), s('variant')].filter(Boolean).join(' ')
    case 'rare_stats':
      return s('method')
    case 'cnv_interpret':
      return s('result')
    case 'report_export':
      return s('report')
    default: {
      for (const k of ['variant', 'query', 'gene', 'symbol', 'text']) if (s(k)) return s(k)
      return ''
    }
  }
  }
}

// ------------------------------------------------------------ the tool's answer

type Envelope = { warnings?: unknown; ledger?: unknown; sources?: unknown; result?: unknown }

/** The JSON envelope register.tsx answers a zebra tool with, from a ToolResult's output, or null. */
export function envelopeOf(output: unknown): Envelope | null {
  let text: string | undefined
  if (typeof output === 'string') text = output
  else if (Array.isArray(output)) {
    text = output.map(b => (b && typeof b === 'object' && typeof (b as { text?: unknown }).text === 'string' ? (b as { text: string }).text : '')).join('')
  } else if (output && typeof output === 'object') {
    if ('ledger' in output || 'sources' in output) return output as Envelope
    const t = (output as { text?: unknown }).text
    if (typeof t === 'string') text = t
  }
  if (!text || text[0] !== '{') return null
  try {
    const parsed = JSON.parse(text) as unknown
    return parsed && typeof parsed === 'object' && ('ledger' in parsed || 'sources' in parsed) ? (parsed as Envelope) : null
  } catch {
    return null
  }
}

type Digest = { dbs: string[]; cached: number; evidence: number; first: number | null; last: number | null; warnings: number }

export function digest(env: Envelope): Digest {
  const sources = Array.isArray(env.sources) ? env.sources : []
  const dbs: string[] = []
  let cached = 0
  for (const s of sources) {
    if (!s || typeof s !== 'object') continue
    const db = (s as { db?: unknown }).db
    if (typeof db === 'string' && db.trim()) dbs.push(db.trim())
    if ((s as { cached?: unknown }).cached === true) cached++
  }
  const ids = (Array.isArray(env.ledger) ? env.ledger : [])
    .map(x => (typeof x === 'string' ? /^E(\d+)$/.exec(x.trim()) : null))
    .filter((m): m is RegExpExecArray => m !== null)
    .map(m => Number(m[1]))
  const warnings = Array.isArray(env.warnings) ? env.warnings.length : 0
  return {
    dbs: [...new Set(dbs)],
    cached,
    evidence: ids.length,
    first: ids.length ? Math.min(...ids) : null,
    last: ids.length ? Math.max(...ids) : null,
    warnings,
  }
}

function range(first: number | null, last: number | null): string {
  if (first === null || last === null) return ''
  return first === last ? `E${first}` : `E${first}–E${last}`
}

function listed(names: string[], zh: boolean, max = 4): string {
  if (names.length <= max) return names.join(zh ? '、' : ', ')
  const shown = names.slice(0, max).join(zh ? '、' : ', ')
  return zh ? `${shown} 等 ${names.length} 个` : `${shown} +${names.length - max} more`
}

/** The one line under a zebra tool's row: evidence, databases, cache, notes. */
export function summaryLine(env: Envelope, zh: boolean): string {
  const d = digest(env)
  const names = sourceNames(d.dbs, zh)
  const parts: string[] = []
  if (d.evidence) parts.push(zh ? `证据 ${d.evidence} 条（${range(d.first, d.last)}）` : `${d.evidence} evidence rows (${range(d.first, d.last)})`)
  if (names.length) parts.push(zh ? `来源 ${listed(names, zh)}` : `from ${listed(names, zh)}`)
  if (d.cached) parts.push(zh ? `${d.cached} 条来自缓存` : `${d.cached} cached`)
  if (d.warnings) parts.push(zh ? `${d.warnings} 条提示` : `${d.warnings} notes`)
  if (!parts.length) return zh ? '已在本机完成' : 'done on this computer'
  return parts.join(' · ')
}

/** The line left at the end of a turn that used zebra-mod. */
export function turnLine(t: ZebraTurn, zh: boolean): string | undefined {
  if (t.calls === 0 && t.blocked === 0) return undefined
  if (t.calls === 0) {
    return zh
      ? `🛡 本轮隐私闸门拦下 ${t.blocked} 次外发调用，相关内容未发出。`
      : `🛡 The privacy gate stopped ${t.blocked} outgoing call${t.blocked > 1 ? 's' : ''} this turn; nothing in them was sent.`
  }
  const names = sourceNames(t.sources, zh).filter(n => n !== (zh ? '本机' : 'this computer'))
  const parts: string[] = [zh ? `调用 ${t.calls} 次` : `${t.calls} call${t.calls > 1 ? 's' : ''}`]
  if (names.length) parts.push(zh ? `查询 ${names.length} 个数据源（${listed(names, zh, 5)}）` : `${names.length} databases (${listed(names, zh, 5)})`)
  if (t.evidence) parts.push(zh ? `新增证据 ${t.evidence} 条（${range(t.firstEvidence, t.lastEvidence)}）` : `${t.evidence} evidence rows (${range(t.firstEvidence, t.lastEvidence)})`)
  if (t.failed) parts.push(zh ? `${t.failed} 次未执行` : `${t.failed} not run`)
  if (t.blocked) parts.push(zh ? `隐私闸门拦截 ${t.blocked} 次` : `${t.blocked} stopped by the privacy gate`)
  return `🦓 ${zh ? '本轮 zebra-mod：' : 'zebra-mod this turn: '}${parts.join(' · ')}`
}

// ------------------------------------------------------------ hooks

// What the animation and the language choice read between render passes; a reload resets them.
let envZh: boolean | undefined
let bandId: string | undefined
let frame = 0
let scroll = 0
let ticker: { cancel: () => void; owner: string } | undefined

async function isZh($: EngineInterface): Promise<boolean> {
  const b = await read($, board)
  if (b?.language) return b.language.toLowerCase().startsWith('zh')
  const l = await read($, lang)
  if (l) return l === 'zh'
  if (envZh === undefined) {
    const env = (await $.env.get('LC_ALL')) || (await $.env.get('LANG')) || ''
    envZh = env.toLowerCase().startsWith('zh')
  }
  return envZh
}

// The gallop runs inside the tool call that started it (that dispatch lasts as long as the query)
// and stops with it.
function startTicker($: EngineInterface, owner: string): void {
  if (ticker) return
  let timer: { cancel: () => void }
  try {
    timer = $.clock.every(FRAME_MS, () => {
    frame = (frame + 1) % FRAME_COUNT
    scroll++
    if (bandId === undefined) return
    void $.ui
      .blit({ requestId: bandId, key: 'zebra', cells: zebraCells(frame, scroll), columns: RASTER_COLUMNS, rows: RASTER_ROWS })
      .catch(() => undefined)
    })
  } catch {
    return // no timer here: the band stays a still picture
  }
  ticker = { cancel: () => timer.cancel(), owner }
}

function stopTicker(owner: string): void {
  if (ticker?.owner !== owner) return
  ticker.cancel()
  ticker = undefined
}

/** One zebra tool call, wrapped: in flight for the spinner and the band, then tallied for the turn. */
async function observeCall<E extends { tool_use_id?: string }>(
  $: EngineInterface,
  e: E,
  next: (e: E) => Promise<ToolCallResult>,
  tool: string,
  mode: UiMode,
): Promise<ToolCallResult> {
  const now = Date.now()
  const id = String(e.tool_use_id ?? `${tool}:${now}`)
  const run: ZebraRun = {
    id,
    tool,
    queryZh: clip(queryOf(tool, e, true), 80),
    queryEn: clip(queryOf(tool, e, false), 80),
    sources: EXPECTED[tool] ?? [],
    startedAt: now,
  }
  await update($, running, list => [...list, run])
  if (mode === 'full' && run.sources.some(s => s !== LOCAL)) startTicker($, id)
  let answer: ToolCallResult | undefined
  try {
    answer = await next(e)
    return answer
  } finally {
    stopTicker(id)
    await update($, running, list => list.filter(r => r.id !== id))
    await update($, turnStats, t => recordCall(t, tool, answer))
  }
}

/** A tool-result row that carries the privacy gate's refusal: a shield until the next turn. */
async function noteRefusals($: EngineInterface, content: unknown, mode: UiMode): Promise<void> {
  const texts = gateRefusals(content)
  const refusals = texts.length
  if (refusals === 0) return
  const at = Date.now()
  await update($, blockedWhy, () => (texts.some(t => t.includes('gate is closed')) ? 'closed' : 'identifier'))
  await update($, blockedAt, () => at)
  await update($, turnStats, t => ({ ...t, blocked: t.blocked + refusals }))
  if (mode === 'quiet') {
    $.ui.toast(
      (await isZh($))
        ? '🛡 隐私闸门拦下一次外发调用：其中含有受保护的身份信息，未发出。'
        : '🛡 The privacy gate stopped an outgoing call that carried a protected identifier; it was not sent.',
    )
  }
}

async function welcomeOnce($: EngineInterface): Promise<void> {
  if ((await $.store.get('uiWelcomed')) === true) return
  await $.store.set('uiWelcomed', true)
  $.ui.toast(
    (await isZh($))
      ? '🦓 zebra-mod 已启用：遇到罕见病问题时会调用它的研究工具。输入 /zebra 查看用法。'
      : '🦓 zebra-mod is on: rare-disease questions get its research tools. Type /zebra for help.',
  )
}

async function closeTurn($: EngineInterface): Promise<void> {
  const line = turnLine(await read($, turnStats), await isZh($))
  if (line) await $.session.append({ message: { type: 'system', content: [{ type: 'text', text: line }] } }).catch(() => undefined)
  await update($, turnStats, () => EMPTY_TURN)
  await update($, running, () => [])
}

export function registerUi(on: On, options: Options): void {
  const mode: UiMode = options.interface === 'quiet' || options.interface === 'off' ? options.interface : 'full'
  if (mode === 'off') return

  on('turn.start', async ($, e, next) => {
    await update($, turnStats, () => EMPTY_TURN)
    await update($, blockedAt, () => null)
    await update($, running, () => []) // nothing of a previous turn is still running
    if (/[㐀-鿿]/.test(e.text)) await update($, lang, () => 'zh')
    await welcomeOnce($)
    return next(e)
  })

  on('turn.complete', async ($, e, next) => {
    const done = await next(e)
    if (!e.agentId) await closeTurn($)
    return done
  })

  // Wraps the tool hook in register.tsx (registered after these): what runs, and what it brought back.
  for (const tool of Object.keys(LABELS)) {
    on('tool.call', { tool: `${PREFIX}${tool}` }, async ($, e, next) => observeCall($, e, next, tool, mode)).catch(($, e, next) => next(e))
  }

  // The privacy gate's refusals, from the tool results they leave (any tool, zebra's or not).
  on('session.append', async ($, e, next) => {
    const stored = await next(e)
    if (e.door === 'tool-result') await noteRefusals($, e.message.content, mode)
    return stored
  }).catch(($, e, next) => next(e))

  // ---------------------------------------------------------- drawing

  // A zebra tool's row: 🦓, what it is, what it is about; while it runs, the databases it queries.
  on('ui.render', { component: 'ToolUse' }, async ($, e, next) => {
    const tool = e.props.tool.startsWith(PREFIX) ? e.props.tool.slice(PREFIX.length) : ''
    if (!LABELS[tool] || e.props.isInterrupted) return next(e)
    const zh = await isZh($)
    const { Box, Text } = $.ui.resolve(e)
    const width = e.viewport?.columns ?? 100
    const label = labelOf(tool, zh)
    const query = clip(queryOf(tool, e.props.input, zh), Math.max(16, width - label.length - 16))
    const head = [<Text color="claude">🦓 </Text>, <Text bold>{label}</Text>]
    if (query) head.push(<Text dimColor>{`  ${query}`}</Text>)
    if (e.props.isRunning) head.push(<Text color="subtle">{zh ? '  查询中…' : '  running…'}</Text>)
    else if (e.props.isErrored) head.push(<Text color="error">{zh ? '  未执行' : '  not run'}</Text>)
    const rows = [<Box flexDirection="row">{head}</Box>]
    const expected = sourceNames(EXPECTED[tool] ?? [], zh)
    if (e.props.isRunning && expected.length) {
      rows.push(<Text dimColor wrap="truncate-end">{`   ${zh ? '数据源：' : 'sources: '}${expected.join(' · ')}`}</Text>)
    }
    return <Box flexDirection="column">{rows}</Box>
  })

  // Its result: evidence, databases, cache and notes in one line, in place of the JSON envelope.
  on('ui.render', { component: 'ToolResult' }, async ($, e, next) => {
    const tool = e.props.tool.startsWith(PREFIX) ? e.props.tool.slice(PREFIX.length) : ''
    if (!LABELS[tool] || e.props.isErrored) return next(e)
    const env = envelopeOf(e.props.output)
    if (!env) return next(e)
    const { Text } = $.ui.resolve(e)
    return <Text dimColor wrap="truncate-end">{`  ⎿ ${summaryLine(env, await isZh($))}`}</Text>
  })

  // The spinner names the databases being queried.
  on('ui.render', { component: 'Spinner' }, async ($, e, next) => {
    if (mode !== 'full' || e.props.message) return next(e)
    const runs = await read($, running)
    if (!runs.length) return next(e)
    const zh = await isZh($)
    const names = sourceNames(runs.flatMap(r => r.sources), zh).filter(n => n !== (zh ? '本机' : 'this computer'))
    const word = names.length
      ? `🦓 ${zh ? '正在查询' : 'Querying'} ${listed(names, zh, 3)}`
      : `🦓 ${labelOf((runs[0] as ZebraRun).tool, zh)}`
    return next({ ...e, props: { ...e.props, word } })
  })

  // A footer label: the mod is here (🦓), with the open case's title, or a shield after a refusal.
  on('ui.render', { component: 'SessionMode' }, async ($, e, next) => {
    const b = await read($, board)
    const at = await read($, blockedAt)
    const label = `${at !== null ? '🛡' : '🦓'} ${b ? clip(b.title, 20) : 'zebra'}`
    return next({ ...e, props: { ...e.props, modes: [...e.props.modes, label] } })
  })

  // The band above the prompt: a galloping zebra while a query runs; a standing one with a shield
  // after the gate refused a call, until the next turn. Nothing otherwise.
  on('ui.render', { component: 'AbovePrompt' }, async ($, e, next) => {
    if (mode !== 'full' || e.props.hasSurvey) return next(e)
    const runs = (await read($, running)).filter(r => r.sources.some(s => s !== LOCAL))
    const at = await read($, blockedAt)
    if (!runs.length && at === null) return next(e)
    const zh = await isZh($)
    bandId = e.requestId
    const lines = runs.length ? runningLines(runs, zh) : shieldLines(zh, (await read($, blockedWhy)) === 'closed')
    const tone = runs.length ? 'claude' : 'success'
    if (e.surface === 'terminal' && e.props.bodyColumns >= RASTER_COLUMNS + 30) {
      const { Box, Text, Raster } = $.ui.resolve(e)
      return (
        <Box flexDirection="row" gap={2}>
          <Raster key="zebra" columns={RASTER_COLUMNS} rows={RASTER_ROWS} cells={zebraCells(frame, scroll, runs.length === 0)} />
          <Box flexDirection="column" justifyContent="center" width={Math.max(20, e.props.bodyColumns - RASTER_COLUMNS - 2)}>
            <Text bold color={tone} wrap="truncate-end">{lines[0]}</Text>
            <Text wrap="truncate-end">{lines[1]}</Text>
            <Text dimColor wrap="truncate-end">{lines[2]}</Text>
          </Box>
        </Box>
      )
    }
    const { Box, Text } = $.ui.resolve(e)
    return (
      <Box flexDirection="column">
        <Text bold color={tone} wrap="truncate-end">{lines[0]}</Text>
        <Text wrap="truncate-end">{lines[1]}</Text>
        <Text dimColor wrap="truncate-end">{lines[2]}</Text>
      </Box>
    )
  })
}

/** The band's three lines while queries run. */
export function runningLines(runs: readonly ZebraRun[], zh: boolean): [string, string, string] {
  const first = runs[0] as ZebraRun
  const query = zh ? first.queryZh : first.queryEn
  const what = runs.length === 1
    ? `${labelOf(first.tool, zh)}${query ? `${zh ? '：' : ': '}${query}` : ''}`
    : (zh ? `${runs.length} 项查询同时进行` : `${runs.length} queries at once`)
  return [
    zh ? '🦓 zebra-mod 正在查询公开数据库' : '🦓 zebra-mod is querying public databases',
    what,
    `${zh ? '数据源：' : 'sources: '}${sourceNames(runs.flatMap(r => r.sources), zh).join(' · ')}`,
  ]
}

/** The band's three lines after the privacy gate refused a call. */
export function shieldLines(zh: boolean, isClosed = false): [string, string, string] {
  if (isClosed) {
    return [
      zh ? '🛡 隐私闸门已关闭外发调用' : '🛡 The privacy gate is holding outgoing calls',
      zh ? '病例的受保护身份信息名单无法读取，调用未发出。' : "The case's list of protected identifiers cannot be read; the call was not sent.",
      zh ? '修复或重建 case.json，或用 /zebra close 关闭病例后再查询。' : 'Fix or re-create case.json, or /zebra close, then query again.',
    ]
  }
  return [
    zh ? '🛡 隐私闸门拦下了一次外发调用' : '🛡 The privacy gate stopped an outgoing call',
    zh ? '调用中含有受保护的身份信息，内容未发出。' : 'It carried a protected identifier; it was not sent.',
    zh ? '查询公开数据库时只使用 HPO 编号、基因、变异和疾病名称。' : 'Public databases are queried with HPO ids, genes, variants and disease names only.',
  ]
}

/** The privacy gate's refusal texts among a tool-result row's blocks (register.tsx words them). */
export function gateRefusals(content: unknown): string[] {
  const out: string[] = []
  const visit = (block: unknown): void => {
    if (typeof block === 'string') {
      if (block.startsWith('zebra-mod privacy gate:') && !block.includes('confirm.')) out.push(block)
      return
    }
    if (Array.isArray(block)) return block.forEach(visit)
    if (!block || typeof block !== 'object') return
    const b = block as { type?: unknown; text?: unknown; content?: unknown; is_error?: unknown }
    if (b.type === 'tool_result' && b.is_error !== true) return
    if (typeof b.text === 'string') visit(b.text)
    if (b.content !== undefined) visit(b.content)
  }
  visit(content)
  return out
}

/** The turn's tally after one zebra call answered `answer`. */
export function recordCall(t: ZebraTurn, tool: string, answer: ToolCallResult | undefined): ZebraTurn {
  const next: ZebraTurn = { ...t, calls: t.calls + 1, sources: [...t.sources] }
  if (!answer || 'deny' in answer && answer.deny !== undefined) {
    next.failed++
    return next
  }
  const env = envelopeOf(answer.result)
  if (!env) {
    if ((EXPECTED[tool] ?? []).includes(LOCAL) && !next.sources.includes(LOCAL)) next.sources.push(LOCAL)
    return next
  }
  const d = digest(env)
  for (const db of d.dbs) if (!next.sources.includes(db)) next.sources.push(db)
  next.cached += d.cached
  next.evidence += d.evidence
  if (d.first !== null) next.firstEvidence = next.firstEvidence === null ? d.first : Math.min(next.firstEvidence, d.first)
  if (d.last !== null) next.lastEvidence = next.lastEvidence === null ? d.last : Math.max(next.lastEvidence, d.last)
  return next
}
