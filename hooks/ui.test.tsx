import { describe, expect, test } from 'claude-code/testing'

import { base64, FRAME_COUNT, RASTER_COLUMNS, RASTER_ROWS, zebraCells } from './ui-sprite'
import { envelopeOf, gateRefusals, queryOf, recordCall, shieldLines, sourceName, summaryLine, turnLine } from './ui'

const P = 'mcp__zebra-mod__'
const ENVELOPE = JSON.stringify({
  warnings: ['2 genes overlap this position'],
  ledger: ['E39', 'E40', 'E41', 'E53'],
  sources: [
    { db: 'Ensembl VEP', cached: false }, { db: 'gnomAD', cached: true }, { db: 'ClinVar esearch' }, { db: 'ClinVar' },
    { db: 'LitVar2 autocomplete' }, { db: 'HPO annotations (phenotype.hpoa)' }, { db: 'China drug approvals (NMPA/CDE/MOF documents)' },
  ],
  result: {},
})
const TURN0 = { calls: 0, failed: 0, blocked: 0, sources: [], cached: 0, evidence: 0, firstEvidence: null, lastEvidence: null }
const SURFACES = ['terminal', 'desktop'] as const

function toolRow(tool: string, input: unknown, over: Record<string, unknown> = {}) {
  return { tool_use_id: 'tu1', tool, input, isRunning: false, isErrored: false, isInterrupted: false, ...over }
}

/** What the engine would draw beneath the plugin: one Text, so a pass-through is visible. */
function engineDraws(on: any, component: string, text: string, seen?: (e: any) => void) {
  on('ui.render', { component }, ($$: any, e: any) => {
    seen?.(e)
    const { Text } = $$.ui.resolve(e)
    return <Text>{text}</Text>
  })
}

function noState(on: any, key: string, value: unknown) {
  on('state.get', { plugin: 'zebra-mod', key }, () => ({ value: { value, version: 1 } }))
}

describe('words', () => {
  test('a source is named by its database, not its endpoint', async () => {
    expect(sourceName('ClinVar esearch', false)).toBe('ClinVar')
    expect(sourceName('HPO annotations (phenotype.hpoa)', false)).toBe('HPO')
    expect(sourceName('gnomAD constraint', false)).toBe('gnomAD')
    expect(sourceName('Ensembl VEP', false)).toBe('Ensembl VEP')
    expect(sourceName('Monarch gene-disease associations', false)).toBe('Monarch')
    expect(sourceName('China drug approvals (NMPA/CDE/MOF documents)', true)).toBe('NMPA/CDE')
    expect(sourceName('国家医保药品目录 (NHSA)', true)).toBe('国家医保目录')
    expect(sourceName('China rare disease list alias layer', true)).toBe('')
  })

  test('a case_update row never shows the identifiers it registers', async () => {
    const input = { identifiers: ['王小雨', 'MZ20240518007', '2022-02-14'], phenotypes: [{ id: 'HP:0001250' }] }
    for (const zh of [true, false]) {
      const q = queryOf('case_update', input, zh)
      expect(q.includes('王小雨')).toBe(false)
      expect(q.includes('MZ20240518007')).toBe(false)
      expect(q.includes('2022')).toBe(false)
      expect(q.includes('3')).toBe(true)
    }
  })

  test('a row says what the call is about even with unusual arguments, never the reserved keys', async () => {
    expect(queryOf('literature_search', { pmids: ['34867351', '41623181'], tool: 'x', tool_use_id: 'y' }, true)).toBe('34867351 41623181')
    expect(queryOf('literature_search', { gene: 'SCN8A' }, false)).toBe('SCN8A')
    expect(queryOf('case_status', { tool: 'mcp__zebra-mod__case_status' }, false)).toBe('')
  })

  test('the result line counts evidence, names databases once, says what was cached', async () => {
    const env = envelopeOf(ENVELOPE)
    expect(env).not.toBe(null)
    const en = summaryLine(env!, false)
    expect(en.includes('4 evidence rows (E39–E53)')).toBe(true)
    expect(en.includes('Ensembl VEP, gnomAD, ClinVar, LitVar2 +2 more')).toBe(true)
    expect(en.includes('1 cached')).toBe(true)
    expect(en.includes('1 notes')).toBe(true)
    const zh = summaryLine(env!, true)
    expect(zh.startsWith('证据 4 条（E39–E53）')).toBe(true)
    expect(summaryLine({ ledger: [], sources: [] }, true)).toBe('已在本机完成')
  })

  test('an envelope is read from a string, text blocks or an object; anything else is not one', async () => {
    expect(envelopeOf([{ type: 'text', text: ENVELOPE }])).not.toBe(null)
    expect(envelopeOf({ text: ENVELOPE })).not.toBe(null)
    expect(envelopeOf('zebra variant_card failed: NoOutput')).toBe(null)
    expect(envelopeOf('{"not": "ours"}')).toBe(null)
  })

  test('the turn line adds up calls, databases, evidence and refusals', async () => {
    let t = recordCall(TURN0, 'variant_card', { result: ENVELOPE })
    t = recordCall(t, 'acmg', { deny: 'acmg was not approved' })
    t = { ...t, blocked: 1 }
    const line = turnLine(t, false) ?? ''
    expect(line.startsWith('🦓 zebra-mod this turn: 2 calls')).toBe(true)
    expect(line.includes('4 evidence rows (E39–E53)')).toBe(true)
    expect(line.includes('1 not run')).toBe(true)
    expect(line.includes('1 stopped by the privacy gate')).toBe(true)
    expect(turnLine(TURN0, true)).toBe(undefined)
    expect((turnLine({ ...TURN0, blocked: 2 }, true) ?? '').startsWith('🛡 本轮隐私闸门拦下 2 次')).toBe(true)
  })

  test('the shield says why: an identifier in the call, or a gate closed on an unreadable case', async () => {
    expect(shieldLines(true)[1]).toBe('调用中含有受保护的身份信息，内容未发出。')
    expect(shieldLines(true, true)[1].includes('无法读取')).toBe(true)
    expect(shieldLines(false, true)[0]).toBe('🛡 The privacy gate is holding outgoing calls')
    const closed = 'zebra-mod privacy gate: /x/case.json cannot be read, so the protected identifiers are unknown and the gate is closed. Fix or re-create the case file, or /zebra close to work without one.'
    expect(gateRefusals([{ type: 'tool_result', is_error: true, content: closed }])).toHaveLength(1)
  })

  test('only the gate’s refusals count, not its questions', async () => {
    const deny = 'zebra-mod privacy gate: this call would send protected identifier #1 from the case off this machine. Query public databases with HPO ids, gene symbols, variants and disease names only.'
    const ask = 'zebra-mod privacy gate: this call would send an email address off this machine. If it belongs to a patient, do not send it; if it is yours or public, confirm.'
    expect(gateRefusals([{ type: 'tool_result', is_error: true, content: [{ type: 'text', text: deny }] }])).toHaveLength(1)
    expect(gateRefusals([{ type: 'tool_result', is_error: true, content: ask }])).toHaveLength(0)
    expect(gateRefusals([{ type: 'tool_result', is_error: false, content: deny }])).toHaveLength(0)
  })
})

describe('the galloping zebra', () => {
  test('every frame is a whole Raster of braille cells in the terminal’s colour', async () => {
    const bytes = RASTER_COLUMNS * RASTER_ROWS * 12
    const decode = (b64: string): Uint32Array => {
      const bin = atob(b64)
      const u8 = new Uint8Array(bin.length)
      for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i)
      return new Uint32Array(u8.buffer)
    }
    for (let f = 0; f < FRAME_COUNT; f++) {
      const cells = zebraCells(f, f * 2)
      expect(cells.length).toBe(Math.ceil(bytes / 3) * 4)
      const words = decode(cells)
      let horse = 0
      for (let i = 0; i < words.length; i += 3) {
        const cp = words[i] as number
        expect(cp === 0x20 || (cp >= 0x2800 && cp <= 0x28ff)).toBe(true)
        expect(words[i + 2]).toBe(0x01000000) // background left to the terminal
        if (words[i + 1] === 0x01000000 && cp !== 0x20) horse++
      }
      expect(horse > 40).toBe(true)
    }
    expect(FRAME_COUNT).toBe(11)
    expect(zebraCells(0, 0)).not.toBe(zebraCells(1, 0))
    expect(zebraCells(0, 0)).not.toBe(zebraCells(0, 2)) // the ground runs
    expect(zebraCells(0, 0, true)).toBe(zebraCells(2, 7, true))
  })

  test('the body stays level from frame to frame: only the legs and head move', async () => {
    const decode = (b64: string): Uint32Array => {
      const bin = atob(b64)
      const u8 = new Uint8Array(bin.length)
      for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i)
      return new Uint32Array(u8.buffer)
    }
    const BITS: Array<[number, number, number]> = [[0, 0, 0x01], [0, 1, 0x02], [0, 2, 0x04], [1, 0, 0x08], [1, 1, 0x10], [1, 2, 0x20], [0, 3, 0x40], [1, 3, 0x80]]
    const backs: number[] = []
    for (let f = 0; f < FRAME_COUNT; f++) {
      const w = decode(zebraCells(f, 0))
      const W = RASTER_COLUMNS * 2
      const H = RASTER_ROWS * 4
      const dot = Array.from({ length: H }, () => new Array<boolean>(W).fill(false))
      for (let r = 0; r < RASTER_ROWS; r++) {
        for (let c = 0; c < RASTER_COLUMNS; c++) {
          const i = (r * RASTER_COLUMNS + c) * 3
          if (w[i + 1] !== 0x01000000 || (w[i] as number) < 0x2800) continue // the zebra, not the ground
          const b = (w[i] as number) - 0x2800
          for (const [dx, dy, bit] of BITS) if (b & bit) dot[r * 4 + dy]![c * 2 + dx] = true
        }
      }
      const cols = [...Array(W).keys()].filter(x => dot.some(row => row[x]))
      const xa = cols[0]!
      const span = cols[cols.length - 1]! - xa + 1
      const tops: number[] = []
      for (let x = xa + Math.floor(span * 0.35); x < xa + Math.floor(span * 0.6); x++) {
        const y = dot.findIndex(row => row[x])
        if (y >= 0) tops.push(y)
      }
      tops.sort((a, b) => a - b)
      backs.push(tops[Math.floor(tops.length / 2)]!)
    }
    expect(Math.max(...backs) - Math.min(...backs)).toBe(0)
  })

  test('base64 matches the standard alphabet and padding', async () => {
    const enc = (s: string) => base64(new TextEncoder().encode(s))
    expect(enc('')).toBe('')
    expect(enc('f')).toBe('Zg==')
    expect(enc('fo')).toBe('Zm8=')
    expect(enc('foo')).toBe('Zm9v')
    expect(enc('zebra-mod')).toBe('emVicmEtbW9k')
  })
})

describe('drawing', () => {
  test('a zebra tool gets its own row on every surface; another tool is left to the engine', async ($, on) => {
    noState(on, 'board', null)
    noState(on, 'lang', 'en')
    engineDraws(on, 'ToolUse', 'engine row')
    for (const surface of SURFACES) {
      const ui = await $.ui.mount({
        plugin: 'zebra-mod', surface, component: 'ToolUse',
        props: toolRow(`${P}variant_card`, { variant: 'NM_001165963.4:c.2134C>T' }, { isRunning: true }),
      })
      expect(await ui.find({ type: 'Text', text: 'Variant card' })).toBeDefined()
      expect(await ui.find({ type: 'Text', text: /NM_001165963\.4:c\.2134C>T/ })).toBeDefined()
      expect(await ui.find({ type: 'Text', text: /sources: Ensembl VEP · gnomAD · ClinVar · LitVar2/ })).toBeDefined()
      await ui.unmount()

      const bash = await $.ui.mount({ plugin: 'zebra-mod', surface, component: 'ToolUse', props: toolRow('Bash', { command: 'ls' }) })
      expect(await bash.find({ type: 'Text', text: 'engine row' })).toBeDefined()
      expect(await bash.find({ type: 'Text', text: /🦓/ })).toBe(undefined)
      await bash.unmount()
    }
  })

  test('a case_update row in Chinese counts identifiers and shows none', async ($, on) => {
    noState(on, 'board', { language: 'zh', title: 't' })
    const ui = await $.ui.mount({
      plugin: 'zebra-mod', surface: 'terminal', component: 'ToolUse',
      props: toolRow(`${P}case_update`, { identifiers: ['王小雨', 'MZ20240518007'] }),
    })
    expect(await ui.find({ type: 'Text', text: /登记受保护身份信息 2 项（内容不显示）/ })).toBeDefined()
    expect(await ui.find({ type: 'Text', text: /王小雨/ })).toBe(undefined)
    await ui.unmount()
  })

  test('a zebra result is one summary line in place of the JSON', async ($, on) => {
    noState(on, 'board', null)
    noState(on, 'lang', 'en')
    for (const surface of SURFACES) {
      const ui = await $.ui.mount({
        plugin: 'zebra-mod', surface, component: 'ToolResult',
        props: { tool_use_id: 'tu1', tool: `${P}variant_card`, output: ENVELOPE, isErrored: false },
      })
      expect(await ui.find({ type: 'Text', text: /⎿ 4 evidence rows \(E39–E53\)/ })).toBeDefined()
      expect(await ui.find({ type: 'Text', text: /"ledger"/ })).toBe(undefined)
      await ui.unmount()
    }
  })

  test('the spinner names the databases while a query runs, and is left alone otherwise', async ($, on) => {
    noState(on, 'board', null)
    noState(on, 'lang', 'en')
    let runs: unknown[] = [{ id: 'a', tool: 'variant_card', queryZh: '', queryEn: '', sources: ['Ensembl VEP', 'gnomAD', 'ClinVar', 'LitVar2'], startedAt: 0 }]
    on('state.get', { plugin: 'zebra-mod', key: 'running' }, () => ({ value: { value: runs, version: runs.length } }))
    let seen = ''
    engineDraws(on, 'Spinner', 'spinner', e => { seen = e.props.word })
    const props = { word: 'Sauteing', message: null, suffix: '…', mode: 'tool-use' as const }
    const ui = await $.ui.mount({ plugin: 'zebra-mod', surface: 'terminal', component: 'Spinner', props })
    expect(seen).toBe('🦓 Querying Ensembl VEP, gnomAD, ClinVar +1 more')
    await ui.unmount()
    runs = []
    const idle = await $.ui.mount({ plugin: 'zebra-mod', surface: 'terminal', component: 'Spinner', props })
    expect(seen).toBe('Sauteing')
    await idle.unmount()
  })

  test('the footer carries 🦓 and the case, or a shield after a refusal', async ($, on) => {
    noState(on, 'board', { language: 'zh', title: 'Lily — seizures since 6 months, long title' })
    let at: number | null = null
    on('state.get', { plugin: 'zebra-mod', key: 'blockedAt' }, () => ({ value: { value: at, version: at ?? 0 } }))
    let modes: readonly string[] = []
    engineDraws(on, 'SessionMode', 'modes', e => { modes = e.props.modes })
    const a = await $.ui.mount({ plugin: 'zebra-mod', surface: 'terminal', component: 'SessionMode', props: { modes: ['focus'] } })
    expect(modes).toEqual(['focus', '🦓 Lily — seizures sin… · HPO 0 · E0'])
    await a.unmount()
    at = 1
    const b = await $.ui.mount({ plugin: 'zebra-mod', surface: 'terminal', component: 'SessionMode', props: { modes: [] } })
    expect(modes[0]?.startsWith('🛡 ')).toBe(true)
    await b.unmount()
  })

  test('the band gallops in the terminal, speaks in words elsewhere, and is empty when idle', async ($, on) => {
    noState(on, 'board', null)
    noState(on, 'lang', 'zh')
    let runs: unknown[] = [{ id: 'a', tool: 'variant_card', queryZh: 'NM_001165963.4:c.2134C>T', queryEn: '', sources: ['Ensembl VEP', 'gnomAD'], startedAt: 0 }]
    on('state.get', { plugin: 'zebra-mod', key: 'running' }, () => ({ value: { value: runs, version: runs.length } }))
    noState(on, 'blockedAt', null)
    engineDraws(on, 'AbovePrompt', 'engine band')
    const band = (surface: 'terminal' | 'desktop') => ({
      plugin: 'zebra-mod', surface, component: 'AbovePrompt' as const,
      props: { hasSurvey: false, isWorking: true, maxRows: 10, bodyColumns: 120, scroll: { offset: 0, bodyRows: 9 } as any, view: {} as any },
    })
    const t = await $.ui.mount(band('terminal'))
    expect(await t.find({ type: 'Raster', key: 'zebra' })).toBeDefined()
    expect(await t.find({ type: 'Text', text: '🦓 zebra-mod 正在查询公开数据库' })).toBeDefined()
    expect(await t.find({ type: 'Text', text: '变异卡：NM_001165963.4:c.2134C>T' })).toBeDefined()
    await t.unmount()
    const d = await $.ui.mount(band('desktop'))
    expect(await d.find({ type: 'Raster' })).toBe(undefined)
    expect(await d.find({ type: 'Text', text: '数据源：Ensembl VEP · gnomAD' })).toBeDefined()
    await d.unmount()
    runs = []
    const idle = await $.ui.mount(band('terminal'))
    expect(await idle.find({ type: 'Text', text: 'engine band' })).toBeDefined()
    expect(await idle.find({ type: 'Text', text: /zebra-mod/ })).toBe(undefined)
    await idle.unmount()
  })
})

describe('interface option', () => {
  test('off leaves every drawing to the engine', { options: { interface: 'off' } }, async ($, on) => {
    noState(on, 'board', null)
    engineDraws(on, 'ToolUse', 'engine row')
    const ui = await $.ui.mount({
      plugin: 'zebra-mod', surface: 'terminal', component: 'ToolUse',
      props: toolRow(`${P}variant_card`, { variant: 'NM_001165963.4:c.2134C>T' }),
    })
    expect(await ui.find({ type: 'Text', text: /🦓/ })).toBe(undefined)
    await ui.unmount()
  })
})

/** The plugin's state, kept in memory beneath it, as a session's host keeps it. */
function memoryState(on: any) {
  const mem = new Map<string, { value: unknown; version: number }>()
  on('state.get', ($$: unknown, e: { key: string }) => ({ value: mem.get(e.key) ?? { value: undefined, version: 0 } }))
  on('state.set', ($$: unknown, e: { key: string; value: unknown; ifVersion?: number }) => {
    const held = mem.get(e.key)
    const version = (held?.version ?? 0) + 1
    if (e.ifVersion !== undefined && e.ifVersion !== (held?.version ?? 0)) return { value: { isSet: false, version: held?.version ?? 0 } }
    mem.set(e.key, { value: e.value, version })
    return { value: { isSet: true, version } }
  })
  return mem
}

describe('the call is watched', () => {
  test('a turn that used zebra tools ends with one line saying what they did', async ($, on) => {
    const mem = memoryState(on)
    on('tool.check', () => ({ decision: 'allow' as const }))
    on('env.get', () => ({ value: '' }))
    on('store.get', () => ({ value: true }))
    on('process.run', () => ({ value: { exitCode: 0, stdout: JSON.stringify({ ok: true, ...JSON.parse(ENVELOPE) }), stderr: '', isStdoutTruncated: false, isStderrTruncated: false } }))
    const appended: string[] = []
    on('session.append', ($$: unknown, e: { message: { content: Array<{ text?: string }> } }, next: (e: unknown) => unknown) => {
      appended.push(e.message.content.map(b => b.text ?? '').join(''))
      return next(e)
    })
    on('turn.start', ($$: unknown, e: { turnId: string }) => ({ turnId: e.turnId }))
    on('turn.complete', () => ({ text: '' }))

    await $.turn.start({ text: 'what does SCN1A c.2134C>T mean?', turnId: 't1' } as never)
    const r = await $.tool.call({ tool: `${P}variant_card`, variant: 'NM_001165963.4:c.2134C>T' } as never)
    expect('result' in r).toBe(true)
    expect((mem.get('running')?.value as unknown[] | undefined) ?? []).toHaveLength(0)
    expect((mem.get('turn')?.value as { calls?: number } | undefined)?.calls).toBe(1)
    await $.turn.complete({ answer: '', durationMs: 1200, isAborted: false, turnId: 't1', reason: 'answer' } as never)
    const line = appended.find(t => t.startsWith('🦓 zebra-mod this turn')) ?? ''
    expect(line.includes('1 call')).toBe(true)
    expect(line.includes('4 evidence rows (E39–E53)')).toBe(true)
    expect((mem.get('turn')?.value as { calls?: number } | undefined)?.calls).toBe(0)
  })

  test('a turn without zebra tools leaves no line', async ($, on) => {
    memoryState(on)
    on('store.get', () => ({ value: true }))
    on('env.get', () => ({ value: '' }))
    const appended: string[] = []
    on('session.append', ($$: unknown, e: { message: unknown }, next: (e: unknown) => unknown) => {
      appended.push(JSON.stringify(e.message))
      return next(e)
    })
    on('turn.start', ($$: unknown, e: { turnId: string }) => ({ turnId: e.turnId }))
    on('turn.complete', () => ({ text: '' }))
    await $.turn.start({ text: 'reverse a linked list', turnId: 't2' } as never)
    await $.turn.complete({ answer: '', durationMs: 5, isAborted: false, turnId: 't2', reason: 'answer' } as never)
    expect(appended).toHaveLength(0)
  })
})
