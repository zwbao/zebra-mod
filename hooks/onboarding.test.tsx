import { describe, expect, test } from 'claude-code/testing'

import { welcomeLines } from './ui'

const DEMO = {
  path: '/Users/test/zebra-cases/demo-xiaoyu',
  reused: false,
  language: 'zh',
  title: '示例：小雨 — 发热抽搐（合成数据）',
  records: ['基因检测报告.txt', '门诊病历.txt'],
  identifiers: 6,
  first_prompt: '我把女儿的门诊病历和基因检测报告放在这个病例的 records 文件夹里了。……下一步该问医生什么？',
  next_prompts: ['Dravet 综合征现在有哪些获批的药和正在招募的临床试验？我们在中国。'],
}
const CASE_JSON = JSON.stringify({ schema: 'zebra.case/1', title: DEMO.title, language: 'zh', privacy: { identifiers: ['王小雨'] } })
const RUN = (stdout: string) => ({ value: { exitCode: 0, stdout, stderr: '', isStdoutTruncated: false, isStderrTruncated: false } })
const COMMAND = { command: 'zebra', origin: { kind: 'composer' as const }, presentation: { isFullscreen: false, columns: 120 } }

/** The plugin's state, kept in memory beneath it, as a session's host keeps it. */
function memoryState(on: any) {
  const mem = new Map<string, { value: unknown; version: number }>()
  on('state.get', ($$: unknown, e: { key: string }) => ({ value: mem.get(e.key) ?? { value: undefined, version: 0 } }))
  on('state.set', ($$: unknown, e: { key: string; value: unknown; ifVersion?: number }) => {
    const held = mem.get(e.key)
    if (e.ifVersion !== undefined && e.ifVersion !== (held?.version ?? 0)) return { value: { isSet: false, version: held?.version ?? 0 } }
    const version = (held?.version ?? 0) + 1
    mem.set(e.key, { value: e.value, version })
    return { value: { isSet: true, version } }
  })
  return mem
}

const BOARD = { path: DEMO.path, id: 'demo-xiaoyu', title: DEMO.title, role: 'family', language: 'zh', updated_at: '', phenotypes: [], variants: [], hypotheses: [], therapy_leads: [], questions: [], evidence_count: 0, identifiers: 6 }
const envelope = (result: unknown) => JSON.stringify({ ok: true, result, sources: [], warnings: [], ledger: [] })

/** A host with LANG, a store, a filesystem holding the demo case, and the CLI: `case demo` answers `stdout`, `case summary` the board. */
function host(on: any, o: { lang?: string; stdout?: string; store?: Map<string, unknown>; argv?: string[][] } = {}) {
  const store = o.store ?? new Map<string, unknown>()
  on('env.get', ($$: unknown, e: { name?: string; key?: string }) => ({ value: (e.name ?? e.key) === 'LANG' ? (o.lang ?? 'zh_CN.UTF-8') : '' }))
  on('env.set', () => ({ value: undefined }))
  on('store.get', ($$: unknown, e: { key: string }) => ({ value: store.get(e.key) }))
  on('store.set', ($$: unknown, e: { key: string; value: unknown }) => {
    store.set(e.key, e.value)
    return { value: undefined }
  })
  on('fs.read', ($$: unknown, e: { path: string }) => (e.path.endsWith('case.json') ? { value: CASE_JSON } : { deny: 'ENOENT' }))
  on('fs.exists', () => ({ value: true }))
  on('fs.stat', ($$: unknown, e: { path: string }) => ({ value: { kind: 'file', size: 1, mtimeMs: 1, isLink: false, realPath: e.path } }))
  on('session.cwd', () => ({ value: '/Users/test' }))
  on('process.run', ($$: unknown, e: { argv: string[] }) => {
    o.argv?.push(e.argv)
    if (e.argv.includes('summary')) return RUN(envelope(BOARD))
    if (e.argv.includes('demo')) return RUN(o.stdout ?? envelope(DEMO))
    return RUN(envelope({}))
  })
  on('ui.open', () => ({ value: { isPlaced: true } }))
  on('ui.status', () => ({ value: undefined }))
  return store
}

describe('/zebra demo', () => {
  test('opens the demo case, makes it active and puts the first question in the prompt box', async ($, on) => {
    const mem = memoryState(on)
    const argv: string[][] = []
    host(on, { argv })
    let filled = ''
    on('prompt.fill', ($$: unknown, e: { text: string }) => {
      filled = e.text
      return { isFilled: true }
    })
    const r = await $.command.run({ ...COMMAND, args: 'demo' } as never)
    const text = String((r as { text?: string }).text)
    expect(argv.some(a => a.includes('case') && a.includes('demo') && a.includes('zh'))).toBe(true)
    expect(filled).toBe(DEMO.first_prompt)
    expect(mem.get('casePath')?.value).toBe(DEMO.path)
    expect(text.includes('已打开示例病例')).toBe(true)
    expect(text.includes('第一个问题已经放进输入框，按回车就开始')).toBe(true)
    expect(text.includes('Dravet 综合征')).toBe(true)
  })

  test('says how to begin when there is no prompt box to fill, and speaks English when asked', async ($, on) => {
    memoryState(on)
    host(on, { lang: 'en_US.UTF-8', stdout: envelope({ ...DEMO, language: 'en', title: 'Demo: Lily' }) })
    on('prompt.fill', () => ({ isFilled: false }))
    const r = await $.command.run({ ...COMMAND, args: 'demo en' } as never)
    const text = String((r as { text?: string }).text)
    expect(text.includes('Demo case open')).toBe(true)
    expect(text.includes('Send Claude this to begin')).toBe(true)
  })

  test('a CLI refusal is said plainly', async ($, on) => {
    memoryState(on)
    host(on, { stdout: JSON.stringify({ ok: false, error: { type: 'UsageError', message: '/x holds a case that is not the demo; give another folder' } }) })
    on('prompt.fill', () => ({ isFilled: true }))
    const r = await $.command.run({ ...COMMAND, args: 'demo /x' } as never)
    expect(String((r as { text?: string }).text).startsWith('示例病例没能准备好：/x holds a case that is not the demo')).toBe(true)
  })
})

describe('/zebra is the guide', () => {
  const guide = async ($: any, on: any, args: string, lang: string): Promise<string> => {
    memoryState(on)
    host(on, { lang })
    const r = await $.command.run({ ...COMMAND, args })
    return String((r as { text?: string }).text)
  }
  test('in Chinese: just ask, the demo, what it does, the commands', async ($, on) => {
    const text = await guide($, on, '', 'zh_CN.UTF-8')
    for (const w of ['直接提问', '/zebra demo', '能做什么', '/zebra new <目录> [标题]', '隐私']) expect(text.includes(w)).toBe(true)
  })
  test('in English from LANG', async ($, on) => {
    const text = await guide($, on, '', 'en_US.UTF-8')
    for (const w of ['just ask', '/zebra demo', 'What it does', '/zebra new <dir> [title]']) expect(text.includes(w)).toBe(true)
  })
  test('in the language asked for', async ($, on) => {
    expect((await guide($, on, 'en', 'zh_CN.UTF-8')).includes('just ask')).toBe(true)
  })
})

describe('onboarding', () => {
  test('the first interactive session shows a toast at once and a card until the first prompt, then never again', async ($, on) => {
    const mem = memoryState(on)
    const store = host(on)
    const toasts: string[] = []
    on('ui.toast', ($$: unknown, e: { text: string }) => {
      toasts.push(e.text)
      return { value: undefined }
    })
    on('session.start', ($$: unknown, e: { cwd: string }) => ({ cwd: e.cwd }))
    on('turn.start', ($$: unknown, e: { turnId: string }) => ({ turnId: e.turnId }))
    on('tool.register', () => ({ value: undefined }))
    on('command.register', () => ({ value: undefined }))

    await $.session.start({ cwd: '/Users/test', surface: 'terminal', isInteractive: true } as never)
    expect(toasts.some(t => t.includes('/zebra demo'))).toBe(true)
    expect(mem.get('onboarding')?.value).toBe(true)

    await $.turn.start({ text: '你好', turnId: 't1' } as never)
    expect(mem.get('onboarding')?.value).toBe(false)
    expect(store.get('onboarded')).toBe(true)

    toasts.length = 0
    await $.session.start({ cwd: '/Users/test', surface: 'terminal', isInteractive: true } as never)
    expect(toasts.some(t => t.includes('/zebra demo'))).toBe(false)
  })

  test('the card shows above the prompt until a case is open', async ($, on) => {
    const mem = memoryState(on)
    host(on)
    mem.set('onboarding', { value: true, version: 1 })
    on('ui.render', { component: 'AbovePrompt' }, ($$: any, e: any) => {
      const { Text } = $$.ui.resolve(e)
      return <Text>engine band</Text>
    })
    const band = { plugin: 'zebra-mod', surface: 'desktop' as const, component: 'AbovePrompt' as const,
      props: { hasSurvey: false, isWorking: false, maxRows: 10, bodyColumns: 120, scroll: { offset: 0, bodyRows: 9 } as any, view: {} as any } }
    const a = await $.ui.mount(band)
    expect(await a.find({ type: 'Text', text: /zebra-mod 已就绪/ })).toBeDefined()
    await a.unmount()
    mem.set('board', { value: { title: 'demo', language: 'zh' }, version: 1 })
    const b = await $.ui.mount(band)
    expect(await b.find({ type: 'Text', text: 'engine band' })).toBeDefined()
    await b.unmount()
  })

  test('the welcome card says to just ask and names the demo', async () => {
    expect(welcomeLines(true)[1].startsWith('直接用中文描述症状')).toBe(true)
    expect(welcomeLines(true)[2].includes('/zebra demo')).toBe(true)
    expect(welcomeLines(false)[2].includes('/zebra demo')).toBe(true)
  })
})
