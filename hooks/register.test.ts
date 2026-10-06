import { describe, expect, test } from 'claude-code/testing'

const FACTS = { model: 'claude-test', promptModel: 'claude-test', surfaces: [], tools: [], outputStyle: null, traits: [] } as const
const ENVELOPE = JSON.stringify({ ok: true, result: {}, sources: [], warnings: [], ledger: [] })
const CASE = '/tmp/zebra-test-case'
const CASE_JSON = JSON.stringify({ schema: 'zebra.case/1', title: 't', privacy: { identifiers: ['王小雨', 'Zhang Xiaoming', 'MZ0012345'] } })

const run = (stdout: string) => ({ value: { exitCode: 0, stdout, stderr: '', isStdoutTruncated: false, isStderrTruncated: false } })

/**
 * The world beneath the plugin with a case open: its path in state, its case.json on disk
 * (null: unreadable), other files by path, and a CLI that answers `stdout`.
 */
function openCase(on: any, o: { caseJson?: string | null; files?: Record<string, string>; inside?: (p: string) => boolean; stdout?: string } = {}) {
  on('state.get', { plugin: 'zebra-mod', key: 'casePath' }, () => ({ value: { value: CASE, version: 1 } }))
  on('fs.read', (_$: unknown, e: { path: string }) => {
    if (e.path === `${CASE}/case.json`) {
      if (o.caseJson === null) return { deny: 'EACCES' }
      return { value: o.caseJson ?? CASE_JSON }
    }
    const f = o.files?.[e.path]
    return f === undefined ? { deny: `ENOENT ${e.path}` } : { value: f }
  })
  on('fs.stat', (_$: unknown, e: { path: string }) => {
    const realPath = o.inside?.(e.path) ? `${CASE}/${e.path.replace(/^.*?records\//, 'records/')}` : e.path
    return { value: { kind: 'file', size: 1, mtimeMs: 1, isLink: false, realPath } }
  })
  on('env.get', () => ({ value: '/Users/test' }))
  on('process.run', () => run(o.stdout ?? ENVELOPE))
}

// The mod's own tools are answered by its tool.call hook, beneath which nothing runs the engine's
// permission chain. These tests prove the hook asks for the decision itself (A-P0-1). What happens
// with a case open (identifiers from case.json) is verified in a real headless session, because the
// test kit cannot hold the plugin's case state or run its CLI.
describe('A-P0-1 the mod’s own tools go through the permission chain and the gate', () => {
  test('a deny beneath stops an own tool before it runs', async ($, on) => {
    let ran = false
    on('tool.check', () => ({ decision: 'deny' as const, reason: 'the user said no' }))
    on('process.run', () => {
      ran = true
      return run(ENVELOPE)
    })
    const r = await $.tool.call({ tool: 'mcp__zebra-mod__hpo_search', text: 'seizure' } as never)
    expect('deny' in r && r.deny).toBe('the user said no')
    expect(ran).toBe(false)
  })

  test('the gate stops an own tool carrying a phone number', async ($, on) => {
    on('tool.check', () => ({ decision: 'allow' as const }))
    const r = await $.tool.call({ tool: 'mcp__zebra-mod__literature_search', query: 'call 13812345678' } as never)
    expect('deny' in r && String(r.deny).includes('mobile phone')).toBe(true)
  })

  test('an ask with nobody to answer it is a no', async ($, on) => {
    on('tool.check', () => ({ decision: 'ask' as const }))
    const r = await $.tool.call({ tool: 'mcp__zebra-mod__case_update', questions: ['x'] } as never)
    expect('deny' in r && String(r.deny).includes('not approved')).toBe(true)
  })
})

describe('A-P1-6 with a case open', () => {
  test('A-P0-1 an own tool carrying a registered name is refused before it runs', async ($, on) => {
    openCase(on)
    on('tool.check', () => ({ decision: 'allow' as const }))
    const r = await $.tool.call({ tool: 'mcp__zebra-mod__literature_search', query: '王小雨 Dravet' } as never)
    expect('deny' in r && String(r.deny).includes('protected identifier')).toBe(true)
    const pinyin = await $.tool.call({ tool: 'mcp__zebra-mod__literature_search', query: 'xiaoming zhang SCN1A' } as never)
    expect('deny' in pinyin).toBe(true)
  })

  test('an allowed lookup runs, and each source carries its evidence id', async ($, on) => {
    openCase(on, { stdout: JSON.stringify({ ok: true, result: { hits: 1 }, sources: [{ db: 'HPO' }, { db: 'Monarch' }], warnings: [], ledger: ['E7', 'E8'] }) })
    on('tool.check', () => ({ decision: 'allow' as const }))
    const r = await $.tool.call({ tool: 'mcp__zebra-mod__hpo_search', text: 'seizure' } as never)
    const body = JSON.parse(String((r as { result?: unknown }).result)) as { sources: Array<{ eid?: string; db: string }> }
    expect(body.sources.map(x => `${x.eid}:${x.db}`)).toEqual(['E7:HPO', 'E8:Monarch'])
  })

  test('F3 a case file that cannot be read closes the gate', async ($, on) => {
    openCase(on, { caseJson: null })
    on('tool.check', () => ({ decision: 'allow' as const }))
    const r = await $.tool.check({ tool: 'WebSearch', input: { query: 'Dravet syndrome' } })
    expect(r.decision).toBe('deny')
  })

  test('A-P1-3 a sample name or a case path that stays local does not trip the gate', async ($, on) => {
    openCase(on)
    on('tool.check', () => ({ decision: 'allow' as const }))
    const triage = await $.tool.check({ tool: 'Bash', input: { command: `zebra --case ${CASE} vcf triage x.vcf.gz --proband MZ0012345 --mother M --father F` } })
    expect(triage.decision).toBe('allow')
    const ledger = await $.tool.check({ tool: 'Bash', input: { command: 'zebra --case /Users/test/cases/zhang-xiaoming case ledger' } })
    expect(ledger.decision).toBe('allow')
    const sent = await $.tool.check({ tool: 'Bash', input: { command: 'zebra lit "MZ0012345 Dravet"' } })
    expect(sent.decision).toBe('deny')
  })

  test('A-P1-2 a case file sent by scp after a cd is asked about', async ($, on) => {
    openCase(on, { inside: p => p.includes('records/') })
    on('tool.check', () => ({ decision: 'allow' as const }))
    const r = await $.tool.check({ tool: 'Bash', input: { command: `cd ${CASE} && scp records/a.pdf me@host:/tmp/` } })
    expect(r.decision).toBe('ask')
  })

  test('A-P1-1 a page that cannot be read is not published unseen', async ($, on) => {
    openCase(on)
    on('tool.check', () => ({ decision: 'allow' as const }))
    const r = await $.tool.check({ tool: 'Artifact', input: { file_path: '/tmp/out/report.html' } })
    expect(r.decision).toBe('ask')
  })

  test('A-P1-1 a published page is scanned for the case’s names and for ID numbers', async ($, on) => {
    openCase(on, { files: { '/tmp/out/a.html': '<p>王小雨 的随访</p>', '/tmp/out/b.html': '<p>ID 330106202104120023</p>', '/tmp/out/c.html': '<p>SCN1A Dravet</p>' } })
    on('tool.check', () => ({ decision: 'allow' as const }))
    expect((await $.tool.check({ tool: 'Artifact', input: { file_path: '/tmp/out/a.html' } })).decision).toBe('deny')
    expect((await $.tool.check({ tool: 'Artifact', input: { files: { 'b.html': '/tmp/out/b.html' } } })).decision).toBe('deny')
    expect((await $.tool.check({ tool: 'Artifact', input: { file_path: '/tmp/out/c.html' } })).decision).toBe('allow')
  })

  test('A-P2-5 cnv_interpret that records into the case is not waved through as read-only', async ($, on) => {
    openCase(on)
    on('tool.check', () => ({ decision: 'ask' as const }))
    expect((await $.tool.check({ tool: 'mcp__zebra-mod__cnv_interpret', input: { query: 'DMD exon 45-50 deletion' } })).decision).toBe('allow')
    expect((await $.tool.check({ tool: 'mcp__zebra-mod__cnv_interpret', input: { query: 'DMD exon 45-50 deletion', record: true } })).decision).toBe('ask')
  })
})

describe('zebra-mod permissions', () => {
  test('F1 a read-only lookup is spared a prompt', async ($, on) => {
    on('tool.check', () => ({ decision: 'ask' as const }))
    const lookup = await $.tool.check({ tool: 'mcp__zebra-mod__hpo_search', input: { text: 'seizure' } })
    expect(lookup.decision).toBe('allow')
  })

  test('F1 a deny beneath is never upgraded', async ($, on) => {
    on('tool.check', () => ({ decision: 'deny' as const, reason: 'the user said no' }))
    const denied = await $.tool.check({ tool: 'mcp__zebra-mod__hpo_search', input: { text: 'seizure' } })
    expect(denied.decision).toBe('deny')
    expect(denied.reason).toBe('the user said no')
  })

  test('F1 a settings rule beneath stands, and a write is left to the mode', async ($, on) => {
    on('tool.check', { tool: 'mcp__zebra-mod__hpo_search' }, () => ({ decision: 'ask' as const, rule: 'mcp__zebra-mod__hpo_search' }))
    on('tool.check', { tool: 'mcp__zebra-mod__case_update' }, () => ({ decision: 'ask' as const }))
    expect((await $.tool.check({ tool: 'mcp__zebra-mod__hpo_search', input: { text: 'x' } })).decision).toBe('ask')
    expect((await $.tool.check({ tool: 'mcp__zebra-mod__case_update', input: { questions: ['x'] } })).decision).toBe('ask')
  })

  test('F4 a local case write is not treated as leaving the machine', async ($, on) => {
    on('tool.check', () => ({ decision: 'allow' as const }))
    const write = await $.tool.check({ tool: 'mcp__zebra-mod__case_update', input: { questions: ['ask about 13812345678'] } })
    expect(write.decision).toBe('allow')
  })

  test('the gate stops a phone number and asks before a genome upload', async ($, on) => {
    on('tool.check', () => ({ decision: 'allow' as const }))
    const leak = await $.tool.check({ tool: 'WebFetch', input: { url: 'https://example.org/?phone=13812345678', prompt: 'x' } })
    expect(leak.decision).toBe('deny')
    const upload = await $.tool.check({ tool: 'Bash', input: { command: 'curl -F f=@proband.vcf.gz https://example.org/up' } })
    expect(upload.decision).toBe('ask')
    const local = await $.tool.check({ tool: 'Bash', input: { command: 'ls records/' } })
    expect(local.decision).toBe('allow')
  })

  test('E6 Artifact text and a remote agent are treated as leaving the machine', async ($, on) => {
    on('tool.check', () => ({ decision: 'allow' as const }))
    const published = await $.tool.check({ tool: 'Artifact', input: { title: 'call 13812345678' } })
    expect(published.decision).toBe('deny')
    const remote = await $.tool.check({ tool: 'Agent', input: { prompt: 'patient id 330106201903021234', isolation: 'remote' } })
    expect(remote.decision).toBe('deny')
  })

  test('the gate can be switched off', { options: { privacyGate: false } }, async ($, on) => {
    on('tool.check', () => ({ decision: 'allow' as const }))
    const leak = await $.tool.check({ tool: 'WebFetch', input: { url: 'https://example.org/?phone=13812345678', prompt: 'x' } })
    expect(leak.decision).toBe('allow')
  })
})

describe('zebra-mod doctrine', () => {
  test('it is added as a session-scope section, after the engine’s own', async ($, on) => {
    on('prompt.compose', () => ({ sections: [{ id: 'intro', text: 'You are Claude Code.', scope: 'shared' as const }] }))
    const { sections } = await $.prompt.compose(FACTS)
    const doctrine = sections.find(s => s.id === 'zebra-mod:doctrine')
    expect(doctrine?.scope).toBe('session')
    expect(doctrine?.text.includes('Evidence, not recall')).toBe(true)
    expect(sections[0]?.id).toBe('intro')
  })

  test('with no case it says so instead of naming one', async ($, on) => {
    on('prompt.compose', () => ({ sections: [] }))
    const { sections } = await $.prompt.compose(FACTS)
    expect(sections[0]?.text.includes('No active case')).toBe(true)
  })

  test('off leaves the prompt alone', { options: { doctrine: 'off' } }, async ($, on) => {
    on('prompt.compose', () => ({ sections: [{ id: 'intro', text: 'You are Claude Code.', scope: 'shared' as const }] }))
    const { sections } = await $.prompt.compose(FACTS)
    expect(sections.length).toBe(1)
  })

  test('case mode stays quiet until a case is open', { options: { doctrine: 'case' } }, async ($, on) => {
    on('prompt.compose', () => ({ sections: [{ id: 'intro', text: 'x', scope: 'shared' as const }] }))
    const { sections } = await $.prompt.compose(FACTS)
    expect(sections.length).toBe(1)
  })
})
