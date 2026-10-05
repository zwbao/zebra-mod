import { describe, expect, test } from 'claude-code/testing'

const FACTS = { model: 'claude-test', promptModel: 'claude-test', surfaces: [], tools: [], outputStyle: null, traits: [] } as const

describe('zebra-mod permissions', () => {
  test('F1 a read-only lookup is spared a prompt, and nothing else is overridden', async ($, on) => {
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
    const write = await $.tool.check({
      tool: 'mcp__zebra-mod__case_update',
      input: { questions: ['ask about 13812345678'] },
    })
    expect(write.decision).toBe('allow')
  })

  test('the gate stops a phone number and a genome upload on the way out', async ($, on) => {
    on('tool.check', () => ({ decision: 'allow' as const }))
    const leak = await $.tool.check({ tool: 'WebFetch', input: { url: 'https://example.org/?phone=13812345678', prompt: 'x' } })
    expect(leak.decision).toBe('deny')
    const upload = await $.tool.check({ tool: 'Bash', input: { command: 'curl -F f=@proband.vcf.gz https://example.org/up' } })
    expect(upload.decision).toBe('ask')
    const local = await $.tool.check({ tool: 'Bash', input: { command: 'ls records/' } })
    expect(local.decision).toBe('allow')
  })

  test('E6 Artifact and a remote agent are treated as leaving the machine', async ($, on) => {
    on('tool.check', () => ({ decision: 'allow' as const }))
    const published = await $.tool.check({ tool: 'Artifact', input: { file_path: '/tmp/x.html', description: 'call 13812345678' } })
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
