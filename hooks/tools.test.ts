import { describe, expect, test } from 'claude-code/testing'

import { TOOLS, toolArgv } from './tools'

const def = (name: string) => {
  const d = TOOLS.find(t => t.name === name)
  if (!d) throw new Error(`no tool ${name}`)
  return d
}

describe('tool → CLI argv', () => {
  test('every tool has a unique, legal name and an object schema', () => {
    const names = TOOLS.map(t => t.name)
    expect(names.length).toBe(15)
    expect(new Set(names).size).toBe(names.length)
    for (const t of TOOLS) {
      expect(/^[A-Za-z0-9_-]{1,64}$/.test(t.name)).toBe(true)
      expect(t.inputSchema.type).toBe('object')
    }
  })

  test('phenotype_rank', async () => {
    const argv = await toolArgv(def('phenotype_rank'), { present: ['HP:0001250', 'HP:0001263'], excluded: ['HP:0001252'], sources: ['local', 'monarch'], top: 10 }, null)
    expect(argv).toEqual(['phenotype', 'rank', '--present', 'HP:0001250', 'HP:0001263', '--exclude', 'HP:0001252', '--sources', 'local,monarch', '--top', '10'])
  })

  test('case tools need an active case', async () => {
    let message = ''
    try {
      await toolArgv(def('case_status'), {}, null)
    } catch (err) {
      message = String(err)
    }
    expect(message.includes('no active case')).toBe(true)
    expect(await toolArgv(def('case_status'), {}, '/x/case')).toEqual(['case', 'summary', '/x/case'])
  })

  test('rare_stats turns params into flags', async () => {
    const argv = await toolArgv(def('rare_stats'), { method: 'maxaf', params: { prevalence: 0.0001, allelic: 0.1, inheritance: 'biallelic' } }, null)
    expect(argv).toEqual(['stats', 'maxaf', '--prevalence', '0.0001', '--allelic', '0.1', '--inheritance', 'biallelic'])
  })

  test('acmg classify and suggest', async () => {
    expect(await toolArgv(def('acmg'), { mode: 'classify', codes: ['PVS1', 'PM2_Supporting'] }, null)).toEqual(['acmg', 'classify', 'PVS1', 'PM2_Supporting'])
    expect(await toolArgv(def('acmg'), { mode: 'suggest', variant: '2-166001-C-T', assembly: 'GRCh38', inheritance: 'AD' }, null)).toEqual(['acmg', 'suggest', '2-166001-C-T', '--assembly', 'GRCh38', '--inheritance', 'AD'])
  })

  test('A-P2-1 only the schema own keys reach the CLI', async () => {
    const argv = await toolArgv(def('case_update'), {
      tool: 'mcp__zebra-mod__case_update', tool_use_id: 'toolu_1', agentId: 'a1',
      questions: ['parental samples?'],
    }, '/x/case')
    expect(argv.slice(0, 3)).toEqual(['case', 'apply', '/x/case'])
    expect(JSON.parse(argv[4] as string)).toEqual({ questions: ['parental samples?'] })
  })

  test('F7 a value that reads as a flag is refused', async () => {
    let message = ''
    try {
      await toolArgv(def('variant_card'), { variant: '--help' }, null)
    } catch (err) {
      message = String(err)
    }
    expect(message.includes('reads as a command-line flag')).toBe(true)
  })

  test('A-P1-7 rare_stats accepts only the flags its method takes', async () => {
    let message = ''
    try {
      await toolArgv(def('rare_stats'), { method: 'maxaf', params: { case: '/other/case', help: true } }, null)
    } catch (err) {
      message = String(err)
    }
    expect(message.includes('does not take')).toBe(true)
    const ok = await toolArgv(def('rare_stats'), { method: 'carrier', params: { prevalence: 0.0004 } }, null)
    expect(ok).toEqual(['stats', 'carrier', '--prevalence', '0.0004'])
  })

  test('an empty required argument is refused', async () => {
    let failed = false
    try {
      await toolArgv(def('gene_card'), { symbol: '  ' }, null)
    } catch {
      failed = true
    }
    expect(failed).toBe(true)
  })
})
