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
    for (const core of ['case_status', 'case_update', 'hpo_search', 'phenotype_rank', 'gene_card', 'variant_card', 'disease_card',
      'acmg', 'cnv_interpret', 's2f_predict', 'therapy_landscape', 'trials_search', 'literature_search', 'rare_stats',
      'edit_check', 'china_rare', 'report_export']) {
      expect(names.includes(core)).toBe(true)
    }
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

  test('CP1-13 report_export resolves a case-relative report and passes formats', async () => {
    const argv = await toolArgv(def('report_export'), { report: 'reports/family-letter.md', formats: ['docx', 'pdf'] }, '/cases/x/')
    expect(argv).toEqual(['report', 'export', '/cases/x/reports/family-letter.md', '--to', 'docx,pdf'])
    const abs = await toolArgv(def('report_export'), { report: '/tmp/r.md' }, '/cases/x')
    expect(abs).toEqual(['report', 'export', '/tmp/r.md'])
  })

  test('W4-W7 new and widened tools build the CLI argv the parser accepts (tests/test_cli_tools.py parses the same lists)', async () => {
    expect(await toolArgv(def('access'), { query: '脊髓性肌萎缩症', province: '浙江', status: 'ANY' }, null))
      .toEqual(['access', '脊髓性肌萎缩症', '--province', '浙江', '--status', 'ANY'])
    expect(await toolArgv(def('expression'), { gene: 'CFTR', top: 5 }, null)).toEqual(['expression', 'CFTR', '--top', '5'])
    expect(await toolArgv(def('cnv_interpret'), { result: 'SMN1 exon 7 deletion', copies: 0, smn2_copies: 3, sex: 'female' }, null))
      .toEqual(['cnv', 'SMN1 exon 7 deletion', '--copies', '0', '--sex', 'female', '--smn2-copies', '3'])
    expect(await toolArgv(def('phenotype_rank'), { present: ['HP:0001250'], excluded_weight: 1, local_method: 'lr' }, null))
      .toEqual(['phenotype', 'rank', '--present', 'HP:0001250', '--local-method', 'lr', '--excluded-weight', '1'])
    expect(await toolArgv(def('case_recheck'), { plan: true }, '/cases/x')).toEqual(['case', 'recheck', '/cases/x', '--plan'])
    expect(await toolArgv(def('aso_screen'), { variant: 'NM_000492.4:c.3718-2477C>T', event: 'pseudoexon', uniqueness: true }, null))
      .toEqual(['aso', 'NM_000492.4:c.3718-2477C>T', '--event', 'pseudoexon', '--uniqueness'])
    expect(await toolArgv(def('china_rare'), { query: 'hospitals', province: '浙江省' }, null)).toEqual(['china', 'hospitals', '--province', '浙江省'])
    expect(await toolArgv(def('variant_card'), { variant: 'm.3243A>G', heteroplasmy: 35 }, null)).toEqual(['variant', 'm.3243A>G', '--heteroplasmy', '35'])
    expect(await toolArgv(def('acmg'), { mode: 'suggest', variant: '7-117559590-ATCT-A', inheritance: 'AR', prevalence: 0.0004, allelic: 0.9, no_splice_lookup: true }, null))
      .toEqual(['acmg', 'suggest', '7-117559590-ATCT-A', '--inheritance', 'AR', '--prevalence', '0.0004', '--allelic', '0.9', '--no-splice-lookup'])
    expect(await toolArgv(def('trials_search'), { condition: 'SMA', country: 'China', keep_unrelated: true }, null))
      .toEqual(['trials', 'SMA', '--country', 'China', '--keep-unrelated'])
  })
})

