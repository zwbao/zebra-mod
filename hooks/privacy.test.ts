import { describe, expect, test } from 'claude-code/testing'

import { guardInput, isOutboundShell, uploadsGenome } from './privacy'

describe('privacy gate', () => {
  test('passes ordinary research queries', () => {
    expect(guardInput({ url: 'https://rest.ensembl.org/vep/human/hgvs/NM_000492.4:c.1521_1523del' }, [])).toBe(undefined)
    expect(guardInput({ variant: '7-117559590-ATCT-A', assembly: 'GRCh38' }, [])).toBe(undefined)
    expect(guardInput({ query: 'chr1:155235002 PMID 31189966 NCT04442295' }, [])).toBe(undefined)
  })

  test('catches the case identifiers, whatever the case of the letters', () => {
    expect(guardInput({ query: 'Dravet syndrome Zhang San' }, ['zhang san'])).toBe('protected identifier #1 from the case')
    expect(guardInput({ query: '张小明 癫痫' }, ['李四', '张小明'])).toBe('protected identifier #2 from the case')
  })

  test('catches ID numbers, phones and emails', () => {
    expect(guardInput({ q: '330106201903021234' }, [])).toBe('what looks like a Chinese resident ID number')
    expect(guardInput({ q: 'call 13812345678' }, [])).toBe('what looks like a mobile phone number')
    expect(guardInput({ q: 'mom@example.com' }, [])).toBe('an email address')
  })

  test('asks before genome files leave the machine', () => {
    expect(uploadsGenome('curl -F file=@proband.vcf.gz https://example.org/upload')).toBe(true)
    expect(uploadsGenome('scp trio.bam me@server:/data/')).toBe(true)
    expect(uploadsGenome('aws s3 cp sample.cram s3://bucket/')).toBe(true)
    expect(uploadsGenome('bcftools view proband.vcf.gz | head')).toBe(false)
    expect(uploadsGenome('curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz')).toBe(false)
  })

  test('knows which shell commands reach the network', () => {
    expect(isOutboundShell('curl -s https://example.org')).toBe(true)
    expect(isOutboundShell('zebra lit "Dravet syndrome"')).toBe(true)
    expect(isOutboundShell('~/zebra-mod/bin/zebra --json variant 2-166042334-G-A')).toBe(true)
    expect(isOutboundShell('zebra --case ~/c case identifiers --add "张三"')).toBe(false)
    expect(isOutboundShell('ls records/ && cat records/a.txt')).toBe(false)
    expect(isOutboundShell('ls ~/zebra-mod/skills')).toBe(false)
  })
})
