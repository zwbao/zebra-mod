// The privacy gate. Three jobs, each best effort and each said plainly where it is not:
//   1. keep a case's registered identifiers out of anything that leaves the machine,
//   2. notice the shell commands and tools that leave the machine at all,
//   3. notice when a file from the case folder is about to be sent somewhere.
// Matching happens on a normalised form, because a name reaches a web service
// percent-encoded, '+'-joined, full-width or re-spaced far more often than verbatim.

const ZERO_WIDTH = /[​-‏‪-‮⁠﻿]/g
const SEPARATORS = /[\s\-_/.,·'"()[\]{}]+/g

/** Percent-decode, '+'-decode, \uXXXX-decode, NFKC-fold, strip zero-width, lowercase. */
function fold(text: string): string {
  let out = text
  for (let i = 0; i < 3; i++) {
    const before = out
    if (/%[0-9a-fA-F]{2}/.test(out)) {
      try {
        out = decodeURIComponent(out.replace(/\+/g, ' '))
      } catch {
        out = out.replace(/%([0-9a-fA-F]{2})/g, (_m, h) => String.fromCharCode(parseInt(h, 16))).replace(/\+/g, ' ')
      }
    }
    if (/\\u[0-9a-fA-F]{4}/.test(out)) {
      out = out.replace(/\\u([0-9a-fA-F]{4})/g, (_m, h) => String.fromCharCode(parseInt(h, 16)))
    }
    // HTML: a page can spell a name as entities (&#24352;) or split it with tags (<b>张</b>小明)
    if (/&#?\w+;/.test(out)) {
      out = out.replace(/&#x([0-9a-f]+);/gi, (_m, h) => String.fromCodePoint(parseInt(h, 16)))
        .replace(/&#(\d+);/g, (_m, d) => String.fromCodePoint(Number(d)))
        .replace(/&(amp|lt|gt|quot|apos|nbsp);/g, (_m, n) => ({ amp: '&', lt: '<', gt: '>', quot: '"', apos: "'", nbsp: ' ' } as Record<string, string>)[n] ?? _m)
    }
    if (/<\/?[a-z][^<>]{0,200}>/i.test(out)) out = out.replace(/<\/?[a-z][^<>]{0,200}>/gi, '')
    if (out === before) break
  }
  // tone marks (Zhāng Xiǎomíng) and other diacritics never hide a name
  return out.normalize('NFKD').replace(/\p{M}+/gu, '').normalize('NFKC').replace(ZERO_WIDTH, '').toLowerCase()
}

// Identifiers of science, not of people: their digits must never be read as a record number
// (HP:0012345 shares its digits with a record number MZ0012345).
const SCIENTIFIC_IDS = new RegExp([
  String.raw`\b(?:hp|omim|mim|orpha|orphanet|mondo|doid|ncit|efo|uberon|go|hgnc|cl|chebi|mp)\s*[:_]\s*\d+`,
  String.raw`\brs\d+`,
  String.raw`\b(?:nm|nr|np|nc|ng|nt|xm|xp|enst|ensg|ensp|ccds|lrg)_?\d+(?:\.\d+)?`,
  String.raw`\bpmid\s*:?\s*\d+`, String.raw`\bpmc\d+`, String.raw`\bnct\d{8}\b`, String.raw`\bchictr-?\w+`,
  String.raw`\bchr(?:\d{1,2}|x|y|mt?)\s*[:-]\s*[\d,]+(?:\s*[-_]\s*[\d,]+)?`,
  String.raw`(?<![\w-])(?:\d{1,2}|x|y|mt?)\s*:\s*[\d,]+(?:\s*[-_]\s*[\d,]+)?`,
  String.raw`\b[cgpnmr]\.\S+`,
  String.raw`\b10\.\d{4,9}/\S+`,
].join('|'), 'g')

/** The folded text with separators kept (for word-boundary checks), removed, and digits only. */
function forms(text: string): { spaced: string; tight: string; digits: string; scrubbed: string } {
  const folded = fold(text)
  // scientific identifiers removed: what a record number, a date or an ID number could hide in
  const scrubbed = folded.replace(SCIENTIFIC_IDS, ' ').replace(/[\s_]+/g, ' ')
  return {
    spaced: folded.replace(/[\s_]+/g, ' '),
    tight: folded.replace(SEPARATORS, ''),
    digits: scrubbed.replace(/\D+/g, ''),
    scrubbed,
  }
}

const HAS_CJK = /[㐀-鿿豈-﫿]/

const MAX_DEPTH = 12
const MAX_LEAVES = 5000

/**
 * Every string in the input: values, and the keys of free-form objects (another server's
 * `{ data: { "Zhang Wei": 1 } }` carries the name in a key). A key is checked as a whole value,
 * never joined to its neighbours, so a schema key such as "variant" matches nothing.
 * Past MAX_DEPTH or MAX_LEAVES the input is too large to read in full: `TOO_DEEP` is returned
 * among the leaves so the caller refuses instead of passing what it did not read.
 */
export const TOO_DEEP = '\u0000zebra:unread'
function leaves(value: unknown, depth = 0, out: string[] = []): string[] {
  if (depth > MAX_DEPTH || out.length > MAX_LEAVES) {
    if (out[out.length - 1] !== TOO_DEEP) out.push(TOO_DEEP)
    return out
  }
  if (typeof value === 'string') out.push(value)
  else if (Array.isArray(value)) for (const v of value) leaves(v, depth + 1, out)
  else if (value && typeof value === 'object') {
    for (const [k, v] of Object.entries(value)) {
      if (!/^[a-z][a-z0-9_]*$/.test(k)) out.push(k) // a key that is not an identifier-shaped field name
      leaves(v, depth + 1, out)
    }
  } else if (typeof value === 'number' || typeof value === 'boolean') out.push(String(value))
  return out
}

const LONG_DIGITS = /\d{6,}/
// written as one number, or in its 6-8-4 / 6-4-2-2-4 groups; never two separate numbers joined
const CN_RESIDENT_ID = /(?<![\d])[1-9]\d{5}[ -]?(?:18|19|20)\d{2}[ -]?(?:0[1-9]|1[0-2])[ -]?(?:0[1-9]|[12]\d|3[01])[ -]?\d{3}[\dx](?![\d])/
// Separators inside are allowed; a letter next to it is not, so rs13812345678 is a variant id
const CN_MOBILE = /(?<![a-z0-9])(?:\+?0{0,2}86[\s-]?)?1[3-9]\d(?:[\s-]?\d){8}(?![\s-]?\d)(?![a-z0-9])/

const MONTHS = ['january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'september', 'october',
  'november', 'december']

/** [year, month, day] readings of a registered date; d/m/y and m/d/y are both kept when ambiguous. */
function dateParts(identifier: string): Array<[string, number, number]> {
  const t = identifier.trim()
  const ymd = /^(\d{4})\s*[-/.年]\s*(\d{1,2})\s*[-/.月]\s*(\d{1,2})\s*日?$/.exec(t)
  if (ymd) return [[ymd[1] as string, Number(ymd[2]), Number(ymd[3])]]
  const xy = /^(\d{1,2})[-/.](\d{1,2})[-/.](\d{4})$/.exec(t)
  if (!xy) return []
  const a = Number(xy[1])
  const b = Number(xy[2])
  return [[xy[3] as string, b, a], [xy[3] as string, a, b]] // day-month-year, month-day-year
}

/**
 * Patterns for a registered date in the ways people write one: 2019-03-02, 2019/3/2, 20190302,
 * 2019年3月2日, 02/03/2019, 3/2/2019, March 2, 2019, 2 Mar 2019. Digit boundaries on both sides,
 * so the same digits inside a coordinate or another date (12 March) never match.
 */
function dateMatchers(identifier: string): RegExp[] {
  const out: RegExp[] = []
  for (const [y, m, d] of dateParts(identifier)) {
    if (m < 1 || m > 12 || d < 1 || d > 31) continue
    const mo = `0?${m}`
    const da = `0?${d}`
    const sep = String.raw`\s*[-/.]\s*`
    const month = MONTHS[m - 1] as string
    const name = `(?:${month}|${month.slice(0, 3)}\\.?)`
    const ord = '(?:st|nd|rd|th)?'
    const yy = y.slice(2)
    const loose = String.raw`[\s\-/.,]*`
    out.push(
      new RegExp(String.raw`(?<!\d)${y}\s*[-/.年]\s*${mo}\s*[-/.月]\s*${da}(?!\d)`),
      new RegExp(String.raw`(?<!\d)${y}${String(m).padStart(2, '0')}${String(d).padStart(2, '0')}(?!\d)`),
      new RegExp(String.raw`(?<!\d)${da}${sep}${mo}${sep}${y}(?!\d)`),
      new RegExp(String.raw`(?<!\d)${mo}${sep}${da}${sep}${y}(?!\d)`),
      new RegExp(String.raw`(?<!\d)${yy}[-/.]${mo}[-/.]${da}(?!\d)`), // 19-03-02
      new RegExp(String.raw`(?<![a-z])${name}${loose}${da}${ord}${loose}${y}(?!\d)`), // March 2, 2019 / Mar-02-2019
      new RegExp(String.raw`(?<!\d)${da}${ord}${loose}(?:of\s*)?${name}${loose}${y}(?!\d)`), // 2 March 2019 / 02-Mar-2019
    )
  }
  return out
}

function escapeRe(text: string): string {
  return text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')
}

/**
 * What in `input` must not leave the machine, described without repeating it.
 * `identifiers` are the case's registered strings (names, dates of birth, record numbers).
 */
export function guardInput(input: unknown, identifiers: readonly string[]): string | undefined {
  const all = leaves(input)
  if (all.includes(TOO_DEEP)) return 'more nested data than the gate can read in full'
  const values = all.filter(v => !isLocalPath(v))
  if (values.length === 0) return undefined
  // Each value is folded on its own: digits or letters from two unrelated values
  // (an MRN and an HPO id next to it) must never join into a match.
  const folded = values.map(forms)

  for (let i = 0; i < identifiers.length; i++) {
    const raw = (identifiers[i] ?? '').trim()
    if (!raw) continue
    const label = `protected identifier #${i + 1} from the case`
    const id = forms(raw)
    // A short number (a year, a floor, an age) identifies nobody and would block
    // ordinary queries; a record number is longer. Dates are handled below.
    if (/^\d+$/.test(id.tight) && id.tight.length < 6) continue
    const dates = dateMatchers(raw)
    const tokens = id.spaced.split(' ').filter(t => t.length >= 2)

    for (const f of folded) {
      // a date of birth, however it is written
      if (dates.some(re => re.test(f.scrubbed))) return label
      // a record or ID number: digits only, so separators cannot hide it
      // (not for a date: its digits joined across a value's other numbers would match by accident)
      if (dates.length === 0 && LONG_DIGITS.test(id.digits) && f.digits.includes(id.digits)) return label
      // a number-only identifier is matched by its digits alone (a PMID or rsID with the same digits is not it)
      if (/^\d+$/.test(id.tight)) continue
      if (HAS_CJK.test(raw)) {
        // Chinese names: two characters identify, and spacing carries no meaning
        if (id.tight.length >= 2 && f.tight.includes(id.tight)) return label
        continue
      }
      // every token present as a whole word, in any order ("Xiaoming Zhang" too)
      if (tokens.length > 1 && tokens.every(t => wordIn(t, f.spaced))) return label
      if (id.tight.length >= 3) {
        if (wordIn(id.tight, f.spaced)) return label
        if (tokens.length <= 1 && id.tight.length >= 6 && f.tight.includes(id.tight)) return label
      }
    }
  }

  // Patterns that identify a person even when no case is open, checked per value.
  for (const f of folded) {
    if (CN_RESIDENT_ID.test(f.scrubbed)) return 'what looks like a Chinese resident ID number'
    if (CN_MOBILE.test(f.scrubbed)) return 'what looks like a mobile phone number'
    if (emailIn(f.spaced)) return 'an email address'
  }
  return undefined
}

/** An email address in `text` that is not a host login (ssh user@host) or a git remote. */
function emailIn(text: string): boolean {
  // anchored at each '@' and bounded on both sides: a 200 kB sequence without one costs nothing
  for (let at = text.indexOf('@'); at >= 0; at = text.indexOf('@', at + 1)) {
    const left = /[a-z0-9._%+-]{1,64}$/.exec(text.slice(Math.max(0, at - 64), at))
    const right = /^[a-z0-9-]+(?:\.[a-z0-9-]+)*\.[a-z]{2,24}/.exec(text.slice(at + 1, at + 256))
    if (!left || !right) continue
    const addr = `${left[0]}@${right[0]}`
    const m = { index: at - left[0].length }
    if (/^git@/.test(addr) || /@(?:github|gitlab|bitbucket)\.(?:com|org)$/.test(addr)) continue
    const before = text.slice(Math.max(0, m.index - 400), m.index)
    if (/\b(?:ssh|scp|sftp|rsync|mosh|ssh-copy-id|autossh)\b[^;|&]*$/.test(before)) continue
    return true
  }
  return false
}

function wordIn(word: string, text: string): boolean {
  return new RegExp(`(?:^|[^a-z0-9])${escapeRe(word)}(?:[^a-z0-9]|$)`).test(text)
}

/** A bare local path (a case folder, a records file): its text never leaves the machine by itself. */
export function isLocalPath(value: string): boolean {
  const v = value.trim()
  if (!v || /\s/.test(v) || v.includes('://') || /[?=&#]/.test(v)) return false
  return /^(?:~|\$HOME|\.{1,2})?\//.test(v) || /^[A-Za-z]:\\/.test(v)
}

// ---------------------------------------------------------------- shell commands

const NETWORK_TOOL =
  /\b(?:curl|wget|nc|ncat|socat|telnet|scp|sftp|rsync|ssh|mosh|ftp|lftp|gh|glab|aws|gsutil|azcopy|rclone|mail|mailx|sendmail|mutt|osascript|dig|nslookup|host|s2f|httpie|http|xh|aria2c)\b|\bgit\s+(?:push|clone|fetch|pull|remote|ls-remote)\b|\b(?:npm|pnpm|yarn)\s+(?:publish|install|i|add)\b|\bnpx\b|\bpip3?\s+install\b|https?:\/\/|\/dev\/(?:tcp|udp)\/|\bopen\b[^|;&]*\b(?:mailto|https?):|\bdocker\s+(?:run|push|exec)\b|\b(?:Invoke-WebRequest|Invoke-RestMethod|iwr|irm)\b/i
const SCRIPT_RUN = /\bpython3?(?:\.\d+)?\b|\bnode\b|\bdeno\b|\bbun\b|\bruby\b|\bperl\b|\bphp\b|\blua\b|\bjulia\b|\bjava\b|\bRscript\b|\bR\s+(?:-e|--vanilla|-f)\b|\buvx?\b|\bpipx\s+run\b|\b(?:pwsh|powershell)\b/
const ZEBRA_CLI = /(?:^|[\s;&|(/])zebra(?![\w-])/
// `zebra case …` reads and writes the local case only (its HPO checks send ids, never text),
// so it stays usable while the gate is closed and its arguments are not outbound content.
// Only when those are the leading words of the segment, as the shell will read them.
function isZebraLocal(segment: string): boolean {
  const w = shellWords(segment)
  let i = 0
  while (i < w.length && /^[A-Za-z_][A-Za-z0-9_]*=/.test(w[i] as string)) i++ // VAR=value prefixes
  if (/(?:^|\/)python3?(?:\.\d+)?$/.test(w[i] ?? '')) i++
  if (!/(?:^|\/)zebra$/.test(w[i] ?? '')) return false
  i++
  while (i < w.length) {
    const a = w[i] as string
    if (a === '--json') i++
    else if (a === '--case') i += 2
    else if (a.startsWith('--case=')) i++
    else break
  }
  // `case <action>` reads and writes the case; `case recheck` asks public databases again.
  // No readable action (`zebra case $(…)`) is not local. `report export` writes files here only.
  const action = w[i + 1] ?? ''
  if (w[i] === 'case') return /^[a-z][a-z-]*$/.test(action) && action !== 'recheck'
  return w[i] === 'report' && action === 'export'
}

/** The words of a zebra invocation in this segment, from the subcommand on; null when it is not zebra. */
function zebraWords(segment: string): string[] | null {
  const w = shellWords(segment)
  const i = w.findIndex(a => /(?:^|\/)zebra$/.test(a))
  return i < 0 ? null : w.slice(i + 1)
}
// flags whose values stay on this machine in any zebra command (sample names, file paths)
const ZEBRA_LOCAL_FLAGS = new Set(['--case', '--proband', '--mother', '--father', '--sibling', '--out', '--genes',
  '--ped', '--csv', '--workspace', '--vcf'])

/** Split a command into the pieces that run on their own, so one exempt piece cannot cover the rest. */
export function shellSegments(command: string): string[] {
  const out: string[] = []
  let depth = 0
  let current = ''
  let quote: '"' | "'" | null = null
  for (let i = 0; i < command.length; i++) {
    const c = command[i] as string
    const two = command.slice(i, i + 2)
    // inside quotes, ; & | and newlines are text (curl "…?a=1&b=2"), not separators
    if (quote) {
      if (c === '\\' && quote === '"' && i + 1 < command.length) {
        current += c + command[++i]
        continue
      }
      if (c === quote) quote = null
      current += c
      continue
    }
    if (c === '"' || c === "'") {
      quote = c
      current += c
      continue
    }
    if (c === '`') {
      // a backtick substitution runs on its own
      if (current.trim()) out.push(current)
      current = ''
      continue
    }
    if (two === '$(' || two === '<(' || two === '>(') {
      depth++
      i++
      if (current.trim()) out.push(current)
      current = ''
      continue
    }
    if (c === ')' && depth > 0) {
      depth--
      if (current.trim()) out.push(current)
      current = ''
      continue
    }
    if (two === '&&' || two === '||') {
      i++
      if (current.trim()) out.push(current)
      current = ''
      continue
    }
    if (c === ';' || c === '|' || c === '&' || c === '\n') {
      if (current.trim()) out.push(current)
      current = ''
      continue
    }
    current += c
  }
  if (current.trim()) out.push(current)
  return out
}

/**
 * A segment as the shell will read it: quotes removed and words rejoined, so `c'u'rl` and
 * `"cu""rl"` read as curl, and the text of `bash -c "…"` or `eval "…"` is seen as the command it is.
 */
function asRun(segment: string): string {
  return shellWords(segment).join(' ')
}

function outboundSegment(segment: string): boolean {
  if (isZebraLocal(segment)) return false // the local case only
  for (const text of [segment, asRun(segment)]) {
    if (NETWORK_TOOL.test(text)) return true
    if (ZEBRA_CLI.test(text)) return true // zebra's other commands query public databases
    if (SCRIPT_RUN.test(text)) return true // a script can reach the network; its source is not read here
  }
  // `bash -c "…"`, `sh -c`, `eval`: judge the inner command the same way
  const w = shellWords(segment)
  const inner = /^(?:ba|z|da|k)?sh$|^eval$/.test(w[0] ?? '') ? w.slice(1).filter(a => a !== '-c').join(' ') : ''
  return inner !== '' && inner !== segment && isOutboundShell(inner)
}

/** True when a shell command may send data off the machine. Judged per segment. */
export function isOutboundShell(command: string): boolean {
  return shellSegments(command).some(outboundSegment)
}

/**
 * The text of a shell command that may reach the network. Once any part of it does, every part is
 * scanned — a pipe (echo … | curl), a heredoc, a variable set before the call or a script's source
 * all carry text to the outbound piece — except what demonstrably stays here: segments that only
 * touch the local case, the values of zebra's local-only flags (sample names, case paths) inside
 * zebra commands, local paths, and the folder given to `cd`.
 */
export function outboundText(command: string): string {
  const segments = shellSegments(command)
  if (!segments.some(outboundSegment)) return ''
  const kept: string[] = []
  for (const segment of segments) {
    if (isZebraLocal(segment)) continue
    const words = shellWords(segment)
    const zebra = zebraWords(segment) !== null
    if (words[0] === 'cd') continue
    for (let i = 0; i < words.length; i++) {
      const w = words[i] as string
      if (zebra && ZEBRA_LOCAL_FLAGS.has(w)) {
        i++
        continue
      }
      if (zebra && /^--[\w-]+=/.test(w) && ZEBRA_LOCAL_FLAGS.has(w.split('=')[0] as string)) continue
      if (isLocalPath(w)) continue
      kept.push(w)
    }
  }
  return kept.join(' ')
}

const GENOME_FILE = /\.(?:g?vcf(?:\.gz|\.bgz)?|bcf|bam|cram|sam|fastq(?:\.gz)?|fq(?:\.gz)?|bed(?:\.gz)?)(?![a-z0-9])/i
const UPLOADER =
  /\b(?:scp|sftp|rclone|gsutil|azcopy)\b|\baws\s+s3\b|\bgh\s+(?:gist|release)\s+\w+|\brsync\b[^|;&]*\S+:\S*|\bcurl\b[^|;&]*(?:\s-T\b|--upload-file|\s-F\b|--form|--data-binary|\s-d\s*@|--data(?:-raw|-urlencode)?\s*@)|\bwget\b[^|;&]*--post-file|\b(?:nc|ncat|socat)\b|\bssh\b[^|;&]*\bcat\s*>|\bmail(?:x)?\b|\bpython3?\s+-m\s+http\.server\b/i

// folders a desktop client uploads on its own: copying a file there sends it
const SYNC_FOLDER = /(?:^|[\s'"=])(?:~|\$HOME|\/Users\/[^/\s]+|\/home\/[^/\s]+)?\/?(?:Dropbox|Google Drive|GoogleDrive|My Drive|OneDrive[^/\s]*|Library\/Mobile Documents|Library\/CloudStorage|iCloud Drive|Nutstore[^/\s]*|坚果云[^/\s]*|BaiduNetdisk[^/\s]*|百度网盘[^/\s]*|WeDrive[^/\s]*|Box Sync|pCloud Drive|MEGA)(?:\/|['"\s]|$)/i
const COPY_INTO = /\b(?:cp|mv|rsync|ditto|ln|install|tee)\b/

/** True when a zebra command would send a whole VCF's variant list to a web service (triage's MyVariant prefilter). */
export function sendsVariantList(command: string): boolean {
  // argparse accepts any unique prefix of a long option: --pref myvariant, --pre=myvariant
  const isFlag = (a: string) => a.length >= 5 && '--prefilter'.startsWith(a)
  return shellSegments(command).some(seg => {
    const rest = zebraWords(seg)
    if (!rest) return false
    return rest.some((a, i) => {
      const [flag, value] = a.includes('=') ? [a.slice(0, a.indexOf('=')), a.slice(a.indexOf('=') + 1)] : [a, rest[i + 1]]
      return isFlag(flag) && value === 'myvariant'
    })
  })
}

/** True when a command looks like it sends a raw genome file to another machine. */
export function uploadsGenome(command: string): boolean {
  if (GENOME_FILE.test(command) && UPLOADER.test(command)) return true
  // copying genome data into a synced folder uploads it as surely as scp
  if (GENOME_FILE.test(command) && COPY_INTO.test(command) && SYNC_FOLDER.test(command)) return true
  // an archive that holds genome data, then sent: zip -r g.zip genome/ && curl -F f=@g.zip …
  const archive = /\b(?:zip|7z|7za|tar|gzip|bgzip|xz|zstd)\b[^|;&]*?([\w.-]+\.(?:zip|7z|tar(?:\.gz|\.bz2|\.xz|\.zst)?|tgz|gz))\b/i.exec(command)
  if (archive && (GENOME_FILE.test(command) || /\b(?:genome|exome|wes|wgs|vcf|bam|cram|fastq)s?\b/i.test(command)) && UPLOADER.test(command)) return true
  const segments = shellSegments(command)
  if (segments.some(s => GENOME_FILE.test(s)) && segments.some(s => UPLOADER.test(s))) return true
  // a whole directory going out (scp -r genome_dir host:, tar czf - dir | ssh)
  return /\bscp\s+-r\b|\btar\b[^|;&]*\bc/.test(command) && /\b(?:ssh|scp|rclone)\b|\baws\s+s3\b/.test(command)
}

const FILE_ARGS: RegExp[] = [
  /--(?:upload-file|post-file)[=\s]+("[^"]+"|'[^']+'|\S+)/gi,
  /\s-T\s*("[^"]+"|'[^']+'|\S+)/gi,
  /(?:--data-binary|--data-urlencode|--data-raw|--data|-d|-F|--form)\s*(?:[^@\s]*@)("[^"]+"|'[^']+'|[^\s;|&]+)/gi,
  /<\s*("[^"]+"|'[^']+'|[^\s;|&]+)/g,
]
// commands whose non-option, non-remote arguments are local files being sent
const COPY_TOOLS = /^(?:scp|sftp|rsync|rclone|gsutil|azcopy)$/
const REMOTE_ARG = /^(?:[\w.-]+@)?[\w.-]+:(?!\/\/)|^(?:s3|gs|az|https?):\/\/|^[\w-]+:$/

/**
 * The local file paths a command would send somewhere, relative paths resolved against a
 * preceding `cd` in the same command, so the caller can check whether they lie in the case folder.
 */
export function uploadedPaths(command: string): string[] {
  const found = new Set<string>()
  // a file copied into a synced folder is uploaded by the desktop client
  for (const segment of shellSegments(command)) {
    const w = shellWords(segment)
    if (COPY_INTO.test(w[0] ?? '') && w.length >= 3 && SYNC_FOLDER.test(` ${w[w.length - 1]}`)) {
      for (const a of w.slice(1, -1)) if (!a.startsWith('-')) found.add(a)
    }
  }
  let cwd = ''
  const segments = shellSegments(command)
  const join = (p: string) => (cwd && !/^(?:\/|~|\$HOME)/.test(p) ? `${cwd.replace(/\/$/, '')}/${p}` : p)
  // files read on the left of a pipe that ends in an uploader reading stdin
  const pipeReaders: string[] = []
  const archived: string[] = []
  const scripted = SCRIPT_RUN.test(command)
  for (const segment of segments) {
    const words = shellWords(segment)
    const head = words[0] ?? ''
    if (head === 'cd' && words[1]) {
      cwd = words[1]
      continue
    }
    if (/^(?:cat|zcat|gzip|bgzip|tar|base64|xxd|head|tail)$/.test(head)) {
      for (const w of words.slice(1)) if (!w.startsWith('-')) pipeReaders.push(join(w))
    }
    for (const re of FILE_ARGS) {
      re.lastIndex = 0
      let m: RegExpExecArray | null
      while ((m = re.exec(segment)) !== null) {
        const capture = m[1] ?? ''
        const parts = /^['"]/.test(capture) ? [capture] : capture.split(/\s+/)
        for (const part of parts) {
          const clean = part.replace(/^['"]|['"]$/g, '')
          if (clean === '-') pipeReaders.forEach(r => found.add(r))
          else if (clean && !/^https?:/.test(clean)) found.add(join(clean))
        }
      }
    }
    if (/^(?:zip|7z|7za|tar)$/.test(head)) {
      // what goes into an archive leaves with it: zip -r out.zip ~/cases/x, tar czf out.tgz -C ~/cases x
      let base = ''
      for (let k = 1; k < words.length; k++) {
        const a = words[k] as string
        if (a === '-C' && words[k + 1]) {
          base = words[++k] as string
          continue
        }
        if (a.startsWith('-') || /\.(?:zip|7z|tar|tgz|gz|bz2|xz|zst)$/.test(a)) continue
        archived.push(base && !/^(?:\/|~|\$HOME)/.test(a) ? `${base.replace(/\/$/, '')}/${a}` : join(a))
      }
    }
    if (scripted && zebraWords(segment) === null) {
      // a script that opens a file from the case folder sends its contents: python -c "…open('…')…",
      // or a heredoc body; only paths that start a word (not the tail of records/x or host:/tmp)
      for (const m of segment.replace(/https?:\/\/\S+/g, ' ').matchAll(/(?<=^|[\s'"=(,])(?:~|\$HOME)?\/[^\s'"<>|;&(),]{2,}/g)) found.add(m[0])
    }
    if (head === 'gh' && words[1] === 'release' && words[2] === 'upload') {
      for (const w of words.slice(4)) if (!w.startsWith('-')) found.add(join(w))
    }
    if (COPY_TOOLS.test(head) || (head === 'aws' && words[1] === 's3') || (head === 'gh' && words[1] === 'gist')) {
      const args = words.slice(head === 'aws' || head === 'gh' ? 3 : 1)
      for (const w of args) {
        if (w.startsWith('-') || REMOTE_ARG.test(w)) continue
        found.add(join(w))
      }
    }
  }
  // the archive's sources count once anything in the command sends something somewhere
  if (archived.length && (UPLOADER.test(command) || segments.some(outboundSegment))) archived.forEach(a => found.add(a))
  return [...found]
}

/** A tiny shell-words splitter, so `/zebra new "my case" A title` works. */
export function shellWords(text: string): string[] {
  const out: string[] = []
  let current = ''
  let quote: '"' | "'" | null = null
  let started = false
  for (let i = 0; i < text.length; i++) {
    const c = text[i] as string
    if (quote) {
      if (c === quote) quote = null
      else if (c === '\\' && quote === '"' && i + 1 < text.length) current += text[++i]
      else current += c
      continue
    }
    if (c === '"' || c === "'") {
      quote = c as '"' | "'"
      started = true
      continue
    }
    if (c === '\\' && i + 1 < text.length) {
      current += text[++i]
      started = true
      continue
    }
    if (/\s/.test(c)) {
      if (current || started) out.push(current)
      current = ''
      started = false
      continue
    }
    current += c
    started = true
  }
  if (current || started) out.push(current)
  return out
}
