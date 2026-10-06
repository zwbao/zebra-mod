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
    if (out === before) break
  }
  return out.normalize('NFKC').replace(ZERO_WIDTH, '').toLowerCase()
}

// Identifiers of science, not of people: their digits must never be read as a record number
// (HP:0012345 shares its digits with a record number MZ0012345).
const SCIENTIFIC_IDS = new RegExp([
  String.raw`\b(?:hp|omim|mim|orpha|orphanet|mondo|doid|ncit|efo|uberon|go|hgnc|cl|chebi|mp)\s*[:_]\s*\d+`,
  String.raw`\brs\d+`,
  String.raw`\b(?:nm|nr|np|nc|ng|nt|xm|xp|enst|ensg|ensp|ccds|lrg)_?\d+(?:\.\d+)?`,
  String.raw`\bpmid\s*:?\s*\d+`, String.raw`\bpmc\d+`, String.raw`\bnct\d{8}\b`, String.raw`\bchictr-?\w+`,
  String.raw`\b(?:chr)?(?:\d{1,2}|x|y|mt?)\s*[:-]\s*\d+(?:\s*[-_]\s*\d+)?`,
  String.raw`\b[cgpnmr]\.\S+`,
  String.raw`\b10\.\d{4,9}/\S+`,
].join('|'), 'g')

/** The folded text with separators kept (for word-boundary checks), removed, and digits only. */
function forms(text: string): { spaced: string; tight: string; digits: string } {
  const folded = fold(text)
  return {
    spaced: folded.replace(/[\s_]+/g, ' '),
    tight: folded.replace(SEPARATORS, ''),
    // the digits a record number could hide in: scientific identifiers removed first
    digits: folded.replace(SCIENTIFIC_IDS, ' ').replace(/\D+/g, ''),
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
const CN_RESIDENT_ID = /(?<![\d])[1-9]\d{5}(?:18|19|20)\d{2}(?:0[1-9]|1[0-2])(?:0[1-9]|[12]\d|3[01])\d{3}[\dx](?![\d])/
// Separators inside are allowed; a letter next to it is not, so rs13812345678 is a variant id
const CN_MOBILE = /(?<![a-z0-9])(?:\+?0{0,2}86[\s-]?)?1[3-9]\d(?:[\s-]?\d){8}(?![\s-]?\d)(?![a-z0-9])/
const EMAIL = /[a-z0-9._%+-]+@[a-z0-9-]+(?:\.[a-z0-9-]+)*\.[a-z]{2,}/

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
    out.push(
      new RegExp(String.raw`(?<!\d)${y}\s*[-/.年]\s*${mo}\s*[-/.月]\s*${da}(?!\d)`),
      new RegExp(String.raw`(?<!\d)${y}${String(m).padStart(2, '0')}${String(d).padStart(2, '0')}(?!\d)`),
      new RegExp(String.raw`(?<!\d)${da}${sep}${mo}${sep}${y}(?!\d)`),
      new RegExp(String.raw`(?<!\d)${mo}${sep}${da}${sep}${y}(?!\d)`),
      new RegExp(String.raw`(?<![a-z])${name}\s*${da}${ord}\s*,?\s*${y}(?!\d)`),
      new RegExp(String.raw`(?<!\d)${da}${ord}\s*(?:of\s*)?${name}\s*,?\s*${y}(?!\d)`),
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
      if (dates.some(re => re.test(f.spaced))) return label
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
    if (CN_RESIDENT_ID.test(f.tight) || CN_RESIDENT_ID.test(f.digits)) return 'what looks like a Chinese resident ID number'
    if (CN_MOBILE.test(f.spaced)) return 'what looks like a mobile phone number'
    if (emailIn(f.spaced)) return 'an email address'
  }
  return undefined
}

/** An email address in `text` that is not a host login (ssh user@host) or a git remote. */
function emailIn(text: string): boolean {
  for (const m of text.matchAll(new RegExp(EMAIL.source, 'g'))) {
    const addr = m[0]
    if (/^git@/.test(addr) || /@(?:github|gitlab|bitbucket)\.(?:com|org)$/.test(addr)) continue
    const before = text.slice(0, m.index)
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
  if (!v || /\s/.test(v) || v.includes('://')) return false
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
  // `case recheck` asks public databases again: not local
  return w[i] === 'case' && w[i + 1] !== 'recheck'
}
// flags whose values stay on this machine in any zebra command (sample names, file paths)
const ZEBRA_LOCAL_FLAGS = new Set(['--case', '--proband', '--mother', '--father', '--sibling', '--out', '--genes',
  '--ped', '--csv', '--workspace', '--vcf'])

/** Split a command into the pieces that run on their own, so one exempt piece cannot cover the rest. */
export function shellSegments(command: string): string[] {
  const out: string[] = []
  let depth = 0
  let current = ''
  for (let i = 0; i < command.length; i++) {
    const c = command[i] as string
    const two = command.slice(i, i + 2)
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
 * The text of a shell command that may reach the network: only its outbound segments, with local
 * paths and the values of zebra's local-only flags removed — a sample name or a case path is not
 * sent anywhere, so it must not trip the gate.
 */
export function outboundText(command: string): string {
  const kept: string[] = []
  for (const segment of shellSegments(command)) {
    if (!outboundSegment(segment)) continue
    const words = shellWords(segment)
    for (let i = 0; i < words.length; i++) {
      const w = words[i] as string
      if (ZEBRA_LOCAL_FLAGS.has(w)) {
        i++
        continue
      }
      if (/^--[\w-]+=/.test(w) && ZEBRA_LOCAL_FLAGS.has(w.split('=')[0] as string)) continue
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
  return shellSegments(command).some(seg => {
    const w = shellWords(seg)
    const i = w.findIndex(a => /(?:^|\/)zebra$/.test(a))
    if (i < 0) return false
    const rest = w.slice(i + 1)
    const at = rest.indexOf('--prefilter')
    return rest.includes('--prefilter=myvariant') || (at >= 0 && rest[at + 1] === 'myvariant')
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
    if (COPY_TOOLS.test(head) || (head === 'aws' && words[1] === 's3') || (head === 'gh' && words[1] === 'gist')) {
      const args = words.slice(head === 'aws' || head === 'gh' ? 3 : 1)
      for (const w of args) {
        if (w.startsWith('-') || REMOTE_ARG.test(w)) continue
        found.add(join(w))
      }
    }
  }
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
