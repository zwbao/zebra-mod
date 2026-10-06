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

/** The folded text with separators kept (for word-boundary checks), removed, and digits only. */
function forms(text: string): { spaced: string; tight: string; digits: string } {
  const folded = fold(text)
  return {
    spaced: folded.replace(/[\s_]+/g, ' '),
    tight: folded.replace(SEPARATORS, ''),
    digits: folded.replace(/\D+/g, ''),
  }
}

const HAS_CJK = /[㐀-鿿豈-﫿]/

/** Every string VALUE in the input, keys never: a key like "variant" must not match an identifier. */
function leaves(value: unknown, depth = 0, out: string[] = []): string[] {
  if (depth > 8 || out.length > 2000) return out
  if (typeof value === 'string') out.push(value)
  else if (Array.isArray(value)) for (const v of value) leaves(v, depth + 1, out)
  else if (value && typeof value === 'object') for (const v of Object.values(value)) leaves(v, depth + 1, out)
  else if (typeof value === 'number' || typeof value === 'boolean') out.push(String(value))
  return out
}

const LONG_DIGITS = /\d{6,}/
const CN_RESIDENT_ID = /(?<![\d])[1-9]\d{5}(?:18|19|20)\d{2}(?:0[1-9]|1[0-2])(?:0[1-9]|[12]\d|3[01])\d{3}[\dx](?![\d])/
// Separators inside are allowed; a letter next to it is not, so rs13812345678 is a variant id
const CN_MOBILE = /(?<![a-z0-9])(?:\+?0{0,2}86[\s-]?)?1[3-9]\d(?:[\s-]?\d){8}(?![\s-]?\d)(?![a-z0-9])/
const EMAIL = /[a-z0-9._%+-]+@[a-z0-9-]+(?:\.[a-z0-9-]+)*\.[a-z]{2,}/
// user@host in an ssh-family command or a git remote is not a person's address
const HOST_LOGIN = /\b(?:ssh|scp|sftp|rsync|git|mosh|ansible|kubectl)\b|git@|@(?:github|gitlab|bitbucket)\.com\b/

function dateRenderings(identifier: string): string[] {
  const m = /^(\d{4})[-/.](\d{1,2})[-/.](\d{1,2})$/.exec(identifier.trim())
  if (!m) return []
  const y = m[1] as string
  const mo = m[2] as string
  const d = m[3] as string
  const mm = mo.padStart(2, '0')
  const dd = d.padStart(2, '0')
  return [`${y}${mm}${dd}`, `${y}-${mm}-${dd}`, `${y}/${mm}/${dd}`, `${y}.${mm}.${dd}`, `${y}年${mo}月${d}日`]
}

function escapeRe(text: string): string {
  return text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')
}

/**
 * What in `input` must not leave the machine, described without repeating it.
 * `identifiers` are the case's registered strings (names, dates of birth, record numbers).
 */
export function guardInput(input: unknown, identifiers: readonly string[]): string | undefined {
  const values = leaves(input).filter(v => !isLocalPath(v))
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
    const renderings = dateRenderings(raw).map(r => forms(r).tight).filter(Boolean)
    const tokens = id.spaced.split(' ').filter(t => t.length >= 2)

    for (const f of folded) {
      // a date of birth, however it is written
      if (renderings.some(r => f.tight.includes(r))) return label
      // a record or ID number: digits only, so separators cannot hide it
      if (LONG_DIGITS.test(id.digits) && f.digits.includes(id.digits)) return label
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
    if (EMAIL.test(f.spaced) && !HOST_LOGIN.test(f.spaced)) return 'an email address'
  }
  return undefined
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
  /\b(?:curl|wget|nc|ncat|socat|telnet|scp|sftp|rsync|ssh|mosh|ftp|lftp|gh|glab|aws|gsutil|azcopy|rclone|mail|mailx|sendmail|osascript|dig|nslookup|host)\b|\bgit\s+(?:push|clone|fetch|pull|remote|ls-remote)\b|\bnpm\s+(?:publish|install)\b|\bnpx\b|\bpip3?\s+install\b|https?:\/\//
const SCRIPT_RUN = /\bpython3?(?:\.\d+)?\b|\bnode\b|\bdeno\b|\bbun\b|\bruby\b|\bperl\b/
const ZEBRA_CLI = /(?:^|[\s;&|(/])zebra(?![\w-])/
// `zebra case …` reads and writes the local case only (its HPO checks send ids, never text),
// so it stays usable while the gate is closed and its arguments are not outbound content
const ZEBRA_LOCAL = /(?:^|[\s;&|(/])zebra(?![\w-])(?:\s+--?(?:json|case)(?:[=\s]+\S+)?)*\s+case\b/
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
    if (two === '$(') {
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

/** True when a shell command may send data off the machine. Judged per segment. */
export function isOutboundShell(command: string): boolean {
  for (const segment of shellSegments(command)) {
    if (ZEBRA_LOCAL.test(segment)) continue // the local case only
    if (NETWORK_TOOL.test(segment)) return true
    if (ZEBRA_CLI.test(segment)) return true // zebra's other commands query public databases
    if (SCRIPT_RUN.test(segment)) return true // a script can reach the network; its source is not read here
  }
  return false
}

/**
 * The text of a shell command that may reach the network: only its outbound segments, with local
 * paths and the values of zebra's local-only flags removed — a sample name or a case path is not
 * sent anywhere, so it must not trip the gate.
 */
export function outboundText(command: string): string {
  const kept: string[] = []
  for (const segment of shellSegments(command)) {
    if (ZEBRA_LOCAL.test(segment)) continue
    if (!NETWORK_TOOL.test(segment) && !ZEBRA_CLI.test(segment) && !SCRIPT_RUN.test(segment)) continue
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

/** True when a command looks like it sends a raw genome file to another machine. */
export function uploadsGenome(command: string): boolean {
  if (GENOME_FILE.test(command) && UPLOADER.test(command)) return true
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
