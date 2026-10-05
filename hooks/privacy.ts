// The privacy gate's two checks. Both read a tool call's input as text; both
// are best effort over spellings, which is why the gate also asks the person
// before raw genome files leave the machine instead of trying to be clever.

const CN_RESIDENT_ID = /(?<![\dA-Za-z])[1-9]\d{5}(?:18|19|20)\d{2}(?:0[1-9]|1[0-2])(?:0[1-9]|[12]\d|3[01])\d{3}[\dXx](?![\dA-Za-z])/
const CN_MOBILE = /(?<![\d])1[3-9]\d{9}(?![\d])/
const EMAIL = /[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)*\.[A-Za-z]{2,}/

/** What in `input` must not leave the machine, described without repeating it; undefined when nothing. */
export function guardInput(input: unknown, identifiers: readonly string[]): string | undefined {
  let text: string
  try {
    text = typeof input === 'string' ? input : JSON.stringify(input) ?? ''
  } catch {
    return undefined
  }
  const lower = text.toLowerCase()
  for (let i = 0; i < identifiers.length; i++) {
    const id = identifiers[i]
    if (id && id.length >= 2 && lower.includes(id.toLowerCase())) return `protected identifier #${i + 1} from the case`
    // JSON.stringify escapes non-ASCII only when asked; also try the escaped form
    if (id && text.includes(JSON.stringify(id).slice(1, -1))) return `protected identifier #${i + 1} from the case`
  }
  if (CN_RESIDENT_ID.test(text)) return 'what looks like a Chinese resident ID number'
  if (CN_MOBILE.test(text)) return 'what looks like a mobile phone number'
  if (EMAIL.test(text) && !/git@github\.com/.test(text)) return 'an email address'
  return undefined
}

const GENOME_FILE = /\.(?:g?vcf(?:\.gz|\.bgz)?|bcf|bam|cram|sam|fastq(?:\.gz)?|fq(?:\.gz)?|bed\.gz)(?![A-Za-z0-9])/i
const UPLOADER =
  /\b(?:scp|sftp|rclone|gsutil|azcopy)\b|\baws\s+s3\b|\bgh\s+release\s+upload\b|\brsync\b[^|;&]*\S+:\S*|\bcurl\b[^|;&]*(?:\s-T\b|--upload-file|\s-F\b|--form|--data-binary|\s-d\s*@|--data(?:-raw|-urlencode)?\s*@)|\bwget\b[^|;&]*--post-file/i

/** True when a shell command looks like it sends a raw genome file to another machine. */
export function uploadsGenome(command: string): boolean {
  return UPLOADER.test(command) && GENOME_FILE.test(command)
}
