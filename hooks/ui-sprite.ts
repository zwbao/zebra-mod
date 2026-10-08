// The galloping zebra drawn in the band above the prompt while zebra-mod queries a database.
// Eleven positions of one gallop stride, traced from Eadweard Muybridge's "The Horse in Motion" (1878,
// public domain): rider removed, frames aligned on the back (croup to withers) so the body stays
// level while the legs move, given zebra stripes, drawn in braille (2 x 4 dots to a terminal cell). Regenerate with
// `python3 tools/zebra_sprite/build.py`. One colour, the terminal's own text colour, so the zebra reads
// on light and dark themes alike; the ground is a dim dotted line that runs backwards.
// Pure functions of (frame, scroll): no state here.

const FRAMES: readonly (readonly string[])[] = [
  [
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⢠⡴⣟⣶⣤⡀⠀',
    '⠀⠀⠀⠀⢀⣀⢠⣤⡖⣤⢠⣤⡆⣶⣾⡟⣵⡏⠉⠉⠁⠀',
    '⠀⠈⠷⠾⠛⠁⢻⡟⣼⣿⢸⣿⡇⣿⢏⣾⠟⠀⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠀⡟⢟⠘⠋⠈⠘⠛⠃⢻⡿⠟⠀⠀⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠴⠃⠀⠱⠀⠀⠒⠦⠤⠟⠀⠀⠀⠀⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀',
  ],
  [
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣰⣾⣳⣦⣀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⢀⣤⣶⣤⢠⣤⡄⣦⣿⡟⣽⠉⠉⠙⠀⠀',
    '⠀⢴⣤⡶⠟⠁⢿⡟⣼⣿⢸⣿⡇⣿⢏⣾⠃⠀⠀⠀⠀⠀',
    '⠀⠀⠉⠀⠀⠀⢀⣼⣿⠏⠘⠛⠃⣧⢿⡋⠀⠀⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠀⠈⢻⢓⣢⣀⡼⠳⠔⠛⠀⠀⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀',
  ],
  [
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣤⡿⢶⣤⡀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⢀⣴⣶⣦⢠⣤⡄⣦⣿⡟⡵⠉⠉⠋⠀⠀',
    '⠀⢠⣤⠾⠛⠃⢿⡟⣼⣿⢸⣿⡇⣿⢏⣾⠃⠀⠀⠀⠀⠀',
    '⠀⠀⠁⠀⠀⠀⠀⢨⣿⣿⠘⠛⠇⣧⠿⠟⣄⠀⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠀⠈⠉⠙⠛⠦⡤⠏⠰⠖⠁⠀⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀',
  ],
  [
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⡷⡶⣤⡀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⢀⣤⣤⣤⢠⣤⡄⣦⣿⡟⣵⠛⠛⠿⠂⠀',
    '⠀⢀⣴⠿⠻⠃⢿⡟⣼⣿⢸⣿⡇⣿⢏⣾⡇⠀⠀⠀⠀⠀',
    '⠀⠙⠁⠀⠀⠀⠈⠘⢿⣿⠸⠿⠇⠧⣿⣟⠔⢲⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠀⠀⠘⢿⡓⠒⠦⠜⠶⠚⠁⠸⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢙⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀',
  ],
  [
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣄⡶⡶⣦⡀⠀⠀',
    '⠀⠀⠀⢀⣀⣀⢀⣤⣤⣤⢠⣤⡄⣦⣿⡟⣵⠋⠛⠻⠂⠀',
    '⠀⢠⣾⠟⠛⠁⣿⡟⣼⣿⢸⣿⡇⣿⢏⣾⡇⠀⠀⠀⠀⠀',
    '⠀⠋⠁⠀⠀⠀⠈⣼⠟⢿⡘⠛⠃⠧⠟⠛⢼⠿⢄⡀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠰⡏⠀⠈⠙⠢⢤⡀⠀⠘⠋⠀⠀⠙⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠀⠳⠤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀',
  ],
  [
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⡦⣣⡷⣄⠀⠀',
    '⠀⠀⠀⣀⣀⣀⢀⣤⡄⣤⢠⣤⡄⣶⣾⡟⣵⠟⠙⠻⠷⠀',
    '⠀⣴⠟⠋⠉⢡⣿⡟⣼⣿⢸⣿⡇⣿⣏⣾⡏⠀⠀⠀⠀⠀',
    '⠀⠁⠀⠀⢀⣠⠟⠘⣿⠉⠘⠛⠃⠛⠛⠛⠸⢿⡇⠀⠀⠀',
    '⠀⠀⠀⢀⡞⠁⠀⠀⠻⡀⠀⠀⠀⠀⠀⠀⠀⠐⠋⠓⠦⠄',
    '⠀⠀⠀⠸⠷⠀⠀⠀⠀⠹⠶⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀',
  ],
  [
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⡰⣠⡶⣄⠀⠀',
    '⠀⠀⠀⠀⣀⣀⢠⣤⣤⣤⢀⣤⡄⣶⣾⡟⣵⠟⠙⠻⠷⠀',
    '⠤⠾⠛⠋⠉⢠⣿⡟⣼⣿⢸⣿⡇⣿⢏⣾⡏⠀⠀⠀⠀⠀',
    '⠀⠀⠀⣠⡠⣾⡟⠈⠉⠉⠘⠛⠃⠛⠛⠛⣜⠛⠒⢤⡀⠀',
    '⠠⠤⠞⠁⢰⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢧⡀⠀⠙⠂',
    '⠀⠀⠀⠀⠘⠓⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⠂⠀⠀',
  ],
  [
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⡴⣣⣶⣄⠀⠀',
    '⠀⠀⠀⠀⣀⡤⢠⣤⣤⣤⢠⣤⡄⣶⣾⡟⣵⡿⠛⠺⠷⠀',
    '⠛⠻⠟⠛⠉⢠⣿⡟⣼⣿⢸⣿⡇⣿⣏⣾⡿⠀⠀⠀⠀⠀',
    '⠀⢀⣠⠴⣢⡿⠛⠈⠀⠉⠘⠛⠃⠛⠻⡟⠘⠷⣤⡀⠀⠀',
    '⠉⠉⢀⡔⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣧⠀⠀⠀⠈⠉⠋',
    '⠀⠀⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠂⠀⠀⠀⠀⠀',
  ],
  [
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⢠⡰⣲⢶⣄⠀⠀',
    '⢀⢀⠀⢀⣠⣤⢠⣤⣤⣤⢠⣤⡄⣶⣾⡟⣵⡿⠛⠻⠷⠀',
    '⠈⠛⠛⠛⠉⢠⣿⡟⣼⣿⢸⣿⡇⣿⢏⣾⡟⠁⠀⠀⠀⠀',
    '⠀⠀⢀⣼⡧⠟⠛⠈⠁⠉⠘⠛⠇⣿⠿⠛⢴⡄⠀⠀⠀⠀',
    '⠀⠛⠋⠁⠀⠀⠀⠀⠀⠀⠀⠀⡰⠋⠀⠀⠀⠙⢢⡀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠃⠀⠀⠀⠀⠀⠀⠉⠋⠀',
  ],
  [
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⢀⡠⣞⣶⣄⠀⠀',
    '⠀⢀⡀⢀⣀⣤⢠⣴⡖⣤⢠⣤⡆⣶⣾⡟⣵⡟⠙⠛⠃⠀',
    '⠀⠀⠻⠟⠋⠁⣿⡟⣼⣿⢸⣿⡇⣿⢏⣾⠏⠀⠀⠀⠀⠀',
    '⠀⠀⠀⢀⠞⢻⠛⠘⠉⠈⠘⠛⠃⣿⣿⠋⠀⠀⠀⠀⠀⠀',
    '⠀⠀⠈⠉⠀⠰⠁⠀⠀⠀⠀⠀⠼⠁⠿⠀⠀⠀⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀',
  ],
  [
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⢠⡴⣞⣶⣄⠀⠀',
    '⠀⠀⡀⣀⣀⣀⢀⣤⡄⣤⢠⣤⡄⣶⣾⡟⣵⡟⠙⠛⠃⠀',
    '⠀⠘⠿⠿⠟⠃⢹⡟⣼⣿⢸⣿⡇⣿⢏⣾⡟⠀⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠀⠐⢞⣼⠟⠉⠘⠛⠃⣿⡿⠟⠀⠀⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠈⡏⠙⢢⠀⣀⠿⠿⠇⠀⠀⠀⠀⠀⠀⠀',
    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀',
  ],
]

const DEFAULT = 0x01000000 // the terminal's own colour
const GROUND = 0x6e6c64
const SPRITE_COLS = 22
const STANDING = 5 // the pose shown still (the welcome card, the shield)

/** Columns and rows of the zebra's Raster: the sprite with a little ground before and behind it. */
export const RASTER_COLUMNS = 28
export const RASTER_ROWS = 6
export const FRAME_COUNT = FRAMES.length
const OFFSET = 3

/** The braille dot bits of the cell at (col, row): the zebra's, then the ground's in the lowest dot row. */
function cell(frame: number, scroll: number, col: number, row: number, isStanding: boolean): { horse: number; ground: number } {
  const rows = FRAMES[isStanding ? STANDING : frame % FRAMES.length] as readonly string[]
  const sc = col - OFFSET
  const ch = sc >= 0 && sc < SPRITE_COLS ? (rows[row] as string).codePointAt(sc) ?? 0x2800 : 0x2800
  const horse = ch - 0x2800
  let ground = 0
  if (row === RASTER_ROWS - 1) {
    // dot 7 (left column, bottom) and dot 8 (right column, bottom), one in every six dots
    const s = isStanding ? 0 : scroll
    if ((col * 2 + s) % 6 === 0) ground |= 0x40
    if ((col * 2 + 1 + s) % 6 === 0) ground |= 0x80
  }
  return { horse, ground }
}

/** The Raster `cells` of one frame: braille dots in the terminal's colour, the ground dimmed. */
export function zebraCells(frame: number, scroll: number, isStanding = false): string {
  const words = new Uint32Array(RASTER_COLUMNS * RASTER_ROWS * 3)
  let i = 0
  for (let row = 0; row < RASTER_ROWS; row++) {
    for (let col = 0; col < RASTER_COLUMNS; col++) {
      const { horse, ground } = cell(frame, scroll, col, row, isStanding)
      const bits = horse | ground
      words[i++] = bits ? 0x2800 + bits : 0x20
      words[i++] = horse ? DEFAULT : ground ? GROUND : DEFAULT
      words[i++] = DEFAULT
    }
  }
  return base64(new Uint8Array(words.buffer))
}

const B64 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'

/** Standard padded base64 (the module runs with no Node and no DOM, so no Buffer or btoa to lean on). */
export function base64(bytes: Uint8Array): string {
  let out = ''
  let i = 0
  for (; i + 2 < bytes.length; i += 3) {
    const n = ((bytes[i] as number) << 16) | ((bytes[i + 1] as number) << 8) | (bytes[i + 2] as number)
    out += B64[(n >> 18) & 63]! + B64[(n >> 12) & 63]! + B64[(n >> 6) & 63]! + B64[n & 63]!
  }
  const rest = bytes.length - i
  if (rest === 1) {
    const n = (bytes[i] as number) << 16
    out += B64[(n >> 18) & 63]! + B64[(n >> 12) & 63]! + '=='
  } else if (rest === 2) {
    const n = ((bytes[i] as number) << 16) | ((bytes[i + 1] as number) << 8)
    out += B64[(n >> 18) & 63]! + B64[(n >> 12) & 63]! + B64[(n >> 6) & 63]! + '='
  }
  return out
}
