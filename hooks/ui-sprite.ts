// The galloping zebra drawn in the band above the prompt while zebra-mod queries a database.
// Four frames of a 28 x 12 pixel sprite, drawn two pixels to a terminal cell with half blocks,
// on a ground of dots that scrolls left. Pure functions of (frame, scroll): no state here.

// . empty  W body  K stripe and mane  E eye  N muzzle  O outline  L leg  H hoof  D dust
const FRAMES: readonly (readonly string[])[] = [
  [
    '.....................KKO....', '....................KWWWOOO.', '...................KWWEWWWWO', '...OO.OO.OO.OO.OO.KKWWWWWWNN',
    '..OWWKWWKWWKWWKWWKWWKWOOOO..', '.OWWWKWWKWWKWWKWWWWWWO......', '.KWWWWKWWKWWKWWKWWWWWO......', 'K.OWWWKWWKWWKWWKWWWWO.......',
    '....L..L........L..L........', '.D.K...K........K...K.......', 'D.L.....L......L.....L......', 'DH.......H....H.......H.....',
  ],
  [
    '.....................KKO....', '....................KWWWOOO.', '...................KWWEWWWWO', '...OO.OO.OO.OO.OO.KKWWWWWWNN',
    '..OWWKWWKWWKWWKWWKWWKWOOOO..', '.OWWWKWWKWWKWWKWWWWWWO......', '.KWWWWKWWKWWKWWKWWWWWO......', 'K.OWWWKWWKWWKWWKWWWWO.......',
    '......L.L.......L.L.........', '.......K.K.....K.K..........', '........L.L...L..L..........', '.........HH..H..H...........',
  ],
  [
    '.....................KKO....', '....................KWWWOOO.', '...................KWWEWWWWO', '...OO.OO.OO.OO.OO.KKWWWWWWNN',
    '..OWWKWWKWWKWWKWWKWWKWOOOO..', '.OWWWKWWKWWKWWKWWWWWWO......', '.KWWWWKWWKWWKWWKWWWWWO......', 'K.OWWWKWWKWWKWWKWWWWO.......',
    '.....L.L.........L.L........', '.....KK...........K.K.......', '...DL.L............L.L......', '..D.HH..............H.H.....',
  ],
  [
    '.....................KKO....', '....................KWWWOOO.', '...................KWWEWWWWO', '...OO.OO.OO.OO.OO.KKWWWWWWNN',
    '..OWWKWWKWWKWWKWWKWWKWOOOO..', '.OWWWKWWKWWKWWKWWWWWWO......', '.KWWWWKWWKWWKWWKWWWWWO......', 'K.OWWWKWWKWWKWWKWWWWO.......',
    '.....L.L.........L.L........', '.....K..K........K.K........', '.....L..L........L..L.......', '......H.H.........H..H......',
  ],
]

// colors read on a light and on a dark terminal alike: the outline carries the shape on white
const PALETTE: Record<string, number> = {
  W: 0xf2efe8, K: 0x222222, E: 0x0a0a0a, N: 0x373737, O: 0x787876, L: 0xc4c0b8, H: 0x3c3c3c, D: 0xa09687, G: 0x6e6c64,
}
const DEFAULT = 0x01000000
const SPRITE_W = 28
const SPRITE_H = 12

/** Columns and rows of the zebra's Raster: the sprite with a little ground before and behind it. */
export const RASTER_COLUMNS = 36
export const RASTER_ROWS = SPRITE_H / 2
export const FRAME_COUNT = FRAMES.length
const OFFSET_X = 3

/** The pixel at (x, y) of `frame` scrolled by `scroll`: a color, or undefined where nothing is drawn. */
function pixel(frame: number, scroll: number, x: number, y: number, isStanding: boolean): number | undefined {
  const rows = FRAMES[isStanding ? 3 : frame % FRAMES.length] as readonly string[]
  const sx = x - OFFSET_X
  const code = sx >= 0 && sx < SPRITE_W ? (rows[y] as string)[sx] : '.'
  if (code !== undefined && code !== '.' && !(isStanding && code === 'D')) return PALETTE[code]
  // the ground: a dotted line under the hooves that moves backwards while the zebra runs
  if (y === SPRITE_H - 1 && (x + (isStanding ? 0 : scroll)) % 5 === 0) return PALETTE.G
  return undefined
}

/** The Raster `cells` of one frame: each cell an upper half block, top pixel as foreground, bottom as background. */
export function zebraCells(frame: number, scroll: number, isStanding = false): string {
  const words = new Uint32Array(RASTER_COLUMNS * RASTER_ROWS * 3)
  let i = 0
  for (let row = 0; row < RASTER_ROWS; row++) {
    for (let col = 0; col < RASTER_COLUMNS; col++) {
      const top = pixel(frame, scroll, col, row * 2, isStanding)
      const bottom = pixel(frame, scroll, col, row * 2 + 1, isStanding)
      if (top === undefined && bottom === undefined) {
        words[i++] = 0x20; words[i++] = DEFAULT; words[i++] = DEFAULT
      } else if (top === undefined) {
        words[i++] = 0x2584; words[i++] = bottom as number; words[i++] = DEFAULT // ▄
      } else {
        words[i++] = 0x2580; words[i++] = top; words[i++] = bottom ?? DEFAULT // ▀
      }
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
