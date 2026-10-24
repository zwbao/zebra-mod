// Render tools/promo/promo.html to PNG frames with headless Chrome over the DevTools protocol.
// No dependencies (Node >= 22 for the global WebSocket). Then encode with tools/promo/build.sh.
//
//   node tools/promo/render.mjs <out_dir> [fps=30] [chrome path]
import { spawn } from 'node:child_process'
import { mkdirSync, writeFileSync, mkdtempSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join, resolve, dirname } from 'node:path'
import { fileURLToPath, pathToFileURL } from 'node:url'

const here = dirname(fileURLToPath(import.meta.url))
const outDir = resolve(process.argv[2] ?? 'promo-frames')
const fps = Number(process.argv[3] ?? 30)
const chrome = process.argv[4] ?? '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
mkdirSync(outDir, { recursive: true })

const port = 9333 + Math.floor(Math.random() * 500)
const profile = mkdtempSync(join(tmpdir(), 'promo-chrome-'))
const proc = spawn(chrome, ['--headless=new', `--remote-debugging-port=${port}`, `--user-data-dir=${profile}`,
  '--hide-scrollbars', '--force-device-scale-factor=1', '--window-size=1920,1080', '--disable-gpu', 'about:blank'],
  { stdio: 'ignore' })

const sleep = ms => new Promise(r => setTimeout(r, ms))
let target
for (let i = 0; i < 400 && !target; i++) {
  try {
    const list = await (await fetch(`http://127.0.0.1:${port}/json`)).json()
    target = list.find(t => t.type === 'page')
  } catch { /* not listening yet */ }
  if (!target) await sleep(150)
}
if (!target) { proc.kill('SIGKILL'); throw new Error('Chrome did not start') }

const ws = new WebSocket(target.webSocketDebuggerUrl)
await new Promise(r => ws.addEventListener('open', r, { once: true }))
let seq = 0
const pending = new Map()
ws.addEventListener('message', ev => {
  const msg = JSON.parse(ev.data)
  if (msg.id && pending.has(msg.id)) { pending.get(msg.id)(msg); pending.delete(msg.id) }
})
const send = (method, params = {}) => new Promise((res, rej) => {
  const id = ++seq
  pending.set(id, m => (m.error ? rej(new Error(`${method}: ${m.error.message}`)) : res(m.result)))
  ws.send(JSON.stringify({ id, method, params }))
})

await send('Page.enable')
await send('Emulation.setDeviceMetricsOverride', { width: 1920, height: 1080, deviceScaleFactor: 1, mobile: false })
const url = pathToFileURL(join(here, 'promo.html')).href + '?capture=1'
await send('Page.navigate', { url })
await sleep(1500)
await send('Runtime.evaluate', { expression: 'document.fonts.ready.then(() => true)', awaitPromise: true })
const { result } = await send('Runtime.evaluate', { expression: 'window.DURATION', returnByValue: true })
const frames = Math.round(result.value * fps)
for (let f = 0; f < frames; f++) {
  await send('Runtime.evaluate', { expression: `render(${(f / fps).toFixed(5)})` })
  const shot = await send('Page.captureScreenshot', { format: 'png', captureBeyondViewport: false })
  writeFileSync(join(outDir, `f${String(f).padStart(5, '0')}.png`), Buffer.from(shot.data, 'base64'))
  if (f % 60 === 0) process.stdout.write(`frame ${f}/${frames}\n`)
}
ws.close()
proc.kill('SIGKILL')
console.log(`done: ${frames} frames in ${outDir}`)
