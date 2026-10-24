# zebra-mod promo — direction

**Brief.** 30 s, 16:9 (1920×1080, 30 fps), English headlines with Chinese subtitles in formal written register.
Placed under the README's one-paragraph introduction: an inline GIF (silent, 800 px, 10 fps) linking to the MP4 with sound.
Scope: the whole product as of 0.3.9; platform: the Claude Code mod (TUI) and its `zebra` CLI. Not covered: nothing
else exists (no desktop, mobile or hosted API).

## Three directions considered

| | Axis | Idea | Verdict |
|---|---|---|---|
| A | Product walkthrough | Terminal mock-ups of a session: records → tool rows → ACMG card | First cut. Accurate but reads like documentation; small text at README width; nothing to stop the scroll. |
| B | **Kinetic manifesto** | Problem in numbers → one question → brand → six one-line claims, each with one big visual, on a beat grid | **Chosen.** A README visitor decides in two seconds; big type and one idea per shot survive being scaled to 800 px. |
| C | Patient story | Follow one synthetic child from first seizure to a care plan | Emotional, but 30 s cannot carry a story and six capabilities; risks sentimentality. |

## Devices derived from the product

1. **Zebra stripes as the transition**, because the product is named for the medical aphorism ("think zebras"): the stripe wipe that reveals the logo and the slow stripes behind the end card.
2. **A heartbeat in the score** (lub-dub sub thumps in the opening and in the breakdown), because every case in the tool is a patient: the music starts and returns with a pulse.
3. **Chips that gather into one diagnosis**, because the mod's core move is phenotype terms → a differential: seven real phenotype chips of the demo case converge into "Dravet syndrome · SCN1A".
4. **A counter that the engine owns**: the ACMG points count to 13 with one tick per point, because zebra computes the points in code (`zebra acmg classify PVS1 PS2 PM2_Supporting` → 13, Pathogenic).
5. **Five model cards lit one by one, never merged**, because S2F models are reported axis by axis and never averaged.
6. **A rain of real source names** behind "22 · 35+", because every number on screen is backed by a named database (docs/DATA-SOURCES.md has 35 rows).

## Picture rules

- Background #070b10 with a slow teal glow (#5eead4 / #22d3ee); the palette of the report export (#0f5c6e) as the deep accent.
- Type: SF Pro Display 800 for headlines (112–230 px), PingFang SC 500 for Chinese (38 px, letter-spacing 2 px); mono SF Mono 46 px.
- One idea per shot; every shot 2.5 s = 6 beats at 144 BPM; every cut on a beat with a 0.2 s, 10% white flash.
- Entrances vary by shot: count-up, punch-in, slide-up, back-ease pop, gather, fan, typewriter.
- Frame 0 is the finished opening ("300,000,000"), so a thumbnail taken from the first frame is a complete picture.

## Shot list

| Time (s) | Shot | Picture | Sound |
|---|---|---|---|
| 0.0–2.5 | Problem | "300,000,000 people live with a rare disease." / 全球罕见病患者约 3 亿 | drone, heartbeat |
| 2.5–5.0 | Delay | "~5 years to get a diagnosis. On average." / 平均确诊历时近 5 年 | small hit, bells |
| 5.0–7.5 | Question | "What if it took one conversation?" / 倘若，只需一次对话 | riser, gap at 7.35 |
| 7.5–10.0 | Brand | stripe wipe → 🦓 zebra-mod, "Your rare-disease research team, inside Claude Code." | impact, pad bloom, bell motif, whoosh |
| 10.0–12.5 | Scale | 22 tools · 35+ live data sources, source names drifting | groove starts (Dm9) |
| 12.5–15.0 | Diagnose | Records → Diagnosis; seven phenotype chips gather into "Dravet syndrome · SCN1A" | pops per chip, resolve |
| 15.0–17.5 | Interpret | Variant → Verdict; PVS1 PS2 PM2_Supporting; 0 → 13 points; "Pathogenic" | clicks, one tick per point, stamp |
| 17.5–20.0 | S2F | 5 genome AI models: SpliceAI, Pangolin, GPN-MSA, Evo 2, AlphaGenome | pops, sweep as they light |
| 20.0–22.5 | Treat | Genotype → Therapy; stiripentol · fenfluramine · cannabidiol · SCN1A antisense in trials; FDA · EMA · NMPA · 医保目录 · Trials in China | pops, clicks |
| 22.5–25.0 | Report | Doctor-ready reports: HTML · Word · PDF fan out | whoosh |
| 25.0–27.5 | Trust | 🛡️ Every claim cited. Every name protected. | breakdown, heartbeat returns, success |
| 27.5–30.0 | End | Think zebras · 🦓 zebra-mod · 让每一位罕见病患者都被看见 · `Install https://github.com/zwbao/zebra-mod` · ★ Star on GitHub | final impact, resolve to F, typing, chime |

## Claims and their evidence

| Claim on screen | Evidence |
|---|---|
| ~300 million people live with a rare disease | Nguengang Wakap et al., Eur J Hum Genet 2020 (263–446 million) |
| ~5 years to a diagnosis on average | EURORDIS Rare Barometer (average diagnostic delay of about 5 years) |
| 22 tools | `hooks/tools.ts` registers 22 |
| 35+ live data sources | `docs/DATA-SOURCES.md`, 35 source rows |
| Dravet syndrome · SCN1A; PVS1 PS2 PM2_Supporting = 13, Pathogenic | synthetic demo case (`examples/demo`), `zebra acmg classify` |
| Five sequence models | `zebra/s2f.py`, `zebra/s2f_backends.py` |
| stiripentol, fenfluramine, cannabidiol approved for Dravet; SCN1A antisense in trials | FDA approvals 2018–2020; zorevunersen phase 3 |
| FDA · EMA · NMPA · NRDL · trials in China | the `access` tool |
| HTML · Word · PDF | `zebra/report_export.py` |
| Names never sent | the privacy gate (`hooks/privacy.ts`) |

## Deviation from the skill's component rule

The product is a terminal plug-in: its interface is Claude Code's own TUI rows, which cannot be mounted as browser
components. The first cut (direction A) rebuilt them as terminal mock-ups; the chosen direction uses the product's real
data (the demo case, the engine's numbers, the real source and model names) as large type instead of UI.

## Sound

Code-original score (`score.py`, seed 1878, 144 BPM, D minor → F major), action effects from the skill's built-in
code-original WAVs (`sfx/`, see `sfx/SOURCE.md`), mixed by `mix.py` with 3 dB ducking under each effect.
Checked by spectrogram, loudness (−12.4 LUFS integrated) and effect-versus-music peak levels; **not listened to by a person yet.**
