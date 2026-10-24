"""Mix the score with action sound effects, ducking the music under each effect, and write the master.

Effects are the guizang-product-video-skill built-in WAVs (code-original, 48 kHz mono, peak ~ -20 dBFS),
copied into tools/promo/sfx/. Each cue: (action time on screen, file, landmark offset, gain dB, pan).
`at` = action time - landmark: clicks and pops align their onset, whooshes and sweeps their peak.

Usage: .venv/bin/python mix.py music.wav master.wav
"""
import sys
from pathlib import Path
import numpy as np
import soundfile as sf

SR = 48000
HERE = Path(__file__).resolve().parent
SFX = HERE / "sfx"
LAND = {  # measured with sfx_landmarks.py: (onset, peak)
    "click-alt": (.001, .007), "click": (.001, .008), "ding-dong": (.002, .206), "pop": (.001, .006),
    "resolve": (.003, .232), "success": (.003, .213), "sweep": (.135, .342), "toggle": (.001, .080),
    "typing": (.001, .496), "whoosh": (.071, .202),
}
PEAK_ALIGNED = {"whoosh", "sweep"}

CUES = []
def cue(t, name, gain_db, pan=0.0):
    CUES.append((t, name, gain_db, pan))

cue(7.45, "whoosh", 16)                                       # stripe wipe into the logo
for i in range(7): cue(12.6 + i * .2083 + .08, "pop", 11, (-1) ** i * .4)   # phenotype chips
cue(14.4, "resolve", 8)                                       # the differential locks
for i in range(3): cue(15.2 + i * .06, "click", 13, (i - 1) * .3)  # ACMG codes
for k in range(1, 14):                                        # the points count up
    x = k / 13; tk = 15.6 + .85 * (1 - (1 - x) ** (1 / 4))
    cue(tk, "click-alt", 11 + 3 * x, .1)
cue(16.55, "toggle", 12)                                      # "Pathogenic" lands
for i in range(5): cue(17.6 + i * .2083 + .1, "pop", 10, (i - 2) * .3)   # five model cards
cue(18.95, "sweep", 12)                                       # the cards light up
for i in range(4): cue(20.2 + i * .2083 + .1, "pop", 10, (i - 1.5) * .3)  # drugs
for i in range(5): cue(21.1 + i * .1, "click", 11, (i - 2) * .3)          # agency badges
cue(22.75, "whoosh", 12)                                      # reports fan out
cue(25.05, "success", 10)                                      # the shield
for k in range(2): cue(27.75 + k * .48, "typing", 16)         # the install command
cue(28.9, "ding-dong", 4)                                     # call to action


def load(name):
    x, sr = sf.read(SFX / f"{name}.wav", dtype="float64")
    assert sr == SR, name
    return x if x.ndim == 1 else x.mean(1)


def main(music_path, out_path):
    music, sr = sf.read(music_path, dtype="float64")
    assert sr == SR
    n = len(music)
    fx = np.zeros((n, 2))
    duck = np.ones(n)
    for t, name, gain_db, pan in CUES:
        x = load(name) * 10 ** (gain_db / 20)
        onset, peak = LAND[name]
        at = t - (peak if name in PEAK_ALIGNED else onset)
        i = int(round(at * SR)); j = min(n, i + len(x))
        if i >= n: continue
        st = np.stack([x * np.sqrt((1 - pan) / 2), x * np.sqrt((1 + pan) / 2)], 1)[: j - i]
        fx[max(0, i):j] += st[max(0, -i):]
        # duck the music 3 dB around each effect, with 60 ms ramps
        a, b = max(0, i - int(.06 * SR)), min(n, j + int(.06 * SR))
        duck[a:b] = np.minimum(duck[a:b], 10 ** (-3 / 20))
    k = int(.06 * SR)
    duck = np.convolve(duck, np.ones(k) / k, mode="same")
    mix = music * duck[:, None] + fx
    mix = np.tanh(mix * 1.05) / np.tanh(1.05)
    mix /= max(1e-9, np.abs(mix).max()) / .9
    sf.write(out_path, mix.astype(np.float32), SR, subtype="FLOAT")
    sf.write(Path(out_path).with_name("sfx-stem.wav"), (fx / max(1e-9, np.abs(fx).max()) * .9).astype(np.float32), SR)
    print(out_path, len(CUES), "cues")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
