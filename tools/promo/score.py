"""zebra-mod promo score — 30 s, 144 BPM (one shot = 6 beats = 2.5 s), D minor resolving to F major.
Code-original synthesis (numpy/scipy/soundfile), deterministic (seed 1878, Muybridge's year).
Techniques adapted from the guizang-product-video-skill cinematic example; the piece, key, tempo
and every section time are derived from this film's shot list (tools/promo/DIRECTION.md):

  0.0–5.0   heartbeat sub pulses over a low drone (the patient); small hit at 2.5 ("~5 years")
  5.0–7.5   riser + reverse swell, a near-silent gap at 7.35
  7.5       impact, pad bloom, bell motif (brand)
 10.0–25.0  montage: four-on-the-floor, clap on 2 and 4, 8th plucks; one chord per shot,
            a light hit on every cut (10, 12.5, 15 … 22.5)
 25.0–27.5  breakdown: drums out, heartbeat returns under the pad ("every name protected")
 27.5–30.0  final impact, resolve to F (add9), bells, tail

Usage: .venv/bin/python score.py music.wav
"""
import sys
import numpy as np
import soundfile as sf
from scipy.signal import butter, sosfilt, fftconvolve

SR = 48000; BPM = 144; BEAT = 60 / BPM; SHOT = 6 * BEAT; DUR = 30.0; N = int(SR * DUR)
rng = np.random.default_rng(1878)
buses = {k: np.zeros((N, 2)) for k in ['drums', 'bass', 'pad', 'arp', 'bell', 'fx', 'heart']}


def hz(n): return 440 * 2 ** ((n - 69) / 12)
def tt(d): return np.arange(int(d * SR)) / SR
def put(bus, t0, sig, pan=0.0, gain=1.0):
    i = int(round(t0 * SR))
    if i >= N: return
    if sig.ndim == 1: sig = np.stack([sig * np.sqrt((1 - pan) / 2), sig * np.sqrt((1 + pan) / 2)], 1)
    j = min(N, i + len(sig)); buses[bus][max(0, i):j] += sig[max(0, -i):j - i] * gain
def env(t, a, d, r_at=None, r=0.3):
    e = np.minimum(1, t / max(a, 1e-4)) * np.exp(-t * d)
    if r_at is not None: e *= np.clip((r_at + r - t) / r, 0, 1)
    return e
def lp(x, fc, order=2): return sosfilt(butter(order, min(fc, SR * .45), 'low', fs=SR, output='sos'), x, axis=0)
def hp(x, fc, order=2): return sosfilt(butter(order, fc, 'high', fs=SR, output='sos'), x, axis=0)
def bp(x, lo, hi): return sosfilt(butter(2, [lo, hi], 'band', fs=SR, output='sos'), x, axis=0)
def sweep_lp(x, f0, f1, block=1024):
    out = np.zeros_like(x); zi = None; nb = int(np.ceil(len(x) / block))
    for b in range(nb):
        fc = f0 * (f1 / f0) ** (b / max(1, nb - 1))
        sos = butter(2, min(fc, SR * .45), 'low', fs=SR, output='sos')
        if zi is None: zi = np.zeros((sos.shape[0], 2) + x.shape[1:])
        y, zi = sosfilt(sos, x[b * block:(b + 1) * block], axis=0, zi=zi); out[b * block:(b + 1) * block] = y
    return out


def saw(f, t, k_max=30, phase=0.0):
    s = np.zeros_like(t)
    for k in range(1, int(min(k_max, (SR * .45) // f)) + 1): s += np.sin(2 * np.pi * k * f * t + phase * k) / k
    return s
def supersaw(note, dur, amp, cut=1800, bright_to=None, attack=.5, release=.7):
    t = tt(dur + release); f = hz(note); L = np.zeros_like(t); R = np.zeros_like(t)
    for i, c in enumerate([-12, -6, 0, 6, 12]):
        v = saw(f * 2 ** (c / 1200), t, 22, rng.uniform(0, 6.28)); p = (i - 2) / 2.2
        L += v * np.sqrt((1 - p) / 2); R += v * np.sqrt((1 + p) / 2)
    x = sweep_lp(np.stack([L, R], 1) / 5, cut, bright_to or cut)
    return x * (np.minimum(1, t / attack) * np.clip((dur + release - t) / release, 0, 1))[:, None] * amp
def pluck(note, amp, dur=.7):
    t = tt(dur); f = hz(note); s = np.zeros_like(t)
    for k in range(1, int(min(24, (SR * .45) // f)) + 1): s += np.sin(2 * np.pi * k * f * t) / k * np.exp(-t * (5 + k * 1.8))
    return s * np.minimum(1, t / .002) * amp
def bell(note, amp, dur=2.4):
    t = tt(dur); f = hz(note)
    return np.sin(2 * np.pi * f * t + 2.0 * np.exp(-t * 3) * np.sin(2 * np.pi * f * 3.5 * t)) * np.exp(-t * 1.7) * np.minimum(1, t / .003) * amp
def sub(note, dur, amp):
    t = tt(dur); f = hz(note); s = np.sin(2 * np.pi * f * t) + .22 * np.sin(4 * np.pi * f * t)
    return np.tanh(1.5 * s) * env(t, .005, 0, dur - .06, .06) * amp
def kick(amp):
    t = tt(.4); ph = 2 * np.pi * (48 * t + (140 - 48) * .04 * (1 - np.exp(-t / .04)))
    return np.tanh(1.4 * (np.sin(ph) * np.exp(-t * 8) + hp(rng.uniform(-1, 1, len(t)), 2500) * np.exp(-t * 200) * .3)) * amp
def heartbeat(amp):
    """Lub-dub: two soft, low thumps — the patient under the whole film."""
    def thump(f0, a):
        t = tt(.32); ph = 2 * np.pi * (f0 * t + 30 * .05 * (1 - np.exp(-t / .05)))
        return lp(np.sin(ph) * np.exp(-t * 11) * a, 220)
    x = np.zeros(int(.6 * SR)); a = thump(52, 1.0); b = thump(46, .7)
    x[:len(a)] += a; o = int(.2 * SR); x[o:o + len(b)] += b
    return x * amp
def clap(amp):
    t = tt(.3); n = rng.uniform(-1, 1, len(t)); e = np.zeros_like(t)
    for o in [0, .01, .021]: e += (t >= o) * np.exp(-np.maximum(0, t - o) * (60 if o < .02 else 15))
    return bp(n * e, 900, 3400) * amp
def hat(amp):
    t = tt(.05); return hp(rng.uniform(-1, 1, len(t)), 7500, 4) * np.exp(-t * 80) * amp
def impact(amp):
    t = tt(3.0); ph = 2 * np.pi * (30 * t + (75 - 30) * .22 * (1 - np.exp(-t / .22)))
    return np.tanh(1.3 * (np.sin(ph) * np.exp(-t * 1.7) + lp(rng.uniform(-1, 1, len(t)), 900) * np.exp(-t * 6) * .8)) * amp
def riser(dur, amp, f0=300, f1=9000):
    t = tt(dur); x = hp(sweep_lp(rng.uniform(-1, 1, (len(t), 2)), f0, f1), 200)
    glide = np.sin(2 * np.pi * (220 * t + (880 - 220) * t ** 2 / (2 * dur)))[:, None] * .1
    return (x + glide) * ((t / dur) ** 2.2)[:, None] * amp
def swell(dur, amp):
    t = tt(dur); return hp(rng.uniform(-1, 1, (len(t), 2)), 3000) * ((t / dur) ** 3)[:, None] * amp


# D minor → F major. (bass root, voicing)
CH = {'Dm9': (38, [53, 57, 60, 64]), 'Bb': (34, [53, 58, 62, 65]), 'F': (41, [53, 57, 60, 65]),
      'C': (36, [52, 55, 60, 64]), 'Bbmaj7': (34, [53, 57, 62, 65]), 'Fadd9': (41, [55, 57, 60, 65, 69])}
MONTAGE = [(10.0, 'Dm9'), (12.5, 'Bb'), (15.0, 'F'), (17.5, 'C'), (20.0, 'Dm9'), (22.5, 'Bb')]
ARP = [0, 2, 1, 3, 2, 1]

# 0–7.5 hook: drone + heartbeat, hit at 2.5, riser into a gap
put('pad', 0.0, supersaw(50, 7.3, .08, cut=260, bright_to=1500, attack=.8, release=.05))
put('bass', 0.0, sub(26, 7.28, .09))
for k in range(int(7.3 / (2 * BEAT))):
    put('heart', k * 2 * BEAT, heartbeat(.55 + .25 * (k * 2 * BEAT) / 7.3))
put('fx', 2.5, impact(.42)); put('fx', 1.9, swell(.6, .14))
for i, n in enumerate([62, 65, 69]): put('bell', 2.52 + i * .09, bell(n, .05), (i - 1) * .4)
put('fx', 4.6, riser(2.75, .34)); put('fx', 6.3, swell(1.05, .22))
# 7.5 brand
put('fx', 7.5, impact(.95))
for i, n in enumerate(CH['F'][1] + [69]): put('pad', 7.5, supersaw(n, 2.4, .05, cut=500, bright_to=3600, attack=.04), (i - 2) * .3)
put('bass', 7.5, sub(29, 2.3, .2))
for i, (t0, n) in enumerate([(7.5, 72), (7.92, 77), (8.33, 81), (8.75, 84), (9.17, 81), (9.58, 77)]):
    put('bell', t0, bell(n, .14), -.4 + .16 * i)
put('fx', 9.0, riser(1.0, .16, 900, 8000))
# 10–25 montage
for t0, name in MONTAGE:
    root, voic = CH[name]
    put('fx', t0, impact(.38)); put('fx', t0 - .55, swell(.55, .12))
    for i, n in enumerate(voic): put('pad', t0, supersaw(n, SHOT, .042, cut=900, bright_to=2400, attack=.12, release=.25), (i - 1.5) * .3)
    for b in range(6):
        tb = t0 + b * BEAT
        put('drums', tb, kick(.78))
        if b % 2 == 1: put('drums', tb, clap(.26), .05)
        for h in range(2): put('drums', tb + h * BEAT / 2, hat(.06 if h else .04), .3 if h else -.3)
        put('bass', tb + BEAT / 2, sub(root + (12 if b % 3 == 2 else 0), BEAT / 2 - .02, .24))
        for s8 in range(2):
            n = voic[ARP[(b * 2 + s8) % 6] % len(voic)] + 12
            put('arp', tb + s8 * BEAT / 2, pluck(n, .075 if s8 == 0 else .055), -.35 if s8 else .35)
put('fx', 24.0, riser(1.0, .14, 600, 6000))
# 25–27.5 breakdown
put('fx', 25.0, impact(.5))
for i, n in enumerate(CH['Bbmaj7'][1]): put('pad', 25.0, supersaw(n, 2.4, .05, cut=1400, bright_to=700, attack=.05, release=.3), (i - 1.5) * .3)
put('bass', 25.0, sub(34, 2.3, .18))
for k in range(3): put('heart', 25.0 + k * 2 * BEAT, heartbeat(.6))
for i, (t0, n) in enumerate([(25.2, 74), (25.83, 77), (26.45, 81)]): put('bell', t0, bell(n, .09), -.3 + .3 * i)
put('fx', 26.7, riser(.8, .2, 700, 9000)); put('fx', 27.0, swell(.5, .2))
# 27.5–30 resolve
put('fx', 27.5, impact(1.0))
for i, n in enumerate(CH['Fadd9'][1]): put('pad', 27.5, supersaw(n, 2.2, .05, cut=2600, bright_to=800, attack=.03, release=.6), (i - 2) * .25)
put('bass', 27.5, sub(29, 2.3, .22))
for t0, n in [(27.5, 77), (27.92, 81), (28.33, 84), (28.75, 89), (29.17, 84)]: put('bell', t0, bell(n, .11), .1)

# mix
t_all = np.arange(N) / SR
def sidechain(depth, rel):
    g = np.ones(N)
    for t0, _ in MONTAGE:
        for b in range(6):
            i = int((t0 + b * BEAT) * SR); j = min(N, i + int(.4 * SR)); x = np.arange(j - i) / SR
            g[i:j] = np.minimum(g[i:j], 1 - depth * np.exp(-x / rel * 2.2))
    return g[:, None]
buses['pad'] *= sidechain(.5, .2); buses['bass'] *= sidechain(.65, .1); buses['arp'] *= sidechain(.3, .12)
def delay(x, time, fb, mix):
    d = int(time * SR); y = x.copy()
    for k in range(1, 6):
        tap = np.zeros_like(x); tap[k * d:] = x[:N - k * d]
        if k % 2: tap = tap[:, ::-1]
        y += tap * mix * fb ** (k - 1)
    return y
buses['arp'] = delay(lp(buses['arp'], 5000), BEAT * .75, .36, .3)
buses['bell'] = delay(buses['bell'], BEAT, .3, .24)
ti = tt(2.6); ir = np.stack([rng.normal(0, 1, len(ti)), rng.normal(0, 1, len(ti))], 1) * np.exp(-ti * 2.5)[:, None]
ir = lp(ir, 6000); ir /= np.sqrt((ir ** 2).sum())
send = buses['pad'] * .35 + buses['bell'] * .6 + buses['arp'] * .3 + buses['fx'] * .3 + hp(buses['drums'], 400) * .1
wet = np.stack([fftconvolve(send[:, c], ir[:, c])[:N] for c in range(2)], 1) * .85
mix = (buses['drums'] * .95 + buses['bass'] * .62 + buses['pad'] * .9 + buses['arp'] * .8 + buses['bell'] * .85
       + buses['fx'] * .85 + buses['heart'] * .9 + wet)
mix = hp(mix, 30)
mix *= (np.clip((DUR - t_all) / 1.6, 0, 1) ** 1.4 * np.minimum(1, t_all / .03))[:, None]
mix = np.tanh(mix * 1.1) / np.tanh(1.1)
mix /= np.abs(mix).max() / .9
out = sys.argv[1] if len(sys.argv) > 1 else 'music.wav'
sf.write(out, mix.astype(np.float32), SR, subtype='FLOAT')
print(out, f'{DUR}s', BPM, 'BPM')
