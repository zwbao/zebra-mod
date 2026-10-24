#!/bin/sh
# Build the README promo from tools/promo/:
#   promo.html --render.mjs--> frames; score.py + mix.py --> master.wav;
#   docs/media/zebra-mod-promo.mp4 (1080p, with sound), .gif (silent README preview), -poster.png
# Needs node >= 22, Google Chrome, ffmpeg, and a venv with numpy scipy soundfile (tools/promo/.venv).
set -eu
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
work=${1:-$(mktemp -d)}
frames="$work/promo-frames"
out="$root/docs/media"
mkdir -p "$out"
[ -f "$frames/f00899.png" ] || node "$here/render.mjs" "$frames" 30
[ -f "$work/music.wav" ] || "$here/.venv/bin/python" "$here/score.py" "$work/music.wav"
"$here/.venv/bin/python" "$here/mix.py" "$work/music.wav" "$work/master.wav"
ffmpeg -loglevel error -y -framerate 30 -i "$frames/f%05d.png" -i "$work/master.wav" \
  -c:v libx264 -preset slow -crf 18 -pix_fmt yuv420p -c:a aac -b:a 192k -shortest -movflags +faststart \
  "$out/zebra-mod-promo.mp4"
ffmpeg -loglevel error -y -framerate 30 -i "$frames/f%05d.png" \
  -vf "fps=10,scale=800:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=80:stats_mode=diff[p];[b][p]paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle" \
  "$out/zebra-mod-promo.gif"
cp "$frames/f00270.png" "$out/zebra-mod-promo-poster.png"
ls -la "$out"
