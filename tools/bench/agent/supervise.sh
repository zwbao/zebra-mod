#!/bin/sh
# Keep a benchmark run going across usage-limit windows: run, and when run.py stops after
# consecutive failures (a usage limit, an outage), wait and run again, until every case of the
# selection has a result.json or MAX_ROUNDS is reached.
#   supervise.sh LOG N_EXPECTED OUT_DIR -- <run.py args...>
LOG=$1; N=$2; OUT=$3; shift 4
WAIT=${WAIT:-900}; MAX_ROUNDS=${MAX_ROUNDS:-60}
HERE=$(cd "$(dirname "$0")" && pwd)
i=0
while [ $i -lt $MAX_ROUNDS ]; do
  i=$((i+1))
  echo "== round $i $(date)" >> "$LOG"
  python3 -I "$HERE/run.py" --out "$OUT" "$@" >> "$LOG" 2>&1
  have=$(ls "$OUT"/cases/*/result.json 2>/dev/null | wc -l | tr -d ' ')
  echo "== round $i ended $(date): $have/$N results" >> "$LOG"
  [ "$have" -ge "$N" ] && { echo "ALL DONE $have/$N"; exit 0; }
  sleep "$WAIT"
done
echo "gave up after $MAX_ROUNDS rounds: $have/$N"
