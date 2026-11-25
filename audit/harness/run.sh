#!/usr/bin/env bash
# Run a harness program (sweep.py or confirm.py) in one process; restart it when the
# watchdog has ended it on a timed-out task (exit 1 with the in-flight file present).
#   bash run.sh sweep.py --targets scipy.special
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${AUDIT_PY:-$HERE/../.venv/Scripts/python.exe}"
PROG="$1"; shift
NAME="${AUDIT_RUNNER_NAME:-$(basename "$PROG" .py)}"
if [ -n "$AUDIT_STATE_DIR" ]; then STATE="$AUDIT_STATE_DIR"; else STATE="$(cygpath -u "${LOCALAPPDATA:-/tmp}" 2>/dev/null || echo /tmp)/certifield_audit"; fi
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
n=0
while true; do
  "$PY" "$HERE/$PROG" "$@"
  code=$?
  # restart only when a task was in flight: the watchdog (exit 1), or a crash while it
  # dumped tracebacks (seen once as exit 139); a real error exits 2 with no in-flight file
  if [ $code -eq 0 ] || [ ! -f "$STATE/$NAME.inflight.json" ]; then
    exit $code
  fi
  n=$((n + 1))
  echo "[run.sh] restart $n (exit $code) after a timed-out task: $(cat "$STATE/$NAME.inflight.json")"
  if [ $n -gt 20000 ]; then echo "[run.sh] too many restarts"; exit 3; fi
done
