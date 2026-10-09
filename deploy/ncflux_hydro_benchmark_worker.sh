#!/usr/bin/env bash
# Detached single-run worker. Install in a NEW prepared benchmark directory.
# No build/retry/physics changes: never reuse out/ or relaunch this worker.
set -euo pipefail
run_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$run_dir"
exec 9>worker.lock
flock -n 9 || { echo 'This benchmark worker is already running.' >&2; exit 1; }
[[ ! -e launched.flag ]] || { echo 'Already launched; preserve results and prepare a new directory.' >&2; exit 1; }
[[ -d out && -z "$(ls -A out)" ]] || { echo 'Requires a new, empty out directory.' >&2; exit 1; }
test -x bin/drutes
test -r drutes.conf/global.conf
test -r drutes.conf/netcdf/mRM_Fluxes_States.nc
test -r drutes.conf/netcdf/dem.nc
date -u +%FT%TZ > launched.flag
date -u +%FT%TZ > started-utc.txt
printf '%s\n' RUNNING > status.txt
printf '%s\n' "$$" > worker-pid.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export GFORTRAN_UNBUFFERED_ALL=y
model_pid=''
interrupted=0
interrupt_worker() {
  interrupted=1
  [[ -z "$model_pid" ]] || kill -TERM "$model_pid" 2>/dev/null || true
}
trap interrupt_worker INT TERM
bin/drutes > terminal.log 2>&1 &
model_pid=$!
printf '%s\n' "$model_pid" > model-pid.txt
printf 'Model PID=%s; directory=%s\n' "$model_pid" "$run_dir"
rc=0
wait "$model_pid" || rc=$?
if [[ "$interrupted" == 1 ]]; then
  wait "$model_pid" 2>/dev/null || true
  printf '%s\n' INTERRUPTED > status.txt
elif [[ "$rc" == 0 ]]; then
  printf '%s\n' FINISHED > status.txt
elif grep -q 'F I N I S H E D' terminal.log; then
  printf '%s\n' NUMERICALLY_FINISHED_WITH_POSTPROCESS_ERROR > status.txt
else
  printf '%s\n' FAILED > status.txt
fi
printf '%s\n' "$rc" > exit-code.txt
date -u +%FT%TZ > finished-utc.txt
printf 'Model ended; exit=%s; status=%s\n' "$rc" "$(<status.txt)"
exit "$rc"
