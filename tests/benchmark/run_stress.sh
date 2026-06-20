#!/usr/bin/env bash
# Stress-test runner for an extreme mesh. Measures load time, peak RSS and
# per-step throughput. Caps virtual memory so a worst case aborts cleanly instead
# of thrashing swap / triggering the OOM killer on the session.
#   tests/benchmark/run_stress.sh [conf]
set -u
REPO="$(cd "$(dirname "$0")/../.." && pwd)"; cd "$REPO"
CONF="${1:-tests/benchmark/stress_hp_outer.conf}"
BIN="$REPO/MC-DC_Simulator"
OUT="tests/benchmark/output"; mkdir -p "$OUT"

# Safety cap: ~53 GB of virtual memory (box has 62 GB). If the mesh needs more,
# allocation throws and the process aborts cleanly rather than swap-thrashing.
ulimit -v 56000000

echo "== MC-DC STRESS test =="
echo "conf : $CONF"
echo "cores: $(grep -E '^num_process' "$CONF" | awk '{print $2}')   N/T: $(grep -E '^N ' "$CONF" | awk '{print $2}')/$(grep -E '^T ' "$CONF" | awk '{print $2}')"
echo "vmem cap: $(ulimit -v) KB   free now: $(free -h | awk '/Mem:/{print $7}')"
echo "starting (mesh load may take minutes)..."

TIMEFILE="$OUT/_stress_time.txt"
/usr/bin/time -v "$BIN" --conf "$CONF" >"$OUT/_stress_run.log" 2>"$TIMEFILE"
rc=$?

WALL=$(grep -i "wall clock" "$TIMEFILE" | sed 's/.*: //')
RSS=$(grep -i "Maximum resident set size" "$TIMEFILE" | sed 's/.*: //')
RSS_GB=$(awk "BEGIN{printf \"%.2f\", ${RSS:-0}/1048576}")
SIM=$(grep -i "ended after" "$OUT/_stress_run.log" | tail -1 | sed 's/.*ended after: //')
echo "---- result ----"
echo "exit code  : $rc   (137=killed/OOM, 134=abort/bad_alloc)"
echo "wall clock : $WALL   (includes mesh load)"
echo "sim time   : ${SIM:-n/a}   (compute only; load ~= wall - sim)"
echo "peak RSS   : ${RSS:-?} kB (${RSS_GB} GB)"
grep -iE "ERROR|bad_alloc|terminate|what\(\)|cannot|Segmentation" "$OUT/_stress_run.log" "$TIMEFILE" 2>/dev/null | head -5
