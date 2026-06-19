#!/usr/bin/env bash
# Performance benchmark runner for the MC-DC simulator.
# Runs the demanding PorusMedia config under /usr/bin/time and reports wall time,
# peak resident memory, and the DWI result hash (so optimizations can be checked
# for both speed and bit-exact results). Run from the repo root.
#
#   tests/benchmark/run_bench.sh [conf]
#
# Notes:
#  - Wall time / peak RSS are MACHINE-SPECIFIC (a relative before/after reference,
#    not a portable gate). The DWI sha256 IS portable for bit-exact changes.
#  - Peak RSS here is dominated by the mesh + AABB grid (roughly N-independent);
#    time scales with N*T. Tier 0/1 optimizations mainly move TIME.
set -u
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
CONF="${1:-tests/benchmark/porous_media_bench.conf}"
BIN="$REPO/MC-DC_Simulator"
OUT="tests/benchmark/output"
mkdir -p "$OUT"
PREFIX="$OUT/porous_media"

echo "== MC-DC benchmark =="
"$BIN" 2>/dev/null | grep -i version | head -1
echo "conf: $CONF"
echo "git:  $(git rev-parse --short HEAD 2>/dev/null)"
echo "cores: $(grep -E '^num_process' "$CONF" | awk '{print $2}')   N/T: $(grep -E '^N ' "$CONF" | awk '{print $2}')/$(grep -E '^T ' "$CONF" | awk '{print $2}')"
echo "running..."

TIMEFILE="$OUT/_time.txt"
/usr/bin/time -v "$BIN" --conf "$CONF" >"$OUT/_run.log" 2>"$TIMEFILE"
rc=$?
if grep -q "\[ERROR\]" "$OUT/_run.log"; then echo "!! ERROR in run:"; grep "\[ERROR\]" "$OUT/_run.log" | head; fi

WALL=$(grep -i "wall clock" "$TIMEFILE" | sed 's/.*: //')
RSS=$(grep -i "Maximum resident set size" "$TIMEFILE" | sed 's/.*: //')
RSS_MB=$(awk "BEGIN{printf \"%.1f\", ${RSS:-0}/1024}")
SIMSEC=$(grep -i "simulations ended after" "$OUT/_run.log" | tail -1 | sed 's/.*ended after: //')
HASH=$(cat ${PREFIX}_DWI.txt ${PREFIX}_DWI_intra.txt ${PREFIX}_DWI_extra.txt 2>/dev/null | sha256sum | cut -c1-16)

echo "---- result ----"
echo "exit code        : $rc"
echo "wall clock       : $WALL"
echo "sim time (avg)   : $SIMSEC"
echo "peak RSS         : ${RSS} kB (${RSS_MB} MB)"
echo "DWI sha256[0:16] : $HASH"
echo "DWI b0 (tot/in/ex): $(sed -n 1p ${PREFIX}_DWI.txt 2>/dev/null) / $(sed -n 1p ${PREFIX}_DWI_intra.txt 2>/dev/null) / $(sed -n 1p ${PREFIX}_DWI_extra.txt 2>/dev/null)"
