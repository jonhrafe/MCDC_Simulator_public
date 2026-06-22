#!/usr/bin/env bash
# =============================================================================
#  GPU stress test -- full tissue model (~10.9k PLY strands, ~49M triangles).
#
#  Validates the GPU PLY pipeline end to end:
#    1. PLY-LIST loading (ply_file_list) of thousands of individual meshes
#    2. crossing sentinel OFF (deportation 0) -- the strands are not watertight
#    3. merged-mesh GPU simulation (one device AABB grid over all triangles) and a
#       signal sanity check with enough particles.
#
#  The mesh paths are machine-specific, so the .list files are generated here (not
#  committed). The voxel was set from the global vertex bbox of all PLYs; walkers
#  spawn in a central dense cube (see strands_stress.conf).
#
#  Usage: tests/benchmark/run_strands_stress.sh [MESH_DIR]
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/../.."                       # repo root
BIN=./MC-DC_Simulator
MESH_DIR="${1:-/home/jonathan/Documents/EPFL/Main_Studies/M_DISCO/decimated_holeless/decimated2}"
SCALE=0.001                                      # PLY file units (um) -> internal mm

[ -x "$BIN" ] || { echo "build $BIN first"; exit 1; }
ls "$MESH_DIR"/*.ply >/dev/null 2>&1 || { echo "no PLYs in $MESH_DIR"; exit 1; }

mkdir -p tests/benchmark/output
echo "Generating PLY lists from $MESH_DIR ..."
ls -1 "$MESH_DIR"/*.ply | awk -v s="$SCALE" '{print $0" "s}' > tests/benchmark/strands.list
head -10 tests/benchmark/strands.list                                  > tests/benchmark/strands_v10.list
echo "  full list: $(wc -l < tests/benchmark/strands.list) PLYs"

echo; echo "=== [1/2] merge validation: 10 strands, GPU vs CPU (within MC noise) ==="
python3 tests/test_gpu_restricted.py tests/benchmark/strands_v10.conf || true

echo; echo "=== [2/2] full stress test: all strands on GPU ==="
/usr/bin/time -v "$BIN" --conf tests/benchmark/strands_stress.conf 2>&1 \
    | grep -iE "GPU path:|triangles|wall clock|Maximum resident|irregular"

f=$(ls -t tests/benchmark/output/strands_stress*_DWI.txt 2>/dev/null | grep -vE "_intra|_extra" | head -1)
N=$(awk '/^N /{print $2; exit}' tests/benchmark/strands_stress.conf)
awk -v N="$N" '{s+=$1;n++; if(NR==1||$1<mn)mn=$1; if($1>mx)mx=$1}
    END{printf "signal: %d measurements, S/N range %.3f .. %.3f (mean %.3f)\n", n, mn/N, mx/N, s/n/N}' "$f"
