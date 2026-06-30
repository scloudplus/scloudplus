#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
usage: Self_Evaluation/benchmark/measure_instance.sh <build_dir> [bench_scloudplus|all] [seconds] [reps]

Examples:
  cmake -S Implementations/Reference_Implementation/Scloudplus-128/kem \
        -B build/ref-128-sm3 \
        -DSCLOUDPLUS_FAMILY=SM3
  Self_Evaluation/benchmark/measure_instance.sh build/ref-128-sm3 all 0.25 3

The build directory must already be configured from a concrete
Implementations/*_Implementation/Scloudplus-*/kem selection.
The selector `all` is accepted as an alias for the leaf-local
`bench_scloudplus` target.
USAGE
}

if [[ $# -lt 1 || $# -gt 4 ]]; then
  usage
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="$1"
SELECTOR="${2:-all}"
SECONDS_PER_REP="${3:-1}"
REPS="${4:-5}"

if [[ "$BUILD_DIR" != /* ]]; then
  BUILD_DIR="$PACKAGE_ROOT/$BUILD_DIR"
fi

if [[ ! -d "$BUILD_DIR" ]]; then
  echo "error: build directory not found: $BUILD_DIR" >&2
  exit 3
fi
BUILD_DIR="$(cd "$BUILD_DIR" && pwd)"
if ! awk -v v="$SECONDS_PER_REP" 'BEGIN { exit !(v ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/ && v + 0 > 0) }'; then
  echo "error: seconds must be a positive decimal number" >&2
  exit 4
fi
case "$REPS" in
  ''|*[!0-9]*) echo "error: reps must be a positive integer" >&2; exit 4 ;;
esac
if [[ "$REPS" -le 0 ]]; then
  echo "error: reps must be a positive integer" >&2
  exit 4
fi
if ! command -v cmake >/dev/null 2>&1; then
  echo "error: cmake is required" >&2
  exit 7
fi

RESULT_ROOT="$PACKAGE_ROOT/Self_Evaluation/results"
STAMP="$(date -u +%Y%m%d_%H%M%S)"
OUTDIR="$RESULT_ROOT/bench_${STAMP}"
CSV="$OUTDIR/benchmark_${SECONDS_PER_REP}s_${REPS}reps_summary.csv"
mkdir -p "$OUTDIR"

TIME_TOOL=()
if [[ -x /usr/bin/time ]] && /usr/bin/time -v true >/dev/null 2>&1; then
  TIME_TOOL=(/usr/bin/time -v)
fi

collect_environment() {
  {
    echo "timestamp_utc_start=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "package_root=$PACKAGE_ROOT"
    echo "build_dir=$BUILD_DIR"
    echo "selector=$SELECTOR"
    echo "seconds_per_repetition=$SECONDS_PER_REP"
    echo "repetitions=$REPS"
    echo "time_tool=${TIME_TOOL[*]:-unavailable}"
    echo "size_tool=$(command -v size 2>/dev/null || echo unavailable)"
    echo "uname=$(uname -a)"
    if command -v cc >/dev/null 2>&1; then
      echo "cc=$(cc --version | head -n 1)"
    else
      echo "cc=unavailable"
    fi
    if command -v cmake >/dev/null 2>&1; then
      echo "cmake=$(cmake --version | head -n 1)"
    fi
    if command -v lscpu >/dev/null 2>&1; then
      lscpu | sed 's/^/lscpu: /'
    fi
  } > "$OUTDIR/environment.txt"
}

collect_build_info() {
  {
    echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "package_root=$PACKAGE_ROOT"
    echo "build_dir=$BUILD_DIR"
    echo "selector=$SELECTOR"
    if [[ -f "$BUILD_DIR/CMakeCache.txt" ]]; then
      echo "-- CMake cache summary --"
      rg '^(CMAKE_C_COMPILER:|CMAKE_C_FLAGS:|CMAKE_BUILD_TYPE:|CMAKE_GENERATOR:|SCLOUDPLUS_)' \
        "$BUILD_DIR/CMakeCache.txt" || true
    else
      echo "cmake_cache=unavailable"
    fi
  } > "$OUTDIR/build_info.txt"
}

build_target() {
  local target="$1"
  cmake --build "$BUILD_DIR" --target "$target" > "$OUTDIR/${target}.build.log" 2>&1
}

kv_colon() {
  local file="$1"
  local key="$2"
  awk -F': ' -v key="$key" '$1 == key { print $2; found=1; exit } END { if (!found) print "NA" }' "$file"
}

static_size_bytes() {
  local file="$1"
  awk 'NR == 2 { print $4; found=1; exit } END { if (!found) print "NA" }' "$file"
}

peak_rss_bytes() {
  local file="$1"
  awk -F: '/Maximum resident set size/ { gsub(/^[ \t]+/, "", $2); print $2 * 1024; found=1; exit } END { if (!found) print "NA" }' "$file"
}

append_csv_header() {
  cat > "$CSV" <<'HEADER'
target,algorithm_instance,operation,repetitions,total_iterations,avg_time_us,avg_cycles,pk_bytes,sk_bytes,ct_bytes,ss_bytes,static_bytes,peak_rss_bytes,output_file
HEADER
}

append_csv_rows() {
  local name="$1"
  local out="$OUTDIR/${name}.out"
  local timef="$OUTDIR/${name}.time"
  local sizef="$OUTDIR/${name}.size"
  local alg pk sk ct ss static_bytes rss out_base

  alg="$(kv_colon "$out" Algorithm)"
  pk="$(kv_colon "$out" pk_bytes)"
  sk="$(kv_colon "$out" sk_bytes)"
  ct="$(kv_colon "$out" ct_bytes)"
  ss="$(kv_colon "$out" ss_bytes)"
  static_bytes="$(static_size_bytes "$sizef")"
  rss="$(peak_rss_bytes "$timef")"
  out_base="$(basename "$out")"

  awk -F, \
    -v target="$name" \
    -v alg="$alg" \
    -v pk="$pk" \
    -v sk="$sk" \
    -v ct="$ct" \
    -v ss="$ss" \
    -v static_bytes="$static_bytes" \
    -v rss="$rss" \
    -v out_base="$out_base" \
    '/^[^,]+,[0-9]+,[0-9]+,/ {
      print target "," alg "," $1 "," $2 "," $3 "," $4 "," $5 "," pk "," sk "," ct "," ss "," static_bytes "," rss "," out_base
    }' "$out" >> "$CSV"
}

run_one() {
  local exe="$1"
  local name
  name="$(basename "$exe")"

  echo "running $name" | tee -a "$OUTDIR/run.log"
  build_target "$name"

  if command -v size >/dev/null 2>&1; then
    size "$exe" > "$OUTDIR/${name}.size"
  else
    echo "size_unavailable" > "$OUTDIR/${name}.size"
  fi

  if [[ ${#TIME_TOOL[@]} -gt 0 ]]; then
    "${TIME_TOOL[@]}" "$exe" --seconds "$SECONDS_PER_REP" --reps "$REPS" \
      > "$OUTDIR/${name}.out" \
      2> "$OUTDIR/${name}.time"
  else
    "$exe" --seconds "$SECONDS_PER_REP" --reps "$REPS" \
      > "$OUTDIR/${name}.out" \
      2> "$OUTDIR/${name}.time"
  fi

  append_csv_rows "$name"
  {
    echo "target=$name"
    echo "output=$OUTDIR/${name}.out"
    echo "time=$OUTDIR/${name}.time"
    echo "size=$OUTDIR/${name}.size"
  } >> "$OUTDIR/summary.txt"
}

collect_environment
collect_build_info
: > "$OUTDIR/run.log"
: > "$OUTDIR/summary.txt"
append_csv_header

case "$SELECTOR" in
  all|bench_scloudplus)
    target=bench_scloudplus
    ;;
  *)
    target="$SELECTOR"
    if [[ "$target" != bench_* ]]; then
      target="bench_${target}"
    fi
    ;;
esac

build_target "$target"
exe="$BUILD_DIR/$target"
if [[ ! -x "$exe" ]]; then
  echo "error: benchmark executable not found after build: $exe" >&2
  exit 6
fi
run_one "$exe"

echo "timestamp_utc_end=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$OUTDIR/environment.txt"
echo "results_dir=$OUTDIR" | tee -a "$OUTDIR/summary.txt"

base_outdir="$(basename "$OUTDIR")"
echo "$base_outdir" > "$RESULT_ROOT/latest_benchmark_dir.txt"
