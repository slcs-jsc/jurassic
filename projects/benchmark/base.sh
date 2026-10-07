# Shared helpers for the benchmark scripts: paths, modules, build, pinned runs.

bench_analyse_topology() {
JR_TOPO_MAP="$JR_RUN_DIR/cpu_topology.csv"
  lscpu -p=CPU,CORE,SOCKET,NODE | grep -v '^#' > "$JR_TOPO_MAP"

  JR_N_LOGICAL=$(wc -l < "$JR_TOPO_MAP")
  JR_N_SOCKETS=$(awk -F, '{print $3}' "$JR_TOPO_MAP" | sort -un | wc -l)
  JR_N_NUMA=$(awk -F, '{print $4}' "$JR_TOPO_MAP" | sort -un | wc -l)
  JR_N_PHYS=$(awk -F, '{print $3"-"$2}' "$JR_TOPO_MAP" | sort -u | wc -l)
  JR_SMT=$(( JR_N_LOGICAL / JR_N_PHYS ))
  JR_PHYS_PER_SOCKET=$(( JR_N_PHYS / JR_N_SOCKETS ))

  {
    echo "logical_cpus=$JR_N_LOGICAL"
    echo "physical_cores=$JR_N_PHYS"
    echo "sockets=$JR_N_SOCKETS"
    echo "numa_domains=$JR_N_NUMA"
    echo "smt_per_core=$JR_SMT"
    echo "phys_cores_per_socket=$JR_PHYS_PER_SOCKET"
  } | tee "$JR_RUN_DIR/topology.txt"
}

find_repo_root() {
    local dir="$1"
    while [ "$dir" != "/" ] && [ -n "$dir" ]; do
        if [ -f "$dir/projects/benchmark/configs/baseline_cases.tsv" ]; then
            printf '%s\n' "$dir"
            return 0
        fi
        dir=$(dirname "$dir")
    done
    return 1
}

bench_init() {

    if [ -n "${JR_SCRIPTS_DIR_OVERRIDE:-}" ] && [ -f "$JR_SCRIPTS_DIR_OVERRIDE/base.sh" ]; then
      JR_SCRIPT_DIR="$JR_SCRIPTS_DIR_OVERRIDE"
    else
        local script_source=${BASH_SOURCE[1]:-$0}
        JR_SCRIPT_DIR=$(cd "$(dirname "$script_source")" && pwd)
    fi

    JR_REPO_ROOT=$(find_repo_root "$JR_SCRIPT_DIR") || {
        echo "Could not locate repo root (no projects/benchmark/configs/baseline_cases.tsv found above $JR_SCRIPT_DIR)" >&2
        exit 1
    }

    if [ ! -f "$JR_REPO_ROOT/projects/benchmark/configs/baseline_cases.tsv" ] \
        && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
        local submit_root
        submit_root=$(cd "$SLURM_SUBMIT_DIR/../../.." && pwd)
        if [ -f "$submit_root/projects/benchmark/configs/baseline_cases.tsv" ]; then
            JR_REPO_ROOT=$submit_root
            JR_SCRIPT_DIR="$JR_REPO_ROOT/projects/benchmark/scripts"
        fi
    fi

    JR_SRC_DIR="$JR_REPO_ROOT/src"
    JR_RUNS_ROOT="$JR_REPO_ROOT/projects/benchmark/runs"
    JR_RUN_ID=${RUN_ID:-${JR_EXPERIMENT:-exp}_${SLURM_JOB_ID:-manual}}
    JR_RUN_DIR="$JR_RUNS_ROOT/$JR_RUN_ID"
    JR_WORK_DIR="$JR_RUN_DIR/work"
    mkdir -p "$JR_RUN_DIR" "$JR_WORK_DIR"

    JR_CASE_NAME=${CASE_NAME:-${GEOMETRY:-zenith}_baseline}
    local cases="$JR_REPO_ROOT/projects/benchmark/configs/baseline_cases.tsv"
    local row
    row=$(awk -F'\t' -v key="$JR_CASE_NAME" \
            'NR > 1 && $1 == key { print; exit }' "$cases")
    if [ -z "$row" ]; then
        echo "Unknown benchmark baseline case: $JR_CASE_NAME" >&2
        exit 1
    fi
    JR_GEOMETRY=$(printf '%s\n' "$row" | awk -F'\t' '{print $2}')
    local ctl_rel
    ctl_rel=$(printf '%s\n' "$row" | awk -F'\t' '{print $3}')
    JR_CTL_TEMPLATE=${CTLFILE:-$JR_REPO_ROOT/$ctl_rel}
    
    case "$JR_GEOMETRY" in
        zenith|nadir|limb) ;;
        *) echo "Unsupported geometry: $JR_GEOMETRY" >&2; exit 1 ;;
    esac
    [ -f "$JR_CTL_TEMPLATE" ] || {
        echo "Control file not found: $JR_CTL_TEMPLATE" >&2; exit 1; }
    
    JR_TBLBASE=${BENCH_TBLBASE:-/p/data1/slmet/model_data/jurassic/tab/tria_1cm/nc_1e-6/tria}
    [ -d "$(dirname "$JR_TBLBASE")" ] || {
        echo "Benchmark LUT directory not found: $(dirname "$JR_TBLBASE")" >&2; exit 1; }
    
    JR_COMPILER=${COMPILER_CPU:-gcc}
    JR_MARCH_NATIVE=${MARCH_NATIVE:-1}
    # CPU target gcc resolves -march=native to (e.g. znver2, skylake-avx512);
    # native binaries only run on that CPU, so caches must be keyed by it
    JR_MARCH_TAG=generic
    if [ "$JR_MARCH_NATIVE" = 1 ]; then
        # awk must read all input: exiting early SIGPIPEs gcc, which fails under pipefail
        JR_MARCH_TAG=$("${JR_COMPILER}" -march=native -Q --help=target 2>/dev/null \
            | awk '$1 == "-march=" && t == "" { t = $2 } END { print t }') || true
        JR_MARCH_TAG=${JR_MARCH_TAG:-native}
    fi
    JR_MPI=${MPI:-0}
    JR_MPICC=${MPICC:-mpicc}

    cd "$JR_WORK_DIR"
    export LANG=C LC_ALL=C

    if command -v ml >/dev/null 2>&1; then
        ml Stages/2026 GCC/14.3.0
        ml ParaStationMPI
        ml CMake/4.0.3
        ml ecBuild
        ml SciPy-bundle/2025.07
        ml netcdf4-python/1.7.2
    fi

    export LD_LIBRARY_PATH="$JR_REPO_ROOT/libs/build/lib:$JR_REPO_ROOT/libs/build/lib64:${LD_LIBRARY_PATH:-}"
    
    # pinning is set per run in bench_run_time
    unset OMP_PLACES OMP_PROC_BIND

    JR_ACTIVE_CTL="$JR_WORK_DIR/${JR_CASE_NAME}.ctl"
    awk -v tblbase="$JR_TBLBASE" \
    '{ if ($1 == "TBLBASE") print "TBLBASE = " tblbase; else print $0; }' \
    "$JR_CTL_TEMPLATE" > "$JR_ACTIVE_CTL"

    command -v lscpu    >/dev/null 2>&1 && lscpu > "$JR_RUN_DIR/lscpu.txt"
    command -v numactl  >/dev/null 2>&1 && numactl --hardware > "$JR_RUN_DIR/numactl.txt"
    ( cd "$JR_REPO_ROOT" && git rev-parse HEAD 2>/dev/null ) \
        > "$JR_RUN_DIR/git_commit.txt" || true
    
    bench_analyse_topology

    {
        echo "experiment=${JR_EXPERIMENT:-unknown}"
        echo "run_id=$JR_RUN_ID"
        echo "case_name=$JR_CASE_NAME"
        echo "geometry=$JR_GEOMETRY"
        echo "ctl_template=$JR_CTL_TEMPLATE"
        echo "tblbase=$JR_TBLBASE"
        echo "compiler=$JR_COMPILER"
        echo "march_native=$JR_MARCH_NATIVE"
        echo "march=$JR_MARCH_TAG"
        echo "partition=${SLURM_JOB_PARTITION:-unknown}"
        echo "nodelist=${SLURM_JOB_NODELIST:-unknown}"
    } > "$JR_RUN_DIR/config.txt"
}

# Build in a private copy of src/
# (so two concurrent jobs on the same repo checkout never race on $JR_SRC_DIR's object files ) 
bench_build_isolated() {
    local build_root="$JR_WORK_DIR/build_src"
    local build_dir="$build_root/src"

    rm -rf "$build_root"
    mkdir -p "$build_root"
    local entry base
    for entry in "$JR_REPO_ROOT"/*; do
        base="$(basename "$entry")"
        if [[ "$base" == "src" ]]; then
            cp -r "$entry" "$build_dir"
        else
            ln -s "$entry" "$build_root/$base"
        fi
    done

    if ! ( cd "$build_dir" && make clean && make -j MARCH_NATIVE="$JR_MARCH_NATIVE" "$@" ) 1>&2; then
        echo "ERROR: isolated build failed -> $build_dir" >&2
        rm -rf "$build_root"
        return 1
    fi

    mkdir -p "$JR_WORK_DIR/bin"
    rm -f "$JR_WORK_DIR/bin"/*
    local bin
    while IFS= read -r -d '' bin; do
        cp "$bin" "$JR_WORK_DIR/bin/"
    done < <(find "$build_dir" -maxdepth 1 -type f -executable -print0)

    rm -rf "$build_root"

    JR_BIN_DIR="$JR_WORK_DIR/bin"
}

# Timed forward run, one thread pinned per CPU of the comma-separated core list.
bench_run_time() {
  local label=$1 threads=$2 batch=$3 rep=$4 cores=$5

  if [ -z "$cores" ]; then
    echo "ERROR: bench_run_time: core list (arg 5) is required" >&2
    return 1
  fi

  local tag="${label}.t${threads}.b${batch}.rep${rep}"
  local txt="$JR_WORK_DIR/out/${tag}.txt"
  local tab="/tmp/jurassic_${JR_RUN_ID//\//_}_${tag}.tab"
  local formod_bin="${JR_FORMOD_BIN:-${JR_BIN_DIR:-$JR_SRC_DIR}/formod}"
  local places
  places=$(sed 's/[0-9]\+/{&}/g' <<< "$cores")
  # MAX_ITER set: time exactly that many batches, time budget out of the way
  local budget=${TIME_BUDGET:-10}
  [ -n "${MAX_ITER:-}" ] && budget=1e9
  mkdir -p "$JR_WORK_DIR/out"

  echo "--- $tag (cores=$cores) ---"

  set +e
  OMP_NUM_THREADS=$threads OMP_PLACES="$places" OMP_PROC_BIND=close \
    JURASSIC_TIME_BUDGET="$budget" JURASSIC_MAX_ITER="${MAX_ITER:-2147483647}" \
    "$formod_bin" "$JR_ACTIVE_CTL" data/obs.tab data/atm.tab "$tab" \
    TASK time BATCH_SIZE "$batch" \
    > "$txt" 2>&1
  local rc=$?
  set -e

  {
    echo "label=$label"
    echo "threads=$threads"
    echo "batch_size=$batch"
    echo "rep=$rep"
    echo "core_list=$cores"
    echo "time_budget=$budget"
    echo "max_iter=${MAX_ITER:-}"
    echo "exit_code=$rc"
  } >> "$txt"

  rm -f "$tab"
  [ "$rc" -ne 0 ] && echo "WARNING: run failed ($tag, exit $rc)" >&2
  return 0
}

# copy results to run directory
bench_finish() {
  cp -a "$JR_WORK_DIR/out" "$JR_RUN_DIR/" 2>/dev/null || true
  cp -a "$JR_ACTIVE_CTL" "$JR_RUN_DIR/" 2>/dev/null || true
  echo
  echo "=== ${JR_EXPERIMENT:-experiment} complete ==="
  echo "Run directory: $JR_RUN_DIR"
  echo "Raw output:    $JR_RUN_DIR/out/<label>.t<N>.b<N>.rep<N>.txt"
}

_cpus_require_count() {
  local wanted=$1 got_list=$2 who=$3
  local got=0
  [ -n "$got_list" ] && got=$(($(grep -o , <<< "$got_list" | wc -l) + 1))
  if [ "$got" -ne "$wanted" ]; then
    echo "FATAL: $who asked for $wanted CPUs but found only $got on this node" \
         "(topology: $JR_TOPO_MAP)." >&2
    exit 1
  fi
  printf '%s' "$got_list"
}

# N physical cores on socket $2 (default socket 0)
cpus_phys() {
  local n=$1 sock=${2:-0}
  local list
  list=$(awk -F, -v s="$sock" '$3 == s { if (!($2 in c)) { c[$2]=1; print $1 } }' "$JR_TOPO_MAP" \
    | sort -n | head -n "$n" | paste -sd,)
  _cpus_require_count "$n" "$list" "cpus_phys $n $sock"
}

# N physical cores, filling socket 0 completely before using socket 1, ...
cpus_phys_compact() {
  local n=$1 left=$1 s take parts=()
  for (( s=0; left>0 && s<JR_N_SOCKETS; s++ )); do
    take=$(( left < JR_PHYS_PER_SOCKET ? left : JR_PHYS_PER_SOCKET ))
    parts+=("$(cpus_phys "$take" "$s")")
    left=$(( left - take ))
  done
  local IFS=,
  _cpus_require_count "$n" "${parts[*]}" "cpus_phys_compact $n"
}

# N physical cores split evenly over all sockets (N must be a multiple of the socket count)
cpus_phys_split() {
  local n=$1 s parts=()
  if (( n % JR_N_SOCKETS != 0 )); then
    echo "FATAL: cpus_phys_split $n: not divisible by $JR_N_SOCKETS sockets." >&2
    exit 1
  fi
  for (( s=0; s<JR_N_SOCKETS; s++ )); do
    parts+=("$(cpus_phys $(( n / JR_N_SOCKETS )) "$s")")
  done
  local IFS=,
  _cpus_require_count "$n" "${parts[*]}" "cpus_phys_split $n"
}
