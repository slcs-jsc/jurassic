# Shared logic for the likwid benchmark profiling runner scripts:
# Resolve paths, select benchmark case, load modules, set env

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

# N physical cores on socket 0 (no SMT)
cpus_phys() {
  awk -F, '$3 == 0 { if (!($2 in s)) { s[$2]=1; print $1 } }' "$JR_TOPO_MAP" \
    | sort -n | head -n "$1" | paste -sd,
}

# All N physical cores across all sockets
cpus_all_phys() {
  awk -F, '{ k=$3"-"$2; if (!(k in s)) { s[k]=1; print $1 } }' "$JR_TOPO_MAP" \
    | sort -n | head -n "$1" | paste -sd,
}

# N logical CPUs on socket 0: physical cores first, then SMT
cpus_smt() {
  awk -F, '$3 == 0 {print $1}' "$JR_TOPO_MAP" | sort -n | head -n "$1" | paste -sd,
}

# N physical cores spread across all sockets
cpus_spread() {
  awk -F, '{ k=$3"-"$2; if (!(k in s)) { s[k]=1; i[$3]++; print i[$3], $3, $1 } }' \
    "$JR_TOPO_MAP" | sort -k1,1n -k2,2n | awk '{print $3}' \
    | head -n "$1" | paste -sd,
}

cpus_all_phys_socket_qualified() {
  local n=$1
  local remaining=$n
  local parts=()
  local sock
  for sock in $(awk -F, '{print $3}' "$JR_TOPO_MAP" | sort -un); do
    [ "$remaining" -le 0 ] && break
 
    local sock_phys_count
    sock_phys_count=$(awk -F, -v s="$sock" \
      '$3 == s { if (!($2 in seen)) { seen[$2]=1; c++ } } END { print c+0 }' \
      "$JR_TOPO_MAP")
 
    local take=$remaining
    [ "$take" -gt "$sock_phys_count" ] && take=$sock_phys_count
    [ "$take" -le 0 ] && continue

    local idx_list
    idx_list=$(seq 0 $((take - 1)) | paste -sd,)
 
    parts+=("S${sock}:${idx_list}")
    remaining=$(( remaining - take ))
  done
  local IFS='@'
  echo "${parts[*]}"
}


bench_init() {
    local script_source=${BASH_SOURCE[1]:-$0}
    JR_SCRIPT_DIR=$(cd "$(dirname "$script_source")" && pwd)
    JR_REPO_ROOT=$(cd "$JR_SCRIPT_DIR/../../.." && pwd)

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
    JR_MPI=${MPI:-0}
    JR_MPICC=${MPICC:-mpicc}

    cd "$JR_WORK_DIR"
    export LANG=C LC_ALL=C
    
    if command -v ml >/dev/null 2>&1; then
        ml Stages/2026 GCC/14.3.0
        ml ParaStationMPI
        ml likwid/5.4.1
        ml CMake/4.0.3
        ml ecBuild
        ml SciPy-bundle/2025.07
        ml netcdf4-python/1.7.2
    fi
    
    command -v likwid-perfctr >/dev/null 2>&1 || {
        echo "likwid-perfctr not found after 'ml likwid'. Try: ml spider likwid" >&2
        exit 1; }
    
    export LD_LIBRARY_PATH="$JR_REPO_ROOT/libs/build/lib:$JR_REPO_ROOT/libs/build/lib64:${LD_LIBRARY_PATH:-}"
    
    # LIKWID pins via -C; a competing OpenMP affinity policy corrupts the mapping.
    # If LIKWID is not used: set export OMP_PLACES=cores, export OMP_PROC_BIND=close
    unset OMP_PLACES OMP_PROC_BIND

    JR_ACTIVE_CTL="$JR_WORK_DIR/${JR_CASE_NAME}.ctl"
    awk -v tblbase="$JR_TBLBASE" \
    '{ if ($1 == "TBLBASE") print "TBLBASE = " tblbase; else print $0; }' \
    "$JR_CTL_TEMPLATE" > "$JR_ACTIVE_CTL"

    command -v lscpu    >/dev/null 2>&1 && lscpu > "$JR_RUN_DIR/lscpu.txt"
    command -v numactl  >/dev/null 2>&1 && numactl --hardware > "$JR_RUN_DIR/numactl.txt"
    likwid-perfctr -a > "$JR_RUN_DIR/likwid_available_groups.txt" 2>&1 || true
    echo "perf_event_paranoid: $(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo unavailable)" \
        > "$JR_RUN_DIR/perf_paranoid_status.txt"
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
        echo "partition=${SLURM_JOB_PARTITION:-unknown}"
        echo "nodelist=${SLURM_JOB_NODELIST:-unknown}"
    } > "$JR_RUN_DIR/config.txt"
}

# build CPU binaries
# nomemset: exclude  memset(los, 0, sizeof(*los)) in jurassic.c
bench_build_forward() {
    local variant=$1
    local extra=""
    case "$variant" in
        base)     extra="" ;;
        flat_array)  extra="FLAT_ARRAYS=1" ;;                        
        *) echo "Unknown build variant: $variant" >&2; return 1 ;;
    esac
    
    echo "=== building variant '$variant' (EXTRA_CFLAGS='$extra') ==="

    if [ -n "$extra" ]; then
    ( cd "$JR_SRC_DIR" && make clean && make -j MPI="$JR_MPI" MPICC="$JR_MPICC" \
        COMPILER="$JR_COMPILER" GPU=0 LIKWID=1 $extra ) || return 1
    else
        ( cd "$JR_SRC_DIR" && make clean && make -j MPI="$JR_MPI" MPICC="$JR_MPICC" \
            COMPILER="$JR_COMPILER" GPU=0 LIKWID=1 ) || return 1
    fi

    echo "$variant" > "$JR_RUN_DIR/build_variant.txt"
    cd "$JR_WORK_DIR"
}

bench_build_retrieval() {
  local variant=$1
    local extra=""
    case "$variant" in
        base)     extra="" ;;
        flat_array)  extra="FLAT_ARRAYS=1" ;;
        *) echo "Unknown build variant: $variant" >&2; return 1 ;;
    esac

    echo "=== building retrieval variant '$variant' (EXTRA_CFLAGS='$extra') ==="

    if [ -n "$extra" ]; then
    ( cd "$JR_SRC_DIR" && make clean && make -j MPI=1 MPICC="$JR_MPICC" \
        COMPILER="$JR_COMPILER" GPU=0 LIKWID=1 $extra ) || return 1
    else
        ( cd "$JR_SRC_DIR" && make clean && make -j MPI=1 MPICC="$JR_MPICC" \
            COMPILER="$JR_COMPILER" GPU=0 LIKWID=1 ) || return 1
    fi

    echo "$variant" > "$JR_RUN_DIR/build_variant.txt"
    cd "$JR_WORK_DIR"
}

# generate data/atm.tab and data/obs.tab
bench_prepare_inputs() {
  cd "$JR_WORK_DIR"
  rm -rf data && mkdir -p data
  "$JR_SRC_DIR/climatology" "$JR_ACTIVE_CTL" data/atm.tab
  "$JR_SRC_DIR/$JR_GEOMETRY" "$JR_ACTIVE_CTL" data/obs.tab
  cp -a data "$JR_RUN_DIR/data" 2>/dev/null || true
}

# run validation
bench_validate() {
  JR_VALIDATION_OK=1
  if [ "${SKIP_VALIDATION:-0}" = "1" ]; then
    echo "exit_code=skipped" > "$JR_RUN_DIR/validation_status.txt"
    return 0
  fi
 
  set +e
  ( cd "$JR_REPO_ROOT/projects/validation" \
    && VALIDATION_TBLBASE="$JR_TBLBASE" scripts/run_validation.py \
       > "$JR_RUN_DIR/validation.log" 2>&1 )
  local rc=$?
  set -e
 
  echo "exit_code=$rc" > "$JR_RUN_DIR/validation_status.txt"
  local latest
  latest=$(ls -td "$JR_REPO_ROOT/projects/validation/runs/validation_"*/ 2>/dev/null | head -n1 || true)
  [ -n "$latest" ] && cp -a "${latest}summary.tsv" "$JR_RUN_DIR/validation_summary.tsv" 2>/dev/null || true
 
  if [ "$rc" -ne 0 ]; then
    echo "Validation FAILED (exit $rc) -- results are not trustworthy." >&2
    JR_VALIDATION_OK=0
  fi
  return 0
}

# filter requested LIKWID groups against the groups available on this node
bench_check_groups() {
  JR_GROUPS=()
  local g
  for g in $1; do
    if grep -q -w "$g" "$JR_RUN_DIR/likwid_available_groups.txt"; then
      JR_GROUPS+=("$g")
    else
      echo "WARNING: LIKWID group '$g' unavailable on this node -> skipping." >&2
    fi
  done
  if [ ${#JR_GROUPS[@]} -eq 0 ]; then
    echo "No requested LIKWID groups available. See likwid_available_groups.txt" >&2
    return 1
  fi
  return 0
}

bench_run_forward() {
  local label=$1 threads=$2 group=$3 batch=$4 rep=$5
  local cores=${6:-$(cpus_all_phys_socket_qualified "$threads")}
  local tag="${label}.t${threads}.${group}.b${batch}.rep${rep}"
  local csv="$JR_WORK_DIR/out/${tag}.csv"
  local txt="$JR_WORK_DIR/out/${tag}.txt"
  local tab="/tmp/jurassic_${JR_RUN_ID}_${tag}.tab"
  mkdir -p "$JR_WORK_DIR/out"

  # Optional NUMA memory policy independent of LIKWID's CPU pinning above. 
  # NUMA_POLICY=interleave=all 
  local numa_wrap=()
  if [ -n "${NUMA_POLICY:-}" ] && command -v numactl >/dev/null 2>&1; then
    numa_wrap=(numactl "--${NUMA_POLICY}")
  fi

  echo "--- $tag (cores=$cores${NUMA_POLICY:+, numa=$NUMA_POLICY}) ---"

  set +e
  OMP_NUM_THREADS=$threads likwid-perfctr -C "$cores" -g "$group" -m \
    -o "$csv" \
    "${numa_wrap[@]}" "$JR_SRC_DIR/formod" "$JR_ACTIVE_CTL" data/obs.tab data/atm.tab "$tab" \
    JURASSIC_TIME_BUDGET=60 TASK time BATCH_SIZE "$batch" \
    > "$txt" 2>&1
  local rc=$?
  set -e

  {
    echo "label=$label"
    echo "threads=$threads"
    echo "group=$group"
    echo "batch_size=$batch"
    echo "rep=$rep"
    echo "core_list=$cores"
    echo "exit_code=$rc"
  } >> "$txt"
 
  rm -f "$tab"
  [ "$rc" -ne 0 ] && echo "WARNING: run failed ($tag, exit $rc)" >&2
  return 0
}

  
bench_run_retrieval() {
  local label=$1 ranks=$2 threads=$3 group=$4 rep=$5 cores=$6
  # core setup needs to be explicitly defined in order to set up MPI ranks
  if [ -z "$cores" ]; then
    echo "bench_run_retrieval: core list (arg 6) is required" >&2
    return 1
  fi
  local tag="${label}.r${ranks}.t${threads}.${group}.rep${rep}"
  local csv="$JR_WORK_DIR/out/${tag}.csv"
  local txt="$JR_WORK_DIR/out/${tag}.txt"
  mkdir -p "$JR_WORK_DIR/out"
 
  echo "--- $tag (cores=$cores) ---"

  set +e
  OMP_NUM_THREADS=$threads mpirun -np "$ranks" \
    likwid-perfctr -C "$cores" -g "$group" -m -o  "$csv" \
    "$JR_SRC_DIR/retrieval" "$JR_ACTIVE_CTL" "$JR_RET_DIRLIST" \
    > "$txt" 2>&1
  local rc=$?
  set -e
 
  {
    echo "label=$label"
    echo "ranks=$ranks"
    echo "threads=$threads"
    echo "group=$group"
    echo "rep=$rep"
    echo "core_list=$cores"
    echo "exit_code=$rc"
  } >> "$txt"

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
  echo "Raw output:    $JR_RUN_DIR/out/<label>.t<N>.<GROUP>.b<N>.rep<N>.{csv,txt}"
}
 