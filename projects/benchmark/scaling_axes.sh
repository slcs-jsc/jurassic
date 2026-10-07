# OpenMP scaling of formod function over geometry, channels (ND) and gas sets (NG), one axis
# at a time around the baselines of CASES. Run via run_scaling_axes_<machine>.sh;
# as an array job, task i runs the i-th case into runs/scaling_axes_<array id>/<case>/.
# Channels are ND evenly spaced rows of CHANNEL_FILE (default configs/channels_alt3.tsv),
# gases come from configs/gas_sets/<set>.txt.
#
# MODES (default: strong t1check batches):
#   cost     only the 1-thread reference runs, to measure the cost per setting
#   strong   STRONG_BATCH scenes over STRONG_CURVE threads (up to one socket);
#            settings not in CURVE_SETTINGS (default: all) only get STRONG_THREADS.
#            T1 comes from a 1-thread run of SCENES_PER_THREAD scenes.
#   t1check  1 thread on the full STRONG_BATCH for T1_CHECK (default: baselines)
#   batches  BATCH_THREADS (default: one socket) over BATCH_LIST for BATCH_SETTINGS
#   weak     SCENES_PER_THREAD scenes per thread over THREAD_LIST
# SPREAD_THREADS adds runs with the threads split over all sockets.
# Every run times MAX_ITER batches (default 3) after an untimed warm-up batch.
#
# Overrides: CASES, AXES, MODES, CHANNEL_FILE, CHANNEL_LIST, GAS_SETS, BASE_GAS, THREAD_LIST,
#            CURVE_SETTINGS, STRONG_CURVE, STRONG_THREADS, STRONG_BATCH, SPREAD_THREADS,
#            SCENES_PER_THREAD, MAX_ITER, REPS, T1_CHECK, BATCH_SETTINGS, BATCH_THREADS,
#            BATCH_LIST

set -euo pipefail

JR_EXPERIMENT=scaling_axes

cases=${CASES:-"zenith_baseline nadir_baseline limb_baseline"}
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  read -r -a case_arr <<< "$cases"
  if [ "$SLURM_ARRAY_TASK_ID" -ge "${#case_arr[@]}" ]; then
    echo "Array task $SLURM_ARRAY_TASK_ID: no case left in CASES='$cases', nothing to do."
    exit 0
  fi
  cases=${case_arr[$SLURM_ARRAY_TASK_ID]}
  RUN_ID=${RUN_ID:-scaling_axes_${SLURM_ARRAY_JOB_ID}}/$cases
fi
RUN_ID=${RUN_ID:-scaling_axes_${SLURM_JOB_ID:-manual}}

export JR_SCRIPTS_DIR_OVERRIDE="$jr_scripts_dir"
source "$jr_scripts_dir/base.sh"

CASE_NAME=${cases%% *}
bench_init

axes=${AXES:-"geometry channels gases"}
modes=${MODES:-"strong t1check batches"}
for m in $modes; do
  case "$m" in cost|strong|weak|t1check|batches) ;;
    *) echo "Unknown mode '$m' (expected cost, strong, weak, t1check, batches)" >&2; exit 1 ;;
  esac
done
has_mode() { [[ " $modes " == *" $1 "* ]]; }
MAX_ITER=${MAX_ITER:-3}
reps=${REPS:-2}
k=${SCENES_PER_THREAD:-4}
channel_list=${CHANNEL_LIST:-"8 16 32 64 128"}
gas_sets=${GAS_SETS:-"ng04 ng08 ng13 ng18"}

config_dir="$JR_REPO_ROOT/projects/benchmark/configs"
cases_tsv="$config_dir/baseline_cases.tsv"
channel_file=${CHANNEL_FILE:-$config_dir/channels_alt3.tsv}
P=$JR_PHYS_PER_SOCKET
N_PHYS=$(( P * JR_N_SOCKETS ))

if [ -z "${THREAD_LIST:-}" ]; then
  THREAD_LIST=""
  for (( t=1; t<P; t*=2 )); do THREAD_LIST="$THREAD_LIST $t"; done
  for (( s=1; s<=JR_N_SOCKETS; s++ )); do THREAD_LIST="$THREAD_LIST $(( s * P ))"; done
fi
spread_threads=${SPREAD_THREADS:-}
[ "$JR_N_SOCKETS" -lt 2 ] && spread_threads=""
# strong curve: up to one socket
strong_curve=${STRONG_CURVE:-}
if [ -z "$strong_curve" ]; then
  for (( t=2; t<P; t*=2 )); do strong_curve="$strong_curve $t"; done
  strong_curve="$strong_curve $P"
fi
strong_threads=${STRONG_THREADS-$P}
curve_settings=${CURVE_SETTINGS:-all}
batch_threads=${BATCH_THREADS:-$P}
batch_list=${BATCH_LIST:-"$(( 2 * batch_threads )) $(( 4 * batch_threads )) $(( 8 * batch_threads )) $(( 16 * batch_threads ))"}

# strong batch: lcm of all strong thread counts, >= 2 scenes per thread
if has_mode strong || has_mode t1check; then
  all_strong=""
  for t in $strong_curve $strong_threads $spread_threads; do
    [ "$t" -gt 1 ] && [ "$t" -le "$N_PHYS" ] && all_strong="$all_strong $t"
  done
  if [ -z "${STRONG_BATCH:-}" ]; then
    gcd() { local a=$1 b=$2; while [ "$b" -gt 0 ]; do set -- "$b" $(( a % b )); a=$1; b=$2; done; echo "$a"; }
    lcm=1 tmax=1
    for t in $all_strong; do
      lcm=$(( lcm * t / $(gcd "$lcm" "$t") )); [ "$t" -gt "$tmax" ] && tmax=$t
    done
    STRONG_BATCH=$lcm
    while [ "$STRONG_BATCH" -lt $(( 2 * tmax )) ]; do STRONG_BATCH=$(( STRONG_BATCH + lcm )); done
  fi
  for t in $all_strong; do
    (( STRONG_BATCH % t == 0 )) || { echo "STRONG_BATCH=$STRONG_BATCH is not divisible by $t threads." >&2; exit 1; }
  done
fi

{
  echo "cases=$cases"; echo "axes=$axes"; echo "threads=$THREAD_LIST"; echo "spread_threads=$spread_threads"
  echo "modes=$modes"; echo "scenes_per_thread=$k"; echo "reps=$reps"; echo "max_iter=$MAX_ITER"
  echo "batch_threads=$batch_threads"; echo "batch_list=$batch_list"
  echo "curve_settings=$curve_settings"; echo "strong_curve=$strong_curve"; echo "strong_threads=$strong_threads"; echo "strong_batch=${STRONG_BATCH:-}"
  echo "channel_file=$channel_file"; echo "channel_list=$channel_list"; echo "gas_sets=$gas_sets"
} >> "$JR_RUN_DIR/config.txt"

case_field() { awk -F'\t' -v c="$1" -v f="$2" 'NR > 1 && $1 == c { print $f; exit }' "$cases_tsv"; }
for c in $cases; do
  [ -n "$(case_field "$c" 1)" ] || { echo "Unknown case '$c' (see $cases_tsv)" >&2; exit 1; }
done
# baseline gas set: BASE_GAS or the case's entry in baseline_cases.tsv
base_gas() { echo "${BASE_GAS:-$(case_field "$1" 6)}"; }

# settings: "label geometry nd gas_set"
settings=()
declare -A setting_axes is_baseline
add_setting() {        # axis geometry nd gas_set
  local label="$2_nd$3_$4"
  [ -z "${setting_axes[$label]:-}" ] && settings+=("$label $2 $3 $4")
  setting_axes[$label]="${setting_axes[$label]:+${setting_axes[$label]},}$1"
}
for axis in $axes; do
  case "$axis" in
    geometry)
      for c in $cases; do add_setting geometry "$(case_field "$c" 2)" "$(case_field "$c" 5)" "$(base_gas "$c")"; done ;;
    channels|gases)
      for c in $cases; do
        g=$(case_field "$c" 2); nd=$(case_field "$c" 5)
        if [ "$axis" = channels ]; then
          for x in $channel_list; do add_setting "channels@$g" "$g" "$x" "$(base_gas "$c")"; done
        else
          for gs in $gas_sets; do add_setting "gases@$g" "$g" "$nd" "$gs"; done
        fi
      done ;;
    *) echo "Unknown axis '$axis' (expected geometry, channels, gases)" >&2; exit 1 ;;
  esac
done

for c in $cases; do is_baseline[$(case_field "$c" 2)_nd$(case_field "$c" 5)_$(base_gas "$c")]=1; done
t1_check=${T1_CHECK-${!is_baseline[*]}}
batch_settings=${BATCH_SETTINGS-${!is_baseline[*]}}
{ echo "t1_check=$t1_check"; echo "batch_settings=$batch_settings"; } >> "$JR_RUN_DIR/config.txt"

# one build per (ND, NG)
declare -A bin_dir
build_variant() {
  local nd=$1 ng=$2 key="nd$1_ng$2"
  [ -n "${bin_dir[$key]:-}" ] && return
  echo "=== building ND=$nd NG=$ng ==="
  bench_build_isolated MPI="$JR_MPI" MPICC="$JR_MPICC" COMPILER="$JR_COMPILER" \
    GPU=0 DEFINES="-DND=$nd -DNG=$ng"
  mv "$JR_WORK_DIR/bin" "$JR_WORK_DIR/bin_$key"
  bin_dir[$key]="$JR_WORK_DIR/bin_$key"
}

# per-setting inputs
settings_tsv="$JR_RUN_DIR/settings.tsv"
printf 'setting\taxes\tgeometry\tnd\tng\tgas_set\tnr\tactive_pairs\n' > "$settings_tsv"
declare -A setting_dir setting_bin

for s in "${settings[@]}"; do
  read -r label geom nd gas_set <<< "$s"
  gas_file="$config_dir/gas_sets/$gas_set.txt"
  ng=$(grep -c '^[^#[:space:]]' "$gas_file")
  row=$(awk -F'\t' -v g="$geom" 'NR > 1 && $2 == g { print; exit }' "$cases_tsv")
  ctl_template="$JR_REPO_ROOT/$(cut -f3 <<< "$row")"
  nr=$(cut -f4 <<< "$row")

  build_variant "$nd" "$ng"
  dir="$JR_WORK_DIR/settings/$label"
  mkdir -p "$dir/data"
  awk -v tblbase="$JR_TBLBASE" '{ if ($1 == "TBLBASE") print "TBLBASE = " tblbase; else print $0 }' \
    "$ctl_template" > "$dir/base.ctl"
  active=$(python3 "$JR_SCRIPT_DIR/generate_ctl.py" --channels "$channel_file" --nd "$nd" \
    --gas-file "$gas_file" "$dir/base.ctl" "$dir/run.ctl")
  ( cd "$dir" && "${bin_dir[nd${nd}_ng${ng}]}/climatology" run.ctl data/atm.tab \
      && "${bin_dir[nd${nd}_ng${ng}]}/$geom" run.ctl data/obs.tab ) > "$dir/prepare.log" 2>&1

  setting_dir[$label]=$dir
  setting_bin[$label]="${bin_dir[nd${nd}_ng${ng}]}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$label" "${setting_axes[$label]}" "$geom" "$nd" "$ng" "$gas_set" "$nr" "$active" >> "$settings_tsv"
done
mkdir -p "$JR_RUN_DIR/ctl"
for label in "${!setting_dir[@]}"; do cp "${setting_dir[$label]}/run.ctl" "$JR_RUN_DIR/ctl/$label.ctl"; done

for rep in $(seq 1 "$reps"); do
  for s in "${settings[@]}"; do
    label=${s%% *}
    cd "${setting_dir[$label]}"
    JR_ACTIVE_CTL=run.ctl
    JR_FORMOD_BIN="${setting_bin[$label]}/formod"
    # 1-thread reference
    bench_run_time "${label}_compact" 1 "$k" "$rep" "$(cpus_phys_compact 1)"
    if has_mode weak; then
      for t in $THREAD_LIST; do
        [ "$t" -gt 1 ] && [ "$t" -le "$N_PHYS" ] || continue
        bench_run_time "${label}_compact" "$t" $(( k * t )) "$rep" "$(cpus_phys_compact "$t")"
      done
      for t in $spread_threads; do
        bench_run_time "${label}_spread" "$t" $(( k * t )) "$rep" "$(cpus_phys_split "$t")"
      done
    fi
    if has_mode strong; then
      read -r _ _ s_nd s_gas <<< "$s"
      ns=$strong_threads
      if [ -n "${is_baseline[$label]:-}" ] || [[ " $curve_settings " == *" all "* ]] \
         || [[ " $curve_settings " == *" nd$s_nd "* ]] || [[ " $curve_settings " == *" $s_gas "* ]]; then
        ns=$strong_curve
      fi
      for t in $ns; do
        [ "$t" -gt 1 ] && [ "$t" -le "$N_PHYS" ] || continue
        bench_run_time "${label}_strongcompact" "$t" "$STRONG_BATCH" "$rep" "$(cpus_phys_compact "$t")"
      done
      for t in $spread_threads; do
        bench_run_time "${label}_strongspread" "$t" "$STRONG_BATCH" "$rep" "$(cpus_phys_split "$t")"
      done
    fi
    if has_mode t1check && [[ " $t1_check " == *" $label "* ]]; then
      bench_run_time "${label}_t1full" 1 "$STRONG_BATCH" "$rep" "$(cpus_phys_compact 1)"
    fi
    if has_mode batches && [[ " $batch_settings " == *" $label "* ]]; then
      for b in $batch_list; do
        bench_run_time "${label}_batch" "$batch_threads" "$b" "$rep" "$(cpus_phys_compact "$batch_threads")"
      done
    fi
  done
done
unset JR_FORMOD_BIN

bench_finish
