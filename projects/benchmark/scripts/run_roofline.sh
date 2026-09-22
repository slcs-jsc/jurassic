#!/bin/bash
#SBATCH --account=slmet
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --disable-perfparanoid
#SBATCH --job-name=e1_roofline

# TODO: single core (cores own bandwith) vs node-level roofline (full memory bandwith ceiling with all cores active)
# cache aware roofline (separate L1/L2/L3/DRAM ceilings?)

# Performs LIKWID profiling (FLOPS_DP, MEM_DP) to measure operational intensity 
# i.e. the ratio of computational work to data movement: 
#     operational intensity = DP [MFLOP/s] / (Memory bandwidth [MBytes/s])
#                           = FLOPs per byte moved
# Goal: analyse whether algorithm is memory- or compute-bound. Cache-residency crossover
# Reference: https://github.com/RRZE-HPC/likwid/wiki/Tutorial:-Empirical-Roofline-Model

# Output: out/size_<case>.t1.<GROUP>.b<N>.rep<N>.csv

set -euo pipefail
 
JR_EXPERIMENT=e1_roofline
RUN_ID=${RUN_ID:-e1_roofline_${SLURM_JOB_ID:-manual}}
 
if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/base.sh" ]; then
  jr_scripts_dir="$SLURM_SUBMIT_DIR"
else
  jr_scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "$jr_scripts_dir/base.sh"
 
reps=${REPS:-3}
threads=${THREADS:-48}
batch=${BATCH_SIZE:-240}
groups=${LIKWID_GROUPS:-"FLOPS_DP MEM_DP CACHE"}
 
# Cases to sweep -- defaults to all three geometries.
# Override with e.g. CASE_LIST="zenith_baseline" for a single case.
case_list=${CASE_LIST:-"zenith_baseline nadir_baseline limb_baseline"}
 
first=1
for case_name in $case_list; do
  export CASE_NAME="$case_name"
 
  bench_init
  if [ "$first" -eq 1 ]; then
    bench_build_forward "${VARIANT:-base}"
    bench_validate
    bench_check_groups "$groups"
 
    # Measure roofline ceilings via likwid-bench
    cores="$(cpus_all_phys "$threads")"
    ceilings_file="$JR_RUN_DIR/ceilings.txt"
    bench_ws=${BENCH_WORKING_SET:-2GB}
    flops_bench=${FLOPS_BENCH:-peakflops}  # Set peakflops_avx_fma if needed
    bw_bench=${BW_BENCH:-stream_mem}

    # Build one -w workgroup per socket, splitting threads evenly across sockets. 
    sockets_needed=$(( (threads + JR_PHYS_PER_SOCKET - 1) / JR_PHYS_PER_SOCKET ))
    [ "$sockets_needed" -gt "$JR_N_SOCKETS" ] && sockets_needed=$JR_N_SOCKETS
    remaining=$threads
    workgroup_args=()
    for (( s=0; s<sockets_needed; s++ )); do
      take=$(( remaining < JR_PHYS_PER_SOCKET ? remaining : JR_PHYS_PER_SOCKET ))
      workgroup_args+=( -w "S${s}:${bench_ws}:${take}" )
      remaining=$(( remaining - take ))
    done

    echo "=== measuring compute ceiling ($flops_bench, threads=$threads, workgroups=${workgroup_args[*]}) ==="
    peak_flops=$(likwid-bench -t "$flops_bench" "${workgroup_args[@]}" 2>&1 \
      | tee "$JR_RUN_DIR/likwid_bench_flops.txt" \
      | awk '/MFlops\/s:/ { print $2; exit }')

    echo "=== measuring bandwidth ceiling ($bw_bench, threads=$threads, workgroups=${workgroup_args[*]}) ==="
    stream_bw=$(likwid-bench -t "$bw_bench" "${workgroup_args[@]}" 2>&1 \
      | tee "$JR_RUN_DIR/likwid_bench_bw.txt" \
      | awk '/MByte\/s:/ { print $2; exit }')
 
    echo "peak_flops_mflops=${peak_flops}" | tee "$ceilings_file"
    echo "stream_bw_mbytes=${stream_bw}"   | tee -a "$ceilings_file"
    echo "threads=${threads}"              | tee -a "$ceilings_file"
 
    first=0
  fi
  bench_prepare_inputs
 
  for rep in $(seq 1 "$reps"); do
    for group in "${JR_GROUPS[@]}"; do
      bench_run_forward "$case_name" "$threads" "$group" "$batch" "$rep"
    done
  done
done
 
bench_finish
