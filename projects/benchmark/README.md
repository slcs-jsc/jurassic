# CPU/OpenMP Benchmark

This directory contains the OpenMP scaling experiment for JURASSIC's forward model (`formod`).

## Scaling over geometry, channels and gas sets (`scaling_axes.sh`)

The experiment varies one axis at a time: geometry (limb, nadir, zenith), number of
channels (ND) and gas set (NG). The other two axes stay at a reference setting of
32 channels and all 18 gases (`ng18`), defined per geometry in `configs/baseline_cases.tsv`.
This setting lies in the middle of the channel range (8–128) and uses the full gas list.
The reference setting is timed on 1 to 64 threads (one socket) and is the only setting
used for the batch-size sweep and the T1 check; all other settings are timed on 1, 4, 16
and 64 threads.

Each setting is built with its own ND/NG and timed with `formod TASK time`.

### Channels and gases

JURASSIC only computes a gas in a channel if a lookup table exists for that gas at the channel's wavenumber. The runtime therefore depends on the number of active (channel, gas) pairs, i.e. the pairs for which a lookup table exists, and not simply on ND × NG.

The file `configs/channels_alt3.tsv` contains a list of 128 channels between 587 and 739 cm⁻¹, each listed with the gases for which a lookup table exists. 

![Active gases per channel in channels_alt3.tsv](figures/channels_alt3_gas_matrix.png)

Blue: lookup table exists (number of channels listed in brackets). 
Orange: the 32 channels of the reference setting. (`python3 plot_channels.py`).

The script picks ND channels spread evenly over the list, including the first and the
last one. The selected channels therefore always cover the whole range from 587 to
739 cm⁻¹, no matter how many are picked. The reference setting has 458 active pairs.

The gas sets in `configs/gas_sets/` are nested (`ng04` ⊂ `ng08` ⊂ `ng13` ⊂ `ng18`).
`ng04` to `ng13` contain only gases active in all 128 channels, so active pairs = NG × ND;
`ng18` adds the 5 partly active gases (+42 pairs at ND = 32). 
The evaluation reports the expected active pairs and the number of tables formod read.

### Modes

* `strong` (default): fixed batch size, 1 thread up to one socket.
* `t1check`: 1 thread on the full batch, to check the extrapolated T1.
* `batches` (default): one socket's threads over several batch sizes.
* `weak`: fixed scenes per thread.

Threads are pinned with `OMP_PLACES`/`OMP_PROC_BIND=close` on physical cores of socket 0.
Overrides are listed at the top of `scaling_axes.sh`.

```sh
cd projects/benchmark
sbatch run_scaling_axes_jureca.sh          # or run_scaling_axes_juwels.sh
python3 eval_scaling_axes.py runs/scaling_axes_<array job id>
```

The evaluation reports runtime, speedup and efficiency for the timed `formod_batch` call and for the whole application (adding table read, reference run and output), plus single-thread cost and efficiency vs batch size.
