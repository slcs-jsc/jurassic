# Baseline source snapshot

Unmodified copy of `src/` from `performance-optimization` at `f1fc09a`, i.e.
the forward model before any of the hermes optimizations. It exists so that
`scripts/run_hermes_compare.sh` can profile the baseline and the optimized
code side by side from a single checkout.

Do not edit these files. The Makefile resolves the bundled libraries through
`../libs/build`, which does not exist here, so build with explicit paths:

```bash
make -C projects/benchmark/baseline/src \
  INCDIR="-I $PWD/libs/build/include" LIBDIR="-L $PWD/libs/build/lib" \
  VERSION=baseline
```

`run_hermes_profile.sh` does this itself when called with
`SRC_DIR=projects/benchmark/baseline/src`.
