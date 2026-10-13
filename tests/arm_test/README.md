# ARM forward-model references

`make arm_test` runs `tests/formod_test/run.sh` with the 16 reference files in
`data.ref/`. It covers limb, nadir, and zenith geometry, radiances, CGA, limb
field of view, and kernels. The regular `make check` still uses the default
references in `tests/formod_test/data.ref/`.

These files were generated on a JUPITER ARM login node (`aarch64`, GCC 14.3.0)
from JURASSIC commit `1fed3c2f2cd92ad288da9b19683c7097b912165a` on
2026-10-13. The same commands and inputs as `formod_test` produced them. The
comparison is exact, so intentional changes to forward-model output require
reviewing and updating both reference sets as appropriate.
