#!/bin/bash

set -euo pipefail

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}

test -d data.ref
cd ../formod_test
JURASSIC_REF_DIR=../arm_test/data.ref bash ./run.sh
