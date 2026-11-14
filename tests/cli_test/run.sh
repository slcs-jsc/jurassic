#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C
set -euo pipefail

# Setup...
jurassic=../../src

# Create directory...
rm -rf data && mkdir -p data

# Test common CLI behavior for every executable listed in src/Makefile.
# Keeping the binary list in the Makefile avoids a second list here that
# would need manual updates whenever a tool is added or removed.
error=0
while IFS= read -r tool ; do
    [ -n "$tool" ] || continue
    echo "checking $tool..."

    # Calling a tool without arguments must fail and print the standard
    # "missing arguments" diagnostic used throughout the command-line tools.
    noarg_out="data/${tool}_noarg.txt"
    if "$jurassic/$tool" > "$noarg_out" 2>&1 ; then
        echo "Missing-argument command succeeded: $tool"
        error=1
    else
        if ! grep -q "Missing or invalid command-line arguments." "$noarg_out"; then
            echo "Missing argument error text: $tool"
            error=1
        fi
    fi

    # Both documented help flags must succeed and print a usage section.
    for flag in -h --help ; do
        echo "  help flag $flag"
        outfile="data/${tool}_${flag#--}.txt"
        if ! "$jurassic/$tool" "$flag" > "$outfile" ; then
            echo "Help command failed: $tool $flag"
            error=1
            continue
        fi
        if ! grep -q "Usage:" "$outfile" ; then
            echo "Missing usage text: $tool $flag"
            error=1
        fi

        # Help should be discoverable even when users add extra arguments.
        # This keeps help behavior forgiving for quick terminal exploration.
        extra_out="data/${tool}_${flag#--}_extra.txt"
        if ! "$jurassic/$tool" "$flag" extra-arg > "$extra_out" ; then
            echo "Help command with extra arguments failed: $tool $flag"
            error=1
            continue
        fi
        if ! grep -q "Usage:" "$extra_out" ; then
            echo "Missing usage text with extra arguments: $tool $flag"
            error=1
        fi
    done

# Read the executable list from the Makefile so this test follows the same
# source of truth as the build and installation rules.
done < <(make --no-print-directory -C ../../src print-exc)

# Spectral windows must identify an existing integer window for each channel.
check_window_error() {
    local name=$1 expected=$2
    shift 2
    if "$jurassic/filter" - "data/${name}.filt" "$@" > "data/${name}.log" 2>&1; then
        echo "Invalid window configuration was accepted: $name"
        error=1
    elif ! grep -Fq "$expected" "data/${name}.log"; then
        echo "Missing window error: $name"
        error=1
    fi
}
check_window_error zero_windows 'NW > 0' ND 1 'NU[0]' 700 NW 0
check_window_error negative_window 'WINDOW[0]' ND 1 'NU[0]' 700 NW 1 'WINDOW[0]' -1
check_window_error high_window 'WINDOW[0]' ND 1 'NU[0]' 700 NW 1 'WINDOW[0]' 1
check_window_error fractional_window 'WINDOW[0]' ND 1 'NU[0]' 700 NW 2 'WINDOW[0]' 1.5
"$jurassic/filter" - data/no_windows.filt NW 0 > data/no_windows.log 2>&1
"$jurassic/filter" - data/last_window.filt ND 1 'NU[0]' 700 NW 2 \
    'WINDOW[0]' 1 > data/last_window.log 2>&1

# A one-point filter cannot provide an interpolation interval.
echo '700 1' > data/one_point.filt
if "$jurassic/tblgen" 1000 285 - data/one_point.filt > data/one_point.log 2>&1; then
    echo 'One-point filter was accepted'
    error=1
elif ! grep -q 'Filter function needs at least two points!' data/one_point.log; then
    echo 'Missing one-point filter error'
    error=1
fi

exit $error
