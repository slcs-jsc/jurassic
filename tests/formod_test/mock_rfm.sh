#!/bin/bash

set -euo pipefail

test -f rfm.drv
test -f rfm.atm
IFS= read -r input_line
test -z "$input_line"

read -r nu0 nu1 < <(awk '$1 == "*SPC" { getline; print $1, $2; exit }' rfm.drv)
printf '%s %s\n' "$PWD" "$nu0" >> "$MOCK_RECORD"
printf 'mock RFM log\n' > rfm.log

if grep -q '^FAIL$' rfm.drv; then
    exit 7
fi

if grep -q '^STALE$' rfm.drv && \
   awk -v nu="$nu0" 'BEGIN { exit !(nu > 800) }'; then
    exit 0
fi

if awk -v nu="$nu0" 'BEGIN { exit !(nu > 800) }'; then
    rad=200000
else
    rad=100000
fi

awk '$1 == "*TAN" { active=1; next } active && /^\*/ { exit }
     active && NF { printf "%05d\n", int($1 * 1000) }' rfm.drv | \
while read -r height; do
    for index in "$((10#$height))" "$((10#$height + 1))"; do
        printf -v label '%05d' "$index"
        for kind in rad tra; do
            value=$rad
            if [ "$kind" = tra ]; then value=0.5; fi
            awk -v lo="$nu0" -v hi="$nu1" -v value="$value" \
                'BEGIN { print "#\n#\n#";
                         printf "3 %.4f %.4f %.4f\n", lo, (hi-lo)/2, hi;
                         print value, value, value }' > "${kind}_${label}.asc"
        done
    done
done
