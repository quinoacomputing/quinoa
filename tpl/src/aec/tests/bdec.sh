#!/bin/sh
set -e
AEC=../src/aec

if [ ! -f bench.rz ]; then
    echo "No encoded file found. Encoding now..."
    path=$(echo $0 | sed -e 's:[^/]*$::')
    "${path}"/benc.sh
fi
rm -f dec.dat
bsize=$(stat -c "%s" bench.dat)
utime=$(../src/utime $AEC -d -n16 -j64 -r256 -m bench.rz dec.dat 2>&1)
perf=$(echo "$bsize/1048576/$utime" | bc)
echo "[0;32m*** Decoding with $perf MiB/s user time ***[0m"
cmp bench.dat dec.dat
