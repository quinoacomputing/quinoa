#!/bin/sh
set -e
AEC=../src/aec

if [ ! -f  typical.dat ]; then
    rm -f typical.rz
    wget https://www.dkrz.de/redmine/attachments/download/442/typical.rz
    $AEC -d -n16 -j64 -r256 -m typical.rz typical.dat
    rm -f bench.dat
fi
if [ ! -f  bench.dat ]; then
    for i in $(seq 0 499);
    do
        cat typical.dat >> bench.dat
    done
fi
rm -f bench.rz
utime=$(../src/utime $AEC -n16 -j64 -r256 -m bench.dat bench.rz 2>&1)
bsize=$(stat -c "%s" bench.dat)
perf=$(echo "$bsize/1048576/$utime" | bc)
echo "[0;32m*** Encoding with $perf MiB/s user time ***[0m"
