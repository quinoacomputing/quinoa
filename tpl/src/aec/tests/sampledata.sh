#!/bin/sh
set -e
AEC=../src/aec
CCSDS_DATA=121B2TestData
ALLO=${CCSDS_DATA}/AllOptions
EXTP=${CCSDS_DATA}/ExtendedParameters
LOWE=${CCSDS_DATA}/LowEntropyOptions
archive=121B2TestData.zip

if [ ! -f $archive ]; then
    wget http://cwe.ccsds.org/sls/docs/SLS-DC/BB121B2TestData/$archive
fi
unzip -oq $archive

decode () {
    $AEC -d $3 $1 test.dat
    cmp -n $(stat -c "%s" $2) $2 test.dat
}

code () {
    $AEC $3 $2 test.rz
    cmp $1 test.rz
}

code_size () {
    $AEC $3 $2 test.rz
    if [ ! $(stat -c "%s" test.rz) -eq $(stat -c "%s" $1) ]; then
        echo "$1 size mismatch"
        exit 1
    fi
}

codec () {
    code "$@"
    decode "$@"
}

cosdec () {
    code_size "$@"
    decode "$@"
}

echo All Options
ln -f ${ALLO}/test_P512n22.dat ${ALLO}/test_p512n22.dat
for i in 01 02 03 04
do
    uf=$ALLO/test_p256n${i}.dat
    codec $ALLO/test_p256n${i}-basic.rz $uf "-n$i -j16 -r16"
    codec $ALLO/test_p256n${i}-restricted.rz $uf "-n$i -j16 -r16 -t"
done
for i in 05 06 07 08 09 10 11 12 13 14 15 16
do
    cosdec $ALLO/test_p256n${i}.rz $ALLO/test_p256n${i}.dat \
        "-n$i -j16 -r16"
done
for i in 17 18 19 20 21 22 23 24
do
    cosdec $ALLO/test_p512n${i}.rz $ALLO/test_p512n${i}.dat \
        "-n$i -j16 -r32"
done

echo Low Entropy Options
for i in 1 2 3
do
    for j in 01 02 03 04
    do
        uf=$LOWE/Lowset${i}_8bit.dat
        codec $LOWE/Lowset${i}_8bit.n${j}-basic.rz $uf "-n$j -j16 -r64"
        codec $LOWE/Lowset${i}_8bit.n${j}-restricted.rz $uf "-n$j -j16 -r64 -t"
    done
    for j in 05 06 07 08
    do
        codec $LOWE/Lowset${i}_8bit.n${j}.rz $LOWE/Lowset${i}_8bit.dat \
            "-n$j -j16 -r64"
    done
done

echo Extended Parameters
decode $EXTP/sar32bit.j16.r256.rz $EXTP/sar32bit.dat "-n32 -j16 -r256 -p"
decode $EXTP/sar32bit.j64.r4096.rz $EXTP/sar32bit.dat "-n32 -j64 -r4096 -p"
