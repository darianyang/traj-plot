#!/bin/bash
# copy.sh

SYSTEMS=(w4f w5f w6f w7f wt)
#SYSTEMS=(w4f)
FF="ipq"
OUT_ROOT=1us_noion
SOURCE=/bgfs/lchong/dty7/19F_ff15ipq/CypA/$FF

for SYS in ${SYSTEMS[@]} ; do
    mkdir -p $FF/$SYS
    cd $FF/$SYS
    for VER in $(seq -f "%02g" 0 1 4) ; do
        mkdir v$VER
        rsync -axhvP dty7@h2p.crc.pitt.edu:$SOURCE/$SYS/v$VER/$OUT_ROOT/*.dat v$VER/$OUT_ROOT/
        #rsync -axhvP dty7@h2p.crc.pitt.edu:$SOURCE/$SYS/v$VER/$OUT_ROOT/qmgbsa/*.dat v$VER/$OUT_ROOT/
        #rsync -axhvP dty7@h2p.crc.pitt.edu:$SOURCE/$SYS/v$VER/$OUT_ROOT/mmpbsa/*INP1.dat v$VER/$OUT_ROOT/
        #scp dty7@h2p.crc.pitt.edu:$SOURCE/$SYS/v$VER/$OUT_ROOT/* v$VER/$OUT_ROOT/
    done
    cd ../..
done
