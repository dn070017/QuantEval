#!/usr/bin/env bash

set -e

THREADS=32
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"
TRIMMOMATIC="$HOME/bin/Trimmomatic-0.36/trimmomatic-0.36.jar"

mkdir -p "$HOME/temp"
export TMP_DIR="$HOME/temp"

for SPECIES in 'mouse' 'dog' 'yeast'
do
    REFDIR="$BASEDIR/reference/$SPECIES"
    if [ $SPECIES == 'mouse' ]
    then
        #DEPTHS=( "25x" "50x" "100x" )
        DEPTHS=( "50x" )
    else
        DEPTHS=( "50x" )
    fi

    for DEPTH in "${DEPTHS[@]}"
    do
        SIMDIR=$BASEDIR/simulation/${SPECIES}_${DEPTH}/simulation/
        READDIR=$BASEDIR/simulation/${SPECIES}_${DEPTH}/reads
        mkdir -p $SIMDIR
        cp -r "$REFDIR/chromosome" $SIMDIR
        cp -r "$REFDIR/flux_simulator_clean.gtf" $SIMDIR/
        
        flux-simulator -p $SIMDIR/flux_simulator.par -l -s -x > $SIMDIR/flux_simulator.out 2> $SIMDIR/flux_simulator.err 
        
        $BASEDIR/scripts/split_interleaved_reads.py $SIMDIR/flux_simulator.fastq
        
        mkdir -p $READDIR
        ln $SIMDIR/flux_simulator_r*.fastq $READDIR
    done
done
