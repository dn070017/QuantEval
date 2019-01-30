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
    SIM_50X_DIR=$BASEDIR/simulation_50x/$SPECIES/simulation/
    SIM_50M_DIR=$BASEDIR/simulation_50M/$SPECIES/simulation/
    READDIR=$BASEDIR/simulation_50M/$SPECIES/reads
    
    mkdir -p $SIM_50M_DIR
    mkdir -p $READDIR

    cp -r "$SIM_50X_DIR/chromosome" $SIM_50M_DIR
    cp -r "$SIM_50X_DIR/flux_simulator_clean.gtf" $SIM_50M_DIR
    cp -r "$SIM_50X_DIR/flux_simulator.par" $SIM_50M_DIR
    cp -r "$SIM_50X_DIR/flux_simulator.pro" $SIM_50M_DIR

    flux-simulator -p $SIM_50M_DIR/flux_simulator.par -l -s > $SIM_50M_DIR/flux_simulator.out 2> $SIM_50M_DIR/flux_simulator.err 
        
    $BASEDIR/scripts/split_interleaved_reads.py $SIM_50M_DIR/flux_simulator.fastq
        
    ln $SIM_50M_DIR/flux_simulator_r*.fastq $READDIR
done
