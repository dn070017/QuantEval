#!/usr/bin/env bash

set -e

THREADS=32
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"
TRIMMOMATIC="$HOME/bin/Trimmomatic-0.36/trimmomatic-0.36.jar"

mkdir -p "$HOME/temp"
export TMP_DIR="$HOME/temp"

for FOLDER in 'simulation_50x' 'simulation_50M'
do
    for SPECIES in 'mouse' 'dog' 'yeast'
    do
        REFDIR="$BASEDIR/reference/$SPECIES"
        SIMDIR=$BASEDIR/$FOLDER/$SPECIES/simulation/
        READDIR=$BASEDIR/$FOLDER/$SPECIES/reads
        
        mkdir -p $SIMDIR
        mkdir -p $READDIR
    
        flux-simulator -p $SIMDIR/flux_simulator.par -l -s > $SIMDIR/flux_simulator.out 2> $SIMDIR/flux_simulator.err 
    done
done
