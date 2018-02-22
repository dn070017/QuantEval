#!/usr/bin/env bash

set -e

THREADS=32
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"

for SPECIES in 'yeast' 'dog' 'mouse'
do
    if [ $SPECIES == 'mouse' ]
    then
        #DEPTHS=( "25x" "50x" "100x" )
        DEPTHS=( "50x" )
    else
        DEPTHS=( "50x" )
    fi

    for DEPTH in "${DEPTHS[@]}"
    do
        for SEQ in 'mRNA' 'rnaspades' 'trinity' 'transabyss'
        do
            $BASEDIR/scripts/extract_features.py $BASEDIR/simulation/${SPECIES}_$DEPTHS $BASEDIR/reference/$SPECIES $SEQ
        done
    done
done
