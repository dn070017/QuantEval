#!/usr/bin/env bash

set -e

THREADS=32
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"

mkdir -p "$HOME/temp"
export TMP_DIR="$HOME/temp"

for SPECIES in 'mouse' 'dog' 'yeast'
do
    REFDIR="$BASEDIR/reference/$SPECIES" 
    $BASEDIR/scripts/prepare_reference.py $REFDIR 500 > $REFDIR/prepare_reference.log 
    awk 'BEGIN{FS="\t";OFS="\t"}{split($NF,a," ");pfx="";s="";for(i=1;i<=length(a);i+=2){if(a[i]=="transcript_id"){pfx=a[i]" "a[i+1]}else{s=s" "a[i]" "a[i+1]}}if(pfx==""){print "[WARN] line "NR" without transcript_id!" > "/dev/stderr"}else{$NF=pfx""s;print$0} }' $REFDIR/flux_simulator.gtf > $REFDIR/flux_simulator_clean.gtf
done
