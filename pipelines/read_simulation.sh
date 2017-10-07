#!/usr/bin/env bash

set -e

THREADS=32
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"
TRIMMOMATIC="$HOME/bin/Trimmomatic-0.36/trimmomatic-0.36.jar"

for SPECIES in 'mouse' 'dog' 'yeast'
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
        SIMDIR=$BASEDIR/simulation/${SPECIES}_${DEPTH}/simulation/
        #$mkdir -p $SIMDIR
        $BASEDIR/scripts/prepare_reference.py $BASEDIR/reference/$SPECIES/ $SIMDIR 500 > $SIMDIR/log 
        awk 'BEGIN{FS="\t";OFS="\t"}{split($NF,a," ");pfx="";s="";for(i=1;i<=length(a);i+=2){if(a[i]=="transcript_id"){pfx=a[i]" "a[i+1]}else{s=s" "a[i]" "a[i+1]}}if(pfx==""){print "[WARN] line "NR" without transcript_id!" > "/dev/stderr"}else{$NF=pfx""s;print$0} }' $SIMDIR/flux_simulator.gtf > $SIMDIR/flux_simulator_clean.gtf
        flux-simulator -p $SIMDIR/flux_simulator.par -l -s -x > $SIMDIR/flux_simulator.out 2> $SIMDIR/flux_simulator.err 
        $BASEDIR/scripts/split_interleaved_reads.py $SIMDIR/flux_simulator.fastq
        
        READDIR=$BASEDIR/simulation/${SPECIES}_${DEPTH}/reads 
        mkdir -p $READDIR/fastqc
        ln $SIMDIR/flux_simulator_r*.fastq $READDIR
        fastqc -t $THREADS -f fastq -o $READDIR/fastqc --nogroup $READDIR/flux_simulator_r*.fastq

        java -jar $TRIMMOMATIC PE -threads $THREADS -phred33 $READDIR/flux_simulator_r1.fastq $READDIR/flux_simulator_r2.fastq $READDIR/read_1.fastq $READDIR/unpaired_read_1.fastq $READDIR/read_2.fastq $READDIR/unpaired_read_2.fastq SLIDINGWINDOW:4:20 MINLEN:30
        fastqc -t $THREADS -f fastq -o $READDIR/fastqc --nogroup $READDIR/read_*.fastq
    done
done
