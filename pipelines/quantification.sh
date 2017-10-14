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
            SEQDIR=$BASEDIR/simulation/${SPECIES}_${DEPTH}/$SEQ
            MRNADIR=$BASEDIR/simulation/${SPECIES}_${DEPTH}/mRNA
            READDIR=$BASEDIR/simulation/${SPECIES}_${DEPTH}/reads
            
            FLMEAN="$(grep 'insert size average:' $MRNADIR/bwa/stats.txt | grep -oP '(\d+\.*\d+)')"
            FLSD="$(grep 'insert size standard deviation:' $MRNADIR/bwa/stats.txt | grep -oP '(\d+\.*\d+)')"
            
            mkdir -p $SEQDIR/kallisto
            kallisto index -i $SEQDIR/kallisto/kallisto.index -k 31 $SEQDIR/${SEQ}.fasta
            kallisto quant -i $SEQDIR/kallisto/kallisto.index -o $SEQDIR/kallisto -t $THREADS $READDIR/read_1.fastq $READDIR/read_2.fastq > $SEQDIR/kallisto/kallisto.out 2> $SEQDIR/kallisto/kallisto.err

            mkdir -p $SEQDIR/rsem
            rsem-prepare-reference -p $THREADS --bowtie2 $SEQDIR/${SEQ}.fasta $SEQDIR/rsem/rsem.index
            rsem-calculate-expression --paired-end --strandedness none -p $THREADS --bowtie2 --time $READDIR/read_1.fastq $READDIR/read_2.fastq $SEQDIR/rsem/rsem.index $SEQDIR/rsem/rsem > $SEQDIR/rsem/rsem.out 2> $SEQDIR/rsem/rsem.err

            mkdir -p $SEQDIR/salmon
            salmon index -i $SEQDIR/salmon/salmon.index -t $SEQDIR/${SEQ}.fasta --type quasi -k 31
            salmon quant -i $SEQDIR/salmon/salmon.index -p $THREADS -l A -1 $READDIR/read_1.fastq -2 $READDIR/read_2.fastq -o $SEQDIR/salmon > $SEQDIR/salmon/salmon.out 2> $SEQDIR/salmon/salmon.err
        done
    done
done
