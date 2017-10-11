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
        READDIR=$BASEDIR/simulation/${SPECIES}_${DEPTH}/reads 
        mkdir -p $READDIR/fastqc
        ln $SIMDIR/flux_simulator_r*.fastq $READDIR
        fastqc -t $THREADS -f fastq -o $READDIR/fastqc --nogroup $READDIR/flux_simulator_r*.fastq

        java -jar $TRIMMOMATIC PE -threads $THREADS -phred33 $READDIR/flux_simulator_r1.fastq $READDIR/flux_simulator_r2.fastq $READDIR/read_1.fastq $READDIR/unpaired_read_1.fastq $READDIR/read_2.fastq $READDIR/unpaired_read_2.fastq SLIDINGWINDOW:4:20 MINLEN:30
        fastqc -t $THREADS -f fastq -o $READDIR/fastqc --nogroup $READDIR/read_*.fastq

        REFDIR=$BASEDIR/reference/$SPECIES
        MRNADIR=$BASEDIR/simulation/${SPECIES}_${DEPTH}/mRNA
        mkdir -p $MRNADIR/bwa
        ln $REFDIR/mRNA.fasta $MRNADIR
        ln $REFDIR/mRNA.fasta $MRNADIR/bwa
        bwa index -a is -p $MRNADIR/bwa/mRNA $MRNADIR/bwa/mRNA.fasta
        bwa mem -t $THREADS $MRNADIR/bwa/mRNA $READDIR/read_1.fastq $READDIR/read_2.fastq | samtools sort -@ $THREADS -O BAM -o $MRNADIR/bwa/bwa.bam - 
        samtools stats $MRNADIR/bwa/bwa.bam > $MRNADIR/bwa/stats.txt 
    done
done
