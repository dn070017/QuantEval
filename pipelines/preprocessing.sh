#!/usr/bin/env bash

set -e

THREADS=32
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"
TRIMMOMATIC="$HOME/bin/Trimmomatic-0.36/trimmomatic-0.36.jar"

#for TYPE in 'real' 'simulation'
for TYPE in 'real'
do
    if [ $TYPE == 'real' ]
    then
        READ="raw"
    else
        READ="flux_simulator"
    fi
    
    #for SPECIES in 'mouse' 'dog' 'yeast'
    for SPECIES in 'yeast'
    do
        if [ $TYPE == 'simulation' ] && [ $SPECIES == 'mouse' ]
        then
            #SURFFIXES=( "_25x" "_50x" "_100x" )
            SURFFIXED=( "_50x" )
        elif [ $TYPE == 'simulation' ]
        then
            SURFFIXES=( "_50x" )
        else
            SURFFIXES=( "" )
        fi
        
        for SURFFIX in "${SURFFIXES[@]}" 
        do
            READDIR=$BASEDIR/$TYPE/$SPECIES$SURFFIX/reads
            mkdir -p $READDIR/fastqc
            fastqc -t $THREADS -f fastq -o $READDIR/fastqc --nogroup $READDIR/${READ}_r*.fastq

            java -jar $TRIMMOMATIC PE -threads $THREADS -phred33 $READDIR/${READ}_r1.fastq $READDIR/${READ}_r2.fastq $READDIR/read_1.fastq $READDIR/unpaired_read_1.fastq $READDIR/read_2.fastq $READDIR/unpaired_read_2.fastq SLIDINGWINDOW:4:20 MINLEN:30
            fastqc -t $THREADS -f fastq -o $READDIR/fastqc --nogroup $READDIR/read_*.fastq

            REFDIR=$BASEDIR/reference/$SPECIES
            MRNADIR=$BASEDIR/simulation/$SPECIES/mRNA
            mkdir -p $MRNADIR/bwa
            ln $REFDIR/mRNA.fasta $MRNADIR
            ln $REFDIR/mRNA.fasta $MRNADIR/bwa
            bwa index -a is -p $MRNADIR/bwa/mRNA $MRNADIR/bwa/mRNA.fasta
            bwa mem -t $THREADS $MRNADIR/bwa/mRNA $READDIR/read_1.fastq $READDIR/read_2.fastq | samtools sort -@ $THREADS -O BAM -o $MRNADIR/bwa/bwa.bam - 
            samtools stats $MRNADIR/bwa/bwa.bam > $MRNADIR/bwa/stats.txt 
        done
    done
done
