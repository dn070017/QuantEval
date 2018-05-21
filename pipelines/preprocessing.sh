#!/usr/bin/env bash

set -e

THREADS=32
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"
TRIMMOMATIC="$HOME/bin/Trimmomatic-0.36/trimmomatic-0.36.jar"

for TYPE in 'real' 'simulation'
do
    if [ $TYPE == 'simulation' ]
    then
        READ="flux_simulator"
    else
        READ="raw"
    fi
    
    for SPECIES in 'yeast' 'dog' 'mouse'
    do
        READDIR=$BASEDIR/$TYPE/$SPECIES/reads
        SIMDIR=$BASEDIR/$TYPE/$SPECIES/simulation
        MRNADIR=$BASEDIR/$TYPE/$SPECIES/mRNA
       
        mkdir -p $MRNADIR
        mkdir -p $READDIR/fastqc
        fastqc -t $THREADS -f fastq -o $READDIR/fastqc --nogroup $READDIR/${READ}_r*.fastq

        java -jar $TRIMMOMATIC PE -threads $THREADS -phred33 $READDIR/${READ}_r1.fastq $READDIR/${READ}_r2.fastq $READDIR/read_1.fastq $READDIR/unpaired_read_1.fastq $READDIR/read_2.fastq $READDIR/unpaired_read_2.fastq SLIDINGWINDOW:4:20 MINLEN:30
        fastqc -t $THREADS -f fastq -o $READDIR/fastqc --nogroup $READDIR/read_*.fastq
            
        if [ $TYPE == 'simulation' ]
        then
            mkdir -p $MRNADIR/answer/
            $BASEDIR/scripts/count_unpaired_read.py $READDIR/unpaired_read_1.fastq $READDIR/unpaired_read_2.fastq $READDIR/unpaired_read_count.picklee
            $BASEDIR/scripts/generate_answer_tpm.py $TYPE $SIMDIR/flux_simulator.pro $SIMDIR/flux_simulator.lib $READDIR/unpaired_read_count.pickle $MRNADIR/answer/answer_tpm.tsv
        fi

        REFDIR=$BASEDIR/reference/$SPECIES
        mkdir -p $MRNADIR/bwa
        ln $REFDIR/mRNA.fasta $MRNADIR
        ln $REFDIR/mRNA.fasta $MRNADIR/bwa
        bwa index -a is -p $MRNADIR/bwa/mRNA $MRNADIR/bwa/mRNA.fasta
        bwa mem -t $THREADS $MRNADIR/bwa/mRNA $READDIR/read_1.fastq $READDIR/read_2.fastq | samtools sort -@ $THREADS -O BAM -o $MRNADIR/bwa/bwa.bam - 
        samtools stats $MRNADIR/bwa/bwa.bam > $MRNADIR/bwa/stats.txt 
    done
done
