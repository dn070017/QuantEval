#!/usr/bin/env bash

set -e

THREADS=32
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"
TRIMMOMATIC="$HOME/bin/Trimmomatic-0.36/trimmomatic-0.36.jar"

#for TYPE in 'real' 'profile' 'simulation'
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
            if [ $TYPE == 'profile' ]
            then
                mkdir -p $READDIR
                #ln $BASEDIR/real/$SPECIES$SURFFIX/mRNA/rsem_sim/rsem_sim_1.fq $READDIR/${READ}_r1.fastq
                #ln $BASEDIR/real/$SPECIES$SURFFIX/mRNA/rsem_sim/rsem_sim_2.fq $READDIR/${READ}_r2.fastq
            fi
            mkdir -p $READDIR/fastqc
            #fastqc -t $THREADS -f fastq -o $READDIR/fastqc --nogroup $READDIR/${READ}_r*.fastq

            #java -jar $TRIMMOMATIC PE -threads $THREADS -phred33 $READDIR/${READ}_r1.fastq $READDIR/${READ}_r2.fastq $READDIR/read_1.fastq $READDIR/unpaired_read_1.fastq $READDIR/read_2.fastq $READDIR/unpaired_read_2.fastq SLIDINGWINDOW:4:20 MINLEN:30
            #fastqc -t $THREADS -f fastq -o $READDIR/fastqc --nogroup $READDIR/read_*.fastq

            $BASEDIR/scripts/count_unpaired_read.py $READDIR/unpaired_read_1.fastq $READDIR/unpaired_read_2.fastq $READDIR/unpaired_read_count.pickle

            REFDIR=$BASEDIR/reference/$SPECIES
            MRNADIR=$BASEDIR/$TYPE/$SPECIES$SURFFIX/mRNA
            mkdir -p $MRNADIR/bwa
            #ln $REFDIR/mRNA.fasta $MRNADIR
            #ln $REFDIR/mRNA.fasta $MRNADIR/bwa
            #bwa index -a is -p $MRNADIR/bwa/mRNA $MRNADIR/bwa/mRNA.fasta
            #bwa mem -t $THREADS $MRNADIR/bwa/mRNA $READDIR/read_1.fastq $READDIR/read_2.fastq | samtools sort -@ $THREADS -O BAM -o $MRNADIR/bwa/bwa.bam - 
            #samtools stats $MRNADIR/bwa/bwa.bam > $MRNADIR/bwa/stats.txt 
            
            if [ $TYPE == 'real' ]
            then
                mkdir -p $MRNADIR/rsem_sim
                #rsem-prepare-reference -p $THREADS --bowtie2 $MRNADIR/mRNA.fasta $MRNADIR/rsem_sim/rsem.index
                if [ $SPECIES == 'mouse' ]
                then
                    STRAND='reverse'
                else
                    STRAND='none'
                fi
                #rsem-calculate-expression --paired-end --strandedness $STRAND -p $THREADS --bowtie2 --time $READDIR/read_1.fastq $READDIR/read_2.fastq $MRNADIR/rsem_sim/rsem.index $MRNADIR/rsem_sim/rsem > $MRNADIR/rsem_sim/rsem.out 2> $MRNADIR/rsem_sim/rsem.err
                
                THETA=`awk 'NR == 3 {print $1}' $MRNADIR/rsem_sim/rsem.stat/rsem.theta`
                N=`awk 'NR == 1 {print $4}' $MRNADIR/rsem_sim/rsem.stat/rsem.cnt`

                #rsem-simulate-reads $MRNADIR/rsem_sim/rsem.index $MRNADIR/rsem_sim/rsem.stat/rsem.model $MRNADIR/rsem_sim/rsem.isoforms.results $THETA $N $MRNADIR/rsem_sim/rsem_sim
            fi
        done
    done
done
