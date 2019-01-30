#!/usr/bin/env bash

set -e

THREADS=32
IDENTITY=95
EVALUE="1e-5"

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"

for TYPE in 'real_low' 'real_high' 'simulation_low' 'simulation_high'
do
    for SPECIES in 'yeast' 'dog' 'mouse'
    do
        #for SEQ in 'mRNA' 'rnaspades' 'trinity' 'transabyss'
        for SEQ in 'rnaspades' 'trinity' 'transabyss'
        do
            SEQDIR=$BASEDIR/$TYPE/$SPECIES/$SEQ
            MRNADIR=$BASEDIR/$TYPE/$SPECIES/mRNA
            READDIR=$BASEDIR/$TYPE/$SPECIES/reads
            
            #FLMEAN="$(grep 'insert size average:' $MRNADIR/bwa/stats.txt | grep -oP '(\d+\.*\d+)')"
            #FLSD="$(grep 'insert size standard deviation:' $MRNADIR/bwa/stats.txt | grep -oP '(\d+\.*\d+)')"
            
            mkdir -p $SEQDIR/blastn
            ln $SEQDIR/${SEQ}.fasta $SEQDIR/blastn
            makeblastdb -in $SEQDIR/blastn/${SEQ}.fasta -dbtype nucl
            blastn -db $SEQDIR/blastn/${SEQ}.fasta -query $SEQDIR/${SEQ}.fasta -outfmt 6 -evalue $EVALUE -perc_identity $IDENTITY -out $SEQDIR/blastn/self.tsv 
            if [ SEQ != 'mRNA' ]
            then
                blastn -db $SEQDIR/blastn/${SEQ}.fasta -query $MRNADIR/mRNA.fasta -outfmt 6 -out $SEQDIR/blastn/mRNA_to_contig.tsv 
                blastn -db $MRNADIR/blastn/mRNA.fasta -query $SEQDIR/${SEQ}.fasta -outfmt 6 -out $SEQDIR/blastn/contig_to_mRNA.tsv 
                mkdir -p $SEQDIR/transrate
                transrate --assembly $SEQDIR/${SEQ}.fasta --threads $THREADS --output $SEQDIR/transrate --left $READDIR/read_1.fastq --right $READDIR/read_2.fastq
            fi
        done
    done
done
