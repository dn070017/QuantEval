#!/usr/bin/env bash

set -e

THREADS=32
IDENTITY=95
EVALUE="1e-5"

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"

for TYPE in 'real'
do
    for SPECIES in 'yeast' 'dog' 'mouse'
    do
        if [ $TYPE == 'simulation' ] && [ $SPECIS == 'mouse' ]
        then
            #DEPTHS=( "25x" "50x" "100x" )
            DEPTHS=( "_50x" )
        elif [ $TYPE == 'simulation' ]
        then
            DEPTHS=( "_50x" )
        else
            DEPTHS=( "" )
        fi

        for DEPTH in "${DEPTHS[@]}"
        do
            for SEQ in 'mRNA' 'rnaspades' 'trinity' 'transabyss'
            do
                SEQDIR=$BASEDIR/$TYPE/${SPECIES}${DEPTH}/$SEQ
                MRNADIR=$BASEDIR/$TYPE/${SPECIES}${DEPTH}/mRNA
                READDIR=$BASEDIR/$TYPE/${SPECIES}${DEPTH}/reads
            
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
done
