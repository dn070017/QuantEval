#!/usr/bin/env bash

set -e

MEMORY=350
THREADS=32
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"

for TYPE in 'real_low' 'real_high' 'simulation_50x' 'simulation_50M'
do
    for SPECIES in 'yeast' 'dog' 'mouse'
    do
        if [ $TYPE != 'simulation_50x' ] && [ $TYPE != 'simulation_50M' ] && [ $SPECIES == 'mouse' ]
        then
            SS_RNASPADES="--ss-rf"
            SS_TRANSABYSS="--SS"
            SS_TRINITY="--SS_lib_type RF"
        else
            SS_RNASPADES=""
            SS_TRANSABYSS=""
            SS_TRINITY=""
        fi

        DATADIR=$BASEDIR/$TYPE/$SPECIES/
        MRNADIR=$BASEDIR/$TYPE/$SPECIES/mRNA
        READDIR=$BASEDIR/$TYPE/$SPECIES/reads
            
        #FLMEAN="$(grep 'insert size average:' $MRNADIR/bwa/stats.txt | grep -oP '(\d+\.*\d+)')"
        #FLSD="$(grep 'insert size standard deviation:' $MRNADIR/bwa/stats.txt | grep -oP '(\d+\.*\d+)')"
            
        mkdir -p $DATADIR/rnaspades/rnaspades
        rnaspades.py -t $THREADS -m $MEMORY $SS_RNASPADES -1 $READDIR/read_1.fastq -2 $READDIR/read_2.fastq -o $DATADIR/rnaspades/rnaspades > $DATADIR/rnaspades/rnaspades/rnaspades.out 2> $DATADIR/rnaspades/rnaspades/rnaspades.err 
        $BASEDIR/scripts/fasta_length_filter.py $DATADIR/rnaspades/rnaspades/transcripts.fasta $DATADIR/rnaspades/rnaspades.fasta 500
        rm -rf $DATADIR/rnaspades/rnaspades

        mkdir -p $DATADIR/transabyss/transabyss
        transabyss --pe $READDIR/read_1.fastq $READDIR/read_2.fastq --outdir $DATADIR/transabyss/transabyss $SS_TRANSABYSS --threads $THREADS --length 500 > $DATADIR/transabyss/transabyss/transabyss.out 2> $DATADIR/transabyss/transabyss/transabyss.err
        ln $DATADIR/transabyss/transabyss/transabyss-final.fa $DATADIR/transabyss/transabyss.fasta
        rm -rf $DATADIR/transabyss/transabyss

        mkdir -p $DATADIR/trinity/trinity
        Trinity --seqType fq --max_memory ${MEMORY}G --left $READDIR/read_1.fastq --right $READDIR/read_2.fastq $SS_TRINITY --CPU $THREADS --min_contig_length 500 --output $DATADIR/trinity/trinity > $DATADIR/trinity/trinity/trinity.out 2> $DATADIR/trinity/trinity/trinity.err
        ln $DATADIR/trinity/trinity/Trinity.fasta $DATADIR/trinity/trinity.fasta
        rm -rf $DATADIR/trinity/trinity
    done
done
