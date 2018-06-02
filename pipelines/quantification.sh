#!/usr/bin/env bash

set -e

THREADS=32
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"

for TYPE in 'real' 'simulation'
do
    for SPECIES in 'yeast' 'dog' 'mouse'
    do
        if [ $TYPE != 'simulation' ] && [ $SPECIES == 'mouse' ]
        then
            SS_KALLISTO="--rf-stranded"
            SS_RSEM="reverse"
            #SS_SALMON="-l ISR"
        else
            SS_KALLISTO=""
            SS_RSEM="none"
            #SS_SALMON=""
        fi

        for SEQ in 'mRNA' 'rnaspades' 'trinity' 'transabyss'
        do
            SEQDIR=$BASEDIR/$TYPE/$SPECIES/$SEQ
            MRNADIR=$BASEDIR/$TYPE/$SPECIES/mRNA
            READDIR=$BASEDIR/$TYPE/$SPECIES/reads
           
            FLMEAN="$(grep 'insert size average:' $MRNADIR/bwa/stats.txt | grep -oP '(\d+\.*\d+)')"
            FLSD="$(grep 'insert size standard deviation:' $MRNADIR/bwa/stats.txt | grep -oP '(\d+\.*\d+)')"
            
            mkdir -p $SEQDIR/kallisto
            kallisto index -i $SEQDIR/kallisto/kallisto.index -k 31 $SEQDIR/${SEQ}.fasta
            kallisto quant -i $SEQDIR/kallisto/kallisto.index -o $SEQDIR/kallisto -t $THREADS $SS_KALLISTO $READDIR/read_1.fastq $READDIR/read_2.fastq > $SEQDIR/kallisto/kallisto.out 2> $SEQDIR/kallisto/kallisto.err

            mkdir -p $SEQDIR/rsem
            rsem-prepare-reference -p $THREADS --bowtie2 $SEQDIR/${SEQ}.fasta $SEQDIR/rsem/rsem.index
            rsem-calculate-expression --paired-end --strandedness $SS_RSEM -p $THREADS --bowtie2 --time $READDIR/read_1.fastq $READDIR/read_2.fastq $SEQDIR/rsem/rsem.index $SEQDIR/rsem/rsem > $SEQDIR/rsem/rsem.out 2> $SEQDIR/rsem/rsem.err

            mkdir -p $SEQDIR/salmon
            salmon index -i $SEQDIR/salmon/salmon.index -t $SEQDIR/${SEQ}.fasta --type quasi -k 31
            salmon quant -i $SEQDIR/salmon/salmon.index -p $THREADS -l A -1 $READDIR/read_1.fastq -2 $READDIR/read_2.fastq -o $SEQDIR/salmon > $SEQDIR/salmon/salmon.out 2> $SEQDIR/salmon/salmon.err

            mkdir -p $SEQDIR/bowtie2
            ln $SEQDIR/${SEQ}.fasta $SEQDIR/bowtie2/${SEQ}.fasta
            bowtie2-build $SEQDIR/bowtie2/${SEQ}.fasta $SEQDIR/bowtie2/bowtie2.index
            bowtie2 -q --phred33 --no-mixed --no-discordant -p $THREADS -k 200 -x $SEQDIR/bowtie2/bowtie2.index -1 $READDIR/read_1.fastq -2 $READDIR/read_2.fastq 2> $SEQDIR/bowtie2/summary.txt | samtools sort -n -o $SEQDIR/bowtie2/bowtie2.bam -
            samtools sort -@ $THREADS $SEQDIR/bowtie2/bowtie2.bam -o $SEQDIR/bowtie2/bowtie2.sorted.bam
            samtools index $SEQDIR/bowtie2/bowtie2.sorted.bam $SEQDIR/bowtie2/bowtie2.sorted.bai
            
            mkdir -p $SEQDIR/corset
            corset -p $SEQDIR/corset/corset $SEQDIR/bowtie2/bowtie2.bam > $SEQDIR/corset/corset.out 2> $SEQDIR/corset/corset.err

            if [ $SEQ == 'mRNA' ] && [ $TYPE == 'real' ]
            then
                mkdir -p $SEQDIR/answer
                $BASEDIR/scripts/generate_answer_xprs.py $TYPE $SEQDIR/kallisto/abundance.tsv $SEQDIR/rsem/rsem.isoforms.results $SEQDIR/salmon/quant.sf $SEQDIR/answer/answer_tpm.tsv
            fi
        done
    done
done
