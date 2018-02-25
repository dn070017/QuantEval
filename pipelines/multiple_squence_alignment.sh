#!/usr/bin/env bash

set -e

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"

for SPECIES in 'yeast' 'dog' 'mouse'
do 
    if [ $SPECIES == 'mouse' ]
    then
        DEPTHS=( "50x" )
    else
        DEPTHS=( "50x" )
    fi
    
    for DEPTH in "${DEPTHS[@]}"
    do
        for SEQ in 'rnaspades' 'trinity' 'transabyss'
        do
            SEQDIR=$BASEDIR/simulation/${SPECIES}_${DEPTH}/$SEQ
            MRNADIR=$BASEDIR/simulation/${SPECIES}_${DEPTH}/mRNA
            REFDIR=$BASEDIR/reference/${SPECIES}
            
            mkdir -p $SEQDIR/mafft/duplication
            mkdir -p $SEQDIR/mafft/family_collapse
            
            MSADIR=$SEQDIR/mafft/duplication
            awk -F$'\t' "\$1 != \"t_name\" { print \$1}" $MSADIR/duplication.tsv | sort | uniq > $MSADIR/t_name.txt
            while read t_name
            do
                #echo $t_name > $MSADIR/_t_name.txt
                #awk -F$'\t' "\$3 == \"$t_name\" {print \$2}" $SEQDIR/features/features.tsv > $MSADIR/_c_name.txt
                #get_seqs_from_id.pl -k $MSADIR/_t_name.txt -t $REFDIR/mRNA.fasta > $MSADIR/t.fasta
                #get_seqs_from_id.pl -k $MSADIR/_c_name.txt -t $SEQDIR/${SEQ}.fasta > $MSADIR/c.fasta
                #cat $MSADIR/t.fasta $MSADIR/c.fasta > $MSADIR/${t_name}.fasta
                #rm $MSADIR/t.fasta $MSADIR/c.fasta $MSADIR/_t_name.txt $MSADIR/_c_name.txt
                
                mafft --thread 24 --adjustdirection --globalpair --clustalout --quiet --maxiterate 1000 $MSADIR/${t_name}.fasta > $MSADIR/${t_name}.clustal.out
            done <$MSADIR/t_name.txt
            
            if [ $SPECIES == 'yeast' ]
            then
                MSADIR=$SEQDIR/mafft/family_collapse
                awk -F$'\t' "\$2 != \"c_name\" { print \$2}" $MSADIR/family_collapse.tsv | sort | uniq > $MSADIR/c_name.txt
                while read c_name
                do
                    #echo $c_name > $MSADIR/_c_name.txt
                    #awk -F$'\t' "\$2 == \"$c_name\" {print \$3}" $SEQDIR/features/features.tsv > $MSADIR/_t_name.txt
                    #get_seqs_from_id.pl -k $MSADIR/_t_name.txt -t $REFDIR/mRNA.fasta > $MSADIR/t.fasta
                    #get_seqs_from_id.pl -k $MSADIR/_c_name.txt -t $SEQDIR/${SEQ}.fasta > $MSADIR/c.fasta
                    #cat $MSADIR/t.fasta $MSADIR/c.fasta > $MSADIR/${c_name}.fasta
                    #rm $MSADIR/t.fasta $MSADIR/c.fasta $MSADIR/_t_name.txt $MSADIR/_c_name.txt
                
                    mafft --thread 24 --adjustdirection --globalpair --clustalout --quiet --maxiterate 1000 $MSADIR/${c_name}.fasta > $MSADIR/${c_name}.clustal.out
                done <$MSADIR/c_name.txt
            fi

        done
    done
done
