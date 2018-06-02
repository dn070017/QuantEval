#!/usr/bin/env bash

set -e

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEDIR="$SCRIPTDIR/../"

for TYPE in 'real' 'simulation'
do 
    for SPECIES in 'yeast' 'dog' 'mouse'
    do
        for SEQ in 'mRNA' 'rnaspades' 'trinity' 'transabyss'
        do
            
            echo "processing $TYPE $SPECIES $SEQ"

            SEQDIR=$BASEDIR/$TYPE/$SPECIES/$SEQ
            REFDIR=$BASEDIR/$TYPE/$SPECIES/mRNA
            
            if [ $SEQ == 'mRNA' ]
            then
                cat > $REFDIR/QuantEval.json << EOF
{
    "ref_fasta": "$REFDIR/mRNA.fasta",
    "ref_blastn": "$REFDIR/blastn/self.tsv",
    "ref_gtf": "$BASEDIR/simulation/$SPECIES/simulation/flux_simulator_clean.gtf",
    "ref_write_pickle": "$REFDIR/QuantEval/ref.pickle",
    "ref_xprs_file": ["$REFDIR/answer/answer_xprs.tsv",
                      "$REFDIR/kallisto/abundance.tsv",
                      "$REFDIR/rsem/rsem.isoforms.results",
                      "$REFDIR/salmon/quant.sf"],
    "ref_xprs_label": ["answer", "kallisto", "rsem", "salmon"],
    "ref_xprs_header": [true, true, true, true],
    "ref_xprs_name_col": [1, 1, 1, 1], 
    "ref_xprs_tpm_col": [2, 5, 6, 4], 
    "ref_xprs_count_col": [3, 4, 5, 5],
    "ref_transrate": "$REFDIR/transrate/mRNA/contigs.csv",
    "output_dir": "$REFDIR/QuantEval/"
}
EOF
               python3 $BASEDIR/scripts/QuantEval.py --input $REFDIR/QuantEval.json --reference
           else
               cat > $SEQDIR/QuantEval.json << EOF
{
    "ref_read_pickle": "$REFDIR/QuantEval/ref.pickle",
    "contig_fasta": "$SEQDIR/${SEQ}.fasta",
    "contig_blastn": "$SEQDIR/blastn/self.tsv",
    "contig_xprs_file": ["$SEQDIR/kallisto/abundance.tsv",
                         "$SEQDIR/rsem/rsem.isoforms.results",
                         "$SEQDIR/salmon/quant.sf"],
    "contig_xprs_label": ["kallisto", "rsem", "salmon"],
    "contig_xprs_header": [true, true, true],
    "contig_xprs_name_col": [1, 1, 1], 
    "contig_xprs_tpm_col": [5, 6, 4], 
    "contig_xprs_count_col": [4, 5, 5],
    "contig_transrate": "$SEQDIR/transrate/$SEQ/contigs.csv",
    "match_blastn": "$SEQDIR/blastn/contig_to_mRNA.tsv",
    "output_dir": "$SEQDIR/QuantEval/"
}
EOF
                python3 $BASEDIR/scripts/QuantEval.py --input $SEQDIR/QuantEval.json --reference --contig --match
            fi
        done
    done
done
