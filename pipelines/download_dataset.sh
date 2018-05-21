#!/usr/bin/env bash

set -e

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$SCRIPTDIR/../reference"
BASEDIR="$SCRIPTDIR/../reference"
REALDIR="$SCRIPTDIR/../real/"

mkdir -p "$BASEDIR/dog/"
cd "$BASEDIR/dog/"
wget ftp://ftp.ensembl.org/pub/release-90/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.dna.toplevel.fa.gz -O genome.fasta.gz
wget ftp://ftp.ensembl.org/pub/release-90/fasta/canis_familiaris/cdna/Canis_familiaris.CanFam3.1.cdna.all.fa.gz -O transcriptome.fasta.gz
wget ftp://ftp.ensembl.org/pub/release-90/gff3/canis_familiaris/Canis_familiaris.CanFam3.1.90.chr.gff3.gz -O annotation.gff.gz
gzip -d *

mkdir -p "$BASEDIR/mouse/"
cd "$BASEDIR/mouse/"
wget ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz -O genome.fasta.gz
wget ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz -O transcriptome.fasta.gz
wget ftp://ftp.ensembl.org/pub/release-90/gff3/mus_musculus/Mus_musculus.GRCm38.90.chr.gff3.gz -O annotation.gff.gz
gzip -d *

mkdir -p "$BASEDIR/yeast/"
cd "$BASEDIR/yeast/"
wget ftp://ftp.ensembl.org/pub/release-90/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz -O genome.fasta.gz
wget ftp://ftp.ensembl.org/pub/release-90/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -O transcriptome.fasta.gz
wget ftp://ftp.ensembl.org/pub/release-90/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.90.gff3.gz -O annotation.gff.gz
gzip -d *

mkdir -p "$REALDIR/dog/reads"
cd "$REALDIR/dog/reads"
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR882109
ln SRR882109_1.fastq raw_r1.fastq 
ln SRR882109_2.fastq raw_r2.fastq

mkdir -p "$REALDIR/mouse/reads"
cd "$REALDIR/mouse/reads"
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR203276
ln SRR203276_1.fastq raw_r1.fastq
ln SRR203276_2.fastq raw_r2.fastq

mkdir -p "$REALDIR/yeast/reads"
cd "$REALDIR/yeast/reads"
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR453566
ln SRR453566_1.fastq raw_r1.fastq
ln SRR453566_2.fastq raw_r2.fastq


