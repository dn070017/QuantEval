#!/usr/bin/env bash

set -e

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$SCRIPTDIR/../reference"
BASEDIR="$SCRIPTDIR/../reference"

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
