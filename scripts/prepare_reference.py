#!/usr/bin/env python3

import os
import re
import sys

def usage():
    if len(sys.argv) != 4:
        print('usage: ', os.path.basename(__file__), ' [reference directory] [simulation directory] [length threshold]') 
        sys.exit(1)
    return

def main():
    usage()

    annotation_pool = list()
    chromosome_pool = set()
    transcript_pool = set()
    gene_pool = set()
    mRNA_pool = set()
    mRNA_parent = dict()

    refdir = os.path.abspath(sys.argv[1])
    outdir = os.path.abspath(sys.argv[2])

    shortest_length = int(sys.argv[3])

    if not os.path.exists(outdir + '/chromosome'):
        os.makedirs(outdir + '/chromosome')

    print('load transcriptome.fasta')
    transcript_length = 0
    total_transcript_length = 0
    with open(refdir + '/transcriptome.fasta', 'r') as transcript_file:
        for i, transcript_line in enumerate(transcript_file):
            if transcript_line[0] == '>':
                if i != 0 and transcript_length >= shortest_length:
                    transcript_pool.add(transcript_name)

                transcript_length = 0
                regex = re.match('>(\S+)', transcript_line)
                transcript_name = regex.group(1)
                regex = re.match('(\S+)\.', transcript_name)
                if regex:
                    transcript_name = regex.group(1)
            else:
                transcript_length += len(transcript_line.rstrip())
                total_transcript_length += len(transcript_line.rstrip())
        
    print('    - number of transcripts: {:>16d}'.format(len(transcript_pool)))
    print('    - total length of transcripts: {:>10d}'.format(total_transcript_length))

    total_transcript_length = 0
    print('extract mRNA in annotation.gff')
    with open(refdir + '/annotation.gff', 'r') as annotation_file:
        for annotation_line in annotation_file:
            if annotation_line[0] == '#':
                continue
            annotation_line = annotation_line.rstrip()
            annotation_data = annotation_line.split('\t')
   
            chromosome = annotation_data[0]
            category = annotation_data[2]
            attributes = annotation_data[8]
           
            if chromosome in ['MT', 'Mito']:
                continue

            if category not in ['gene', 'mRNA', 'exon']:
                continue

            annotation_pool.append(annotation_data)

            if category == 'mRNA':
                regex = re.match('ID=transcript:(\S+?);', attributes)
                transcript_name = regex.group(1)
                regex = re.search('Parent=gene:(\S+?);', attributes)
                gene_name = regex.group(1)
                mRNA_parent[transcript_name] = gene_name

            elif category == 'exon':
                regex = re.match('Parent=transcript:(\S+?);', attributes)
                transcript_name = regex.group(1)
                if transcript_name in transcript_pool and transcript_name in mRNA_parent:
                    total_transcript_length += int(annotation_data[4]) - int(annotation_data[3]) + 1
                    chromosome_pool.add(chromosome)
                    gene_pool.add(mRNA_parent[transcript_name])
                    mRNA_pool.add(transcript_name)

    print('    - number of chromosomes: {:>16d}'.format(len(chromosome_pool)))
    print('    - number of genes: {:>22d}'.format(len(gene_pool)))
    print('    - number of transcripts: {:>16d}'.format(len(mRNA_pool)))
    print('    - total length of transcripts: {:>10d}'.format(total_transcript_length))

    print('output flux_simulator.gtf and mRNA.gtf')
    mRNA_count = 0
    total_transcript_length = 0
    annotation_out = open(refdir + '/mRNA.gtf', 'w')
    flux_annotation_out = open(outdir + '/flux_simulator.gtf', 'w')
    for annotation_data in annotation_pool:
        category = annotation_data[2]
        attributes = annotation_data[8]
        
        if category == 'gene':
            regex = re.match('ID=gene:(\S+?);', attributes)
            gene_name = regex.group(1)
            if gene_name in gene_pool:
                new_attributes = 'gene_id "' + gene_name + '";'
                print('\t'.join(annotation_data[:-1]) + '\t' + new_attributes, file=annotation_out)

        elif category == 'mRNA':
            regex = re.match('ID=transcript:(\S+?);', attributes)
            transcript_name = regex.group(1)
            regex = re.search('Parent=gene:(\S+?);', attributes)
            gene_name = regex.group(1)
            if transcript_name in mRNA_pool and gene_name in gene_pool:
                mRNA_count += 1
                new_attributes = 'gene_id "' + gene_name + '"; transcript_id "' + transcript_name + '";'
                print('\t'.join(annotation_data[:-1]) + '\t' + new_attributes, file=annotation_out)

        elif category == 'exon':
            regex = re.match('Parent=transcript:(\S+?);', attributes)
            transcript_name = regex.group(1)
            if transcript_name in mRNA_pool and transcript_name in mRNA_parent: 
                gene_name = mRNA_parent[transcript_name]
                if gene_name not in gene_pool:
                    continue
                total_transcript_length += int(annotation_data[4]) - int(annotation_data[3]) + 1
                new_attributes = 'gene_id "' + gene_name + '"; transcript_id "' + transcript_name + '";'
                print('\t'.join(annotation_data[:-1]) + '\t' + new_attributes, file=annotation_out)
                print('\t'.join(annotation_data[:-1]) + '\t' + new_attributes, file=flux_annotation_out)
    
    print('    - number of chromosomes: {:>16d}'.format(len(chromosome_pool)))
    print('    - number of genes: {:>22d}'.format(len(gene_pool)))
    print('    - number of transcripts: {:>16d}'.format(mRNA_count))
    print('    - total length of transcripts: {:>10d}'.format(total_transcript_length)) 
    annotation_out.close()
    flux_annotation_out.close()

    print('output genome.fasta and chromosome.fasta')
    chromosome_count = 0
    retain_chromosome = False
    genome_in = open(refdir + '/genome.fasta', 'r')
    genome_out = open(refdir + '/chromosome.fasta', 'w')
    for i, genome_line in enumerate(genome_in):
        if genome_line[0] == '>':
            try:
                chromosome_out.close()
            except:
                pass
            regex = re.match('^>(\S+)', genome_line)
            genome_name = regex.group(1)
            if genome_name in chromosome_pool:
                chromosome_count += 1
                retain_chromosome = True
                chromosome_out = open(outdir + '/chromosome/' + genome_name + '.fa', 'w')
            else:
                retain_chromosome = False        
        if retain_chromosome:
            print(genome_line, file=chromosome_out, end='')
            print(genome_line, file=genome_out, end='')
    genome_out.close()
    print('    - number of chromosomes: {:>16d}'.format(chromosome_count))
    
    print('output mRNA.fasta')
    mRNA_count = 0
    total_transcript_length = 0
    retain_transcript = False
    transcript_in = open(refdir + '/transcriptome.fasta', 'r')
    transcript_out = open(refdir + '/mRNA.fasta', 'w')
    for transcript_line in transcript_in:
        if transcript_line[0] == '>':
            regex = re.match('^>(\S+)', transcript_line)
            transcript_name = regex.group(1)
            regex = re.match('(\S+)\.', transcript_name)
            if regex:
                transcript_name = regex.group(1)
            if transcript_name in mRNA_pool: 
                mRNA_count += 1
                retain_transcript = True
                print('>' + transcript_name, file=transcript_out)
            else:
               retain_transcript = False
        elif retain_transcript == True:
            total_transcript_length += len(transcript_line.strip())
            print(transcript_line, file=transcript_out, end='')
    print('    - number of transcripts: {:>16d}'.format(mRNA_count))
    print('    - total length of transcripts: {:>10d}'.format(total_transcript_length)) 
    transcript_in.close()
    transcript_out.close()
    return

if __name__ == '__main__':
    main()
    sys.exit(0)
