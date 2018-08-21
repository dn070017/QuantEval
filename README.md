 # QuantEval
 ## About
 QuantEval was released for two purposes: 1. a user can use the scripts in QuantEval to reproduce all the analyses in the following study (Hsieh et al.); 2. a user can follow the example in QuantEval to conduct the same analyses on his own study. For the first purpose, there are three modes in the QuantEval main program. (1) <b>Reference Mode</b>, (2) <b>Contig Mode</b> and (3) <b>Match Mode</b>. The first two modes read the quantification results and build a ambiguity cluster based on connected components for the reference transcripts and contig sequences. The match mode builds relations between contigs and reference transcripts. For the second purpose, the users are encouraged to follow the provided example to conduct new analyses on his own data.
 ## Reference
 > Ping-Han Hsieh, Yen-Jen Oyang and Chien-Yu Chen. Effect of de novo transcriptome assembly on transcript quantification. 2018, bioRxiv 380998.
 ## Requirement
 - QuantEval main program:
    - Python3 (3.5.2)
    - Python packages: pandas (0.20.3), numpy (1.12.1)
 - Generate figures and table:
    - R (3.3.0)
    - R pacakges: gridExtra, grid, stats, tidyverse, plyr, ggplot2, reshape2
 - Utilities:
    - Bowtie2 (2.3.0), BLASTn (2.5.0), Flux Simulator (1.2.1), RSEM (1.2.31), Kallisto (0.43.0), rnaSPAdes (3.11.1), Salmon (0.8.2), Trans-ABySS (1.5.5), TransRate (1.0.3), Trinity (2.4.0)

 ## Manual
 - Run QuantEval individually:
 ```shell
 python3 ./scripts/QuantEval.py --reference --contig --match --input input.json
 ```
 the first three parameters (<b>--reference, --contig, --match</b>) indicate which mode to run (the three mode can be run independantly, but one has to run both reference mode and contig mode <b>before</b> running match mode, it is recommended to run three mode together) and the <b>input.json</b> file specify the input parameters for the QuantEval main program. Because the main program of QuantEval <b>does not</b> include a wrapper for quantification/sequence alignment/contig evaluation, which are essesntial steps for QuantEval main program, one might need to run quantification algorithms (i.e. RSEM/Kallisto/Salmon), sequence alignment (BLASTn) and contig evaluation (Transrate) by themself in order to get similar analysis result in the reference research.
 - Run pairwise BLASTn for reference/contig mode:
 ```shell
 # reference mode
 blastn -db ref.fasta -query ref.fasta -outfmt 6 -evalue 1e-5 -perc_identity 95 -out ./blastn/ref.self.tsv 
 
 # contig mode
 blastn -db contig.fasta -query contig.fasta -outfmt 6 -evalue 1e-5 -perc_identity 95 -out ./blastn/contig.self.tsv
 ```
 - Run BLASTn for the mapping of reference and contig sequence (match mode)
 ```shell
 blastn -db ref.fasta -query contig.fasta -outfmt 6 -out ./blastn/contig_to_ref.tsv 
 ```
 - Run quantification for reference/contig mode with default parameters:
 ```shell
 # RSEM
 rsem-prepare-reference --bowtie2 ref.fasta ./rsem/rsem.index
 rsem-calculate-expression --paired-end --strandedness none --bowtie2 --time ref_read/read_1.fastq ref_read/read_2.fastq ./rsem/rsem.index ./rsem/rsem 

 # Kallisto
 kallisto index -i ./kallisto/kallisto.index -k 31 ref.fasta
 kallisto quant -i ./kallisto/kallisto.index -o ./kallisto ref_read/read_1.fastq ref_read/read_2.fastq

 # Salmon
 salmon index -i ./salmon/salmon.index -t ref.fasta --type quasi -k 31
 salmon quant -i ./salmon/salmon.index -l A -1 ref_read/read_1.fastq -2 ref_read/read_2.fastq -o ./salmon 
 ```
 - Run TransRate for reference/contig mode: 
 ```
 # reference mode
 transrate --assembly ref.fasta --output ./transrate/ref --left ref_read/read_1.fastq --right ref_read/read_2.fastq

 # contig mode
 transrate --assembly contig.fasta --output ./transrate/ref --left contig_read/read_1.fastq --right contig_read/read_2.fastq
 ```
 - Example of input.json:
 ```json
{
    "ref_fasta": "ref.fasta",
    "ref_blastn": "./blastn/ref.self.tsv",
    "ref_gtf": "ref.gtf",
    "ref_xprs_file": ["./answer/answer_xprs.tsv",
                      "./kallisto/ref/abundance.tsv",
                      "./rsem/ref/rsem.isoforms.results",
                      "./salmon/ref/quant.sf"],
    "ref_xprs_label": ["answer", "kallisto", "rsem", "salmon"],
    "ref_xprs_header": [true, true, true, true],
    "ref_xprs_name_col": [1, 1, 1, 1], 
    "ref_xprs_tpm_col": [2, 5, 6, 4], 
    "ref_xprs_count_col": [3, 4, 5, 5],
    "ref_transrate": "./transrate/ref/contigs.csv",
    "contig_fasta": "contig.fasta",
    "contig_blastn": "./blastn/contig.self.tsv",
    "contig_xprs_file": ["./kallisto/contig/abundance.tsv",
                         "./rsem/contig/rsem.isoforms.results",
                         "./salmon/contig/quant.sf"],
    "contig_xprs_label": ["kallisto", "rsem", "salmon"],
    "contig_xprs_header": [true, true, true],
    "contig_xprs_name_col": [1, 1, 1], 
    "contig_xprs_tpm_col": [5, 6, 4], 
    "contig_xprs_count_col": [4, 5, 5],
    "contig_transrate": "./transrate/contig/contigs.csv",
    "match_blastn": "./blastn/contig_to_ref.tsv",
    "output_dir": "./QuantEval/"
}
 ```
- Run example:
```
cd example
python3 ../scripts/QuantEval.py --reference --contig --match --input ./example.json
```
One can also import the functions in utilities.py to built their own analysis pipeline.
- Construct connected component for reference only:
```python
from utilities import construct_sequence, filter_blastn, intersect_match, construct_grap
import copy

input_file = dict()
input_file["contig_ref_file"] = ["./answer/answer_xprs.tsv",
                                 "./kallisto/abundance.tsv",
                                 "./rsem/rsem.isoforms.results",
                                 "./salmon/quant.sf"]
input_file["ref_xprs_label"]: ["answer", "kallisto", "rsem", "salmon"],
input_file["ref_xprs_header"]: [true, true, true, true],
input_file["ref_xprs_name_col"]: [1, 1, 1, 1], 
input_file["ref_xprs_tpm_col"]: [2, 5, 6, 4], 
input_file["ref_xprs_count_col"]: [3, 4, 5, 5],   
ref_seq_dict = construct_sequence('ref.fasta')
ref_self_blastn = filter_blastn('./blastn/ref.self.tsv')
read_expression(input_file, ref_seq_dict, 'ref')
ref_self_match_dict = intersect_match(ref_self_blastn, ref_seq_dict, copy.deepcopy(ref_seq_dict))
ref_uf, ref_component_dict = construct_graph(ref_seq_dict, ref_self_match_dict)

print(ref_uf.component_label)
print(ref_uf.parent)
```
- Run all the analysis in the study (<b>time consuming</b>):
```shell
./pipelines/run_analysis.sh
```
- Output format

| column | description |
|--------|-------------|
| match_name | alignment (ref.contig.strand) |
| contig_name | target contig name |
| ref_name | target ref name |
| accuracy | accuracy of the alignment |
| recovery | recovery of the alignment |
| ***contig/ref***\_length | length of ***contig/ref*** |
| ***contig/ref***\_tr\_***transrate_score*** | ***transrate score*** of ***contig/ref*** |
| ***contig/ref***\_xprs\_***tpm/count***_***quantifier*** | quantification result of ***contig/ref*** |
| ***contig/ref***\_component | label of connected component of ***contig/ref*** |
| ***contig/ref***\_component_size | number of sequences in the connected component of ***contig/ref*** |
| ***contig/ref***\_component_contribute_xprs_***tpm/count***\_***quantifier*** | proportion of ***TPM/read count (RPEA)*** in the connected component of ***contig/ref*** |
| ***contig/ref***\_component_relative_xprs_***tpm/count***\_***quantifier*** | ***TPM/count*** of ***contig/ref*** / highest ***TPM/count*** in the same connected component |
| ***contig/ref***\_component_max_xprs_***tpm/count***\_***quantifier*** | highest ***TPM/count*** of ***contig/ref*** in the same connected component |
| ***contig/ref***\_component_avg_xprs_***tpm/count***\_***quantifier*** | average ***TPM/count*** of ***contig/ref*** in the same connected component |
| ***contig/ref***\_component_tot_xprs_***tpm/count***\_***quantifier*** | total ***TPM/count*** of ***contig/ref*** in the same connected component |
| ref_gene_contribute_xprs_***tpm/count***\_***quantifier*** | proportion of ***TPM/read count*** in the gene of ref |
| ref_gene_relative_xprs_***tpm/count***\_***quantifier*** | ***TPM/count*** of ref / highest ***TPM/count*** in the same gene |
| ref_gene_max_xprs_***tpm/count***\_***quantifier*** | highest ***TPM/count*** of ref in the same gene |
| ref_gene_avg_xprs_***tpm/count***\_***quantifier*** | average ***TPM/count*** of ref in the same gene |
| ref_gene_tot_xprs_***tpm/count***\_***quantifier*** | total ***TPM/count*** of ref in the same gene |
| length_difference | the difference of length between contig and reference |
| xprs_***tpm/count***\_error_***quantifier*** | quantificaion error for the estimated abundance of contig | 
