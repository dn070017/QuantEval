 # QuantEval
 ## About
 QuantEval is an analsis pipeline which evaluate the reliability of quantification tools. There are three modes in the QuantEval main program. (1) <b>Reference Mode</b>, (2) <b>Contig Mode</b> and (3) <b>Match Mode</b>. The first two modes read the quantification results and build a ambiguity cluster based on connected components for the reference transcripts and contig sequences. The match mode built relations between contigs and reference transcripts.
 ## Reference
 > Ping-Han Hsieh, Yen-Jen Oyang and Chien-Yu Chen. Effect of de novo transcriptome assembly on quality of read mapping and transcript quantification. 2018, bioRxiv 380998.
 ## Manual
 - Run all the analysis in the study:
 ```shell
 ./pipelines/run_analysis.sh
 ```
 - Run QuantEval indivisually:
 ```shell
 ./scripts/QuantEval.py --reference --contig --match --input input.json
 ```
 the first three parameters indicate which mode to run, the three mode can be run independantly and the input.json file specify the input parameters for the QuantEval main program. Here's an example of input.json for reference mode:
 ```json
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
 ```
 and here's an example for contig and match mode:
```json
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
```
note that one has to run reference mode and contig mode before running match mode. One can also import the functions in utilities.py to built their own analysis pipeline.

Because the main program of QuantEval does not include a wrapper for quantification/sequence alignment/contig evaluation, which are essesntial steps for QuantEval main program, one might need to use the same parameters in the modules located in pipelines (assembly.sh, quantification.sh, postprocessing.sh) and run quantification algorithms (i.e. RSEM/Kallisto/Salmon) and sequence alignment (BLASTn) and contig evaluation (Transrate) by themself in order to get similar analysis result in the reference research.