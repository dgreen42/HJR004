# HJR004 pipeline

### Required software:
 - minimap2
 - samtools
 - R
 - R packages:
    - Bambu
    - edgeR
    - stageR

### Plotting functions

There are extra plotting functions in the utility script (bin/utils.R) that can be used to visualize results from the pipeline. Inputs for these functions come from the output of the pipeline. The plotting functions use the following inputs from the listed outputs:

 - propTable: analysis/dex_isoform_proportions.csv
 - gene: cts/bambu_out/HJR4_counts_transcript.txt; field: GENEID
 - annotation: annotation from resources
 - acronym_list: acronym list from resources

### Data Sources

 - Genome and Annotation: https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/downloads/1.9/20220708_MtrunA17r5.0-ANR-EGN-r1.9_Annotation.zip
 - Gene Acronyms: https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/downloads/1.9/Mt_gene_Acronyms_IDs_latest.xlsx
 - RNA seq data: N/A (will exist after publication)

### Example run

```
nextflow run main.nf \
    --fastq 'fastq_dir' \
    --ref_annotation 'annotation.gtf' \
    --ref_genome 'genome.fasta' \
    --sample_sheet 'samples.csv'
```
