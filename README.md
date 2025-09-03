# HJR004 pipeline

### Required software:
 - minimap2
 - samtools
 - R
 - R package Bambu

### Example run

```
nextflow run main.nf \
    --fastq 'fastq_dir' \
    --ref_annotation 'annotation.gtf' \
    --ref_genome 'genome.fasta' \
    --sample_sheet 'samples.csv'
```
