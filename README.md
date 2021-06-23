# t2g SNP [Version 1.0.0]

Script to convert transcriptome positioned SNPs to genomically positioned. 

**t2g** will convert transcriptome SNP positions to genomic positions without the need for remapping raw reads and variant calling. In a nutshell the script takes transcriptome SNP positions and BLAST aligns 100bp right & left flanking sequences to the genome to identify the corresponding genomic positions. As the genome reference may contain haplotigs (2+ scaffolds from the same region that correspond to the different haplotypes), the complete transcriptome contig is BLAST aligned to the genome to identify the relevant haplotigs where there is >=80% coverage (of the transcriptome contig). Flanking sequences matching are only considered within these scaffolds/contigs. Where multiple haplotigs exist the position within both are reported. Positions where the left and right flanking seq do not directly flank the SNP position in the genome are filtered out. Positions where there are 3 or more haplotigs are deemed unreliable and are also filtered out.

## Getting Started

### Requirements:
- BEDTools
- BLAST+
- R
- R packages:
	- tidyverse
	- stringr
	- Biostrings
	- argparse

### Installing

```shell
git archive --remote=git@github.com:lpembleton/t2g.git master | tar -x
```

## Workflow
The Rscript script requires three input files:
- [SNP_LIST] SNP list as 1-based position with 4 columns (CHROM POS REF,ALT SNP_Name)
- [TRANSCRIPT_REF] Transcriptome reference as fasta
- [GENOME_REF] Genome reference as fasta

```shell
Rscript t2g-snp.R -s SNP_LIST -tr TRANSCRIPT_REF -gr GENOME_REF
```

This will output a text table with columns: 

|Column|Description|
|---|---|
|SNP_Name|name of SNP from the input SNP list|
|g_CHROM|genome scaffold where the SNP has been positioned|
|bed_start|bed start (0-based) position of the SNP|
|bed_end|bed end position of the SNP|
|Comment|flanking seq info i.e. dual, left or right|
|g_hits|number of quality genomic hits|
|t_CHROM|transcriptome contig where the SNP is positioned|
|t_POS|transcriptome SNP position (1-based)|
|t_REF,t_ALT|transcriptome REF and ALT allele|

An R workspace image will also be saved to allow easier manual inspection and editing of results

## Authors
- L.W. Pembleton


