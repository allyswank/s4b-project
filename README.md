# Analysis of Gene Expression Differentiation Following Environmental Stressors

Next generation sequencing has provided opportunities to understand the molecular and genomic consequences of climate change induced environmental stressors on organisms. Specifically, RNA-seq has enabled us to understand what genes are being up- or down-regulated in response to these stressors. By analyzing gene expression data, we have a much greater understanding of the physiological consequences that have been previously described in many model species. 


## Project Objectives

In this project, I intend to develop scripts that will:
* Clean RNA-seq reads and align them to a reference genome (bash)
* Find differentially expressed genes (DEGs) (R)

These analyses will consist of RNA-seq data extracted from the common mouse species, *Mus musculus*, that was exposed to two separate treatment groups. 

# All Output Files

Throughout these scripts, many outputs are generated. The following section contains descriptions of what will be contained in each of these files. 

## FastQC (Fast Quality Checks)
* **FILENAME_fastqc.html** : Contains all quality check information. Having good read quality is essential prior to further analysis.
* FILENAME_fastqc.zip : a compress version of the .html file 

## TrimGalore
* lessl

## STAR (Spliced Transcripts Alignement to a Reference)
* **Aligned.out.sam** : Alignments to the genome in standard SAM format.
* **Log.out** : Main log file with a lot of detailed information about the run. This file is most useful for troubleshooting and debugging.
* **Log.progress.out** : Reports job progress statistics, such as the number of processed reads, % of mapped reads etc.
* **Log.final.out** :  Summary mapping statistics after mapping job is complete, very useful for quality control. The statistics are calculated for each read (single- or paired-end) and then summed or averaged over all reads.
* **ReadsPerGene.out.tab** : 4 columns - 1) gene ID, 2) counts for unstranded RNA-seq, 3) counts for the 1st read strand aligned with RNA, 4) counts for the 2nd read strand aligned with RNA.
* SJ.out.tab : Contains high confidence collapsed splice junctions in tab-delimited format. 
* chrLength.txt : 1 column - length in base pairs of each chromosome
* chrNameLength.txt : 2 columns - 1) chromosome name and 2) chromosome length
* **chrName.txt** : 1 column - chromosome names. After indexing the genome, you can change the chromosome names in this file while keeping the same order. The names in this file will be used in all output alignment files (i.e. Alignment.out.sam).
* chrStart.txt : 1 column - starting point of each chromosome in order
* exonGeTrInfo.tab : position of exons (based on transcript file [.gtf])
* exonInfo.tab : length of exons
* geneInfo.tab : empty
* **Genome** : Comprise binary genome sequence, suffix arrays, text chromosome names/lengths, splice junctions coordinates, and transcripts/genes information. Most of these files use internal STAR format and are not intended to be utilized by the end user.
* genomeParameters.txt : User input parameters for STAR to run, including defaults
* **SA** : Suffix Array (a sorted array of all suffixes of a string)
* **SAindex** : Indexed representation of Genome
* sjdbInfo.txt : transcript length (set by user, based on sequencing adapters, Illumina is 99)
* sjdbList.fromGTF.out.tab : empty
* sjdbList.out.tab : empty
* transctiptInfo.tab : all transcripts 


## Box Folder

Here is the link to the Box Folder with our data: https://auburn.box.com/s/x17154cah63h3yp4wd14ak6p0fa3t5jz
