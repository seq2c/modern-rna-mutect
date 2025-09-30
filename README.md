# Modern RNA-MuTect Pipeline (GATK4 + Mutect2)

## Introduction

This repository provides a robust Bash script that modernizes and automates the workflow from the original [RNA_MUTECT_1.0-1](https://github.com/broadinstitute/RNA_MUTECT_1.0-1/tree/master) by [Yizhak et al 2019](https://www.science.org/doi/10.1126/science.aaw0726). The original pipeline is a powerful method for somatic variant detection in RNA-Seq, but it is challenging to use today as it relies on outdated tools (GATK3, MuTect1), an old reference genome (hg19), and unmaintained MATLAB executables.

This script replicates the core logic of the pipeline—initial discovery, annotation, targeted re-alignment, and final re-calling—using modern, maintained tools like **GATK4** (including **Mutect2**,**Funcotator**, etc.), and **HISAT2** on the **GRCh38/hg38** human reference genome. The output is a high-confidence VCF file, ready for the further filtering steps.

Modifications:
* Allows for tumor-only RNA-seq inputs in addition to matched tumor/normal pairs
* Runs computationally intensive steps in parallel by splitting work across chromosomes
* No more hardcoded hg19

## 1\. Tool installation

Tools:
- gatk
- samtools
- hisat2
- gsutil

**GATK**

Download the prebuilt gatk from [github](https://github.com/broadinstitute/gatk/releases) and unzip it. For example here it is the newest version of gatk 4.6.2
```
wget https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip
unzip gatk-4.6.2.0.zip
````
then add the `/path/to/gatk-4.6.2.0/` to $PATH in ~/.bashrc

Next we use [miniconda](https://docs.conda.io/en/latest/miniconda.html) to install required packages as this is straightforward.
```
# create a new environment
conda create -n mutect python=3.13 -y
conda activate mutect

# configure the necessary channels 
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# install
conda install samtools hisat2 google-cloud-sdk -y
```


## 2\. Reference file preparation

This pipeline requires several reference files. It is easier that all files are based on the same reference genome. Mixing different builds (e.g., UCSC vs. Ensembl) would cause contig naming errors.

Here we use GRCh38 as an example.

To run RNA-mutect we need:
* Reference genome fasta
* Reference gene annotation
* [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial#1.1.2) somatic data source from gatk
* Germline variant source (gnomad) from gatk
* Panel of normal (1000 genome project) from gatk


**Download reference files**

```bash
cd /path/to/reference
# reference genome
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.fa.gz
# gtf
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.ensGene.gtf.gz

```

For data from gatk we can simply use google cloud sdk, but you gotta set up a google cloud account first (usually with free trial)

```bash
# add account to gsutil
gcloud init
gcloud auth login

gsutil -u [project id] cp gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz .

gsutil -u [project id] cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz .

gsutil -u [project id] cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz .

```

**Prepare and index files**

```bash
gunzip hg38.fa.gz
samtools faidx hg38.fa
gatk CreateSequenceDictionary -R hg38.fa
gunzip hg38.ensGene.gtf.gz

tabix -p vcf af-only-gnomad.hg38.vcf.gz
tabix -p vcf 1000g_pon.hg38.vcf.gz

tar -zxvf funcotator_dataSources.v1.8.hg38.20230908s.tar.gz

# build HISAT2 index (This could be a long step)
mkdir -p hisat2_index
# 1. extract splice sites and exons from the GTF
hisat2_extract_splice_sites.py hg38.ensGene.gtf > hg38.ss
hisat2_extract_exons.py hg38.ensGene.gtf > hg38.exon

# 2. build the index
hisat2-build -p 16 --ss hg38.ss --exon hg38.exon hg38.fa ./hisat2_index/hg38
```

## 3\. Run the pipeline
In the original RNA mutect paper, they recommended align the RNA-seq data using STAR. They filtered duplicated reads in a more downstream step, however, we can simply deduplicate RNA-seq reads before running variant calling. Here is an example with `umi_tools`, or with any other deduplicate tool.

First extract umis
```bash
conda install umi_tools

sample=SAMPLE_ID
# umi_tools extract
umi_tools extract -I ${read1}.fastq.gz \
  --read2-in=${read2}.fastq.gz \
  -S ${sample}.umi_extract_1.fastq.gz \
  --read2-out=${sample}.umi_extract_2.fastq.gz \
  --bc-pattern="NNNCC" --bc-pattern2="NNNCC" # you can figure out if your reads have umis and if so, what the patterns are.
```

then run star
```bash
# build STAR index
STAR --runMode genomeGenerate \
     --runThreadN 16 \
     --genomeDir /path/to/star_index/ \
     --genomeFastaFiles /path/to/hg38.fa \
     --sjdbGTFfile /path/to/hg38.ensGene.gtf \
     --sjdbOverhang 99  # Should be ReadLength - 1. 99 is a common value for 100bp reads.

# run star
STAR --runThreadN 16 \
     --genomeDir /path/to/star_index/ \
     --readFilesIn /path/to/${sample}.umi_extract_1.fastq.gz /path/to/${sample}.umi_extract_2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMattrRGline ID:${sample} SM:${sample} PL:Illumina \
     --outFileNamePrefix /path/to/rna_ \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes Standard \
     --outFilterMultimapNmax 50 \
     --peOverlapNbasesMin 10 \
  	 --alignSplicedMateMapLminOverLmate 0.5 \
  	 --alignSJstitchMismatchNmax 5 -1 5 5 \
     --twopassMode Basic
```
then deduplicate with umi_tools
```bash
# umi_tools dedup
umi_tools dedup -I ${sample}.Aligned.sortedByCoord.out.bam \
  -S ${sample}.dedup.bam -L ${sample}.dedup.log \
  --output-stats ${sample}.dedup --paired --random-seed=100
```

### Tumor-only mode
```bash
# rna.bam could be ${sample}.dedup.bam
./run_modern_rna_mutect.sh \
    --rna-bam /path/to/rna.bam \
    --genome-fasta /path/to/hg38.fa \
    --pon /path/to/1000g_pon.hg38.vcf.gz \
    --germline-resource /path/to/af-only-gnomad.hg38.vcf.gz \
    --funcotator-sources /path/to/funcotator_dataSources.v1.8.hg38.20230908s \
    --hisat2-index /path/to/hisat2_index/hg38 \
    --output-dir /path/to/output \
    --threads 16
```

### Matched-normal mode
```bash
./run_modern_rna_mutect.sh \
    --rna-bam /path/to/rna.bam \
    --normal-bam /path/to/dna_normal.bam \
    --normal-sample-name NORMAL_SAMPLE_NAME_IN_BAM \
    --genome-fasta /path/to/hg38.fa \
    --pon /path/to/1000g_pon.hg38.vcf.gz \
    --germline-resource /path/to/af-only-gnomad.hg38.vcf.gz \
    --funcotator-sources /path/to/funcotator_dataSources.v1.8.hg38.20230908s \
    --hisat2-index /path/to/hisat2_index/hg38 \
    --output-dir /path/to/output \
    --threads 16
```

If you look at the script, it actually includes the following:

* Runs `SplitNCigarReads`, `Mutect2`, and `Funcotator` in parallel for speed
* Follows the paper’s approach: extract site reads and re-align with HISAT2 for targeted checks
* Supports both tumor-only and matched-normal modes


## To do
* Rewrite SNV filtering scripts
* Allows BAMs from various aligners
