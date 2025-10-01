#!/bin/bash

# =============================================================================
# Modern RNA-MuTect Pipeline (GATK4 + HISAT2)
# Version 1.0
#
# This script modernizes the RNA-MuTect v1.0 pipeline. It takes an
# RNA BAM (STAR-aligned) and an optional matched-normal DNA BAM and produces
# a re-aligned, re-called VCF ready for the further filtering steps.
# =============================================================================

set -e
set -o pipefail

# args
usage() {
    echo "Usage: $0 -r <RNA_BAM> [-n <NORMAL_DNA_BAM> -s <NORMAL_SAMPLE_NAME>] -g <GENOME_FASTA> -p <PANEL_OF_NORMALS> -a <GERMLINE_RESOURCE> -f <FUNCOTATOR_SOURCES> -h <HISAT2_INDEX> -o <OUTPUT_DIR> [-t <THREADS>]"
    echo ""
    echo "Required:"
    echo "  -r, --rna-bam                Path to the input RNA BAM file (STAR-aligned)."
    echo "  -g, --genome-fasta           Path to the reference genome FASTA file (hg38/GRCh38 with 'chr' prefixes)."
    echo "  -p, --pon                    Path to the Panel of Normals VCF file."
    echo "  -a, --germline-resource      Path to the Germline Resource VCF (e.g., af-only-gnomad)."
    echo "  -f, --funcotator-sources     Path to the Funcotator data sources directory."
    echo "  -h, --hisat2-index           Path and prefix for the HISAT2 index (must match the genome)."
    echo "  -o, --output-dir             Path to the main output directory for all results."
    echo ""
    echo "Matched-Normal Mode (optional):"
    echo "  -n, --normal-bam             Path to the matched-normal DNA BAM file."
    echo "  -s, --normal-sample-name     Sample name (SM tag) of the normal sample in the BAM header."
    echo ""
    echo "Optional:"
    echo "  -t, --threads                Number of parallel jobs to run (default: 8)."
    echo ""
    exit 1
}

NUM_PARALLEL_JOBS=8
RNA_BAM_PATH=""; DNA_BAM_PATH=""; NORMAL_SAMPLE_NAME=""
REFERENCE_FASTA=""; PANEL_OF_NORMALS=""; GERMLINE_RESOURCE=""
FUNCOTATOR_DATA_SOURCES=""; HISAT2_INDEX=""; OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--rna-bam) RNA_BAM_PATH="$2"; shift;;
    -n|--normal-bam) DNA_BAM_PATH="$2"; shift;;
    -s|--normal-sample-name) NORMAL_SAMPLE_NAME="$2"; shift;;
    -g|--genome-fasta) REFERENCE_FASTA="$2"; shift;;
    -p|--pon) PANEL_OF_NORMALS="$2"; shift;;
    -a|--germline-resource) GERMLINE_RESOURCE="$2"; shift;;
    -f|--funcotator-sources) FUNCOTATOR_DATA_SOURCES="$2"; shift;;
    -h|--hisat2-index) HISAT2_INDEX="$2"; shift;;
    -o|--output-dir) OUTPUT_DIR="$2"; shift;;
    -t|--threads) NUM_PARALLEL_JOBS="$2"; shift;;
    -H|--help) usage;;
    *) echo "Unknown param: $1"; usage;;
  esac
  shift
done

# validate inputs
if [ -z "$RNA_BAM_PATH" ] || [ -z "$REFERENCE_FASTA" ] || [ -z "$PANEL_OF_NORMALS" ] || [ -z "$GERMLINE_RESOURCE" ] || [ -z "$FUNCOTATOR_DATA_SOURCES" ] || [ -z "$HISAT2_INDEX" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: One or more required arguments are missing."
    usage
fi

if [[ -n "$DNA_BAM_PATH" && -z "$NORMAL_SAMPLE_NAME" ]] || [[ -z "$DNA_BAM_PATH" && -n "$NORMAL_SAMPLE_NAME" ]]; then
    echo "Error: If using a matched normal, both -n/--normal-bam and -s/--normal-sample-name must be provided."
    usage
fi

echo "Starting RNA-MuTect pipeline..."
date

# 1. setup
echo "[STAGE 1/5] Setting up directories and variables..."
BASENAME=$(basename ${RNA_BAM_PATH} .bam)
mkdir -p ${OUTPUT_DIR}
REALIGN_DIR=${OUTPUT_DIR}/hisat2
mkdir -p ${OUTPUT_DIR} ${REALIGN_DIR}

GN=${REFERENCE_FASTA}
DICT=${GN%.*}.dict
if [[ ! -f "$DICT" ]]; then
  echo " -> Creating sequence dictionary..."
  gatk CreateSequenceDictionary -R ${GN}
fi

RNA_BAM_SPLIT=${OUTPUT_DIR}/${BASENAME}.split.bam
MERGED_VCF=${OUTPUT_DIR}/${BASENAME}.merged.vcf.gz
FUNCOMAF=${OUTPUT_DIR}/${BASENAME}.funcotated.maf
CONTIG_LIST=$(cut -f1 "${REFERENCE_FASTA}.fai")


# Export variables to be available in xargs subshells
export GN PANEL_OF_NORMALS GERMLINE_RESOURCE OUTPUT_DIR BASENAME 
export RNA_BAM_PATH RNA_BAM_SPLIT MERGED_VCF FUNCOMAF FUNCOTATOR_DATA_SOURCES


# 2. INITIAL VARIANT CALLING
echo "[STAGE 2/5] Performing initial variant discovery with Mutect2..."

# run SplitNCigarReads on each contig in parallel
printf "%s\n" ${CONTIG_LIST} | xargs -I {} -P ${NUM_PARALLEL_JOBS} bash -c '
    chr="{}"
    gatk SplitNCigarReads \
        -R "${GN}" \
        -I "${RNA_BAM_PATH}" \
        -L "${chr}" \
        -O "${OUTPUT_DIR}/${BASENAME}.${chr}.split.bam"
'
# GATHER splited bams
find ${OUTPUT_DIR} -name ${BASENAME}.*.split.bam > ${OUTPUT_DIR}/split_bam_list.txt
gatk GatherBamFiles \
    -I ${OUTPUT_DIR}/split_bam_list.txt \
    -O ${RNA_BAM_SPLIT} \
    -R ${GN}
samtools index -@ ${NUM_PARALLEL_JOBS} ${RNA_BAM_SPLIT}

# clean up
while read -r b; do 
   rm -f ${b} ${b%.bam}.bai
done < ${OUTPUT_DIR}/split_bam_list.txt
rm -f ${OUTPUT_DIR}/split_bam_list.txt


# parallel Mutect2
MUTECT2_CMD_BASE='gatk Mutect2 -R "${GN}" -I "${RNA_BAM_SPLIT}" --panel-of-normals "${PANEL_OF_NORMALS}" --germline-resource "${GERMLINE_RESOURCE}" -L "${chr}" -O "${OUTPUT_DIR}/${BASENAME}.${chr}.vcf.gz"'
if [[ -n "${DNA_BAM_PATH}" ]]; then
  echo " -> Matched-normal mode"
  export DNA_BAM_N="${DNA_BAM_PATH}" NORMAL_SAMPLE_NAME
  MUTECT2_CMD="${MUTECT2_CMD_BASE} -I \"\${DNA_BAM_N}\" -normal \"\${NORMAL_SAMPLE_NAME}\""
else
  echo " -> Tumor-only mode"
  MUTECT2_CMD="${MUTECT2_CMD_BASE}"
fi
export MUTECT2_CMD
printf "%s\n" ${CONTIG_LIST} | xargs -I {} -P "${NUM_PARALLEL_JOBS}" bash -c 'chr="{}"; eval "$MUTECT2_CMD"'

# merge VCFs
echo "--> Merging parallel VCF results..."
find ${OUTPUT_DIR} -name ${BASENAME}.*.vcf.gz > ${OUTPUT_DIR}/vcf_list.txt
gatk MergeVcfs -I ${OUTPUT_DIR}/vcf_list.txt -O ${MERGED_VCF}

# clean up
while read -r vcf_file; do 
    rm -f ${vcf_file} ${vcf_file}.tbi ${vcf_file}.stats
done < ${OUTPUT_DIR}/vcf_list.txt
rm -f ${OUTPUT_DIR}/vcf_list.txt


# annotation (NOTE: hg38 is hardcoded here)
echo "[STAGE 3/5] Annotating variants with Funcotator..."
# parallel Funcotator
printf "%s\n" ${CONTIG_LIST} | xargs -I {} -P ${NUM_PARALLEL_JOBS} bash -c '
  chr="{}"
  echo "    -> Annotating contig: ${chr}"
  gatk Funcotator \
    --variant "${MERGED_VCF}" --reference "${GN}" --ref-version hg38 \
    -L "${chr}" --data-sources-path "${FUNCOTATOR_DATA_SOURCES}" \
    --output "${OUTPUT_DIR}/${BASENAME}.${chr}.maf" --output-file-format MAF
'
# merge MAFs
echo "--> Merging parallel MAF results..."
MAF_FILES=$(find "${OUTPUT_DIR}" -name "${BASENAME}.*.maf" -size +0)
if [ -z "${MAF_FILES}" ]; then
    echo "FATAL: No MAF files were generated by Funcotator." && exit 1;
fi
FIRST_MAF=$(echo "${MAF_FILES}" | head -n 1)
grep "^Hugo_Symbol" "${FIRST_MAF}" > "${FUNCOMAF}"
echo ${MAF_FILES} | xargs -n 1 -I {} grep -hv "^Hugo_Symbol" {} >> "${FUNCOMAF}"
# cleanup
echo ${MAF_FILES} | xargs rm

# targeted re alignment
echo "[STAGE 4/5] Performing targeted re-alignment with HISAT2..."

echo "--> Creating BED file from MAF..."
awk 'BEGIN {OFS="\t"}
  # Condition: Skip header, comments, AND lines where columns 6 or 7 are not numbers
  NR > 1 && !/^#/ && $6 ~ /^[0-9]+$/ && $7 ~ /^[0-9]+$/ {
    chrom = $5;
     
    # Logic to correctly handle all contig name variations in your reference
    if (chrom == "MT") {
      chrom = "chrM";
    } else if (chrom !~ /^KN/ && chrom !~ /^JTFH/) {
      if (chrom !~ /^chr/) {
        chrom = "chr" chrom;
      }
    }
     
    start = $6 - 1; # Convert to 0-based start
    end = $7;
    print chrom, start, end;
  }' ${FUNCOMAF} > ${REALIGN_DIR}/variants.bed

gatk BedToIntervalList -I ${REALIGN_DIR}/variants.bed \
  -O ${REALIGN_DIR}/variants.interval_list -SD ${GN}

# process RNA sample
echo "--> Re-aligning RNA reads..."
samtools view -@ ${NUM_PARALLEL_JOBS} \
  -L ${REALIGN_DIR}/variants.bed ${RNA_BAM_PATH} | \
  cut -f1 | sort -u > ${REALIGN_DIR}/${BASENAME}_read_names.txt

gatk FilterSamReads -I ${RNA_BAM_PATH} \
  -O ${REALIGN_DIR}/${BASENAME}.filtered.bam \
  --READ_LIST_FILE ${REALIGN_DIR}/${BASENAME}_read_names.txt \
  --FILTER includeReadList

gatk SamToFastq -I ${REALIGN_DIR}/${BASENAME}.filtered.bam \
  -F ${REALIGN_DIR}/${BASENAME}_1.fastq.gz -F2 ${REALIGN_DIR}/${BASENAME}_2.fastq.gz

hisat2 -p ${NUM_PARALLEL_JOBS} \
  -x ${HISAT2_INDEX} \
  -1 ${REALIGN_DIR}/${BASENAME}_1.fastq.gz \
  -2 ${REALIGN_DIR}/${BASENAME}_2.fastq.gz \
  --summary-file ${REALIGN_DIR}/${BASENAME}.hisat2.summary.txt \
  --rg-id ${BASENAME} --rg "SM:${BASENAME}" | \
  samtools sort -o ${REALIGN_DIR}/${BASENAME}.realigned.bam
samtools index ${REALIGN_DIR}/${BASENAME}.realigned.bam

# if input DNA bam
if [ -n "${DNA_BAM_N}" ]; then
  echo "--> Re-aligning DNA reads..."
  samtools view -@ ${NUM_PARALLEL_JOBS} \
    -L ${REALIGN_DIR}/variants.bed ${DNA_BAM_N} | \
    cut -f1 | sort -u > ${REALIGN_DIR}/${NORMAL_SAMPLE_NAME}_read_names.txt

  echo "  -> Filtering BAM by read name to include pairs..."
  gatk FilterSamReads \
      -I ${DNA_BAM_N} \
      -O ${REALIGN_DIR}/${NORMAL_SAMPLE_NAME}.filtered.bam \
      --READ_LIST_FILE ${REALIGN_DIR}/${NORMAL_SAMPLE_NAME}_read_names.txt \
      --FILTER includeReadList

  echo "  -> Converting filtered BAM to FASTQ..."
  gatk SamToFastq \
      -I ${REALIGN_DIR}/${NORMAL_SAMPLE_NAME}.filtered.bam \
      -F ${REALIGN_DIR}/${NORMAL_SAMPLE_NAME}_1.fastq.gz \
      -F2 ${REALIGN_DIR}/${NORMAL_SAMPLE_NAME}_2.fastq.gz

  echo "  -> Re-aligning Control reads with HISAT2..."
  hisat2 -p ${NUM_PARALLEL_JOBS} \
      -x ${HISAT2_INDEX} \
      -1 ${REALIGN_DIR}/${NORMAL_SAMPLE_NAME}_1.fastq.gz \
      -2 ${REALIGN_DIR}/${NORMAL_SAMPLE_NAME}_2.fastq.gz \
      --summary-file ${REALIGN_DIR}/${NORMAL_SAMPLE_NAME}.hisat2.summary.txt \
      --rg-id ${NORMAL_SAMPLE_NAME} \
      --rg "SM:${NORMAL_SAMPLE_NAME}" \
      | samtools sort -o ${REALIGN_DIR}/${NORMAL_SAMPLE_NAME}.realigned.bam
  samtools index ${REALIGN_DIR}/${NORMAL_SAMPLE_NAME}.realigned.bam
fi

rm -f ${REALIGN_DIR}/*.filtered.bam ${REALIGN_DIR}/*_read_names.txt

# final variant calling
echo "[STAGE 5/5] Final variant re-calling..."
# The final Mutect2 call must also be conditional
FINAL_VCF_OUT=${REALIGN_DIR}/${BASENAME}.realigned.vcf.gz
FINAL_MUTECT2_CMD_BASE="gatk Mutect2 -R '${GN}' -I '${REALIGN_DIR}/${BASENAME}.realigned.bam' -L '${REALIGN_DIR}/variants.interval_list' --germline-resource '${GERMLINE_RESOURCE}' --panel-of-normals '${PANEL_OF_NORMALS}' -O '${FINAL_VCF_OUT}'"

# The variable MUST be quoted inside the [ ] test to handle the empty case correctly.
if [ -n "${DNA_BAM_N}" ]; then
    FINAL_MUTECT2_CMD="${FINAL_MUTECT2_CMD_BASE} -I '${REALIGN_DIR}/${NORMAL_SAMPLE_NAME}.realigned.bam' -normal '${NORMAL_SAMPLE_NAME}'"
else
    FINAL_MUTECT2_CMD="${FINAL_MUTECT2_CMD_BASE}"
fi
eval $FINAL_MUTECT2_CMD

echo "--> Complete."
date
