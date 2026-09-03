#!/usr/bin/env bash
# =============================================================================
# RNA-seq preprocessing: SRA -> FASTQ -> QC -> trim -> align (GRCh38) -> counts
# Data: NCBI BioProject PRJNA961804
# Tools: SRA Toolkit, FastQC, MultiQC, Trimmomatic, Bowtie2, SAMtools, featureCounts (Subread)
#
# Edit the ACCESSIONS list and the reference/annotation paths, then run:
#   bash preprocessing.sh
# =============================================================================
set -euo pipefail

# ---- configure -------------------------------------------------------------
ACCESSIONS=(SRRxxxxxxx1 SRRxxxxxxx2 SRRxxxxxxx3 SRRxxxxxxx4)   # <-- fill in PRJNA961804 runs
THREADS=8
REF_FASTA="GRCh38.fa"           # GRCh38 genome FASTA
GTF="GRCh38.gtf"                # matching gene annotation (GTF)
INDEX="GRCh38_bt2"              # Bowtie2 index prefix
ADAPTER_TRIM="SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36"   # no ILLUMINACLIP

mkdir -p raw trimmed qc bam counts

# ---- 0. build Bowtie2 index (once) -----------------------------------------
if [ ! -f "${INDEX}.1.bt2" ]; then
  bowtie2-build --threads "$THREADS" "$REF_FASTA" "$INDEX"
fi

for SRR in "${ACCESSIONS[@]}"; do
  # ---- 1. download & convert to FASTQ (paired-end) -------------------------
  prefetch "$SRR"
  fasterq-dump "$SRR" --split-files --threads "$THREADS" -O raw/

  # ---- 2. pre-trim QC ------------------------------------------------------
  fastqc -t "$THREADS" -o qc/ raw/${SRR}_1.fastq raw/${SRR}_2.fastq

  # ---- 3. adapter/quality trimming (Trimmomatic, paired-end) ---------------
  trimmomatic PE -threads "$THREADS" \
    raw/${SRR}_1.fastq raw/${SRR}_2.fastq \
    trimmed/${SRR}_1P.fastq trimmed/${SRR}_1U.fastq \
    trimmed/${SRR}_2P.fastq trimmed/${SRR}_2U.fastq \
    $ADAPTER_TRIM

  # ---- 4. post-trim QC -----------------------------------------------------
  fastqc -t "$THREADS" -o qc/ trimmed/${SRR}_1P.fastq trimmed/${SRR}_2P.fastq

  # ---- 5. alignment to GRCh38 (Bowtie2) ------------------------------------
  bowtie2 -p "$THREADS" -x "$INDEX" \
    -1 trimmed/${SRR}_1P.fastq -2 trimmed/${SRR}_2P.fastq \
    -S bam/${SRR}.sam

  # ---- 6. SAM -> sorted, indexed BAM (SAMtools) ----------------------------
  samtools view -@ "$THREADS" -bS bam/${SRR}.sam | \
    samtools sort -@ "$THREADS" -o bam/${SRR}.sorted.bam
  samtools index bam/${SRR}.sorted.bam
  rm bam/${SRR}.sam
done

# ---- 7. aggregate QC -------------------------------------------------------
multiqc qc/ -o qc/

# ---- 8. gene-level read counting (featureCounts) ---------------------------
featureCounts -T "$THREADS" -p --countReadPairs -a "$GTF" -g gene_name \
  -o counts/counts_raw.txt bam/*.sorted.bam

# ---- 9. tidy to a genes x samples matrix -> cleaned_counts.csv --------------
# Keep the Geneid column and the per-sample count columns (drop columns 2-6:
# Chr, Start, End, Strand, Length), then use as input to deseq2_analysis.R.
tail -n +2 counts/counts_raw.txt | cut -f1,7- > counts/cleaned_counts.tsv
echo "featureCounts done -> counts/cleaned_counts.tsv"
echo "Convert to cleaned_counts.csv (comma-separated; gene IDs as row names) for DESeq2."
