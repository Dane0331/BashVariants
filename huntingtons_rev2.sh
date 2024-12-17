#!/bin/bash

# Display a friendly introduction message
echo "Welcome to the Genomic Analysis Pipeline Script!"

# Prompt for user inputs
echo "Choose the source for sequencing data:"
echo "1. SRA Accession"
echo "2. NCBI Nucleotide/Gene Accession"
read -p "Enter your choice (1 or 2): " SEQ_SOURCE
if [[ $SEQ_SOURCE -eq 1 ]]; then
    read -p "Enter the SRA Accession Number (e.g., SRR1972739): " SEQ_ACCESSION
elif [[ $SEQ_SOURCE -eq 2 ]]; then
    read -p "Enter the NCBI Gene/Nucleotide Accession Numbers (separated by spaces, e.g., NM_000001 NM_000002): " ACCESSION_LIST
else
    echo "Invalid choice. Exiting."
    exit 1
fi

read -p "Enter the Reference Genome Accession Number (e.g., NC_000001): " REF_GENOME_ID
read -p "Enter the target Chromosome and Gene Region, if any, for the Reference Sequence (e.g., 4:3074000-3240000): " GENE_REGION
read -p "Enter the maximum number of reads to process (e.g., 10000): " MAX_READS
read -p "Enter the project directory path (e.g., /path/to/project): " PROJECT_DIR

# Create the project directory structure
mkdir -p "$PROJECT_DIR/data" "$PROJECT_DIR/results" "$PROJECT_DIR/logs" "$PROJECT_DIR/temp" "$PROJECT_DIR/clinvar" "$PROJECT_DIR/variants"

# Define directory paths
DATA_DIR="$PROJECT_DIR/data"
RESULTS_DIR="$PROJECT_DIR/results"
LOGS_DIR="$PROJECT_DIR/logs"
TEMP_DIR="$PROJECT_DIR/temp"
ANNOTATION_DIR="$PROJECT_DIR/clinvar"
VARIANT_DIR="$PROJECT_DIR/variants"

# Step 1: Download Reference Genome in FASTA format
echo "Downloading reference genome..."
efetch -db nucleotide -id $REF_GENOME_ID -format fasta > "$DATA_DIR/reference.fasta" 2>>"$LOGS_DIR/error.log"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to download reference genome. Check the logs for details."
    exit 1
fi

# Step 2: Download Sequence Data
if [[ $SEQ_SOURCE -eq 1 ]]; then
    # SRA workflow
    SRA_SUBFOLDER="$TEMP_DIR/$SEQ_ACCESSION"
    SRA_FILE="$SRA_SUBFOLDER/$SEQ_ACCESSION.sra"

    if [[ -f "$SRA_FILE" ]]; then
        echo "SRA file for accession '$SEQ_ACCESSION' already exists. Skipping prefetch."
    else
        echo "Prefetching SRA data..."
        prefetch $SEQ_ACCESSION --output-directory "$TEMP_DIR/" 2>>"$LOGS_DIR/error.log"
        if [[ $? -ne 0 ]]; then
            echo "Error: Failed to prefetch SRA data. Check the logs for details."
            exit 1
        fi
    fi

    # Convert SRA to FASTQ
    echo "Converting SRA file to FASTQ format and limiting reads to $MAX_READS..."
    fastq-dump -X $MAX_READS --split-3 "$SRA_FILE" -O "$DATA_DIR/" 2>>"$LOGS_DIR/error.log"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to convert SRA to FASTQ. Check the logs for details."
        exit 1
    fi
elif [[ $SEQ_SOURCE -eq 2 ]]; then
    # NCBI nucleotide/gene accession workflow
    for ACCESSION in $ACCESSION_LIST; do
        echo "Downloading sequence data for $ACCESSION from NCBI..."
        efetch -db nucleotide -id $ACCESSION -format fasta > "$DATA_DIR/${ACCESSION}.fasta" 2>>"$LOGS_DIR/error.log"
        if [[ $? -ne 0 ]]; then
            echo "Error: Failed to download sequence data for $ACCESSION. Check the logs for details."
            exit 1
        fi
        echo "Sequence data for $ACCESSION downloaded successfully."
    done
    
    # Combine all FASTA files into one
    COMBINED_FASTA="$DATA_DIR/combined_sequences.fasta"
    cat "$DATA_DIR"/*.fasta > "$COMBINED_FASTA"
    echo "All sequences combined into $COMBINED_FASTA."
fi

# Step 3: Perform Quality Control
echo "Performing quality control with fastp..."
if [[ -f "$DATA_DIR/${SEQ_ACCESSION}_1.fastq" && -f "$DATA_DIR/${SEQ_ACCESSION}_2.fastq" ]]; then
    # Paired-end data
    fastp \
        -i "$DATA_DIR/${SEQ_ACCESSION}_1.fastq" \
        -I "$DATA_DIR/${SEQ_ACCESSION}_2.fastq" \
        -o "$DATA_DIR/${SEQ_ACCESSION}_1.cleaned.fastq" \
        -O "$DATA_DIR/${SEQ_ACCESSION}_2.cleaned.fastq" \
        --html "$RESULTS_DIR/fastp_report.html" \
        --json "$RESULTS_DIR/fastp_report.json" \
        --thread 4 \
        2>>"$LOGS_DIR/error.log"
    if [[ $? -ne 0 ]]; then
        echo "Error: Quality control failed for paired-end data. Check the logs for details."
        exit 1
    fi
    echo "Paired-end quality control completed successfully."
elif [[ -f "$DATA_DIR/${SEQ_ACCESSION}.fastq" ]]; then
    # Single-end data
    fastp \
        -i "$DATA_DIR/${SEQ_ACCESSION}.fastq" \
        -o "$DATA_DIR/${SEQ_ACCESSION}.cleaned.fastq" \
        --html "$RESULTS_DIR/fastp_report.html" \
        --json "$RESULTS_DIR/fastp_report.json" \
        --thread 4 \
        --qualified_quality_phred 15 \
        --unqualified_percent_limit 50 \
        --n_base_limit 5 \
        --disable_length_filtering \
        --verbose
        2>>"$LOGS_DIR/error.log"
    if [[ $? -ne 0 ]]; then
        echo "Error: Quality control failed for single-end data. Check the logs for details."
        exit 1
    fi
    echo "Single-end quality control completed successfully."
else
    echo "Error: No valid FASTQ files found. Exiting."
    exit 1
fi

# Index the reference genome
bwa index "$DATA_DIR/reference.fasta" 2>>"$LOGS_DIR/error.log"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to index the reference genome. Check the logs for details."
    exit 1
fi
echo "Reference genome indexed successfully."

# Step 4: Prepare the Reference Genome
echo "Preparing the reference genome..."
samtools faidx "$DATA_DIR/reference.fasta" 2>>"$LOGS_DIR/error.log"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to index the reference genome. Check the logs for details."
    exit 1
fi

echo "Creating sequence dictionary..."
picard CreateSequenceDictionary R="$DATA_DIR/reference.fasta" O="$DATA_DIR/reference.dict"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to create sequence dictionary. Check the logs for details."
    exit 1
fi

# Step 5: Align Reads
echo "Aligning reads with BWA..."
if [[ -f "$DATA_DIR/${SEQ_ACCESSION}_1.cleaned.fastq" && -f "$DATA_DIR/${SEQ_ACCESSION}_2.cleaned.fastq" ]]; then
    # Paired-end alignment
    bwa mem "$DATA_DIR/reference.fasta" \
        "$DATA_DIR/${SEQ_ACCESSION}_1.cleaned.fastq" \
        "$DATA_DIR/${SEQ_ACCESSION}_2.cleaned.fastq" \
        > "$RESULTS_DIR/aligned.sam" 2>>"$LOGS_DIR/error.log"
elif [[ -f "$DATA_DIR/${SEQ_ACCESSION}.cleaned.fastq" ]]; then
    # Single-end alignment
    bwa mem "$DATA_DIR/reference.fasta" \
        "$DATA_DIR/${SEQ_ACCESSION}.cleaned.fastq" \
        > "$RESULTS_DIR/aligned.sam" 2>>"$LOGS_DIR/error.log"
else
    echo "Error: No cleaned FASTQ files found for alignment. Exiting."
    exit 1
fi

if [[ $? -ne 0 ]]; then
    echo "Error: Failed to align reads with BWA. Check the logs for details."
    exit 1
fi
echo "Alignment completed successfully."


# Step 6: Convert SAM to BAM and Sort
echo "Converting SAM to BAM and sorting..."
samtools view -bS "$RESULTS_DIR/aligned.sam" | samtools sort -o "$RESULTS_DIR/aligned_sorted.bam" 2>>"$LOGS_DIR/error.log"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to convert and sort BAM file. Check the logs for details."
    exit 1
fi

# Step 6.1: Add Read Groups # Reminder to check and experiment with the parameters here
echo "Adding read groups to BAM file..."
gatk AddOrReplaceReadGroups \
    -I "$RESULTS_DIR/aligned_sorted.bam" \
    -O "$RESULTS_DIR/aligned_sorted_rg.bam" \
    --RGID 1 \
    --RGLB lib1 \
    --RGPL illumina \
    --RGPU unit1 \
    --RGSM sample1 \
    2>>"$LOGS_DIR/error.log"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to add read groups. Check the logs for details."
    exit 1
fi
echo "Read groups added successfully."

# Step 7: Mark Duplicates
echo "Marking duplicates..."
gatk MarkDuplicates \
    -I "$RESULTS_DIR/aligned_sorted_rg.bam" \
    -O "$RESULTS_DIR/dedup.bam" \
    -M "$RESULTS_DIR/marked_dup_metrics.txt" \
    2>>"$LOGS_DIR/error.log"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to mark duplicates. Check the logs for details."
    exit 1
fi

# Index the deduplicated BAM file
echo "Indexing deduplicated BAM file..."
samtools index "$RESULTS_DIR/dedup.bam" 2>>"$LOGS_DIR/error.log"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to index deduplicated BAM file. Check the logs for details."
    exit 1
fi
echo "Deduplicated BAM file indexed successfully."

# Step 8: Variant Calling
echo "Calling variants with GATK HaplotypeCaller..."
gatk HaplotypeCaller \
    -R "$DATA_DIR/reference.fasta" \
    -I "$RESULTS_DIR/dedup.bam" \
    -O "$RESULTS_DIR/raw_variants.vcf" \
    2>>"$LOGS_DIR/error.log"
if [[ $? -ne 0 ]]; then
    echo "Error: Variant calling failed. Check the logs for details."
    exit 1
fi

# Step 9: Filter Variants using GATK
echo "Filtering variants with GATK VariantFiltration..."
gatk VariantFiltration \
    -R "$DATA_DIR/reference.fasta" \
    -V "$RESULTS_DIR/raw_variants.vcf" \
    -O "$RESULTS_DIR/filtered_variants.vcf" \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "LowQual" \
    2>>"$LOGS_DIR/error.log"
if [[ $? -ne 0 ]]; then
    echo "Error: Variant filtering failed. Check the logs for details."
    exit 1
fi
echo "Variant filtering completed successfully. Filtered variants saved to $RESULTS_DIR/filtered_variants.vcf."

# Step 10: Download ClinVar Variants for the Gene
echo "Step 10: Downloading ClinVar variants for $REF_GENOME_ID..."
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
wget -P "$ANNOTATION_DIR" "$CLINVAR_URL"
wget -P "$ANNOTATION_DIR" "${CLINVAR_URL}.tbi"

# Step 11: Filter ClinVar Variants for the Gene Region
echo "Step 11: Filtering ClinVar variants for $REF_GENOME_ID..."
if [[ -z "$GENE_REGION" ]]; then
    echo "Error: Gene region is not specified. Please provide a valid region." >&2
    exit 1
fi

bcftools view -r "$GENE_REGION" "$ANNOTATION_DIR/clinvar.vcf.gz" > "$VARIANT_DIR/${REF_GENOME_ID}_clinvar_variants.vcf"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to filter ClinVar variants. Check input files and region." >&2
    exit 1
fi

# Step 11.1: Check if variants are found
if [[ $(grep -v "^#" "$VARIANT_DIR/${REF_GENOME_ID}_clinvar_variants.vcf" | wc -l) -eq 0 ]]; then
    echo "No variants found in region $GENE_REGION. Exiting." >&2
    exit 1
fi

# Step 12: Functional Annotation with SnpEff 
echo "Step 12: Annotating variants with SnpEff..." # Reminder to try and configure this to use other snpEff databases
snpEff GRCh38.99 \ 
    -v "$VARIANT_DIR/${REF_GENOME_ID}_clinvar_variants.vcf" \
    -stats "$ANNOTATION_DIR/${REF_GENOME_ID}_annotated_variants.html" \
    > "$ANNOTATION_DIR/${REF_GENOME_ID}_annotated_variants.vcf"

bgzip -c "$VARIANT_DIR/${REF_GENOME_ID}_clinvar_variants.vcf" > "$VARIANT_DIR/${REF_GENOME_ID}_clinvar_variants.vcf.gz"
tabix -p vcf "$VARIANT_DIR/${REF_GENOME_ID}_clinvar_variants.vcf.gz"

bgzip -c "$ANNOTATION_DIR/${REF_GENOME_ID}_annotated_variants.vcf" > "$ANNOTATION_DIR/${REF_GENOME_ID}_annotated_variants.vcf.gz"
tabix -p vcf "$ANNOTATION_DIR/${REF_GENOME_ID}_annotated_variants.vcf.gz"


# Display output file summary
echo "Pipeline completed successfully!"
echo "Results Summary:"
echo "- BAM file with Read Groups: $RESULTS_DIR/aligned_sorted_rg.bam"
echo "- Deduplicated BAM file: $RESULTS_DIR/dedup.bam"
echo "- Filtered Variants (VCF): $RESULTS_DIR/filtered_variants.vcf"

# Farewell message
echo "Goodbye! You can visualize the results using IGV. Logs are saved in the '$LOGS_DIR' folder."

