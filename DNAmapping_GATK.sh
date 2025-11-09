#!/bin/bash

GEO Data Download

#prefetch --max-size 0 --option-file GEOdownload.txt
#--max-size 30G # Default is 20G, changed to 30G to download large data; 0 means no size limit

Convert to FASTQ format

:<<!
OUTPUT_DIR="/Users/lixia/Data/data/DNAseq/fastq"
mkdir -p $OUTPUT_DIR
sample="SRR33264153 SRR33264154"
for s in $sample
do
fasterq-dump -e 8 --split-files -O $OUTPUT_DIR $s
pigz -p 8 -k "${OUTPUT_DIR}/${s}"*.fastq

done
!

Reads Quality Control (QC)

:<<!
Input_dir="/Users/lixia/Data/data/DNAseq/fastq"
Out_fastQC_dir="/Users/lixia/Data/data/DNAseq/fastq/fastQC"
mkdir -p "$Out_fastQC_dir"
find "$Input_dir" -name "*.fastq" | 

parallel -j 8 "fastqc {} -o $Out_fastQC_dir"
!

Remove Adaptors with fastp

:<<!
Input_fastq_dir="/Users/lixia/Data/data/DNAseq/fastq"
Out_trimmed_dir="/Users/lixia/Data/data/DNAseq/fastq/trimmed"
mkdir -p "$Out_trimmed_dir"
TN5_ADAPTER="CTGTCTCTTATACACATCT"

sample="R"
for s in $sample
do
Reads1="${Input_fastq_dir}/${s}.1.fastq.gz"
Reads2="${Input_fastq_dir}/${s}.2.fastq.gz"

fastp 

        -i "$Reads1" -I "$Reads2" \
        -o "${Out_trimmed_dir}/${s}_1.trimmed.fastq.gz" 

        -O "${Out_trimmed_dir}/${s}_2.trimmed.fastq.gz" 

        --adapter_sequence "$TN5_ADAPTER" \
        --adapter_sequence_r2 "$TN5_ADAPTER" \
        --length_required 30 \
        --thread 8 \
        --html "${Out_trimmed_dir}/${s}.fastp.html" 

        --json "${Out_trimmed_dir}/${s}.fastp.json"
done
!
######### BWA-MEM Mapping, Samtools Sort and Index ##########

Generate test data

#zmore XXX.trimmed.fastq.gz | head -n 40000 > BWAtest.fastq
:<<!
Input_fastq_dir="/Users/lixia/Data/data/DNAseq/fastq"
Out_trimmed_dir="/Users/lixia/Data/data/DNAseq/fastq/trimmed"
Ref_genome="/Users/lixia/Data/database/ref_genome/GRCh38.GENCODE/BWA_index/GRCh38.primary_assembly.genome.fa"

sample="R"
for s in $sample
do

bwa mem -t 8 

  -R '@RG\tID:'"$s"'\tSM:'"$s"'\tPL:ILLUMINA' \
  -M \
  "$Ref_genome" \
  "${Out_trimmed_dir}/${s}_1.trimmed.fastq.gz" \
  "${Out_trimmed_dir}/${s}_2.trimmed.fastq.gz" \
  | samtools view -@ 8 -Sb - \
  | samtools sort -@ 8 -o "${s}_sorted.bam" - 

  && samtools index "${s}_sorted.bam"

done

######### GATK Mark Duplicates #########
#echo "--- 5.1 GATK Mark Duplicates ---"
sample="R"
for s in $sample
do

gatk MarkDuplicates 

    -I "${s}_sorted.bam" \
    -O "${s}_marked_duplicates.bam" 

    -M "${s}_marked_dup_metrics.txt"
!Ref_GATK="/Users/lixia/Data/database/GATK/GRCh38.primary_assembly.genome.fa"

GATK Base Quality Score Recalibration (BaseRecalibrator)

Ref_GATK="/Users/lixia/Data/database/GATK/GRCh38.primary_assembly.genome.fa"
GATK_SNP="/Users/lixia/Data/database/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz"
GATK_Indel="/Users/lixia/Data/database/GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
Chromsome="/Users/lixia/Data/database/GATK/chroms.list"

:<<!
gatk BaseRecalibrator 

  -I "${s}_marked_duplicates.bam" 

  -R "$Ref_GATK" 

  --known-sites "$GATK_SNP" 

  --known-sites "$GATK_Indel" \
  -O "${s}_recal_data.table"
!

GATK Apply Base Quality Score Recalibration (ApplyBQSR)

:<<!
gatk ApplyBQSR 

  -I "${s}_marked_duplicates.bam" \
  -R "$Ref_GATK" \
  --bqsr-recal-file "${s}_recal_data.table" \
  -O "${s}_recalibrated.bam" 

  --create-output-bam-index true
 # --compression-level 5
!

GATK Variant Calling (HaplotypeCaller)

Memory limit: gatk --java-options "-Xmx30g"

#CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

Process chromosomes in batches: for CHR in {1..22} X Y; do

:<<!
gatk HaplotypeCaller 

  -R "$Ref_GATK" \
  -I "${s}_recalibrated.bam" 

  -O "${s}_raw_variants.vcf" 

  -ERC GVCF 

  --native-pair-hmm-threads 8
!

Run GATK HaplotypeCaller for each chromosome

:<<!
cat "$Chromsome" | parallel -j $(sysctl -n hw.ncpu) \
  "gatk HaplotypeCaller \
  -I "${s}recalibrated.bam" 

  -R "$Ref_GATK" 

  -L {} 

  --dbsnp ${GATK_SNP} \
  -O "${s}{}.g.vcf.gz" 

  -ERC GVCF"
!

After running HaplotypeCaller for all chromosomes, combine all GVCF files

First, create a list file containing all GVCF files generated for each chromosome

:<<!
find "$(pwd)" -name "${s}_*.g.vcf.gz" | sort > "${s}.gvcfs.list"

Run CombineGVCFs

    gatk CombineGVCFs 

        -R "${Ref_GATK}" 

        $(while read vcf; do echo "-V $vcf"; done < "${s}.gvcfs.list") \
        -O "${s}.combined.g.vcf.gz"
!

#rm -rf ${s}.gvcfs.list
:<<!
gatk GenotypeGVCFs \
  -R "${Ref_GATK}" 

  -V "${s}.combined.g.vcf.gz" \
  --dbsnp ${GATK_SNP} \
  -O "${s}.final.vcf.gz"
done
    # Select SNPs
    gatk SelectVariants \
        -R "${Ref_GATK}" 

        -V "${s}.final.vcf.gz" \
        -select-type SNP \
        -O "${s}.raw_snps.vcf.gz"
done
    # Filter SNPs
    gatk VariantFiltration 

        -R "${Ref_GATK}" \
        -V "${s}.raw_snps.vcf.gz" 

        --filter-expression "QD < 2.0" --filter-name "QD2" 

        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" 

        --filter-expression "SOR > 3.0" --filter-name "SOR3" 

        --filter-expression "FS > 60.0" --filter-name "FS60" 

        --filter-expression "MQ < 40.0" --filter-name "MQ40" 

        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" 

        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" 

        -O "${s}.filtered_snps.vcf.gz"
    done
    # Select INDELs
    gatk SelectVariants \
        -R "${Ref_GATK}" 

        -V "${s}.final.vcf.gz" \
        -select-type INDEL \
        -O "${s}.raw_indels.vcf.gz"
    
    # Filter INDELs
    gatk VariantFiltration 

        -R "${Ref_GATK}" \
        -V "${s}.raw_indels.vcf.gz" 

        --filter-expression "QD < 2.0" --filter-name "QD2" 

        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" 

        --filter-expression "FS > 200.0" --filter-name "FS200" 

        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" 

        -O "${s}.filtered_indels.vcf.gz"

    gatk MergeVcfs 

        -I "${s}.filtered_snps.vcf.gz" \
        -I "${s}.filtered_indels.vcf.gz" 

        -O "${s}.filtered.vcf.gz"
done

:<<!    
    # Select variants that passed filtering
    gatk SelectVariants 

        -R "${Ref_GATK}" \
        -V "${s}.filtered.vcf.gz" 

        --exclude-filtered 

        -O "${s}.results.vcf.gz"

done

Extract mutation information

AD{1} represents the second position, which is the count of the alternate (mutant) allele

#bcftools query -f '%CHROM\t%POS\t%REF->%ALT\t[%AD{1}]\t[%DP]\t%QUAL\n' R.results.vcf.gz > mutation_info.txt
!
####### Calculate GATK DepthOfCoverage ####
Ref_GATK="/Users/lixia/Data/database/GATK/GRCh38.primary_assembly.genome.fa"
sample="R"
for s in $sample
do

gatk DepthOfCoverage 

    -R "${Ref_GATK}" \
    -I "${s}_sorted.bam" 

    -O coverage_output 

    --summary-coverage-threshold 1 

    --summary-coverage-threshold 5 

    --summary-coverage-threshold 10 

    --summary-coverage-threshold 20
