#!/bin/bash
#######计算GATK DepthOfCoverage####
Ref_GATK="/Users/lixia/Data/database/GATK/GRCh38.primary_assembly.genome.fa"
sample="R"
for s in $sample
do

gatk DepthOfCoverage \
    -R "${Ref_GATK}" \
    -I "${s}_sorted.bam" \
    -O coverage_output \
    --summary-coverage-threshold 1 \
    --summary-coverage-threshold 5 \
    --summary-coverage-threshold 10 \
    --summary-coverage-threshold 20
done
