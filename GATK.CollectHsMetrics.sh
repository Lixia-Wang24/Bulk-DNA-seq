#!/bin/bash
Ref_GATK="/Users/lixia/Data/database/GATK/GRCh38.primary_assembly.genome.fa"
sample="R"
for s in $sample
do

gatk CollectHsMetrics \
  -R"${Ref_GATK}" \
  -I "${s}_sorted.bam" \
  -O "${s}".metrics.txt
done
