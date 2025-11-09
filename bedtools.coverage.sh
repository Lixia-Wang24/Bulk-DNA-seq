#!/bin/bash
# 最基本的coverage计算
Ref_GATK="/Users/lixia/Data/database/GATK/GRCh38.primary_assembly.genome.fa"
#bedtools genomecov -ibam R_sorted.bam  -g "$Ref_GATK" > coverage1.txt

#Input 是bam文件，则不用给genome reference, 即下面的命名即可
bedtools genomecov -ibam R_sorted.bam > coverage2.txt


#好像没有-hist命令
#bedtools genomecov -ibam R_sorted.bam -hist > genome_coverage_hist1.txt


# 或者输出bedgraph格式（更常用）
bedtools genomecov -ibam R_sorted.bam  -bg > coverage1.bedgraph

