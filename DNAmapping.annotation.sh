#!/bin/bash:
###### 1-(1) GEO 下载数据################################
prefetch --max-size 0 --option-file GEOdownload.txt
#--max-size 30G  #默认20G，改成30G，可以下载大体积数据； 0表示取消大小限制

###### 1-(2)转换成fastq格式#################################
OUTPUT_DIR="/Users/lixia/Data/data/DNAseq/fastq"
mkdir -p $OUTPUT_DIR
sample="SRR33264153 SRR33264154"
for s in $sample
do
fasterq-dump -e 8 --split-files -O $OUTPUT_DIR $s
pigz -p 8 -k "${OUTPUT_DIR}/${s}"*.fastq
done

##### 1-(3)reads 质控QC###########
Input_dir="/Users/lixia/Data/data/DNAseq/fastq"
Out_fastQC_dir="/Users/lixia/Data/data/DNAseq/fastq/fastQC"
mkdir -p "$Out_fastQC_dir"
find "$Input_dir" -name "*.fastq" | \
parallel -j 8 "fastqc {} -o $Out_fastQC_dir"

##### 2- remove adaptor,fastp###############
Input_fastq_dir="/Users/lixia/Data/data/DNAseq/fastq"
Out_trimmed_dir="/Users/lixia/Data/data/DNAseq/fastq/trimmed"
mkdir -p "$Out_trimmed_dir"
TN5_ADAPTER="CTGTCTCTTATACACATCT"

sample="R"  #SRRxxxxx数据太大，换成了小数据量的R样品
for s in $sample
do

Reads1="${Input_fastq_dir}/${s}.1.fastq.gz"
Reads2="${Input_fastq_dir}/${s}.2.fastq.gz"

fastp \
        -i "$Reads1" -I "$Reads2" \
        -o "${Out_trimmed_dir}/${s}_1.trimmed.fastq.gz" \
        -O "${Out_trimmed_dir}/${s}_2.trimmed.fastq.gz" \
        --adapter_sequence "$TN5_ADAPTER" \
        --adapter_sequence_r2 "$TN5_ADAPTER" \
        --length_required 30 \
        --thread 8 \
        --html "${Out_trimmed_dir}/${s}.fastp.html" \
        --json "${Out_trimmed_dir}/${s}.fastp.json"
done
#参数注释
#        --adapter_sequence # Read 1 的3'端查找并去除该序列
#        --adapter_sequence_r2  #在 Read 2 的3'端查找并去除该序列
#        --trim_poly_g #从3'端修剪掉连续的 G (默认 10)
#        --length_required 30 #保留读段的最小长度
#        --thread 8 #CPU 线程数为 8
#        --html  #HTML 格式质控报告
#        --json #SON 格式质控报告
#        --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \ #3端滑动,去除低质量碱基
#        --correction # 启用碱基校正
#        --trim_poly_x #修剪PolyX（A/T/C）尾巴--trim_poly_G  --trim_poly_A
#        --low_complexity_filter #低复杂度序列过滤，默认0.3
#        --detect_adapter_for_pe # 让 fastp 自动检测双端数据的接头
#        -Q # 质量值阈值 (默认 Q15)
#        -u # 允许-Q的最大百分比 (默认 40%)
#        -n # 允许的 N 碱基最大数目 (默认 5)
#        -e # 整条 read 的平均质量值阈值
#        --length_limit # 设定 reads 允许的最大长度
#        -D / --dedup：#全局去重
#        --correction: 启用碱基校正。双端数据，重叠区域的碱基不匹配，fastp会根据质量分数进行校正
#        -5 -3 -W 4 \  #使用窗口大小为 4 的滑动窗口进行 5' 和 3' 端质量修剪 (默认阈值 Q20)
#        -n 5 \ # 丢弃包含超过 5 个 N 碱基的 reads

######### 3- BWA-MEM Mapping, Samtools 排序和索引##########
# 生成测试数据
#zmore x.trimmed.fastq.gz | head -n 40000 > BWAtest.fastq
Input_fastq_dir="/Users/lixia/Data/data/DNAseq/fastq"
Out_trimmed_dir="/Users/lixia/Data/data/DNAseq/fastq/trimmed"
Ref_genome="/Users/lixia/Data/database/ref_genome/GRCh38.GENCODE/BWA_index/GRCh38.primary_assembly.genome.fa"

sample="R"
for s in $sample
do

bwa mem -t 8 \ #线程数
  -R '@RG\tID:'"$s"'\tSM:'"$s"'\tPL:ILLUMINA' \ #将完整的读取组信息添加到 SAM 头中
  -M \ #将较短的比对标记为辅助比对。这对于某些下游工具（如 GATK）是推荐的
  "$Ref_genome" \
  "${Out_trimmed_dir}/${s}_1.trimmed.fastqgz" \
  "${Out_trimmed_dir}/${s}_2.trimmed.fastq.gz" \
  | samtools view -@ 8 -Sb - \ #将SAM文件变成bam文件
  | samtools sort -@ 8 -o "${s}_sorted.bam" - \ #排序
  && samtools index "${s}_sorted.bam" #索引
done
######### 4-(1) GATK 标记重复序列#########
sample="R"
for s in $sample
do

gatk MarkDuplicates \
    -I "${s}_sorted.bam" \
    -O "${s}_marked_duplicates.bam" \
    -M "${s}_marked_dup_metrics.txt"
done
####### 4-(2) GATK 碱基质量分数重校准 (BaseRecalibrator)#########
Ref_GATK="/Users/lixia/Data/database/GATK/GRCh38.primary_assembly.genome.fa"
GATK_SNP="/Users/lixia/Data/database/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz"
GATK_Indel="/Users/lixia/Data/database/GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
Chromsome="/Users/lixia/Data/database/GATK/chroms.list"

sample="R"
for s in $sample
do

gatk BaseRecalibrator \
  -I "${s}_marked_duplicates.bam" \
  -R "$Ref_GATK" \
  --known-sites "$GATK_SNP" \
  --known-sites "$GATK_Indel" \
  -O "${s}_recal_data.table"

####### 4-(3) GATK 应用碱基校准 (ApplyBQSR)##########
gatk ApplyBQSR \
  -I "${s}_marked_duplicates.bam" \
  -R "$Ref_GATK" \
  --bqsr-recal-file "${s}_recal_data.table" \
  -O "${s}_recalibrated.bam" \
  --create-output-bam-index true
 # --compression-level 5

###### 4-(4) GATK 变异检测 (HaplotypeCaller),查找单倍体基因型#########
#此步骤很耗内存，花时间，可以限制内存gatk --java-options "-Xmx30g"，或分染色体
###### 4-(4-1)不分染色体如下############
#gatk --java-options "-Xmx30g" #限制内存
gatk HaplotypeCaller \
  -R "$Ref_GATK" \
  -I "${s}_recalibrated.bam" \
  -O "${s}_raw_variants.vcf" \
  -ERC GVCF \
  --native-pair-hmm-threads 8

###### 4-(4-2-1)分染色体， 为每个染色体运行GATK HaplotypeCaller###########
#"$Chromsome" 只包含了主染色体，没有存contig
cat "$Chromsome" | parallel -j $(sysctl -n hw.ncpu) \ #paralle根据系统自动分派线程
  "gatk HaplotypeCaller \
  -I "${s}_recalibrated.bam" \
  -R "$Ref_GATK" \
  -L {} \ #占位符
  --dbsnp ${GATK_SNP} \
  -O "${s}_{}.g.vcf.gz" \
  -ERC GVCF" #指定输出模式为 GVCF,这个是GATK分析的格式

###### 4-（4-2-2）合并染色体，合并所有GVCF文件##########
# 首先创建一个列表文件，包含所有染色体生成的GVCF文件
find "$(pwd)" -name "${s}_*.g.vcf.gz" | sort > "${s}.gvcfs.list"

    gatk CombineGVCFs \
        -R "${Ref_GATK}" \
        $(while read vcf; do echo "-V $vcf"; done < "${s}.gvcfs.list") \
        -O "${s}.combined.g.vcf.gz"

###### 4-（4-2-2）合并染色体的常用方案，但是此命令不work.
# -V后面的语法有错，下次尝试一次输入每一个"${s}_*.g.vcf.gz"试下
       #gatk GenomicsDBImport \
       # -V "vcfs://$(pwd)/${s}.gvcfs.list" \
       # --genomicsdb-workspace-path "${s}.genomics_db" \
       # --tmp-dir ./tmp

###### 4-(4-3) 联合基因分型，整合多倍体基因型做变硬calling####
gatk GenotypeGVCFs \
  -R "${Ref_GATK}" \
  -V "${s}.combined.g.vcf.gz" \
  --dbsnp ${GATK_SNP} \
  -O "${s}.final.vcf.gz"

##### 4-（4-4）变异过滤，硬过滤（Hard Filtering），适用于单个样本分析######
# 提取SNP
gatk SelectVariants \
   -R "${Ref_GATK}" \
   -V "${s}.final.vcf.gz" \
    -select-type SNP \
   -O "${s}.raw_snps.vcf.gz"

 # 过滤SNP
 gatk VariantFiltration \
        -R "${Ref_GATK}" \
        -V "${s}.raw_snps.vcf.gz" \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        --filter-expression "SOR > 3.0" --filter-name "SOR3" \
        --filter-expression "FS > 60.0" --filter-name "FS60" \
        --filter-expression "MQ < 40.0" --filter-name "MQ40" \
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O "${s}.filtered_snps.vcf.gz"

# 提取INDEL
    gatk SelectVariants \
        -R "${Ref_GATK}" \
        -V "${s}.final.vcf.gz" \
        -select-type INDEL \
        -O "${s}.raw_indels.vcf.gz"

# 过滤INDEL
    gatk VariantFiltration \
        -R "${Ref_GATK}" \
        -V "${s}.raw_indels.vcf.gz" \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        --filter-expression "FS > 200.0" --filter-name "FS200" \
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O "${s}.filtered_indels.vcf.gz"

##### 4-（5） 合并过滤后的变异（如果使用硬过滤）
   gatk MergeVcfs \
        -I "${s}.filtered_snps.vcf.gz" \
        -I "${s}.filtered_indels.vcf.gz" \
        -O "${s}.filtered.vcf.gz"

##### 4-（6）移除被过滤的变异#######
    gatk SelectVariants \
        -R "${Ref_GATK}" \
        -V "${s}.filtered.vcf.gz" \
        --exclude-filtered \
        -O "${s}.results.vcf.gz"
##### 生成统计信息，这步不work,所以不用，vcf结果和Ref_GATK文件头不一致，因为vcf结果缺少contig染色体

##### 4-（7）提取突变信息##########
# AD{1} 表示第二位即突变的count数量
bcftools query -f '%CHROM\t%POS\t%REF->%ALT\t[%AD{1}]\t[%DP]\t%QUAL\n' R.results.vcf.gz > mutation_info.txt
