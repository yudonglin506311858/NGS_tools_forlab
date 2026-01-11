#计算在不同基因组区域的结合特异性


#!/bin/bash
set -e

### 参数区（需修改）
sgRNA_seq="TGGCACTGGCTTAGGAG"
chr="chr11"
start=5227214
end=5227234
up_down=250


#hbb
chr="chr11"
start=5227056
end=5227596
up_down=100

#ifng
chr="chr12"
start=68159591
end=68160195
up_down=0


#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
genome_size="/data8/yudonglin/HIV_project/chip-seq/PAMless所有/hg38.chrom.sizes"
#bam_list=("WT-NT-Chip-2_sorted_rmdup.bam" "WT-NT-chip-3_sorted_rmdup.bam"
#          "PmY-NT-Chip-2_sorted_rmdup.bam" "PmY-NT-chip-3_sorted_rmdup.bam"
#          "WT-sg30-32-chip_sorted_rmdup.bam" "WT-sg30-33_sorted_rmdup.bam"
#          "PmY-sg30-32-chip_sorted_rmdup.bam" "PmY-sg30-33_sorted_rmdup.bam")

bam_dir="/data8/yudonglin/PAMless_project/汇总/bowtie2/bam"
bam_list=($(find "$bam_dir" -type f -name "*.bam" -printf "%f\n"))
          
### 输出目录
outdir="./results-hbb"
#outdir="./results-ifng"
mkdir -p $outdir

### Step1: 定义目标区域
target_start=$((start - up_down))
target_end=$((end + up_down))
echo -e "${chr}\t${target_start}\t${target_end}" > target.bed

### Step2: 非目标区域
bedtools complement -i target.bed -g $genome_size > regions_without_seed_pam.bed

### Step3: 对每个BAM执行循环计算
for bam in "${bam_list[@]}"; do
  sample=${bam%%_sorted_rmdup.bam}

  echo "==== Processing $sample ===="

  ## (1) 总reads数
  samtools idxstats $bam | awk '{s+=$3} END {print s}' > ${outdir}/${sample}_total_reads.txt
  total_reads=$(cat ${outdir}/${sample}_total_reads.txt)

  ## (2) 含seed/PAM区域
  bedtools coverage -a target.bed -b $bam -counts > ${outdir}/${sample}_with_seed_pam_counts.bed

  ## (3) 不含seed/PAM区域
  bedtools coverage -a regions_without_seed_pam.bed -b $bam -counts > ${outdir}/${sample}_without_seed_pam_counts.bed

  ## (4) 合并非目标区域
  total_length=$(awk '{sum += $3-$2} END {print sum}' ${outdir}/${sample}_without_seed_pam_counts.bed)
  total_reads_nt=$(awk '{sum += $4} END {print sum}' ${outdir}/${sample}_without_seed_pam_counts.bed)
  echo -e "non-target\t0\t${total_length}\t${total_reads_nt}" > ${outdir}/${sample}_non-target.bed

  ## (5) 计算 RPKM
  awk -v total=$total_reads 'BEGIN{OFS="\t"} {
    len = ($3 - $2) / 1000.0;
    rpkm = ($4 + 0.0) / len / (total / 1e6);
    print $1, $2, $3, $4, rpkm
  }' ${outdir}/${sample}_with_seed_pam_counts.bed > ${outdir}/${sample}_with_seed_pam_rpkm.txt

  awk -v total=$total_reads 'BEGIN{OFS="\t"} {
    len = ($3 - $2) / 1000.0;
    rpkm = ($4 + 0.0) / len / (total / 1e6);
    print $1, $2, $3, $4, rpkm
  }' ${outdir}/${sample}_non-target.bed > ${outdir}/${sample}_without_seed_pam_rpkm.txt

  ## (6) 计算 CPM
  awk -v total=$total_reads 'BEGIN{OFS="\t"} {
    cpm = ($4 + 0.0) / total * 1e6;
    print $1, $2, $3, $4, cpm
  }' ${outdir}/${sample}_with_seed_pam_counts.bed > ${outdir}/${sample}_with_seed_pam_cpm.txt

  awk -v total=$total_reads 'BEGIN{OFS="\t"} {
    cpm = ($4 + 0.0) / total * 1e6;
    print $1, $2, $3, $4, cpm
  }' ${outdir}/${sample}_non-target.bed > ${outdir}/${sample}_without_seed_pam_cpm.txt

done

echo " All samples processed. Outputs in $outdir/"



