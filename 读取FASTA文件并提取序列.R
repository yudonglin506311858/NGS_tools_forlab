# ======================
# 读取FASTA文件并提取序列
# ======================

# # 方法1：使用Biostrings包（专业生物信息学方法）
# if (!require("Biostrings", quietly = TRUE)) {
#   if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install("Biostrings")
# }
library(Biostrings)

# 读取FASTA文件
pcsk9_seq <- readDNAStringSet("/data/yudonglin/reference/virus/PCSK9/PCSK9.fa")  # 请确保文件路径正确

# 提取序列信息
seq_name <- names(pcsk9_seq)[1]  # 获取序列名称
seq_length <- width(pcsk9_seq)[1]  # 获取序列长度
cat("序列名称:", seq_name, "\n")
cat("序列长度:", seq_length, "\n")

# 提取2955-2975位序列
target_seq <- subseq(pcsk9_seq[[1]], start = 2955, end = 2975)
cat("\n提取的序列(2955-2975):\n")
print(target_seq)
cat("序列长度:", width(target_seq), "个碱基\n")

# 提取2955-2975位序列
target_seq <- subseq(pcsk9_seq[[1]], start = 1641, end = 1664)
cat("\n提取的序列(2955-2975):\n")
print(target_seq)
cat("序列长度:", width(target_seq), "个碱基\n")


#CTTCCTTCCTCCCCCACCTCCCTC
# 创建目标序列文件（用于精确匹配）
echo "CTTCCTTCCTCCCCCACCTCCCTC" > target_seq.txt

# 批量统计R1文件中的出现次数
for r1_file in *.fq.gz; do
sample=${r1_file%.fq.gz}
echo -n "处理 ${sample}... "

# 统计序列出现次数（精确匹配，不分大小写）
count=$(zgrep -o -i "CTTCCTTCCTCCCCCACCTCCCTC" "$r1_file" | wc -l)

# 计算总reads数（FASTQ文件每4行一个read）
total_reads=$(zcat "$r1_file" | wc -l | awk '{print $1/4}')

# 计算频率
if [ $total_reads -gt 0 ]; then
frequency=$(echo "scale=6; $count / $total_reads * 100" | bc)
else
  frequency=0
fi

echo "出现次数: $count, 总reads: $total_reads, 频率: ${frequency}%"
done




