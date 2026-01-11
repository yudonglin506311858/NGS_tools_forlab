setwd("/data8/yudonglin/PAMless_project/汇总/bowtie2/bam/results-ifng")
setwd("/data8/yudonglin/PAMless_project/汇总/bowtie2/bam/results-hbb")

library(dplyr)
library(stringr)
library(readr)

# 自动获取所有RPKM文件
files_with <- list.files(pattern = "_with_seed_pam_rpkm.txt$")
files_without <- list.files(pattern = "_without_seed_pam_rpkm.txt$")

# 定义函数：读取文件并标注信息
read_rpkm <- function(file, condition_label) {
  # 提取样本名，例如 WT_sg30_1
  sample_name <- str_replace(file, "_with_seed_pam_rpkm.txt|_without_seed_pam_rpkm.txt", "")
  
  # 读取数据
  dat <- read.table(file, sep="\t", header=FALSE)
  
  # 确定group（WT或PmY）
  group_prefix <- ifelse(grepl("^PmCas12m", sample_name), "PmCas12m", "PmY")
  
  # 返回统一格式
  data.frame(
    sample = sample_name,
    condition = condition_label,
    group = paste0(group_prefix, "_", condition_label),
    FPKM = dat$V5
  )
}

# 读取两类数据
data_with <- lapply(files_with, read_rpkm, condition_label="Seed&PAM") %>% bind_rows()
data_without <- lapply(files_without, read_rpkm, condition_label="-Seed&PAM") %>% bind_rows()

# 合并
data_all <- bind_rows(data_with, data_without)

# 输出到文件
write.csv(data_all, "fpkm_ifng-20260102.csv")

# 查看汇总结果
print(data_all)


