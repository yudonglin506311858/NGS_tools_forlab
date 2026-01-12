setwd("/data1/yudonglin/CRISPRscreen-microbe/Legionella_pneumophila")




data<-read_excel("小鼠 CRISPR-Cas9 文库.xlsx",col_names = F)
head(data)
data<-as.data.frame(data)
# 直接使用基础R语法，避免dplyr的rename问题
data_transformed <- data.frame(
  id = data$`...1`,
  gRNA.sequence = data$`...2`,
  gene = sapply(strsplit(data$`...1`, "_"), function(x) x[2]),
  stringsAsFactors = FALSE
)

# 查看结果
head(data_transformed)
dim(data_transformed)  # 应该是 242,961 × 3

write.csv(data_transformed,"library_mouse.csv",quote = F,row.names = F)





data<-read_excel("人类 CRISPR-Cas9 文库.xlsx",col_names = F)
head(data)
data<-as.data.frame(data)
# 直接使用基础R语法，避免dplyr的rename问题
data_transformed <- data.frame(
  id = data$`...1`,
  gRNA.sequence = data$`...2`,
  gene = sapply(strsplit(data$`...1`, "_"), function(x) x[2]),
  stringsAsFactors = FALSE
)

# 查看结果
head(data_transformed)
dim(data_transformed)  # 应该是 242,961 × 3

write.csv(data_transformed,"library_human.csv",quote = F,row.names = F)
