# ADMIXTURE结果可视化R脚本

library(ggplot2)
library(reshape2)
library(RColorBrewer)

# 读取Q矩阵文件 (以K=3为例)
k_value <- 3
q_file <- paste0("admixture_k", k_value, ".Q")
sample_file <- "samples.txt"  # 样本信息文件

# 读取Q矩阵数据
q_data <- read.table(q_file, header=FALSE)
colnames(q_data) <- paste0("Ancestry", 1:k_value)

# 添加样本ID (需要根据实际情况调整)
samples <- read.table(sample_file, header=TRUE)  # 假设有样本信息文件
q_data$Sample <- samples$Sample_ID

# 转换数据格式用于绘图
q_melt <- melt(q_data, id.vars="Sample", variable.name="Ancestry", value.name="Proportion")

# 创建群体结构图
p <- ggplot(q_melt, aes(x=Sample, y=Proportion, fill=Ancestry)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_brewer(palette="Set2") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title=paste("Population Structure (K =", k_value, ")"),
         x="Samples", y="Ancestry Proportion")

print(p)

# 保存图形
ggsave("admixture_structure.png", width=12, height=6)
