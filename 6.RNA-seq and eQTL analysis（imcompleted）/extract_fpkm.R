# 创建R脚本提取FPKM值   
# Ballgown提取FPKM值的R脚本  

# 加载必要的库  
library(ballgown)  

# 设置工作目录和输入路径  
args <- commandArgs(trailingOnly = TRUE)  
if(length(args) < 1) {  
    stop("请提供Ballgown输入目录")  
}  

input_dir <- args[1]  

# 检查输入目录是否存在  
if(!dir.exists(input_dir)) {  
    stop(paste("输入目录不存在:", input_dir))  
}  

# 创建Ballgown对象  
cat("创建Ballgown对象...\n")  

# 假设所有GTF文件都在input_dir目录下  
sample_dirs <- list.dirs(input_dir, recursive = FALSE)  

if(length(sample_dirs) == 0) {  
    # 如果没有子目录，直接从当前目录读取  
    bg <- ballgown(dataDir = input_dir, meas = 'all')  
} else {  
    # 如果有子目录，从子目录读取  
    bg <- ballgown(samples = sample_dirs, meas = 'all')  
}  

cat("Ballgown对象创建完成\n")  

# 提取基因表达量  
cat("提取基因表达量...\n")  

# 获取基因ID和表达量  
gene_expression <- gexpr(bg)  

# 计算每个基因的中位数表达量  
median_expression <- apply(gene_expression, 1, median)  

# 计算每个基因的平均表达量  
mean_expression <- apply(gene_expression, 1, mean)  

# 计算每个基因的FPKM值 (Ballgown已经提供了FPKM)  
fpkm_values <- fpkm(bg)  

# 获取样本信息  
sample_names <- sampleNames(bg)  

# 创建表达量矩阵  
expression_matrix <- data.frame(  
    GeneID = rownames(gene_expression),  
    Median_FPKM = median_expression,  
    Mean_FPKM = mean_expression,  
    gene_expression,  
    stringsAsFactors = FALSE,  
    check.names = FALSE  
)  

# 添加FPKM矩阵 (如果可用)  
if(!is.null(fpkm_values)) {  
    fpkm_matrix <- data.frame(  
        GeneID = rownames(fpkm_values),  
        fpkm_values,  
        stringsAsFactors = FALSE,  
        check.names = FALSE  
    )  
    
    write.table(fpkm_matrix, "gene_fpkm_matrix.txt",   
                row.names = FALSE, quote = FALSE, sep = "\t")  
    
    cat("FPKM矩阵已保存到: gene_fpkm_matrix.txt\n")  
}  

# 保存表达量矩阵  
write.table(expression_matrix, "gene_expression_matrix.txt",   
            row.names = FALSE, quote = FALSE, sep = "\t")  

cat("表达量矩阵已保存到: gene_expression_matrix.txt\n")  

# 统计信息  
cat("\n表达量统计:\n")  
cat("基因总数:", nrow(expression_matrix), "\n")  
cat("样本数量:", ncol(gene_expression), "\n")  
cat("中位数FPKM > 1的基因数:", sum(median_expression > 1), "\n")  

# 创建统计摘要文件  
sink("expression_summary.txt")  
cat("基因表达量提取摘要\n")  
cat("====================\n\n")  

cat("分析时间:", date(), "\n")  
cat("样本数量:", length(sample_names), "\n")  
cat("基因总数:", nrow(expression_matrix), "\n\n")  

cat("样本列表:\n")  
for(sample in sample_names) {  
    cat("- ", sample, "\n")  
}  

cat("\n表达量统计:\n")  
cat("中位数FPKM > 1的基因数:", sum(median_expression > 1), "\n")  
cat("平均FPKM > 1的基因数:", sum(mean_expression > 1), "\n")  

sink()  

cat("表达量提取完成\n")  
   