# 创建R脚本进行表达量标准化和eQTL分析准备  
# eQTL分析准备的R脚本  

# 加载必要的库  
library(qqman)  

# 读取基因表达量数据  
cat("读取基因表达量数据...\n")  

if(file.exists("gene_expression_matrix.txt")) {  
    expr_data <- read.table("gene_expression_matrix.txt", header = TRUE, sep = "\t", check.names = FALSE)  
    
    # 提取基因ID和表达量列  
    gene_ids <- expr_data$GeneID  
    
    # 提取数值型表达量列 (排除GeneID列)  
    expr_matrix <- expr_data[, !(names(expr_data) %in% c("GeneID", "Median_FPKM", "Mean_FPKM"))]  
    
    cat("表达量数据维度:", dim(expr_matrix), "\n")  
    
} else {  
    stop("未找到基因表达量矩阵文件")  
}  

# 表达量标准化 (正态分位数变换)  
cat("进行正态分位数变换...\n")  

normalize_quantile <- function(x) {  
    # 移除NA值并进行正态分位数变换  
    na_idx <- is.na(x)  
    if(sum(!na_idx) > 0) {  
        x_norm <- x  
        x_norm[!na_idx] <- qqnorm(x[!na_idx], plot.it = FALSE)$x  
        return(x_norm)  
    } else {  
        return(x)  
        }  
}  

# 对每个基因进行正态分位数变换  
expr_normalized <- apply(expr_matrix, 1, normalize_quantile)  

# 转置矩阵使其与原始格式一致  
expr_normalized <- t(expr_normalized)  

# 添加基因ID列  
expr_normalized_df <- data.frame(  
    GeneID = gene_ids,  
    expr_normalized,  
    stringsAsFactors = FALSE,  
    check.names = FALSE  
)  

# 保存标准化后的表达量数据  
write.table(expr_normalized_df, "normalized_expression_matrix.txt",   
            row.names = FALSE, quote = FALSE, sep = "\t")  

cat("标准化表达量数据已保存到: normalized_expression_matrix.txt\n")  

# 筛选高表达基因 (中位数FPKM > 1)  
if("Median_FPKM" %in% names(expr_data)) {  
    high_expr_genes <- expr_data[expr_data$Median_FPKM > 1, ]  
    
    cat("高表达基因数 (FPKM > 1):", nrow(high_expr_genes), "\n")  
    
    # 对高表达基因进行标准化处理  
    if(nrow(high_expr_genes) > 0) {  
        high_expr_matrix <- high_expr_genes[, !(names(high_expr_genes) %in% c("GeneID", "Median_FPKM", "Mean_FPKM"))]  
        
        high_expr_normalized <- apply(high_expr_matrix, 1, normalize_quantile)  
        high_expr_normalized <- t(high_expr_normalized)  
        
        high_expr_normalized_df <- data.frame(  
            GeneID = high_expr_genes$GeneID,  
            high_expr_normalized,  
            stringsAsFactors = FALSE,  
            check.names = FALSE  
        )  
        
        write.table(high_expr_normalized_df, "high_expr_normalized_matrix.txt",   
                    row.names = FALSE, quote = FALSE, sep = "\t")  
        
        cat("高表达基因标准化数据已保存到: high_expr_normalized_matrix.txt\n")  
        
        # 创建eQTL分析准备报告  
        sink("eqtl_preparation_summary.txt")  
        cat("eQTL分析准备摘要\n")  
        cat("=================\n\n")  

        cat("分析时间:", date(), "\n\n")  

        cat("原始基因数:", nrow(expr_data), "\n")  
        cat("高表达基因数 (FPKM > 1):", nrow(high_expr_genes), "\n\n")  

        cat("输出文件:\n")  
        cat("- 标准化表达量矩阵: normalized_expression_matrix.txt\n")  
        cat("- 高表达基因标准化矩阵: high_expr_normalized_matrix.txt\n")  

        sink()  

        cat("eQTL分析准备完成\n")  

} 