# 准备GCTA格式表型文件的R脚本

# 读取表型数据
if(file.exists("phenotype_data.txt")) {
    pheno_data <- read.table("phenotype_data.txt", header = TRUE, sep = "\t")
} else {
    stop("未找到表型数据文件")
}

# GCTA表型文件格式: FID IID Trait1 Trait2 ...
# 假设SampleID列包含个体ID

# 创建GCTA格式表型文件
gcta_pheno <- data.frame(
    FID = pheno_data$SampleID,
    IID = pheno_data$SampleID,
    pheno_data[, !(names(pheno_data) %in% c("SampleID"))],
    stringsAsFactors = FALSE
)

write.table(gcta_pheno, "gcta_phenotype.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("GCTA格式表型文件已创建: gcta_phenotype.txt\n")
