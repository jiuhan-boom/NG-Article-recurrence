# 创建R脚本进行广义遗传力分析
# 广义遗传力分析R脚本

# 加载必要的库
library(lme4)
library(nlme)

# 读取表型数据
cat("读取表型数据...\n")

# 假设表型数据格式: SampleID, Trait1, Trait2, ..., Genotype, Year
# 这里创建示例数据结构
if(file.exists("phenotype_data.txt")) {
    pheno_data <- read.table("phenotype_data.txt", header = TRUE, sep = "\t")
} else {
    # 创建示例表型数据
    set.seed(123)
    n_samples <- 200
    n_genotypes <- 50
    n_years <- 3
    
    pheno_data <- data.frame(
        SampleID = paste0("Sample", 1:n_samples),
        Genotype = rep(paste0("Geno", 1:n_genotypes), each = n_samples/n_genotypes),
        Year = rep(2020:2022, length.out = n_samples),
        Trait1 = rnorm(n_samples, 50, 10),
        Trait2 = rnorm(n_samples, 30, 5),
        stringsAsFactors = FALSE
    )
    
    write.table(pheno_data, "phenotype_data.txt", row.names = FALSE, quote = FALSE, sep = "\t")
}

cat("表型数据维度:", dim(pheno_data), "\n")

# 计算广义遗传力 H2 = σg² / (σg² + σy²/ny)
# 其中 σg² 是基因型方差分量，σy² 是年份方差分量，ny 是年份数量

# 对每个性状进行分析
traits <- names(pheno_data)[!(names(pheno_data) %in% c("SampleID", "Genotype", "Year"))]
heritability_results <- data.frame()

for(trait in traits) {
    cat("分析性状:", trait, "\n")
    
    # 构建混合效应模型
    # 假设模型: Trait ~ (1|Genotype) + (1|Year)
    formula_str <- paste(trait, "~ (1|Genotype) + (1|Year)")
    
    tryCatch({
        model <- lmer(as.formula(formula_str), data = pheno_data, REML = TRUE)
        
        # 提取方差分量
        varcomp <- VarCorr(model)
        geno_var <- as.numeric(VarCorr(model)$Genotype[1,1])
        year_var <- as.numeric(VarCorr(model)$Year[1,1])
        resid_var <- attr(VarCorr(model), "sc")^2
        
        # 计算年份数量
        n_years <- length(unique(pheno_data$Year))
        
        # 计算广义遗传力
        H2 <- geno_var / (geno_var + year_var/n_years + resid_var)
        
        # 计算标准误 (简化方法)
        H2_se <- sqrt(H2 * (1 - H2) / nrow(pheno_data))
        
        # 保存结果
        result_row <- data.frame(
            Trait = trait,
            Genotype_Variance = geno_var,
            Year_Variance = year_var,
            Residual_Variance = resid_var,
            Number_Years = n_years,
            Heritability_H2 = H2,
            H2_SE = H2_se,
            stringsAsFactors = FALSE
        )
        
        heritability_results <- rbind(heritability_results, result_row)
        
        cat("性状", trait, "广义遗传力:", round(H2, 4), "±", round(H2_se, 4), "\n")
        
    }, error = function(e) {
        cat("性状", trait, "分析失败:", e$message, "\n")
        result_row <- data.frame(
            Trait = trait,
            Genotype_Variance = NA,
            Year_Variance = NA,
            Residual_Variance = NA,
            Number_Years = NA,
            Heritability_H2 = NA,
            H2_SE = NA,
            stringsAsFactors = FALSE
        )
        heritability_results <- rbind(heritability_results, result_row)
    })
}

# 保存结果
write.table(heritability_results, "heritability_results.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# 创建可视化图
if(require(ggplot2)) {
    heritability_plot <- ggplot(heritability_results, aes(x = Trait, y = Heritability_H2)) +
        geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
        geom_errorbar(aes(ymin = Heritability_H2 - H2_SE, ymax = Heritability_H2 + H2_SE), 
                      width = 0.2) +
        labs(title = "广义遗传力估计", x = "性状", y = "广义遗传力 (H²)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave("heritability_plot.png", width = 8, height = 6)
}

# 创建分析报告
sink("heritability_summary.txt")
cat("广义遗传力分析摘要\n")
cat("==================\n\n")

cat("分析时间:", date(), "\n")
cat("样本数量:", nrow(pheno_data), "\n")
cat("基因型数量:", length(unique(pheno_data$Genotype)), "\n")
cat("年份数量:", length(unique(pheno_data$Year)), "\n\n")

cat("性状遗传力估计:\n")
for(i in 1:nrow(heritability_results)) {
    if(!is.na(heritability_results$Heritability_H2[i])) {
        cat(sprintf("%s: %.4f ± %.4f\n", 
            heritability_results$Trait[i],
            heritability_results$Heritability_H2[i],
            heritability_results$H2_SE[i]))
    }
}

cat("\n输出文件:\n")
cat("- 遗传力结果: heritability_results.txt\n")
cat("- 可视化图: heritability_plot.png\n")

sink()

cat("广义遗传力分析完成\n")
    