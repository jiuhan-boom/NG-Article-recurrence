# 创建R脚本进行eQTL热点识别和置换检验
# eQTL热点识别和置换检验R脚本

library(dplyr)

# 读取eQTL结果数据 (假设包含染色体、位置、P值等信息)
if(file.exists("eqtl_results.txt")) {
    eqtl_data <- read.table("eqtl_results.txt", header = TRUE, sep = "\t")
} else {
    # 创建示例eQTL数据用于演示
    set.seed(456)
    n_eqtls <- 10000
    
    eqtl_data <- data.frame(
        SNP = paste0("SNP", 1:n_eqtls),
        Gene = paste0("Gene", sample(1:5000, n_eqtls, replace = TRUE)),
        Chr = sample(1:10, n_eqtls, replace = TRUE),
        Pos = sample(1:100000000, n_eqtls),
        Pvalue = runif(n_eqtls),
        Beta = rnorm(n_eqtls, 0, 0.5),
        stringsAsFactors = FALSE
    )
    
    write.table(eqtl_data, "eqtl_results.txt", row.names = FALSE, quote = FALSE, sep = "\t")
}

cat("eQTL数据维度:", dim(eqtl_data), "\n")

# eQTL热点识别 (按染色体区域统计eQTL密度)
cat("识别eQTL热点区域...\n")

# 定义窗口大小 (例如1Mb窗口)
window_size <- 1000000

# 按染色体分组统计eQTL密度
eqtl_density <- eqtl_data %>%
    mutate(Window_Start = (Pos %/% window_size) * window_size) %>%
    group_by(Chr, Window_Start) %>%
    summarise(
        Window_End = Window_Start + window_size,
        EQTL_Count = n(),
        Mean_Pvalue = mean(Pvalue),
        .groups = 'drop'
    ) %>%
    arrange(Chr, Window_Start)

# 识别显著热点 (top 5%密度区域)
threshold_count <- quantile(eqtl_density$EQTL_Count, 0.95)
hotspot_regions <- eqtl_density[eqtl_density$EQTL_Count >= threshold_count, ]

cat("识别到", nrow(hotspot_regions), "个eQTL热点区域\n")

# 置换检验评估观察数据分布偏离均匀分布的统计显著性
cat("进行置换检验...\n")

perform_permutation_test <- function(observed_data, n_permutations = 1000) {
    observed_statistic <- length(hotspot_regions$EQTL_Count)
    
    permutation_statistics <- numeric(n_permutations)
    
    for(i in 1:n_permutations) {
        # 随机打乱染色体位置分配
        permuted_data <- observed_data
        permuted_data$Chr <- sample(permuted_data$Chr)
        
        # 计算置换后的统计量 (热点区域数)
        permuted_density <- permuted_data %>%
            mutate(Window_Start = (Pos %/% window_size) * window_size) %>%
            group_by(Chr, Window_Start) %>%
            summarise(EQTL_Count = n(), .groups = 'drop') %>%
            arrange(Chr, Window_Start)
        
        threshold_perm <- quantile(permuted_density$EQTL_Count, 0.95)
        permuted_hotspots <- permuted_density[permuted_density$EQTL_Count >= threshold_perm, ]
        
        permutation_statistics[i] <- nrow(permuted_hotspots)
        
        if(i %% 100 == 0) {
            cat("完成置换", i, "/", n_permutations, "\n")
        }
    }
    
    # 计算p值 (观察值大于等于置换值的比例)
    p_value <- sum(permutation_statistics >= observed_statistic) / n_permutations
    
    return(list(
        observed_statistic = observed_statistic,
        permutation_statistics = permutation_statistics,
        p_value = p_value
    ))
}

# 执行置换检验 (使用较少的置换次数以节省时间)
permutation_result <- perform_permutation_test(eqtl_data, n_permutations = 100)

# 保存结果
write.table(eqtl_density, "eqtl_density_by_region.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(hotspot_regions, "eqtl_hotspot_regions.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# 创建可视化图
if(require(ggplot2)) {
    # eQTL密度图
    density_plot <- ggplot(eqtl_density, aes(x = Window_Start/1000000, y = EQTL_Count)) +
        geom_point(alpha = 0.6) +
        facet_wrap(~Chr, scales = "free_x") +
        labs(title = "eQTL密度分布", x = "染色体位置 (Mb)", y = "eQTL数量") +
        theme_minimal()
    
    ggsave("eqtl_density_plot.png", width = 12, height = 8)
    
    # 置换检验结果图
    if(exists("permutation_result")) {
        hist_plot <- ggplot(data.frame(statistics = permutation_result$permutation_statistics), 
                           aes(x = statistics)) +
            geom_histogram(bins = 30, fill = "lightblue", alpha = 0.7) +
            geom_vline(xintercept = permutation_result$observed_statistic, 
                      color = "red", linetype = "dashed", size = 1) +
            labs(title = "置换检验结果", 
                 x = "热点区域数", 
                 y = "频率",
                 subtitle = paste("观察值:", permutation_result$observed_statistic,
                                ", P值:", round(permutation_result$p_value, 4))) +
            theme_minimal()
        
        ggsave("permutation_test_plot.png", width = 8, height = 6)
    }
}

# 创建分析报告
sink("eqtl_hotspot_summary.txt")
cat("eQTL热点识别和置换检验摘要\n")
cat("============================\n\n")

cat("分析时间:", date(), "\n")
cat("eQTL总数:", nrow(eqtl_data), "\n")
cat("窗口大小:", window_size/1000000, "Mb\n\n")

cat("热点识别结果:\n")
cat("- 热点区域数:", nrow(hotspot_regions), "\n")
cat("- 密度阈值:", threshold_count, "\n\n")

if(exists("permutation_result")) {
    cat("置换检验结果:\n")
    cat("- 置换次数: 100\n")
    cat("- 观察统计量:", permutation_result$observed_statistic, "\n")
    cat("- P值:", round(permutation_result$p_value, 4), "\n\n")
}

cat("输出文件:\n")
cat("- 区域eQTL密度: eqtl_density_by_region.txt\n")
cat("- 热点区域: eqtl_hotspot_regions.txt\n")
cat("- 密度图: eqtl_density_plot.png\n")
if(exists("permutation_result")) {
    cat("- 置换检验图: permutation_test_plot.png\n")
}

sink()

cat("eQTL热点识别和置换检验完成\n")