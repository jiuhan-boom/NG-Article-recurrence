library(ggplot2)  
library(dplyr)  

# 读取XP-CLR结果  
if(file.exists("merged_candidate_sweeps_xpclr.txt")) {  
    xpclr_results <- read.table("merged_candidate_sweeps_xpclr.txt", header=TRUE)  
    xpclr_results$Method <- "XP-CLR"  
    
    cat("XP-CLR候选区域数:", nrow(xpclr_results), "\n")  
}  

# 读取XP-EHH结果  
if(file.exists("merged_significant_xpehh_regions.txt")) {  
    xpehh_results <- read.table("merged_significant_xpehh_regions.txt", header=TRUE)  
    colnames(xpehh_results)[1:3] <- c("Chr", "Start", "End")  
    xpehh_results$Method <- "XP-EHH"  
    
    cat("XP-EHH显著区域数:", nrow(xpehh_results), "\n")  
}  

# 合并两种方法的结果  
if(exists("xpclr_results") && exists("xpehh_results")) {  
    combined_results <- rbind(xpclr_results[,c("Chr", "Start", "End", "Method")],   
                             xpehh_results[,c("Chr", "Start", "End", "Method")])  
    
    # 按染色体统计每种方法的区域数  
    method_summary <- combined_results %>%  
        group_by(Chr, Method) %>%  
        summarise(RegionCount = n(), .groups = 'drop') %>%  
        spread(Method, RegionCount, fill = 0)  
    
    print(method_summary)  
    
    # 寻找两种方法共同检测到的区域 (重叠区域)  
    find_overlaps <- function(xpclr_df, xpehh_df) {  
        overlaps <- data.frame()  
        
        for(chr in unique(c(xpclr_df$Chr, xpehh_df$Chr))) {  
            xpclr_chr <- xpclr_df[xpclr_df$Chr == chr, ]  
            xpehh_chr <- xpehh_df[xpehh_df$Chr == chr, ]  
            
            for(i in 1:nrow(xpclr_chr)) {  
                for(j in 1:nrow(xpehh_chr)) {  
                    # 检查区间是否重叠  
                    if(xpclr_chr$Start[i] <= xpehh_chr$End[j] &&   
                       xpclr_chr$End[i] >= xpehh_chr$Start[j]) {  
                        
                        overlap_start <- max(xpclr_chr$Start[i], xpehh_chr$Start[j])  
                        overlap_end <- min(xpclr_chr$End[i], xpehh_chr$End[j])  
                        
                        overlap_row <- data.frame(  
                            Chr = chr,  
                            Start = overlap_start,  
                            End = overlap_end,  
                            XPCLR_Method = TRUE,  
                            XPEHH_Method = TRUE,  
                            stringsAsFactors = FALSE  
                        )  
                        
                        overlaps <- rbind(overlaps, overlap_row)  
                    }  
                }  
            }  
        }  
        
        return(overlaps)  
    }  
    
    common_regions <- find_overlaps(xpclr_results, xpehh_results)  
    
    cat("两种方法共同检测到的区域数:", nrow(common_regions), "\n")  
    
    # 保存整合结果  
    write.table(combined_results, "combined_selective_sweeps.txt",   
                row.names=FALSE, quote=FALSE, sep="\t")  
    
    if(nrow(common_regions) > 0) {  
        write.table(common_regions, "common_selective_sweeps.txt",   
                    row.names=FALSE, quote=FALSE, sep="\t")  
    }  
    
}  

# 创建结果汇总报告  
sink("selective_sweep_summary.txt")  
cat("全基因组选择清除扫描结果汇总\n")  
cat("================================\n\n")  

cat("分析时间:", date(), "\n\n")  

if(exists("xpclr_results")) {  
    cat("XP-CLR分析:\n")  
    cat("- 候选区域数:", nrow(xpclr_results), "\n")  
}  

if(exists("xpehh_results")) {  
    cat("XP-EHH分析:\n")  
    cat("- 显著区域数:", nrow(xpehh_results), "\n")  
}  

if(exists("common_regions")) {  
    cat("共同检测区域数:", nrow(common_regions), "\n")  
}  

cat("\n输出文件:\n")  
cat("- 整合结果: combined_selective_sweeps.txt\n")  

if(exists("common_regions") && nrow(common_regions) > 0) {  
    cat("- 共同区域: common_selective_sweeps.txt\n")  
}  

sink()  

cat("结果整合完成\n")  
    