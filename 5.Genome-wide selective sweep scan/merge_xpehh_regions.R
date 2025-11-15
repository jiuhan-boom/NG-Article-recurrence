# 合并相邻且无间隔的区域形成连续区域
# 合并XP-EHH相邻区域的R脚本

# 读取显著选择区域 (标准化后的结果)
regions <- read.table("significant_xpehh_regions.txt", header=FALSE)
colnames(regions) <- c("Chr", "Start", "End", "iHH1", "iHH2", "XPEHH", "NormXPEHH", "Pvalue")

# 按染色体和位置排序
regions <- regions[order(regions$Chr, regions$Start), ]

# 合并相邻区域的函数  
merge_regions <- function(reg_regions) {
    if(nrow(reg_regions) <= 1) return(reg_regions)
    
    merged <- reg_regions[1, ]
    
    for(i in 2:nrow(reg_regions)) {
        current <- reg_regions[i, ]
        last_merged <- merged[nrow(merged), ]
        
        # 如果当前区域与上一个合并区域相邻或重叠，则合并  
        if(current$Chr == last_merged$Chr && current$Start <= last_merged$End) {
            merged[nrow(merged), "End"] <- max(last_merged$End, current$End)
            merged[nrow(merged), "XPEHH"] <- max(last_merged$XPEHH, current$XPEHH)
            merged[nrow(merged), "NormXPEHH"] <- max(last_merged$NormXPEHH, current$NormXPEHH)
            merged[nrow(merged), "Pvalue"] <- min(last_merged$Pvalue, current$Pvalue)
        } else {
            merged <- rbind(merged, current)
        }
    }
    
    return(merged)
}

# 按染色体分组处理  
library(dplyr)
merged_regions <- regions %>%
    group_by(Chr) %>%
    do(merge_regions(.)) %>%
    ungroup()

# 保存合并后的结果  
write.table(merged_regions, "merged_significant_xpehh_regions.txt", 
            row.names=FALSE, quote=FALSE, sep="\t")

# 统计结果  
cat("XP-EHH显著选择区域统计:\n")
cat("原始显著区域数:", nrow(regions), "\n")
cat("合并后区域数:", nrow(merged_regions), "\n")
cat("平均区域大小:", mean(merged_regions$End - merged_regions$Start), "bp\n")