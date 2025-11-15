# 合并相邻且无间隔的区域形成连续区域
# 合并相邻区域的R脚本

# 读取候选选择清除区域
sweeps <- read.table("candidate_sweeps_xpclr.txt", header=FALSE)
colnames(sweeps) <- c("Chr", "Start", "End", "XPCLR")

# 按染色体和位置排序
sweeps <- sweeps[order(sweeps$Chr, sweeps$Start), ]

# 合并相邻区域的函数
merge_regions <- function(regions) {
    if(nrow(regions) <= 1) return(regions)
    
    merged <- regions[1, ]
    
    for(i in 2:nrow(regions)) {
        current <- regions[i, ]
        last_merged <- merged[nrow(merged), ]
        
        # 如果当前区域与上一个合并区域相邻或重叠，则合并
        if(current$Chr == last_merged$Chr && current$Start <= last_merged$End) {
            merged[nrow(merged), "End"] <- max(last_merged$End, current$End)
            merged[nrow(merged), "XPCLR"] <- max(last_merged$XPCLR, current$XPCLR)
        } else {
            merged <- rbind(merged, current)
        }
    }
    
    return(merged)
}

# 按染色体分组处理
library(dplyr)
merged_sweeps <- sweeps %>%
    group_by(Chr) %>%
    do(merge_regions(.)) %>%
    ungroup()

# 保存合并后的结果
write.table(merged_sweeps, "merged_candidate_sweeps_xpclr.txt", 
            row.names=FALSE, quote=FALSE, sep="\t")

# 统计结果
cat("XP-CLR候选选择清除区域统计:\n")
cat("原始候选区域数:", nrow(sweeps), "\n")
cat("合并后区域数:", nrow(merged_sweeps), "\n")
cat("平均区域大小:", mean(merged_sweeps$End - merged_sweeps$Start), "bp\n")
    