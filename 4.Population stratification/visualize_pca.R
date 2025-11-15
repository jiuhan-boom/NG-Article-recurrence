library(ggplot2)  

# 读取PCA结果文件 (特征向量)  
pca_data <- read.table("pca_results.eigenvec", header=FALSE)  

# 假设前两列是家族ID和个体ID，后面是主成分值  
colnames(pca_data) <- c("Family", "Individual", paste0("PC", 1:(ncol(pca_data)-2)))  

# 创建PCA图 (PC1 vs PC2)  
p <- ggplot(pca_data, aes(x=PC1, y=PC2)) +  
    geom_point(size=2, alpha=0.7) +  
    labs(title="Principal Component Analysis",  
         x=paste0("PC1 (", round(summary(lm(PC2~PC1, data=pca_data))$r.squared*100, 2), "% variance)"),   
         y=paste0("PC2 (", round(summary(lm(PC3~PC2, data=pca_data))$r.squared*100, 2), "% variance)")) +  
    theme_minimal()  

print(p)  

# 保存图形  
ggsave("pca_plot.png", width=8, height=6)  

# 创建群体着色的PCA图 (如果有群体信息)  
if(file.exists("population_info.txt")) {  
    pop_info <- read.table("population_info.txt", header=TRUE)  
    pca_data$Population <- pop_info$Population[match(pca_data$Individual, pop_info$Sample)]  
    
    p2 <- ggplot(pca_data, aes(x=PC1, y=PC2, color=Population)) +  
        geom_point(size=2, alpha=0.7) +  
        labs(title="PCA by Population", x="PC1", y="PC2") +  
        theme_minimal() +  
        scale_color_brewer(palette="Set1")  
    
    print(p2)  
    ggsave("pca_by_population.png", width=8, height=6)  
}  
