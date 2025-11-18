# ==========================================
# 统计分析完整流程
# Statistics analysis pipeline
# ==========================================

# 设置工作目录和环境变量
WORK_DIR=$(pwd)
STATISTICS_DIR="$WORK_DIR/statistics_analysis"
HERITABILITY_DIR="$STATISTICS_DIR/heritability"
SNP_HERITABILITY_DIR="$STATISTICS_DIR/snp_heritability"
EQTL_DIR="$STATISTICS_DIR/eqtl_analysis"
REGRESSION_DIR="$STATISTICS_DIR/regression_analysis"
ENRICHMENT_DIR="$STATISTICS_DIR/enrichment_analysis"
LOG_FILE="$WORK_DIR/statistics_analysis.log"

# 创建目录结构
mkdir -p $STATISTICS_DIR $HERITABILITY_DIR $SNP_HERITABILITY_DIR $EQTL_DIR $REGRESSION_DIR $ENRICHMENT_DIR
mkdir -p $WORK_DIR/logs $WORK_DIR/scripts

echo "==========================================" | tee $LOG_FILE
echo "统计分析流程开始" | tee -a $LOG_FILE
echo "工作目录: $WORK_DIR" | tee -a $LOG_FILE
echo "==========================================" | tee -a $LOG_FILE

# ==========================================
# 第一步：广义遗传力分析 (Broad-sense heritability)
# ==========================================
echo "[$(date)] 开始广义遗传力分析..." | tee -a $LOG_FILE

run_heritability_analysis() {
    local phenotype_file=$1      # 表型数据文件
    local genotype_file=$2       # 基因型数据文件 (可选)
    local output_dir=$3
    
    cd $output_dir
    
    # 运行R脚本
    Rscript heritability_analysis.R >> $WORK_DIR/logs/heritability_analysis.log 2>&1
    
    if [ $? -eq 0 ]; then
        echo "广义遗传力分析完成" | tee -a $LOG_FILE
    else
        echo "广义遗传力分析失败，请检查日志文件" | tee -a $LOG_FILE
    fi
    
}

# ==========================================
# 第二步：SNP遗传力分析 (使用GCTA)
# ==========================================
echo "[$(date)] 开始SNP遗传力分析..." | tee -a $LOG_FILE

run_snp_heritability_analysis() {
    local genotype_file=$1       # VCF格式基因型文件
    local phenotype_file=$2      # 表型数据文件
    local output_dir=$3
    
    cd $output_dir
    
    # 检查GCTA是否已安装
    if ! command -v gcta64 &> /dev/null; then
        echo "警告：未找到GCTA，请先安装GCTA v1.94.1" | tee -a $LOG_FILE
        
          
        chmod +x install_gcta.sh
        echo "已创建GCTA安装脚本: install_gcta.sh" | tee -a $LOG_FILE
        return 1
    fi
    
    # 检查VCFtools是否已安装 (用于格式转换)
    if ! command -v vcftools &> /dev/null; then
        echo "警告：未找到VCFtools，请先安装VCFtools" | tee -a $LOG_FILE
        return 1
    fi
    
    # 转换VCF到PLINK格式以满足GCTA要求
    echo "转换VCF到PLINK格式..." | tee -a $LOG_FILE
    
    vcftools \
        --gzvcf "$genotype_file" \
        --plink \
        --out genotype_plink \
        >> $WORK_DIR/logs/vcftools_conversion.log 2>&1
    
    if [ $? -eq 0 ]; then
        echo "格式转换完成" | tee -a $LOG_FILE
        
        # 准备表型文件 (GCTA格式)
        echo "准备GCTA格式表型文件..." | tee -a $LOG_FILE
            
        Rscript prepare_phenotype.R >> $WORK_DIR/logs/prepare_phenotype.log 2>&1
        
        # 运行GCTA SNP遗传力分析
        echo "运行GCTA SNP遗传力分析..." | tee -a $LOG_FILE
        
        # GCTA REML分析估计SNP遗传力
        gcta64 \
            --bfile genotype_plink \
            --pheno gcta_phenotype.txt \
            --reml \
            --out snp_heritability_results \
            >> $WORK_DIR/logs/gcta_analysis.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo "SNP遗传力分析完成" | tee -a $LOG_FILE
            
            # 创建SNP遗传力分析报告
            cat > snp_heritability_summary.txt << EOF
# SNP遗传力分析摘要

分析时间: $(date)
基因型文件: $(basename "$genotype_file")
表型文件: $(basename "$phenotype_file")

## 分析参数:
- 使用GCTA v1.94.1进行REML分析
- 转换VCF到PLINK格式满足GCTA要求

## 输出文件:
- PLINK格式基因型: genotype_plink.{bed,bim,fam}
- GCTA表型文件: gcta_phenotype.txt
- SNP遗传力结果: snp_heritability_results.*

EOF
            
            echo "SNP遗传力分析完成，摘要报告保存到: snp_heritability_summary.txt" | tee -a $LOG_FILE
            
        else
            echo "SNP遗传力分析失败，请检查日志文件" | tee -a $LOG_FILE
        fi
        
    else
        echo "VCF到PLINK格式转换失败" | tee -a $LOG_FILE
    fi
    
}

# ==========================================
# 第三步：eQTL热点识别和置换检验
# ==========================================
echo "[$(date)] 开始eQTL热点识别和置换检验..." | tee -a $LOG_FILE

run_eqtl_hotspot_analysis() {
    local eqtl_results=$1        # eQTL结果文件
    local output_dir=$2
}

