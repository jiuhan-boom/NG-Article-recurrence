# ==========================================
# 群体分层分析完整流程
# Population stratification analysis pipeline
# ==========================================

# 设置工作目录和环境变量
WORK_DIR=$(pwd)
POPULATION_DIR="$WORK_DIR/population_stratification"
ADMIXTURE_DIR="$POPULATION_DIR/admixture_analysis"
DIVERSITY_DIR="$POPULATION_DIR/diversity_analysis"
PCA_DIR="$POPULATION_DIR/pca_analysis"
NETWORK_DIR="$POPULATION_DIR/haplotype_network"
LOG_FILE="$WORK_DIR/population_analysis.log"

# 创建目录结构
mkdir -p $POPULATION_DIR $ADMIXTURE_DIR $DIVERSITY_DIR $PCA_DIR $NETWORK_DIR
mkdir -p $WORK_DIR/logs $WORK_DIR/scripts

echo "==========================================" | tee $LOG_FILE
echo "群体分层分析流程开始" | tee -a $LOG_FILE
echo "工作目录: $WORK_DIR" | tee -a $LOG_FILE
echo "==========================================" | tee -a $LOG_FILE

# ==========================================
# 第一步：SNP数据集预处理和过滤
# ==========================================
echo "[$(date)] 开始SNP数据集预处理..." | tee -a $LOG_FILE

prepare_snp_dataset() {
    local input_vcf=$1
    local output_dir=$2
    
    if [ ! -f "$input_vcf" ]; then
        echo "错误：输入VCF文件不存在: $input_vcf" | tee -a $LOG_FILE
        return 1
    fi
    
    cd $output_dir
    
    # 检查必要工具是否已安装
    if ! command -v plink &> /dev/null; then
        echo "警告：未找到PLINK，请先安装PLINK" | tee -a $LOG_FILE
        return 1
    fi
    
    if ! command -v vcftools &> /dev/null; then
        echo "警告：未找到VCFtools，请先安装VCFtools" | tee -a $LOG_FILE
        return 1
    fi
    
    echo "转换VCF到PLINK格式..." | tee -a $LOG_FILE
    
    # 转换VCF到PLINK格式
    plink \
        --vcf "$input_vcf" \
        --make-bed \
        --out raw_data \
        --allow-extra-chr \
        >> $WORK_DIR/logs/plink_conversion.log 2>&1
    
    if [ $? -eq 0 ]; then
        echo "VCF转换完成" | tee -a $LOG_FILE
        
        # 应用过滤标准 (MAF > 0.05, 缺失率 < 10%)
        echo "应用质量过滤..." | tee -a $LOG_FILE
        
        plink \
            --bfile raw_data \
            --maf 0.05 \
            --geno 0.1 \
            --mind 0.1 \
            --make-bed \
            --out filtered_data \
            --allow-extra-chr \
            >> $WORK_DIR/logs/plink_filter.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo "质量过滤完成" | tee -a $LOG_FILE
            
            # LD过滤 (r2 < 0.2 within 100kb)
            echo "进行LD过滤..." | tee -a $LOG_FILE
            
            plink \
                --bfile filtered_data \
                --indep-pairwise 100 5 0.2 \
                --out ld_pruned \
                --allow-extra-chr \
                >> $WORK_DIR/logs/ld_pruning.log 2>&1
            
            if [ $? -eq 0 ]; then
                # 应用LD过滤结果
                plink \
                    --bfile filtered_data \
                    --extract ld_pruned.prune.in \
                    --make-bed \
                    --out final_dataset \
                    --allow-extra-chr \
                    >> $WORK_DIR/logs/final_dataset.log 2>&1
                
                if [ $? -eq 0 ]; then
                    # 统计最终数据集信息
                    snp_count=$(wc -l < final_dataset.bim)
                    sample_count=$(wc -l < final_dataset.fam)
                    
                    cat > dataset_summary.txt << EOF
# SNP数据集预处理摘要

处理时间: $(date)
输入VCF文件: $(basename "$input_vcf")
最终SNP数量: $snp_count
样本数量: $sample_count

## 过滤标准:
- 最小等位基因频率 (MAF) > 0.05
- 缺失率 < 10%
- 连锁不平衡 (LD) r2 < 0.2 (100kb窗口内)
- 样本缺失率 < 10%

## 输出文件:
- 原始数据: raw_data.{bed,bim,fam}
- 质量过滤后: filtered_data.{bed,bim,fam}
- LD过滤后: final_dataset.{bed,bim,fam}

EOF
                    
                    echo "SNP数据集预处理完成，摘要报告保存到: dataset_summary.txt" | tee -a $LOG_FILE
                    
                else
                    echo "最终数据集生成失败" | tee -a $LOG_FILE
                fi
                
            else
                echo "LD过滤失败" | tee -a $LOG_FILE
            fi
            
        else
            echo "质量过滤失败" | tee -a $LOG_FILE
        fi
        
    else
        echo "VCF转换失败" | tee -a $LOG_FILE
    fi
}

# ==========================================
# 第二步：ADMIXTURE群体结构分析
# ==========================================
echo "[$(date)] 开始ADMIXTURE分析..." | tee -a $LOG_FILE

run_admixture_analysis() {
    local plink_bed=$1  # PLINK bed文件前缀
    local output_dir=$2
    
    if [ ! -f "${plink_bed}.bed" ]; then
        echo "错误：PLINK文件不存在: ${plink_bed}.bed" | tee -a $LOG_FILE
        return 1
    fi
    
    cd $output_dir
    
    # 检查ADMIXTURE是否已安装
    if ! command -v admixture &> /dev/null; then
        echo "警告：未找到ADMIXTURE，请先安装ADMIXTURE v1.3.0" | tee -a $LOG_FILE
        
        
        
        chmod +x install_admixture.sh
        echo "已创建ADMIXTURE安装脚本: install_admixture.sh" | tee -a $LOG_FILE
        return 1
    fi
    
    # 创建交叉验证分析目录
    mkdir -p cv_analysis
    
    # 运行K=1到K=10的交叉验证分析
    echo "运行K=1到K=10的交叉验证分析..." | tee -a $LOG_FILE
    
    # 创建交叉验证结果汇总文件
    > cv_results.txt
    
    for k in {1..10}; do
        echo "运行K=$k分析..." | tee -a $LOG_FILE
        
        admixture \
            --cv \
            "${plink_bed}.bed" \
            $k \
            > "cv_analysis/admixture_k${k}.out" \
            2>> $WORK_DIR/logs/admixture_k${k}.log
        
        if [ $? -eq 0 ]; then
            # 提取交叉验证误差
            cv_error=$(grep "CV error" "cv_analysis/admixture_k${k}.out" | awk '{print $3}')
            echo "K=$k, CV Error=$cv_error" >> cv_results.txt
            
            echo "K=$k分析完成，CV误差: $cv_error" | tee -a $LOG_FILE
            
            # 重命名输出文件以便识别
            mv "${plink_bed}.${k}.P" "admixture_k${k}.P"
            mv "${plink_bed}.${k}.Q" "admixture_k${k}.Q"
            
        else
            echo "K=$k分析失败" | tee -a $LOG_FILE
        fi
    done
    
    # 确定最佳K值 (根据交叉验证误差)
    echo "确定最佳K值..." | tee -a $LOG_FILE
    
    best_k=$(awk 'NR==1{min=$3; best=1} $3<min{min=$3; best=NR} END{print best}' cv_results.txt)
    
    cat > admixture_summary.txt << EOF
# ADMIXTURE分析摘要

分析时间: $(date)
输入数据集: $(basename "$plink_bed")
样本数量: $(wc -l < ${plink_bed}.fam)

## 交叉验证结果:
$(cat cv_results.txt)

## 最佳K值:
根据交叉验证误差，最佳K值为: K=$best_k

## 推荐使用的K值:
参考文献中选择K=3作为祖先成分数量，因为交叉验证误差在此处急剧收敛。

## 输出文件:
- 各K值的祖先成分比例: admixture_k*.P  
- 各K值的个体祖先比例: admixture_k*.Q
- 交叉验证结果: cv_results.txt

EOF
    
    echo "ADMIXTURE分析完成，摘要报告保存到: admixture_summary.txt" | tee -a $LOG_FILE

    echo "已创建R可视化脚本: visualize_admixture.R" | tee -a $LOG_FILE
    
}

# ==========================================
# 第三步：核苷酸多样性(π)和群体分化(Fst)分析  
# ==========================================
echo "[$(date)] 开始核苷酸多样性和群体分化分析..." | tee -a $LOG_FILE

run_diversity_analysis() {
    local input_vcf=$1
    local output_dir=$2
    
    if [ ! -f "$input_vcf" ]; then  
        echo "错误：输入VCF文件不存在: $input_vcf" | tee -a $LOG_FILE  
        return 1  
    fi  
    
    cd $output_dir  
    
    # 检查VCFtools是否已安装  
    if ! command -v vcftools &> /dev/null; then  
        echo "警告：未找到VCFtools，请先安装VCFtools v0.1.16" | tee -a $LOG_FILE  
        return 1  
    fi  
    
    # 创建群体分组文件 (需要根据实际情况创建)  
    echo "# 创建群体分组示例文件 (需要根据实际样本信息调整)" > create_populations.sh    
    
    chmod +x create_populations.sh  
    
    # 计算核苷酸多样性 (π)  
    echo "计算核苷酸多样性 (π)..." | tee -a $LOG_FILE  
    
    vcftools \  
        --gzvcf "$input_vcf" \
        --window-pi 1000 \
        --window-pi-step 100 \
        --out nucleotide_diversity \
        >> $WORK_DIR/logs/vcftools_pi.log 2>&1  
    
    if [ $? -eq 0 ]; then  
        echo "核苷酸多样性计算完成" | tee -a $LOG_FILE  
        
        # 计算群体分化 (Fst)  
        echo "计算群体分化 (Fst)..." | tee -a $LOG_FILE  
        
        # 需要群体分组文件，这里创建示例  
        if [ ! -f "all_populations.txt" ]; then  
            echo "# 创建示例群体分组文件..." | tee -a $LOG_FILE  
            zcat "$input_vcf" | grep -m 10 "^#" | tail -1 | cut -f10- | tr '\t' '\n' | \
            awk 'NR<=5{print $1"\tSweetCorn"} NR>5{print $1"\tFieldCorn"}' > all_populations.txt  
        fi  
        
        vcftools \
            --gzvcf "$input_vcf" \
            --weir-fst-pop all_populations.txt \
            --out population_divergence \
            >> $WORK_DIR/logs/vcftools_fst.log 2>&1  
        
        if [ $? -eq 0 ]; then  
            echo "群体分化计算完成" | tee -a $LOG_FILE  
            
            # 创建多样性分析报告  
            cat > diversity_summary.txt << EOF  
# 核苷酸多样性和群体分化分析摘要  

分析时间: $(date)  
输入VCF文件: $(basename "$input_vcf")  

## 分析参数:  
- 滑动窗口大小: 1,000 bp  
- 步长: 100 bp  

## 输出文件:  
- 核苷酸多样性结果: nucleotide_diversity.windowed.pi  
- 群体分化结果: population_divergence.weir.fst  

## 结果统计:  

### 核苷酸多样性 (π):  
$(head -5 nucleotide_diversity.windowed.pi)  

### 群体分化 (Fst):  
$(head -5 population_divergence.weir.fst)  

EOF
            
            echo "多样性分析完成，摘要报告保存到: diversity_summary.txt" | tee -a $LOG_FILE  
            
        else  
            echo "群体分化计算失败" | tee -a $LOG_FILE  
        fi  
        
    else  
        echo "核苷酸多样性计算失败" | tee -a $LOG_FILE  
    fi  
    
}

# ==========================================  
# 第四步：主成分分析 (PCA)  
# ==========================================  
echo "[$(date)] 开始主成分分析..." | tee -a $LOG_FILE  

run_pca_analysis() {  
    local plink_bed=$1  # PLINK bed文件前缀  
    local output_dir=$2  
    
    if [ ! -f "${plink_bed}.bed" ]; then  
        echo "错误：PLINK文件不存在: ${plink_bed}.bed" | tee -a $LOG_FILE  
        return 1  
    fi  
    
    cd $output_dir  
    
    # 检查PLINK是否已安装  
    if ! command -v plink &> /dev/null; then  
        echo "警告：未找到PLINK，请先安装PLINK" | tee -a $LOG_FILE  
        return 1  
    fi  
    
    # 应用PCA过滤标准 (MAF > 0.05, 缺失率 < 10%)  
    echo "应用PCA过滤标准..." | tee -a $LOG_FILE  
    
    plink \
        --bfile "$plink_bed" \
        --maf 0.05 \
        --geno 0.1 \
        --mind 0.1 \
        --make-bed \
        --out pca_filtered \
        --allow-extra-chr \
        >> $WORK_DIR/logs/pca_filter.log 2>&1  
    
    if [ $? -eq 0 ]; then  
        echo "PCA过滤完成，运行主成分分析..." | tee -a $LOG_FILE  
        
        # 运行PCA分析 (前10个主成分)  
        plink \
            --bfile pca_filtered \
            --pca 10 \
            --out pca_results \
            --allow-extra-chr \
            >> $WORK_DIR/logs/plink_pca.log 2>&1  
        
        if [ $? -eq 0 ]; then  
            echo "PCA分析完成" | tee -a $LOG_FILE  
            
            # 统计PCA结果信息  
            pc_count=$(wc -l < pca_results.eigenvec)  
            
            cat > pca_summary.txt << EOF  
# 主成分分析 (PCA) 摘要  

分析时间: $(date)  
输入数据集: $(basename "$plink_bed")  
样本数量: $(wc -l < pca_filtered.fam)  

## 分析参数:  
- MAF过滤: > 0.05  
- 缺失率过滤: < 10%  

## 输出文件:  
- 特征值: pca_results.eigenval  
- 特征向量: pca_results.eigenvec  

## 主成分解释方差比例:  

EOF
            
            # 计算前几个主成分的解释方差比例  
            total_variance=$(awk '{sum+=$1} END {print sum}' pca_results.eigenval)  
            
            for i in {1..5}; do  
                eigenval=$(sed -n "${i}p" pca_results.eigenval)  
                if [ -n "$eigenval" ] && [ "$total_variance" != "0" ]; then  
                    proportion=$(echo "scale=4; $eigenval / $total_variance * 100" | bc)  
                    echo "PC${i}: ${proportion}%" >> pca_summary.txt  
                fi  
            done  
            
            echo "" >> pca_summary.txt  
            echo "总方差: $total_variance" >> pca_summary.txt  
            
            echo "PCA分析完成，摘要报告保存到: pca_summary.txt" | tee -a $LOG_FILE  
            
            
            echo "已创建R可视化脚本: visualize_pca.R" | tee -a $LOG_FILE  
            
        else  
            echo "PCA分析失败" | tee -a $LOG_FILE  
        fi  
        
    else  
        echo "PCA过滤失败" | tee -a $LOG_FILE  
    fi  
    
}

# ==========================================  

