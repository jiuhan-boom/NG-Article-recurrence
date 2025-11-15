# ==========================================
# 全基因组选择清除扫描完整流程
# Genome-wide selective sweep scan pipeline
# ==========================================

# 设置工作目录和环境变量
WORK_DIR=$(pwd)
SELECTIVE_DIR="$WORK_DIR/selective_sweep"
XPCLR_DIR="$SELECTIVE_DIR/xpclr_analysis"
XPEHH_DIR="$SELECTIVE_DIR/xpehh_analysis"
RESULTS_DIR="$SELECTIVE_DIR/results"
LOG_FILE="$WORK_DIR/selective_sweep_analysis.log"

# 创建目录结构
mkdir -p $SELECTIVE_DIR $XPCLR_DIR $XPEHH_DIR $RESULTS_DIR
mkdir -p $WORK_DIR/logs $WORK_DIR/scripts

echo "==========================================" | tee $LOG_FILE
echo "全基因组选择清除扫描流程开始" | tee -a $LOG_FILE
echo "工作目录: $WORK_DIR" | tee -a $LOG_FILE
echo "==========================================" | tee -a $LOG_FILE

# ==========================================
# 第一步：XP-CLR分析准备和运行
# ==========================================
echo "[$(date)] 开始XP-CLR分析..." | tee -a $LOG_FILE

run_xpclr_analysis() {
    local sweet_corn_vcf=$1      # 甜玉米VCF文件 (测试群体)
    local field_corn_vcf=$2      # 田间玉米VCF文件 (参考群体)
    local genetic_map=$3         # 遗传图谱文件 (可选)
    local output_dir=$4
    
    if [ ! -f "$sweet_corn_vcf" ] || [ ! -f "$field_corn_vcf" ]; then
        echo "错误：输入VCF文件不存在" | tee -a $LOG_FILE
        return 1
    fi
    
    cd $output_dir
    
    # 检查XP-CLR是否已安装
    if ! command -v xpclr &> /dev/null; then
        echo "警告：未找到XP-CLR，请先安装XP-CLR v1.0" | tee -a $LOG_FILE
        
        chmod +x install_xpclr.sh
        echo "已创建XP-CLR安装脚本: install_xpclr.sh" | tee -a $LOG_FILE
        return 1
    fi
    
    # 准备群体文件 (需要根据实际数据创建)
    echo "准备群体文件..." | tee -a $LOG_FILE
    
    # 创建甜玉米群体文件 (假设前半部分样本是甜玉米)
    zcat "$sweet_corn_vcf" | grep -m 1 "^#" | cut -f10- | tr '\t' '\n' > sweet_corn_samples.txt
    
    # 创建田间玉米群体文件 (假设后半部分样本是田间玉米)
    zcat "$field_corn_vcf" | grep -m 1 "^#" | cut -f10- | tr '\t' '\n' > field_corn_samples.txt
    
    # 转换VCF到XP-CLR格式
    echo "转换VCF到XP-CLR格式..." | tee -a $LOG_FILE
    
    # 转换两个群体的VCF文件
    python3 convert_vcp_xpclr.py "$sweet_corn_vcf" "sweet_corn" "sweet_corn_samples.txt"
    python3 convert_vcp_xpclr.py "$field_corn_vcf" "field_corn" "field_corn_samples.txt"
    
    # 创建遗传图谱文件 (如果没有提供)
    if [ ! -f "$genetic_map" ]; then
        echo "创建示例遗传图谱文件..." | tee -a $LOG_FILE
        
        # 从VCF文件提取染色体和位置信息创建遗传图谱
        zcat "$sweet_corn_vcf" | grep -v "^#" | head -1000 | \
        awk '{print $1"\t"$2"\t"$2/1000000}' > genetic_map.txt
        
        genetic_map="genetic_map.txt"
    fi
    
    # 运行XP-CLR分析
    echo "运行XP-CLR分析..." | tee -a $LOG_FILE
    
    # XP-CLR参数设置：
    # - 滑动窗口大小: 0.5 cM
    # - 网格大小: 20 kb  
    # - 窗口内最大SNP数: 200
    # - 相关系数截断值: 0.70
    
    xpclr \
        -xpclr sweet_corn.xpclr.txt \
        -ref field_corn.xpclr.txt \
        -map "$genetic_map" \
        -out xpclr_results.txt \
        -size 0.5 \
        -grid 20000 \
        -maxsnps 200 \
        -corr 0.70 \
        >> $WORK_DIR/logs/xpclr_analysis.log 2>&1
    
    if [ $? -eq 0 ]; then
        echo "XP-CLR分析完成" | tee -a $LOG_FILE
        
        # 识别候选选择清除区域 (top 5% highest XP-CLR values)
        echo "识别候选选择清除区域..." | tee -a $LOG_FILE
        
        # 计算XP-CLR值的95百分位数作为阈值
        threshold=$(awk 'NR>1 {print $4}' xpclr_results.txt | sort -n | awk 'NR==int(NR*0.95)' RS='\n' FS='\n')
        
        # 提取top 5%的区域
        awk -v thresh="$threshold" 'NR>1 && $4>=thresh {print}' xpclr_results.txt > candidate_sweeps_xpclr.txt
        
        Rscript merge_adjacent_regions.R >> $WORK_DIR/logs/merge_regions.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo "区域合并完成" | tee -a $LOG_FILE
            
            # 创建XP-CLR分析报告
            original_count=$(wc -l < candidate_sweeps_xpclr.txt)
            merged_count=$(wc -l < merged_candidate_sweeps_xpclr.txt)
            
            cat > xpclr_summary.txt << EOF
# XP-CLR分析摘要

分析时间: $(date)
测试群体 (甜玉米): $(basename "$sweet_corn_vcf")
参考群体 (田间玉米): $(basename "$field_corn_vcf")

## 分析参数:
- 滑动窗口大小: 0.5 cM
- 网格大小: 20 kb
- 窗口内最大SNP数: 200
- 相关系数截断值: 0.70
- 显著性阈值: top 5% (XP-CLR值 >= $(printf "%.3f" $threshold))

## 结果统计:
- 原始候选区域数: $original_count
- 合并后连续区域数: $merged_count

## 输出文件:
- XP-CLR结果: xpclr_results.txt
- 候选选择清除区域: candidate_sweeps_xpclr.txt
- 合并后区域: merged_candidate_sweeps_xpclr.txt

EOF
            
            echo "XP-CLR分析完成，摘要报告保存到: xpclr_summary.txt" | tee -a $LOG_FILE
            
        else
            echo "区域合并失败" | tee -a $LOG_FILE
        fi
        
    else
        echo "XP-CLR分析失败，请检查日志文件" | tee -a $LOG_FILE
    fi
    
}

# ==========================================
# 第二步：XP-EHH分析准备和运行
# ==========================================
echo "[$(date)] 开始XP-EHH分析..." | tee -a $LOG_FILE

run_xpehh_analysis() {
    local sweet_corn_vcf=$1      # 甜玉米VCF文件 (群体1)
    local field_corn_vcf=$2      # 田间玉米VCF文件 (群体2)
    local output_dir=$3
    
    if [ ! -f "$sweet_corn_vcf" ] || [ ! -f "$field_corn_vcf" ]; then
        echo "错误：输入VCF文件不存在" | tee -a $LOG_FILE
        return 1
    fi
    
    cd $output_dir
    
    # 检查selscan是否已安装
    if ! command -v selscan &> /dev/null; then
        echo "警告：未找到selscan，请先安装selscan v1.3.0" | tee -a $LOG_FILE
        
        chmod +x install_selscan.sh
        echo "已创建selscan安装脚本: install_selscan.sh" | tee -a $LOG_FILE
        return 1
    fi
    
    # 准备群体文件和转换格式
    echo "准备群体文件和格式转换..." | tee -a $LOG_FILE
    
    # 转换两个群体的VCF文件到hapmap格式
    python3 convert_to_hapmap.py "$sweet_corn_vcf" "sweet_corn" "SweetCorn"
    python3 convert_to_hapmap.py "$field_corn_vcf" "field_corn" "FieldCorn"
    
    # 运行XP-EHH分析 (使用20kb非重叠窗口)
    echo "运行XP-EHH分析..." | tee -a $LOG_FILE
    
    selscan \
        --xpehh \
        --hap sweet_corn.hapmap \
        --ref field_corn.hapmap \
        --out xpehh_results \
        --window 20000 \
        --keep-low-freq \
        >> $WORK_DIR/logs/selscan_xpehh.log 2>&1
    
    if [ $? -eq 0 ]; then
        echo "XP-EHH分析完成" | tee -a $LOG_FILE
        
        # 标准化XP-EHH值 (使用norm工具)
        echo "标准化XP-EHH值..." | tee -a $LOG_FILE
        
        norm \
            --xpehh \
            --files xpehh_results.xpehh.out \
            >> $WORK_DIR/logs/norm_xpehh.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo "XP-EHH标准化完成" | tee -a $LOG_FILE
            
            # 识别显著选择区域 (P < 0.05)
            echo "识别显著选择区域..." | tee -a $LOG_FILE
            
            # 提取P值小于0.05的区域
            awk '$8 < 0.05 {print}' xpehh_results.xpehh.out.norm > significant_xpehh_regions.txt
                    
            Rscript merge_xpehh_regions.R >> $WORK_DIR/logs/merge_xpehh.log 2>&1
            
            if [ $? -eq 0 ]; then
                echo "XP-EHH区域合并完成" | tee -a $LOG_FILE
                
                # 创建XP-EHH分析报告  
                original_xpehh_count=$(wc -l < significant_xpehh_regions.txt)  
                merged_xpehh_count=$(wc -l < merged_significant_xpehh_regions.txt)  
                
                cat > xpehh_summary.txt << EOF  
# XP-EHH分析摘要  

分析时间: $(date)  
群体1 (甜玉米): $(basename "$sweet_corn_vcf")  
群体2 (田间玉米): $(basename "$field_corn_vcf")  

## 分析参数:  
- 窗口大小: 20 kb 非重叠窗口  
- 显著性阈值: P < 0.05  

## 结果统计:  
- 原始显著区域数: $original_xpehh_count  
- 合并后连续区域数: $merged_xpehh_count  

## 输出文件:  
- XP-EHH结果: xpehh_results.xpehh.out  
- 标准化结果: xpehh_results.xpehh.out.norm  
- 显著选择区域: significant_xpehh_regions.txt  
- 合并后区域: merged_significant_xpehh_regions.txt  

EOF
                
                echo "XP-EHH分析完成，摘要报告保存到: xpehh_summary.txt" | tee -a $LOG_FILE  
                
            else  
                echo "XP-EHH区域合并失败" | tee -a $LOG_FILE  
            fi  
            
        else  
            echo "XP-EHH标准化失败" | tee -a $LOG_FILE  
        fi  
        
    else  
        echo "XP-EHH分析失败，请检查日志文件" | tee -a $LOG_FILE  
    fi  
    
}

# ==========================================  
# 第三步：结果整合和可视化  
# ==========================================  
echo "[$(date)] 开始结果整合和可视化..." | tee -a $LOG_FILE  

integrate_results() {  
    local xpclr_results=$1  
    local xpehh_results=$2  
    local output_dir=$3  
    
    cd $output_dir  
    
    # 整合两种方法的结果  
    echo "整合XP-CLR和XP-EHH结果..." | tee -a $LOG_FILE  
}

