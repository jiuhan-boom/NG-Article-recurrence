# ==========================================
# 重测序和变异检测完整流程
# Resequencing and variant calling pipeline
# ==========================================

# 设置工作目录和环境变量
WORK_DIR=$(pwd)
RESEQ_DIR="$WORK_DIR/resequencing"
VARIANT_DIR="$WORK_DIR/variant_calling"
POPULATION_DIR="$WORK_DIR/population_analysis"
LOG_FILE="$WORK_DIR/resequencing_analysis.log"

# 创建目录结构
mkdir -p $RESEQ_DIR $VARIANT_DIR $POPULATION_DIR
mkdir -p $RESEQ_DIR/raw_data $RESEQ_DIR/quality_control
mkdir -p $VARIANT_DIR/alignment $VARIANT_DIR/calling $VARIANT_DIR/filtering
mkdir -p $POPULATION_DIR/pca $POPULATION_DIR/tree $POPULATION_DIR/admixture
mkdir -p $WORK_DIR/logs $WORK_DIR/scripts

echo "==========================================" | tee $LOG_FILE
echo "重测序和变异检测流程开始" | tee -a $LOG_FILE
echo "工作目录: $WORK_DIR" | tee -a $LOG_FILE
echo "==========================================" | tee -a $LOG_FILE

# ==========================================
# 第一步：数据预处理和质量控制
# ==========================================
echo "[$(date)] 开始数据预处理和质量控制..." | tee -a $LOG_FILE

run_quality_control() {
    local input_dir=$1
    local output_dir=$2
    
    cd $output_dir
    
    # 检查fastp是否已安装
    if ! command -v fastp &> /dev/null; then
        echo "警告：未找到fastp，请先安装fastp" | tee -a $LOG_FILE
        return 1
    fi
    
    # 处理所有FASTQ文件
    echo "处理FASTQ文件进行质量控制..." | tee -a $LOG_FILE
    
    # 创建质量控制报告
    cat > qc_summary.txt << EOF
# 重测序数据质量控制报告

分析时间: $(date)
输入目录: $input_dir
输出目录: $output_dir

## 处理的样本列表:
EOF
    
    # 处理配对末端数据
    for fastq1 in $input_dir/*_1.fastq.gz $input_dir/*_R1.fastq.gz $input_dir/*_R1_001.fastq.gz; do
        if [ -f "$fastq1" ]; then
            # 获取对应的配对文件
            fastq2=${fastq1/_1.fastq.gz/_2.fastq.gz}
            fastq2=${fastq2/_R1.fastq.gz/_R2.fastq.gz}
            fastq2=${fastq2/_R1_001.fastq.gz/_R2_001.fastq.gz}
            
            if [ -f "$fastq2" ]; then
                sample_name=$(basename "$fastq1" | sed 's/_1\.fastq\.gz$//' | sed 's/_R1\.fastq\.gz$//' | sed 's/_R1_001\.fastq\.gz$//')
                
                echo "处理样本: $sample_name" | tee -a $LOG_FILE
                echo "- $sample_name" >> qc_summary.txt
                
                # 运行fastp质量控制
                fastp \
                    --in1 "$fastq1" \
                    --in2 "$fastq2" \
                    --out1 "${sample_name}_clean_1.fastq.gz" \
                    --out2 "${sample_name}_clean_2.fastq.gz" \
                    --json "${sample_name}_fastp.json" \
                    --html "${sample_name}_fastp.html" \
                    --thread 8 \
                    --average_qual 20 \
                    --length_required 100 \
                    --cut_front \
                    --cut_tail \
                    >> $WORK_DIR/logs/fastp_${sample_name}.log 2>&1
                
                if [ $? -eq 0 ]; then
                    echo "样本 $sample_name 质量控制完成" | tee -a $LOG_FILE
                else
                    echo "样本 $sample_name 质量控制失败" | tee -a $LOG_FILE
                fi
            fi
        fi
    done
    
    echo "质量控制完成，报告保存到: qc_summary.txt" | tee -a $LOG_FILE
}

# ==========================================
# 第二步：参考基因组准备和索引构建
# ==========================================
echo "[$(date)] 开始参考基因组准备..." | tee -a $LOG_FILE

prepare_reference_genome() {
    local ref_genome=$1
    local output_dir=$2
    
    if [ ! -f "$ref_genome" ]; then
        echo "错误：参考基因组文件不存在: $ref_genome" | tee -a $LOG_FILE
        return 1
    fi
    
    cd $output_dir
    
    # 检查BWA是否已安装
    if ! command -v bwa &> /dev/null; then
        echo "警告：未找到BWA" | tee -a $LOG_FILE
        return 1
    fi
    
    # 检查samtools是否已安装
    if ! command -v samtools &> /dev/null; then
        echo "警告：未找到samtools" | tee -a $LOG_FILE
        return 1
    fi
    
    # 检查picard是否已安装
    if ! command -v picard &> /dev/null; then
        echo "警告：未找到picard" | tee -a $LOG_FILE
    fi
    
    echo "构建BWA索引..." | tee -a $LOG_FILE
    bwa index "$ref_genome" > bwa_index.log 2>&1 &
    
    echo "构建samtools索引..." | tee -a $LOG_FILE
    samtools faidx "$ref_genome" > samtools_faidx.log 2>&1 &
    
    echo "创建序列字典..." | tee -a $LOG_FILE
    if command -v picard &> /dev/null; then
        picard CreateSequenceDictionary \
            R="$ref_genome" \
            O="${ref_genome%.fasta}.dict" \
            > picard_dict.log 2>&1 &
    fi
    
    # 等待所有后台进程完成
    wait
    
    echo "参考基因组准备完成" | tee -a $LOG_FILE
    
    # 创建参考基因组信息文件
    cat > reference_info.txt << EOF
# 参考基因组信息

基因组文件: $(basename "$ref_genome")
基因组大小: $(ls -lh "$ref_genome" | awk '{print $5}')
染色体数量: $(grep -c "^>" "$ref_genome")
总长度: $(awk '/^>/ {next} {sum+=length($0)} END {print sum}' "$ref_genome")

## 索引文件:
- BWA索引: ${ref_genome}.amb, ${ref_genome}.ann, ${ref_genome}.bwt, ${ref_genome}.pac, ${ref_genome}.sa
- samtools索引: ${ref_genome}.fai
- 序列字典: ${ref_genome%.fasta}.dict

EOF
    
    echo "参考基因组信息已保存到: reference_info.txt" | tee -a $LOG_FILE
}

# ==========================================
# 第三步：序列比对 (BWA-MEM)
# ==========================================
echo "[$(date)] 开始序列比对..." | tee -a $LOG_FILE

run_sequence_alignment() {
    local clean_data_dir=$1
    local ref_genome=$2
    local output_dir=$3
    
    if [ ! -f "$ref_genome" ]; then
        echo "错误：参考基因组文件不存在: $ref_genome" | tee -a $LOG_FILE
        return 1
    fi
    
    cd $output_dir
    
    # 检查BWA索引是否存在
    if [ ! -f "${ref_genome}.bwt" ]; then
        echo "错误：BWA索引不存在" | tee -a $LOG_FILE
        return 1
    fi
    
    # 处理所有清理后的数据
    echo "开始序列比对..." | tee -a $LOG_FILE
    
    for fastq1 in $clean_data_dir/*_clean_1.fastq.gz; do
        if [ -f "$fastq1" ]; then
            fastq2=${fastq1/_clean_1.fastq.gz/_clean_2.fastq.gz}
            sample_name=$(basename "$fastq1" | sed 's/_clean_1\.fastq\.gz$//')
            
            if [ -f "$fastq2" ]; then
                echo "比对样本: $sample_name" | tee -a $LOG_FILE
                
                # BWA-MEM比对
                bwa mem \
                    -t 16 \
                    -M \
                    -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA" \
                    "$ref_genome" \
                    "$fastq1" "$fastq2" \
                    > "${sample_name}.sam" \
                    2>> $WORK_DIR/logs/bwa_${sample_name}.log
                
                if [ $? -eq 0 ]; then
                    echo "样本 $sample_name 比对完成，生成SAM文件" | tee -a $LOG_FILE
                    
                    # SAM转BAM并排序
                    samtools view -bS "${sample_name}.sam" | \
                    samtools sort -@ 8 -o "${sample_name}.sorted.bam" -
                    
                    if [ $? -eq 0 ]; then
                        echo "样本 $sample_name BAM文件生成完成" | tee -a $LOG_FILE
                        
                        # 删除临时SAM文件以节省空间
                        rm "${sample_name}.sam"
                        
                        # 索引BAM文件
                        samtools index "${sample_name}.sorted.bam"
                        
                        # 统计比对结果
                        samtools flagstat "${sample_name}.sorted.bam" > "${sample_name}_alignment_stats.txt"
                        
                        echo "样本 $sample_name 比对统计完成" | tee -a $LOG_FILE
                        
                    else
                        echo "样本 $sample_name BAM转换失败" | tee -a $LOG_FILE
                    fi
                    
                else
                    echo "样本 $sample_name 比对失败" | tee -a $LOG_FILE
                fi
            fi
        fi
    done
    
    echo "序列比对完成" | tee -a $LOG_FILE
    
    # 创建比对统计汇总报告
    cat > alignment_summary.txt << EOF
# 序列比对统计汇总

分析时间: $(date)
参考基因组: $(basename "$ref_genome")
处理样本数: $(ls *_alignment_stats.txt 2>/dev/null | wc -l)

## 各样本比对统计:
EOF
    
    for stats_file in *_alignment_stats.txt; do
        if [ -f "$stats_file" ]; then
            sample=$(basename "$stats_file" _alignment_stats.txt)
            total_reads=$(head -1 "$stats_file" | awk '{print $1}')
            mapped_reads=$(grep "mapped (" "$stats_file" | head -1 | awk '{print $1}')
            mapping_rate=$(grep "mapped (" "$stats_file" | head -1 | awk '{print $5}' | sed 's/[()%]//g')
            
            cat >> alignment_summary.txt << EOF

样本: $sample
总读取数: $total_reads  
映射读取数: $mapped_reads  
映射率: ${mapping_rate}%
EOF
            
        fi
    done
    
    echo "比对统计汇总报告已保存到: alignment_summary.txt" | tee -a $LOG_FILE
}

# ==========================================
# 第四步：变异检测 (GATK HaplotypeCaller)
# ==========================================
echo "[$(date)] 开始变异检测..." | tee -a $LOG_FILE

run_variant_calling() {
    local bam_dir=$1
    local ref_genome=$2
    local output_dir=$3
    
    if [ ! -f "$ref_genome" ]; then
        echo "错误：参考基因组文件不存在: $ref_genome" | tee -a $LOG_FILE
        return 1
    fi
    
    cd $output_dir
    
    # 检查GATK是否已安装  
    if ! command -v gatk &> /dev/null; then  
        echo "警告：未找到GATK，请先安装GATK" | tee -a $LOG_FILE  
        return 1  
    fi  
    
    # 创建输出目录  
    mkdir -p gvcf individual_vcf joint_vcf  
    
    # 1. 对每个样本运行HaplotypeCaller生成gVCF  
    echo "生成个体gVCF文件..." | tee -a $LOG_FILE  
    
    for bam_file in $bam_dir/*.sorted.bam; do  
        if [ -f "$bam_file" ]; then  
            sample_name=$(basename "$bam_file" .sorted.bam)  
            echo "处理样本: $sample_name" | tee -a $LOG_FILE  
            
            gatk --java-options "-Xmx8g" HaplotypeCaller \
                -R "$ref_genome" \
                -I "$bam_file" \
                -O "gvcf/${sample_name}.g.vcf.gz" \
                -ERC GVCF \
                --native-pair-hmm-threads 4 \
                >> $WORK_DIR/logs/hc_${sample_name}.log 2>&1  
            
            if [ $? -eq 0 ]; then  
                echo "样本 $sample_name gVCF生成完成" | tee -a $LOG_FILE  
            else  
                echo "样本 $sample_name gVCF生成失败" | tee -a $LOG_FILE  
            fi  
        fi  
    done  
    
    # 2. 创建样本列表文件用于联合基因型分析  
    echo "创建样本列表文件..." | tee -a $LOG_FILE  
    ls gvcf/*.g.vcf.gz > sample_list.txt  
    
    # 3. 联合基因型分析 (GenotypeGVCFs)  
    echo "运行联合基因型分析..." | tee -a $LOG_FILE  
    
    gatk --java-options "-Xmx16g" GenotypeGVCFs \
        -R "$ref_genome" \
        -V sample_list.txt \
        -O "joint_vcf/all_samples.vcf.gz" \
        >> $WORK_DIR/logs/genotype_gvcfs.log 2>&1  
    
    if [ $? -eq 0 ]; then  
        echo "联合基因型分析完成，生成VCF文件: joint_vcf/all_samples.vcf.gz" | tee -a $LOG_FILE  
        
        # 统计变异位点数量  
        zcat joint_vcf/all_samples.vcf.gz | grep -v "^#" | wc -l > variant_count.txt  
        variant_count=$(cat variant_count.txt)  
        
        cat > variant_calling_summary.txt << EOF  
# 变异检测汇总报告  

分析时间: $(date)  
参考基因组: $(basename "$ref_genome")  
样本数量: $(wc -l < sample_list.txt)  
总变异位点数: $variant_count  

## 处理步骤:  
1. HaplotypeCaller生成gVCF (每个样本)  
2. GenotypeGVCFs联合基因型分析  

## 输出文件:  
- 个体gVCF文件: gvcf/  
- 联合VCF文件: joint_vcf/all_samples.vcf.gz  

EOF
        
        echo "变异检测汇总报告已保存到: variant_calling_summary.txt" | tee -a $LOG_FILE  
        
    else  
        echo "联合基因型分析失败，请检查日志文件" | tee -a $LOG_FILE  
    fi  
}

# ==========================================  
# 第五步：变异过滤和质量控制  
# ==========================================  
echo "[$(date)] 开始变异过滤..." | tee -a $LOG_FILE  

run_variant_filtering() {  
    local vcf_file=$1  
    local ref_genome=$2  
    local output_dir=$3  
    
    if [ ! -f "$vcf_file" ]; then  
        echo "错误：VCF文件不存在: $vcf_file" | tee -a $LOG_FILE  
        return 1  
    fi  
    
    cd $output_dir  
    
    # 检查bcftools是否已安装  
    if ! command -v bcftools &> /dev/null; then  
        echo "警告：未找到bcftools，请先安装bcftools" | tee -a $LOG_FILE  
        return 1  
    fi  
    
    # 1. 基本统计信息  
    echo "生成VCF统计信息..." | tee -a $LOG_FILE  
    bcftools stats "$vcf_file" > vcf_stats.txt  
    
    # 2. 变异过滤 (基于质量分数)  
    echo "运行变异过滤..." | tee -a $LOG_FILE  
    
    # SNP过滤  
    bcftools filter \  
        -i 'TYPE="snp" && QUAL>30 && INFO/DP>10' \  
        "$vcf_file" \  
        > snps_filtered.vcf.gz \  
        2>> $WORK_DIR/logs/filter.log  
    
    # InDel过滤  
    bcftools filter \  
        -i 'TYPE="indel" && QUAL>30 && INFO/DP>10' \  
        "$vcf_file" \  
        > indels_filtered.vcf.gz \  
        2>> $WORK_DIR/logs/filter.log  
    
    # 综合过滤 (保留高质量变异)  
    bcftools filter \  
        -i 'QUAL>30 && INFO/DP>10 && (TYPE="snp" || TYPE="indel")' \  
        "$vcf_file" \  
        > high_quality_variants.vcf.gz \  
        2>> $WORK_DIR/logs/filter.log  
    
    if [ $? -eq 0 ]; then  
        echo "变异过滤完成" | tee -a $LOG_FILE  
        
        # 统计过滤结果  
        snp_count=$(zcat snps_filtered.vcf.gz | grep -v "^#" | wc -l)  
        indel_count=$(zcat indels_filtered.vcf.gz | grep -v "^#" | wc -l)  
        high_qual_count=$(zcat high_quality_variants.vcf.gz | grep -v "^#" | wc -l)  
        
        cat > filtering_summary.txt << EOF  
# 变异过滤汇总报告  

分析时间: $(date)  

## 过滤统计:  
原始变异位点数: $(zcat "$vcf_file" | grep -v "^#" | wc -l)  
高质量SNP数: $snp_count  
高质量InDel数: $indel_count  
总高质量变异数: $high_qual_count  

## 过滤标准:  
- 质量分数 (QUAL) > 30  
- 覆盖度 (DP) > 10  
- SNP和InDel分别处理  

## 输出文件:  
- SNP过滤结果: snps_filtered.vcf.gz  
- InDel过滤结果: indels_filtered.vcf.gz  
- 高质量变异: high_quality_variants.vcf.gz  

EOF
        
        echo "过滤汇总报告已保存到: filtering_summary.txt" | tee -a $LOG_FILE  
        
    else  
        echo "变异过滤失败，请检查日志文件" | tee -a $LOG_FILE  
    fi  
    
}
