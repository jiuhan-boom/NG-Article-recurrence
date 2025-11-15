# ==========================================
# RNA-seq和eQTL分析完整流程
# RNA-seq and eQTL analysis pipeline
# ==========================================

# 设置工作目录和环境变量
WORK_DIR=$(pwd)
RNA_DIR="$WORK_DIR/rna_seq_analysis"
EQTL_DIR="$WORK_DIR/eqtl_analysis"
LOG_FILE="$WORK_DIR/rna_eqtl_analysis.log"

# 创建目录结构
mkdir -p $RNA_DIR $EQTL_DIR
mkdir -p $RNA_DIR/raw_data $RNA_DIR/quality_control $RNA_DIR/alignment $RNA_DIR/quantification
mkdir -p $EQTL_DIR/expression_data $EQTL_DIR/genotype_data $EQTL_DIR/eqtl_mapping $EQTL_DIR/results
mkdir -p $WORK_DIR/logs $WORK_DIR/scripts

echo "==========================================" | tee $LOG_FILE
echo "RNA-seq和eQTL分析流程开始" | tee -a $LOG_FILE
echo "工作目录: $WORK_DIR" | tee -a $LOG_FILE
echo "==========================================" | tee -a $LOG_FILE

# ==========================================
# 第一步：RNA-seq数据质量控制和预处理
# ==========================================
echo "[$(date)] 开始RNA-seq数据质量控制..." | tee -a $LOG_FILE

run_rna_quality_control() {
    local input_dir=$1
    local output_dir=$2
    
    cd $output_dir
    
    # 检查Trimmomatic是否已安装
    if ! command -v trimmomatic &> /dev/null; then
        echo "警告：未找到Trimmomatic，请先安装Trimmomatic" | tee -a $LOG_FILE
        return 1
    fi
    
    # 创建质量控制报告
    cat > qc_summary.txt << EOF
# RNA-seq数据质量控制报告

分析时间: $(date)
输入目录: $input_dir
输出目录: $output_dir

## 处理的样本列表:
EOF
    
    # 处理配对末端FASTQ文件
    for fastq1 in $input_dir/*_1.fastq.gz $input_dir/*_R1.fastq.gz; do
        if [ -f "$fastq1" ]; then
            # 获取对应的配对文件
            fastq2=${fastq1/_1.fastq.gz/_2.fastq.gz}
            fastq2=${fastq2/_R1.fastq.gz/_R2.fastq.gz}
            
            if [ -f "$fastq2" ]; then
                sample_name=$(basename "$fastq1" | sed 's/_1\.fastq\.gz$//' | sed 's/_R1\.fastq\.gz$//')
                
                echo "处理样本: $sample_name" | tee -a $LOG_FILE
                echo "- $sample_name" >> qc_summary.txt
                
                # 运行Trimmomatic质量控制
                trimmomatic PE \
                    -threads 8 \
                    -phred33 \
                    "$fastq1" "$fastq2" \
                    "${sample_name}_clean_1.fastq.gz" \
                    "${sample_name}_unpaired_1.fastq.gz" \
                    "${sample_name}_clean_2.fastq.gz" \
                    "${sample_name}_unpaired_2.fastq.gz" \
                    ILLUMINACLIP:/path/to/adapters/TruSeq3-PE.fa:2:30:10 \
                    LEADING:3 \
                    TRAILING:3 \
                    SLIDINGWINDOW:4:15 \
                    MINLEN:36 \
                    >> $WORK_DIR/logs/trimmomatic_${sample_name}.log 2>&1
                
                if [ $? -eq 0 ]; then
                    echo "样本 $sample_name 质量控制完成" | tee -a $LOG_FILE
                    
                    # 统计清理后的读取数
                    clean_reads_1=$(zcat "${sample_name}_clean_1.fastq.gz" | wc -l)
                    clean_reads_2=$(zcat "${sample_name}_clean_2.fastq.gz" | wc -l)
                    
                    echo "  清理后读取数 (R1): $((clean_reads_1/4))" >> qc_summary.txt
                    echo "  清理后读取数 (R2): $((clean_reads_2/4))" >> qc_summary.txt
                    
                else
                    echo "样本 $sample_name 质量控制失败" | tee -a $LOG_FILE
                fi
            fi
        fi
    done
    
    echo "RNA-seq质量控制完成，报告保存到: qc_summary.txt" | tee -a $LOG_FILE
}

# ==========================================
# 第二步：序列比对到参考基因组 (Hisat2)
# ==========================================
echo "[$(date)] 开始序列比对..." | tee -a $LOG_FILE

run_rna_alignment() {
    local clean_data_dir=$1
    local ref_genome=$2
    local output_dir=$3
    
    if [ ! -f "$ref_genome" ]; then
        echo "错误：参考基因组文件不存在: $ref_genome" | tee -a $LOG_FILE
        return 1
    fi
    
    cd $output_dir
    
    # 检查Hisat2是否已安装
    if ! command -v hisat2 &> /dev/null; then
        echo "警告：未找到Hisat2，请先安装Hisat2 v2.1.0" | tee -a $LOG_FILE
        
        
        chmod +x install_hisat2.sh
        echo "已创建Hisat2安装脚本: install_hisat2.sh" | tee -a $LOG_FILE
        return 1
    fi
    
    # 检查samtools是否已安装
    if ! command -v samtools &> /dev/null; then
        echo "警告：未找到samtools，请先安装samtools" | tee -a $LOG_FILE
        return 1
    fi
    
    # 构建Hisat2索引 (如果不存在)
    if [ ! -f "${ref_genome}.1.ht2" ]; then
        echo "构建Hisat2索引..." | tee -a $LOG_FILE
        hisat2-build "$ref_genome" "$ref_genome" > hisat2_build.log 2>&1 &
        wait
    fi
    
    # 处理所有清理后的数据进行比对
    echo "开始序列比对..." | tee -a $LOG_FILE
    
    for fastq1 in $clean_data_dir/*_clean_1.fastq.gz; do
        if [ -f "$fastq1" ]; then
            fastq2=${fastq1/_clean_1.fastq.gz/_clean_2.fastq.gz}
            sample_name=$(basename "$fastq1" | sed 's/_clean_1\.fastq\.gz$//')
            
            if [ -f "$fastq2" ]; then
                echo "比对样本: $sample_name" | tee -a $LOG_FILE
                
                # Hisat2比对 (使用stranded参数)
                hisat2 \
                    -p 16 \
                    --dta \
                    --rna-strandness RF \
                    -x "$ref_genome" \
                    -1 "$fastq1" \
                    -2 "$fastq2" \
                    2>> $WORK_DIR/logs/hisat2_${sample_name}.log | \
                samtools view -bS - | \
                samtools sort -@ 8 -o "${sample_name}.sorted.bam" -
                
                if [ $? -eq 0 ]; then
                    echo "样本 $sample_name 比对完成" | tee -a $LOG_FILE
                    
                    # 索引BAM文件
                    samtools index "${sample_name}.sorted.bam"
                    
                    # 统计比对结果
                    samtools flagstat "${sample_name}.sorted.bam" > "${sample_name}_alignment_stats.txt"
                    
                    # 过滤高比对质量的读取 (MAPQ > 30)
                    samtools view -b -q 30 "${sample_name}.sorted.bam" > "${sample_name}_hq.sorted.bam"
                    samtools index "${sample_name}_hq.sorted.bam"
                    
                    echo "样本 $sample_name 高质量比对完成" | tee -a $LOG_FILE
                    
                else
                    echo "样本 $sample_name 比对失败" | tee -a $LOG_FILE
                fi
            fi
        fi
    done
    
    echo "序列比对完成" | tee -a $LOG_FILE
    
    # 创建比对统计汇总报告
    cat > alignment_summary.txt << EOF
# RNA-seq序列比对统计汇总

分析时间: $(date)
参考基因组: $(basename "$ref_genome")
处理样本数: $(ls *_alignment_stats.txt 2>/dev/null | wc -l)

## 各样本比对统计:
EOF
    
    total_reads=0
    mapped_reads=0
    
    for stats_file in *_alignment_stats.txt; do
        if [ -f "$stats_file" ]; then
            sample=$(basename "$stats_file" _alignment_stats.txt)
            sample_total=$(head -1 "$stats_file" | awk '{print $1}')
            sample_mapped=$(grep "mapped (" "$stats_file" | head -1 | awk '{print $1}')
            mapping_rate=$(grep "mapped (" "$stats_file" | head -1 | awk '{print $5}' | sed 's/[()%]//g')
            
            echo "" >> alignment_summary.txt
            echo "样本: $sample" >> alignment_summary.txt
            echo "总读取数: $sample_total" >> alignment_summary.txt
            echo "映射读取数: $sample_mapped" >> alignment_summary.txt
            echo "映射率: ${mapping_rate}%" >> alignment_summary.txt
            
            total_reads=$((total_reads + sample_total))
            mapped_reads=$((mapped_reads + sample_mapped))
        fi
    done
    
    avg_reads=$((total_reads / $(ls *_alignment_stats.txt 2>/dev/null | wc -l)))
    avg_mapped=$((mapped_reads / $(ls *_alignment_stats.txt 2>/dev/null | wc -l)))
    
    echo "" >> alignment_summary.txt
    echo "平均统计:" >> alignment_summary.txt
    echo "平均总读取数: $avg_reads" >> alignment_summary.txt
    echo "平均映射读取数: $avg_mapped" >> alignment_summary.txt
    
    echo "比对统计汇总报告已保存到: alignment_summary.txt" | tee -a $LOG_FILE
    
}

# ==========================================
# 第三步：转录本组装和表达量估计 (StringTie)
# ==========================================
echo "[$(date)] 开始转录本组装和表达量估计..." | tee -a $LOG_FILE

run_transcript_assembly() {
    local bam_dir=$1
    local ref_genome=$2
    local ref_annotation=$3  # GTF注释文件 (可选)
    local output_dir=$4
    
    cd $output_dir
    
    # 检查StringTie是否已安装
    if ! command -v stringtie &> /dev/null; then
        echo "警告：未找到StringTie，请先安装StringTie v2.20" | tee -a $LOG_FILE
        
        chmod +x install_stringtie.sh
        echo "已创建StringTie安装脚本: install_stringtie.sh" | tee -a $LOG_FILE
        return 1
    fi
    
    # 创建目录存储每个样本的结果
    mkdir -p individual_assemblies merged_assembly
    
    # 1. 对每个样本进行转录本组装和表达量估计
    echo "进行个体样本转录本组装..." | tee -a $LOG_FILE
    
    for bam_file in $bam_dir/*_hq.sorted.bam; do
        if [ -f "$bam_file" ]; then
            sample_name=$(basename "$bam_file" _hq.sorted.bam)
            echo "处理样本: $sample_name" | tee -a $LOG_FILE
            
            # StringTie组装和表达量估计
            if [ -f "$ref_annotation" ]; then
                # 使用参考注释进行引导组装
                stringtie \
                    "$bam_file" \
                    -p 8 \
                    -G "$ref_annotation" \
                    -o "individual_assemblies/${sample_name}.gtf" \
                    -A "individual_assemblies/${sample_name}_gene_abund.tab" \
                    -B \
                    >> $WORK_DIR/logs/stringtie_${sample_name}.log 2>&1
            else
                # 从头组装
                stringtie \
                    "$bam_file" \
                    -p 8 \
                    -o "individual_assemblies/${sample_name}.gtf" \
                    -A "individual_assemblies/${sample_name}_gene_abund.tab" \
                    >> $WORK_DIR/logs/stringtie_${sample_name}.log 2>&1
            fi
            
            if [ $? -eq 0 ]; then
                echo "样本 $sample_name 转录本组装完成" | tee -a $LOG_FILE
            else
                echo "样本 $sample_name 转录本组装失败" | tee -a $LOG_FILE
            fi
        fi
    done
    
    # 2. 创建样本列表文件用于合并组装 (可选)
    echo "创建样本列表文件..." | tee -a $LOG_FILE
    ls individual_assemblies/*.gtf > sample_list.txt
    
    # 3. 合并所有样本的转录本 (StringTie --merge)
    echo "合并所有样本转录本..." | tee -a $LOG_FILE
    
    if [ -f "$ref_annotation" ]; then
        stringtie \
            --merge \
            -p 8 \
            -G "$ref_annotation" \
            -o merged_assembly/merged_transcripts.gtf \
            sample_list.txt \
            >> $WORK_DIR/logs/stringtie_merge.log 2>&1
    else
        stringtie \
            --merge \
            -p 8 \
            -o merged_assembly/merged_transcripts.gtf \
            sample_list.txt \
            >> $WORK_DIR/logs/stringtie_merge.log 2>&1
    fi
    
    if [ $? -eq 0 ]; then
        echo "转录本合并完成" | tee -a $LOG_FILE
        
        # 4. 对每个样本重新估算相对于合并转录本的表达量
        echo "重新估算表达量..." | tee -a $LOG_FILE
        
        for bam_file in $bam_dir/*_hq.sorted.bam; do
            if [ -f "$bam_file" ]; then
                sample_name=$(basename "$bam_file" _hq.sorted.bam)
                echo "重新估算样本: $sample_name 表达量" | tee -a $LOG_FILE
                
                stringtie \
                    "$bam_file" \
                    -e \
                    -B \
                    -p 8 \
                    -G merged_assembly/merged_transcripts.gtf \
                    -o "merged_assembly/${sample_name}_merged.gtf" \
                    -A "merged_assembly/${sample_name}_merged_gene_abund.tab" \
                    >> $WORK_DIR/logs/stringtie_reestimate_${sample_name}.log 2>&1
                
                if [ $? -eq 0 ]; then
                    echo "样本 $sample_name 表达量重新估算完成" | tee -a $LOG_FILE
                else
                    echo "样本 $sample_name 表达量重新估算失败" | tee -a $LOG_FILE
                fi
                
            fi
        done
        
        echo "转录本组装和表达量估计完成" | tee -a $LOG_FILE
        
        # 创建组装统计报告
        cat > assembly_summary.txt << EOF
# 转录本组装统计摘要

分析时间: $(date)
参考基因组: $(basename "$ref_genome")
样本数量: $(ls individual_assemblies/*.gtf 2>/dev/null | wc -l)

## 组装结果:
- 合并转录本文件: merged_assembly/merged_transcripts.gtf  
- 各样本表达量文件: merged_assembly/*_merged_gene_abund.tab

## 统计信息:
$(stringtie --version)

EOF
        
        echo "转录本组装摘要报告已保存到: assembly_summary.txt" | tee -a $LOG_FILE
        
    else  
        echo "转录本合并失败，请检查日志文件" | tee -a $LOG_FILE  
    fi  
    
}

# ==========================================  
# 第四步：表达量提取和FPKM计算 (Ballgown)  
# ==========================================  
echo "[$(date)] 开始表达量提取和FPKM计算..." | tee -a $LOG_FILE  

extract_expression_values() {  
    local assembly_dir=$1  
    local output_dir=$2  
    
    cd $output_dir  
    
    # 检查R和Ballgown是否已安装  
    if ! R --version &> /dev/null; then  
        echo "警告：未找到R，请先安装R" | tee -a $LOG_FILE  
        return 1  
    fi   
}