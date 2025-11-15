#!/bin/bash

# ==========================================
# 基因组数据批量下载和预处理流程
# 参考：Genome assembly pipeline
# ==========================================

# 设置工作目录
WORK_DIR=$(pwd)
LOG_FILE="$WORK_DIR/download_pipeline.log"
DOWNLOAD_DIR="$WORK_DIR/raw_data"
ASSEMBLY_DIR="$WORK_DIR/assembly"
HIC_DIR="$WORK_DIR/hic_analysis"

# 创建目录结构
mkdir -p $DOWNLOAD_DIR $ASSEMBLY_DIR $HIC_DIR/logs

echo "==========================================" | tee $LOG_FILE
echo "基因组数据批量下载和处理流程开始" | tee -a $LOG_FILE
echo "工作目录: $WORK_DIR" | tee -a $LOG_FILE
echo "==========================================" | tee -a $LOG_FILE

# ==========================================
# 第一步：批量下载原始数据
# ==========================================
echo "[$(date)] 开始批量下载原始数据..." | tee -a $LOG_FILE

cd $DOWNLOAD_DIR

# 检查FTP链接文件,以其中一个为例
if [ ! -f "$WORK_DIR/data_download_links_CNP0003213_ftp.txt" ]; then
    echo "错误：FTP链接文件不存在" | tee -a $LOG_FILE
    exit 1
fi

# 创建下载状态文件
DOWNLOAD_STATUS="$DOWNLOAD_DIR/download_status.txt"
> $DOWNLOAD_STATUS

# 逐行读取FTP链接并下载
TOTAL_LINKS=$(grep -v "^#" "$WORK_DIR/data_download_links_CNP0003213_ftp.txt" | grep -v "^$" | wc -l)
CURRENT=0

while IFS= read -r url || [ -n "$url" ]; do
    # 跳过空行和注释行
    if [[ -n "$url" && ! "$url" =~ ^[[:space:]]*# ]]; then
        CURRENT=$((CURRENT + 1))
        echo "[$(date)] [$CURRENT/$TOTAL_LINKS] 正在下载: $url" | tee -a $LOG_FILE
        
        # 提取文件名
        FILENAME=$(basename "$url")
        
        # 使用wget下载，支持断点续传
        if wget -c --timeout=30 --tries=3 --show-progress "$url" -O "$FILENAME" 2>>$LOG_FILE; then
            FILE_SIZE=$(ls -lh "$FILENAME" | awk '{print $5}')
            echo "SUCCESS: $FILENAME (大小: $FILE_SIZE)" >> $DOWNLOAD_STATUS
            echo "[$(date)] 下载成功: $FILENAME" | tee -a $LOG_FILE
        else
            echo "FAILED: $FILENAME" >> $DOWNLOAD_STATUS
            echo "[$(date)] 下载失败: $FILENAME" | tee -a $LOG_FILE
        fi
    fi
done < "$WORK_DIR/data_download_links_CNP0003213_ftp.txt"

# 统计下载结果
SUCCESS_COUNT=$(grep "SUCCESS" $DOWNLOAD_STATUS | wc -l)
FAILED_COUNT=$(grep "FAILED" $DOWNLOAD_STATUS | wc -l)

echo "==========================================" | tee -a $LOG_FILE
echo "下载完成统计:" | tee -a $LOG_FILE
echo "成功: $SUCCESS_COUNT 个文件" | tee -a $LOG_FILE
echo "失败: $FAILED_COUNT 个文件" | tee -a $LOG_FILE
echo "==========================================" | tee -a $LOG_FILE

if [ $FAILED_COUNT -gt 0 ]; then
    echo "警告：有 $FAILED_COUNT 个文件下载失败，请检查日志" | tee -a $LOG_FILE
    grep "FAILED" $DOWNLOAD_STATUS | tee -a $LOG_FILE
fi

# ==========================================
# 第二步：数据质量控制和过滤 (参考fastp处理)
# ==========================================
echo "[$(date)] 开始数据质量控制..." | tee -a $LOG_FILE

QC_DIR="$DOWNLOAD_DIR/qc_results"
mkdir -p $QC_DIR

# 处理所有FASTQ文件
for fastq_file in *.fastq *.fq *.fastq.gz *.fq.gz; do
    if [ -f "$fastq_file" ]; then
        echo "[$(date)] 处理文件: $fastq_file" | tee -a $LOG_FILE
        
        BASENAME=$(basename "$fastq_file" | sed 's/\.[^.]*$//' | sed 's/\.[^.]*$//')
        
        # 检查是否已安装fastp
        if ! command -v fastp &> /dev/null; then
            echo "警告：未找到fastp，请先安装fastp" | tee -a $LOG_FILE
            continue
        fi
        
        # 使用fastp进行质量控制 (参考参数: --average_qual 15 -l 150)
        fastp \
            --in1 "$fastq_file" \
            --out1 "$QC_DIR/${BASENAME}_clean.fastq.gz" \
            --json "$QC_DIR/${BASENAME}_fastp.json" \
            --html "$QC_DIR/${BASENAME}_fastp.html" \
            --average_qual 15 \
            --length_required 150 \
            --thread 4 \
            >> $QC_DIR/fastp.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo "[$(date)] 质量控制完成: $fastq_file" | tee -a $LOG_FILE
        else
            echo "[$(date)] 质量控制失败: $fastq_file" | tee -a $LOG_FILE
        fi
    fi
done

# ==========================================
# 第三步：基因组组装准备 (Necat流程模拟)
# ==========================================
echo "[$(date)] 开始基因组组装准备..." | tee -a $LOG_FILE

cd $ASSEMBLY_DIR

# 创建Necat配置文件模板
cat > necat_config.txt << EOF
PROJECT=RCAsssembly
ONT_READ_LIST=read_list.txt
GENOME_SIZE=3000000000
THREADS=16
MIN_READ_LENGTH=1000
CNS_OVLP_OPTIONS=cns_ovlp_options.txt
CNS_OPTIONS=cns_options.txt
CNS_OUTPUT_COVERAGE=30
CLEANUP=1
POLISH_CONTIGS=1
EOF

echo "已创建Necat配置文件模板" | tee -a $LOG_FILE

# 创建read列表文件 (需要根据实际数据调整)
echo "# Necat read list (请根据实际文件调整)" > read_list.txt
for clean_fastq in $QC_DIR/*_clean.fastq.gz; do
    if [ -f "$clean_fastq" ]; then
        echo "$clean_fastq" >> read_list.txt
    fi
done

# ==========================================
# 第四步：Hi-C数据分析准备 (Juicer流程)
# ==========================================
echo "[$(date)] 开始Hi-C数据分析准备..." | tee -a $LOG_FILE

cd $HIC_DIR

# 创建Juicer目录结构
mkdir -p fastq references aligned hic

# 创建参考基因组索引 (需要实际的组装基因组)
echo "# Juicer参考基因组准备脚本" > prepare_juicer.sh
cat >> prepare_juicer.sh << 'EOF'
#!/bin/bash
# 假设有组装好的基因组文件
REFERENCE_GENOME="../assembly/contigs.fasta"

if [ -f "$REFERENCE_GENOME" ]; then
    echo "创建BWA索引..."
    bwa index $REFERENCE_GENOME
    
    echo "创建Juicer限制性酶切位点文件..."
    # 假设使用MboI酶切 (GATC)
    python -c "
import sys
with open('$REFERENCE_GENOME', 'r') as f:
    chrom = ''
    seq = ''
    for line in f:
        if line.startswith('>'):
            if chrom and seq:
                print(f'{chrom}\t{len(seq)}')
            chrom = line.strip()[1:].split()[0]
            seq = ''
        else:
            seq += line.strip()
    if chrom and seq:
        print(f'{chrom}\t{len(seq)}')
" > references/chrom_sizes.txt
    
    echo "Juicer准备完成"
else
    echo "警告：未找到参考基因组文件，请先完成基因组组装"
fi
EOF

chmod +x prepare_juicer.sh

# ==========================================
# 第五步：3D-DNA分析准备
# ==========================================
echo "[$(date)] 开始3D-DNA分析准备..." | tee -a $LOG_FILE

# 创建3D-DNA分析脚本模板
cat > run_3d_dna.sh << 'EOF'
#!/bin/bash
# 3D-DNA染色体构建脚本

INPUT_ASM="../assembly/contigs.fasta"
INPUT_HIC="aligned/aligned_reads.bam"

if [ ! -f "$INPUT_ASM" ]; then
    echo "错误：未找到输入组装文件"
    exit 1
fi

if [ ! -f "$INPUT_HIC" ]; then
    echo "错误：未找到Hi-C比对文件"
    exit 1
fi

echo "运行3D-DNA染色体构建..."
# run-asm-pipeline.sh $INPUT_ASM $INPUT_HIC

echo "可视化结果..."
# run-asm-visualizer.sh output.final.assembly output.hic

echo "3D-DNA分析完成"
EOF

chmod +x run_3d_dna.sh

# ==========================================
# 第六步：生成处理报告和热图可视化准备
# ==========================================
echo "[$(date)] 生成处理报告..." | tee -a $LOG_FILE

REPORT_DIR="$WORK_DIR/reports"
mkdir -p $REPORT_DIR

# 创建处理总结报告
cat > $REPORT_DIR/pipeline_summary.txt << EOF
==========================================
基因组组装和Hi-C分析处理报告
生成时间: $(date)
==========================================

1. 数据下载统计:
   - 总链接数: $TOTAL_LINKS
   - 成功下载: $SUCCESS_COUNT
   - 下载失败: $FAILED_COUNT

2. 数据质量控制:
   - 处理目录: $QC_DIR
   - 质控工具: fastp v0.23.2
   - 质控参数: --average_qual 15 -l 150

3. 基因组组装:
   - 组装工具: Necat v0.01 + Pilon v1.22
   - 预期流程:
     * Nanopore读取校正 (Necat)
     * 初始contig组装 (Necat)
     * contig校正 (Necat)
     * indel和SNP错误校正 (Pilon)

4. Hi-C分析:
   - 分析工具: Juicer v1.6 + 3D-DNA v170123
   - 预期结果:
     * 877.40 million reads (263.2Gb clean data)
     * 65.79% read pairs唯一比对到基因组
     * 47.96%有效交互对用于Hi-C组装

5. 后续步骤:
   - 运行基因组组装: cd $ASSEMBLY_DIR && necat.pl config.txt
   - 运行Hi-C分析: cd $HIC_DIR && ./prepare_juicer.sh && juicer.sh
   - 运行3D-DNA: cd $HIC_DIR && ./run_3d_dna.sh

==========================================
EOF

echo "处理报告已生成: $REPORT_DIR/pipeline_summary.txt" | tee -a $LOG_FILE

# ==========================================
# 第七步：创建运行脚本和说明文档
# ==========================================

# 创建主运行脚本
cat > run_pipeline.sh << 'EOF'
#!/bin/bash

echo "选择要运行的步骤:"
echo "1. 数据质量控制 (fastp)"
echo "2. 基因组组装 (Necat)"
echo "3. Hi-C数据处理 (Juicer)"
echo "4. 染色体构建 (3D-DNA)"
echo "5. 运行所有步骤"
echo "请输入选项 (1-5): "

read choice

case $choice in
    1)
        echo "运行数据质量控制..."
        cd qc_analysis && ./run_fastp.sh
        ;;
    2)
        echo "运行基因组组装..."
        cd assembly && necat.pl config.txt
        ;;
    3)
        echo "运行Hi-C数据处理..."
        cd hic_analysis && ./prepare_juicer.sh && juicer.sh
        ;;
    4)
        echo "运行染色体构建..."
        cd hic_analysis && ./run_3d_dna.sh
        ;;
    5)
        echo "运行所有步骤..."
        cd qc_analysis && ./run_fastp.sh
        cd ../assembly && necat.pl config.txt  
        cd ../hic_analysis && ./prepare_juicer.sh && juicer.sh && ./run_3d_dna.sh
        ;;
    *)
        echo "无效选项"
        ;;
esac
EOF

chmod +x run_pipeline.sh

echo "==========================================" | tee -a $LOG_FILE
echo "流程设置完成！" | tee -a $LOG_FILE
echo "" | tee -a $LOG_FILE
echo "主要目录结构:" | tee -a $LOG_FILE
echo "- 原始数据: $DOWNLOAD_DIR" | tee -a $LOG_FILE  
echo "- 质量控制: $QC_DIR" | tee -a $LOG_FILE
echo "- 基因组组装: $ASSEMBLY_DIR" | tee -a $LOG_FILE
echo "- Hi-C分析: $HIC_DIR" | tee -a $LOG_FILE
echo "- 处理报告: $REPORT_DIR" | tee -a $LOG_FILE
echo "" | tee -a $LOG_FILE
echo "运行步骤:" | tee -a $LOG_FILE
echo "1. 首先确保所有数据已成功下载" | tee -a $LOG_FILE
echo "2. 运行质量控制: cd qc_analysis && ./run_fastp.sh" | tee -a $LOG_FILE  
echo "3. 运行基因组组装: cd assembly && necat.pl config.txt" | tee -a $LOG_FILE
echo "4. 运行Hi-C分析: cd hic_analysis && ./prepare_juicer.sh && juicer.sh" | tee -a $LOG_FILE
echo "" | tee -a $LOG_FILE  
echo "查看详细日志: cat $LOG_FILE" | tee -a $LOG_FILE
echo "查看下载状态: cat $DOWNLOAD_STATUS" | tee -a $LOG_FILE  
echo "查看处理报告: cat $REPORT_DIR/pipeline_summary.txt" | tee -a $LOG_FILE
echo "==========================================" | tee -a $LOG_FILE

echo ""
echo "流程脚本已创建完成！"
echo "请运行以下命令开始处理："
echo "bash run_pipeline.sh"
