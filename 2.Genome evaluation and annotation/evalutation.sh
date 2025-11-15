# ==========================================
# 基因组评估和注释完整流程
# Genome evaluation and annotation pipeline
# ==========================================

# 设置工作目录和环境变量
WORK_DIR=$(pwd)
EVAL_DIR="$WORK_DIR/genome_evaluation"
ANNOTATION_DIR="$WORK_DIR/gene_annotation"
LOG_FILE="$WORK_DIR/genome_analysis.log"

# 创建目录结构
mkdir -p $EVAL_DIR $ANNOTATION_DIR
mkdir -p $EVAL_DIR/busco $EVAL_DIR/repeats $EVAL_DIR/kmer_analysis $EVAL_DIR/ltr_analysis
mkdir -p $ANNOTATION_DIR/homolog $ANNOTATION_DIR/denovo $ANNOTATION_DIR/maker
mkdir -p $WORK_DIR/logs

echo "==========================================" | tee $LOG_FILE
echo "基因组评估和注释流程开始" | tee -a $LOG_FILE
echo "工作目录: $WORK_DIR" | tee -a $LOG_FILE
echo "==========================================" | tee -a $LOG_FILE

# ==========================================
# 第一步：基因组完整性评估 (BUSCO)
# ==========================================
echo "[$(date)] 开始BUSCO基因组完整性评估..." | tee -a $LOG_FILE

run_busco_analysis() {
    local genome_file=$1
    local output_dir=$2
    
    if [ ! -f "$genome_file" ]; then
        echo "错误：基因组文件不存在: $genome_file" | tee -a $LOG_FILE
        return 1
    fi

    # 下载embryophyta_odb10数据库（如果不存在）
    if [ ! -d "embryophyta_odb10" ]; then
        echo "下载embryophyta_odb10数据库..." | tee -a $LOG_FILE
        busco --download embryophyta_odb10 >> $WORK_DIR/logs/busco_download.log 2>&1
    fi
    
    # 运行BUSCO分析
    echo "运行BUSCO分析..." | tee -a $LOG_FILE
    busco \
        --in "$genome_file" \
        --out busco_embryophyta \
        --lineage_dataset embryophyta_odb10 \
        --mode genome \
        --cpu 16 \
        --force \
        >> $WORK_DIR/logs/busco_analysis.log 2>&1
    
    if [ $? -eq 0 ]; then
        echo "[$(date)] BUSCO分析完成" | tee -a $LOG_FILE
        echo "结果目录: $output_dir/run_busco_embryophyta" | tee -a $LOG_FILE
        
        # 生成摘要报告
        if [ -f "run_busco_embryophyta/short_summary.specific.embryophyta_odb10.busco_embryophyta.txt" ]; then
            cat "run_busco_embryophyta/short_summary.specific.embryophyta_odb10.busco_embryophyta.txt" > busco_summary.txt
            echo "BUSCO摘要已保存到: busco_summary.txt" | tee -a $LOG_FILE
        fi
    else
        echo "[$(date)] BUSCO分析失败" | tee -a $LOG_FILE
    fi
}

# ==========================================
# 第二步：重复序列注释
# ==========================================
echo "[$(date)] 开始重复序列注释..." | tee -a $LOG_FILE

run_repeat_annotation() {
    local genome_file=$1
    local output_dir=$2
    
    if [ ! -f "$genome_file" ]; then
        echo "错误：基因组文件不存在: $genome_file" | tee -a $LOG_FILE
        return 1
    fi
    
    cd $output_dir
    
    # 创建重复序列分析报告
    echo "# 重复序列注释分析报告" > repeat_annotation_report.txt
    echo "分析时间: $(date)" >> repeat_annotation_report.txt
    echo "基因组文件: $genome_file" >> repeat_annotation_report.txt
    echo "" >> repeat_annotation_report.txt
    
    # 1. de novo预测 (RepeatModeler)
    echo "1. 运行RepeatModeler de novo预测..." | tee -a $LOG_FILE
    
    if command -v BuildDatabase &> /dev/null && command -v RepeatModeler &> /dev/null; then
        # 构建数据库
        BuildDatabase -name genome_db -engine ncbi "$genome_file" >> $WORK_DIR/logs/repeatmodeler.log 2>&1
        
        # 运行RepeatModeler
        RepeatModeler \
            -database genome_db \
            -engine ncbi \
            -pa 16 \
            >> $WORK_DIR/logs/repeatmodeler.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo "RepeatModeler分析完成" >> repeat_annotation_report.txt
            echo "结果文件: genome_db-families.fa" >> repeat_annotation_report.txt
        else
            echo "RepeatModeler分析失败" >> repeat_annotation_report.txt
        fi
    else
        echo "警告：未找到RepeatModeler" >> repeat_annotation_report.txt
    fi
    
    # 2. LTR预测 (LTR_FINDER)
    echo "2. 运行LTR_FINDER预测..." | tee -a $LOG_FILE
    
    if command -v LTR_FINDER_parallel &> /dev/null; then
        LTR_FINDER_parallel \
            -seq "$genome_file" \
            -threads 16 \
            -harvest_out \
            > ltr_finder_output.scn 2>> $WORK_DIR/logs/ltr_finder.log
        
        if [ $? -eq 0 ]; then
            echo "LTR_FINDER分析完成" >> repeat_annotation_report.txt
            echo "结果文件: ltr_finder_output.scn" >> repeat_annotation_report.txt
        else
            echo "LTR_FINDER分析失败" >> repeat_annotation_report.txt
        fi
    else
        echo "警告：未找到LTR_FINDER" >> repeat_annotation_report.txt
    fi
    
    # 3. 同源搜索 (RepeatMasker)
    echo "3. 运行RepeatMasker同源搜索..." | tee -a $LOG_FILE
    
    if command -v RepeatMasker &> /dev/null; then
        # 检查RepBase数据库是否存在
        if [ ! -d "/usr/local/lib/RepeatMasker/Libraries" ]; then
            echo "警告：未找到RepBase数据库" >> repeat_annotation_report.txt
        fi
        
        RepeatMasker \
            -pa 16 \
            -gff \
            -lib genome_db-families.fa \
            "$genome_file" \
            >> $WORK_DIR/logs/repeatmasker.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo "RepeatMasker分析完成" >> repeat_annotation_report.txt
            echo "结果文件: $(basename ${genome_file}.out)" >> repeat_annotation_report.txt
            echo "GFF文件: $(basename ${genome_file}.out.gff)" >> repeat_annotation_report.txt
        else
            echo "RepeatMasker分析失败" >> repeat_annotation_report.txt
        fi
    else
        echo "警告：未找到RepeatMasker" >> repeat_annotation_report.txt
    fi
    
    # 4. 蛋白质重复序列掩蔽 (RepeatProteinMask)
    echo "4. 运行RepeatProteinMask..." | tee -a $LOG_FILE
    
    if command -v RepeatProteinMask &> /dev/null; then
        RepeatProteinMask \
            -engine ncbi \
            "$genome_file" \
            >> $WORK_DIR/logs/repeatproteinmask.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo "RepeatProteinMask分析完成" >> repeat_annotation_report.txt
        else
            echo "RepeatProteinMask分析失败" >> repeat_annotation_report.txt
        fi
    else
        echo "警告：未找到RepeatProteinMask" >> repeat_annotation_report.txt
    fi
    
    echo "重复序列注释完成，报告保存到: $output_dir/repeat_annotation_report.txt" | tee -a $LOG_FILE
}

# ==========================================
# 第三步：k-mer分析 (Meryl + Merqury)
# ==========================================
echo "[$(date)] 开始k-mer分析..." | tee -a $LOG_FILE

run_kmer_analysis() {
    local reads_files=$1  # 读取文件列表，用逗号分隔
    local genome_file=$2
    local output_dir=$3
    
    if [ ! -f "$genome_file" ]; then
        echo "错误：基因组文件不存在: $genome_file" | tee -a $LOG_FILE
        return 1
    fi
    
    cd $output_dir
    
    # 检查Meryl是否已安装
    if ! command -v meryl &> /dev/null; then
        echo "警告：未找到Meryl" | tee -a $LOG_FILE
        return 1
    fi
    
    # 检查Merqury是否已安装  
    if ! command -v merqury.sh &> /dev/null; then
        echo "警告：未找到Merqury" | tee -a $LOG_FILE
        return 1
    fi
    
    # 构建k-mer数据库 (k=20)
    echo "构建k-mer数据库 (k=20)..." | tee -a $LOG_FILE
    
    # 处理读取文件列表
    IFS=',' read -ra READ_FILES <<< "$reads_files"
    
    # 合并所有读取文件的k-mer计数
    meryl_commands=""
    for i in "${!READ_FILES[@]}"; do
        read_file="${READ_FILES[$i]}"
        if [ -f "$read_file" ]; then
            meryl count k=20 output read_${i}.meryl "$read_file" >> $WORK_DIR/logs/meryl.log 2>&1 &
            meryl_commands="$meryl_commands read_${i}.meryl"
        fi
    done
    
    # 等待所有后台进程完成
    wait
    
    # 合并k-mer数据库
    if [ -n "$meryl_commands" ]; then
        meryl union-sum output genome_reads.meryl $meryl_commands >> $WORK_DIR/logs/meryl.log 2>&1
        
        # 运行Merqury进行基因组质量评估
        echo "运行Merqury进行基因组质量评估..." | tee -a $LOG_FILE
        merqury.sh genome_reads.meryl "$genome_file" merqury_analysis >> $WORK_DIR/logs/merqury.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo "[$(date)] k-mer分析完成" | tee -a $LOG_FILE
            
            # 生成质量报告摘要
            cat > kmer_quality_summary.txt << EOF
# k-mer质量分析摘要报告

分析时间: $(date)
基因组文件: $genome_file

## Merqury分析结果:
$(cat merqury_analysis*.qv 2>/dev/null || echo "暂无QV结果")

## k-mer统计:
$(cat merqury_analysis*.spectra-cn.fl.summary 2>/dev/null || echo "暂无统计结果")

## 基因组完整性评估:
- Consensus k-mer QV score: 查看merqury_analysis*.qv文件获取详细信息
- k-mer大小: 20

EOF
            
            echo "k-mer质量报告已保存到: kmer_quality_summary.txt" | tee -a $LOG_FILE
            
        else
            echo "[$(date)] Merqury分析失败" | tee -a $LOG_FILE
        fi
        
    else
        echo "[$(date)] 未找到有效的读取文件进行k-mer分析" | tee -a $LOG_FILE
    fi
}

# ==========================================
# 第四步：LTR分析 (EDTA + LTR_retriever)
# ==========================================
echo "[$(date)] 开始LTR分析..." | tee -a $LOG_FILE

run_ltr_analysis() {
    local genome_file=$1
    local output_dir=$2
    
    if [ ! -f "$genome_file" ]; then
        echo "错误：基因组文件不存在: $genome_file" | tee -a $LOG_FILE
        return 1
    fi
    
    cd $output_dir
    
    # 检查EDTA是否已安装  
    if ! command -v EDTA.pl &> /dev/null; then
        echo "警告：未找到EDTA，请先安装EDTA v2.2" | tee -a $LOG_FILE
        
        # 创建替代的LTR_retriever分析脚本（如果EDTA不可用）
        cat > run_ltr_retriever_only.sh << 'EOF'
#!/bin/bash

GENOME_FILE=$1

if [ ! -f "$GENOME_FILE" ]; then
    echo "错误：基因组文件不存在"
    exit 1  
fi

# 检查LTR_retriever是否已安装  
if ! command -v LTR_retriever &> /dev/null; then
    echo "错误：未找到LTR_retriever，请先安装LTR_retriever v3.0.1"
    exit 1  
fi

echo "运行LTR_retriever..."
LTR_retriever \
    -genome "$GENOME_FILE" \
    -inhouse "$GENOME_FILE" \
    -threads 16 \
    > ltr_retriever.log 2>&1

if [ $? -eq 0 ]; then  
    echo "LTR_retriever分析完成"
    
    # 计算LAI (LTR Assembly Index)
    if [ -f "${GENOME_FILE}.mod.out" ]; then  
        LAI=$(LTR_retriever -lai "${GENOME_FILE}.mod.out")
        echo "LAI值: $LAI"
        echo "LAI值: $LAI" > lai_score.txt  
    fi
    
else  
    echo "LTR_retriever分析失败"
fi

EOF
        
        chmod +x run_ltr_retriever_only.sh  
        
        # 运行替代脚本  
        ./run_ltr_retriever_only.sh "$genome_file"
        
        return 0  
    fi

    # 运行EDTA完整分析  
    echo "运行EDTA LTR预测..." | tee -a $LOG_FILE  
    EDTA.pl \
        --genome "$genome_file" \
        --species others \
        --step all \
        --threads 16 \
        > edta_analysis.log 2>&1

    if [ $? -eq 0 ]; then  
        echo "[$(date)] EDTA LTR分析完成" | tee -a $LOG_FILE
        
        # 运行LTR_retriever计算LAI  
        if command -v LTR_retriever &> /dev/null; then  
            echo "计算LTR Assembly Index (LAI)..." | tee -a $LOG_FILE
            
            # 使用EDTA输出作为输入  
            edta_output="${genome_file}.EDTA.raw/LTR/${genome_file%.fasta}.ltr.fa"
            
            if [ -f "$edta_output" ]; then  
                LTR_retriever \
                    -genome "$genome_file" \
                    -inhouse "$edta_output" \
                    -threads 16 \
                    > ltr_retriever.log 2>&1
                
                # 计算LAI  
                if [ -f "${genome_file}.mod.out" ]; then  
                    LAI=$(LTR_retriever -lai "${genome_file}.mod.out")
                    echo "LAI值: $LAI"
                    echo "LAI值: $LAI" > lai_score.txt
                    
                    # 创建LTR分析报告  
                    cat > ltr_analysis_report.txt << EOF  
# LTR分析报告  

分析时间: $(date)  
基因组文件: $genome_file  

## EDTA分析结果:  
- LTR预测完成  
- 结果目录: ${genome_file}.EDTA.raw  

## LTR Assembly Index (LAI):  
- LAI值: $LAI  

## 分析说明:  
LAI值用于评估组装的重复序列连续性：  
- LAI > 10: 高质量组装  
- LAI 5-10: 中等质量组装  
- LAI < 5: 需要改进的组装  

EOF
                    
                    echo "LTR分析报告已保存到: ltr_analysis_report.txt" | tee -a $LOG_FILE
                    
                fi
                
            else  
                echo "警告：未找到EDTA输出文件用于LAI计算" | tee -a $LOG_FILE  
            fi
            
        else  
            echo "警告：未找到LTR_retriever，跳过LAI计算" | tee -a $LOG_FILE  
        fi
        
    else  
        echo "[$(date)] EDTA LTR分析失败" | tee -a $LOG_FILE  
    fi  
}

# ==========================================  
# 第五步：同源分析和基因注释准备  
# ==========================================  
echo "[$(date)] 开始同源分析和基因注释准备..." | tee -a $LOG_FILE  

run_homolog_analysis() {  
    local genome_file=$1  
    local output_dir=$2  
    
    if [ ! -f "$genome_file" ]; then  
        echo "错误：基因组文件不存在: $genome_file" | tee -a $LOG_FILE  
        return 1  
    fi  
    
    cd $output_dir  
    
    # 创建同源分析目录结构  
    mkdir -p {setaria,oryza,brachypodium,sorghum,arabidopsis,maize_b73,maize_mo17,maize_sk}  
    
    # 创建同源分析配置文件  
    cat > homolog_analysis_config.txt << EOF  
# 同源分析配置文件  

## 模型物种列表:  
1. Setaria italica (谷子)  
2. Oryza sativa (水稻)  
3. Brachypodium distachyon (短柄草)  
4. Sorghum bicolor (高粱)  
5. Arabidopsis thaliana (拟南芥)  

## 玉米基因组:  
6. B73 (玉米B73参考基因组)  
7. Mo17 (玉米Mo17基因组)  
8. SK (玉米SK基因组)  

## 分析流程:  
- 蛋白质序列比对 (BLASTP)  
- 基因家族聚类分析  
- 系统发育分析准备  

EOF  
    
    # 创建BLAST数据库构建脚本模板  
    cat > build_blast_db.sh << 'EOF'  
#!/bin/bash  

# 构建BLAST数据库的示例脚本  

GENOME_FASTA=$1  

if [ ! -f "$GENOME_FASTA" ]; then  
    echo "错误：基因组FASTA文件不存在"  
    exit 1  
fi  

echo "构建蛋白质BLAST数据库..."  
makeblastdb \  
    -in "$GENOME_FASTA" \  
    -dbtype nucl \  
    -parse_seqids \  
    > blast_db.log 2>&1  

if [ $? -eq 0 ]; then  
    echo "BLAST数据库构建完成"  
else  
    echo "BLAST数据库构建失败"  
fi  

EOF  
    
    chmod +x build_blast_db.sh  
    
    # 创建同源基因预测脚本模板  
    cat > predict_homologs.sh << 'EOF'  
#!/bin/bash  

# 同源基因预测脚本模板  

GENOME_FASTA=$1  
PROTEIN_DB=$2  

if [ ! -f "$GENOME_FASTA" ] || [ ! -f "$PROTEIN_DB" ]; then  
    echo "错误：输入文件不存在"  
    exit 1  
fi  

echo "运行蛋白质比对..."  

# tBLASTn搜索同源基因  
tblastn \  
    -query "$PROTEIN_DB" \  
    -db "$GENOME_FASTA" \  
    -evalue 1e-5 \  
    -outfmt 6 \  
    -num_threads 16 \  
    > homolog_blast_results.txt \  
    2> tblastn.log  

if [ $? -eq 0 ]; then  
    echo "同源基因预测完成，结果保存到: homolog_blast_results.txt"  
    
    # 统计结果  
    HITS=$(wc -l < homolog_blast_results.txt)  
    UNIQUE_QUERIES=$(cut -f1 homolog_blast_results.txt | sort | uniq | wc -l)  
    
    cat > homolog_summary.txt << SUMMARY_EOF  
# 同源基因预测摘要  

总比对结果数: $HITS  
唯一查询序列数: $UNIQUE_QUERIES  

SUMMARY_EOF  
    
else  
    echo "同源基因预测失败，请检查日志文件: tblastn.log"  
fi  

EOF  
    
    chmod +x predict_homologs.sh  
    
    echo "同源分析准备完成，相关脚本已创建:" | tee -a $LOG_FILE  
    echo "- 配置文件: homolog_analysis_config.txt" | tee -a $LOG_FILE  
    echo "- BLAST数据库构建脚本: build_blast_db.sh" | tee -a $LOG_FILE  
    echo "- 同源基因预测脚本: predict_homologs.sh" | tee -a $LOG_FILE  
    
}
