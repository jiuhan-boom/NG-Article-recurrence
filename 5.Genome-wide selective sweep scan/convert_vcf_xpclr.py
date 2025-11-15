import sys
import gzip

def vcf_to_xpclr(vcf_file, output_prefix, sample_file):
    """将VCF文件转换为XP-CLR格式"""
    
    # 读取样本列表
    with open(sample_file, 'r') as f:
        target_samples = [line.strip() for line in f]
    
    # 打开VCF文件
    if vcf_file.endswith('.gz'):
        vcf_handle = gzip.open(vcf_file, 'rt')
    else:
        vcf_handle = open(vcf_file, 'r')
    
    # 读取样本头信息
    for line in vcf_handle:
        if line.startswith('#CHROM'):
            header = line.strip().split('\t')
            sample_names = header[9:]
            break
    
    # 找到目标样本的索引
    sample_indices = []
    for sample in target_samples:
        if sample in sample_names:
            sample_indices.append(sample_names.index(sample))
    
    print(f"找到 {len(sample_indices)} 个目标样本", file=sys.stderr)
    
    # 处理变异位点
    snp_count = 0
    with open(f"{output_prefix}.xpclr.txt", 'w') as out_file:
        vcf_handle.seek(0)  # 重新开始读取
        for line in vcf_handle:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            
            # 只处理双等位基因位点
            if ',' in alt:
                continue
            
            # 提取目标样本的基因型
            genotypes = []
            for idx in sample_indices:
                gt_field = fields[9 + idx].split(':')[0]
                if gt_field == '0/0':
                    genotypes.extend([ref, ref])
                elif gt_field == '0/1':
                    genotypes.extend([ref, alt])
                elif gt_field == '1/1':
                    genotypes.extend([alt, alt])
                else:
                    genotypes.extend(['N', 'N'])  # 缺失数据
            
            # 写入XP-CLR格式
            out_file.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t" + '\t'.join(genotypes) + '\n')
            snp_count += 1
    
    print(f"转换完成，共处理 {snp_count} 个SNP位点", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("使用方法: python vcf_to_xpclr.py <vcf_file> <output_prefix> <sample_file>")
        sys.exit(1)
    
    vcf_to_xpclr(sys.argv[1], sys.argv[2], sys.argv[3])
    