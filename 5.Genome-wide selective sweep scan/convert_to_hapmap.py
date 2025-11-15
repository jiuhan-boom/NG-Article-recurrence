import sys
import gzip

def vcf_to_hapmap(vcf_file, output_prefix, population_name):
    """将VCF转换为hapmap格式用于selscan"""
    
    # 打开VCF文件
    if vcf_file.endswith('.gz'):
        vcf_handle = gzip.open(vcf_file, 'rt')
    else:
        vcf_handle = open(vcf_file, 'r')
    
    # 读取样本头信息
    sample_names = []
    for line in vcf_handle:
        if line.startswith('#CHROM'):
            sample_names = line.strip().split('\t')[9:]
            break
    
    print(f"样本数量: {len(sample_names)}", file=sys.stderr)
    
    # 创建hapmap输出文件
    hapmap_file = f"{output_prefix}.hapmap"
    
    with open(hapmap_file, 'w') as out_file:
        # 写入hapmap头信息
        header = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#", 
                  "center", "protLSID", "assayLSID", "panelLSID", "QCcode"] + sample_names
        out_file.write('\t'.join(header) + '\n')
        
        # 处理每个SNP位点
        snp_count = 0
        vcf_handle.seek(0)  # 重新开始读取
        
        for line in vcf_handle:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            id_ = f"{chrom}_{pos}"
            ref = fields[3]
            alt = fields[4]
            
            # 只处理双等位基因位点
            if ',' in alt:
                continue
            
            alleles = f"{ref}/{alt}"
            
            # 提取基因型信息
            genotypes = []
            for i in range(9, len(fields)):
                gt_field = fields[i].split(':')[0]
                if gt_field == '0/0':
                    genotypes.append(f"{ref}{ref}")
                elif gt_field == '0/1':
                    genotypes.append(f"{ref}{alt}")
                elif gt_field == '1/1':
                    genotypes.append(f"{alt}{alt}")
                else:
                    genotypes.append("NN")  # 缺失数据
            
            # 写入hapmap格式行
            row_data = [id_, alleles, chrom, pos, "+", "NA", "NA", "NA", "NA", "NA", "NA"] + genotypes
            out_file.write('\t'.join(row_data) + '\n')
            
            snp_count += 1
            
            if snp_count % 10000 == 0:
                print(f"已处理 {snp_count} 个SNP位点", file=sys.stderr)
    
    print(f"转换完成，共处理 {snp_count} 个SNP位点", file=sys.stderr)
    
    return hapmap_file

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("使用方法: python convert_to_hapmap.py <vcf_file> <output_prefix> <population_name>")
        sys.exit(1)
    
    vcf_to_hapmap(sys.argv[1], sys.argv[2], sys.argv[3])