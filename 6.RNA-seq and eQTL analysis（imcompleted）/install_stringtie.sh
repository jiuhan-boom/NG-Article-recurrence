#!/bin/bash

# StringTie安装脚本

echo "下载StringTie v2.20..."
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.0.Linux_x86_64.tar.gz

echo "解压文件..."
tar -xzf stringtie-2.2.0.Linux_x86_64.tar.gz

echo "安装到系统路径..."
sudo cp stringtie-2.2.0.Linux_x86_64/stringtie /usr/local/bin/

echo "验证安装..."
stringtie --version
        