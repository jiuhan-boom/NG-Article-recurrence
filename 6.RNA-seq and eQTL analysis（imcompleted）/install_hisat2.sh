#!/bin/bash

# Hisat2安装脚本

echo "下载Hisat2 v2.1.0..."
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip

echo "解压文件..."
unzip hisat2-2.1.0-Linux_x86_64.zip

echo "安装到系统路径..."
sudo cp hisat2-2.1.0/hisat2* /usr/local/bin/

echo "验证安装..."
hisat2 --version 