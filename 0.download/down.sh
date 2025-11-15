#!/bin/bash

# 检查文件是否存在
if [ ! -f "data_download_links_CNP0003213_ftp.txt" ]; then
    echo "错误：文件 data_download_links_CNP0003213_ftp.txt 不存在"
    exit 1
fi

# 创建下载目录
mkdir -p downloads
cd downloads

# 逐行读取FTP链接并下载
while IFS= read -r url; do
    # 跳过空行和注释行
    if [[ -n "$url" && ! "$url" =~ ^[[:space:]]*# ]]; then
        echo "正在下载: $url"
        wget -c "$url" --show-progress
        if [ $? -eq 0 ]; then
            echo "下载成功: $url"
        else
            echo "下载失败: $url"
        fi
    fi
done < "../data_download_links_CNP0003213_ftp.txt"

echo "批量下载完成！"
