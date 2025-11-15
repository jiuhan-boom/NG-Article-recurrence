# GCTA安装脚本

echo "下载GCTA v1.94.1..."
wget https://cnsgenomics.com/software/gcta/bin/gcta_1.94.1beta_linux_kernel_4_x86_64.zip

echo "解压文件..."
unzip gcta_1.94.1beta_linux_kernel_4_x86_64.zip

echo "安装到系统路径..."
sudo cp gcta_1.94.1beta/gcta64 /usr/local/bin/

echo "验证安装..."
gcta64 --help
      