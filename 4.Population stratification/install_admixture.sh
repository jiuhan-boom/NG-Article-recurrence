# ADMIXTURE安装脚本

echo "下载ADMIXTURE v1.3.0..."
wget http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz

echo "解压文件..."
tar -xzf admixture_linux-1.3.0.tar.gz

echo "安装到系统路径..."
sudo cp admixture_linux-1.3.0/admixture /usr/local/bin/

echo "验证安装..."
admixture --help