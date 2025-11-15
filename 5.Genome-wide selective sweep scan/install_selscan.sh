# selscan安装脚本

echo "下载selscan v1.3.0..."
wget https://github.com/szpiech/selscan/archive/v1.3.0.tar.gz

echo "解压文件..."
tar -xzf v1.3.0.tar.gz

echo "编译selscan..."
cd selscan-1.3.0
make

echo "安装到系统路径..."
sudo cp selscan /usr/local/bin/
sudo cp norm /usr/local/bin/

echo "验证安装..."
selscan --help