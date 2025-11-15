# XP-CLR安装脚本

echo "克隆XP-CLR仓库..."
git clone https://github.com/hardingnj/XPEHH.git

echo "编译XP-CLR..."
cd XPEHH
make

echo "安装到系统路径..."
sudo cp xpclr /usr/local/bin/

echo "验证安装..."
xpclr --help