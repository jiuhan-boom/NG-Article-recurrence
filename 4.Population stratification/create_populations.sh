# 根据样本信息创建群体分组文件  
# 格式: SampleID\tPopulation  

echo "Sample1	SweetCorn" > sweet_corn.pop  
echo "Sample2	SweetCorn" >> sweet_corn.pop  
echo "Sample3	FieldCorn" >> field_corn.pop  
echo "Sample4	FieldCorn" >> field_corn.pop  

# 合并为一个群体文件  
cat sweet_corn.pop field_corn.pop > all_populations.txt  
