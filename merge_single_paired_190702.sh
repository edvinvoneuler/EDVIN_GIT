#!/bin/bash -l
#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J merge single end - paired end
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

for file in /home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/data_unmapped_SE_human/*.gz
do
	base=$(basename $file)
	sample=${base%%_*.gz}
	outfile=`find /home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/data_unmapped_human/ -regex ".*/${sample}.*"`
	zcat $file >> $outfile 
	echo $file "concatenated into" $outfile
done
echo "Done"
