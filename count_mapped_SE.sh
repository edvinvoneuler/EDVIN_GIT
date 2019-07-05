#!/bin/bash -l
#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH -J count_mapped_SE_reads
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

module load bioinfo-tools samtools 

outdir=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/data_mapped_SE_human
echo -e "Sample\tReads" > $outdir/p6_Mapped_Jena_readcount.txt

for file in $outdir/*.bam
do
	readcount=`samtools view -c $file`
	base=$(basename $file)
	name=${base%.bam}
	echo -e "${name}\t${readcount}" >> $outdir/p6_Mapped_Jena_readcount.txt
done


outdir=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/data_unmapped_SE_human
echo -e "Sample\tReads" > $outdir/p6_Unmapped_Jena_readcount.txt

for file in $outdir/*.bam
do
	readcount=`samtools view -c $file`
	base=$(basename $file)
	name=${base%.bam}
	echo -e "${name}\t${readcount}" >> $outdir/p6_Unmapped_Jena_readcount.txt

done
