#!/bin/bash -l
#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH -J fastQC_Jena
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

module load bioinfo-tools
module load FastQC

DATDIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/RAW_FQ/ftp.shh.mpg.de/private/gorilla_calculus/
OUTDIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/FASTQC


for file in ${DATDIR}/*/*
do	
	sample=$(basename $(dirname $file))
	read=$(basename $file)
	read=${read%.fastq.gz}
	mkdir -p $OUTDIR/$sample/$read
    	fastqc $file -o $OUTDIR/$sample/$read
	#echo $OUTDIR/$sample/$read
done
