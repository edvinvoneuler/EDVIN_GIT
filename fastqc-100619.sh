#!/bin/bash -l
#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 1:00:00
#SBATCH -J fastQC_Jena
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

#module load bioinfo-tools
#module load fastQC

DATDIR=/home/edvo1850/DENTAL_CALC/gorilla_calculus
OUTDIR=/home/edv01850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/FASTQC


for file in ${DATDIR}/*/*
do	
	sample=$(basename $(dirname $file))
	read=$(basename $file)
	mkdir -p $OUTDIR/$sample/$read
    	fastqc $file -o $OUTDIR/$sample/$read
	#echo $OUTDIR/$sample/$read
done
