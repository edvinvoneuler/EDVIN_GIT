#!/bin/bash -l
#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH -J adrm_Jena 
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

module load bioinfo-tools
module load AdapterRemoval 

OUTDIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P2_adrm_singleend/
mkdir -p $OUTDIR
cd $OUTDIR
for file in /home/edvo1850/DENTAL_CALC/gorilla_calculus/*/*_R1_*.fastq.gz
do
    pair1=$file
    pair2=$(echo $file | sed 's/_R1_/_R2_/')
    basename=$(echo $file | sed 's/_R1_/_/')
    basename=${basename##*/}
    basename=${basename%.fastq.gz}
    if [ ! -f "$pair2" ]
    then
    	AdapterRemoval --trimns --trimqualities --minquality 30 --minlength 30 --mm 3 --minalignmentlength 11 --file1 $pair1 --basename $basename 
    	#echo $pair1
    	#echo $pair2
    	#echo $basename
    fi
done

