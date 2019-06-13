#!/bin/bash -l

#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J concatenate_collapsed
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

for file in /home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P2_adrm_190612/\w{6}\.\w{5}_S\d_L\d{3}_\d{3}.collapsed$
do
    #cat $file >> ${file}_m_concat.fastq.gz
    #cat ${file}.truncated >> ${file}_m_concat.fastq.gz
    echo $file
done