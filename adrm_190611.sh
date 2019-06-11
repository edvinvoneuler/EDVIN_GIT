#!/bin/bash -l
#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH -J fastQC_Jena
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

module load bioinfo-tools
module load AdapterRemoval 

DATDIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/RAW_FQ/ftp.shh.mpg.de/private/gorilla_calculus/
OUTDIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P2_adrm_unmerged/

AdapterRemoval --minlength 30 /home/edvo1850/DENTAL_CALC/gorilla_calculus/MTM010.A010/MTM010.A0101_S0_L006_R1_001.fastq.gz 
 



