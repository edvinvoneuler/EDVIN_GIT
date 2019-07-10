#!/bin/bash -l

#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 5:00:00
#SBATCH -J prinsq_jena_uncollapsed
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

module load bioinfo-tools
module load prinseq

# Preprocessing pipeline

# Step 3: Quality filter reads with prinseq
# remove reads with mean quality < 30

DATADIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P2_adrm_190612
OUTDIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P3_prinseq_uncollapsed_190710/

mkdir -p $OUTDIR
for forward_read in $DATADIR/*pair1*
do
	reverse_read=`echo $forward_read | sed 's/pair1/pair2/'`
	forward_base=$(basename $forward_read)
	reverse_base=$(basename $reverse_read)
  	echo -e "forward read: $forward_base, reverse read: $reverse_base"
	prinseq -fastq $forward_read -fastq2 $reverse_read -min_qual_mean 30 -out_good $OUTDIR/${forward_base%%_*}_passed -out_bad $OUTDIR/${forward_base%%_*}_failed -log $OUTDIR/log_prinseq.log
    
done
pigz -p 5 ${OUTDIR}/*

