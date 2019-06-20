#!/bin/bash -l

#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 5:00:00
#SBATCH -J prinsq_jena
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

module load bioinfo-tools
module load prinseq

# Preprocessing pipeline
# Step 3: Quality filter reads with prinseq
# remove reads with mean quality < 30

DATADIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P2_adrm_singleend/
OUTDIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P3_prinseq_SE_190619/

mkdir -p $OUTDIR
cd $DATADIR
for i in *.truncated;
do
    cat $i | prinseq -fastq stdin -min_qual_mean 30 -out_good $OUTDIR/${i%.truncated}_passed -out_bad $OUTDIR/${i%.truncated}_failed -log $OUTDIR/log_prinseq.log
    # output is not compressed, so compress with gzip in parallel
    pigz -p 5 ${OUTDIR}/${i%.truncated}_passed.fastq
    pigz -p 5 ${OUTDIR}/${i%.truncated}_failed.fastq
done
