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

DATADIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P2_adrm_merged_190613/
OUTDIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P3_prinseq_190613/

#mkdir $OUTDIR
cd $DATADIR
for i in *.gz;
do
    zcat $i | prinseq -fastq stdin -min_qual_mean 30 -out_good $OUTDIR/${i%_m_concat.gz}_passed -out_bad $OUTDIR/${i%_m_concat.gz}_failed -log $OUTDIR/log_prinseq.log
    
    # output is not compressed, so compress with gzip in parallel
    pigz -p 5 ${OUTDIR}/${i%_m_concat.gz}_passed.fastq
    pigz -p 5 ${OUTDIR}/${i%_m_concat.gz}_failed.fastq
done
