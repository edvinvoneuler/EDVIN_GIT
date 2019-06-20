#!/bin/bash -l

#SBATCH -A snic2018-3-449 
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 4:00:00
#SBATCH -J rmdup
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

# Preprocessing pipeline
# Step 4: Remove duplicates with Tom's inhouse script

# may need to increase memory depending on input file size

DATADIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P3_prinseq_SE_190619
OUTDIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P4_dedup_SE_190620

rmdup_script=/proj/sllstore2017021/nobackup/EXAMPLES_for_new_users/remove_duplicates_single_end.py

mkdir -p $OUTDIR

for i in $DATADIR/*_passed.fastq.gz;
do
    cd $OUTDIR
    python $rmdup_script $i ${i%_passed.fastq.gz};
done

# output is not compressed, compress with gzip in parallel (pigz)

for i in $OUTDIR/*.fastq;
do
    pigz -p 5 $i;
done
