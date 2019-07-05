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

DATADIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P3_prinseq_190613
OUTDIR=/proj/sllstore2017020/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P4_dedup_190614

rmdup_script=/proj/sllstore2017021/nobackup/EXAMPLES_for_new_users/remove_duplicates_single_end.py

mkdir $OUTDIR
cd $DATADIR

for i in *passed.fastq.gz;
do
    cd $OUTDIR
    python $rmdup_script $DATADIR/$i ${i%_filt30_passed.fastq.gz};
done

# output is not compressed, compress with gzip in parallel (pigz)
cd $OUTDIR
for i in *.fastq;
do
    pigz -p 5 $i;
done
