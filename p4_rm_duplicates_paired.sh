#!/bin/bash -l

#SBATCH -A snic2018-3-449 
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 2:00:00
#SBATCH -J rmdup_paired_ends
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

# may need to increase memory depending on input file size

DATADIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P3_prinseq_uncollapsed_190710
OUTDIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P4_dedup_uncollapsed_190718

rmdup_script=/crex/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/SCRIPTS/EDVIN_GIT/rmdup_paired_end.py

mkdir -p $OUTDIR
cd $DATADIR

for file in *passed_1.fastq.gz;
do
    rev_read=${file/_1./_2.}
    echo -e "fwd_read: $file\nrev_read: $rev_read\nout_file: $OUTDIR/${file%_passed_1.fastq.gz}_dedup.fastq"
    python $rmdup_script $DATADIR/$file $DATADIR/$rev_read $OUTDIR/${file%_passed_1.fastq.gz}
    pigz -p 5 *.fastq
done