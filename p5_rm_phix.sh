#!/bin/bash -l

#SBATCH -A snic2018-3-449 
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 6:00:00
#SBATCH -J rmphix
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

# Preprocessing pipeline
# Step 5: Map to PhiX, remove mapped reads (samtools view -f 4)

module load bioinfo-tools samtools bwa BEDTools

INDEXDIR=/proj/sllstore2017021/nobackup/HENRIQUE/phix_ref
DATADIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P4_rmdup_190614
OUTDIR=/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P5_rmphix_190614

mkdir -p $OUTDIR
cd $DATADIR
for i in *fastq.gz;
do
    bwa mem -t 5 $INDEXDIR/phix $i | samtools view -Sb -f 4 - > $OUTDIR/${i%dedup.fastq.gz}rmphix.bam
    bedtools bamtofastq -i $OUTDIR/${i%dedup.fastq.gz}rmphix.bam -fq $OUTDIR/${i%dedup.fastq.gz}rmphix.fastq
    pigz -p 5 $OUTDIR/${i%dedup.fastq.gz}rmphix.fastq;
done
