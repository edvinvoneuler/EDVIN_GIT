#!/bin/bash -l

#SBATCH -A snic2018-3-449 
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 48:00:00
#SBATCH -J maphuman_Jena
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

# Preprocessing pipeline
# Step 6: Map to human genome
# Save mapped reads for future human contamination analysis
# Save unmapped reads for future microbial analysis

module load bioinfo-tools samtools bwa BEDTools

INDEXDIR=/proj/sllstore2017021/nobackup/JAELLE/REFERENCES/gorilla_gorilla_human/calculushost
DATADIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P5_rmphix_SE_190620
UNMAPDIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/data_unmapped_SE_human
MAPDIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/data_mapped_SE_human

mkdir $MAPDIR
mkdir $UNMAPDIR
cd $DATADIR
for i in *fastq.gz;
do
    # mapped
    bwa mem -t 10 $INDEXDIR $i | samtools view -Sb -F 4 -@ 10 - > $MAPDIR/${i%rmphix.fastq.gz}human_mapped.bam;
done

for i in *fastq.gz;
do
    # unmapped
    bwa mem -t 10 $INDEXDIR $i | samtools view -Sb -f 4 -@ 10 - > $UNMAPDIR/${i%rmphix.fastq.gz}unmapped.bam
    bedtools bamtofastq -i $UNMAPDIR/${i%rmphix.fastq.gz}unmapped.bam -fq $UNMAPDIR/${i%rmphix.fastq.gz}unmapped.fastq
    pigz -p 10 $UNMAPDIR/${i%rmphix.fastq.gz}unmapped.fastq;
done
