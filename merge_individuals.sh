#!/bin/bash -l
#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 4:00:00
#SBATCH -J merge_individuals
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

outdir=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P6_merged_samples_unmapped
for file in /home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P6_data_unmapped_human/*.fastq
do
    basename=$(basename $file)
    ind_name=$(echo $basename | cut -c1-6)
    
    echo "Input file:" $file "Concatenated into" $outdir/$ind_name.fastq
    cat $file >> $outdir/$ind_name.fastq
done