#!/bin/bash -l

#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J concatenate_collapsed
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

OUTDIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P2_adrm_merged_190613
i=0
for file in `find /home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P2_adrm_190612/ -regex '.+collapsed'`
do
    basename=${file%_S0_L00[0-9]_00[0-9].collapsed}
    outfile=$(basename $basename)_m_concat.fastq.gz
    #cat $file >> $OUTDIR/$outfile
    #cat ${file}.truncated >> $OUTDIR/$outfile
    echo "infile0" $(basename $file) 
    echo "infile2" $(basename $file.truncated)
    echo "concatenated into" $outfile
    echo $i
    i=$((i+1))
done
echo "Number of concatenated files: $i"
