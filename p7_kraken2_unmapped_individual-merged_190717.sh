#!/bin/bash -l

#SBATCH -A snic2018-3-449
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 48:00:00
#SBATCH -J kraken_jena_unmapped-host_individuals-merged

# Microbial assignments using kraken
# Standard kraken database (RefSeq bacteria+archaea+viral+human)
# build by UPPMAX $KRAKEN_DB (updated started of each month)

# Unmerged and merged

module load bioinfo-tools Kraken2
echo $(module list)
echo $(readlink -f $KRAKEN2_DEFAULT_DB)

DATDIR=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P6_merged_samples_unmapped
OUTDIR=/crex/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P7_kraken2_merged-indiviuals_190717/

# Set up database
MY_DB_DIR=$SNIC_TMP/Kraken2
MY_DB=$MY_DB_DIR/${KRAKEN2_DEFAULT_DB##*/}
mkdir -p $MY_DB
cp -av $KRAKEN2_DEFAULT_DB/* $MY_DB/


mkdir $OUTDIR
cd $DATDIR

for i in *.fastq
do
	echo "Sample: ${i%.fastq}"
    	kraken2 --db $MY_DB $i --threads 20 --fastq-input --report $OUTDIR/${i%.fastq}kraken2_report.txt --report-zero-counts --output $OUTDIR/${i%.fastq}kraken2_output.txt;
done
