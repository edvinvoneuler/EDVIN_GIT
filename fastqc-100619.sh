#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 1:00:00
#SBATCH -J fastQC
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

load bioinfo-tools fastQC

DATDIR=/home/edvo1850/DENTAL_CALC/gorilla_calculus
OUTDIR=/home/edv01850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/FASTQC


for file in ${DATDIR}/*/*
do	
	name=${file##*/}
        base=${name%.fastq.gz}
	mkdir ${OUTDIR}/$base 	
    	fastqc $file -o ${OUTDIR}/$base
	#echo  $file -o ${OUTDIR}/$base
done
