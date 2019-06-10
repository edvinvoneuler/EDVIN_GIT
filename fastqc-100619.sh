#SBATCH -A snic2018-3-449
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 1:00:00
#SBATCH -J fastQC
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com

load bioinfo-tools fastQC

DATDIR=/home/edvo1850/DENTAL_CALC/gorilla_calculus

for file in ${DATDIR}/*/*
    print file