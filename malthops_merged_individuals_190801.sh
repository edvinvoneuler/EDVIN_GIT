#!/bin/bash -l
#SBATCH -A snic2018-3-449
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem1TB
#SBATCH -t 7-00:00:00
#SBATCH -J malthops_trial_fatnode
#SBATCH --mail-type=ALL
#SBATCH --mail-user edvinvoneuler@gmail.com


HOPS=/proj/sllstore2017021/nobackup/EDVIN/MALTHOPS/HOPS/AMPS/Resources

echo "Samples:"  $1 $2 $3 $4 $5

/proj/sllstore2017021/nobackup/EDVIN/jdk-13/bin/java -jar $HOPS/hops0.2.jar -input $1 $2 $3 $4 $5 -output /home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/malt_full_190801 -m malt -c /proj/sllstore2017021/nobackup/EDVIN/MALTHOPS/config.txt 
