
datadir=/crex/proj/sllstore2017021/nobackup/JAELLE/DENTAL_CALCULUS_JENA_GORILLA_181009/P6_merged_individuals_unmapped


for i in `seq 1 5 21`
do
	input=`find $datadir  -type f | sed -n "$i,$((i+4))p" | sed ':a;N;$!ba;s/\n/ /g'`
	echo -e "Input files:" $input "\n"
	sbatch malthops_merged_individuals_190801.sh $input
done
