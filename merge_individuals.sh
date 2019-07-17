#!/bin/bash

outdir=/home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P6_merged_samples_unmapped
for file in /home/edvo1850/DENTAL_CALC/DENTAL_CALCULUS_JENA_GORILLA_181009/P6_data_unmapped_human
do
    basename=$(basename $file)
    ind_name=$(echo $basename | cut -c1-6)

    cat $file >> $outdir/$ind_name