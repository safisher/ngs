#!/bin/bash
# Jamie Shallcross, 2016

#set -e
cd src/1* #Make sure to change sample 


if [[ ! -z $1 ]]; then
    mask="--use-bases-mask $1"
    echo $mask
fi 

bcl2fastq --ignore-missing-bcl --ignore-missing-control --input-dir Data/Intensities/BaseCalls --barcode-mismatches 1 $mask -o ../../raw 

for s in `grep ",,,," SampleSheet.csv | cut -d',' -f1`; do 
    mkdir -v ../../raw/Sample_$s
    mv -v ../../raw/${s}_S*.gz ../../raw/Sample_$s/
done

chgrp -R repo-admin ../../raw/*

