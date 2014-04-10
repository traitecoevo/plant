#!/bin/bash

mkdir output
for d in 7.5 10.6 15 21.2 30 42.4 60 84.8 120 169.7 240; do
 nohup Rscript -e "source('scripts/batch-fun.R'); out<-find_equilibirum($d); saveRDS(out,paste0('output/equil-',$d,'.rds'));" >> output/out-$d.txt 2>&1 &
done
