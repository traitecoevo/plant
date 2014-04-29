#!/bin/bash

for d in 7.5 10.6 20 30 60 120; do
 nohup Rscript -e "source('scripts/assembly_lma.R'); assemble(time.disturbance=$d, prod=1.0, n.steps=1000, output.dir='output');" >> output/out-$d.txt 2>&1 &
done
