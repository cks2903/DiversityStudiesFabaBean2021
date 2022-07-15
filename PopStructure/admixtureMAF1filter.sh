#!/bin/bash
#SBATCH -A norfab
#SBATCH -c 16
#SBATCH -t 48:00:00

for K in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; \
	do admixture --cv=10 genotypesADMIXTURE_MAFfiltered.bed $K | tee logMAFFilter${K}.out; done