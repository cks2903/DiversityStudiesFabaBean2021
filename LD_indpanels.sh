#!/bin/sh
#SBATCH -A norfab


plink --file $1 --allow-extra-chr --chr $2 --maf 0.05 --r2 --ld-window 2000000000 --ld-window-kb 2000000 --ld-window-r2 0 --out $3

