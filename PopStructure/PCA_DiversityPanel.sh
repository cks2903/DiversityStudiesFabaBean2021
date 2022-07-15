#!/bin/sh
#SBATCH -A norfab


plink --file $1 --allow-extra-chr --maf 0.01 --pca 685 --out pca

