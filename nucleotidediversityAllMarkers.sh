#!/bin/bash
#SBATCH -A norfab
#SBATCH --mem 200g
#SBATCH -c 16


vcftools --vcf Mapped_HighQualSNPs2022_RedFiltered20220126.vcf --keep $1 --site-pi --out $2
