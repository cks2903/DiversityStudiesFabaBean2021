#!/bin/bash
#SBATCH -A norfab
#SBATCH --mem 200g
vcftools --vcf Mapped_HighQualSNPs2022-01-06_imputed.vcf --SNPdensity 1000000 --chr chr1S --out SNPdensities_chr1S # calculate the average number of SNPs pr. 1 gbp
vcftools --vcf Mapped_HighQualSNPs2022-01-06_imputed.vcf --SNPdensity 1000000 --chr chr1L --out SNPdensities_chr1L # calculate the average number of SNPs pr. 1 gbp
vcftools --vcf Mapped_HighQualSNPs2022-01-06_imputed.vcf --SNPdensity 1000000 --chr chr2 --out SNPdensities_chr2 # calculate the average number of SNPs pr. 1 gbp
vcftools --vcf Mapped_HighQualSNPs2022-01-06_imputed.vcf --SNPdensity 1000000 --chr chr3 --out SNPdensities_chr3 # calculate the average number of SNPs pr. 1 gbp
vcftools --vcf Mapped_HighQualSNPs2022-01-06_imputed.vcf --SNPdensity 1000000 --chr chr4 --out SNPdensities_chr4 # calculate the average number of SNPs pr. 1 gbp
vcftools --vcf Mapped_HighQualSNPs2022-01-06_imputed.vcf --SNPdensity 1000000 --chr chr5 --out SNPdensities_chr5 # calculate the average number of SNPs pr. 1 gbp
vcftools --vcf Mapped_HighQualSNPs2022-01-06_imputed.vcf --SNPdensity 1000000 --chr chr6 --out SNPdensities_chr6 # calculate the average number of SNPs pr. 1 gbp