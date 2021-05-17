# PGS_construct
Pipeline of Polygenic Score (PGS) construction.

This repo contains codes for automation of PGS construction using various tools (e.g. PRSice-2, PRS-CS)
as well as logs documenting results from runs.

Note that this repo does not contain the PGS construction source code, i.e. PRSice-2 software.

## Standard polygenic score with PRSice
_Run_PRSice_cmd.R_

is an R wrapper to construt a simple, plain vanilla PGS (weighted average of all the SNPs, without clumping or anything) both for all the SNPs and genomewide significant SNPs only


Example call (example produces several PGS scores, one for every GWAS summary statistics ending in *.sumstat in the rGE folder, for the MCS data
$ head ~/biroli/geighei/data/GWAS_sumstats/clean/rGE/educYears.sumstats      #check the column names
$  nohup Rscript Run_PRSice_cmd.R
         --prsice ~/tools/PRSice-2/
         --geno ~/MCS/genomeclean/cleaned_MCS_genome
         --geno_name MCS
         --gwas_dir ~/biroli/geighei/data/GWAS_sumstats/clean/rGE/
         --gwas_names bmi.sumstats
         --out_dir ~/MCS/PGS/
         --chr CHR
         --bp POS
         --A1 A1
         --A2 A2
         --snp SNP
         --stat BETA
         --pvalue P
         &
