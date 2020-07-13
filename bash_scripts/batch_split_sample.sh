#!/usr/bin/env bash
shopt -s extglob

SUMSTATS="/home/ubuntu/biroli/geighei/data/GWAS_sumstats/clean/rGE/split_sample"
PGS_SCRIPT="/home/ubuntu/biroli/geighei/code/PGS_construct/Run_PRSice_Cmd.R"
GWAS_PATTERN="*.gz"

## gunzip all .gz files
gunzip $SUMSTATS/*.gz


## ELSA
Rscript $PGS_SCRIPT --prsice ~/tools/PRSice-2/ \
--geno /home/ubuntu/biroli/geighei/data/ELSA/genomeclean/cleaned_ega-box-163_ForwardStrand_excREL \
--geno_name ELSA --gwas_dir /home/ubuntu/biroli/geighei/data/GWAS_sumstats/clean/rGE/split_sample/ \
--gwas_pattern $GWAS_PATTERN --out_dir /home/ubuntu/biroli/geighei/data/ELSA/PGS/PRSice/split_sample/ \
--combine

## HRS


## WLS