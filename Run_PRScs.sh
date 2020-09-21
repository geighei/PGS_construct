#!/usr/bin/env bash
shopt -s extglob

##  --------------  File Paths ------------
# run inputs
n=575705
phenotype="smokeInit"
dataset="HRS"
bim="cleaned_phg000515-v1-HRS-phase123-bestGuessImp"

# file paths
prscs="/home/ubuntu/tools/PRScs"
plink="/home/ubuntu/tools/plink/plink"
geno="/home/ubuntu/biroli/geighei/data/$dataset/genomeclean/$bim"
sumstat="/home/ubuntu/biroli/geighei/data/GWAS_sumstats/clean/rGE/$phenotype.sumstats"
out_prscs="/home/ubuntu/biroli/geighei/data/$dataset/PGS/PRScs/$phenotype"
out_prscs_name="$out_prscs/$phenotype"
appended_prscs=$out_prscs_name"_PRScs_estimates_ALL.txt"
pgs=$out_prscs_name"_PGS"

##  --------------  Run PRS-CS (estimate effect size) ------------
python $prscs/PRScs.py \
--ref_dir=$prscs/ldblk_1kg_eur \
--bim_prefix=$geno \
--sst_file=$sumstat \
--n_gwas=$n \
--out_dir=$out_prscs_name

##  --------------  Append PRS-CS output files ------------
find "$out_prscs_name"*chr*.txt | xargs cat > $appended_prscs

##  --------------  Run Plink (polygenic risk score) ------------

$plink \
--bfile $geno \
--score $appended_prscs 2 4 6 \
--out $pgs
