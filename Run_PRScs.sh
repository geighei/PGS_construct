#!/usr/bin/env bash
shopt -s extglob

##  -------------- Read in options ------------
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -n|--n) n=$2; shift ;;
        -s|--sumstat) sumstat="$2"; shift ;;
        -g|--geno) geno="$2"; shift ;;
        -o|--out) out_prscs="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

#  --------------  File Paths ------------
# other file paths
prscs="/home/ubuntu/tools/PRScs"
plink="/home/ubuntu/tools/plink/plink2"
appended_prscs=$out_prscs"_PRScs_estimates_ALL.txt"
pgs=$out_prscs"_PGS"

##  --------------  Run PRS-CS (estimate effect size) ------------
python $prscs/PRScs.py \
--ref_dir=$prscs/ldblk_1kg_eur \
--bim_prefix=$geno \
--sst_file=$sumstat \
--n_gwas=$n \
--out_dir=$out_prscs

##  --------------  Append PRS-CS output files ------------
find "$out_prscs"*chr*.txt | xargs cat > $appended_prscs

##  --------------  Run Plink (polygenic risk score) ------------
$plink \
--bfile $geno \
--score $appended_prscs 2 4 6 \
--out $pgs
