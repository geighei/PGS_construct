#!/usr/bin/env bash
shopt -s extglob

SUMSTATS="/home/ubuntu/biroli/geighei/data/GWAS_sumstats/clean/rGE/split_sample"
PGS_CODE="/home/ubuntu/biroli/geighei/code/PGS_construct/"
GWAS_PATTERN="*"

## gunzip all .gz files
gunzip $SUMSTATS/*.gz

pushd $PGS_CODE
for file in $SUMSTATS/$GWAS_PATTERN
do
## ELSA
	Rscript "/home/ubuntu/tools/PRSice-2/PRSice.R"  \
	--dir "/home/ubuntu/biroli/geighei/data/ELSA/PGS/PRSice/split_sample" \
	--prsice "/home/ubuntu/tools/PRSice-2/PRSice_linux" \
	--base "$file" \
	--snp SNP --A1 EFFECT_ALLELE --A2 OTHER_ALLELE --stat BETA --pvalue PVAL \
	--target "/home/ubuntu/biroli/geighei/data/ELSA/genomeclean/cleaned_ega-box-163_ForwardStrand_excREL" \
	--beta --no-regress T \
	--out "/home/ubuntu/biroli/geighei/data/ELSA/PGS/PRSice/split_sample/$(basename "$file")" --bar-levels 1e-2 --fastscore --no-clump
	# $(basename xxx) convention isolates filename from path
	
## HRS
	Rscript "/home/ubuntu/tools/PRSice-2/PRSice.R"  \
	--dir "/home/ubuntu/biroli/geighei/data/HRS/PGS/PRSice/split_sample" \
	--prsice "/home/ubuntu/tools/PRSice-2/PRSice_linux" \
	--base "$file" \
	--snp SNP --A1 EFFECT_ALLELE --A2 OTHER_ALLELE --stat BETA --pvalue PVAL \
	--target "/home/ubuntu/biroli/geighei/data/HRS/genomeclean/cleaned_phg000515-v1-HRS-phase123-bestGuessImp" \
	--beta --no-regress T \
	--out "/home/ubuntu/biroli/geighei/data/HRS/PGS/PRSice/split_sample/$(basename "$file")" --bar-levels 1e-2 --fastscore --no-clump


## WLS
	Rscript "/home/ubuntu/tools/PRSice-2/PRSice.R"  \
	--dir "/home/ubuntu/biroli/geighei/data/WLS/PGS/PRSice/split_sample" \
	--prsice "/home/ubuntu/tools/PRSice-2/PRSice_linux" \
	--base "$file" \
	--snp SNP --A1 EFFECT_ALLELE --A2 OTHER_ALLELE --stat BETA --pvalue PVAL \
	--target "/home/ubuntu/biroli/geighei/data/WLS/genomeclean/cleanedWLS_bestGuessImp" \
	--beta --no-regress T \
	--out "/home/ubuntu/biroli/geighei/data/WLS/PGS/PRSice/split_sample/$(basename "$file")" --bar-levels 1e-2 --fastscore --no-clump
done
popd