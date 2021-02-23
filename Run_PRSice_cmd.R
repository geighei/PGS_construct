##################################
######
###### Automated PRS Calculation
###### Author: Jeremy Vollen
###### Date: 29.10.2019
######
###### Runs PRSice for two p-values: 5e-8 and 1
###### for arbitrarily long list of target datasets and arbitrarily long list of phenotype/GWAS
###### 
##################################

## ----------------------------- Import Libraries, Macros, Sources -------------------------------

# Clear workspace and load necessary libraries
rm(list = ls())
gc()
library(tidyverse)
library(lubridate)
library(rstudioapi)
library(optparse)

## Import combining functions
source("combine_prsice.R")

## Define inputs; these might be user-inputted in the future
pval_thresholds <- as.character(c(5e-8, 1))


## ----------------------------- Import Libraries, Macros, Sources -------------------------------

# Read Command-line arguments
option_list = list(
  make_option(c("--prsice"), type="character", default=NULL, 
              help="PRSice path, e.g. ~/plink/PRSice-2/", metavar="character"),
  make_option(c("--geno"), type="character", default=NULL, 
              help="Cleaned genome data, e.g. ~/biroli/geighei/data/ELSA/genomeclean/cleaned_ega-box-163_ForwardStrand_excREL", metavar="character"),
  make_option(c("--geno_name"), type="character", default="", 
              help="Name of genome dataset, e.g. ELSA", metavar="character"),
  make_option(c("--gwas_dir"), type="character", default=NULL, 
              help="Path to sumstats, e.g. ~/biroli/geighei/data/GWAS_sumstats/clean/rGE/", metavar="character"),
  make_option(c("--gwas_names"), type="character", default=NULL, 
              help="Sumstats files separated by comma, e.g. bmi.sumstats,ea.txt", metavar="character"),
  make_option(c("--gwas_pattern"), type="character", default=".sumstats",
              help="Regex pattern to match sumstats to generate PGS for. Defaults to all with .sumstats suffix.", metavar="character"),
  make_option(c("--out_dir"), type="character", default=NULL, 
              help="Path for PRS output, e.g. ~/biroli/geighei/data/ELSA/PGS/", metavar="character"),
  make_option(c("--combine"), action="store_true", type="character", default=FALSE, 
              help="Use this flag if you want individual phenotype results combined into one file. Be careful this doesn't overwrite an existing file."),
  make_option(c("--snp"), type="character", default="SNP", help="Name of SNP column in sumstat. Default: 'SNP'.", metavar="character"),
  make_option(c("--chr"), type="character", default="CHR", help="Name of chromosome column in sumstat. Default: 'CHR'.", metavar="character"),
  make_option(c("--bp"), type="character", default="POS", help="Name of position column in sumstat. Default: 'POS'.", metavar="character"),
  make_option(c("--A1"), type="character", default="A1", help="Name of effect allele column in sumstat. Default: 'A1'.", metavar="character"),
  make_option(c("--A2"), type="character", default="A2", help="Name of other allele column in sumstat. Default: 'A2'.", metavar="character"),
  make_option(c("--stat"), type="character", default="BETA", help="Name of effect size column in sumstat. Default: 'BETA'.", metavar="character"),
  make_option(c("--pvalue"), type="character", default="P", help="Name of p-value column in sumstat. Default: 'P'.", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check for missing arguments and return error statement
if (is.null(opt$prsice)) {
  print_help(opt_parser)
  stop("Use --prsice option to point to PRSice software directory.", call.=FALSE)
}
if (is.null(opt$geno)) {
  print_help(opt_parser)
  stop("Use --geno option to point to cleaned genome data.", call.=FALSE)
}
if (is.null(opt$geno_name)) {
  print_help(opt_parser)
  stop("Use --geno_name option to point to cleaned genome data.", call.=FALSE)
}
if (is.null(opt$gwas_dir)) {
  print_help(opt_parser)
  stop("Use --gwas_dir option to point to munged summary statistics files.", call.=FALSE)
}
if (is.null(opt$out_dir)) {
  print_help(opt_parser)
  stop("Use --out_dir option to specify the path where PRSice should write results.")
}

# If GWAS names is null, we list all files ending in '.sumstats'; if not null, 
# we split the string on commas or semicolons and form in list of file names
if (is.null(opt$gwas_names)) {
  writeLines("Argument --gwas_names, used to specify summary statistic files to run PRSice on, not detected. 
             By default, we use --gwas_pattern argument to match all sumstats matching pattern, which defaults to all files with .sumstats ending if not provided.")
  # We find all files in the given GWAS directory recursively and filter on those matching the pattern given or default '.sumstats'
  gwas_names <- list.files(opt$gwas_dir, pattern = opt$gwas_pattern, recursive = TRUE)
} else {
  gwas_names <- str_split(opt$gwas_names, pattern = ",|;")[[1]]
}

# Ascertain the PRSice executable from the path provided with the --prsice argument.
# Note: I may need to add the PRSice executable used on Windows.. I'm not sure at the moment what that's called as I don't ever use Windows
prsice_exe_list <- list.files(opt$prsice, "mac|linux")
if (length(prsice_exe_list) > 1) {
  prsice_exe <- prsice_exe_list[[1]]
  warning("More than one PRSice executable detected. Using ", prsice_exe, ".")
} else if (length(prsice_exe_list) == 0) {
  stop("No PRSice executable found. Please make sure a PRSice executable exists in the PRSice path provided with the --prsice option.")
} else {
  prsice_exe <- prsice_exe_list[[1]]
}

# set working directory to output directory since PRSice sometimes doesn't heed the --dir argument 
# and instead writes results to the working directory
setwd(opt$out_dir)
## Get path of where this file is saved to ascertain what machine we're running it on
current_dir_vec <- str_split(getwd(), "/")[[1]]

# These paths differ by machine and are not hosted in geighei so we define them here
if(current_dir_vec[[2]] == "Volumes"){
  prsice_exe <- "PRSice_mac"
}else if(current_dir_vec[[2]] == "home"){
  prsice_exe <- "PRSice_linux"
}else {
  stop('Cannot recognize file structure. This program is set up for Linux and OSX.')
}

## ----------------------------- Calculate PRS -------------------------------

# Consider using mapply or map2 if we end up looping over both GWAS and target data
i <- 0
output <- list()
# Command line call to PRSice
for(gwas_name in gwas_names){
  i <- i + 1
  gwas_prefix <- str_split(gwas_name, pattern = "\\.")[[1]][1]
  output_file_name <- str_c(opt$geno_name, gwas_prefix, sep="_")
  cmd <- paste("Rscript", str_c(opt$prsice, "PRSice.R"),
               # point to directory to write PRS output
               "--dir", opt$out_dir,
               # point to PRSice executable
               "--prsice", str_c(opt$prsice, prsice_exe),
               # point to summary statistic to use (we are iterating over a list of these)
               "--base", str_c(opt$gwas_dir, gwas_name),
               # point to column names corresponding to each flag, this should be standardized to fit the below values
               "--snp", opt$snp, "--chr", opt$chr, "--bp", opt$bp, "--A1", opt$A1, "--A2", opt$A2, 
               "--stat", opt$stat, "--pvalue", opt$pvalue,
               # point to reference/target genome data, post-QC
               "--target", opt$geno,
               # --beta denotes continuous phenotype, --no-regress is necessary since we don't supply phenotype
               "--beta --no-regress T --out", output_file_name, 
               "--bar-levels 5e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1", 
               "--seed 112794 --fastscore --no-clump")
  
  # Run command line in terminal
  system(cmd)
  
  if(opt$combine) {
    output[[i]] <- cleanPRSiceOutput(str_c(output_file_name, ".all.score"), gwas_prefix, pval_thresholds)
  }
}


## ----------------------------- Combine PRS output -------------------------------
if(opt$combine){
  # Combine PRS from all phenotypes at selected p-value thresholds
  all_pheno <- output %>%
    reduce(full_join, by = c("FID", "IID")) 
  
  # Write combined output to output directory defined at beginning of code
  write_tsv(all_pheno, path = str_c(opt$out_dir, "/", today(), "_", opt$geno_name, "_Combined-PGS.txt"))  
}
