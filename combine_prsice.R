##################################
######
###### Combine PRSice output
###### Author: Jeremy Vollen
###### Date: 18.12.2019
######
###### Combines individual PRSice .all.score outputs into a merged,
###### fully-labeled dataset
###### 
##################################

library(tidyverse)
library(lubridate)

## Combine PRSice outputs (individual for each phenotype) into a combined data frame,
## allowing selection of subset of p-values and automatically selecting all "all.score" files
## 
## Args:
##  dir - directory where PGS are located and output is written 
##    (default: working directory)
##  file_list - list of file names (relative) with PRSice output 
##    (default: all 'all.score' suffixes in dir provided)
##  thresholds - p-value columns to select from each PRSice output
##  
## Usage: combinePRSiceOutput(dir = "/home/ubuntu/PGS_dir/", file_list = c("bmi.txt", "ea.txt"), 
##                            thresholds = c("5e-8", "1"), output_prefix = "ELSA")
combinePRSiceOutput <- function(dir = "./", 
                                file_list = list.files(path = dir, pattern = "all.score"), 
                                thresholds, output_prefix){
  # extract phenotype/trait string for each from provided file list
  phenotypes <- map(file_list, 
                    function(x) 
                      tail(str_split(
                        str_split(x, pattern = "\\.")[[1]][1],
                        pattern = "_")[[1]], n=1))
  
  # read and clean PRSice output for each file specified by user
  all_outputs <- map2(file_list, phenotypes, 
                      ~ cleanPRSiceOutput(str_c(dir, .x), .y, thresholds)) 
  
  # Combine PRS from all phenotypes at selected p-value thresholds
  combined <- all_outputs %>%
    reduce(full_join, by = c("FID", "IID"))
  
  # Write combined output to output directory defined at beginning of code
  write_tsv(combined, path = str_c(dir, today(), "_", output_prefix, "_Combined-PGS.txt")) 
  
  return(combined)
}


## Given a file path to PRSice output, select only the columns we care about and rename them
cleanPRSiceOutput <- function(file, gwas_prefix, thresholds){
  # Read table outputted by PRSice, path provided by user
  prs_output <- read_table2(file)
  
  # Select only the PRS for the p-value thresholds we care about
  prs_clean <- prs_output %>%
    select(FID, IID, one_of(thresholds))
  
  # For p-value columns, change name to be of form {prefix}_pval_{threshold}
  names(prs_clean) <- map(names(prs_clean), 
                          function(x) 
                            if(x %in% c("FID", "IID")){x}
                          else{str_c(gwas_prefix, 
                                     str_replace_all(x, 
                                                     pattern = "-0",
                                                     replacement = ""),
                                     sep = "_")})
  
  return(prs_clean)
}
