# build model fitting scripts for each ebv ----------------------------------------------------
library(glue)

# data ----------------------------------------------------------------------------------------


full_M <- readRDS("R/dardel/full_M.rds")
model_script <- read_lines("R/model_script.txt") |> paste(collapse = "\n")

loc_list    <- names(full_M)
target_list <- names(full_M[[1]])


# Build each model fitting script -------------------------------------------------------------

# Loop over locations and EBVs 
for(i in loc_list){
  for(j in target_list){
  loc_var <- i # Location
  tar_var <- j # Target
  # Build model script
  mScript <-  glue("
      # libraries ----------------------------------------------------------------
        
      library(tidyverse)
      library(xgboost)
        
      # load feature & target matrices --------------------------------------------------------------
      ft_M <- readRDS('full_M.rds')${loc_var}${tar_var} |> drop_na() 
      
      ### Model fitting code ###  
      {model_script}
      #########################
      
      # save ----------------------------------------------------------------------------------------
        
       saveRDS(saveList , 'results/results_{loc_var}_{tar_var}.rds')
        
      ")
  # Save script
writeLines(mScript , glue("R/dardel/fit_model_{loc_var}_{tar_var}.R"))
  }
}


# Build shell scripts -------------------------------------------------------------------------


for(i in loc_list){
  for(j in target_list){
    loc_var <- i # Location
    tar_var <- j # Target
    
   shScript <-  glue("
#!/bin/bash -l
#SBATCH -A naiss2025-22-41
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 05:00:00
#SBATCH --mem 20000
#SBATCH -e logfiles/error_{loc_var}_{tar_var}.e
#SBATCH --mail-type=FAIL
#SBATCH -J fit_{loc_var}_{tar_var}

module load PDC/23.12    
module load R/4.4.0
Rscript fit_model_{loc_var}_{tar_var}.R > logfiles/out_{loc_var}_{tar_var}.out") 
    
   writeLines(shScript , glue("R/dardel/submit_model_{loc_var}_{tar_var}.sh"))
   
  }
  }

