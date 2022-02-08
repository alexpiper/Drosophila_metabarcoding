#!/bin/bash
#SBATCH --job-name MVAbund
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=10G
#SBATCH --time=240:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.piper@agriculture.vic.gov.au


#Load modules and run R scripts
module load R/4.1.0-foss-2021a
Rscript run_mvabund.R