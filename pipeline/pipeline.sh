#!/bin/bash

# Job name:
#SBATCH --job-name=Pathogen_pipeline
#
# Project:
#SBATCH --account=nn9383k
#
# Wall clock limit:
#SBATCH --time=96:00:00
#

# Number of cores:
#SBATCH --mem-per-cpu=8G
## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors


## Do some work:
module load sdag

#run complete pipeline
sdag pipeline.sdag
