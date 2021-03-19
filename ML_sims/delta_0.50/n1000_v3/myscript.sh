#!/bin/bash
#!
#! Example SLURM job script
#!

#!#######################################
#!#### SLURM CONFIGURATION OPTIONS ######
#!#######################################

#! Name of the job: change this to anything you like
#! How many cores do you need in total?
#SBATCH -N 1
#SBATCH -c 20
#! How much wallclock time will be required? Use format DD-HH:MM:SS
#SBATCH -t 0-03:00  

#SBATCH -p short                        
#SBATCH --mem=50000                 
#SBATCH -o hostname_%j.out                
#SBATCH -e hostname_%j.err                
#SBATCH --mail-type=FAIL               
#SBATCH --mail-user=lanwen@hsph.harvard.edu  

#!#########################
#!#### COMMAND TO RUN #####
#!#########################

module load gcc/6.2.0 R/4.0.1
R CMD BATCH --quiet --no-restore --no-save dr_ss_param_haz.R 



