#!/bin/bash

#SBATCH -c 2                  		# number of cores
#SBATCH --mem=12GB              	# specify additional memory to accomodate bbmap
#SBATCH -t 2-00:00:00			# wall time (D-HH:MM)
#SBATCH -p public			# Uses the htc partition
#SBATCH -q public			# Submits job queue
#SBATCH -o scinkd.%j.out       	 	# STDOUT (%j = JobId)
#SBATCH -e scinkd.%j.err       	 	# STDERR (%j = JobId)
#SBATCH --mail-type=ALL             # Send a notification when a job starts, stops, or fails
#SBATCH --mail-user=bjp004@morningside.edu # send-to address
#SBATCH --export=NONE

module purge

source activate scimr

time snakemake --use-conda --rerun-incomplete --nolock --cores 2 -s scimr_v0.2.py
