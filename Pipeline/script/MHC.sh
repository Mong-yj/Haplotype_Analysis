#!/bin/bash
script_path="./Pipeline/script"

# set umask to avoid locking each other out of directories
umask 002

# get variables from workflow/variables.env
source $script_path/env/variables_MHC.env

# execute snakemake
snakemake \
    --profile $script_path/profiles/sge \
	--conda-frontend conda \
	--nolock \
    --snakefile ./smk/MHC.smk \
    #--singularity-args "--bind "

