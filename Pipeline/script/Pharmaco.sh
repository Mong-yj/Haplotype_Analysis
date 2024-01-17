#!/bin/bash
script_path="/ess/dlstibm/Workspace/workspace.ryj/Haplotype/Pipeline/script"

# set umask to avoid locking each other out of directories
umask 002

# get variables from workflow/variables.env
source $script_path/env/variables.Pharmaco.env

# execute snakemake
# PharmCAT
snakemake \
    --profile $script_path/profiles/sge \
	--conda-frontend conda \
	--nolock \
    --snakefile ./smk/Pharmaco.smk \
    #--singularity-args "--bind "
