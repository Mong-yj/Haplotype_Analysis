#!/bin/bash
#$ -A 100humans
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -pe smp 4
#$ -o ./cluster_logs/sge-$JOB_NAME-$JOB_ID-$HOSTNAME.out
#$ -e ./cluster_logs/sge-$JOB_NAME-$JOB_ID-$HOSTNAME.err

# USAGE: qsub workflow/process_smrtcells.sge.sh

# set umask to avoid locking each other out of directories
umask 002

# get variables from workflow/variables.env
source script/env/variables.env

# execute snakemake
snakemake \
    --profile script/profiles/sge \
	--conda-frontend conda \
	--nolock \
    --snakefile script/smk/haplogrep.smk

