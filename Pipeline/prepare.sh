#!/bin/bash
if [ $# -ne 1 ]
then
        echo ""
        echo -e "       prepare.sh \$Tool_name ( whatshap / yleaf / MHC / Pharmaco )"
        echo -e ""
else
    tool=$1
    if [ ${tool} == "whatshap" ]
    then
        echo "You need bam and vcf input files"
    elif [ ${tool} == "yleaf" ]
    then 
        echo "You need bam input files"
    elif [ ${tool} == "MHC" ]
    then
        echo "You need fastq and bam input files"
    elif [ ${tool} == "Pharmaco" ]
    then
        echo "You need bam and vcf input files"
    else
        echo ${tool}" is an unavailable tool"
        echo "choose one of whatshap / yleaf / MHC / Pharmaco "
        exit 1
    fi
        
    pipeline_path="./Pipeline"
    script_path="./Pipeline/script"

    # source script copy
    mkdir -p smk/config/
    mkdir -p cluster_logs/
    cp $script_path/${tool}.sh ./
    cp $script_path/smk/${tool}.smk smk/
    cp $script_path/smk/config/Haplotype.yaml smk/config/

    if [ ${tool} != "Pharmaco" ]
    then
        cp $script_path/smk/config/${tool}.yaml smk/config/
    else
        cp $pipeline_path/PGx_target ./
    fi

    if [ ${tool} != "yleaf" ]
    then
        cp $script_path/smk/script/${tool}*.py smk/
    fi

fi

