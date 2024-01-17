import glob

with open("Data/Sample_list.txt",'r') as f:
    SAMPLES = [sam.strip() for sam in f]

rule all:
	input:
		expand(["Haplotype/Yleaf/{sample}/hg_prediction.hg",
                "Haplotype/Yleaf/{sample}/{sample}.recal/{sample}.recal.chr",
                "Haplotype/Yleaf/{sample}/{sample}.recal/{sample}.recal.fmf",
                "Haplotype/Yleaf/{sample}/{sample}.recal/{sample}.recal.info",
                "Haplotype/Yleaf/{sample}/{sample}.recal/{sample}.recal.out",
                "Haplotype/Yleaf/{sample}/{sample}.recal/{sample}.recal.pmu",
                "Haplotype/Yleaf/Summary/Summary.txt"],
				sample=SAMPLES)

rule yleaf_run:
    input:
        bam="Data/bam/{sample}.recal.bam",
        bai="Data/bam/{sample}.recal.bam.bai",
    output: 
        hg = "Haplotype/Yleaf/{sample}/hg_prediction.hg",
        chr = "Haplotype/Yleaf/{sample}/{sample}.recal/{sample}.recal.chr",
        fmf = "Haplotype/Yleaf/{sample}/{sample}.recal/{sample}.recal.fmf",
        info = "Haplotype/Yleaf/{sample}/{sample}.recal/{sample}.recal.info",
        out = "Haplotype/Yleaf/{sample}/{sample}.recal/{sample}.recal.out",
        pmu = "Haplotype/Yleaf/{sample}/{sample}.recal/{sample}.recal.pmu",
    log: "Haplotype/Yleaf/log/{sample}/{sample}.smk.log"
    benchmark: "Haplotype/Yleaf/benchmarks/{sample}/{sample}.smk.benchmarks.tsv"
    params: 
        base_bam_name=lambda wildcards: f"{wildcards.sample}.recal.bam",
        sam=lambda wildcards: f"{wildcards.sample}"
    conda: "config/yleaf.yaml",
    shell:
        """
        (
        #git clone https://github.com/genid/Yleaf.git
        #cd /ess/dlstibm/Workspace/workspace.ryj/Haplotype/Step3.Apply/Yleaf/
        #pip install -U pip setuptools
        #pip install -e .
        #cd -
        
        Yleaf \
        -bam Data/bam/{params.base_bam_name} \
        --reference_genome hg38 \
        -p \
        -dh \
        -force \
        -o Haplotype/Yleaf/{params.sam}
        
        cp Haplotype/Yleaf/{params.sam}/hg_prediction.hg Haplotype/Yleaf/{params.sam}/{params.sam}.recal) > {log} 2>&1
        """



rule summary:
    input:
        hg = expand(["Haplotype/Yleaf/{sample}/hg_prediction.hg"], sample=SAMPLES),
    output: 
        Summary = "Haplotype/Yleaf/Summary/Summary.txt",
    log: "Haplotype/Yleaf/log/Summary.log"
    conda: "config/yleaf.yaml",
    shell:
        """
        cat {input.hg} | cut -f2 | sort | uniq -c | sort -nrk1 | grep -v Hg |\
        awk '{{print $2"\t"$1}}' | sed 1iAllele"\t"Count > {output.Summary}

        # cat Male.result | cut -f2 | sed 1d | awk -F '(' '{{print $1}}' | sed -e 's/*//g' |\
        # sort | uniq -c | sort -nrk1 | awk -F ' ' '{{print $2"\t"$1}}' | sed 1iAllele"\t"Count > Summary.male.txt
        """
