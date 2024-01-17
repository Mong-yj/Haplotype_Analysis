import glob

with open("Data/Sample_list.txt",'r') as f:
    SAMPLES = [sam.strip() for sam in f]

rule all:
	input:
		expand(["Haplotype/whatshap/{sample}/{sample}.phased.vcf",
				"Haplotype/whatshap/{sample}/{sample}.phased.gtf",
                "Haplotype/whatshap/{sample}/{sample}.phased.tsv",
                "Haplotype/whatshap/{sample}/{sample}.phased.blocklist",
                "Haplotype/whatshap/Summary/All.Summary"],
				sample=SAMPLES)

rule whatshap_phase:
    input:
        reference = "/ess/dlstibm/Workspace/workspace.ryj/Haplotype/Data/Homo_sapiens_assembly38.fasta",
        vcf = "Data/vcf/HaplotypeCaller_{sample}.vcf.gz",
        tbi = "Data/vcf/HaplotypeCaller_{sample}.vcf.gz.tbi",
        bam = "Data/bam/{sample}.recal.bam",
        bamindex = "Data/bam/{sample}.recal.bam.bai",
    output: "Haplotype/whatshap/{sample}/{sample}.phased.vcf"
    log: "Haplotype/whatshap/{sample}/log/{sample}.phase.log"
    benchmark: "Haplotype/whatshap/{sample}/benchmarks/{sample}.phase.tsv"
    params: 
        base_vcf_name=lambda wildcards: f"HaplotypeCaller_{wildcards.sample}.vcf.gz",
        base_bam_name=lambda wildcards: f"{wildcards.sample}.recal.bam"
    conda: "config/whatshap.yaml"
    resources:
        load=100
    shell:
        """
        (whatshap phase \
		 	--output {output} \
            --reference {input.reference} \
            Data/vcf/{params.base_vcf_name} \
            Data/bam/{params.base_bam_name}) > {log} 2>&1
        """

rule whatshap_stats:
    input:
        vcf = "Haplotype/whatshap/{sample}/{sample}.phased.vcf",
    output:
        gtf = "Haplotype/whatshap/{sample}/{sample}.phased.gtf",
        tsv = "Haplotype/whatshap/{sample}/{sample}.phased.tsv",
        blocklist = "Haplotype/whatshap/{sample}/{sample}.phased.blocklist"
    log: "Haplotype/whatshap/{sample}/log/{sample}.stats.log"
    benchmark: "Haplotype/whatshap/{sample}/benchmarks/{sample}.stats.tsv"
    params: 
        vcf_name=lambda wildcards: f"{wildcards.sample}/{wildcards.sample}.phased.vcf"
    conda: "config/whatshap.yaml"
    shell:
        """
        (whatshap stats \
            --gtf {output.gtf} \
            --tsv {output.tsv} \
            --block-list {output.blocklist} \
            Haplotype/whatshap/{params.vcf_name}) > {log} 2>&1
        """



rule whatshap_Summary:
    input:
        tsv = expand(["Haplotype/whatshap/{sample}/{sample}.phased.tsv"], sample=SAMPLES),
    output:
        merge = "Haplotype/whatshap/Summary/All.merge",
        Summary = "Haplotype/whatshap/Summary/All.Summary",
    log: "Haplotype/whatshap/Summary/Summary.log"
    conda: "config/Haplotype.yaml"
    script:
        "whatshap_Result.py"