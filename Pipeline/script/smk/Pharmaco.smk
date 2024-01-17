import glob

with open("Data/Sample_list.txt",'r') as f:
    SAMPLES = [sam.strip() for sam in f]

TARGETS = []
with open('./PGx_target','r') as f:
	for i in f:
		TARGETS.append(i.strip())

rule all:
	input:
		expand(["PGx/PharmCAT/{sample}/{sample}.match.json",
				"PGx/PharmCAT/{sample}/{sample}.phenotype.json",
				"PGx/PharmCAT/{sample}/{sample}.report.html",
				"PGx/Stargazer/{sample}/{target}/report.tsv",
				"PGx/Aldy/{sample}/{target}/result.tsv",
				"PGx/Summary/{sample}/Summary.txt",
				"PGx/Summary/All/Summary.2.log"],
				sample=SAMPLES, target=[t.lower() for t in TARGETS])

rule pharmcat:
	input:
		vcf = "Data/vcf/HaplotypeCaller_{sample}.vcf.gz",
		tbi = "Data/vcf/HaplotypeCaller_{sample}.vcf.gz.tbi",
	output:
		json1 = "PGx/PharmCAT/{sample}/{sample}.match.json",
		json2 = "PGx/PharmCAT/{sample}/{sample}.phenotype.json",
		html = "PGx/PharmCAT/{sample}/{sample}.report.html",
	log: "PGx/PharmCAT/{sample}/log/{sample}.log"
	benchmark: "PGx/PharmCAT/{sample}/benchmarks/{sample}.tsv"
	params: 
		base_vcf_name=lambda wildcards: f"HaplotypeCaller_{wildcards.sample}.vcf.gz",
		sam=lambda wildcards: f"{wildcards.sample}"
	singularity: f"docker://pgkb/pharmcat:latest",
	shell:
		"""
			(
			# 변수설정
			JAVA_HOME=""
			PARH=$PATH:/pharmcat
			TMP="/tmp"
			TMPDIR="/tmp"
			
			# run
			/pharmcat/pharmcat_pipeline \
			-o PGx/PharmCAT/{params.sam}/ \
			-bf {params.sam} \
			Data/vcf/{params.base_vcf_name}) > {log} 2>&1
		"""

rule stargazer:
	input:
		vcf = "Data/vcf/HaplotypeCaller_{sample}.vcf.gz",
		tbi = "Data/vcf/HaplotypeCaller_{sample}.vcf.gz.tbi",
	output:
		tsv = "PGx/Stargazer/{sample}/{target}/report.tsv",
	log: "PGx/Stargazer/{sample}/log/{sample}.{target}.log"
	benchmark: "PGx/Stargazer/{sample}/benchmarks/{sample}.{target}.tsv"
	params: 
		base_vcf_name=lambda wildcards: f"HaplotypeCaller_{wildcards.sample}.vcf.gz",
		sam=lambda wildcards: f"{wildcards.sample}",
		target=lambda wildcards: f"{wildcards.target}"
	singularity: f"docker://yuyae/stargazer:v2.0.2",
	shell:
		"""
			(python \
			/home/Tool/Stargazer/stargazer-grc38-v.2.0.2/stargazer \
			-i Data/vcf/{params.base_vcf_name} \
			-a grc38 \
			-o PGx/Stargazer/{params.sam}/{params.target} \
			-d wgs \
			-t {params.target}) > {log} 2>&1
		"""

rule aldy:
	input:
		bam = "Data/bam/{sample}.recal.bam",
		bamindex = "Data/bam/{sample}.recal.bam.bai",
	output:
		tsv = "PGx/Aldy/{sample}/{target}/result.tsv",
	log: "PGx/Aldy/{sample}/log/{sample}.{target}.log"
	benchmark: "PGx/Aldy/{sample}/benchmarks/{sample}.{target}.tsv"
	params: 
		base_bam_name=lambda wildcards: f"{wildcards.sample}.recal.bam",
		sam=lambda wildcards: f"{wildcards.sample}",
		target=lambda wildcards: f"{wildcards.target}",
	conda: "config/Haplotype.yaml",
	shell:
		"""
			(aldy genotype \
			-p wgs \
			-g {params.target} \
			Data/bam/{params.base_bam_name} \
			--log PGx/Aldy/{params.sam}/{params.target}/aldy.log \
			--output PGx/Aldy/{params.sam}/{params.target}/result.tsv \
			--genome hg38) > {log} 2>&1
		"""

rule Summary:
	input:
		pharmcat = "PGx/PharmCAT/{sample}/{sample}.match.json",
	output:
		summary = "PGx/Summary/{sample}/Summary.txt"
	log: "PGx/Summary/{sample}/log/Summary.1.log"
	params:
		sam=lambda wildcards: f"{wildcards.sample}"
	conda: "config/Haplotype.yaml"
	script:
		"Pharmaco_Result.py"

rule Summary_All:
	input:
		summary = expand(["PGx/Summary/{sample}/Summary.txt"], sample=SAMPLES),
	output:
		All_summary = "PGx/Summary/All/Summary.2.log"
	conda: "config/Haplotype.yaml"
	script:
		"Pharmaco_Summary.py"
