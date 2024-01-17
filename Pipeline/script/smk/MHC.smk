import glob

with open("Data/Sample_list.txt",'r') as f:
    SAMPLES = [sam.strip() for sam in f]


rule all:
	input:
		expand(["MHC/{sample}/Optitype/{sample}_result.tsv",
				"MHC/{sample}/Optitype/{sample}_coverage_plot.pdf",
				"MHC/{sample}/xHLA/_SUCCESS",
				"MHC/{sample}/xHLA/report-{sample}-hla.json",
				"MHC/{sample}/HLA-LA/{sample}/hla/R1_bestguess_G.txt",
				"MHC/{sample}/HLA-LA/{sample}/hla/R1_bestguess.txt",
				"MHC/{sample}/HLA-LA/{sample}/hla/histogram_matchesPerRead.txt",
				"MHC/{sample}/kourami/{sample}.result",
				"MHC/Summary/{sample}/Summary.txt",
				"MHC/Summary/All/Summary.2.log"],
				sample=SAMPLES)
				

rule optitype:
	input:
		fq1="Data/fastq/{sample}-x1_1.fastq.gz",
		fq2="Data/fastq/{sample}-x1_2.fastq.gz",
	output:
		tsv = "MHC/{sample}/Optitype/{sample}_result.tsv",
		pdf = "MHC/{sample}/Optitype/{sample}_coverage_plot.pdf",
	log: "MHC/{sample}/Optitype/log/{sample}.log"
	benchmark: "MHC/{sample}/Optitype/benchmarks/{sample}.tsv"
	params:
		fastq1=lambda wildcards: f"{wildcards.sample}-x1_1.fastq.gz",
		fastq2=lambda wildcards: f"{wildcards.sample}-x1_2.fastq.gz",
		sam=lambda wildcards: f"{wildcards.sample}"
	singularity: f"docker://fred2/optitype:latest",
	shell:
		"""
			(OptiTypePipeline.py \
			-i Data/fastq/{params.fastq1} -i Data/fastq/{params.fastq2} \
			-d -o "MHC/{params.sam}/Optitype/" \
			-p {params.sam}) > {log} 2>&1
		"""


rule xHLA:
	input:
		bam="Data/bam/{sample}.recal.bam",
		bai="Data/bam/{sample}.recal.bam.bai",
	output:
		sucess = "MHC/{sample}/xHLA/_SUCCESS",
		json = "MHC/{sample}/xHLA/report-{sample}-hla.json",
	log: "MHC/{sample}/xHLA/log/{sample}.log"
	benchmark: "MHC/{sample}/xHLA/benchmarks/{sample}.tsv"
	params:
		base_bam_name=lambda wildcards: f"{wildcards.sample}.recal.bam",
		sam=lambda wildcards: f"{wildcards.sample}"
	singularity: f"docker://humanlongevity/hla:latest",
	shell:
		"""
			(python /opt/bin/run.py \
			--sample_id {params.sam} \
			--input_bam_path Data/bam/{params.base_bam_name} \
			--output_path MHC/{params.sam}/xHLA/
			mv hla-{params.sam} MHC/{params.sam}/xHLA/) > {log} 2>&1
		"""

rule hlala: 
	input:
		bam="Data/bam/{sample}.recal.bam",
		bai="Data/bam/{sample}.recal.bam.bai",
	output:
		bestguess_G = "MHC/{sample}/HLA-LA/{sample}/hla/R1_bestguess_G.txt",
		bestguess = "MHC/{sample}/HLA-LA/{sample}/hla/R1_bestguess.txt",
		histogram = "MHC/{sample}/HLA-LA/{sample}/hla/histogram_matchesPerRead.txt"
	log: "MHC/{sample}/HLA-LA/{sample}/log/{sample}.log"
	benchmark: "MHC/{sample}/HLA-LA/{sample}/benchmarks/{sample}.tsv"
	params: 
		base_bam_name=lambda wildcards: f"{wildcards.sample}.recal.bam",
		sam=lambda wildcards: f"{wildcards.sample}"
	conda: "config/MHC.yaml"
	#singularity: f"docker://yuyae/hla-la:latest",
	shell:
		"""
		(~/anaconda3/envs/hla-la/bin/HLA-LA.pl \
		--BAM Data/bam/{params.base_bam_name} \
		--graph PRG_MHC_GRCh38_withIMGT \
		--sampleID {params.sam} \
		--workingDir MHC/{params.sam}/HLA-LA/ \
		--maxThreads 8
		) > {log} 2>&1
		"""

rule kourami_step1: 
	input:
		bam="Data/bam/{sample}.recal.bam",
		bai="Data/bam/{sample}.recal.bam.bai",
	output:
		bam = "MHC/{sample}/kourami/{sample}_on_KouramiPanel.bam",
		fq1 = "MHC/{sample}/kourami/{sample}_extract_1.fq.gz",
		fq2 = "MHC/{sample}/kourami/{sample}_extract_2.fq.gz",
	log: "MHC/{sample}/kourami/log/{sample}.1.log"
	benchmark: "MHC/{sample}/kourami/benchmarks/{sample}.1.tsv"
	params: 
		base_bam_name=lambda wildcards: f"{wildcards.sample}.recal.bam",
		sam=lambda wildcards: f"{wildcards.sample}"
	singularity: f"docker://yuyae/kourami:latest",
	shell:
		"""
		(/kourami/scripts/alignAndExtract_hs38DH.sh \
		{params.sam} \
		Data/bam/{params.base_bam_name}
		
		mv {params.sam}_extract*.fq.gz MHC/{params.sam}/kourami/
		mv {params.sam}_on_KouramiPanel.bam MHC/{params.sam}/kourami/) > {log} 2>&1
		"""

rule kourami_step2: 
	input:
		bam = "MHC/{sample}/kourami/{sample}_on_KouramiPanel.bam",
	output:
		result = "MHC/{sample}/kourami/{sample}.result",
	log: "MHC/{sample}/kourami/log/{sample}.2.log"
	benchmark: "MHC/{sample}/kourami/benchmarks/{sample}.2.tsv"
	params: 
		base_bam_name=lambda wildcards: f"{wildcards.sample}_on_KouramiPanel.bam",
		sam=lambda wildcards: f"{wildcards.sample}"
	singularity: f"docker://yuyae/kourami:latest",
	shell:
		"""
		(java -jar /kourami/target/Kourami.jar \
		-d /kourami/db/ \
		-o MHC/{params.sam}/kourami/{params.sam} \
		MHC/{params.sam}/kourami/{params.base_bam_name}) > {log} 2>&1
		"""


rule Summary:
	input:
		Optitype = "MHC/{sample}/Optitype/{sample}_result.tsv",
		xHLA = "MHC/{sample}/xHLA/report-{sample}-hla.json",
		HLA_LA = "MHC/{sample}/HLA-LA/{sample}/hla/R1_bestguess_G.txt",
		kourami = "MHC/{sample}/kourami/{sample}.result",
	output:
		summary = "MHC/Summary/{sample}/Summary.txt",
	log: "MHC/Summary/{sample}/Summary.1.log"
	params: 
		sample=lambda wildcards: f"{wildcards.sample}"
	conda: "config/Haplotype.yaml"
	script:
		"MHC_Result.py"


rule Summary_All:
	input:
		summary = expand(["MHC/Summary/{sample}/Summary.txt"], sample=SAMPLES),
	output:
		All_summary = "MHC/Summary/All/Summary.2.log"
	conda: "config/Haplotype.yaml"
	script:
		"MHC_Summary.py"