import glob

input_bam = glob.glob("Data/bam/*.bam")
SAMPLES = [f.split('/')[-1].split('.')[0] for f in input_bam[:2]]
#SAMPLES = glob_wildcards("Data/test//{sample}.vcf.gz")


rule all:
	input:
		expand(["Haplotype/Haplogrep3/{sample}/{sample}.haplogrep",
                "Haplotype/Haplogrep3/{sample}/{sample}.fasta",
                "Haplotype/Haplogrep3/{sample}/{sample}.qc.txt",],
				sample=SAMPLES)


rule MT_call:
    shell:
        """
        java -jar /home/vcpa/pkg/bin/picard.jar CollectWgsMetrics \
        I=$bam \
        O=$prefix.collect_wgs_metrics.txt \
        R=/ess/dlstibm/Workspace/workspace.ryj/Haplotype/Step3.Apply/Data/Homo_sapiens_assembly38.fasta
        
        samtools view -bh $bam chrM > $prefix.chrM.bam
        samtools index $prefix.chrM.bam
        
        /home/vcpa/pkg/bin/gatk --java-options -Xmx16g Mutect2 \
        -R /ess/dlstibm/Workspace/workspace.ryj/Haplotype/Step3.Apply/Data/Homo_sapiens_assembly38.fasta \
        -I $prefix.chrM.bam \
        -L chrM \
        --mitochondria-mode \
        -O $prefix.chrM.vcf.gz

        """

rule Haplogrep3_run:
    input:
        vcf = "Haplotype/Haplogrep3/{sample}/MT/{sample}.chrM.vcf.gz",
        tbi = "Haplotype/Haplogrep3/{sample}/MT/{sample}.chrM.vcf.gz.tbi",
    output: 
        haplogrep = "Haplotype/Haplogrep3/{sample}/{sample}.haplogrep",
        fasta = "Haplotype/Haplogrep3/{sample}/{sample}.fasta",
        qc = "Haplotype/Haplogrep3/{sample}/{sample}.qc.txt",
    log: "Haplotype/Haplogrep3/{sample}/log/{sample}.MT.log"
    benchmark: "Haplotype/Haplogrep3/{sample}/benchmarks/{sample}.MT.benchmarks.tsv"
    params: 
        base_vcf_name=lambda wildcards: f"{wildcards.sample}.chrM.vcf.gz",
        sam=lambda wildcards: f"{wildcards.sample}"
    singularity: f"docker://yuyae/haplogrep3:latest",
    shell:
        """
        (haplogrep3 classify \
        --in Haplotype/Haplogrep3/{params.sam}/MT/{params.base_vcf_name} \
        --out {output.haplogrep} \
        --tree phylotree-rcrs@17.2 \
        --extend-report --write-qc --write-fasta) > {log} 2>&1
        """
