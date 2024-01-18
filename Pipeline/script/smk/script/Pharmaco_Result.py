import json
import pandas as pd
import sys
import os

# os.chdir("/ess/dlstibm/Workspace/workspace.ryj/Haplotype/Step3.Apply")

p = snakemake.input.pharmcat
result = snakemake.output.summary
sample = p.split("/")[-2]

# target gene
with open("./Pipeline/Overlap_gene", "r") as f:
    target_gene = [i.strip() for i in f]

Result = pd.DataFrame(columns=["Sample", "Gene", "Stargazer", "Pharmcat", "Aldy"])

for i, gene in enumerate(target_gene):
    print(i)
    print(gene)
    # load output
    stargazer_info = pd.read_table(
        "PGx/Stargazer/%s/%s/report.tsv" % (sample, gene.lower()), sep="\t", header=0
    )
    with open("PGx/Aldy/%s/%s/result.tsv" % (sample, gene.lower()), "r") as f:
        aldy_info = [i.strip() for i in f]
    with open(p, "r") as file3:
        json_content = json.load(file3)

    # Extracting information
    pharmcat_info = [
        gene_info
        for gene_info in json_content["results"]
        if (gene_info["gene"] == gene.lower()) or (gene_info["gene"] == gene)
    ]
    diplotypes = pharmcat_info[0]["diplotypes"] if pharmcat_info else []
    diplotypes_name = diplotypes[0]["name"] if diplotypes else "NA"

    sample_report = stargazer_info["Sample"][0]
    gene_report = stargazer_info["Gene"][0].upper()
    diplotype_report = stargazer_info["Diplotype"][0]
    try:
        solution1_result = aldy_info[-1].split("\t")[3]
    except IndexError:
        solution1_result = "NA"

    # Creating a DataFrame using pandas
    df = pd.DataFrame(
        {
            "Sample": [sample_report],
            "Gene": [gene_report],
            "Stargazer": [diplotype_report],
            "Pharmcat": [diplotypes_name],
            "Aldy": [solution1_result],
        }
    )
    Result = pd.concat([Result, df], axis=0).reset_index(drop=True)

Result.to_csv(result, sep="\t")
