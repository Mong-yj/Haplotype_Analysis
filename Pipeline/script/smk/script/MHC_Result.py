import json
import pandas as pd
import numpy as np
import sys
import os
import re


def set_index(input_list):
    counts = {}
    result = []
    for item in input_list:
        if item in counts:
            counts[item] += 1
            result.append(f"{item}.{counts[item]}")
        else:
            counts[item] = 1
            result.append(f"{item}.1")
    return result


def most_frequent(lst):
    lst_head = lst.apply(
        lambda x: ":".join(x.split(":")[0:2]) if x is not np.NaN else x
    )
    if len(list(lst)) == 1:
        return lst.value_counts()
    else:
        return lst_head.value_counts()


def select_allele(row):
    # If Optitype is not NaN, return Optitype, otherwise return xHLA
    if pd.notna(row["Optitype"]):
        return row["Optitype"]
    elif pd.notna(row["xHLA"]):
        return row["xHLA"]
    elif pd.notna(row["HLA_LA"]):
        return row["HLA_LA"]
    else:
        return row["kourami"]


def select_allele_with_most_frequent(row):
    # Use most_frequent function to find the most common allele, ignoring NaN
    most_common = most_frequent(row.dropna()).idxmax()
    most_count = most_frequent(row).max()
    unique_values = list(row.dropna().unique())

    if most_count != 1 or len(unique_values) == 1:
        return most_common
    elif len(unique_values) > 1:
        return select_allele(row)
    else:
        return np.NaN


def main():
    # os.chdir("/ess/dlstibm/Workspace/workspace.ryj/Haplotype/Step3.Apply")

    kourami = snakemake.input.kourami  # MHC/{sample}/kourami/{sample}.result
    HLA_LA = (
        snakemake.input.HLA_LA
    )  # MHC/{sample}/HLA-LA/{sample}/hla/R1_bestguess_G.txt
    Optitype = snakemake.input.Optitype  # MHC/{sample}/Optitype/{sample}_result.tsv
    xHLA = snakemake.input.xHLA  # MHC/{sample}/xHLA/report-{sample}-hla.json
    Result = snakemake.output.summary

    # Reading the contents of each file
    kourami_content = pd.DataFrame(pd.read_table(kourami, sep="\t", header=None)[0])
    kourami_content.columns = ["kourami"]
    kourami_content = kourami_content.sort_values("kourami")
    kourami_content.index = set_index(
        [i[0] for i in kourami_content["kourami"].str.split("*")]
    )

    HLA_LA_df = pd.read_table(HLA_LA, sep="\t", header=0)
    HLA_LA_content = pd.DataFrame(HLA_LA_df["Allele"].tolist(), columns=["HLA_LA"])
    HLA_LA_content = HLA_LA_content.sort_values("HLA_LA")
    HLA_LA_content.index = set_index(
        [i[0] for i in HLA_LA_content["HLA_LA"].str.split("*")]
    )

    Optitype_content = pd.read_table(Optitype, sep="\t", header=0, index_col=0).T.iloc[
        :-2, :
    ]
    Optitype_content.columns = ["Optitype"]
    Optitype_content = Optitype_content.sort_values("Optitype")
    Optitype_content.index = set_index(
        [i[0] for i in Optitype_content["Optitype"].str.split("*")]
    )

    with open(xHLA, "r") as f:
        xHLA_json = json.load(f)
    xHLA_content = pd.DataFrame(xHLA_json["hla"]["alleles"], columns=["xHLA"])
    xHLA_content = xHLA_content.sort_values("xHLA")
    xHLA_content.index = set_index([i[0] for i in xHLA_content["xHLA"].str.split("*")])

    # Result joining
    locus = HLA_LA_df["Locus"].tolist()
    Result_df = pd.DataFrame(index=set_index(locus))
    Result_df = pd.concat(
        [Result_df, kourami_content, HLA_LA_content, Optitype_content, xHLA_content],
        axis=1,
    )

    # Select Allele
    Result_df["Allele"] = Result_df.apply(select_allele_with_most_frequent, axis=1)
    Result_df = Result_df.fillna("NA")
    Result_df.to_csv(Result, sep="\t")


if __name__ == "__main__":
    main()
