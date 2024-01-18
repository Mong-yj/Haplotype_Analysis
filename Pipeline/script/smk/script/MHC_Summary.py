import json
import pandas as pd
import numpy as np
import sys
import os
import re


def main():
    input_summary = snakemake.input.summary
    All_summary = snakemake.output.All_summary

    with open("./Pipeline/MHC.allele", "r") as f:
        MHC_allele = [i.strip() for i in f]

    result = {}
    for allele in MHC_allele:
        result[allele] = pd.DataFrame(index=[allele])

    for Summary in input_summary:
        df = pd.read_table(Summary, sep="\t", header=0, index_col=0).fillna("NA")
        for i, allele in enumerate(MHC_allele):
            Genotype = df[df.index.str.startswith(allele)]["Allele"].tolist()
            for G in Genotype:
                if G == "NA":
                    continue
                G_re = ":".join(G.split(":")[:2])
                if G_re in list(result[allele].columns):
                    result[allele].loc[allele, G_re] += 1
                else:
                    result[allele].loc[allele, G_re] = 1

    for df in result.keys():
        result[df] = result[df].astype(int)
        result[df] = result[df].T.sort_values(df, ascending=False)
        result[df].index.name = "Allele"
        result[df].to_csv("./MHC/Summary/All/%s.Summary" % df, sep="\t")

    with open(All_summary, "w") as f:
        f.write("Done!\n")


if __name__ == "__main__":
    main()
