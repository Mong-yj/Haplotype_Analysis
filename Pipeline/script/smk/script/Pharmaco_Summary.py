import pandas as pd
import sys
import os

# os.chdir("/ess/dlstibm/Workspace/workspace.ryj/Haplotype/Step3.Apply")


def process_summary_file(file_path):
    df = pd.read_csv(file_path, delimiter="\t", index_col=0)

    # Function to determine the consensus or concatenate the values
    def get_consensus_value(row):
        values = [
            str(val).replace("*WT", "*1")
            for val in [row["Stargazer"], row["Pharmcat"], row["Aldy"]]
            if str(val) != "nan" and str(val) != "NA"
        ]
        if len(set(values)) == 1:
            return values[0]
        return "|".join(sorted(set(values)))

    df["Consensus"] = df.apply(get_consensus_value, axis=1)

    return df


def main():
    with open(
        "/ess/dlstibm/Workspace/workspace.ryj/Haplotype/Pipeline/Overlap_gene", "r"
    ) as f:
        target_gene = [i.strip() for i in f]

    result = {}
    for g in target_gene:
        result[g] = pd.DataFrame(index=[g])

    input_summary = snakemake.input.summary
    All_summary = snakemake.output.All_summary

    for Summary in input_summary:
        df = process_summary_file(Summary)
        for i, g in enumerate(target_gene):
            Genotype = df[df["Gene"] == g]["Consensus"].iloc[0]
            if Genotype in list(result[g].columns):
                result[g].loc[g, Genotype] = result[g].loc[g, Genotype] + 1
            else:
                result[g].loc[g, Genotype] = 1

    for df in result.keys():
        result[df] = result[df].astype(int)
        result[df] = result[df].T.sort_values(df, ascending=False)
        result[df].index.name = "Allele"
        result[df].to_csv("./PGx/Summary/All/%s.Summary" % df, sep="\t")

    with open(All_summary, "w") as f:
        f.write("Done!\n")


if __name__ == "__main__":
    main()
