import pandas as pd
import sys

output_merge = snakemake.output.merge
output_summary = snakemake.output.Summary
input_tsvs = snakemake.input.tsv

all_data = []

for idx, tsv_file in enumerate(input_tsvs):
    df = pd.read_csv(tsv_file, sep='\t', header=0)
    if idx == 0:  # 첫 번째 파일에서 헤더를 가져오기
        header = df.columns.tolist()
    all_data.append(df[df[header[1]].str.contains('ALL')])

# 모든 데이터 병합
merged_data = pd.concat(all_data)
merged_data.to_csv(output_merge, index=False, sep='\t', float_format='%.2f')

# 통계 계산
summary_data = merged_data.iloc[:, 3:8] 
summary_stats = summary_data.describe().loc[['count','mean', 'min', '50%', 'max']].T
summary_stats.columns = ['count','mean', 'min', 'median', 'max']

# Summary 저장
with open(output_summary, 'w') as final_summary:
    final_summary.write('Index\tCount\tMean\tMin\tMedian\tMax\n')
    final_summary.write(summary_stats.to_csv(sep='\t', header=False, float_format='%.2f'))