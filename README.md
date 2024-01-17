# Haplotype_Analysis

## 실행 방법

1. 분석을 진행하고자 하는 위치로 이동
2. Data 폴더 생성하여 **input** 으로 사용 될 bam, vcf, fastq 파일 등을 복사 
    
    & **Sample_list.txt** 파일 생성하여 진행 할 Sample ID 기입 
    
    - Data 폴더 내 구조 예시
        
        ![Untitled](https://prod-files-secure.s3.us-west-2.amazonaws.com/7ae9e35a-94e9-410c-a0b9-512409c3cf25/204f6aec-7d2a-49c5-8980-e038664e6274/Untitled.png)
        
3. prepare 코드 실행
    
    ```bash
    sh /ess/dlstibm/Workspace/workspace.ryj/Haplotype/Pipeline/prepare.sh [**whatshap / yleaf / MHC / Pharmaco**] 
    ```
    
    → **whatshap / yleaf / MHC / Pharmaco** 중 진행하고자 하는 Tool 이름 넣어 해당 스크립트 실행
    
    **만약, Data 폴더 내 input 파일이 symbolic link 라면 main script 에 원본 파일 위치를 bind 해줘야 함**
    
    ex) 위 스크립트 실행 후 복사 된 "whatshap.sh" 파일 내에 bind 주소 넣어주기
    
    ```bash
    # 원본
    snakemake \
        --profile $script_path/profiles/sge \
        --conda-frontend conda \
        --nolock \
    		--snakefile ./smk/whatshap.smk \
    
    # 수정
    snakemake \
        --profile $script_path/profiles/sge \
        --conda-frontend conda \
        --nolock \
    		--snakefile ./smk/whatshap.smk \
    		**--singularity-args** "--bind /ess/dlstibm/PROJECT/KDNA-DATA/submit/ --bind /ess/dlstibm/Analysis/20210915_KRIBB_PJH_WGS_1/results/Preprocessing/"
    ```
    
4. conda 환경 실행 후 run


