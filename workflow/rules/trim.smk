rule sort_raw_reads:
    output:
        "reads/raw/{sample}_{unit}_{fq}_raw.fastq.gz",
    log:
        "logs/sort/{sample}_{unit}_{fq}_raw.log",
    resources:
        mem="6G",
        rmem="4G",
    params:
        input=lambda w: samples.loc[w.sample].loc[w.unit].loc[w.fq],
        gzipped=lambda w, input: 1 if samples.loc[w.sample].loc[w.unit].loc[w.fq].endswith(".gz") else 0,
        unzipped_out = "reads/raw/{sample}_{unit}_{fq}_raw.fastq",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 1
    shell:
        """
        if [[ -f "{params.input}" ]] ;
          then
            if [[ {params.gzipped} -eq 1 ]] ;
              then             
              ln -s $(readlink -f {params.input}) {output} 2> {log}
              else
              gzip -c {params.input} > {output} 2>{log}
            fi
          else
          if [[ {params.gzipped} -eq 1 ]] ;
            then
            wget {params.input} -O {output} 2> {log}
            else
            wget {params.input} -O {params.unzipped_out} 2> {log} &&
            gzip {params.unzipped_out} 2>>{log}
          fi
        fi
        """

rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "0.77.0/bio/sra-tools/fasterq-dump"

rule umi_extract_se:
    input:
        "reads/raw/{sample}_{unit}_fq1_raw.fastq.gz",
    output:
        fq="reads/umi_extracted/{sample}_{unit}_se_raw.fastq.gz",
    params:
        bc_pattern=lambda w: samples.loc[w.sample].loc[w.unit, "umi_bc"],
        no_umi=lambda w: 1 if str(samples.loc[w.sample].loc[w.unit, "umi_bc"])=="nan" else 0,
        begin_umi="workflow/scripts/awk/begin_umi.awk",
    conda:
        "../envs/umi_tools.yaml"
    resources:
         mem="8G",
         rmem="6G",
    log:
        "logs/umi/{sample}_{unit}_se_extract.log"   
    shell:
        """
        if [[ {params.no_umi} -eq 1 ]] ;
          then
          zcat {input} |
          awk -f {params.begin_umi} - |
          gzip > {output}
          else
          umi_tools extract \
          --stdin={input} \
          --bc-pattern={params.bc_pattern} \
          --log={log} \
          --stdout {output}
        fi  
        """ 
    
rule umi_extract_pe:
    input:
        fq1="reads/raw/{sample}_{unit}_fq1_raw.fastq.gz",
        fq2="reads/raw/{sample}_{unit}_fq2_raw.fastq.gz",
    output:
        fq1="reads/umi_extracted/{sample}_{unit}_fq1_raw.fastq.gz",
        fq2="reads/umi_extracted/{sample}_{unit}_fq2_raw.fastq.gz",
    params:
        bc_pattern=lambda w: samples.loc[w.sample].loc[w.unit, "umi_bc"],
        no_umi=lambda w: 1 if str(samples.loc[w.sample].loc[w.unit, "umi_bc"])=="nan" else 0,
        begin_umi="workflow/scripts/awk/begin_umi.awk",
    log:
        "logs/umi/{sample}_{unit}_pe_extract.log"
    conda:
        "../envs/umi_tools.yaml"
    resources:
         mem="8G",
         rmem="6G",
    shell:
        """
        if [[ {params.no_umi} -eq 1 ]] ;
          then
          zcat {input.fq1} |        
          awk -f {params.begin_umi} - |
          gzip > {output.fq1} &&
          zcat {input.fq2} |
          awk -f {params.begin_umi} - |
          gzip > {output.fq2}
          else  
          umi_tools extract \
          -I {input.fq1} \
          --bc-pattern={params.bc_pattern} \ 
          --read2-in={input.fq2} \
          --stdout={output.fq1} \
          --read2-out={output.fq2}
        fi
        """

rule cutadapt_pe:
    input:
        get_cutadapt_input,
    output:
        fastq1="reads/trimmed/{sample}_{unit}_fq1_trimmed.fastq.gz",
        fastq2="reads/trimmed/{sample}_{unit}_fq2_trimmed.fastq.gz",
        qc="qc/{sample}-{unit}.paired.qc.txt",
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    params:
        extra= ("--minimum-length " + str(config["trimming"]["minimum_read_length"]) + " -q " + str(config["trimming"]["minimum_read_quality"])) ,
        adapters=lambda w: str(samples.loc[w.sample].loc[w.unit].squeeze(axis=0)["adapters"]),
    threads: 4
    resources:
        mem="8G",
        rmem="6G",
    wrapper:
        "0.80.2/bio/cutadapt/pe"


rule cutadapt_se:
    input:
        get_cutadapt_input,
    output:
        fastq="reads/trimmed/{sample}_{unit}_single.fastq.gz",
        qc="qc/{sample}_{unit}_trimmed.qc.txt",
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    params:
        extra= ("--minimum-length " + str(config["trimming"]["minimum_read_length"]) + " -q " + str(config["trimming"]["minimum_read_quality"])) ,
        adapters_r1=lambda w: str(samples.loc[w.sample].loc[w.unit].squeeze(axis=0)[ "adapters"]),
    threads: 4
    resources:
        mem="16G",
        rmem="12G",
    wrapper:
        "0.80.2/bio/cutadapt/se"

