rule sort_raw_reads:
    output:
        "reads/raw/{sample}_{unit}_{fq}_raw.fastq.gz",
    log:
        "logs/sort/{sample}_{unit}_{fq}_raw.log",
    resources:
        mem="6G",
        rmem="4G",
    params:
        gzipped=lambda w, input: 1 if input[0].endswith(".gz") else 0,
        unzipped_out = "reads/raw/{sample}_{unit}_{fq}_raw.fastq",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 1
    shell:
        """
        if [[ -f "{input}" ]] ;
          then
            if [[ {params.gzipped} -eq 1 ]] ;
              then             
              ln -s $(readlink -f {input}) {output} 2> {log}
              else
              gzip -c {input} > {output} 2>{log}
            fi
          else
          if [[ {params.gzipped} -eq 1 ]] ;
            then
            wget {input} -O {output} 2> {log}
            else
            wget {input} -O {params.unzipped_out} 2> {log} &&
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
       fq="reads/umi_extracted/{sample}_{unit}_raw.fastq.gz",
    params:
       bc_pattern=lambda w: samples.loc[w.sample].loc[w.unit, "umi_bc"],
    log:
       "logs/sort/{sample}_{unit}_{fq}_raw.log"   
    shell:
       """
       umi_tools extract \
       --stdin={input} \
       --bc-pattern={params.bc_pattern} \
       --log={log} \
       --stdout {output}  
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
    log:
       "logs/sort/{sample}_{unit}_{fq}_raw.log"
    shell:
       """
       umi_tools extract \
       -I {input.fq1} \
       --bc-pattern={params.bc_pattern} \ 
       --read2-in={input.fq2} \
       --stdout={output.fq1} \
       --read2-out={output.fq2}
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
        fastq="reads/trimmed/{sample}_{unit}_trimmed.fastq.gz",
        qc="qc/{sample}_{unit}_trimmed.qc.txt",
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    params:
        extra= ("--minimum-length " + str(config["trimming"]["minimum_read_length"]) + " -q " + str(config["trimming"]["minimum_read_quality"])) ,
        adapters_r1=lambda w: str(samples.loc[w.sample].loc[w.unit].squeeze(axis=0)[ "adapters"]),
    threads: 4
    resources:
        mem="8G",
        rmem="6G",
    wrapper:
        "0.80.2/bio/cutadapt/se"

