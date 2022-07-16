rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "0.77.0/bio/sra-tools/fasterq-dump"


rule cutadapt_pe:
    input:
        get_cutadapt_input,
    output:
        fastq1="{experiment}/trimmed_reads/{sample}_{unit}_fq1_trimmed.fastq.gz",
        fastq2="{experiment}/trimmed_reads/{sample}_{unit}_fq2_trimmed.fastq.gz",
        qc="{experiment}/qc/{sample}-{unit}.paired.qc.txt",
    log:
        "logs/cutadapt/{experiment}-{sample}-{unit}.log",
    params:
        extra= ("--minimum-length " + str(config["trimming"]["minimum_read_length"]) + " -q " + str(config["trimming"]["minimum_read_quality"])) ,
        adapters=lambda w: str(samples.loc[w.experiment].loc[w.sample].loc[w.unit].squeeze(axis=0)["adapters"]),
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
        fastq="{experiment}/trimmed_reads/{sample}_{unit}_trimmed.fastq.gz",
        qc="{experiment}/qc/{sample}_{unit}_trimmed.qc.txt",
    log:
        "logs/cutadapt/{experiment}-{sample}-{unit}.log",
    params:
        extra= ("--minimum-length " + str(config["trimming"]["minimum_read_length"]) + " -q " + str(config["trimming"]["minimum_read_quality"])) ,
        adapters_r1=lambda w: str(samples.loc[w.experiment].loc[w.sample].loc[w.unit].squeeze(axis=0)[ "adapters"]),
    threads: 4
    resources:
        mem="8G",
        rmem="6G",
    wrapper:
        "0.80.2/bio/cutadapt/se"
