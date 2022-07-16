##FASTQC on input reads and RSEQC on STAR-aligned reads, generates a MULTIQC html report

rule sort_raw_reads:
    input:
        sort_raw_reads,
    output:
        "{experiment}/raw_reads/{sample}_{unit}_{fq}_raw.{ext}",
    log:
        "logs/sort/{experiment}_{sample}_{unit}_{fq}_raw.{ext}.log",
    resources:
        mem="6G",
        rmem="4G",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 1
    shell:
        "cat {input} > {output} 2> {log}"


rule raw_fastqc:
    input:
        get_raw_fastqc_input,
    output:
        html="{experiment}/qc/fastqc/{sample}-{unit}-{fq}-raw.html",
        zip="{experiment}/qc/fastqc/{sample}-{unit}-{fq}-raw_fastqc.zip"
    params: "--quiet"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/fastqc/{experiment}_{sample}_{unit}_{fq}_raw.log"
    threads: 1
    wrapper:
        "0.80.2/bio/fastqc"


rule trimmed_fastqc:
    input:
        "{experiment}/trimmed_reads/{sample}_{unit}_{fq}_trimmed.fastq.gz",
    output:
        html="{experiment}/qc/fastqc/{sample}-{unit}-{fq}-trimmed.html",
        zip="{experiment}/qc/fastqc/{sample}-{unit}-{fq}-trimmed_fastqc.zip"
    params: "--quiet"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/fastqc/{experiment}_{sample}_{unit}_{fq}.log"
    threads: 1
    wrapper:
        "0.80.2/bio/fastqc"


rule rseqc_gtf2bed:
    input:
        gtf= lambda wildcards: "resources/annotations/{source}genome.gtf".format(source=get_source(wildcards)),
    output:
        bed="{experiment}/qc/rseqc/annotation.bed",
        db=temp("{experiment}/qc/rseqc/annotation.db"),
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/{experiment}_rseqc_gtf2bed.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/py/gtf2bed.py"


rule rseqc:
    input:
        bam="{experiment}/star/{sample}-{unit}/AllAligned.sortedByCoord.out.bam",
        bed="{experiment}/qc/rseqc/annotation.bed",
    output:
        "{experiment}/qc/rseqc/{sample}-{unit}.junctionanno.junction.bed",
        "{experiment}/qc/rseqc/{sample}-{unit}.junctionsat.junctionSaturation_plot.pdf",
        "{experiment}/qc/rseqc/{sample}-{unit}.stats.txt",
        "{experiment}/qc/rseqc/{sample}-{unit}.infer_experiment.txt",
        "{experiment}/qc/rseqc/{sample}-{unit}.inner_distance_freq.inner_distance.txt",
        "{experiment}/qc/rseqc/{sample}-{unit}.readdistribution.txt",
        "{experiment}/qc/rseqc/{sample}-{unit}.readdup.DupRate_plot.pdf",
        "{experiment}/qc/rseqc/{sample}-{unit}.readgc.GC_plot.pdf",
    priority: 1
    threads: 2
    resources:
        mem="20G",
        rmem="16G",
    log:
        "logs/rseqc/{experiment}-{sample}-{unit}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        junction_prefix=lambda w, output: output[0].replace(".junction.bed", ""),
        junction_sat_prefix=lambda w, output: output[1].replace(".junctionSaturation_plot.pdf", ""),
        innerdis_prefix=lambda w, output: output[4].replace(".inner_distance.txt", ""),
        readdup_prefix=lambda w, output: output[6].replace(".DupRate_plot.pdf", ""),
        readgc_prefix=lambda w, output: output[7].replace(".GC_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.junction_prefix} "
        "> {log} 2>&1 ; "
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.junction_sat_prefix} "
        "> {log} 2>&1 ; "
        "bam_stat.py -i {input.bam} > {output[2]} 2> {log} ; "
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output[3]} 2> {log} ; "
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.innerdis_prefix} > {log} 2>&1 ; "
        "read_distribution.py -r {input.bed} -i {input.bam} > {output[5]} 2> {log} ; "
        "read_duplication.py -i {input.bam} -o {params.readdup_prefix} > {log} 2>&1 ; "
        "read_GC.py -i {input.bam} -o {params.readgc_prefix} > {log} 2>&1"


rule multiqc:
    input:
       expand(
            "{sample.experiment}/qc/fastqc/{sample.sample_name}-{sample.unit_name}-{fq}-{trim}_fastqc.zip",
            sample=samples.itertuples(),
            fq=["fq1","fq2"],
            trim=["raw","trimmed"],
        ),
        expand(
            "{sample.experiment}/star/{sample.sample_name}-{sample.unit_name}/Aligned.sortedByCoord.out.bam",
            sample=samples.itertuples(),
        ),
        expand(
            "{sample.experiment}/qc/rseqc/{sample.sample_name}-{sample.unit_name}.{file}",
            sample=samples.itertuples(),
            file=[
                "junctionanno.junction.bed",
                "junctionsat.junctionSaturation_plot.pdf",
                "infer_experiment.txt",
                "stats.txt",
                "inner_distance_freq.inner_distance.txt",
                "readdistribution.txt",
                "readdup.DupRate_plot.pdf",
                "readgc.GC_plot.pdf",
                ],
        ),
        expand(
            "logs/rseqc/{sample.experiment}-{sample.sample_name}-{sample.unit_name}.log",
            sample=samples.itertuples(),
        ),
    output:
        "results/qc/multiqc_report.html",
    threads: 1 
    resources:
        mem="20G",
        rmem="12G",
    log:
        "logs/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    script:
        "../scripts/py/multiqc.py"
