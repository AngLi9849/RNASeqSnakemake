##FASTQC on input reads and RSEQC on STAR-aligned reads, generates a MULTIQC html report


rule raw_fastqc:
    input:
        get_raw_fastqc_input,
    output:
        html="qc/fastqc/{sample}-{unit}-{fq}-raw.html",
        zip="qc/fastqc/{sample}-{unit}-{fq}-raw_fastqc.zip"
    params: "--quiet"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/fastqc/{sample}_{unit}_{fq}_raw.log"
    threads: 1
    wrapper:
        "0.80.2/bio/fastqc"


rule trimmed_fastqc:
    input:
        "trimmed_reads/{sample}_{unit}_{fq}_trimmed.fastq.gz",
    output:
        html="qc/fastqc/{sample}-{unit}-{fq}-trimmed.html",
        zip="qc/fastqc/{sample}-{unit}-{fq}-trimmed_fastqc.zip"
    params: "--quiet"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/fastqc/{sample}_{unit}_{fq}.log"
    threads: 1
    wrapper:
        "0.80.2/bio/fastqc"


rule rseqc_gtf2bed:
    input:
        gtf= "resources/annotations/{reference}/genome.gtf",
    output:
        bed="resources/qc/rseqc/{reference}/annotation.bed",
        db=temp("resources/qc/rseqc/{reference}/annotation.db"),
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/{reference}_rseqc_gtf2bed.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/py/gtf2bed.py"


rule rseqc:
    input:
        bam="star/{sample}/{unit}/{reference}/AllAligned{prefix}.sortedByCoord.out.bam",
        bed="resources/qc/rseqc/{reference}/annotation.bed",
    output:
        junc_anno="qc/rseqc/{sample}.{unit}.{reference}.{prefix}.junctionanno.junction.bed",
        junc_sat="qc/rseqc/{sample}.{unit}.{reference}.{prefix}.junctionsat.junctionSaturation_plot.pdf",
        stats="qc/rseqc/{sample}.{unit}.{reference}.{prefix}.stats.txt",
        infer="qc/rseqc/{sample}.{unit}.{reference}.{prefix}.infer_experiment.txt",
        dist="qc/rseqc/{sample}.{unit}.{reference}.{prefix}.inner_distance_freq.inner_distance.txt",
        distribute="qc/rseqc/{sample}.{unit}.{reference}.{prefix}.readdistribution.txt",
        duprate="qc/rseqc/{sample}.{unit}.{reference}.{prefix}.readdup.DupRate_plot.pdf",
        gc="qc/rseqc/{sample}.{unit}.{reference}.{prefix}.readgc.GC_plot.pdf",
    priority: 1
    threads: 2
    resources:
        mem="48G",
        rmem="24G",
    log:
        "logs/rseqc/{sample}.{unit}.{prefix}.{reference}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        junction_prefix=lambda w, output: output.junc_anno.replace(".junction.bed", ""),
        junction_sat_prefix=lambda w, output: output.junc_sat.replace(".junctionSaturation_plot.pdf", ""),
        innerdis_prefix=lambda w, output: output.dist.replace(".inner_distance.txt", ""),
        readdup_prefix=lambda w, output: output.duprate.replace(".DupRate_plot.pdf", ""),
        readgc_prefix=lambda w, output: output.gc.replace(".GC_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.junction_prefix} > {log} 2>&1 ;
        junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.junction_sat_prefix} > {log} 2>&1 ; "
        bam_stat.py -i {input.bam} > {output[2]} 2> {log} ;
        infer_experiment.py -r {input.bed} -i {input.bam} > {output.infer} 2> {log} ;
        inner_distance.py -r {input.bed} -i {input.bam} -o {params.innerdis_prefix} > {log} 2>&1 ;
        read_distribution.py -r {input.bed} -i {input.bam} > {output.distribute} 2> {log} ;
        read_duplication.py -i {input.bam} -o {params.readdup_prefix} > {log} 2>&1 ;
        read_GC.py -i {input.bam} -o {params.readgc_prefix} > {log} 2>&1 
        """


rule multiqc:
    input:
        lambda wildcards: expand(
            "qc/fastqc/{sample.sample_name}-{sample.unit_name}-{sample.fq}-{trim}_fastqc.zip",
            sample=get_experiment_qc_samples(wildcards.experiment),
            trim=["raw","trimmed"],
        ),
        lambda wildcards: expand(
            "star/{sample.sample_name}-{sample.unit_name}/AllAligned.sortedByCoord.out.bam",
            sample=get_experiment_samples(wildcards.experiment),
        ),
#        expand(
#            "{sample.experiment}/qc/rseqc/{sample.sample_name}-{sample.unit_name}.{file}",
#            sample=samples.itertuples(),
#            file=[
#                "junctionanno.junction.bed",
#                "junctionsat.junctionSaturation_plot.pdf",
#                "infer_experiment.txt",
#                "stats.txt",
#                "inner_distance_freq.inner_distance.txt",
#                "readdistribution.txt",
#                "readdup.DupRate_plot.pdf",
#                "readgc.GC_plot.pdf",
#                ],
#        ),
#        expand(
#            "logs/rseqc/{sample.experiment}-{sample.sample_name}-{sample.unit_name}.log",
#            sample=samples.itertuples(),
#        ),
    output:
        "qc/multiqc/{experiment}.multiqc.html",
    threads: 1 
    resources:
        mem="36G",
        rmem="24G",
    log:
        "logs/{experiment}/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    script:
        "../scripts/py/multiqc.py"
