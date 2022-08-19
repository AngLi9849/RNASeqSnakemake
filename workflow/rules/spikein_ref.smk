
rule prefix_spikein_genome:
    input:
        fasta = lambda wildcards: "resources/genomes/{s}_{{species}}_genome.fasta".format(
            s=get_ref_source(wildcards.species)
        ),
    output:
        "resources/genomes/spikein_{species}_genome.fasta",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/ref/prefix_spikein_{species}_genome.log",
    shell:
        "sed 's/>/>spikein_/' {input} > {output} 2> {log}"

rule spikein_genome2bed:
    input:
        "resources/genomes/spikein_{source}_genome.fasta",
    output:
        "resources/genomes/spikein_{source}_genome.fasta.bed",
    params:
        fasta2bed_awk="workflow/scripts/awk/fasta2bed.awk"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/ref/spikein_{source}_genome2bed.log",
    shell:
        "awk -f {params.fasta2bed_awk} {input} > {output}"



rule prefix_spikein_annotation:
    input:
        lambda wildcards: "resources/annotations/{s}_{{species}}_genome.gtf".format(s=get_ref_source(wildcards.species)),
        awk_spikein_prefix="workflow/scripts/awk/spikein_prefix.awk"
    output:
        "resources/annotations/spikein_{species}_genome.gtf",
    cache: True
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/set_spikein_{species}_annotation.log",
    shell:
        "awk -F'\\t' -v OFS='\\t' -f {input.awk_spikein_prefix} {input[0]} > {output} 2> {log}"


rule spikein_faidx:
    input:
        "resources/genomes/spikein_{source}_genome.fasta",
    output:
        "resources/genomes/spikein_{source}_genome.fasta.fai",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/spikein_{source}_faidx.log",
    cache: True
    wrapper:
        "0.77.0/bio/samtools/faidx"


rule spikein_bwa_index:
    input:
        "resources/genomes/spikein_{source}_genome.fasta",
    output:
        multiext("resources/spikein/spikein_{source}_genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/spikein_{source}_bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "0.77.0/bio/bwa/index"


rule spikein_star_index:
    input:
        fasta="resources/genomes/spikein_{source}_genome.fasta",
        gtf="resources/annotations/spikein_{source}_genome.gtf",
    output:
        directory("resources/star/spikein_{source}_genome"),
    threads: 4
    params:
        extra="--genomeSAindexNbases 10"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/star_index_spikein_{source}_genome.log",
    cache: True
    wrapper:
        "0.80.2/bio/star/index"
