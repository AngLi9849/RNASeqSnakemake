rule star_index:
    input:
        fasta="resources/genomes/{prefix}.fasta",
        gtf="resources/annotations/{prefix}.gtf",
    output:
        directory("resources/star/{prefix}"),
    threads: 4
    resources:
        mem="24G",
        rmem="16G",
    params:
        extra="",
    log:
        "logs/{prefix}_star_index_genome.log",
    wrapper:
        "0.80.2/bio/star/index"


rule star_align:
    input:
        unpack(get_fq),
        index=lambda wildcards:  ("resources/star/" + str(wildcards.reference) + "_genome" ),
        gtf= lambda wildcards: ("resources/annotations/" + str(wildcards.reference) + "_genome.gtf"),
    output:
        "star/{sample}/{unit}/{reference}/Aligned.sortedByCoord.out.bam",
        "star/{sample}/{unit}/{reference}/ReadsPerGene.out.tab",
        "star/{sample}/{unit}/{reference}/Aligned.toTranscriptome.out.bam",
        "star/{sample}/{unit}/{reference}/SJ.out.tab",
    log:
        "logs/star/{reference}-{sample}-{unit}.log",
    params:
        index=lambda wc, input: input.index,
        extra=lambda wc, input: "--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterType BySJout --outSAMattrRGline ID:{} PU:{} SM:{} LB:unknown PL:{} --sjdbGTFfile {} {}".format(
           "{sample}_{unit}","{sample}_{unit}","{sample}_{unit}","illumina", input.gtf, config["params"]["star"],
        ),
    threads: 12
    resources:
        mem="16G",
        rmem="12G",
    wrapper:
        "0.80.2/bio/star/align"
