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
    cache: True
    wrapper:
        "0.80.2/bio/star/index"


rule star_align:
    input:
        unpack(get_fq),
        index= ("resources/star/" + str(get_source) + "genome" ),
        gtf=("resources/annotations/" + str(get_source) + "genome.gtf"),
    output:
        "{experiment}/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "{experiment}/star/{sample}-{unit}/ReadsPerGene.out.tab",
        "{experiment}/star/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
        "{experiment}/star/{sample}-{unit}/SJ.out.tab",
    log:
        "logs/star/{experiment}-{sample}-{unit}.log",
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
