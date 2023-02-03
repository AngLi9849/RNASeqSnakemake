rule star_index:
    input:
        fasta="resources/genomes/{prefix}_genome.fasta",
        gtf="resources/annotations/{prefix}/genome.gtf",
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
        index=lambda wildcards:  ("resources/star/" + str(wildcards.reference)),
        gtf= lambda wildcards: ("resources/annotations/" + str(wildcards.reference) + "/genome.gtf"),
    output:
        "star/{sample}/{unit}/{reference}/Aligned.sortedByCoord.out.bam",
        "star/{sample}/{unit}/{reference}/ReadsPerGene.out.tab",
        "star/{sample}/{unit}/{reference}/Aligned.toTranscriptome.out.bam",
        "star/{sample}/{unit}/{reference}/SJ.out.tab",
        "star/{sample}/{unit}/{reference}/Log.final.out",
    log:
        "logs/star/{reference}-{sample}-{unit}.log",
    params:
        index=lambda wc, input: input.index,
        extra=lambda wc, input: "--outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax {max_multimap} --quantMode TranscriptomeSAM GeneCounts --outFilterType BySJout --outSAMattrRGline ID:{id} PU:{pu} SM:{sm} LB:unknown PL:{pl} --sjdbGTFfile {gtf} {config}".format(
           max_multimap = samples.loc[wc.sample].loc[wc.unit,"max_multimap"],
           id = wc.sample + "_" + wc.unit, 
           pu = wc.sample , 
           sm = wc.unit,
           pl = "illumina", 
           gtf = input.gtf, 
           config = config["params"]["star"],
        ),
    threads: 12
    resources:
        mem="16G",
        rmem="12G",
    wrapper:
        "0.80.2/bio/star/align"
