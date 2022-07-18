rule mv_bam:
    input:
        "star/{sample}/{unit}/{reference}/Aligned.sortedByCoord.out.bam",
    output:
        "star/{sample}/{unit}/{reference}/AllAligned.sortedByCoord.out.bam",
    resources:
        mem="8G",
        rmem="6G",
    shell:
        """
        mv {input} {output}
        """

rule samtools_read_count:
    input:
        bam="star/{sample}/{unit}/{reference}/{prefix}.sortedByCoord.out.bam",
        bai="star/{sample}/{unit}/{reference}/{prefix}.sortedByCoord.out.bam.bai",
    output:
        tab="star/{sample}/{unit}/{reference}/{prefix}.total_read_count.tab",
    threads: 1 
    resources:
        mem="8G",
        rmem="6G",
    log:
        "logs/samtools/{sample}/{unit}/{reference}/{prefix}_read_count.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools idxstat {input.bam} |
        awk -F'\\t' -v OFS='\\t' '
          BEGIN {{
            sum=0 ; spikein=0
          }}
          $1 !~ "spikein_" {{
            sum += $3
          }}
          $1 ~ "spikein_" {{
            spikein += $3
          }}
          END {{
            print "{wildcards.sample}", sum, spikein >> "{output.tab}"
          }} ' -
        """

rule samtools_seperate_splice:
    input:
        "star/{sample}/{unit}/{reference}/AllAligned.sortedByCoord.out.bam",
    output:
        "star/{sample}/{unit}/{reference}/UnsplicedAligned.sortedByCoord.out.bam",
        "star/{sample}/{unit}/{reference}/SplicedAligned.sortedByCoord.out.bam",
    params:
        max=config["unspliced_max_length"],
        min=config["unspliced_min_length"]
    threads: 4
    resources:
        mem="8G",
        rmem="6G",
    log:
        "logs/samtools/{sample}/{unit}/{reference}/separate-splice.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view --threads 3 -h {input}   |
        awk '$0 ~ /^@/ || ($6 !~ /N/ && $9 >= {params.min} && $9 <= {params.max}) || ($6 !~ /N/ && $9 >= (0-{params.max}) && $9 <= (0-{params.min}))'  |
        samtools view --threads 3 -bo {output[0]} &&
        samtools view --threads 3 -h {input}   |
        awk '$0 ~ /^@/ || ($6 ~ /N/) || ($6 !~ /N/ && $9 >= {params.max}) || ($6 !~ /N/ && $9 <= (0-{params.max}))' |
        samtools view --threads 3 -bo {output[1]}
        """

rule samtools_demultimap:
    input:
        "star/{sample}/{unit}/{reference}/{prefix}.sortedByCoord.out.bam",
    output:
        "star/{sample}/{unit}/{reference}/{prefix}Demultimapped.sortedByCoord.out.bam",
    params:
        extra="-b -F 256"
    threads: 4
    resources:
        mem="8G",
        rmem="6G",
    log:
        "logs/samtools/{sample}/{unit}/{reference}/{prefix}-demultimap.log"
    wrapper:
        "0.80.2/bio/samtools/view"


rule samtools_deduplicate:
    input:
        "star/{sample}/{unit}/{reference}/{prefix}.sortedByCoord.out.bam",
    output:
        "star/{sample}/{unit}/{reference}/{prefix}Deduplicated.sortedByCoord.out.bam",
    params:
        extra="-b -F 1024"
    threads: 4 
    resources:
        mem="8G",
        rmem="6G",
    log:
        "logs/samtools/{sample}/{unit}/{reference}/{prefix}-deduplicate.log"
    wrapper:
        "0.80.2/bio/samtools/view"

rule featurecounts:
    input:
        bam="star/{sample}/{unit}/{reference}/{prefix}.sortedByCoord.out.bam",
        bai="star/{sample}/{unit}/{reference}/{prefix}.sortedByCoord.out.bam.bai",
        saf=lambda w: "resources/annotations/{{reference}}_{{lineage}}.{type}.{{valid}}_{{tag}}.{{feature}}.bed.saf".format(
            type="custom" if (w.feature in features["feature_name"].tolist()) else "gtf",
        ), 
    output:
        tab = "star/{sample}/{unit}/{reference}/{prefix}.{lineage}_{valid}.{tag}.{feature}Reads.featurecounts.tab",
    threads: 6
    resources:
        mem=lambda wildcards, input: (str((input.size//3000000000)+4) + "G"),
        rmem=lambda wildcards, input: (str((input.size//6000000000)+4) + "G"),
    log:
        "logs/featurecounts/{sample}/{unit}/{reference}/{prefix}_{lineage}_{valid}.{tag}_{feature}Reads.log"
    conda:
        "../envs/subread.yaml",
    params:
        strand=get_sample_strandedness,
        paired=lambda wildcards:("" if not is_paired_end(wildcards.experiment,wildcards.sample) else "-p"),
        overlap="-O" if config["counting"]["count_every_overlap"] else "",
    shell:
        """
        featureCounts -s {params.strand} {params.paired} -M {params.overlap} -T {threads} -F SAF --verbose -a {input.saf} -o {output.tab} {input.bam}
        """

rule count_matrix:
    input:
        featurecounts=lambda wildcards: expand(
            "star/{sample.sample_name}/{sample.unit_name}/{{reference}}/{{prefix}}.featurecounts.tab",
            sample= get_experiment_samples(wildcards),
        ),
    output:
        counts="featurecounts/{reference}/{experiment}/{prefix}.counts.tsv",
        length="featurecounts/{reference}/{experiment}/{prefix}.lengths.tsv",
    log:
        "logs/{experiment}/{reference}/featurecounts/{prefix}_count_matrix.log",
    params:
        names=lambda wildcards: "\t".join(map(str,samples.loc[wildcards.experiment].sample_name.tolist())),
    resources:
        mem="8G",
        rmem="6G",
    shell:
        """
        paste {input.featurecounts} |
        awk -F'\\t' -v OFS='\\t' '
        FNR==1 {{
          print "gene", "{params.names}"
        }}
        FNR>=3 {{ 
          printf "%s\\t",$1 ; 
          for (i=7;i<=NF;i+=7) {{
            if ((NF-i)>6) {{
              printf "%s\\t", $i
            }}
            else {{
              printf "%s\\n", $i
            }}
          }}
        }}' - > {output.counts} &&
        paste {input.featurecounts} |
        awk -F'\\t' -v OFS='\\t' '
        FNR==1 {{
          print "gene", "Length"
        }}
        FNR>=3 {{ 
          print $1, $6 
        }}' - > {output.length}
        """
