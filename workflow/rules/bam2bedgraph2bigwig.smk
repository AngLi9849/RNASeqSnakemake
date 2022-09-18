rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/samtools_index/{prefix}.log"
    threads: 4 
    resources:
        mem="6G",
        rmem="4G",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index -b -@ 3 {input} {output}
        """

rule stranded_bam:
    input:
        "{prefix}.out.bam",
    output:
        fwd="{prefix}.fwd.bam",
        rev="{prefix}.rev.bam",
    resources:
        mem="12G",
        rmem="8G",
    conda:
        "../envs/samtools.yaml",
    shell:
        """
        samtools view -h -b -f 128 -F 16 {input} > {wildcards.prefix}.a.bam && 
        samtools view -h -b -f 80 {input} > {wildcards.prefix}.b.bam && 
        samtools merge -f {output.fwd} {wildcards.prefix}.a.bam {wildcards.prefix}.b.bam &&
        samtools view -h -b -f 144 {input} > {wildcards.prefix}.c.bam &&
        samtools view -h -b -f 64 -F 16 {input} > {wildcards.prefix}.d.bam &&
        samtools merge -f {output.rev} {wildcards.prefix}.c.bam {wildcards.prefix}.d.bam &&
        rm -r {wildcards.prefix}.["a","b","c","d"].bam
        """


rule unstranded_genomecov:
    input:
        bam="star/{sample}/{unit}/{reference}/{prefix}.sortedByCoord.out.bam",
        bai="star/{sample}/{unit}/{reference}/{prefix}.sortedByCoord.out.bam.bai", 
        chr_size = lambda wildcards: ("resources/genomes/" + str(wildcards.reference) + ".fasta.chrom.sizes")
    output:
        bg="bedgraph/{sample}/{unit}/{reference}/{prefix}.unstranded.bedgraph",
        txt="bedgraph/{sample}/{unit}/{reference}/{prefix}.BaseCoverage.txt",
        bw="bigwigs/{sample}/{unit}/{reference}/{prefix}.unstranded.raw.bigwig",
    params:
        bin_size=config["bigwig_bin_size"]
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/deeptools/{sample}/{unit}/{reference}/{prefix}_unstranded_bamcoverage.log",
    conda:
        "../envs/bedtools.yaml",
    threads: 2
    shell:
        """
        bedtools genomecov -ibam {input.bam} -bga -split |

        awk -v OFS='\\t' -F'\\t' '
          BEGIN {{ i = 0 ; s = 0 }}
          $1 !~ "spikein_" {{i += $4 }}
          $1 ~  "spikein_" {{s += $4 }}
          END {{
            print "{wildcards.sample}", i, s >> "{output.txt}"
          }}' |
        LC_COLLATE=C sort -k1,1 -k2,2n - > {output.bg} &&
        bedGraphToBigWig {output.bg} {input.chr_size} {output.bw} ;
        """
        
rule stranded_genomecov:
    input:
        bam="star/{sample}/{unit}/{reference}/{prefix}.sortedByCoord.{strand}.bam",
        bai="star/{sample}/{unit}/{reference}/{prefix}.sortedByCoord.{strand}.bam.bai",
        chr_size = lambda wildcards: ("resources/genomes/" + str(wildcards.reference) + "_genome.fasta.chrom.sizes")
    output:
        bg="bedgraph/{sample}/{unit}/{reference}/{prefix}.{strand}.bedgraph",
        bw="raw_bw/{sample}/{unit}/{reference}/{prefix}.{strand}.raw.bigwig",
    params:
        bin_size=config["bigwig_bin_size"]
    resources:
        mem="6G",
        rmem="4G",
    wildcard_constraints:
        strand=r"((?!unstranded).)*"
    log:
        "logs/deeptools/{sample}/{unit}/{reference}/{prefix}.{strand}_bamcoverage.log",
    conda:
        "../envs/bedtools.yaml",
    threads: 2
    shell:
        """
        bedtools genomecov -ibam {input.bam} -bga -split |
        LC_COLLATE=C sort -k1,1 -k2,2n - > {output.bg} &&
        bedGraphToBigWig {output.bg} {input.chr_size} {output.bw} ;
        """
 
rule scale_bedgraph2bigwig:
    input:
       scale = lambda w: "deseq2/{{norm_group}}/{{reference}}/All{{prefix}}.{lineage}_{valid}.{norm_type}.{{normaliser}}ReadCount.{{spikein}}_{{pair}}.scale_factors.tsv".format(
            lineage=results.loc[w.sample].loc[w.unit,"diff_lineage"][0],
            valid=VALID,
            norm_type= ("custom-" + str(features.loc[w.normaliser,"prefix_md5"])) if (w.normaliser in features["feature_name"].tolist()) else "gtf",
        ),
       bedgraph = "bedgraph/{sample}/{unit}/{reference}/{splice}{prefix}.{strand}.bedgraph",
       chr_size = lambda wildcards: ("resources/genomes/" + str(wildcards.reference) + "_genome.fasta.chrom.sizes")
    output:
       bg = "bedgraph/{norm_group}/{reference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}{prefix}/{sample}_{unit}.{strand}_{splice}.norm.bedgraph",
       bw = "norm_bw/{norm_group}/{reference}/bigwigs/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}{prefix}/{sample}_{unit}.{strand}_{splice}.normalised.bigwig",
    conda:
       "../envs/bedtools.yaml"
    threads: 1
    resources:
        mem="12G",
        rmem="8G",
    log:
       "logs/{norm_group}/bg2bw/{sample}_{unit}/{reference}/{strand}_by_{spikein}_{pair}_{normaliser}_{prefix}_{splice}_{strand}_bg2bw.log"
    shell:
       """
       awk -F'\\t' -v OFS='\\t' '
         FNR==NR{{scalefactor[$1]=$3}} 
         FNR < NR {{
           print $1,$2,$3,$4*scalefactor["{wildcards.sample}"]
         }}' {input.scale} {input.bedgraph} |

       LC_COLLATE=C sort -k1,1 -k2,2n - > {output.bg} &&
       bedGraphToBigWig {output.bg} {input.chr_size} {output.bw}
       """

