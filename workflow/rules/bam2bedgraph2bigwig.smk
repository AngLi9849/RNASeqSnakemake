rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/samtools_index/{prefix}.log"
    threads: 2 
    resources:
        mem="6G",
        rmem="4G",
    wrapper:
        "0.80.2/bio/samtools/index"


rule stranded_bam:
    input:
        "{prefix}.out.bam",
    output:
        fwd="{prefix}.fwd.bam",
        rev="{prefix}.rev.bam",
    resources:
        mem="6G",
        rmem="4G",
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
        bam="{experiment}/star/{sample}-{unit}/{prefix}.sortedByCoord.out.bam",
        bai="{experiment}/star/{sample}-{unit}/{prefix}.sortedByCoord.out.bam.bai", 
    output:
        "{experiment}/bedgraph/{sample}_{unit}_{prefix}.unstranded.bedgraph",
        "{experiment}/bedgraph/{sample}_{unit}_{prefix}.BaseCoverage.txt",
    params:
        bin_size=config["bigwig_bin_size"]
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/{experiment}/deeptools/{sample}_{unit}_{prefix}_unstranded_bamcoverage.log",
    conda:
        "../envs/bedtools.yaml",
    threads: 2
    shell:
        """
        bedtools genomecov -ibam {input[0]} -bga -split > {output[0]} &&

        awk -v OFS='\\t' -F'\\t' '
          BEGIN {{ i = 0 ; s = 0 }}
          $1 !~ "spikein_" {{i += $4 }}
          $1 ~  "spikein_" {{s += $4 }}
          END {{
            print "{wildcards.sample}", i, s >> "{output[1]}"
          }}' {output[0]}
        """
        
rule stranded_genomecov:
    input:
        fwdbam="{experiment}/star/{sample}-{unit}/{prefix}.sortedByCoord.fwd.bam",
        fwdbai="{experiment}/star/{sample}-{unit}/{prefix}.sortedByCoord.fwd.bam.bai",
        revbam="{experiment}/star/{sample}-{unit}/{prefix}.sortedByCoord.rev.bam",
        revbai="{experiment}/star/{sample}-{unit}/{prefix}.sortedByCoord.rev.bam.bai",
    output:
        bg_fwd="{experiment}/bedgraph/{sample}_{unit}_{prefix}.fwd.bedgraph",
        bg_rev="{experiment}/bedgraph/{sample}_{unit}_{prefix}.rev.bedgraph",
    params:
        bin_size=config["bigwig_bin_size"]
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/{experiment}/deeptools/{sample}_{unit}_{prefix}_stranded_bamcoverage.log",
    conda:
        "../envs/bedtools.yaml",
    threads: 2
    shell:
        """
        bedtools genomecov -ibam {input[0]} -bga -split > {output[0]} &&
        bedtools genomecov -ibam {input[2]} -bga -split > {output[1]}
        """
 
rule scale_bedgraph2bigwig:
    input:
       "{experiment}/deseq2/All{prefix}_{counts}_{normaliser}_scale_factors.tsv",
       "{experiment}/bedgraph/{sample}_{unit}_{splice}{prefix}.{strand}.bedgraph",
       lambda wildcards: ("resources/genomes/" + str(get_source(wildcards)) + ".chrom.sizes")
    output:
       "{experiment}/bedgraph/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{sample}_{unit}.{strand}_{splice}.norm.bedgraph",
       "results/{experiment}/bigwig/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{sample}_{unit}.{strand}_{splice}.bigwig",
    conda:
       "../envs/bedgraphtobigwig.yaml"
    threads: 1
    resources:
        mem="12G",
        rmem="8G",
    log:
       "logs/{experiment}/bg2bw/{sample}_{unit}_{strand}_by_{normaliser}_{counts}_{prefix}_{splice}_bg2bw.log"
    shell:
       """
       awk -F'\\t' -v OFS='\\t' 'FNR==NR{{scalefactor[$1]=$3; next}} {{print $1,$2,$3,$4*scalefactor["{wildcards.sample}"]}}' {input[0]} {input[1]} |
       LC_COLLATE=C sort -k1,1 -k2,2n - > {output[0]} &&
       bedGraphToBigWig {output[0]} {input[2]} {output[1]}
       """

