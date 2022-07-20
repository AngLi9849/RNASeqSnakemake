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
        bam="star/{prefix}.sortedByCoord.out.bam",
        bai="star/{prefix}.sortedByCoord.out.bam.bai", 
    output:
        "bedgraph/{prefix}.unstranded.bedgraph",
        "bedgraph/{prefix}.BaseCoverage.txt",
    params:
        bin_size=config["bigwig_bin_size"]
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/deeptools/{prefix}_unstranded_bamcoverage.log",
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
        fwdbam="star/{prefix}.sortedByCoord.fwd.bam",
        fwdbai="star/{prefix}.sortedByCoord.fwd.bam.bai",
        revbam="star/{prefix}.sortedByCoord.rev.bam",
        revbai="star/{prefix}.sortedByCoord.rev.bam.bai",
    output:
        bg_fwd="bedgraph/{prefix}.fwd.bedgraph",
        bg_rev="bedgraph/{prefix}.rev.bedgraph",
    params:
        bin_size=config["bigwig_bin_size"]
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/deeptools/{prefix}_stranded_bamcoverage.log",
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
       scale = "deseq2/{experiment}/All{prefix}_{counts}_{normaliser}_scale_factors.tsv",
       bedgraph = "bedgraph/{sample}/{unit}/{reference}/{splice}{prefix}.{strand}.bedgraph",
       chr_size = lambda wildcards: ("resources/genomes/" + str(wildcards.reference) + ".fasta.chrom.sizes")
    output:
       bg = "bedgraph/{experiment}/{reference}/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{sample}_{unit}.{strand}_{splice}.norm.bedgraph",
       bw = "results/{experiment}/{reference}/bigwigs/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{sample}_{unit}.{strand}_{splice}.coverage.bigwig",
    conda:
       "../envs/bedgraphtobigwig.yaml"
    threads: 1
    resources:
        mem="12G",
        rmem="8G",
    log:
       "logs/{experiment}/bg2bw/{sample}_{unit}/{reference}/{strand}_by_{normaliser}_{counts}_{prefix}_{splice}_bg2bw.log"
    shell:
       """
       awk -F'\\t' -v OFS='\\t' '
         FNR==NR{{scalefactor[$1]=$3}} 
         FNR < NR {{
           print $1,$2,$3,$4*scalefactor["{wildcards.sample}"]
         }}' {input.scale} {input.bedgraph} |

       LC_COLLATE=C sort -k1,1 -k2,2n - > {output.bg} &&
       bedGraphToBigWig {output.bw} {input.chr_size} {output.bg}
       """

