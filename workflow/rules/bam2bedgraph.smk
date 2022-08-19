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
        "{experiment}/bedgraph/{sample}_{unit}_{prefix}.unstranded.internal_sum.txt",
        "{experiment}/bedgraph/{sample}_{unit}_{prefix}.unstranded.spikein_sum.txt",
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
            print "{wildcards.sample}", i" >> "{output[1]}"
            print "{wildcards.sample}", s" >> "{output[2]}"
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
 

