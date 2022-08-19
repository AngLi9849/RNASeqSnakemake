rule base_coverage_scale_factors:
    input:
        lambda wildcards: expand(
            "{{experiment}}/bedgraph/{sample.sample_name}_{sample.unit_name}_{{prefix}}.unstranded.{{normaliser}}_sum.txt",sample=samples.loc[wildcards.experiment].itertuples(),
        ),
    output:
        "{experiment}/deseq2/{prefix}_BaseCoverage_{normaliser}_scale_factors.tsv",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/{experiment}/bg2bw/{prefix}_scale_{normaliser}_coverage.log",
    shell:
        """
        cat {input} |
        awk -F'\\t' -v OFS='\\t' '
          BEGIN {{
            sum=0
          }}
          FNR == NR {{
            sum += $2 ;
            if ($1 != s) {{
               n += 1 ; s=$1
            }}
          }}
          FNR<NR && ($1 != sample) {{
            if FNR==1{{
              mean=(sum/n) ;
              print "sample_name", "counts", "scale_factor" >> "{output}" ;
            }}
            else {{
              print sample, count, (count/mean) >> "{output}" ;
            }} ;
            sample = $1 ; count = $2 ;
          }}
          FNR < NR && ($1 == sample) {{
            count += $2
          }}
          END {{
            print sample, count, (count/mean) >> "{output}" ;
          }}
        """  

rule read_count_scale_factors:
    input:
        reads = lambda wildcards: expand(
            "{{experiment}}/star/{sample.sample_name}-{sample.unit_name}/total_read_count.tab",sample=samples.loc[wildcards.experiment].itertuples(),
        ),
    output:
        internal = "{experiment}/counts/{prefix}_TotalReadCount_internal_scale_factors.tsv"
        spikein = "{experiment}/deseq2/{prefix}_TotalReadCount_spikein_scale_factors.tsv",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/{experiment}/bg2bw/{prefix}_scale_{normaliser}_CPM.log",
    shell:
        """
        cat {input} |
        awk -F'\\t' -v OFS='\\t' '
          BEGIN {{
            internal_sum=0 ; spikein_sum=0
          }}
          FNR == NR {{
            internal_sum += $2 ; spikein_sum += $3 ;
            if ($1 != s) {{
               n += 1 ; s=$1
            }}
          }}
          FNR<NR && ($1 != sample) {{
            if FNR==1{{
              internal_mean=(internal_sum/n) ; 
              spikein_mean=(spikein_sum/n) ; 
              print "sample_name", "counts", "scale_factor" >> "{output.internal}" ;
              print "sample_name", "counts", "scale_factor" >> "{output.spikein}" ;
            }}
            else {{
              print sample, internal, (internal/internal_mean) >> "{output.internal}" ;
              print sample, spikein, (spikein/spikein_mean) >> "{output.spikein}" ;
            }} ;
            sample = $1 ; internal = $2 ; spikein = $3 ;
          }}
          FNR < NR && ($1 == sample) {{
            internal += $2 ; spikein += $3
          }}
          END {{
            print sample, internal, (internal/internal_mean) >> "{output.internal}" ;
            print sample, spikein, (spikein/spikein_mean) >> "{output.spikein}" ;
          }}
        """


rule scale_bedgraph2bigwig:
    input:
       "{experiment}/deseq2/All{prefix}_{counts}_{normaliser}_scale_factors.tsv",
       "{experiment}/bedgraph/{sample}_{unit}_{splice}{prefix}.{strand}.bedgraph",
       lambda wildcards: (str(select_genome(wildcards)) + ".chrom.sizes")
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
    

