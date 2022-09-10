rule base_coverage_scale_factors:
    input:
        lambda wildcards: expand(
            "bedgraph/{sample.sample_name}/{sample.unit_name}/{{reference}}/{{prefix}}.BaseCoverage.txt",sample=get_norm_group_samples(wildcards.norm_group),
        ),
    output:
        summary = "deseq2/{norm_group}/{reference}/{prefix}.BaseCoverage.summary.tsv",
        internal = "deseq2/{norm_group}/{reference}/{prefix}.BaseCoverage.internal_scale_factors.tsv",
        spikein = "deseq2/{norm_group}/{reference}/{prefix}.BaseCoverage.spikein_scale_factors.tsv",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/scale_factors/{norm_group}/{reference}/{prefix}_base_coverage.log",
    shell:
        """
        cat {input} > {output.summary} &&
        awk -F'\\t' -v OFS='\\t' '
          BEGIN {{
            internal_sum=0 ; spikein_sum=0 ; n=0 ; s="" ; sample=""
          }}
          FNR == NR {{
            internal_sum += $2 ; spikein_sum += $3 ;
            if (s != 1) {{
               n += 1 ; s=$1
            }}
          }}
          FNR<NR {{
            if (sample != $1) {{
              if (FNR==1) {{
                internal_mean=(internal_sum/n) ;
                spikein_mean=(spikein_sum/n) ;
                print "sample_name", "counts", "size_factor" >> "{output.internal}" ;
                print "sample_name", "counts", "size_factor" >> "{output.spikein}" ;
              }}
              else {{
                internal_size = (internal_mean==0)?0:(internal/internal_mean) ;
                spikein_size = (spikein_mean==0)?0:(spikein/spikein_mean) ;
                print sample, internal, internal_size >> "{output.internal}" ;
                print sample, spikein, spikein_size >> "{output.spikein}" ;
              }} ;
              sample = $1 ; internal = $2 ; spikein = $3 ;
            }}
            else {{
              internal += $2 ; spikein += $3
            }}
          }}
          END {{
            internal_size = (internal_mean==0)?0:(internal/internal_mean) ;
            spikein_size = (spikein_mean==0)?0:(spikein/spikein_mean) ;
            print sample, internal, internal_size >> "{output.internal}" ;
            print sample, spikein, spikein_size >> "{output.spikein}" ;
          }}' {output.summary} {output.summary}
        """  

rule total_read_count_size_factors:
    input:
        reads = lambda wildcards: expand(
            "star/{sample.sample_name}/{sample.unit_name}/{{reference}}/{{prefix}}.total_read_count.tab",
            sample=get_norm_group_samples(wildcards.norm_group).itertuples(),
        ),
    output:
        summary = "deseq2/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.gtf.TotalReadCount.{pair}.summary.tsv",
        internal = "deseq2/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.gtf.TotalReadCount.internal_{pair}.scale_factors.tsv",
        spikein = "deseq2/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.gtf.TotalReadCount.spikein_{pair}.scale_factors.tsv",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/scale_factors/{norm_group}/{reference}/{prefix}_{lineage}_{valid}_{pair}_totalreadcount_scale.log",
    shell:
        """
        cat {input} > {output.summary} &&
        awk -F'\\t' -v OFS='\\t' '
          BEGIN {{
            internal_sum=0 ; spikein_sum=0 ; n=0 ; s="" ; sample=""
          }}
          FNR == NR {{
            internal_sum += $2 ; spikein_sum += $3 ;
            if (s != 1) {{
               n += 1 ; s=$1
            }}
          }}
          FNR<NR {{ 
            if (sample != $1) {{
              if (FNR==1) {{
                internal_mean=(internal_sum/n) ; 
                spikein_mean=(spikein_sum/n) ; 
                print "sample_name", "counts", "size_factor" >> "{output.internal}" ;
                print "sample_name", "counts", "size_factor" >> "{output.spikein}" ;
              }}
              else {{
                internal_size = (internal_mean==0)?0:(internal/internal_mean) ;
                spikein_size = (spikein_mean==0)?0:(spikein/spikein_mean) ;
                print sample, internal, internal_size >> "{output.internal}" ;
                print sample, spikein, spikein_size >> "{output.spikein}" ;
              }} ;
              sample = $1 ; internal = $2 ; spikein = $3 ;
            }}
            else {{
              internal += $2 ; spikein += $3
            }}
          }}
          END {{
            internal_size = (internal_mean==0)?0:(internal/internal_mean) ;
            spikein_size = (spikein_mean==0)?0:(spikein/spikein_mean) ;
            print sample, internal, internal_size >> "{output.internal}" ;
            print sample, spikein, spikein_size >> "{output.spikein}" ;
          }}' {output.summary} {output.summary} 
        """

rule feature_count_scale_factors:
    input:
        counts="featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.basic.{feature}Reads.counts.tsv",
        bed="resources/annotations/{reference}/genome.{type}.annotated_basic.{feature}.bed"
    output:
        paired = "deseq2/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{feature}ReadCount.{spikein}_paired.scale_factors.tsv",
        unpaired = "deseq2/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{feature}ReadCount.{spikein}_unpaired.scale_factors.tsv",
    params:
        sample_table="config/samples.tsv",
    resources:
        mem="8G",
        rmem="6G",
    wildcard_constraints:
       feature=r"((?!Total).)*"
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{norm_group}/{reference}/{prefix}_{spikein}.{lineage}_{valid}.{type}.{feature}Reads_scale.log",
    script:
        "../scripts/R/deseq2_feature_scale.R"

