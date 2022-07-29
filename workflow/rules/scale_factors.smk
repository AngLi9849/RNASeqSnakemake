rule base_coverage_scale_factors:
    input:
        lambda wildcards: expand(
            "bedgraph/{sample.sample_name}/{sample.unit_name}/{{reference}}/{{prefix}}.BaseCoverage.txt",sample=get_experiment_samples(wildcards),
        ),
    output:
        summary = "deseq2/{experiment}/{reference}/{prefix}.BaseCoverage.summary.tsv",
        internal = "deseq2/{experiment}/{reference}/{prefix}.BaseCoverage.internal_scale_factors.tsv",
        spikein = "deseq2/{experiment}/{reference}/{prefix}.BaseCoverage.spikein_scale_factors.tsv",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/scale_factors/{experiment}/{reference}/{prefix}_base_coverage.log",
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
                print "sample_name", "counts", "scale_factor" >> "{output.internal}" ;
                print "sample_name", "counts", "scale_factor" >> "{output.spikein}" ;
              }}
              else {{
                internal_scale = (internal_mean==0)?0:(internal/internal_mean) ;
                spikein_scale = (spikein_mean==0)?0:(spikein/spikein_mean) ;
                print sample, internal, internal_scale >> "{output.internal}" ;
                print sample, spikein, spikein_scale >> "{output.spikein}" ;
              }} ;
              sample = $1 ; internal = $2 ; spikein = $3 ;
            }}
            else {{
              internal += $2 ; spikein += $3
            }}
          }}
          END {{
            internal_scale = (internal_mean==0)?0:(internal/internal_mean) ;
            spikein_scale = (spikein_mean==0)?0:(spikein/spikein_mean) ;
            print sample, internal, internal_scale >> "{output.internal}" ;
            print sample, spikein, spikein_scale >> "{output.spikein}" ;
          }}' {output.summary} {output.summary}
        """  

rule read_count_scale_factors:
    input:
        reads = lambda wildcards: expand(
            "star/{sample.sample_name}/{sample.unit_name}/{{reference}}/{{prefix}}.total_read_count.tab",sample=get_experiment_samples(wildcards),
        ),
    output:
        summary = "deseq2/{experiment}/{reference}/{prefix}.TotalAlignedReadsCount.summary.tsv",
        internal = "deseq2/{experiment}/{reference}/{prefix}.TotalAlignedReadsCount.internal_scale_factors.tsv",
        spikein = "deseq2/{experiment}/{reference}/{prefix}.TotalAlignedReadsCount.spikein_scale_factors.tsv",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/scale_factors/{experiment}/{reference}/{prefix}_readcount_scale.log",
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
                print "sample_name", "counts", "scale_factor" >> "{output.internal}" ;
                print "sample_name", "counts", "scale_factor" >> "{output.spikein}" ;
              }}
              else {{
                internal_scale = (internal_mean==0)?0:(internal/internal_mean) ;
                spikein_scale = (spikein_mean==0)?0:(spikein/spikein_mean) ;
                print sample, internal, internal_scale >> "{output.internal}" ;
                print sample, spikein, spikein_scale >> "{output.spikein}" ;
              }} ;
              sample = $1 ; internal = $2 ; spikein = $3 ;
            }}
            else {{
              internal += $2 ; spikein += $3
            }}
          }}
          END {{
            internal_scale = (internal_mean==0)?0:(internal/internal_mean) ;
            spikein_scale = (spikein_mean==0)?0:(spikein/spikein_mean) ;
            print sample, internal, internal_scale >> "{output.internal}" ;
            print sample, spikein, spikein_scale >> "{output.spikein}" ;
          }}' {output.summary} {output.summary} 
        """

rule feature_count_scale_factors:
    input:
        counts="featurecounts/{experiment}/{reference}/{prefix}.genome_annotated.{type}.basic.{feature}Reads.counts.tsv",
        bed="resources/annotations/{reference}_genome.{type}.annotated_basic.{feature}.bed"
    output:
        "deseq2/{experiment}/{reference}/{prefix}.{type}.{feature}ReadCount.{spikein}_scale_factors.tsv",
    params:
        paired=lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["pairRep"]),
        sample_table="config/samples.tsv",
    resources:
        mem="8G",
        rmem="6G",
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{experiment}/{reference}/{prefix}_{spikein}.{type}.{feature}Reads_scale.log",
    script:
        "../scripts/R/deseq2_feature_scale.R"

