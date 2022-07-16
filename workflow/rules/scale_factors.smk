rule base_coverage_scale_factors:
    input:
        lambda wildcards: expand(
            "{{experiment}}/bedgraph/{sample.sample_name}_{sample.unit_name}_{{prefix}}.BaseCoverage.txt",sample=samples.loc[wildcards.experiment].itertuples(),
        ),
    output:
        summary = "{experiment}/deseq2/{prefix}.BaseCoverage.summary.tsv",
        internal = "{experiment}/deseq2/{prefix}.BaseCoverage.internal_scale_factors.tsv",
        spikein = "{experiment}/deseq2/{prefix}.BaseCoverage.spikein_scale_factors.tsv",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/{experiment}/bg2bw/{prefix}_base_coverage.log",
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
            "{{experiment}}/star/{sample.sample_name}-{sample.unit_name}/{{prefix}}.total_read_count.tab",sample=samples.loc[wildcards.experiment].itertuples(),
        ),
    output:
        summary = "{experiment}/deseq2/{prefix}.TotalAlignedReadsCount.summary.tsv",
        internal = "{experiment}/deseq2/{prefix}.TotalAlignedReadsCount.internal_scale_factors.tsv",
        spikein = "{experiment}/deseq2/{prefix}.TotalAlignedReadsCount.spikein_scale_factors.tsv",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/{experiment}/bg2bw/{prefix}_readcount_scale.log",
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
        counts="{experiment}/feature_counts/{prefix}.annotated_basic.{feature}.counts.tsv",
        bed=lambda w: "resources/annotations/{source}.{type}.annotated_basic.{{feature}}.bed".format(
            source = ( str( get_source(w) ) + "genome" ),
            type = "custom" if (w.feature in features["feature_name"].tolist()) else "gtf",
        ),
    output:
        "{experiment}/deseq2/{prefix}.{feature}ReadCount.{spikein}_scale_factors.tsv",
    params:
        paired=lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["paired_analysis"]),
        sample_table="config/samples.tsv",
    resources:
        mem="8G",
        rmem="6G",
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/{experiment}/deseq2/{prefix}_{spikein}_{feature}Reads_scale.log",
    script:
        "../scripts/R/deseq2_feature_scale.R"

