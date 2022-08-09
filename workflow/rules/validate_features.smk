rule star_detected_splice_junctions:
    input:
        sj=lambda wildcards: expand(
                "star/{sample.sample_name}/{sample.unit_name}/{sample.reference}/SJ.out.tab",
                sample=lineage.loc[wildcards.species].loc[wildcards.lineage].itertuples(),
        ),
    output:
        sj="resources/annotations/{species}.{lineage}.star.splice_junctions.bed",
    threads: 1
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/bedtools.yaml",
    log:
        "logs/awk/{species}_{lineage}_splice_junctions.log",
    shell:
        """
        cat {input.sj} |

        awk -F'\\t' -v OFS='\\t' '
          BEGIN {{
            split("nan|GT:AG|CT:AC|GC:AG|CT:GC|AT:AC|GT:AT",m,"|") ;
          }}
          FNR==NR && $6>=1 {{
            s=$4 ;
            $2=$2-1 ; 
            $4=($4==1)?"+":(($4==2)?"-":"") ;
            name=$1":"$2"-"$3":"$4 ;
            motif=m[($5+1)] ;
            print $1, $2, $3, name, motif, $4, $7, $8, $9
          }}' - |

        sort -k4,4 |

        awk -F'\\t' -v OFS='\\t' '
          name != $4 {{
            save=$0 ; $0=load ; $7=n ; $8=m ; $9=o ;
            if (FNR>1) {{
              print
            }} ;
            $0=save ; name=$4 ; n=$7 ; m=$8 ; o=$9 ; load=$0
          }}
          name==$4 {{
            n+=$7 ; m+=$8 ; o=(o>$9)?o:$9 ;
          }}
          END {{
            $0=load ; $7=n ; $8=m ; $9=o ; print
          }}' - > {output.sj}
        """

rule salmon_transcr


rule validate_features:
    input:
        transcripts="resources/annotations/{prefix}.gtf.{tag}_transcripts.bed",
        gene_tab="resources/annotations/{prefix}.gtf.{tag}_gene_info.tab",
        inconfident="resources/annotations/{prefix}.gtf.{tag}_inconfident.bed",
        confident="resources/annotations/{prefix}.gtf.{tag}_confident.bed",
        sj="resources/annotations/{lineage}.star.splice_junctions.bed",
    output:
        valid_features="resources/annotations/{prefix}_{lineage}.gtf.{tag}_validated.bed"
    params:
        min_overhang=config["feature_validation"]["introns"]["minimum_overhang"],
        min_int_uniq=config["feature_validation"]["introns"]["minimum_unique_splice_reads"],
        min_splice = config["feature_validation"]["introns"]["minimum_multimap_splice_reads"],
        min_ret_cov=config["features"]["minimum_retained_intron_coverage"],
        intron_min=config["features"]["minimum_intron_length"],
        feature_fwd="workflow/scripts/awk/feature_index_fwd.awk",
        feature_rev="workflow/scripts/awk/feature_index_rev.awk",
    threads: 1
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/bedtools.yaml",
    log:
        "logs/awk/{prefix}_{lineage}/validate_{tag}_features.log",
    shell:
        """
# Define genebody range using confident (MANE entries or correct biotype and top tsl transcripts) transcript ranges
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $10=={{ 
            
        cat {input.sj} |

        """ 

