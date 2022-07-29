rule extract_annotated_feature:
    input:
        "resources/annotations/{prefix}_{lineage}.gtf.{tag}_{valid}.bed",
    output:
        "resources/annotations/{prefix}_{lineage}.gtf.{valid}_{tag}.{feature}.bed"
    log:
        "logs/features/extract_{prefix}_{lineage}_{valid}_{tag}_{feature}.log"
    threads: 1
    resources:
        mem="4G",
        rmem="6G",
    shell:
        """
        awk -F'\\t' '$7=="{wildcards.feature}"  {{print}}' {input} > {output}
        """

rule provided_feature:
    input:
        bed = lambda wildcards: expand(
           features.loc[wildcards.feature].squeeze(axis=0)["annotation_bed"]
        ),
    output:
        bed = "resources/annotations/{prefix}.custom-{md5}.provided_basic.{feature}.bed" 
    log:
        "logs/features/provided-{md5}_{prefix}_{feature}.log"
    threads: 1
    resources:
        mem="4G",
        rmem="6G",
    shell:
        """
        cat {input.bed} |  
        sort -k1,1 -k2,2n -k3,3n |

        awk -F'\\t' -v OFS='\\t' '
          BEGIN {{ 
            length(FNR)=digit ;
            n=1 ;
          }}
          FNR==NR {{ 
            if (NR<4 || $4=="" || $4=="." || $4=="-") {{
              for (d=1 ; d<=digits ; d+=1) index=index"0" ;
              index = index n ;
              index=substr(index, 1 + length(index) - digits) ;
              $4="{wildcards.feature}"index ;
              n+=1 ;
            }} ;
            $5=$3-$2 ;          
            if (NR<6 || $6=="") {{
              print $1, $2, $3, $4, $5, "+", $4, $4, $4, 0, $4, 1, 1, 1 ;
              print $1, $2, $3, $4, $5, "-", $4, $4, $4, 0, $4, 1, 1, 1 ;
            }} 
            else {{
              print $0, $4, $4, $4, 0, $4, 1, 1, 1 ;
            }}
          }}' -  > {output.bed} 
         """ 
           
rule feature_nuc_info:
    input:
        fasta="resources/genomes/{source}_genome.fasta",
        bed="resources/annotations/{source}_{lineage}.{type}.{valid}_{tag}.{feature}.bed",
    output:
        info="resources/annotations/{source}_{lineage}.{type}.{valid}_{tag}.{feature}.bed.nuc.tab",
    log:
        "logs/features/nuc_info/{source}_{lineage}.{type}.{valid}_{tag}.{feature}.log"
    threads: 1
    resources:
        mem="4G",
        rmem="6G",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools nuc -C -s -fi {input.fasta} -bed {input.bed} |
        awk -F'\\t' -v OFS='\\t' '
        FNR > 1 {{
          print $4, $15, $16, $17, $18, $19, $20, $21, $22
        }}' - |
        sort - |
        uniq - |
        sed  '1 i\\featureID\\tAT\\tGC\\tA\\tC\\tG\\tT\\tN\\tOther' > {output.info}
        """

rule custom_feature:
    input:
        region=lambda wildcards: ("{{prefix}}.{s}.{{type}}.{f}.bed".format(
            s="custom" if (features.loc[wildcards.feature,"region"] in features["feature_name"].tolist()) else "gtf",
            f=features.loc[wildcards.feature,"region"]
        )),
        exclude=lambda wildcards: ("{{prefix}}.{s}.{{type}}.{f}.bed".format(
            s="custom" if (features.loc[wildcards.feature,"exclude"] in features["feature_name"].tolist()) else "gtf",
            f=features.loc[wildcards.feature,"exclude"]
        )),
        group=lambda wildcards: ("{{prefix}}.{s}.{{type}}.{f}.bed".format(
            s="custom" if (features.loc[wildcards.feature,"group"] in features["feature_name"].tolist()) else "gtf",
            f=features.loc[wildcards.feature,"group"]
        )),
        feature=lambda wildcards: ("{{prefix}}.{s}.{{type}}.{f}.bed".format(
            s="custom" if (features.loc[wildcards.feature,"feature"] in features["feature_name"].tolist()) else "gtf",
            f=features.loc[wildcards.feature,"feature"]
        )),
    threads: 1
    resources:
        mem="4G",
        rmem="6G",
    output:
        bed="{prefix}.custom-{md5}.{type}.{feature}.bed"
    params:
        feat=lambda wildcards: features.loc[wildcards.feature,"feature"],
        group=lambda wildcards: features.loc[wildcards.feature,"group"],
        sect=lambda wildcards: features.loc[wildcards.feature,"section"],
        sense=lambda wildcards: features.loc[wildcards.feature,"sense"],
        no_first=lambda wildcards: features.loc[wildcards.feature,"no_frst"],
        no_last=lambda wildcards: features.loc[wildcards.feature,"no_last"],
        min_len=lambda wildcards: features.loc[wildcards.feature,"min_len"],
        tsl=lambda wildcards: features.loc[wildcards.feature,"tsl"],
        before=lambda wildcards: features.loc[wildcards.feature,"len_bef"],
        after=lambda wildcards: features.loc[wildcards.feature,"len_aft"],
        annotation=lambda wildcards: features.loc[wildcards.feature,"annotation_bed"],
    log:
        "logs/features/{prefix}/custom-{md5}_{type}_{feature}.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools subtract -s -a {input.region} -b {input.exclude} |
        sort -k1,1 -k2,2n |
        bedtools merge -s -i - -c 4,5,6 -o collapse,collapse,distinct |
        bedtools intersect -s -a {input.group} -b stdin |
        bedtools intersect -wa -wb -s -a {input.feature} -b stdin |
        awk -F'\\t' -v OFS='\\t' '
          $4==$18 && $10 <= {params.tsl} && $5 >= {params.min_len} {{
            $11=$22; print 
          }}' - |
        cut -f1-14 |
        sort -k7,7 -k4,4 -k2,2n -k3,3n |
        uniq |
        awk -F'\\t' -v OFS='\\t' '{{
          if ($6=="+") {{ 
            $7="{wildcards.feature}" ; a=$2 ; b=$3 ;
            $2 = ( \
              ("{params.sect}"=="body") ? (a - {params.before}) : ( \
                ("{params.sect}"=="start") ? (a - {params.before}) : ( \
                  ("{params.sect}"=="end") ? (b - {params.before}) : ( \
                    ("{params.sect}"=="centre") ? ( int( (b+a)/2 + 0.5 ) - {params.before} ) : (a - {params.before}) \
                  ) \
                ) \
              ) \
            ) ;
            $3 = ( \
              ("{params.sect}"=="body") ? (b + {params.after}) : ( \
                ("{params.sect}"=="start") ? (a + {params.after}) : ( \
                  ("{params.sect}"=="end") ? (b + {params.after}) : ( \
                    ("{params.sect}"=="centre") ? ( int( (b+a)/2 + 0.5 ) + {params.after} ) : (b + {params.after})  \
                  ) \
                ) \
              ) \
            ) ;
            $5 = $3 - $2 ;
            if ( ("{params.sense}" ~ "+") || ("{params.sense}" == "nan") ) {{
              print
            }} ;
            if ( ("{params.sense}" ~ "-") || ("{params.sense}" == "nan") ) {{
              $6 = "-" ; print
            }} ;
          }} 
          else if ($6=="-") {{
            $7="{wildcards.feature}" ; a=$2 ; b=$3 ;
            $2 = ( \
              ("{params.sect}"=="body") ? (a - {params.after}) : ( \
                ("{params.sect}"=="start") ? (a - {params.after}) : ( \
                  ("{params.sect}"=="end") ? (b - {params.after}) : ( \
                    ("{params.sect}"=="centre") ? ( int( (b+a)/2 - 0.5 ) - {params.after} ) : (a - {params.after}) \
                  ) \
                ) \
              ) \
            ) ;
            $3 = ( \
              ("{params.sect}"=="body") ? (b + {params.before}) : ( \
                ("{params.sect}"=="start") ? (a + {params.before}) : ( \
                  ("{params.sect}"=="end") ? (b + {params.before}) : ( \
                    ("{params.sect}"=="centre") ? ( int( (b+a)/2 - 0.5 ) + {params.before} ) : (b + {params.before}) \
                  ) \
                ) \
              ) \
            ) ;
            $5 = $3 - $2 ;
            if ( ("{params.sense}" ~ "+") || ("{params.sense}" == "nan") ) {{
              print
            }} ;
            if ( ("{params.sense}" ~ "-") || ("{params.sense}" == "nan") ) {{
              $6 = "+" ; print
            }} ;
          }}
        }}' - |
        sort -k7,7 -k4,4 -k2,2n -k3,3n |

        awk -F'\\t' -v OFS='\\t' '
          $2>=0 {{ print }} 
        ' - |
        uniq > {output.bed}
        """ 

rule star_detected_splice_junctions:
    input:
        sj=lambda wildcards: expand(
                "star/{sample.sample_name}/{sample.unit_name}/{reference}/SJ.out.tab",
                sample=get_lineage_sj_samples(wildcards).itertuples(),
        ),
    output:
        sj="resources/annotations/{species}.{lineage}.star.splice_junctions.bed",
    threads: 1
    resources:
        mem="16G",
        rmem="12G",
    params:
        min_overhang = config["feature_validation"]["introns"]["minimum_overhang"],
        min_unique = config["feature_validation"]["introns"]["minimum_unique_splice_reads"],
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
            $4=($4=1)?"+":(($4=2)?"-":"") ;
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
          }}' - | 

        awk -F'\\t' -v OFS='\\t' '
          $9 >= {params.min_overhang} && $7 >= {params.min_unique} {{ 
            print $1, $2, $3, $4, $5, $6
          }}
        ' - > {output.sj}
        """

rule validate_features:
    input:
        bed="resources/annotations/{prefix}_genome.gtf.{tag}_annotated.bed",
        sj="resources/annotations/{lineage}.star.splice_junctions.bed",
    output:
        valid_features="resources/annotations/{prefix}_{lineage}.gtf.{tag}_validated.bed"
    params:
        min_overhang=config["feature_validation"]["introns"]["minimum_overhang"],
        min_int_uniq=config["feature_validation"]["introns"]["minimum_unique_splice_reads"],
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
        cat {input.sj} |

        """ 

