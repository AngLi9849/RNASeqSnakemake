rule extract_annotated_feature:
    input:
        "resources/annotations/{prefix}/{lineage}.gtf.{tag}_{valid}.bed",
    output:
        "resources/annotations/{prefix}/{lineage}.gtf.{valid}_{tag}.{feature}.bed"
    log:
        "logs/features/extract_{prefix}/{lineage}_{valid}_{tag}_{feature}.log"
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
        bed="resources/annotations/{source}/{lineage}.{type}.{valid}_{tag}.{feature}.bed",
    output:
        info="resources/annotations/{source}/{lineage}.{type}.{valid}_{tag}.{feature}.bed.nuc.tab",
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
          a[$4] += $17 ; 
          c[$4] += $18 ; 
          g[$4] += $19 ;
          t[$4] += $20 ; 
          n[$4] += $21 ;
          other[$4] += $22 ;
          len[$4] += $23 ; 
        }}
        END {{
          for (i in len) {{
            print i, (a[i]+t[i])/len[i], (g[i]+c[i])/len[i], a[i], c[i], g[i], t[i], n[i], other[i]
          }}
        }}' - |
        sort -k1,1 - |
        sed  '1 i\\featureID\\tAT\\tGC\\tA\\tC\\tG\\tT\\tN\\tOther' > {output.info}
        """

rule custom_feature:
    input:
        region=lambda wildcards: ("{{prefix}}.{s}.{{type}}.{f}.bed".format(
            s=features.loc[features.loc[wildcards.feature,"feature"],"type"] if (features.loc[wildcards.feature,"region"] in features["feature_name"].tolist()) else "gtf",
            f=features.loc[wildcards.feature,"region"]
        )),
        exclude=lambda wildcards: ("{{prefix}}.{s}.{{type}}.{f}.bed".format(
            s=features.loc[features.loc[wildcards.feature,"feature"],"type"] if (features.loc[wildcards.feature,"exclude"] in features["feature_name"].tolist()) else "gtf",
            f=features.loc[wildcards.feature,"exclude"]
        )),
        group=lambda wildcards: ("{{prefix}}.{s}.{{type}}.{f}.bed".format(
            s=features.loc[features.loc[wildcards.feature,"feature"],"type"] if (features.loc[wildcards.feature,"group"] in features["feature_name"].tolist()) else "gtf",
            f=features.loc[wildcards.feature,"group"]
        )),
        feature=lambda wildcards: ("{{prefix}}.{s}.{{type}}.{f}.bed".format(
            s=features.loc[features.loc[wildcards.feature,"feature"],"type"] if (features.loc[wildcards.feature,"feature"] in features["feature_name"].tolist()) else "gtf",
            f=features.loc[wildcards.feature,"feature"]
        )),
    threads: 1
    resources:
        mem="4G",
        rmem="6G",
    output:
        sense = "{prefix}.custom-{md5}.{type}.{feature}.sense.bed",
        antisense = "{prefix}.custom-{md5}.{type}.{feature}.antisense.bed",
        bed="{prefix}.custom-{md5}.{type}.{feature}.bed"
    params:
        feat=lambda wildcards: features.loc[wildcards.feature,"feature"],
        group=lambda wildcards: features.loc[wildcards.feature,"group"],
        sect=lambda wildcards: features.loc[wildcards.feature,"section"],
        sense=lambda wildcards: str(features.loc[wildcards.feature,"sense"]),
        no_frst=lambda wildcards: features.loc[wildcards.feature,"no_frst"],
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
        bedtools intersect -s -a {input.group} -b - |
        bedtools intersect -wb -s -a {input.feature} -b - |
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
              print ; 
              print ("{params.sense}" == "-+")?$0:"" >> "{output.antisense}" ;
              print (("{params.sense}" ~ /^+/) || ("{params.sense}" == "nan") )?$0:"" >> "{output.sense}" ; 
            }} ;            
            if ( ("{params.sense}" ~ "-") || ("{params.sense}" == "nan") ) {{
              $6 = "-" ; print ;
              print ("{params.sense}" ~ /^-/)?$0:"" >> "{output.sense}" ;
              print ("{params.sense}" == "+-")?$0:"" >> "{output.antisense}" ;
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
              print ;
              print ("{params.sense}" == "-+")?$0:"" >> "{output.antisense}" ;
              print (("{params.sense}" ~ /^+/) || ("{params.sense}" == "nan") )?$0:"" >> "{output.sense}" ;
            }} ;
            if ( ("{params.sense}" ~ "-") || ("{params.sense}" == "nan") ) {{
              $6 = "+" ; print ;
              print ("{params.sense}" ~ /^-/)?$0:"" >> "{output.sense}" ;
              print ("{params.sense}" == "+-")?$0:"" >> "{output.antisense}" ;
            }} ;
          }}
        }}' - |
        sort -k7,7 -k4,4 -k2,2n -k3,3n |

        awk -F'\\t' -v OFS='\\t' '
          $2>=0 {{ print }} 
        ' - |
        uniq > {output.bed}
        """ 


rule feature_rpk:
    input:
        counts = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}Reads.counts.tsv",
        bed = "resources/annotations/{reference}/{lineage}.{type}.{valid}_{tag}.{feature}.bed",
    output:
        bed = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.rpk.bed",
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    resources:
        mem="6G",
        rmem="4G",
    shell:
        """
        sort -k6,6 -k8,8 -k2,2n {input.bed} |
        awk -F'\\t' -v OFS='\\t' -v OFMT='%f' '
          BEGIN {{
            id=""
          }}
          FNR == NR {{
            len[$8,$6] += $5 ;
            five[$8,$6]=(five[$8,$6]<=1)?$2:((five[$8,$6] <= $2)?five[$8,$6]:$2) ;
            three[$8,$6]=(three[$8,$6]>=$3)?three[$8,$6]:$3 ;
          }}
          FNR < NR {{
            if ( id != $8":"$6 ) {{
              $2=five[$8,$6] ;
              $3=three[$8,$6] ;
              $5=len[$8,$6] ;
              print
            }} ;
            id=$8":"$6
          }}' {input.bed} - |         
        awk -F'\\t' -v OFS='\\t' -v OFMT='%f' '
          BEGIN {{
            total=0
          }}
          FNR==NR && FNR > 1 {{
            sum[$1]=0 ;
            for (i=2 ; i <=NF ; i ++ ) {{
              sum[$1] += $i ;
              total += $i ;
            }} 
          }}
          FNR < NR {{
# chr start end root length strand rpkm feature_id parent
            print $1, $2, $3, $4, $5, $6,sum[$4]*1000000000/$5, $8, $11
          }}' {input.counts} - > {output.bed}
        """

