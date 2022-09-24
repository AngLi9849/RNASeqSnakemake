rule feature_metagene_annotations:
    input:
        bed=lambda w: "{{prefix}}.custom-{id}.{{type}}.{{feature}}.{{sense}}.bed".format(
            id = features.loc[w.feature,"prefix_md5"],
        ),
    output:
        main="{prefix}.plot-{md5}.{type}.{feature}.{sense}_main.bed",
        before="{prefix}.plot-{md5}.{type}.{feature}.{sense}_plotbef.bed",
        after="{prefix}.plot-{md5}.{type}.{feature}.{sense}_plotaft.bed",
        range = "{prefix}.plot-{md5}.{type}.{feature}.{sense}_range.bed",
    params:
        before=lambda w: features.loc[w.feature,"plotbef"],
        after=lambda w: features.loc[w.feature,"plotaft"],
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/metagene/{prefix}_plot-{md5}.{type}.{feature}.{sense}.plot.log",
    shell:
        """
        sort -k6,6 -k1,1 -k8,8 -k2,2n -k3,3n {input.bed} | 
        awk -F'\\t' -v OFS='\\t' -v id='' '
          FNR==NR {{
            len[$8] += $5 ; 
            five[$8]=(five[$8]<=1)?$2:((five[$8] <= $2)?five[$8]:$2) ;
            three[$8]=(three[$8]>=$3)?three[$8]:$3 ;
            print $1, $2, $3, $8, $5, $6 >> "{output.main}" ;
          }}
          FNR < NR {{
            if ( id != $8 ) {{
              l=len[$8] ;
              bef=(match("{params.before}","([^x]*)x",b))? int(l*b) : {params.before} ;
              aft=(match("{params.after}","([^x]*)x",a))? int(l*a) : {params.after} ;
              if (($6=="+" && "{wildcards.sense}" == "sense") || ($6=="-" && "{wildcards.sense}" == "antisense")) {{ 
                $2 = ( five[$8] >= bef ) ? ( five[$8] - bef ) : 0 ;
                $3 = five[$8] ;
                print $1, $2, $3, $8, l, $6 >> "{output.before}" ;
                $2 = three[$8] ;
                $3 = three[$8] + aft ;
                print $1, $2, $3, $8, l, $6 >> "{output.after}" ;
                $2 = ( five[$8] >= bef ) ? ( five[$8] - bef ) : 0 ;
                $3 = three[$8] + aft ;
                print $1, $2, $3, $8, l, $6 >> "{output.range}" ;
              }} else {{
                $2 = (five[$8] >= aft)? ( five[$8] - aft ) : 0   ;
                $3 = five[$8] ;
                print $1, $2, $3, $8, l, $6 >> "{output.after}" ;
                $2 = three[$8] ; 
                $3 = three[$8] + bef
                print $1, $2, $3, $8, l, $6 >> "{output.before}" ;
                $2 = (five[$8] >= aft)? ( five[$8] - aft ) : 0   ;
                $3 = three[$8] + bef
                print $1, $2, $3, $8, l, $6 >> "{output.range}" ;
              }} ;
              id = $8 ;
            }}
          }}' {input.bed} -
        """


rule compute_raw_matrix:
    input:        
        bed= "resources/annotations/{reference}_{lineage}.plot-{md5}.{valid}_{tag}.{feature}.{sense}_{part}.bed",
        bigwig = "bigwigs/{sample}/{unit}/{reference}/{prefix}.{strand}.raw.bigwig",
    output:
        matrix="matrices/{sample}/{unit}/{reference}/{prefix}.{strand}/{lineage}_{valid}.plot-{md5}.{tag}.{feature}.{sense}_{part}.{bin}bins.matrix.gz",
    log:
        "matrices/{sample}/{unit}/{reference}/{prefix}.{strand}/{lineage}_{valid}.plot-{md5}.{tag}.{feature}.{sense}_{part}.{bin}bins.matrix.log",
    params:
        strand = lambda wildcards:  "+" if (wildcards.strand == "fwd") else "-" if (wildcards.strand == "rev") else "+-",
        temp = "temp/{sample}/{unit}/{reference}/{prefix}.{strand}/{lineage}_{valid}.plot-{md5}.{tag}.{feature}.{sense}_{part}.{strand}.bed",
    threads: 4
    resources:
        mem="10G",
        rmem="6G",
    conda:
        "../envs/deeptools.yaml",
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '
          "{params.strand}" ~ $6 {{
            print
          }}
        ' {input.bed} > {params.temp} &&
        computeMatrix scale-regions -S {input.bigwig} -R {params.temp} -p {threads} --skipZeros --metagene --binSize 1 --averageTypeBins mean --regionBodyLength {wildcards.bin} --sortRegions descend --sortUsing region_length -o {output.matrix} &&
        rm {params.temp}
        """


rule expressed_non_overlapping_feature:
    input:
        rpkm = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.rpkm.bed",
        sense = "resources/annotations/{reference}/{lineage}.{type}.{valid}_{tag}.{feature}.sense.bed",
        genetab = "resources/annotations/{reference}/genome.gtf.{tag}_gene_info.tab",
        noise = lambda wildcards: expand("featurecounts/{{norm_group}}/{{reference}}/{{prefix}}.{{lineage}}_{feat.valid}.{feat.type}.{{tag}}.{feat.feature_name}.rpkm.bed",
            feat=features.loc[str(features.loc[wildcards.feature,"noise"]).split(",")].itertuples()
        ) ,
    output:
        bed = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.plot-{md5}.non_overlap.bed",
    params:
        compat_bt=lambda wildcards: features.loc[wildcards.feature,"comp_bt"],
        bef= lambda wildcards: features.loc[wildcards.feature,"plotbef"],
        aft= lambda wildcards: features.loc[wildcards.feature,"plotaft"],
        sig_min= lambda wildcards: features.loc[wildcards.feature,"s2n_min"],
        sig_max= lambda wildcards: features.loc[wildcards.feature,"s2n_max"],        
        range = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.plot-{md5}.range.bed",        
        temp="featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.plot-{md5}.temp.bed"
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    resources:
        mem="6G",
        rmem="4G",
    shell:
        """
        cat {input.sense} |

        awk -F'\\t' -v OFS='\\t' ' 
          FNR==NR {{
            strand[$8]=$6
          }}
          FNR<NR {{
            print $0, "range" ;
            start=$2 ; 
            end=$3 ;
            if ( ( {params.bef}!=0 && strand[$8]=="+" ) || ({params.aft}!=0 && strand[$8]=="-") ) {{
              before=(strand[$8]=="+")?{params.bef}:{params.aft} ; 
              $2=(start>=before)?start-before:0 ; 
              $3=start ;
              print $0, "outer"
            }} ; 
            if ( ( {params.aft}!=0 && strand[$8]=="+" ) || ({params.bef}!=0 && strand[$8]=="-") ) {{
              $2=end ;
              $3=end + ( (strand[$8]=="+")?{params.aft}:{params.bef} ) ;
              print $0, "outer"
            }} ;
          }}' - {input.rpkm} > {params.range} &&

        cat {input.noise} |        

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            biotype[$1]=$3 ; 
          }}
          FNR < NR {{
            if (biotype[$4] != "") {{
              print $0, biotype[$4] ; 
            }} else {{
              print $0, "NA"
            }}
          }}
        ' {input.genetab} - |
 
        bedtools intersect -a {params.range} -b - -s -wao > {params.temp} &&

        awk -F'\\t' -v OFS='\\t' '
          BEGIN {{
            print "__UNLIST__"
          }}
          FNR ==NR && ($10=="range") && ($21 != 0) && ($9 != $19) && (("{params.compat_bt}" != "nan" && "{params.compat_bt}" !~ $20) || "{params.compat_bt}" == "nan" ) {{
              noise=(($17*$15)-($7*$21))/$15 ; 
              print noise ; 
              if ( ({params.sig_max} > {params.sig_min} && noise*{params.sig_max} < $7 ) || noise*{params.sig_min} > $7) {{
                print $8;
                checked[$18]=1 ;
              }}
            }}
          FNR < NR && ($10=="outer") && ($21 != 0) && ($9 != $19) && (("{params.compat_bt}" != "nan" && "{params.compat_bt}" !~ $20) || "{params.compat_bt}" == "nan" ) {{
            if ( (checked[$18] != 1) && (({params.sig_max} > {params.sig_min} && $17*{params.sig_max}<$7) || $17*{params.sig_min} > $7) ) {{
                print $8 ;
            }}
          }}' {params.temp} {params.temp} > "test.list.txt" && 
  
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{ 
            unlist[$1]=1 ;
          }}
          FNR < NR {{
            if (unlist[$8]==1) {{
              next
            }} else {{
              print $0
            }}
          }}' "test.list.txt" {input.rpkm} > {output.bed}
        """ 


rule sort_raw_matrices:         
    input:
        matrix = lambda wildcards : expand(
          "matrices/{{sample}}/{unit.unit}/{{reference}}/{{prefix}}.{strand}/{{lineage}}_{{valid}}.plot-{{md5}}.{{tag}}.{{feature}}.{{sense}}_{{part}}.{{bin}}bins.matrix.gz",
          unit = samples.loc[wildcards.sample_name].itertuples(),
          strand = ["fwd","rev"] if wildcards.strand=="stranded" else "unstranded",
        ),
    output:
        matrix = "matrices/{sample}/{unit}/{reference}/{prefix}.{strand}/{lineage}_{valid}.plot-{md5}.{tag}.{feature}.{sense}_{part}.{bin}bins.matrix.gz",
    log:
        "logs/matrices/{sample}/{unit}/{reference}/{prefix}.{strand}/{lineage}_{valid}.plot-{md5}.{tag}.{feature}.{sense}_{part}.{bin}bins.matrix.gz",
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    resources:
        mem="6G",
        rmem="4G",
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '
          {{print}}
        """  
