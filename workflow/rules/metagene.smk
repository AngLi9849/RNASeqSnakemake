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
        bed= "resources/annotations/{reference}/{lineage}.plot-{md5}.{valid}_{tag}.{feature}.{sense}_{part}.bed",
        bigwig = "raw_bw/{sample}/{unit}/{reference}/{prefix}.{strand}.raw.bigwig",
    output:
        matrix="matrices/{sample}/{unit}/{reference}/{prefix}.{strand}/{lineage}_{valid}.plot-{md5}.{tag}.{feature}.{sense}_{part}.{bin}bins.matrix.gz",
    log:
        "matrices/{sample}/{unit}/{reference}/{prefix}.{strand}/{lineage}_{valid}.plot-{md5}.{tag}.{feature}.{sense}_{part}.{bin}bins.matrix.log",
    params:
        strand = lambda wildcards:  "+" if (wildcards.strand == "fwd") else "-" if (wildcards.strand == "rev") else "+-",
        temp = "resources/annotations/{reference}/{lineage}.plot-{md5}.{valid}_{tag}.{feature}.{sense}_{part}.{strand}.bed",
        temp_gz = "matrices/{sample}/{unit}/{reference}/{prefix}.{strand}/{lineage}_{valid}.plot-{md5}.{tag}.{feature}.{sense}_{part}.{bin}bins.matrix.temp.gz",
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
        computeMatrix scale-regions -S {input.bigwig} -R {params.temp} -p {threads} --metagene --binSize 1 --averageTypeBins mean --regionBodyLength {wildcards.bin} --sortRegions descend --sortUsing region_length -o {params.temp_gz} &&
        zcat {params.temp_gz} |
        tail -n +2 |
        cut -f4,7- |
        sort -k1,1 |
        gzip - > {output.matrix} &&
        rm {params.temp_gz} && 
        rm {params.temp}
        """


rule feature_signal2background:
    input:
        rpk = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.rpk.bed",
        sense = "resources/annotations/{reference}/{lineage}.{type}.{valid}_{tag}.{feature}.sense.bed",
        genetab = "resources/annotations/{reference}/genome.gtf.{tag}_gene_info.tab",
        background = lambda wildcards: expand("featurecounts/{{norm_group}}/{{reference}}/{{prefix}}.{{lineage}}_{feat.valid}.{feat.type}.{{tag}}.{feat.feature_name}.rpk.bed",
            feat=features.loc[str(features.loc[wildcards.feature,"backgrd"]).split(",")].itertuples()
        ) ,
    output:
        tab = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.plot-{md5}.sig2bg.tab",
    params:
        compat_bt=lambda wildcards: features.loc[wildcards.feature,"comp_bt"],
        bef= lambda wildcards: features.loc[wildcards.feature,"plotbef"],
        aft= lambda wildcards: features.loc[wildcards.feature,"plotaft"],
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
          }}' - {input.rpk} > {params.range} &&

        cat {input.background} |        

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            biotype[$1]=$3 ; 
          }}
          FNR < NR && ($7 >0 ) {{
            if (biotype[$4] != "") {{
              print $0, biotype[$4] ; 
            }} else {{
              print $0, "NA"
            }}
          }}
        ' {input.genetab} - |
 
        bedtools intersect -a {params.range} -b - -s -wao > {params.temp} &&

        awk -F'\\t' -v OFS='\\t' -v OFMT='%f' '
          BEGIN {{
            print "featureID", "sig2bg", "bg2noi"
          }}
          {{
            if ( ($7!=0) && ($21 != 0) && ($9 != $19) && (("{params.compat_bt}" != "nan" && "{params.compat_bt}" !~ $20) || "{params.compat_bt}" == "nan" ) ) {{
              sig_bg=$7/$17 ;
              sb[$8]=((sb[$8]-0)>=sig_bg || sb[$8]=="Inf")?sig_bg:sb[$8] ;
              bs[$8]=($7-0>0)? ( ( (bs[$8]-0) >= ($17/$7) )?(bs[$8]-0):$17/$7 ) :( (bs[$8]-0>0)? bs[$8] : "Inf" ) ;
            }} else {{
              sb[$8]=((sb[$8]-0)>0)?sb[$8]:( ($7==0)?0:"Inf" ) ; 
            }}
          }}
          END {{
            for ( i in sb ) {{
              print i, sb[i], (bs[i]=="Inf")?"Inf":bs[i]-0
            }}
          }}' {params.temp} > {output.tab} 
        """ 


rule sort_raw_matrices:         
    input:
        bef = lambda wildcards : "workflow/null.txt" if int(features.loc[wildcards.feature,"plotbef"])==0 else expand(
          "matrices/{{sample}}/{{unit}}/{{reference}}/{{prefix}}.{strand}/{{lineage}}_{{valid}}.plot-{{md5}}.{{tag}}.{{feature}}.{{sense}}_plotbef.{bin}bins.matrix.gz",
          strand = ["fwd","rev"] if wildcards.strand=="stranded" else "unstranded",
          bin= features.loc[wildcards.feature,"plotbef_bin"],
        ), 
        body = lambda wildcards : expand(
          "matrices/{{sample}}/{{unit}}/{{reference}}/{{prefix}}.{strand}/{{lineage}}_{{valid}}.plot-{{md5}}.{{tag}}.{{feature}}.{{sense}}_main.{bin}bins.matrix.gz",
          strand = ["fwd","rev"] if wildcards.strand=="stranded" else "unstranded",
          bin = features.loc[wildcards.feature,"bin_n"],
        ),
        aft = lambda wildcards : "workflow/null.txt" if int(features.loc[wildcards.feature,"plotaft"])==0 else expand(
          "matrices/{{sample}}/{{unit}}/{{reference}}/{{prefix}}.{strand}/{{lineage}}_{{valid}}.plot-{{md5}}.{{tag}}.{{feature}}.{{sense}}_plotaft.{bin}bins.matrix.gz",
          strand = ["fwd","rev"] if wildcards.strand=="stranded" else "unstranded",
          bin= features.loc[wildcards.feature,"plotaft_bin"],
        ),
    output:
        matrix = "matrices/{sample}/{unit}/{reference}/{prefix}.{strand}/{lineage}_{valid}.plot-{md5}.{tag}.{feature}.{sense}.matrix.gz",
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    resources:
        mem="6G",
        rmem="4G",
    shell:
        """
        zcat {input.bef} {input.body} {input.aft} | gzip - > {output.matrix}
        """  
