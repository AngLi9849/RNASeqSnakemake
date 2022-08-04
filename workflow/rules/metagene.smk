rule feature_metagene_annotations:
    input:
        bed=lambda w: "{{prefix}}.custom-{id}.{{type}}.{{feature}}.{{sense}}.bed".format(
            id = features.loc[w.feature,"prefix_md5"],
        ),
    output:
        main="{prefix}.plot-{md5}.{type}.{feature}.{sense}_main.bed",
        before="{prefix}.plot-{md5}.{type}.{feature}.{sense}_plotbef.bed",
        after="{prefix}.plot-{md5}.{type}.{feature}.{sense}_plotaft.bed",
    threads: 1 
    params:
        before=lambda w: features.loc[w.feature,"plotbef"],
        after=lambda w: features.loc[w.feature,"plotaft"],
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
              }} else {{
                $2 = (five[$8] >= aft)? ( five[$8] - aft ) : 0   ;
                $3 = five[$8] ;
                print $1, $2, $3, $8, l, $6 >> "{output.after}" ;
                $2 = three[$8] ; 
                $3 = three[$8] + bef
                print $1, $2, $3, $8, l, $6 >> "{output.before}" ;
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

rule filter_expressed_feature:
    input:
        counts = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}Reads.counts.tsv",
        matrix = lambda wildcards : expand(
          "matrices/{{sample}}/{unit.unit}/{{reference}}/{{prefix}}.{strand}/{{lineage}}_{{valid}}.plot-{{md5}}.{{tag}}.{{feature}}.{{sense}}_{{part}}.{{bin}}bins.matrix.gz",
          unit = samples.loc[wildcards.sample_name].itertuples(),
          strand = ["fwd","rev"] if wildcards.strand=="stranded" else "unstranded"
        ),
    output:
        matrix = "matrices/{sample}/{unit}/{reference}/{prefix}.{strand}/{lineage}_{valid}.plot-{md5}.{tag}.{feature}.{sense}_{part}.{bin}bins.min{min}reads.matrix.gz",
    log:
        "logs/matrices/{sample}/{unit}/{reference}/{prefix}.{strand}/{lineage}_{valid}.plot-{md5}.{tag}.{feature}.{sense}_{part}.{bin}bins.min{min}reads.matrix.gz",
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    resources:
        mem="6G",
        rmem="4G",
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '
          {{
