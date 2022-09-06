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
        range = "resources/annotations/{reference}_{lineage}.plot-{md5}.{valid}_{tag}.{feature}.{sense}_range.bed",
        genetab = "resources/annotations/{reference}_genome.gtf.{tag}_gene_info.tab",
        background = lambda wildcards: "featurecounts/{{norm_group}}/{{reference}}/{{prefix}}.{{lineage}}_{valid}.{type}.{{tag}}.{base}.rpkm.bed".format(
            type = get_feature_type(features.loc[wildcards.feature,"feature"]),
            valid = get_feature_validity(features.loc[wildcards.feature,"feature"]),
            base = features.loc[wildcards.feature,"feature"],
        ),
    output:
        biotype_bed = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{sense}.{sig}sig2noise.biotype_non_overlap.bed",
        all_bed = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{sense}.{sig}sig2noise.all_non_overlap.bed",
    params:
        sig_noi=lambda wildcards: features.loc[wildcards.feature,"sig2noi"]
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    resources:
        mem="6G",
        rmem="4G",
    shell:
        """
        bedtools intersect -a {input.rpkm} -b {input.background} -s -wa -wb |
        awk -F'\\t' -v OFS='\\t' '{{
          if ($5 >= {wildcards.sig}*$13) {{
            print $1, $2, $3, $4, $5, $6, $7, $8 

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            biotype[$1] = $3
          }}
          FNR < NR {{
            $9 = biotype[$4] ;
            print
          }}' {input.genetab} {input.rpkm} |
        sort -k1,1 -k2,2n | 
        bedtools merge -s -i - -c 4, 6,7,8,9 -o count, distinct,collapse,collapse,distinct |
         
        bedtools intersect -a {input.rpkm} -b {input.background} -s - 
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{ 
            
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            biotype[$1] = $2 ;
            mean[$1] = $3 ;
          }}
          FNR < NR && mean[$4] >= 1 {{
            print $1, $2, $3, $8, mean[$4], $6, biotype[$4]
          }}' - {input.bed} |
        bedtools merge -s -i - -c 5,6,7 -o collapse,collapse,distinct |
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
