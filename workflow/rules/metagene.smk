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
          FNR < NR && $6=="+" {{
            if ( id != $8 ) {{
              l=len[$8] ;
              bef=("{params.before}" ~ /\\..*[1-9]/)? int(l*{params.before}) : {params.before} ;
              aft=("{params.after}" ~ /\\..*[1-9]/)? int(l*{params.after}) : {params.after} ;
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
          }}' - {input.bed} 
        """


rule compute_matrix:
    input:        
        bed=resources/annotations/{reference}_{lineage}.plot-{md5}.{type}.{feature}.{sense}_{part}.bed,
        bigwig = "results/{norm_group}/{reference}/bigwigs/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}{prefix}/{sample}_{unit}.{strand}_{splice}.coverage.bigwig",
    output:
        matrix="compute_matrix/{norm_group}/{reference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}{prefix}/{sample}_{unit}.{strand}_{splice}/{lineage}.{type}.{feature}.plot-{md5}.{sense}_{part}.matrix.gz",
    log:
        "logs/metagene/by_{normaliser}_{counts}_{splice}{prefix}/{biotype}_promptTSS/{sample}_{unit}.matrix.log",
    params:
        bin_num = lambda wildcards: get_part_number(wildcards)
    threads: 4
    resources:
        mem="10G",
        rmem="6G",
    conda:
        "../envs/deeptools.yaml",
    shell:
        """
        computeMatrix scale-regions -S {input.bigwig} -R {input.bed} -p {threads} --metagene --binSize 1 --averageTypeBins mean --regionBodyLength {params.bin_num} --sortRegions descend --sortUsing region_length -o {output.matrix}
        """
    
