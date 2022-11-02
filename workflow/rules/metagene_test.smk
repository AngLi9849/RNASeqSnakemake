rule feature_metagene_annotations:
    input:
        bed=lambda w: "{{prefix}}.custom-{id}.{{type}}.{{feature}}.{{sense}}.main.bed".format(
            id = features.loc[w.feature,"prefix_md5"],
        ),
    output:
        before="{prefix}.plot-{md5}.{type}.{feature}.{sense}.before.bed",
        after="{prefix}.plot-{md5}.{type}.{feature}.{sense}.after.bed",
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
          }}
          FNR < NR && $6=="+" {{
            if ( id != $8 ) {{
              l=len[$8] ;
              bef=("{params.before}" ~ /\\..*[1-9]/)? int(l*{params.before}) : {params.before} ;
              aft=("{params.after}" ~ /\\..*[1-9]/)? int(l*{params.after}) : {params.after} ;
              if (($6=="+" && "{wildcards.sense}" == "sense") || ($6=="-" && "{wildcards.sense}" == "antisense")) {{
                $2 = ( five[$8] >= bef ) ? ( five[$8] - bef ) : 0 ;
                $3 = five[$8] ;
                print >> "{output.before}" ;
                $2 = three[$8] ;
                $3 = three[$8] + aft ;
                print >> "{output.after}" ;
              }} else {{
                $2 = (five[$8] >= aft)? ( five[$8] - aft ) : 0   ;
                $3 = five[$8] ;
                print >> "{output.after}" ;
                $2 = three[$8] ;
                $3 = three[$8] + bef
                print >> "{output.before}" ;
              }} ;
              id = $8 ;
            }}
          }}' - {input.bed}
        """




