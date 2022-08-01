rule non_overlap_feature_annotations:
    input:
        bed=lambda w: "{{prefix}}.custom-{id}.{{type}}.{{feature}}.{{sense}}.bed".format(
            id = features.loc[w.feature,"prefix_md5"],
        ),
    output:
        bed="{prefix}.plot-{md5}.{type}.{feature}.{sense}.non-overlap.{strand}.bed",
    threads: 1
    params:
        strand=lambda w: "-s" if w.strand == "stranded" else "",
        distinct= lambda w: "distinct" if w.strand == "stranded" else "collapse",
        before=lambda w: features.loc[w.feature,"plotbef"],
        after=lambda w: features.loc[w.feature,"plotaft"],
    conda:
        "../envs/bedtools.yaml"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/metagene/{prefix}_plot-{md5}.{type}.{feature}.{sense}.non-overlap.{strand}.log",
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '
          $6=="+" && $2 >= {params.before} {{
            a=$2 ; b=$3 ;
            $2=(a-{params.before}) ;
            $3=(b+{params.after}) ;
            print
          }}
          $6=="-" && $2 >= {params.after} {{
            a=$2 ; b=$3 ;
            $2=(a-{params.after}) ;
            $3=(b+{params.before}) ;
            print
          }}' {input.bed} |
        sort -k1,1 -k2,2n - |
        bedtools merge {params.strand} -i - -c 4,5,6,8 -o count,collapse,{params.distinct},collapse |
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $4==1 {{
            n[$7]=1
          }}
          FNR < NR && n[$8] == 1 {{
            print
          }}' - {input.bed} > {output.bed}
        """


