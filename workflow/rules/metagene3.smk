rule expressed_non_overlapping_feature:
    input:
        bed = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.rpkm.bed",
        range = resources/annotations/{reference}_{lineage}.plot-{md5}.{valid}_{tag}.{feature}.{sense}_range.bed",
        genetab = "resources/annotations/{reference}_genome.gtf.{tag}_gene_info.tab",
    output:
        biotype_bed = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.min{min}reads.{sense}.biotype_non_overlap.bed",
        all_bed = "featurecounts/{norm_group}/{reference}/{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.min{min}reads.{sense}.all_non_overlap.bed",
    threads: 1
    conda:
        "../envs/bedtools.yaml"
    resources:
        mem="6G",
        rmem="4G",
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '
          FNR < NR {{
            print $1, $3, mean[$1]
          }}' {input.counts} {input.genetab} |
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            biotype[$1] = $2 ;
            mean[$1] = $3 ;
          }}
          FNR < NR && mean[$4] >= {wildcards.min} {{
            print $1, $2, $3, $8, mean[$4], $6, biotype[$4]
          }}' - {input.bed} |
        bedtools merge -s -i - -c 5,6,7 -o collapse,collapse,distinct 
        """  
