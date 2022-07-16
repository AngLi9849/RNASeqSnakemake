rule mv_spliced_bam:
    input:
        "{prefix}Aligned{midfix}counts{suffix}",
    output:
        "{prefix}AllAligned{midfix}counts{suffix}",
    shell:
        "mv {input} {output}"

rule mv_Unspliced_bam:
    input:
        "{prefix}pliced{midfix}counts{suffix}",
    output:
        "{prefix}plicedAligned{midfix}counts{suffix}",
    shell:
        "mv {input} {output}"

rule test:
    input:
        "test.a.bed",
    output:
        "test.tsv",
    params:
        a=expand(
            "{s.sample_name}",s=samples.itertuples()
        ), 
        b=["A","B","C"],
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '{{
          print {params.a} ; 
        }}' {input} > {output}
        """
