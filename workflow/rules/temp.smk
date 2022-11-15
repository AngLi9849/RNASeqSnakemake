rule mv_Unspliced_bam:
    input:
        "{prefix}pliced{midfix}counts{suffix}",
    output:
        "{prefix}plicedAligned{midfix}counts{suffix}",
    shell:
        "mv {input} {output}"

rule test:
    output:
        "test.tsv",
    params:
        a=expand(
            "{s.sample_name}",s=samples.itertuples()
        ), 
        b=["A","B","C"],
    shell:
        """
        awk -F'\\t' -v OFS='\\t' ' BEGIN{{
          print {params.a};
          print "\\"" ; 
        }}' > {output}
        """

rule test2:
    input:
        fq=get_salmon_fq, 
    output:
        b="test.{sample}/{unit}.txt"
    params:
        fqs=lambda wildcards, input: get_salmon_input(wildcards, input)
    shell:
        """
        echo "{params.fqs}" >> {output}
        """  

rule validate_test:
    input:
        bed="resources/annotations/{reference}/genome.gtf.bed",
        transcripts="resources/annotations/{reference}/genome.gtf.{tag}_transcripts.bed",
        trs_idx="resources/annotations/{reference}/genome.gtf.{tag}_transcripts.indexed.bed",
        trs_sj="resources/annotations/{reference}/genome.gtf.{tag}_trs_sj.bed",
        salmon_all = lambda wildcards: expand(
            "salmon/{sample.sample_name}/{sample.unit_name}/{sample.reference}/{{tag}}_transcripts/quant.sf",
            sample=lineage[lineage.trs_val.tolist()].loc[get_reference_species(wildcards.reference)].loc[wildcards.lineage].itertuples(),
        ),
        salmon_confident = lambda wildcards: expand(
            "salmon/{sample.sample_name}/{sample.unit_name}/{sample.reference}/{{tag}}_confident/quant.sf",
            sample=lineage[lineage.trs_val.tolist()].loc[get_reference_species(wildcards.reference)].loc[wildcards.lineage].itertuples(),
        ),
        gene_tab="resources/annotations/{reference}/genome.gtf.{tag}_gene_info.tab",
        sj=lambda wildcards: "resources/annotations/{species}.{{lineage}}.star.splice_junctions.bed".format(
            species=get_reference_species(wildcards.reference),
        ),
        confident="resources/annotations/{reference}/genome.gtf.{tag}_confident.bed",
    threads: 1
    resources:
        mem="16G",
        rmem="12G",
    output:
        form="resources/annotations/{reference}/{lineage}.gtf.{tag}_form.tab",
    shell:
        """
        cat {input.salmon_confident} |

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $1!="Name" {{
            rpk[$1]+=($5/$3) ;
            nreads[$1]+=$5 ;
            rpksum+=($5/$3) ;
          }}
          FNR < NR {{
            if (nreads[$8]>0) {{
              form[$4][$15]= 1 ;
              form_rpk[$4][$15]+=rpk[$8] ;
              form_reads[$4][$15]+=nreads[$8] ;
            }} ;
          }}
          END {{
            for (i in form) {{
              for (j in form[i]) {{
                print i, j, form_rpk[i][j], form_reads[i][j] ;
              }}
            }}
          }}
        ' - {input.transcripts} |

        sort -k1,1 -k4,4nr > {output.form}

        """
