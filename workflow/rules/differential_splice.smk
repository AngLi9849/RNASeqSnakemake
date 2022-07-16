rule feature_splice_sites:
    input:
        gene_tab = "resources/annotations/{species}_genome.gtf.{tag}_gene_info.tab",
        bed= lambda wildcards: "resources/annotations/{{species}}_{g}.{{type}}.{{valid}}_{{tag}}.{{feature}}.bed".format(
            g=wildcards.lineage if wildcards.valid=="validated" else "genome",
        ),
        sj = "resources/annotations/{species}.{lineage}.{source}.splice_junctions.bed"
    output:
        bed="resources/annotations/{species}.{lineage}.{type}.{valid}_{tag}.{feature}.{source}_splice_sites.bed",
        saf="resources/annotations/{species}.{lineage}.{type}.{valid}_{tag}.{feature}.{source}_splice_sites.saf",
    params:
        overhang = config["splicing"]["splice_site_overhang"],
    threads: 1
    resources:
        mem="8G",
        rmem="6G",
    log:
        "logs/bedtools/{species}.{lineage}_{valid}.{type}.{tag}.{feature}.{source}_splice_sites.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            multi[$1]=($4>1)?1:0
          }}
          FNR<NR && multi[$4]==1 {{
            print
          }}
        ' {input.gene_tab} {input.bed} |

        bedtools intersect -f 1 -wa -wb -a {input.sj} -b stdin |

        awk -F'\\t' -v OFS='\\t' '
          {{
             a = $2 ; b = $3 ; s = $6 ; $4 = $14 ;
             $2 = a-{params.overhang};
             $3 = a-{params.overhang};
             $7 = (s==1)? 5 : 3 ;
             print $1, $2, $3, $4, $5, $6, $7;
             $2 = b-{params.overhang} ;
             $3 = b+{params.overhang} ;
             $7 = (s==1)? 3 : 5 ;
             print $1, $2, $3, $4, $5, $6, $7;
          }}
        ' - |

        sort -k1,1 -k2,2n - |
        uniq - > {output.bed} &&

        awk -F'\\t' -v OFS='\\t' '{{print $4,$1,$2+1,$3,$6}}' {output.bed} |
        sort -k2,2 -k3,3n - |
        uniq - |
        sed '1 i\\GeneID\\tChr\\tStart\\tEnd\\tStrand' > {output.saf}
        """

rule splice_site_featurecount:
    input:
        splice_bam="{experiment}/star/{sample}-{unit}/Spliced{prefix}.sortedByCoord.out.bam",
        splice_bai="{experiment}/star/{sample}-{unit}/Spliced{prefix}.sortedByCoord.out.bam.bai",
        unsplice_bam="{experiment}/star/{sample}-{unit}/Unspliced{prefix}.sortedByCoord.out.bam",
        unsplice_bai="{experiment}/star/{sample}-{unit}/Unspliced{prefix}.sortedByCoord.out.bam.bai",
        saf=lambda w: "resources/annotations/{species}.{{lineage}}.{type}.{{valid}}_{{tag}}.{{feature}}.{{source}}_splice_sites.saf".format(
            species= str( get_sample_source(w) ),
            type="custom" if (w.feature in features["feature_name"].tolist()) else "gtf",
        ),
    output:
        splice="{experiment}/star/{sample}-{unit}/Spliced{prefix}.{lineage}_{valid}.{tag}.{feature}.{source}.SpliceSite.featurecounts.tab",
        unsplice="{experiment}/star/{sample}-{unit}/Unspliced{prefix}.{lineage}_{valid}.{tag}.{feature}.{source}.SpliceSite.featurecounts.tab",
    threads: 6
    resources:
        mem=lambda wildcards, input: (str((input.size//4000000000)+4) + "G"),
        rmem=lambda wildcards, input: (str((input.size//8000000000)+4) + "G"),
    log:
        "logs/feature_counts/{experiment}/{sample}-{unit}/Spliced{prefix}.{lineage}_{valid}.{tag}.{feature}.{source}.SpliceSiteReads.log"
    conda:
        "../envs/subread.yaml",
    params:
        unspl_overlap= config["splicing"]["splice_site_overhang"] + 1 ,
        strand=get_sample_strandedness,
        paired=lambda wildcards:("" if not is_paired_end(wildcards.experiment,wildcards.sample) else "-p")
    shell:
        """
        featureCounts -s {params.strand} {params.paired} --minOverlap 1 -M -O -T {threads} -F SAF --verbose -a {input.saf} -o {output.splice} {input.splice_bam} &&
        featureCounts -s {params.strand} {params.paired} --minOverlap {params.unspl_overlap} -M -T {threads} -F SAF --verbose -a {input.saf} -o {output.unsplice} {input.unsplice_bam}
        """

rule dexseq_splice_ratio:
    input:
        spliced="{experiment}/counts/Spliced{prefix}.{lineage}_{valid}.{tag}.{feature}.{source}.SpliceSite.counts.tsv",
        unspliced="{experiment}/counts/Unspliced{prefix}.{lineage}_{valid}.{tag}.{feature}.{source}.SpliceSite.counts.tsv",
        genetab=lambda wildcards: (str(get_annotation(experiments.loc[wildcards.experiment,"sample_source"])) + ".{tag}_gene_info.tab"),
    output:
        dir=directory("results/{experiment}/differential_splicing/All{prefix}")
    params:
        plot_script="workflow/scripts/R/differential_plots.R",
        sample_table="config/samples.tsv",
        biotypes=config["biotypes"],
        goi=config["GOI"],
        control=lambda wildcards: experiments.loc[wildcards.experiment,"control_condition"],
        paired=lambda wildcards: str(experiments.loc[wildcards.experiment,"paired_analysis"]),
        dir="results/{experiment}/differential_splicing/All{prefix}",
        ma_number=config["differential_plots"]["ma_gene_name_numbers"],
        volc_number=config["differential_plots"]["volcano_gene_name_numbers"],
        up_col=config["differential_plots"]["up_colour"],
        down_col=config["differential_plots"]["down_colour"],
        p_threshold=config["differential_plots"]["p_value_threshold"],
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/dexseq.yaml"
    log:
        "logs/dexseq/{experiment}/Aligned{prefix}.diffsplice.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/R/dexseq_splice.R"

