rule feature_splice_sites:
    input:
        gene_tab = "resources/annotations/{reference}/genome.gtf.{tag}_gene_info.tab",
        bed= lambda wildcards: "resources/annotations/{reference}/{lineage}.{type}.{valid}_{tag}.{feature}.bed",
        sj="resources/annotations/{reference}/{lineage}.{source}.splice_junctions.bed",
    output:
        bed="resources/annotations/{reference}/{lineage}.{type}.{valid}_{tag}.{feature}.{source}_splice_sites.bed",
        saf="resources/annotations/{reference}/{lineage}.{type}.{valid}_{tag}.{feature}.{source}_splice_sites.saf",
    params:
        overhang = config["splicing"]["splice_site_overhang"],
    threads: 1
    resources:
        mem="8G",
        rmem="6G",
    log:
        "logs/bedtools/{reference}.{lineage}_{valid}.{type}.{tag}.{feature}.{source}_splice_sites.log"
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
             a = $2 ; b = $3 ; s = $6 ; $4 = $17 ;
             $2 = a-{params.overhang};
             $3 = a+{params.overhang};
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
        splice_bam="star/{sample}/{unit}/{reference}/Spliced{prefix}.sortedByCoord.out.bam",
        splice_bai="star/{sample}/{unit}/{reference}/Spliced{prefix}.sortedByCoord.out.bam.bai",
        unsplice_bam="star/{sample}/{unit}/{reference}/Unspliced{prefix}.sortedByCoord.out.bam",
        unsplice_bai="star/{sample}/{unit}/{reference}/Unspliced{prefix}.sortedByCoord.out.bam.bai",
        saf="resources/annotations/{reference}/{lineage}.{type}.{valid}_{tag}.{feature}.{source}_splice_sites.saf",
    output:
        splice="featurecounts/{sample}/{unit}/{reference}/Spliced{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{source}.SpliceSite.featurecounts.tab",
        unsplice="featurecounts/{sample}/{unit}/{reference}/Unspliced{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{source}.SpliceSite.featurecounts.tab",
    threads: 6
    resources:
        mem=lambda wildcards, input: (str((input.size//4000000000)+4) + "G"),
        rmem=lambda wildcards, input: (str((input.size//8000000000)+4) + "G"),
    log:
        "logs/featurecounts/{sample}/{unit}/{reference}/Spliced{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{source}.SpliceSite.featurecounts.log"
    conda:
        "../envs/subread.yaml",
    params:
        unspl_overlap= config["splicing"]["splice_site_overhang"] + 1 ,
        strand=get_sample_strandedness,
        paired=lambda wildcards:("-p" if is_paired_end(wildcards.sample) else ""),
    shell:
        """
        featureCounts -s {params.strand} {params.paired} --minOverlap 1 -M -O -T {threads} -F SAF --verbose -a {input.saf} -o {output.splice} {input.splice_bam} &&
        featureCounts -s {params.strand} {params.paired} --minOverlap {params.unspl_overlap} -M -O -T {threads} -F SAF --verbose -a {input.saf} -o {output.unsplice} {input.unsplice_bam}
        """

rule dexseq_splice_ratio:
    input:
        spliced= lambda w: "featurecounts/{norm_group}/{{reference}}/Spliced{{prefix}}.{{lineage}}_{{valid}}.{{type}}.{{tag}}.{{feature}}.star.SpliceSite.counts.tsv".format(
            norm_group=experiments.loc[w.experiment,"group_name"],
        ),
        unspliced= lambda w: "featurecounts/{norm_group}/{{reference}}/Unspliced{{prefix}}.{{lineage}}_{{valid}}.{{type}}.{{tag}}.{{feature}}.star.SpliceSite.counts.tsv".format(
            norm_group=experiments.loc[w.experiment,"group_name"],
        ),
        genetab=lambda w: "resources/annotations/{source}/genome.gtf.{{tag}}_gene_info.tab".format(
            source= str( get_sample_source(w.experiment) ),
        ),
        nuc=lambda w: "resources/annotations/{source}/{{lineage}}.{{type}}.{{valid}}_{{tag}}.{{feature}}.bed.nuc.tab".format(
            source=  str( get_sample_source(w.experiment) ),
        ),
        bed=lambda w: "resources/annotations/{source}/{{lineage}}.{{type}}.{{valid}}_{{tag}}.{{feature}}.bed".format(
            source=  str( get_sample_source(w.experiment) ),
        ),
        express="differential/{experiment}/{reference}/differential_expression/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.lfc.tab",
    output:
        lfc="differential/{experiment}/{reference}/differential_splicing_ratio/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.lfc.tab",
        levels="differential/{experiment}/{reference}/differential_splicing_ratio/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.levels.tab",
        counts="differential/{experiment}/{reference}/differential_splicing_ratio/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.counts.tab",
    params:
        control=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["control"],
        treat=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["treatment"],
        paired=lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["pairRep"]),
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/dexseq.yaml"
    log:
        "logs/dexseq/{experiment}/{reference}/differential_splicing_ratio/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/R/dexseq_splice.R"

