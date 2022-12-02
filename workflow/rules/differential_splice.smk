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
        ram=int(0.75 * config["max_ram_gb"] * 1000000),
        unspl_overlap= config["splicing"]["splice_site_overhang"] + 1 ,
        strand=get_sample_strandedness,
        paired=lambda wildcards:("-p" if is_paired_end(wildcards.sample) else ""),
    shell:
        """
        featureCounts -s {params.strand} {params.paired} --minOverlap 1 -M -O -T {threads} -F SAF --verbose -a {input.saf} -o {output.splice} {input.splice_bam} &&
        if [[ $(du {input.splice_bam} | cut -f1) -gt {params.ram} ]] ;
        then
          for i in $(cut -f2 {input.saf} | sort | uniq) ; do
            samtools view -bh -@ 5 {input.unsplice_bam} $i > {input.unsplice_bam}"$i".bam &&
            samtools index -b -@ 5 {input.unsplice_bam}"$i".bam {input.unsplice_bam}"$i".bam.bai &&
            featureCounts -s {params.strand} {params.paired} --minOverlap {params.unspl_overlap} -M -O -T {threads} -F SAF --verbose -a {input.saf} -o {output.unsplice}"$i".chr.tab {input.unsplice_bam}"$i".bam ;
          done &&

          for i in $(cut -f2 {input.saf} | sort | uniq) ; do
            cat {output.unsplice}"$i".chr.tab ;
          done |

          awk -F'\\t' -v OFS='\\t' '
            FNR<3 {{
              print >> "{output.unsplice}"
            }}
            FNR>2 && $1!="Geneid" && NF>5 {{
              entry[$1]=$0 ;
              count[$1]+=$7 ;
            }}
            END {{
              for (i in entry) {{
                $0=entry[i] ;
                $7=count[i] ;
                print ;
              }}
            }}
           ' - |
           sort -k1,1 >> {output.unsplice}

        else
        featureCounts -s {params.strand} {params.paired} --minOverlap {params.unspl_overlap} -M -O -T {threads} -F SAF --verbose -a {input.saf} -o {output.unsplice} {input.unsplice_bam}
        fi
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
        express="differential/{experiment}/{reference}/differential_expression/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.lfc.tab",
        base_bed=lambda wildcards: ("resources/annotations/{source}/{{lineage}}.{s}.{{valid}}_{{tag}}.{f}.bed".format(
            source=  str( get_sample_source(wildcards.experiment) ),
            s=features.loc[features.loc[wildcards.feature,"group"],"type"] if (features.loc[wildcards.feature,"group"] in features["feature_name"].tolist()) else "gtf",
            f=features.loc[wildcards.feature,"group"]
        ) ),
    output:
        lfc="differential/{experiment}/{reference}/differential_splicing_ratio/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.lfc.tab",
        levels="differential/{experiment}/{reference}/differential_splicing_ratio/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.levels.tab",
        counts="differential/{experiment}/{reference}/differential_splicing_ratio/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.counts.tab",
    params:
        control=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["control"],
        treat=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["treatment"],
        paired=lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["pairRep"]),
        section=lambda wildcards : features.loc[wildcards.feature,"section"],
        main_int = lambda wildcards : str(features.loc[wildcards.feature,"is_main_int"]),
    resources:
        mem="24G",
        rmem="16G",
    conda:
        "../envs/dexseq.yaml"
    log:
        "logs/dexseq/{experiment}/{reference}/differential_splicing_ratio/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/R/dexseq_splice.R"

