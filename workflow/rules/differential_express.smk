rule pca:
    input:
        "{experiment}/deseq2/{prefix}_{counts}_{normaliser}_deseq2.rds",
    output:
        report("{experiment}/{prefix}_{counts}_{normaliser}_pca.svg", "../report/{prefix}_{counts}_{normaliser}_pca.rst"),
    params:
        pca_labels=config["pca"]["labels"],
    resources:
        mem="8G",
        rmem="6G",
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/{experiment}/{prefix}_{counts}_{normaliser}_pca.log",
    script:
        "../scripts/R/plot-pca.R"


rule deseq2_expression:
    input:
        size_table=lambda w: "deseq2/{{experiment}}/{{reference}}/All{{prefix}}.{norm_type}.{{normaliser}}ReadCount.{{spikein}}_scale_factors.tsv".format(
            norm_type= ("custom-" + str(features.loc[w.normaliser,"prefix_md5"])) if (w.normaliser in features["feature_name"].tolist()) else "gtf",
        ),
        length = "featurecounts/{experiment}/{reference}/{splice}{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}Reads.lengths.tsv",
        counts="featurecounts/{experiment}/{reference}/{splice}{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}Reads.counts.tsv",
        genetab=lambda w: "resources/annotations/{source}_genome.gtf.{{tag}}_gene_info.tab".format(
            source= str( get_sample_source(w.experiment) ),
        ),
        bed=lambda w: "resources/annotations/{source}_{{lineage}}.{{type}}.{{valid}}_{{tag}}.{{feature}}.bed".format(
            source=  str( get_sample_source(w.experiment) ),
        ),
        nuc=lambda w: "resources/annotations/{source}_{{lineage}}.{{type}}.{{valid}}_{{tag}}.{{feature}}.bed.nuc.tab".format(
            source=  str( get_sample_source(w.experiment) ),
        ),
        pptx="resources/templates/{w}cm_wide.{h}cm_tall.pptx".format(
            w=str(config["differential_plots"]["powerpoint"]["width"]),
            h=str(config["differential_plots"]["powerpoint"]["height"]),
        ),
        docx="resources/templates/{c}_{f}.docx".format(
            f=str(config["differential_plots"]["word_docx"]["font_name"]),
            c=str(config["differential_plots"]["word_docx"]["font_colour"]),
        ),
    output:
        normcounts="results/{experiment}/{reference}/differential_expression/{spikein}_{normaliser}ReadCount_normalised.{splice}{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.tsv",
        docx="results/{experiment}/{reference}/differential_expression/{spikein}_{normaliser}ReadCount_normalised.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.docx",
        rpkm="results/{experiment}/{reference}/differential_expression/{spikein}_{normaliser}ReadCount_normalised.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.rpkm.tsv",
    params:
        sample_table="config/samples.tsv",
        control=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["control"],
        treat=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["treatment"],
        paired=lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["pairRep"]),
        descript= lambda wildcards: feature_descript(wildcards),
        dir="results/{experiment}/{reference}/differential_expression/{spikein}_{normaliser}ReadCount_normalised.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}",        
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/differential.yaml"
    log:
        "logs/{experiment}/{reference}/deseq2/{splice}ed{prefix}.{lineage}_{valid}.{type}.{tag}_{spikein}_{feature}_{normaliser}ReadCount.diffexp.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/R/deseq2_express.R"
