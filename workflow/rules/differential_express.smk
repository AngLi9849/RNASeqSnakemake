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
        size_table="deseq2/{experiment}/{reference}/All{prefix}.{normaliser}.{spikein}_scale_factors.tsv",
        length = "featurecounts/{experiment}/{reference}/{splice}{prefix}.{lineage}_{valid}.{tag}.{feature}.lengths.tsv",
        counts="featurecounts/{experiment}/{reference}/{splice}{prefix}.{lineage}_{valid}.{tag}.{feature}Reads.counts.tsv",
        genetab=lambda w: "resources/annotations/{source}_genome.gtf.{{tag}}_gene_info.tab".format(
            source= str( get_sample_source(w) ),
        ),
        bed=lambda w: "resources/annotations/{source}_{{lineage}}.{type}.{{valid}}_{{tag}}.{{feature}}.bed".format(
            source=  str( get_sample_source(w) ),
            type="custom" if (w.feature in features["feature_name"].tolist()) else "gtf",
        ),
        nuc=lambda w: "resources/annotations/{source}_{{lineage}}.{type}.{{valid}}_{{tag}}.{{feature}}.bed.nuc.tab".format(
            source=  str( get_sample_source(w) ),
            type="custom" if (w.feature in features["feature_name"].tolist()) else "gtf",
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
        normcounts="results/{experiment}/{reference}/differential_expression/{spikein}_{normaliser}_normalised.{splice}{prefix}.{lineage}_{valid}.{tag}.{feature}.tsv",
        docx="results/{experiment}/{reference}/differential_expression/{spikein}_{normaliser}_normalised.{splice}_{prefix}.{lineage}_{valid}.{tag}.{feature}.docx",
        rpkm="results/{experiment}/{reference}/differential_expression/{spikein}_{normaliser}_normalised.{splice}_{prefix}.{lineage}_{valid}.{tag}.{feature}.rpkm.tsv",
    params:
        sample_table="config/samples.tsv",
        control=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["control_condition"],
        paired=lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["pairRep"]),
        descript= lambda wildcards: feature_descript(wildcards),
        dir="results/{experiment}/{reference}/differential_expression/{spikein}_{normaliser}_normalised.{splice}_{prefix}.{lineage}_{valid}.{tag}.{feature}",        
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/differential.yaml"
    log:
        "logs/{experiment}/{reference}/deseq2/{splice}ed{prefix}.{lineage}_{valid}.{tag}_{spikein}_{feature}_{normaliser}.diffexp.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/R/deseq2_express.R"
