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
        size_table=lambda w: "deseq2/{norm_group}/{{reference}}/All{{prefix}}.{{lineage}}_{{valid}}.{norm_type}.{{normaliser}}.{{norm_read}}.Count.{{spikein}}_{{pair}}.scale_factors.tsv".format(
            norm_type= ("custom-" + str(features.loc[w.normaliser,"prefix_md5"])) if (w.normaliser in features["feature_name"].tolist()) else "gtf",
            norm_group=experiments.loc[w.experiment,"group_name"],
        ),
        length = lambda w: "featurecounts/{norm_group}/{{reference}}/{{splice}}{{prefix}}.{{lineage}}_{{valid}}.{{type}}.{{tag}}.{{feature}}.{{read}}.lengths.tsv".format(
            norm_group=experiments.loc[w.experiment,"group_name"],
        ),
        counts= lambda w: "featurecounts/{norm_group}/{{reference}}/{{splice}}{{prefix}}.{{lineage}}_{{valid}}.{{type}}.{{tag}}.{{feature}}.{{read}}.counts.tsv".format(
            norm_group=experiments.loc[w.experiment,"group_name"],
        ),
        total_sum = lambda w: "deseq2/{norm_group}/{{reference}}/{{splice}}{{prefix}}.{{lineage}}_{{valid}}.gtf.TotalReadCount.{{pair}}.summary.tsv".format(
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
        base_bed=lambda wildcards: ("resources/annotations/{source}/{{lineage}}.{s}.{{valid}}_{{tag}}.{f}.bed".format(
            source=  str( get_sample_source(wildcards.experiment) ),
            s=features.loc[features.loc[wildcards.feature,"group"],"type"] if (features.loc[wildcards.feature,"group"] in features["feature_name"].tolist()) else "gtf",
            f=features.loc[wildcards.feature,"group"]
        ) ),
    output:
        counts="differential/{experiment}/{reference}/differential_expression/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.counts.tab",
        lfc="differential/{experiment}/{reference}/differential_expression/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.lfc.tab",
        levels="differential/{experiment}/{reference}/differential_expression/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.levels.tab",
        rpkm="differential/{experiment}/{reference}/differential_expression/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.rpkm.tsv",
        toptable="differential/{experiment}/{reference}/differential_expression/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.toptable.tsv",
    params:
        control=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["control"],
        treat=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["treatment"],
        paired=lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["pairRep"]),
        section=lambda wildcards : features.loc[wildcards.feature,"section"],
        main_int = lambda wildcards : str(features.loc[wildcards.feature,"is_main_int"]),
        dir="results/{experiment}/{reference}/differential_expression/{spikein}_{normaliser}.{norm_read}.Count_normalised.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}",
    resources:
        mem="24G",
        rmem="16G",
    conda:
        "../envs/differential.yaml"
    log:
        "logs/{experiment}/{reference}/deseq2/{splice}ed{prefix}.{lineage}_{valid}.{type}.{tag}_{spikein}_{feature}_{read}.{pair}_{normaliser}.{norm_read}_normalised.Count.diffexp.log",
    threads: 1
    script:
        "../scripts/R/deseq2_express.R"

