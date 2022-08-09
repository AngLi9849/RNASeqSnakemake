rule differential_plots:
    input:
        lfc="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.lfc.tab",
        levels="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.levels.tab",
        bed=lambda w: "resources/annotations/{source}/{{lineage}}.{{type}}.{{valid}}_{{tag}}.{{feature}}.bed".format(
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
        docx="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.docx",
    params:
        protocol = lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["protocol"],
        control=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["control"],
        treat=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["treatment"],
        paired=lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["pairRep"]),
        descript= lambda wildcards: feature_descript(wildcards),
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/differential.yaml"
    log:
        "logs/differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.plots.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/R/differential_plots.R"
