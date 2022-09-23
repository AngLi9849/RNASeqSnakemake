rule differential_plots:
    input:
        lfc="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.lfc.tab",
        levels="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.levels.tab",
        counts="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.counts.tab",
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
        docx="diff_plots/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.docx",
    params:
        genesets=lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["gene_sets"]).split(","),
        goi= lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["GOI"]).split(","),
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
    threads: 1
    script:
        "../scripts/R/differential_plots.R"

rule feature_reports:
    input:
        docx="resources/templates/{c}_{f}.docx".format(
            f=str(config["differential_plots"]["word_docx"]["font_name"]),
            c=str(config["differential_plots"]["word_docx"]["font_colour"]),
        ),        
        plots=lambda wildcards: expand(
            "diff_plots/{{experiment}}/{{reference}}/differential_{feat.diff}/{{pair}}.{{spikein}}_{{normaliser}}ReadCount_normalised/{{experiment}}.{{splice}}_{{prefix}}.{{lineage}}_{{valid}}.{{type}}.{{tag}}.{{feature}}.docx",
        feat=feat_res[feat_res.feature_name==wildcards.feature].itertuples(),
        valid=VALID,
        )
    output:
        docx="diff_reports/features/{experiment}/{reference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.docx",
    params:
        protocol = lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["protocol"],
        control=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["control"],
        treat=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["treatment"],
        paired=lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["pairRep"]),
        descript= lambda wildcards: feature_descript(wildcards),
    threads: 1
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/differential.yaml"
    log:
        "feature_reports/{experiment}/{reference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.log",
    script:
        "../scripts/R/feature_report.R"


rule differential_report:
    input:
        docx="resources/templates/{c}_{f}.docx".format(
            f=str(config["differential_plots"]["word_docx"]["font_name"]),
            c=str(config["differential_plots"]["word_docx"]["font_colour"]),
        ),
        plots=lambda wildcards: list(dict.fromkeys(expand(
            "diff_reports/features/{{experiment}}/{exp.reference}/{{pair}}.{{spikein}}_{{normaliser}}ReadCount_normalised/{{experiment}}.{{splice}}_{{prefix}}.{{lineage}}_{valid}.custom-{feat.prefix_md5}.{{tag}}.{feat.feature_name}.docx",
            exp=results[results.experiment==wildcards.experiment].itertuples(),
            feat=feat_res.itertuples(),
            valid=VALID,
        )))
    output:
        docx="diff_reports/experiment_reports/{experiment}/{experiment}.{lineage}.{tag}.{pair}.{spikein}.{normaliser}_normalised.{splice}_{prefix}.differential_report.docx"
    threads: 1
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/differential.yaml"
    log:
        "diff_reports/{splice}_{prefix}/{experiment}/{lineage}.{tag}.{pair}.{spikein}.{normaliser}_normalised/{experiment}.{splice}_{prefix}.differential_report.log"
    script:
        "../scripts/R/differential_report.R"
            
