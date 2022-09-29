rule meta_plot_data:
    input:
        lfc="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.lfc.tab",
        levels="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.levels.tab",
        counts="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.counts.tab",
        sense_mx = lambda wildcards: list(dict.fromkeys(expand("norm_mx/{sample.sample_name}/{sample.unit_name}/{{reference}}/{{splice}}{{prefix}}.{sample.stranded}/{{lineage}}_{{valid}}.plot-{plot_md5}.{{tag}}.{{feature}}.sense.{norm}_matrix.gz",
            sample=results[results.experiment==wildcards.experiment].itertuples(),
            plot_md5=features.loc[wildcards.feature,"plot_md5"],
            norm="norm" if config["metagene"]["norm_per_gene"] else "sum",
        ))),
        antisense_mx=lambda wildcards: [] if not features.loc[wildcards.feature,"antisense"] else list(dict.fromkeys(expand("norm_mx/{sample.sample_name}/{sample.unit_name}/{{reference}}/{{splice}}{{prefix}}.{sample.stranded}/{{lineage}}_{{valid}}.plot-{plot_md5}.{{tag}}.{{feature}}.antisense.{norm}_matrix.gz",
            sample=results[results.experiment==wildcards.experiment].itertuples(),
            plot_md5=features.loc[wildcards.feature,"plot_md5"],
            norm="norm" if config["metagene"]["norm_per_gene"] else "sum",
        ))),
        size_table=lambda w: "deseq2/{norm_group}/{{reference}}/All{{prefix}}.{{lineage}}_{{valid}}.{norm_type}.{{normaliser}}ReadCount.{{spikein}}_{{pair}}.scale_factors.tsv".format(
            norm_type= ("custom-" + str(features.loc[w.normaliser,"prefix_md5"])) if (w.normaliser in features["feature_name"].tolist()) else "gtf",
            norm_group=experiments.loc[w.experiment,"group_name"],
        ),
        sig_bg=lambda wildcards : "featurecounts/{norm_group}/{{reference}}/{{splice}}{{prefix}}.{{lineage}}_{{valid}}.{{type}}.{{tag}}.{{feature}}.plot-{plot_md5}.sig2bg.tab".format(
            norm_group=experiments.loc[wildcards.experiment,"group_name"],
            plot_md5=features.loc[wildcards.feature,"plot_md5"],
        ),
    output:
        mx_data="meta_data/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.mx_data.tab",
    params:
        genesets=lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["gene_sets"]).split(","),
        goi= lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["GOI"]).split(","),
        protocol = lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["protocol"],
        control=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["control"],
        treat=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["treatment"],
        paired=lambda wildcards: str(experiments.loc[wildcards.experiment].squeeze(axis=0)["pairRep"]),
        descript= lambda wildcards: feature_descript(wildcards),
        samples= lambda wildcards: list(dict.fromkeys(expand("{sample.sample_name}",
            sample=results[results.experiment==wildcards.experiment].itertuples(),
        ))),
        sig=lambda wildcards : features.loc[wildcards.feature,"s2b_min"],
        bg=lambda wildcards : features.loc[wildcards.feature,"b2s_min"],
        section=lambda wildcards : features.loc[wildcards.feature,"section"],
        len_bef=lambda wildcards : features.loc[wildcards.feature,"len_bef"],
        len_aft=lambda wildcards : features.loc[wildcards.feature,"len_aft"],
        bef_bin=lambda wildcards : features.loc[wildcards.feature,"bef_bin"],
        main_bin=lambda wildcards : features.loc[wildcards.feature,"bin_n"],
        plotbef_bin=lambda wildcards : features.loc[wildcards.feature,"plotbef_bin"],
        plotaft_bin=lambda wildcards : features.loc[wildcards.feature,"plotaft_bin"],
        base=lambda wildcards : features.loc[wildcards.feature,"feature"],
    resources:
        mem="48G",
        rmem="32G",
    conda:
        "../envs/differential.yaml"
    log:
        "logs/meta_plot_data/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.log",
    threads: 1
    script:
        "../scripts/R/meta_plot_data.R"



rule differential_plots:
    input:
        lfc="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.lfc.tab",
        levels="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.levels.tab",
        counts="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.counts.tab",
        mx_data=lambda wildcards: [] if not wildcards.difference=="expression" else "meta_data/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.mx_data.tab",
        size_table=lambda w: "deseq2/{norm_group}/{{reference}}/All{{prefix}}.{{lineage}}_{{valid}}.{norm_type}.{{normaliser}}ReadCount.{{spikein}}_{{pair}}.scale_factors.tsv".format(
            norm_type= ("custom-" + str(features.loc[w.normaliser,"prefix_md5"])) if (w.normaliser in features["feature_name"].tolist()) else "gtf",
            norm_group=experiments.loc[w.experiment,"group_name"],
        ),
        bed=lambda w: "resources/annotations/{source}/{{lineage}}.{{type}}.{{valid}}_{{tag}}.{{feature}}.bed".format(
            source=  str( get_sample_source(w.experiment) ),
        ),
        sig_bg=lambda wildcards : "featurecounts/{norm_group}/{{reference}}/{{splice}}{{prefix}}.{{lineage}}_{{valid}}.{{type}}.{{tag}}.{{feature}}.plot-{plot_md5}.sig2bg.tab".format(
            norm_group=experiments.loc[wildcards.experiment,"group_name"],
            plot_md5=features.loc[wildcards.feature,"plot_md5"],
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
        samples= lambda wildcards: list(dict.fromkeys(expand("{sample.sample_name}",
            sample=results[results.experiment==wildcards.experiment].itertuples(),
        ))), 
        sig=lambda wildcards : features.loc[wildcards.feature,"s2b_min"],
        bg=lambda wildcards : features.loc[wildcards.feature,"b2s_min"],
        section=lambda wildcards : features.loc[wildcards.feature,"section"],
        len_bef=lambda wildcards : features.loc[wildcards.feature,"len_bef"],
        len_aft=lambda wildcards : features.loc[wildcards.feature,"len_aft"],
        bef_bin=lambda wildcards : features.loc[wildcards.feature,"bef_bin"],
        main_bin=lambda wildcards : features.loc[wildcards.feature,"bin_n"],
        plotbef_bin=lambda wildcards : features.loc[wildcards.feature,"plotbef_bin"],
        plotaft_bin=lambda wildcards : features.loc[wildcards.feature,"plotaft_bin"],        
        base=lambda wildcards : features.loc[wildcards.feature,"feature"],
        is_antisense=lambda wildcards : features.loc[wildcards.feature,"sense_dir"]
    resources:
        mem="486G",
        rmem="32G",
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
            
