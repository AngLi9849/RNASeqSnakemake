rule group_data:
    input:
        lfc=lambda wildcards: list(dict.fromkeys(expand(
            "differential/{exp.experiment}/{exp.reference}/differential_{exp.difference}/{exp.paired}.{exp.normaliser}_{exp.norm_feat}ReadCount_normalised/{exp.splice}_{exp.demulti}{exp.dedup}.{exp.diff_lineage}_{exp.valid}.{exp.type}.{{tag}}.{exp.feature}.lfc.tab",
            exp=groups.loc[wildcards.md5].itertuples(),
        ))),
    params:
        titles=lambda wildcards: groups.loc[wildcars.md5,"title"].tolist()
        group_title=lambda wildcards: str(group_config.loc[wildcards.md5,"group_title"])
    output:
        group_lfc="differential/group_reports/Group.{group}.{title}.{md5}.{tag}.docx",
    log:
        logs="differential/Group.{group}.{title}.{md5}.{tag}.log",
    conda:
        "../envs/differential.yaml"
    script:
        "../scripts/R/group_data.R"

rule group_plots:
    input:
        lfc="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.lfc.tab",
        size_table=lambda w: "deseq2/{norm_group}/{{reference}}/All{{prefix}}.{{lineage}}_{{valid}}.{norm_type}.{{normaliser}}ReadCount.{{spikein}}_{{pair}}.scale_factors.tsv".format(
        base_bed=lambda wildcards: ("resources/annotations/{source}/{{lineage}}.{s}.{{valid}}_{{tag}}.{f}.bed".format(
            source=  str( get_sample_source(wildcards.experiment) ),
            s=features.loc[features.loc[wildcards.feature,"group"],"type"] if (features.loc[wildcards.feature,"group"] in features["feature_name"].tolist()) else "gtf",
            f=features.loc[wildcards.feature,"group"]
        ) ), 
        genetab=lambda w: "resources/annotations/{source}/genome.gtf.{{tag}}_gene_info.tab".format(
            source= str( get_sample_source(w.experiment) ),
        ),
        bed=lambda w: "resources/annotations/{source}/{{lineage}}.{{type}}.{{valid}}_{{tag}}.{{feature}}.bed".format(
            source=  str( get_sample_source(w.experiment) ),
        ),
        docx="resources/templates/{c}_{f}.docx".format(
            f=str(config["differential_plots"]["word_docx"]["font_name"]),
            c=str(config["differential_plots"]["word_docx"]["font_colour"]),
        ),
    output:
        docx="diff_plots/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.docx",
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
        plotbef_len=lambda wildcards : features.loc[wildcards.feature,"plotbef"],
        plotaft_len=lambda wildcards : features.loc[wildcards.feature,"plotaft"],
        plotbef_bin=lambda wildcards : features.loc[wildcards.feature,"plotbef_bin"],
        plotaft_bin=lambda wildcards : features.loc[wildcards.feature,"plotaft_bin"],        
        base = lambda wildcards : features.loc[wildcards.feature,"group"],
        base_feat = lambda wildcards : features.loc[wildcards.feature,"feature"],
        is_antisense=lambda wildcards : features.loc[wildcards.feature,"sense_dir"],
        main_int = lambda wildcards : str(features.loc[wildcards.feature,"is_main_int"]),
        start_name = lambda wildcards : str(features.loc[wildcards.feature,"strt_nm"]),
        end_name = lambda wildcards : str(features.loc[wildcards.feature,"end_nm"]),
        bar_data="diff_plots/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.bar.tab",
        pie_data="diff_plots/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.pie.tab",
        violin_data="diff_plots/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.violin.tab",
    resources:
        mem="48G",
        rmem="32G",
    conda:
        "../envs/differential.yaml"
    log:
        "logs/differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}ReadCount_normalised.{mean}_{norm}/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.plots.log",
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
            "diff_plots/{{experiment}}/{{reference}}/differential_{feat.diff}/{{pair}}.{{spikein}}_{{normaliser}}ReadCount_normalised.{{mean}}_{{norm}}/{{experiment}}.{{splice}}_{{prefix}}.{{lineage}}_{{valid}}.{{type}}.{{tag}}.{{feature}}.docx",
        feat=feat_res[feat_res.feature_name==wildcards.feature].itertuples(),
        valid=VALID,
        )
    output:
        docx="diff_reports/features/{experiment}/{reference}/{pair}.{spikein}_{normaliser}ReadCount_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.docx",
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
        "feature_reports/{experiment}/{reference}/{pair}.{spikein}_{normaliser}ReadCount_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.log",
    script:
        "../scripts/R/feature_report.R"


rule differential_report:
    input:
        docx="resources/templates/{c}_{f}.docx".format(
            f=str(config["differential_plots"]["word_docx"]["font_name"]),
            c=str(config["differential_plots"]["word_docx"]["font_colour"]),
        ),
        plots=lambda wildcards: list(dict.fromkeys(expand(
            "diff_reports/features/{{experiment}}/{exp.reference}/{{pair}}.{{spikein}}_{{normaliser}}ReadCount_normalised.{{mean}}_{{norm}}/{{experiment}}.{{splice}}_{{prefix}}.{{lineage}}_{valid}.custom-{feat.prefix_md5}.{{tag}}.{feat.feature_name}.docx",
            exp=results[results.experiment==wildcards.experiment].itertuples(),
            feat=feat_res.itertuples(),
            valid=VALID,
        )))
    output:
        docx="diff_reports/experiment_reports/{experiment}.{lineage}.{tag}.{pair}.{spikein}.{normaliser}_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.differential_report.docx"
    threads: 1
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/differential.yaml"
    log:
        "diff_reports/{splice}_{prefix}/{experiment}/{lineage}.{tag}.{pair}.{spikein}.{normaliser}_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.differential_report.log"
    script:
        "../scripts/R/differential_report.R"
            
