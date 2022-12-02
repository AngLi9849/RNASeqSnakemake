rule heatmap_data:
    input:
        sense_mx = lambda wildcards: list(dict.fromkeys(expand("norm_mx/{sample.sample_name}/{sample.unit_name}/{{reference}}/{{splice}}{{prefix}}.{sample.stranded}/{{lineage}}_{{valid}}.plot-{plot_md5}.{{tag}}.{{feature}}.{{read}}.sense.{{mean}}_sum_matrix.gz",
            sample=results[results.experiment==wildcards.experiment].itertuples(),
            plot_md5=features.loc[wildcards.feature,"plot_md5"],
        ))),
        antisense_mx=lambda wildcards: [] if not features.loc[wildcards.feature,"antisense"] else list(dict.fromkeys(expand("norm_mx/{sample.sample_name}/{sample.unit_name}/{{reference}}/{{splice}}{{prefix}}.{sample.stranded}/{{lineage}}_{{valid}}.plot-{plot_md5}.{{tag}}.{{feature}}.{{read}}.antisense.{{mean}}_sum_matrix.gz",
            sample=results[results.experiment==wildcards.experiment].itertuples(),
            plot_md5=features.loc[wildcards.feature,"plot_md5"],
        ))),
        size_table=lambda w: "deseq2/{norm_group}/{{reference}}/All{{prefix}}.{{lineage}}_{{valid}}.{norm_type}.{{normaliser}}.{norm_read}.Count.{{spikein}}_{{pair}}.scale_factors.tsv".format(
            norm_type= ("custom-" + str(features.loc[w.normaliser,"prefix_md5"])) if (w.normaliser in features["feature_name"].tolist()) else "gtf",
            norm_group=experiments.loc[w.experiment,"group_name"],
        ),
    output:
        heat_data="heat_data/{experiment}/{reference}/differential_{read}.{count}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{mean}_heat_data.tab",
    resources:
        mem="48G",
        rmem="32G",
    params:
        control=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["control"],
        treat=lambda wildcards: experiments.loc[wildcards.experiment].squeeze(axis=0)["treatment"],
        sample_table=samples,
        bef_bin=lambda wildcards : features.loc[wildcards.feature,"bef_bin"],
        main_bin=lambda wildcards : features.loc[wildcards.feature,"bin_n"],
        plotbef_bin=lambda wildcards : features.loc[wildcards.feature,"plotbef_bin"],
        plotaft_bin=lambda wildcards : features.loc[wildcards.feature,"plotaft_bin"],
    log:
        "logs/heatmap_data/{experiment}/{reference}/differential_{read}.{count}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised_{mean}/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.log",
    threads: 1
    script:
        "../scripts/py/heatmap_data.py"


rule meta_plot_data:
    input:
        lfc="differential/{experiment}/{reference}/differential_{read}.{count}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.lfc.tab",
        levels="differential/{experiment}/{reference}/differential_{read}.{count}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.levels.tab",
        counts="differential/{experiment}/{reference}/differential_{read}.{count}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.counts.tab",
        sense_mx = lambda wildcards: list(dict.fromkeys(expand("norm_mx/{sample.sample_name}/{sample.unit_name}/{{reference}}/{{splice}}{{prefix}}.{sample.stranded}/{{lineage}}_{{valid}}.plot-{plot_md5}.{{tag}}.{{feature}}.{{read}}.sense.{{mean}}_sum_matrix.gz",
            sample=results[results.experiment==wildcards.experiment].itertuples(),
            plot_md5=features.loc[wildcards.feature,"plot_md5"],
        ))),
        antisense_mx=lambda wildcards: [] if not features.loc[wildcards.feature,"antisense"] else list(dict.fromkeys(expand("norm_mx/{sample.sample_name}/{sample.unit_name}/{{reference}}/{{splice}}{{prefix}}.{sample.stranded}/{{lineage}}_{{valid}}.plot-{plot_md5}.{{tag}}.{{feature}}.{{read}}.antisense.{{mean}}_sum_matrix.gz",
            sample=results[results.experiment==wildcards.experiment].itertuples(),
            plot_md5=features.loc[wildcards.feature,"plot_md5"],
        ))),
        size_table=lambda w: "deseq2/{norm_group}/{{reference}}/All{{prefix}}.{{lineage}}_{{valid}}.{norm_type}.{{normaliser}}.{norm_read}.Count.{{spikein}}_{{pair}}.scale_factors.tsv".format(
            norm_type= ("custom-" + str(features.loc[w.normaliser,"prefix_md5"])) if (w.normaliser in features["feature_name"].tolist()) else "gtf",
            norm_group=experiments.loc[w.experiment,"group_name"],
        ),
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
        sig_bg=lambda wildcards : "featurecounts/{norm_group}/{{reference}}/{{splice}}{{prefix}}.{{lineage}}_{{valid}}.{{type}}.{{tag}}.{{feature}}.{{read}}.plot-{plot_md5}.sig2bg.tab".format(
            norm_group=experiments.loc[wildcards.experiment,"group_name"],
            plot_md5=features.loc[wildcards.feature,"plot_md5"],
        ),
    output:
        mx_data="meta_data/{experiment}/{reference}/differential_{read}.{count}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.{mean}_{norm}.mx_data.tab",
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
        base=lambda wildcards : features.loc[wildcards.feature,"group"],
        start_name = lambda wildcards : str(features.loc[wildcards.feature,"strt_nm"]),
        end_name = lambda wildcards : str(features.loc[wildcards.feature,"end_nm"]),
    resources:
        mem="48G",
        rmem="32G",
    conda:
        "../envs/differential.yaml"
    log:
        "logs/meta_plot_data/{experiment}/{reference}/differential_{read}.{count}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised.{mean}_{norm}/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.log",
    threads: 1
    script:
        "../scripts/R/meta_plot_data.R"



rule differential_plots:
    input:
        lfc="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.lfc.tab",
        levels="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.levels.tab",
        counts="differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.counts.tab",
        mx_data=lambda wildcards: [] if wildcards.difference=="splicing_ratio" else "meta_data/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.{mean}_{norm}.mx_data.tab",
        heat_data = "heat_data/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.{mean}_heat_data.tab",
        size_table=lambda w: "deseq2/{norm_group}/{{reference}}/All{{prefix}}.{{lineage}}_{{valid}}.{norm_type}.{{normaliser}}.{norm_read}.Count.{{spikein}}_{{pair}}.scale_factors.tsv".format(
            norm_type= ("custom-" + str(features.loc[w.normaliser,"prefix_md5"])) if (w.normaliser in features["feature_name"].tolist()) else "gtf",
            norm_group=experiments.loc[w.experiment,"group_name"],
        ),
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
        sig_bg=lambda wildcards : "featurecounts/{norm_group}/{{reference}}/{{splice}}{{prefix}}.{{lineage}}_{{valid}}.{{type}}.{{tag}}.{{feature}}.{{read}}.plot-{plot_md5}.sig2bg.tab".format(
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
        docx="diff_plots/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.docx",
        rdata="diff_plots/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.Rdata",
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
        bar_data="diff_plots/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.bar.tab",
        pie_data="diff_plots/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.pie.tab",
        violin_data="diff_plots/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.violin.tab",
    resources:
        mem="48G",
        rmem="32G",
    conda:
        "../envs/differential.yaml"
    log:
        "logs/differential/{experiment}/{reference}/differential_{difference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised.{mean}_{norm}/{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.plots.log",
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
            "diff_plots/{{experiment}}/{{reference}}/differential_{feat.diff}/{{pair}}.{{spikein}}_{{normaliser}}.{norm_read}.Count_normalised.{{mean}}_{{norm}}/{{experiment}}.{{splice}}_{{prefix}}.{{lineage}}_{{valid}}.{{type}}.{{tag}}.{{feature}}.docx",
        feat=feat_res[feat_res.feature_name==wildcards.feature].itertuples(),
        valid=VALID,
        )
    output:
        docx="diff_reports/features/{experiment}/{reference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.docx",
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
        "feature_reports/{experiment}/{reference}/{pair}.{spikein}_{normaliser}.{norm_read}.Count_normalised.{mean}_{norm}/{experiment}.{splice}_{prefix}.{lineage}_{valid}.{type}.{tag}.{feature}.{read}.log",
    script:
        "../scripts/R/feature_report.R"


rule differential_report:
    input:
        docx="resources/templates/{c}_{f}.docx".format(
            f=str(config["differential_plots"]["word_docx"]["font_name"]),
            c=str(config["differential_plots"]["word_docx"]["font_colour"]),
        ),
        plots=lambda wildcards: list(dict.fromkeys(expand(
            "diff_reports/features/{{experiment}}/{exp.reference}/{{pair}}.{{spikein}}_{{normaliser}}.{norm_read}.Count_normalised.{{mean}}_{{norm}}/{{experiment}}.{{splice}}_{{prefix}}.{{lineage}}_{valid}.custom-{feat.prefix_md5}.{{tag}}.{feat.feature_name}.docx",
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
            
