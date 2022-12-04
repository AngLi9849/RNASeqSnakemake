rule group_plots:
    input:
        lfc=lambda wildcards: expand(
            "differential/{exp.experiment}/{exp.reference}/differential_{exp.difference}/{exp.paired}.{exp.normaliser}_{exp.norm_feat}.{exp.norm_read}.Count_normalised/{exp.splice}_Aligned{exp.demulti}{exp.dedup}.{exp.diff_lineage}_{exp.valid}.{exp.type}.{{tag}}.{exp.feature}.lfc.tab",
            exp=groups.loc[wildcards.md5].itertuples(),
        ),
        docx="resources/templates/{c}_{f}.docx".format(
            f=str(config["differential_plots"]["word_docx"]["font_name"]),
            c=str(config["differential_plots"]["word_docx"]["font_colour"]),
        ),
    resources:
        mem="24G",
        rmem="16G",
    threads: 1
    params:
        titles=lambda wildcards: groups.loc[wildcards.md5,"title"].tolist(),
        descripts=lambda wildcards: groups.loc[wildcards.md5,"descript"].tolist(),
        group_title=lambda wildcards: str(group_config.loc[wildcards.md5,"group_title"]),
        genesets=lambda wildcards: str(group_config.loc[wildcards.md5,"gene_sets"]).split(","),
    output:
        docx="group_reports/Group.{group}.{title}.{md5}.{tag}.docx",
    log:
        logs="logs/Group.{group}.{title}.{md5}.{tag}.log",
    conda:
        "../envs/differential.yaml"
    script:
        "../scripts/R/group_data.R"

