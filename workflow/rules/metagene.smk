rule feature_metagene_annotations:
    input:
        bed=lambda w: "{{prefix}}.custom-{id}.{{type}}.{{feature}}.{{sense}}.main.bed".format(
            id = features.loc[w.feature,"prefix_md5"],
        ),
    output:
        before="{prefix}.plot-{md5}.{type}.{feature}.{sense}.before.bed",
        after="{prefix}.plot-{md5}.{type}.{feature}.{sense}.after.bed",
    threads: 1 
    params:
        strand=lambda w: "-s" if w.strand == "stranded" else "",
        distinct= lambda w: "distinct" if w.strand == "stranded" else "collapse",
        before=lambda w: features.loc[w.feature,"plotbef"],
        after=lambda w: features.loc[w.feature,"plotaft"],
    conda:
        "../envs/bedtools.yaml"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/metagene/{prefix}_plot-{md5}.{type}.{feature}.{sense}.non-overlap.{strand}.log",
    shell:
        """
        sort -k6,6 -k1,1 -k8,8 -k2,2n -k3,3n {input.bed} 
        awk -F'\\t' -v OFS='\\t' -v id='' '
          FNR==NR {{
            length[$8] += $5 ; 
            five[$8]=(five[$8]<=$2)?five[$8]:$2 ;
            three[$8]=(three[$8]>=$3)?three[$8]:$3 ;
          }}
          FNR < NR && $6=="+" {{
            if ( id != $8 ) {{
              l=length[$8] ;
              bef=("{params.before}" ~ /\\..*[1-9]/)? int(l*{params.before}) : {params.before} ;
              aft=("{params.after}" ~ /\\..*[1-9]/)? int(l*{params.after}) : {params.after} ;
              if (($6=="+" && "{wildcards.sense}" == "sense") || ($6=="-" && "{wildcards.sense}" == "antisense")) {{ 
                $2 = ( five[$8] >= bef ) ? ( five[$8] - bef ) : 0 ;
                $3 = five[$8] ;
                print >> "{output.before}" ;
                $2 = three[$8] ;
                $3 = three[$8] + aft ;
                print >> "{output.after}" ;
              }} else {{
                $2 = (five[$8] >= aft)? ( five[$8] - aft ) : 0   ;
                $3 = five[$8] ;
                print >> "{output.after}" ;
                $2 = three[$8] ; 
                $3 = three[$8] + bef
                print >> "{output.before}" ;
              }} ;
              id = $8 ;
            }}
          }}' - {input.bed} 
        """


rule compute_matrix_promptTSS_stranded:
    input:        
        fwd_bed=lambda wildcards: (str(get_annotation(experiments.loc[wildcards.experiment,"sample_source"]) + ".{wildcards.biotype}_promptTSS.non_overlap.fwd.bed"),
        bigwig=expand(
            "results/{experiment}/bigwig/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{sample}_{unit}.fwd_{splice}.bigwig",
    output:
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_promptTSS/{sample}_{unit}.sense.fwd.matrix.gz",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_promptTSS/{sample}_{unit}.sense.rev.matrix.gz",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_promptTSS/{sample}_{unit}.antisense.fwd.matrix.gz",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_promptTSS/{sample}_{unit}.antisense.rev.matrix.gz",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_promptTSS/{sample}_{unit}_sum.sense.matrix.tab",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_promptTSS/{sample}_{unit}_per_gene.sense.matrix.tab",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_promptTSS/{sample}_{unit}_sum.antisense.matrix.tab",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_promptTSS/{sample}_{unit}_per_gene.antisense.matrix.tab",
    log:
        "logs/metagene/by_{normaliser}_{counts}_{splice}{prefix}/{biotype}_promptTSS/{sample}_{unit}.matrix.log",
    params:
        bin_size=config["metagene"]["promptTSS"]["bin_size"],
    threads: 4
    resources:
        mem="10G",
        rmem="6G",
    conda:
        "../envs/deeptools.yaml",
    shell:
        """
        computeMatrix reference-point -S {input.fwd_bigwig} -R {input.fwd_bed} -p {threads} -b {params.prompt} -a {params.tss} --referencePoint TSS --binSize {params.bin_size} --averageTypeBins mean --sortRegions descend --sortUsing region_length -o {output[0]} && 
        computeMatrix reference-point -S {input.rev_bigwig} -R {input.rev_bed} -p {threads} -b {params.prompt} -a {params.tss} --referencePoint TSS --binSize {params.bin_size} --averageTypeBins mean --sortRegions descend --sortUsing region_length -o {output[1]} &&
        computeMatrix reference-point -S {input.rev_bigwig} -R {input.fwd_bed} -p {threads} -b {params.prompt} -a {params.tss} --referencePoint TSS --binSize {params.bin_size} --averageTypeBins mean --sortRegions descend --sortUsing region_length -o {output[2]} &&
        computeMatrix reference-point -S {input.fwd_bigwig} -R {input.rev_bed} -p {threads} -b {params.prompt} -a {params.tss} --referencePoint TSS --binSize {params.bin_size} --averageTypeBins mean --sortRegions descend --sortUsing region_length -o {output[3]} &&
        zcat {output[0]} {output[1]} | 
        awk -F'\\t'   -v OFS='\\t' 'NF>1&& $0 !~ "nan" {{$2=sqrt(($3-$2)^2);print}}' - |
        sort -rnk2,2 - |
        cut -f 1,2,4,7- - |
        sed '1 i\\{wildcards.sample}_sense' - > {output[4]} &&
        zcat {output[2]} {output[3]} |
        awk -F'\\t' -v OFS='\\t' 'NF>1&& $0 !~ "nan" {{$2=sqrt(($3-$2)^2);print}}' - |
        sort -rnk2,2 - |
        cut -f 1,2,4,7- - |
        sed '1 i\\{wildcards.sample}_antisense' - > {output[6]} &&
        paste {output[4]} {output[6]} | 
        awk -v OFS='\\t' 'FNR==NR&&FNR>1{{for (i=4; i<=NF; i++) sum[FNR]+=$i;size[FNR]=sum[FNR]/(NF-3);total+=(sum[FNR]/(FNR-1)); next}} FNR<NR&&FNR==1{{print $0;next}} FNR<NR&&FNR>1{{for(i=4;i<=NF;i++) $i=(size[FNR]>0?$i/size[FNR]*total:0);print}}' - {output[4]} > {output[5]} &&
        paste {output[4]} {output[6]} |
        awk -v OFS='\\t' 'FNR==NR&&FNR>1{{for (i=4; i<=NF; i++) sum[FNR]+=$i;size[FNR]=sum[FNR]/(NF-3);total+=(sum[FNR]/(FNR-1)); next}} FNR<NR&&FNR==1{{print $0;next}} FNR<NR&&FNR>1{{for(i=4;i<=NF;i++) $i=(size[FNR]>0?$i/size[FNR]*total:0);print}}' - {output[6]} > {output[7]}
        """
    

rule compute_matrix_TTSpostgene_stranded:
    input:
        fwd_bed=lambda wildcards: (str(get_annotation(experiments.loc[wildcards.experiment,"sample_source"]) + ".{wildcards.biotype}_TTSpostgene.non_overlap.fwd.bed"),
        rev_bed=lambda wildcards: (str(get_annotation(experiments.loc[wildcards.experiment,"sample_source"]) + ".{wildcards.biotype}_TTSpostgene.non_overlap.rev.bed"),
        fwd_bigwig="{experiment}/bigwig/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{sample}_{unit}.fwd_{splice}.bigwig",
        rev_bigwig="{experiment}/bigwig/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{sample}_{unit}.rev_{splice}.bigwig",
    output:
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_TTSpostgene/{sample}_{unit}.fwd.matrix.gz",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_TTSpostgene/{sample}_{unit}.rev.matrix.gz",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_TTSpostgene/{sample}_{unit}_sum.stranded.matrix.tab",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_TTSpostgene/{sample}_{unit}_per_gene.stranded.matrix.tab",
    log:
        "logs/metagene/by_{normaliser}_{counts}_{splice}{prefix}/{biotype}_TTSpostgene/{sample}_{unit}.matrix.log",
    params:
        tts=config["metagene"]["TTSpostgene"]["TTS_length"],
        postgene=config["metagene"]["TTSpostgene"]["postgene_length"],
        bin_size=config["metagene"]["TTSpostgene"]["bin_size"],
    threads: 4
    resources:
        mem="10G",
        rmem="6G",
    conda:
        "../envs/deeptools.yaml",
    shell:
        """
        computeMatrix reference-point -S {input.fwd_bigwig} -R {input.fwd_bed} -p {threads} -b {params.tts} -a {params.postgene} --referencePoint TES --binSize {params.bin_size} --averageTypeBins mean --sortRegions descend --sortUsing region_length -o {output[0]} &&
        computeMatrix reference-point -S {input.rev_bigwig} -R {input.rev_bed} -p {threads} -b {params.tts} -a {params.postgene} --referencePoint TES --binSize {params.bin_size} --averageTypeBins mean --sortRegions descend --sortUsing region_length -o {output[1]} &&
        zcat {output[0]} {output[1]} | 
        awk -F'\\t' -v OFS='\\t' 'NF>1&& $0 !~ "nan" {{$2=sqrt(($3-$2)^2);print}}' - |
        sort -rnk2,2 - |
        cut -f 1,2,4,7- - | 
        sed '1 i\\{wildcards.sample}' - > {output[2]} &&
        awk -v OFS='\\t' 'FNR==NR&&FNR>1{{for (i=4; i<=NF; i++) sum[FNR]+=$i;size[FNR]=sum[FNR]/(NF-3);total+=(sum[FNR]/(FNR-1)); next}} FNR<NR&&FNR==1{{print $0;next}} FNR<NR&&FNR>1{{for(i=4;i<=NF;i++) $i=(size[FNR]>0?$i/size[FNR]*total:0);print}}' {output[2]} {output[2]} > {output[3]}
        """


 
rule compute_matrix_genebody:
    input:
        bed=lambda wildcards: (str(get_annotation(experiments.loc[wildcards.experiment,"sample_source"]) + ".{wildcards.biotype}_genebody.non_overlap.{wildcards.strand}.bed"),
        bigwig="{experiment}/bigwig/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{sample}_{unit}.{strand}_{splice}.bigwig",
    output:
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}/{sample}_{unit}_genebody.{strand}.matrix.gz",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}/{sample}_{unit}_genebody_sum.{strand}.matrix.tab",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}/{sample}_{unit}_genebody_per_gene.{strand}.matrix.tab",
    params:
        before=config["metagene"]["genebody"]["before_TSS"],
        after=config["metagene"]["genebody"]["after_TTS"],
        start=config["metagene"]["genebody"]["TSS_unscaled"],
        end=config["metagene"]["genebody"]["TTS_unscaled"],
        bin_size=config["metagene"]["genebody"]["bin_size"],
        body_length=config["metagene"]["genebody"]["gene_body"],
    wildcard_constraints:
        strand=r"fwd|rev|unstranded",
    resources:
        mem="12G",
        rmem="8G",
    log:
        "logs/metagene/by_{normaliser}_{counts}_{splice}_{prefix}/{sample}_{unit}_{biotype}_genebody.{strand}.matrix.log",
    conda:
        "../envs/deeptools.yaml",
    threads: 2 
    shell:
        """
        computeMatrix scale-regions -S {input.bigwig} -R {input.bed} -p {threads} -b {params.before} -a {params.after} --unscaled5prime {params.start} --unscaled3prime {params.end} --binSize {params.bin_size} --averageTypeBins mean --regionBodyLength {params.body_length} --sortRegions descend --sortUsing region_length -o {output[0]} && 
        zcat {output[0]} |
        awk -F'\\t' -v OFS='\\t' 'NF>1&& $0 !~ "nan" {{$2=sqrt(($3-$2)^2);print}}' - |
        cut -f 1,2,4,7- - |
        sed '1 i\\{wildcards.sample}' - > {output[1]} &&
        awk -v OFS='\\t' 'FNR==NR&&FNR>1{{for (i=(({params.before}/{params.bin_size})+4); i<=((({params.before}+{params.body_length})/{params.bin_size})+3); i++) sum[FNR]+=$i;size[FNR]=sum[FNR]*{params.bin_size}/{params.body_length};total+=(sum[FNR]/(FNR-1)); next}} FNR<NR&&FNR==1{{print $0;next}} FNR<NR&&FNR>1{{for(i=4;i<=NF;i++) $i=(size[FNR]>0?$i/size[FNR]*total:0);print}}' {output[1]} {output[1]} > {output[2]}
        """
       

rule combine_stranded_genebody_matrices:
    input:
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}/{sample}_{unit}_genebody_{score}.fwd.matrix.tab",
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}/{sample}_{unit}_genebody_{score}.rev.matrix.tab",
    output:
        "{experiment}/meta_matrices/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/{biotype}/{sample}_{unit}_genebody_{score}.stranded.matrix.tab",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/meta_matrices/{biotype}_by_{normaliser}_{counts}_{splice}_{prefix}/{sample}_{unit}_genebody_{score}_combine.log"
    shell:
        """
        awk -v OFS='\t' 'FNR>1{{print}}' {input[0]} {input[1]} |  
        sort -rnk2,2 - |
        sed '1 i\\{wildcards.sample}' - > {output[0]} 2>{log}
        """ 

rule plot_promptTSS_stranded_meta_profile:
    input:
        sense_matrices=expand(
            "{{experiment}}/meta_matrices/{{splice}}_{{prefix}}_normalised_by_{{normaliser}}_{{counts}}/{{biotype}}_promptTSS/{sample.sample_name}_{sample.unit_name}_{{score}}.sense.matrix.tab",sample=samples.itertuples(),
        ),
        antisense_matrices=expand(
            "{{experiment}}/meta_matrices/{{splice}}_{{prefix}}_normalised_by_{{normaliser}}_{{counts}}/{{biotype}}_promptTSS/{sample.sample_name}_{sample.unit_name}_{{score}}.antisense.matrix.tab",sample=samples.itertuples(),
        ),
        size_table="{experiment}/meta_profiles/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_genebody/{splice}_{biotype}_genebody_normalised_{score}.stranded.sizes.tab",
        sample_table="config/samples.tsv",
    output:
        meta_profile_png="{experiment}/meta_profiles/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_promptTSS/{splice}_{biotype}_promptTSS_normalised_{score}.stranded.meta_profile.png",
        meta_profile_pdf="{experiment}/meta_profiles/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_promptTSS/{splice}_{biotype}_promptTSS_normalised_{score}.stranded.meta_profile.pdf",
    log:
        "logs/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/meta_graphs/{biotype}_promptTSS/{score}.stranded.meta_profile.log"
    threads: 1
    resources:
        mem=lambda wildcards, input: (str((input.size//10000000)+8) + "G"),
        rmem=lambda wildcards, input: (str((input.size//20000000)+8) + "G"),
    params:
        prompt=config["metagene"]["promptTSS"]["prompt_length"],
        tss=config["metagene"]["promptTSS"]["TSS_length"],
        bin_size=config["metagene"]["promptTSS"]["bin_size"],
        trim=config["metagene"]["promptTSS"]["anomaly_trim"],
        dir="{experiment}/meta_profiles/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_promptTSS",
        control=config["control_condition"],
        paired=config["paired_analysis"],
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/R/stranded_prompt_meta_profile.R"


rule plot_TTSpostgene_meta_profile:
    input:
        sense_matrices=expand(
            "{{experiment}}/meta_matrices/{{splice}}_{{prefix}}_normalised_by_{{normaliser}}_{{counts}}/{{biotype}}_TTSpostgene/{sample.sample_name}_{sample.unit_name}_{{score}}.{{strand}}.matrix.tab",sample=samples.itertuples(),
        ),
        size_table="{{experiment}}/meta_profiles/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_genebody/{splice}_{biotype}_genebody_normalised_{score}.{strand}.sizes.tab",
        sample_table="config/samples.tsv",
    output:
        meta_profile_png="{experiment}/meta_profiles/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_TTSpostgene/{splice}_{biotype}_TTSpostgene_normalised_{score}.{strand}.meta_profile.png",
        meta_profile_pdf="{experiment}/meta_profiles/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_TTSpostgene/{splice}_{biotype}_TTSpostgene_normalised_{score}.{strand}.meta_profile.pdf",
    log:
        "logs/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/meta_graphs/{biotype}_TTSpostgene/{score}.{strand}.meta_profile.log"
    threads: 1
    resources:
        mem=lambda wildcards, input: (str((input.size//60000000)+8) + "G"),
        rmem=lambda wildcards, input: (str((input.size//125000000)+8) + "G"),
    params:
        tts=config["metagene"]["TTSpostgene"]["TTS_length"],
        postgene=config["metagene"]["TTSpostgene"]["postgene_length"],
        bin_size=config["metagene"]["TTSpostgene"]["bin_size"],
        trim=config["metagene"]["TTSpostgene"]["anomaly_trim"],
        dir="{experiment}/meta_profiles/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_TTSpostgene",
        control=config["control_condition"],
        paired=config["paired_analysis"],
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/R/stranded_postgene_meta_profile.R"



rule plot_genebody_meta_profile:
    input:
        matrices=expand(
            "{{experiment}}/meta_matrices/{{splice}}_{{prefix}}_normalised_by_{{normaliser}}_{{counts}}/{{biotype}}/{sample.sample_name}_{sample.unit_name}_genebody_{{score}}.{{strand}}.matrix.tab",sample=samples.itertuples(),
        ),
        sample_table="config/samples.tsv",
    output:
        meta_profile_png="{experiment}/meta_profiles/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_genebody/{splice}_{biotype}_genebody_normalised_{score}.{strand}.meta_profile.png",
        meta_profile_pdf="{experiment}/meta_profiles/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_genebody/{splice}_{biotype}_genebody_normalised_{score}.{strand}.meta_profile.pdf",
        size_table="{experiment}/meta_profiles/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_genebody/{splice}_{biotype}_genebody_normalised_{score}.{strand}.sizes.tab",
    log:
        "logs/{splice}_{prefix}_normalised_by_{normaliser}_{counts}/meta_graphs/{biotype}/genebody_{score}.{strand}.meta_profile.log"
    threads: 1 
    resources:
        mem=lambda wildcards, input: (str((input.size//50000000)+10) + "G"),
        rmem=lambda wildcards, input: (str((input.size//100000000)+10) + "G"),
    params:
        before=config["metagene"]["genebody"]["before_TSS"],
        after=config["metagene"]["genebody"]["after_TTS"],
        start=config["metagene"]["genebody"]["TSS_unscaled"],
        end=config["metagene"]["genebody"]["TTS_unscaled"],
        body_length=config["metagene"]["genebody"]["gene_body"],
        bin_size=config["metagene"]["genebody"]["bin_size"],
        trim=config["metagene"]["genebody"]["anomaly_trim"],
        dir="{experiment}/meta_profiles/{splice}{prefix}_normalised_by_{normaliser}_{counts}/{biotype}_genebody",
        control=config["control_condition"],
        paired=config["paired_analysis"],
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/R/genebody_meta_profile.R"


