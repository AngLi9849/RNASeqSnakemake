rule rmats_prep:
    input:
        bam="star/{sample}/{unit}/{reference}/AllAligned.sortedByCoord.out.bam",
        gtf="resources/annotations/{reference}/genome.gtf",
    output:
        sample="rmats/prep/{sample}/{unit}/{reference}/sample.txt",
        prep=directory("rmats/prep/{sample}/{unit}/{reference}/")
    params:
        read_length=lambda wildcards: samples.loc[wildcards.sample].loc[wildcards.unit,"readlen"],
        read_paired=lambda wildcards: "paired" if is_paired_end(wildcards.sample) else "single",
        temp="rmats/prep/{sample}/{unit}/{reference}/temp",
    threads: 4
    resources:
        mem=lambda wildcards, input: (str((input.size//3000000000)+4) + "G"),
        rmem=lambda wildcards, input: (str((input.size//6000000000)+4) + "G"),
    conda: 
        "../envs/rmats.yaml"
    shell:
        """
        awk 'BEGIN {{ print "{input.bam}"}}' > {output.sample} &&
        python $(conda info | awk 'match($0,/active env location : ([^ ]*)/, a) {{print a[1]"/rMATS/rmats.py" }}') \
          --b1 {output.sample} \
          --gtf {input.gtf} \
          -t {params.read_paired} \
          --readLength {params.read_length} \
          --variable-read-length \
          --nthread {threads} \
          --od {output.prep} \
          --tmp {params.temp} \
          --task prep
        """ 

rule rmats_post:
    input:
        control_bam=lambda wildcards: expand(
          "star/{control.sample_name}/{control.unit_name}/{{reference}}/AllAligned.sortedByCoord.out.bam",
          control=get_experiment_controls(wildcards.experiment).itertuples(),
        ),
        treat_bam=lambda wildcards: expand(
          "star/{treat.sample_name}/{treat.unit_name}/{{reference}}/AllAligned.sortedByCoord.out.bam",
          treat=get_experiment_treatments(wildcards.experiment).itertuples(),
        ),
        prep=lambda wildcards: expand(
          "rmats/prep/{sample.sample_name}/{sample.unit_name}/{{reference}}/temp",
          sample=get_experiment_samples(wildcards.experiment).itertuples(),
        ),
        gtf="resources/annotations/{reference}/genome.gtf",
    output:
        control="rmats/post/{experiment}/{reference}/temp/control.txt",
        treat="rmats/post/{experiment}/{reference}/temp/treat.txt",
        temp=directory("rmats/post/{experiment}/{reference}/temp/"),
#        post=directory("rmats/prep/{experiment}/{reference}/")
    params:
#        read_length=lambda wildcards: samples.loc[wildcards.sample].loc[wildcards.unit,"readlen"],
#        read_paired=lambda wildcards: "paired" if is_paired_end(wildcards.sample) else "single",
        temp="rmats/post/{experiment}/{reference}/temp",
    threads: 6
    resources:
        mem=lambda wildcards, input: (str((input.size//6000000000)+4) + "G"),
        rmem=lambda wildcards, input: (str((input.size//12000000000)+4) + "G"),
    conda:
        "../envs/rmats.yaml"
    shell:
        """
        for i in {input.control_bam} ; do echo "$i" >> {output.control}; done &&
        for i in {input.treat_bam} ; do echo "$i" >> {output.treat}; done &&
        for i in {input.prep} ;
          do for f in "$i"/* ; 
             do cp "$f" {params.temp}/$(awk -v f="$f" 'BEGIN{{ match(f,/rmats\/prep\/([^\/]*)\/([^\/]*)\//,a) ; print a[1]"_"a[2] }}')_$(awk -v f="$f" 'BEGIN{{ match(f,/([^\/]*)$/,a) ; print a[1] }}') ;
          done ;
        done    
#        do python $(conda info | awk 'match($0,/active env location : ([^ ]*)/, a) {{print a[1]"/rMATS/cp_with_prefix.py" }}') \
#         $(awk -v i="$i" 'BEGIN{{ match(i,/rmats\/prep\/([^\/]*)\/([^\/]*)\//,a) ; print a[1]"_"a[2] }}') \
#         $(pwd)/{params.temp} \
#         $(ls -d "$i"/*) \
#        ; done
        """
