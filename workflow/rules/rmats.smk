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
        libtype=lambda wildcards: get_rmats_sample_libtype(wildcards.sample,wildcards.unit), 
    threads: 2
    resources:
        mem=lambda wildcards, input: (str((input.size//1500000000)+4) + "G"),
        rmem=lambda wildcards, input: (str((input.size//3000000000)+4) + "G"),
    conda: 
        "../envs/rmats.yaml"
    shell:
        """
        readlink -f "{input.bam}" >> {output.sample} &&
        python $(conda info | awk 'match($0,/active env location : ([^ ]*)/, a) {{print a[1]"/rMATS/rmats.py" }}') \
          --b1 {output.sample} \
          --gtf {input.gtf} \
          -t {params.read_paired} \
          --readLength {params.read_length} \
          --libType {params.libtype} \
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
          "rmats/prep/{sample.sample_name}/{sample.unit_name}/{{reference}}",
          sample=get_experiment_samples(wildcards.experiment).itertuples(),
        ),
        gtf="resources/annotations/{reference}/genome.gtf",
    output:
        control="rmats/post/{experiment}/{reference}/temp/control.txt",
        treat="rmats/post/{experiment}/{reference}/temp/treat.txt",
        temp=directory("rmats/post/{experiment}/{reference}/temp/"),
        post=directory("rmats/post/{experiment}/{reference}/")
    params:
        read_length=lambda wildcards: get_experiment_readlen(wildcards.experiment),
        read_paired=lambda wildcards: "paired" if is_experiment_readpaired(wildcards.experiment) else "single",
        temp="rmats/post/{experiment}/{reference}/temp",
        rep_paired=lambda wildcards: "--paired-stats" if experiments.loc[wildcards.experiment,"pairRep"] else "",
        libtype=lambda wildcards: get_rmats_exp_libtype(wildcards.experiment),
    threads: 6
    resources:
        mem=lambda wildcards, input: (str((input.size//6000000000)+4) + "G"),
        rmem=lambda wildcards, input: (str((input.size//12000000000)+4) + "G"),
    conda:
        "../envs/rmats.yaml"
    shell:
        """
        for i in {input.control_bam} ; do readlink -f "$i" ; done |
        awk '{{printf ("%s,", $0) >> "{output.control}" }}' && 
        for i in {input.treat_bam} ; do readlink -f "$i" ; done |
        awk '{{printf ("%s,", $0) >> "{output.treat}" }}' &&
        for i in {input.prep} ;
          do for f in "$i"/temp/* ; 
             do cp "$f" {params.temp}/$(awk -v f="$f" 'BEGIN{{ match(f,/rmats\/prep\/([^\/]*)\/([^\/]*)\//,a) ; print a[1]"_"a[2] }}')_$(awk -v f="$f" 'BEGIN{{ match(f,/([^\/]*)$/,a) ; print a[1] }}') ;
          done ;
        done   

#        do python $(conda info | awk 'match($0,/active env location : ([^ ]*)/, a) {{print a[1]"/rMATS/cp_with_prefix.py" }}') \
#         $(awk -v i="$i" 'BEGIN{{ match(i,/rmats\/prep\/([^\/]*)\/([^\/]*)\//,a) ; print a[1]"_"a[2] }}') \
#         $(pwd)/{params.temp} \
#         $(ls -d "$i"/*) \
#        ; done &&

        python $(conda info | awk 'match($0,/active env location : ([^ ]*)/, a) {{print a[1]"/rMATS/rmats.py" }}') \
          --b1 {output.control} \
          --b2 {output.treat} \
          --gtf {input.gtf} \
          -t {params.read_paired} {params.rep_paired} \
          --readLength {params.read_length} \
          --variable-read-length \
          --libType {params.libtype} \
          --nthread {threads} \
          --od {output.post} \
          --tmp {params.temp} \
          --task post        
        """
