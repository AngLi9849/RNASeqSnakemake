rule gtf_exon_intron_annotation:
    input:
        gene_type="{prefix}.gtf.gene_name_biotype.tab",
        bed="{prefix}.gtf.bed",
    output:
        exon="{prefix}.gtf.Exon.bed",
        exon_saf="{prefix}.gtf.Exon.saf",
        exon_var="{prefix}.gtf.ExonVariants.bed",
        exon_var_saf="{prefix}.gtf.ExonVariants.saf",
        intron="{prefix}.gtf.Intron.bed",
        intron_saf="{prefix}.gtf.Intron.saf",
    threads: 1
    resources:
        mem="16G",
        rmem="12G",
    log:
        "logs/awk/gtf_splice_sites.log",
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$8=="exon"&&match($0,/transcript_support_level "([^"]).*".*/,a){{if(a[1]==1){{print $1,$2,$3,$4,$5,$6}}}}' {input.bed} - |
        sort -k1,1 -k4,4 -k2,2n -k3,3n - |
        uniq - |
        awk -F'\\t' -v OFS='\\t' '
          {{
            if (gene!=$4) {{ 
              gene=$4 ; a=$2 ; b=$3 ; n=1 ; v=1 ; print $0, "exon", n, v
            }} 
            else if (gene==$4) {{
              gene==$4 ; a=$2 ;
              if (a>(b+2)) {{ 
                n+=1 ; b=$3 ; v=1 ; print $0, "exon", n, v
              }}
              else if (a<=(b+2)&&a>b) {{
                next 
              }}
              else if (a<=b) {{ 
                n+=0 ; b=b ; v+=1 ; print $0, "exon", n, v
              }} ;
            }}
          }}' - |
        sort -k1,1r -k4,4r -k3,3nr -k2,2nr - |
        uniq - |
        awk -F'\\t' -v OFS='\\t' '
          $6=="-"{{
            if (gene!=$4) {{
              gene=$4 ; a=$3 ; b=$2 ;  n=1 ; v=1 ; k=$8 ; print $1, $2, $3, $4, $5, $6, $7, n, v
            }}
            else if (gene==$4) {{
            gene==$4 ; a=$3 ;
              if (a<(b-2)) {{
                n+=1 ; b=$2 ; v=1 ; k=$8 ; print $1,$2,$3,$4,$5,$6, $7, n, v
              }}
              else if (a>=(b-2)) {{
                n+=0 ; b=b  ; k=k  ; 
                if (k>$8) {{
                  next
                }}
                else if (k==$8) {{
                  v+=1 ; print $1,$2,$3,$4,$5,$6, $7, n, v
                }}
              }} 
            }}
          }}
          $6=="+" {{if (gene!=$4) {{
            gene=$4 ; c=$3 ; d=$2 ; w=$8 ; print 
            }}
            else if (gene==$4) {{
            gene==$4 ; c=$3
              if (c<(d-2)) {{
                w=$8 ; d=$2 ; print
              }}
              else if (c>=(d-2)) {{
                w=w ;
                if (w>$8) {{
                  next
                }}
                else if (w==$8) {{
                  w=$8 ; print
                }}
              }}
            }}
          }}' - |
          sort -k1,1 -k2,2n - |
          awk -F'\\t' -v OFS='\\t' 'FNR==NR{{name[$1]=$2; type[$1]=$3 ; num[$1]=$4>1?"Multiexonic":"Monoexonic"}} FNR<NR{{printf "%s\\t%s\\t%s\\t%s\\t%s %s %s\\n", $0,num[$4],type[$4],name[$4],name[$4],"Exon",$8}}' {input.gene_type} - > {output.exon} &&
          sort -k1,1 -k4,4 -k2,2n -k3,3n {output.exon} |
          awk -F'\\t' -v OFS='\\t' '{{
            if (gene!=$4) {{
              gene=$4 ; a=$2 ; b=$3 
             }}
           else if (gene==$4) {{
             gene==$4 ; a=$2 ; 
              if (a>b) {{
                n=($6=="+")?($8-1):$8 ; printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s %s %s\\n", $1, b+1, a-1, $4, $5, $6, "intron", n, $9, $10, $11, $12, $12, "Intron", n ; b=$3 ; 
              }}
              else if (a<=b) {{
                b=$3 ; next
              }}
            }}
          }}' - |
        sort -k1,1 -k2,2n > {output.intron} &&
        awk -F'\\t' -v OFS='\\t' 'FNR==NR{{ n[$13]=($9>n[$13])?$9:n[$13] }} FNR<NR&&n[$13]>1{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s %s%s%s%s\\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$13,$12,"Ex",$8,"var",$9}}' {output.exon} {output.exon} > {output.exon_var} &&
        awk -F'\\t' -v OFS='\\t' '{{print $13, $1, $2, $3, $6}}' {output.exon} |
        sed '1 i\\GeneID\\tChr\\tStart\\tEnd\\tStrand' - > {output.exon_saf} &&
        awk -F'\\t' -v OFS='\\t' '{{print $13, $1, $2, $3, $6}}' {output.exon_var} |
        sed '1 i\\GeneID\\tChr\\tStart\\tEnd\\tStrand' - > {output.exon_var_saf} &&
        awk -F'\\t' -v OFS='\\t' '{{print $13, $1, $2, $3, $6}}' {output.intron} |
        sed '1 i\\GeneID\\tChr\\tStart\\tEnd\\tStrand' - > {output.intron_saf}        
        """

rule featurecounts_counts:
    input:
        bam="results/star/{sample}-{unit}/{prefix}.sortedByCoord.out.bam",
        bai="results/star/{sample}-{unit}/{prefix}.sortedByCoord.out.bam.bai",
        saf="resources/sample/sample_genome.{feature}.saf",
    output:
        "results/star/{sample}-{unit}/{prefix}_{feature}Usage.featurecounts.tab",
        "results/feature_counts/{sample}-{unit}/{prefix}_{feature}Usage.counts.tab",
    threads: 6
    resources:
        mem=lambda wildcards, input: (str((input.size//3000000000)+4) + "G"),
        rmem=lambda wildcards, input: (str((input.size//6000000000)+4) + "G"),
    log:
        "logs/feature_counts/{sample}-{unit}-{prefix}.{feature}Usage.log"
    conda:
        "../envs/subread.yaml",
    params:
        strand=get_sample_strandedness,
        paired=lambda wildcards:("" if not is_paired_end(wildcards.sample) else "-p")
    shell:
        """
        featureCounts -s {params.strand} {params.paired} -M -O -T {threads} -F SAF --verbose -a {input.saf} -o {output[0]} {input.bam} &&
        awk -F'\\t' -v OFS='\\t' 'FNR>2{{print $1,$7}}' {output[0]} |
        sort -k1,1 - |
        sed '1 i\\gene\\t{wildcards.sample}' > {output[1]}
        """


rule dexseq_exon_intron_usage:
    input:
        counts="results/counts/{prefix}_{feature}Usage_counts.tsv",
        sample_table="config/samples.tsv",
        bed="resources/sample/sample_genome.{feature}.bed",
    output:
        dir=directory("results/differential_usage/{feature}_Usage/{prefix}")
    params:
        plot_script="workflow/scripts/R/differential_plots.R",
        biotypes=config["biotypes"],
        exp=config["seq_experiment_name"],
        goi=config["GOI"],
        control=config["control_condition"],
        paired=config["paired_analysis"],
        dir="results/differential_usage/{feature}_Usage/{prefix}",
        ma_number=config["differential_plots"]["ma_gene_name_numbers"],
        volc_number=config["differential_plots"]["volcano_gene_name_numbers"],
        up_col=config["differential_plots"]["up_colour"],
        down_col=config["differential_plots"]["down_colour"],
        p_threshold=config["differential_plots"]["p_value_threshold"],
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/dexseq.yaml"
    log:
        "logs/dexseq/{prefix}.diff_{feature}_usage.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/R/dexseq_usage.R"

