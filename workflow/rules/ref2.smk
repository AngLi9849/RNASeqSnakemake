rule get_ensembl_genome:
    output:
        "resources/genomes/ensembl_{species}_genome.fasta",
    log:
        "logs/{species}_genome.log",
    params:
        datatype="dna",
        species=lambda wildcards: wildcards.species,
        build=lambda wildcards: references.loc[wildcards.species,"ensembl_build"],
        release=lambda wildcards: references.loc[wildcards.species,"ensembl_release"],
    resources:
        mem="6G",
        rmem="4G",
    cache: True
    wrapper:
        "0.77.0/bio/reference/ensembl-sequence"

rule get_ensembl_annotation:
    output:
        "resources/annotations/ensembl_{species}_genome.gtf",
    params:
        fmt="gtf",
        species=lambda wildcards: wildcards.species,
        build=lambda wildcards: references.loc[wildcards.species,"ensembl_build"],
        release=lambda wildcards: references.loc[wildcards.species,"ensembl_release"],
        flavor="",
    cache: True
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/{species}_annotation.log",
    wrapper:
        "0.77.0/bio/reference/ensembl-annotation"

rule get_local_genome:
    input:
        lambda wildcards: str(references.loc[wildcards.species,"genome_dir"]),
    output:
        "resources/genomes/local_{species}_genome.fasta",
    log:
        "logs/local_{species}_genome.log",
    resources:
        mem="6G",
        rmem="4G",
    cache: True
    shell:
        "cat {input} > {output} 2> {log}"

rule get_local_annotation:
    input:
        lambda wildcards: str(references.loc[wildcards.species,"annotation_dir"]),
    output:
        "resources/annotations/local_{species}_annotations.gtf",
    resources:
        mem="6G",
        rmem="4G",
    cache: True
    log:
        "logs/get_local_{species}_annotation.log",
    shell:
        "cat {input} > {output} 2> {log}"

rule gtf2bed:
    input:
        "{prefix}.gtf",
    output:
        "{prefix}.gtf.bed",
    log:
        "logs/bedops/{prefix}_gtf2bed.log",
    params:
        bedops_fix="workflow/scripts/awk/bedops_fix.awk"
    resources:
        mem="6G",
        rmem="4G",
    conda:
        "../envs/bedops.yaml",  
    shell:
        "awk -f {params.bedops_fix} {input} | " 
        "gtf2bed - > {output}"

rule bed2saf:
    input:
        "{prefix}.bed",
    output:
        "{prefix}.bed.saf",
    log:
        "logs/awk/{prefix}_bed2saf.log",
    resources:
        mem="6G",
        rmem="4G",
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '{{print $4,$1,$2+1,$3,$6}}' {input} |
        sed '1 i\\GeneID\\tChr\\tStart\\tEnd\\tStrand' > {output}
        """

rule gene_name_biotype:
    input:
        "{prefix}.gtf.bed",
    output:
        "{prefix}.gtf.gene_name_biotype.tab"
    log:
        "logs/awk/{prefix}_gene_name_biotype.log",
    resources:
        mem="6G",
        rmem="4G",
    shell:
        """
	awk -F'\\t' -v OFS='\\t' '
          FNR==NR&&$8=="exon"&&match($0,/exon_number "([^"]*)".*/,a) {{
             num[$4]=(num[$4]>=1)?((a[1]>=num[$4])?a[1]:num[$4]):1
          }} 
          FNR<NR&&$8=="gene"&&$0~"gene_name"&&match($0,/gene_name "(.+)"; gene_source.+gene_biotype "(.+)"; .*/,a) {{
             num[$4]=(num[$4]==""?1:num[$4]);print $4,a[1],a[2],num[$4]
          }} 
          FNR<NR&&$8=="gene"&&$0!~"gene_name"&&match($0,/gene_biotype "(.+)"; .*/,a) {{
            num[$4]=(num[$4]==""?1:num[$4]);print $4,$4,a[1],num[$4]
          }}
        ' {input} {input} > {output}
	"""

rule chrom_sizes:
    input:
        "{prefix}.fasta",
    output:
        "{prefix}.fasta.chrom.sizes",
    params:
        chrom_sizes_awk="workflow/scripts/awk/chrom_sizes.awk"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/ref/{prefix}_chrom_sizes.log",
    shell:
        "awk -f {params.chrom_sizes_awk} {input} > {output}"

rule combine_genome:
    input:
        sample_genome=lambda wildcards: get_genome(wildcards.species),
        spikein_genome="resources/genomes/spikein_{spikein_species}_genome.fasta",
    output:
       "resources/genomes/combined_{species}_and_{spikein_species}_genome.fasta",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/ref/combine_{species}_and_{spikein_species}_genome.log",
    shell:
        "cat {input.sample_genome} {input.spikein_genome} > {output} 2> {log}"

rule combine_annotation:
    input:
        sample_gtf=lambda wildcards: get_annotation(wildcards.species),
        spikein_gtf="resources/annotations/spikein_{spikein_species}_genome.gtf",
    output:
       "resources/annotations/combined_{species}_and_{spikein_species}_genome.gtf",
    params:
        col2check_awk="workflow/scripts/awk/col2check.awk"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/ref/combine_{species}_and_{spikein_species}_genome_gtf.log",
    shell:
        "cat {input.sample_gtf} > {output} ; "
        "awk -F'\\t' -v OFS='\\t' -f {params.col2check_awk} {input.spikein_gtf} | tail -n +1 - >> {output}"

rule genome_faidx:
    input:
        "{prefix}.fasta",
    output:
        "{prefix}.fasta.fai",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/{prefix}-faidx.log",
    cache: True
    wrapper:
        "0.77.0/bio/samtools/faidx"

rule bwa_index:
    input:
        select_genome,
    output:
        multiext("{experiment}/bwa/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/bwa_{experiment}_genome_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "0.77.0/bio/bwa/index"

rule gtf_feature_annotation:
    input:
        gene_type="{prefix}.gtf.gene_name_biotype.tab",
        bed="{prefix}.gtf.bed",
    output:
        genes="{prefix}.gtf.genes.bed",
        genes_fwd="{prefix}.gtf.genes_fwd.bed",
        features="{prefix}.gtf.features.bed",
#        long_intron="{prefix}.gtf.long_intron.bed",
#        long_exon="{prefix}.gtf.long_exon.bed",
    params:
        min_tsl=config["features"]["minimum_transcript_support_level"],
        min_ret_tsl=config["features"]["minimum_retained_intron_support_level"],
        intron_min=config["features"]["minimum_intron_length"],
        feature_fwd="workflow/scripts/awk/feature_index_fwd.awk",
        feature_rev="workflow/scripts/awk/feature_index_rev.awk",
    threads: 1
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/bedtools.yaml",
    log:
        "logs/awk/{prefix}/gtf_features.log",
    shell:
        """        
        awk -F'\\t' -v OFS='\\t' -v min_tsl={params.min_tsl} -v min_ret_tsl={params.min_ret_tsl} '
          match($0,/transcript_id "([^"]*)".*transcript_biotype "([^"]*)".*transcript_support_level "([^"]).*".*/,a) {{
            if ((a[2] != "retained_intron" && a[3] <= min_tsl) || (a[2]=="retained_intron" && a[3] <= min_ret_tsl)) {{
              if ($0~"transcript_name") {{
                match($0,/transcript_name "([^"]*)".*/,b)
              }}
              else {{
                b[1]=a[1]
              }} ; 
              if ($8=="five_prime_utr") {{
                f="5_UTR"
              }}
              else if ($8=="three_prime_utr") {{
                f="3_UTR"
              }}
              else {{
                f=$8 
              }} ;
              print $1, $2, $3, $4, $5, $6, a[2], a[1], f, a[1], b[1]
            }}
          }}
        ' {input.bed} |

        sort -k9,9 -k1,1 -k8,8 -k2,2n -k3,3n |
       
        awk -F'\\t' -v OFS='\\t' '{{
          if (id != $8) {{ 
            n=$0 ; t=$8 ;
            $0 = m ; $2 = a ; $3 = b ; print ; 
            m=n ; $0=n ; a=$2 ; b=$3 ; id=t
          }}
          else if ($9 != "exon") {{
            b=$3
          }}
          else {{
            c = $2 ; 
            if ( (c-b) < {params.intron_min} ) {{
              b = $3
            }}
            else {{ 
              d=$2 ; e=$3 ; n=$0 ;
              $0=m ; $2=a ; $3=b ; print ;
              $9 = "intron" ; $2 = b ; $3 = c ; print ; 
              a=d ; b=e ; m=n ;
            }}
          }}
        }} END {{
          $0=m ; $2=a ; $3=b ; print ;
        }}' - |

        sort -k9,9 -k1,1 -k4,4 -k2,2n -k3,3n - |
        uniq - |

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR&&$8=="gene"&&$0~"gene_name"&&match($0,/gene_name "(.+)"; gene_source.+gene_biotype "(.+)"; .*/,a) {{
             name[$4]=a[1] ; type[$4]=a[2] ;
             print $1, $2, $3, $4, $3-$2, $6, type[$4], name[$4], "gene", $4, name[$4]   
          }}
          FNR==NR&&$8=="gene"&&$0!~"gene_name"&&match($0,/gene_biotype "(.+)"; .*/,a) {{
             name[$4]=$4 ; type[$4]=a[1] ;
             print $1, $2, $3, $4, $3-$2, $6, type[$4], name[$4], "gene", $4, name[$4]
          }}
          FNR < NR && $9!="transcript" {{
            $7=type[$4] ; $8=name[$4] ; $5=$3-$2 ; 
            if ($9=="exon") {{ 
              $9="transcript" ; print $0 ;
              $9="exon" ; $10="" ; $11="" ; print ;  
            }} 
            else {{
              $10="" ; $11="" ; print ;
            }}
          }}
        ' {input.bed} - > {output.genes} &&

        sort -k9,9 -k1,1 -k4,4 -k2,2n -k3,3n {output.genes} |
        uniq - |
        awk -F'\\t' -v OFS='\\t' -f {params.feature_fwd} - > {output.genes_fwd} &&
        sort -k9,9r -k1,1r -k4,4r -k3,3nr -k2,2nr {output.genes_fwd} |
        awk -F'\\t' -v OFS='\\t' -f {params.feature_rev} - |
        sort -k1,1 -k2,2n - > {output.features}
        """  
