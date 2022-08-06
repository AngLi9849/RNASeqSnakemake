rule get_ensembl_genome:
    output:
        "resources/genomes/ensembl_{species}.{build}.{release}_genome.fasta",
    log:
        "logs/{species}_{build}_{release}_genome.log",
    params:
        datatype="dna",
        species=lambda wildcards: wildcards.species,
        build=lambda wildcards: wildcards.build,
        release=lambda wildcards: wildcards.release,
    resources:
        mem="6G",
        rmem="4G",
    wrapper:
        "0.77.0/bio/reference/ensembl-sequence"

rule get_ensembl_annotation:
    output:
        "resources/annotations/ensembl_{species}.{build}.{release}_genome.gtf",
    params:
        fmt="gtf",
        species=lambda wildcards: wildcards.species,
        build=lambda wildcards: wildcards.build,
        release=lambda wildcards: wildcards.release,
        flavor="",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/{species}_{build}_{release}_annotation.log",
    wrapper:
        "0.77.0/bio/reference/ensembl-annotation"

rule get_refseq_genome:
    output:
        "resources/genomes/refseq_{species}.{build}.{release}_genome.fasta",
    log:
        "logs/{species}_{build}_{release}_genome.log",
    params:
        datatype="dna",
        species=lambda wildcards: wildcards.species,
        build=lambda wildcards: wildcards.build,
        release=lambda wildcards: wildcards.release,
    resources:
        mem="6G",
        rmem="4G",
    conda:
        "../envs/ncbi.yaml"
    script:
        "../scripts/py/refseq_fasta.py"

rule get_refseq_annotation:
    output:
        "resources/annotations/refseq_{species}.{build}.{release}_genome.gtf",
    params:
        fmt="gtf",
        species=lambda wildcards: wildcards.species,
        build=lambda wildcards: wildcards.build,
        release=lambda wildcards: wildcards.release,
        flavor="",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/{species}_{build}_{release}_annotation.log",
    wrapper:
        "../scripts/py/refseq_gtf.py"

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
    log:
        "logs/get_local_{species}_annotation.log",
    shell:
        "cat {input} > {output} 2> {log}"

rule prefix_spikein_genome:
    input:
        fasta = "resources/genomes/{species}_genome.fasta",
    output:
        "resources/genomes/spikein_{species}_genome.fasta",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/ref/prefix_spikein_{species}_genome.log",
    shell:
        "sed 's/>/>spikein_/' {input} > {output} 2> {log}"

rule prefix_spikein_annotation:
    input:
        gtf="resources/annotations/{species}_genome.gtf",
    output:
        gtf="resources/annotations/spikein_{species}_genome.gtf",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/set_spikein_{species}_annotation.log",
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$2!=""{{ $1 = "spikein_"$1 ; print }}' {input.gtf} > {output.gtf} 
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

rule combine_genome_and_annotation:
    input:
        sample_genome=lambda wildcards: "resources/genomes/{species}_genome.fasta",
        spikein_genome="resources/genomes/spikein_{spikein}_genome.fasta",
        sample_gtf=lambda wildcards: "resources/annotations/{species}_genome.gtf",
        spikein_gtf="resources/annotations/spikein_{spikein}_genome.gtf",
    output:
       fasta="resources/genomes/combined_{species}_and_{spikein}_genome.fasta",
       gtf="resources/annotations/combined_{species}_and_{spikein}_genome.gtf",
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/ref/combine_{species}_and_{spikein}_genome.log",
    shell:
        """
        cat {input.sample_genome} {input.spikein_genome} > {output.fasta} &&
        cat {input.sample_gtf} {input.spikein_gtf} > {output.gtf} 
        """

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
    wrapper:
        "0.77.0/bio/samtools/faidx"

rule bwa_index:
    input:
        lambda wildcards: "resources/genomes/{prefix}_genome.fasta"
    output:
        multiext("resources/bwa/{prefix}_genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/bwa_{prefix}_genome_index.log",
    resources:
        mem_mb=369000,
    wrapper:
        "0.77.0/bio/bwa/index"

rule gtf_features:
    input:
        bed="resources/annotations/{prefix}.gtf.bed",
        chr="resources/genomes/{prefix}.fasta.chrom.sizes",
    output:
        gene_tab="resources/annotations/{prefix}.gtf.{tag}_gene_info.tab",
        inconfident="resources/annotations/{prefix}.gtf.{tag}_inconfident.bed",
        confident="resources/annotations/{prefix}.gtf.{tag}_confident.bed",
        transcripts="resources/annotations/{prefix}.gtf.{tag}_transcripts.bed",
        features="resources/annotations/{prefix}.gtf.{tag}_annotated.bed",
        feature_list="resources/annotations/{prefix}.gtf.{tag}_feature_list.tab",
        long_intron="resources/annotations/{prefix}.gtf.{tag}_long_intron.bed",
        long_exon="resources/annotations/{prefix}.gtf.{tag}_long_exon.bed",
    params:
        min_ret_cov=config["features"]["minimum_retained_intron_coverage"],
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
        "logs/awk/{prefix}/gtf_{tag}_features.log",
    shell:
        """        
# Find gene names, biotypes and exon numbers
        awk -F'\\t' -v OFS='\\t' ' 
          FNR==NR {{
            if ($8=="exon") {{
              match($0,/exon_number "([^"]*)".*/,a);
              num[$4]=(num[$4]=="")?a[1]:((a[1]>=num[$4])?a[1]:num[$4])
            }} ;
            if (match($0,/tag "{wildcards.tag}"/)) {{
              t[$4]+=1
            }}
            else if ("{wildcards.tag}"=="basic"){{
              t[$4]+=1
            }}
            else {{
              t[$4]+=0
            }}
          }}
          FNR<NR&&$8=="gene"&&$0~"gene_name"&&match($0,/gene_name "([^"]*)".*gene_biotype "([^"]*)"; .*/,a) && t[$4]>=1{{
             num[$4]=(num[$4]==""?1:num[$4]);print $4,a[1],a[2],num[$4]
          }}
          FNR<NR&&$8=="gene"&&$0!~"gene_name"&&match($0,/gene_biotype "([^"]*)".*/,a) && t[$4]>=1 {{
            num[$4]=(num[$4]==""?1:num[$4]);print $4,$4,a[1],num[$4]
          }}
        ' {input.bed} {input.bed} > {output.gene_tab} &&

# Sort transcript features to a tab-demilited bed file with transcript id, name, transcript support levels, parental gene id and transcript biotype 
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            exon[$1]=$4 ; presence[$1]=1 ; biotype[$1] = $3 ;
          }}
          FNR<NR && match($0,/transcript_id "([^"]*)".*/,a) && presence[$4]=1 && $8 != "gene" {{  
            if ($0~"transcript_name") {{
              match($0,/transcript_name "([^"]*)".*/,b)
            }}
            else {{
              b[1]=a[1]
            }} ;
             if ($0 ~ "transcript_biotype") {{
              match($0, /transcript_biotype "([^"]).*"/, tbt)
            }}
            else {{
              tbt[1]=biotype[$4]
            }} ;
            if ($0 ~ "transcript_support_level") {{
              match($0, /transcript_support_level "([^"]).*"/, tsl) 
            }}
            else {{ 
              tsl[1]=6
            }} ;
            tsl[1]=((tsl[1]-0)>=1)?tsl[1]:6 ; 
            if ($8=="five_prime_utr") {{
              $8="5UTR"
            }}
            else if ($8=="three_prime_utr") {{
              $8="3UTR"
            }} ;
            if ($0 ~ "tag") {{
              match($0,/tag "([^"]*)"/,tag)
            }} 
            else {{ 
               tag[1]="basic" 
            }}; 
            tag[1]=(tag[1]=="")?"basic":tag[1] ;
            if ( ( tag[1] == "{wildcards.tag}" ) || ( $0 ~ "tag "{wildcards.tag}"" ) ) {{
            print $1, $2, $3, $4, $5, $6, $8, a[1], b[1], tsl[1], $4, tbt[1]
            }}
          }}
        ' {output.gene_tab} {input.bed} |

# Sort through each transcript in to generate introns
        sort -k7,7 -k4,4 -k8,8 -k2,2n -k3,3n |

        awk -F'\\t' -v OFS='\\t' '{{
          if (id != $8) {{
            n=$0 ; t=$8 ;
            $0 = m ; $2 = a ; $3 = b ; print ;
            m=n ; $0=n ; a=$2 ; b=$3 ; id=t
          }}
          else if ($7 != "exon") {{
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
              $7 = "intron" ; $2 = b ; $3 = c ; print ; 
              a=d ; b=e ; m=n ;
            }}
          }}
        }} END {{
          $0=m ; $2=a ; $3=b ; print ;
        }}' - > {output.transcripts} &&

# Filter out transcript features with biotypes mismatching the parental gene biotype to inconfident entries
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            biotype[$1] = $3 ;
          }}
          FNR < NR {{
            if ($12==biotype[$4]) {{
              print
            }} else {{
              print >> "{output.inconfident}"
            }}
          }}' {output.gene_tab} {output.transcripts} |

# Find highest tsl of the gene and filter out transcript features with less than top tsl of the gene
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $7=="transcript" {{
            tsl[$4]= (tsl[$4]==0) ? $10 : ( \ 
              (tsl[$4] <= $10) ? tsl[$4] : $10 \ 
            ) ;
          }} 
          FNR < NR {{ 
            if ($10==tsl[$4]) {{
              print
            }} else {{
              print >> "{output.inconfident}"
            }}
          }}' {output.gene_tab} - > {output.confident}
        """ 
