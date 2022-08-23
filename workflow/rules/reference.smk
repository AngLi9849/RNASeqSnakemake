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

rule mane_genome:
    input:
        "resources/genomes/ensembl_{prefix}.fasta",
    output:
        "resources/genomes/MANE_{prefix}.fasta",
    resources:
        mem="6G",
        rmem="4G",
    shell:
        """
        cat {input} > {output}
        """

rule get_ensembl_annotation:
    output:
        "resources/annotations/ensembl_{species}.{build}.{release}/genome.gtf",
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

rule incorporate_MANE_annotations:
    input:
        ensembl="resources/annotations/ensembl_{prefix}.gtf",
    output:
        gtf="resources/annotations/MANE_{prefix}.gtf",
    params:
        mane_link=config["MANE_annotation"],
        mane = "resources/annotations/MANE.gtf",
    threads: 1
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/awk/{prefix}_MANE_gtf.log",
    shell:
        """
        curl {params.mane_link} -o {params.mane} &&
        awk -F'\\t' -v OFS='\\t' '{{
          $2="MANE" ; 
          print ;
        }}' {params.mane} |
        cat - {input.ensembl} |
        sort -k1,1 -k4,4n > {output.gtf}
        """

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
        "resources/annotations/local_{species}_genome.gtf",
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
        gtf="resources/annotations/spikein_{species}/genome.gtf",
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
        sample_gtf=lambda wildcards: "resources/annotations/{species}/genome.gtf",
        spikein_gtf="resources/annotations/spikein_{spikein}/genome.gtf",
    output:
       fasta="resources/genomes/combined_{species}_and_{spikein}_genome.fasta",
       gtf="resources/annotations/combined_{species}_and_{spikein}/genome.gtf",
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

rule gtf_transcripts:
    input:
        bed="resources/annotations/{prefix}/genome.gtf.bed",
        chr="resources/genomes/{prefix}_genome.fasta.chrom.sizes",
    output:
        gene_tab="resources/annotations/{prefix}/genome.gtf.{tag}_gene_info.tab",
        trs_tab="resources/annotations/{prefix}/genome.gtf.{tag}_trs_info.tab",
        transcripts="resources/annotations/{prefix}/genome.gtf.{tag}_transcripts.bed",
        inconfident="resources/annotations/{prefix}/genome.gtf.{tag}_inconfident.bed",
        confident="resources/annotations/{prefix}/genome.gtf.{tag}_confident.bed",
        biotyped="resources/annotations/{prefix}/genome.gtf.{tag}_biotyped.bed",
        indexed="resources/annotations/{prefix}/genome.gtf.{tag}_transcripts.indexed.bed",
    params:
        intron_min=config["features"]["minimum_intron_length"],
        tsl_tol=config["ensembl_annotations"]["transcript_support_level_tolerance"],
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
              match($0, /transcript_biotype "([^"]*).*"/, tbt)
            }}
            else {{
              tbt[1]=biotype[$4]
            }} ;
            if ($0 ~ "tag") {{
              match($0,/tag "([^"]*)"/,tag)
            }}
            else {{
               tag[1]="basic"
            }};
            tag[1]=(tag[1]=="")?"basic":tag[1] ;
            if ($0 ~ "transcript_support_level") {{
              match($0, /transcript_support_level "([^"]).*"/, tsl) 
            }}
            else {{ 
              tsl[1]=6
            }} ;
            tsl[1]=(tag[1]=="CCDS")?1:((tsl[1]-0)>=1)?tsl[1]:6 ; 
            if ($8=="five_prime_utr") {{
              $8="5UTR"
            }}
            else if ($8=="three_prime_utr") {{
              $8="3UTR"
            }} ;
            if ($7 == "MANE") {{
              match($4,/([^\.]*)\..*/,g) ;
              $4 = g[1] ;
              tag[1]="basic" ;
              tsl[1]=-1 ;
              tbt[1]=biotype[$4] ;
            }} ; 
            if ( ( tag[1] == "{wildcards.tag}" ) || ( $0 ~ "tag \\"{wildcards.tag}\\"" ) ) {{
            print $1, $2, $3, $4, $5, $6, $8, a[1], b[1], tsl[1], $4, tbt[1]
            }}
          }}
        ' {output.gene_tab} {input.bed} |

# Sort through each transcript in to generate introns
        sort -k7,7 -k4,4 -k8,8 -k2,2n -k3,3n |

        awk -F'\\t' -v OFS='\\t' '{{
          if (id != $8) {{
            n=$0 ; t=$8 ;
            $0 = m ; $2 = a ; $3 = b ; 
            if (FNR>1) {{ print }} ;
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
          }}' {output.gene_tab} {output.transcripts} > {output.biotyped} &&

# Find highest tsl of the gene and filter out transcript features with less tsl than tolerance of the gene
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $7=="transcript" && $10 > 0{{
            tsl[$4] = (tsl[$4]==0 || tsl[$4] > $10)?$10:tsl[$4] ;
          }}
          FNR < NR {{
            if ($10<=(tsl[$4]+{params.tsl_tol})) {{
              $12="confident" ;
              print $0, tsl[$4]
            }} else {{
              print >> "{output.inconfident}"
            }}
          }}' {output.biotyped} {output.biotyped} > {output.confident} && 

# Transcript names
        awk -F'\\t' -v OFS='\\\t' '{{ 
          print $8, $9
        }}' {output.transcripts} |
        sort | uniq > {output.trs_tab} &&

# Index features in each transcript by their order
        awk -F'\\t' -v OFS='\\t' ' {{
          $4=$8 ; print
        }}' {output.transcripts} |

        cut -f1-11 | 
        sort -k7,7 -k4,4 -k2,2n -k3,3n |
        awk -F'\\t' -v OFS='\\t' -f {params.feature_fwd} {output.trs_tab} - |
        sort -k7,7r -k4,4r -k3,3nr -k2,2nr |
        awk -F'\\t' -v OFS='\\t' -f {params.feature_rev} - |
        sort -k1,1 -k2,2n > {output.indexed}
        """

rule annotated_features:
    input:
        bed="resources/annotations/{prefix}.gtf.bed",
        chr="resources/genomes/{prefix}.fasta.chrom.sizes",
        transcripts="resources/annotations/{prefix}.gtf.{tag}_transcripts.bed",
        gene_tab="resources/annotations/{prefix}.gtf.{tag}_gene_info.tab",
        inconfident="resources/annotations/{prefix}.gtf.{tag}_inconfident.bed",
        confident="resources/annotations/{prefix}.gtf.{tag}_confident.bed",
    output:
        features="resources/annotations/{prefix}.gtf.{tag}_annotated.bed",
        feature_list="resources/annotations/{prefix}.gtf.{tag}_feature_list.tab",
        long_intron="resources/annotations/{prefix}.gtf.{tag}_long_intron.bed",
        long_exon="resources/annotations/{prefix}.gtf.{tag}_long_exon.bed",
    params:
        min_ret_cov=config["features"]["minimum_retained_intron_coverage"],
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
# Assign each confident transcript features a index
        
        sort -k7,7 -k4,4 -k2,2n -k3,3n - |
        uniq - |



        awk -F'\\t' -v OFS='\\t' '
#          FNR==NR && $8=="gene" {{
#            if ($0~"gene_name") {{
#              match($0,/gene_name "([^"]*)".*/,a)
#            }}
#            else {{
#              a[1]=$4
#            }} ;
#            name[$4]=a[1] ;
#            print $1, $2, $3, $4, $3-$2, $6, "gene", $4, name[$4], 0, $4
#          }}
          FNR < NR && $7!="transcript" {{
            $5=$3-$2 ; 
            if ($7=="exon") {{ 
              $7="trscrpt" ; print ;
              $7="exon" ; $8="" ; $9="" ; print ;  
            }} 
            else {{
              $8="" ; $9="" ; print ;
            }}
          }}
        ' {input.bed} - |

        sort -k7,7 -k4,4 -k2,2n -k3,3n - |

        awk -F'\\t' -v OFS='\\t' '
          {{
            name[FNR]=$1":"$2"-"$3":"$7;
            if (n != name[FNR]) {{
              print s ;
              t=$10 ; n = name[FNR] ; s=$0
            }} 
            else {{
              t = ($10<=t)? $10 : t ; $10=t ; s=$0
            }} ;
          }}
          END {{
            print s
          }}
        ' - |

        sort -k7,7 -k4,4 -k2,2n -k3,3n - |
        uniq - |
        awk -F'\\t' -v OFS='\\t' -f {params.feature_fwd} {input.gene_tab} - |
        sort -k7,7r -k4,4r -k3,3nr -k2,2nr - |
        awk -F'\\t' -v OFS='\\t' -f {params.feature_rev} - |
        awk -F'\\t' -v OFS='\\t' '
          $7=="long_intron"{{print >> "{output.long_intron}"}}
          $7=="long_exon" {{ print >> "{output.long_exon}"}}
          $7!~/^long_/ {{ print }}
        ' - |
        sort -k7,7 -k4,4 -k2,2n -k3,3n - > {output.features} &&

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{ 
            v[$8]=(v[$8]>$12)?v[$8]:$12 
          }}
          FNR < NR {{
            if (($7!="gene" || $7!="trscrpt") && (v[$8]>1)) {{
              print ; $9=($9" var "$12) ; $8=($8"var"$12) ; $7=$7"_var";
              start[$8]=(start[$8]==0 || start[$8]>$2)?$2:start[$8] ;
              end[$8]=(end[$8]==0 || end[$8]<$3)?$3:end[$8] ;
            }}
            else {{ 
             print
            }}
          }}' {output.features} {output.features} |

        sort -k7,7 -k4,4 -k2,2n -k3,3n -o {output.features} - &&      

        awk -F'\\t' -v OFS='\\t' '
          $7=="exon" || $7=="intron" {{print}}
        ' {output.features} |

        bedtools intersect -s -f 1 -wa -wb -a stdin -b {output.long_intron} |

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $7=="exon" {{ $7="skip_ex" ; print >> "{output.features}" }}
          FNR==NR && $7=="intron" {{
            
        ' && 



        awk -F'\\t' -v OFS='\\t' '
          $7=="long_exon"{{print}}
        ' {output.features} > {output.long_exon} &&

        awk -F'\\t' -v OFS='\\t' '
          $7=="intron"{{print}}
        ' {output.features} |

        bedtools intersect -s -f {params.min_ret_cov} -wa -a stdin -b {output.long_exon} |

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{ $7="ret_int" ; print }}
          FNR< NR {{ print }}
        ' - {output.features} |
        
        sort  -k7,7 -k4,4 -k2,2n -k3,3n - |
        uniq - |
        
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            print $1, 0, $2, $1, $2, "+", "chr", $1, $1, 0, $1, 1, 1, 1 ;
            print $1, 0, $2, $1, $2, "-", "chr", $1, $1, 0, $1, 1, 1, 1 ;
            print $1, 0, $2, "{wildcards.prefix}", $2, "+", "genome", "{wildcards.prefix}", "{wildcards.prefix}", 0, "{wildcards.prefix}", 1, 1, 1 ;
            print $1, 0, $2, "{wildcards.prefix}", $2, "+", "genome", "{wildcards.prefix}", "{wildcards.prefix}", 0, "{wildcards.prefix}", 1, 1, 1 ;
          }}
          FNR < NR {{ print }}' {input.chr} - |

        sort -o {output.features} -k1,1 -k2,2n -  &&
        cut -f7 {output.features} |
        sort - |
        uniq - > {output.feature_list}
        """ 
