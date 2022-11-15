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
             num[$4]=(num[$4]=="")?a[1]:((a[1]>=num[$4])?a[1]:num[$4])
          }} 
          FNR<NR&&$8=="gene"&&$0~"gene_name"&&match($0,/gene_name "([^"]*)".*gene_biotype "([^"]*)"; .*/,a) {{
             num[$4]=(num[$4]==""?1:num[$4]);print $4,a[1],a[2],num[$4]
          }} 
          FNR<NR&&$8=="gene"&&$0!~"gene_name"&&match($0,/gene_biotype "([^"]*)".*/,a) {{
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

rule gtf_features:
    input:
        bed="{prefix}.gtf.bed",
    output:
        gene_tab="{prefix}.gtf.gene_info.tab",
        features="{prefix}.gtf.features.bed",
        feature_list="{prefix}.gtf.feature_list.tab",
        long_intron="{prefix}.gtf.long_intron.bed",
        long_exon="{prefix}.gtf.long_exon.bed",
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
        "logs/awk/{prefix}/gtf_features.log",
    shell:
        """        
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR&&$8=="exon"&&match($0,/exon_number "([^"]*)".*/,a) {{
             num[$4]=(num[$4]=="")?a[1]:((a[1]>=num[$4])?a[1]:num[$4])
          }}
          FNR<NR&&$8=="gene"&&$0~"gene_name"&&match($0,/gene_name "([^"]*)".*gene_biotype "([^"]*)"; .*/,a) {{
             num[$4]=(num[$4]==""?1:num[$4]);print $4,a[1],a[2],num[$4]
          }}
          FNR<NR&&$8=="gene"&&$0!~"gene_name"&&match($0,/gene_biotype "([^"]*)".*/,a) {{
            num[$4]=(num[$4]==""?1:num[$4]);print $4,$4,a[1],num[$4]
          }}
        ' {input.bed} {input.bed} > {output.gene_tab}

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            exon[1]=$4
          }}
          FNR<NR && match($0,/transcript_id "([^"]*)".*/,a) {{  
            if ($0~"transcript_name") {{
              match($0,/transcript_name "([^"]*)".*/,b)
            }}
            else {{
              b[1]=a[1]
            }} ; 
            if ($0 ~ "transcript_support_level") {{
              match($0, /transcript_support_level "([^"]).*"/, tsl) 
            }}
            else {{ 
              tsl[1]=0
            }} ;
            tsl[1]=((tsl[1]-0)>=1)?tsl[1]:0 ;  
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
            if (tag[1]=="basic") {{
              print $1, $2, $3, $4, $5, $6, $8, a[1], b[1], tsl[1], tag[1]
            }}
          }}
        ' {output.gene_tab} {input.bed} |

        sort -k7,7 -k1,1 -k8,8 -k2,2n -k3,3n |

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
        }}' - |

        sort -k7,7 -k1,1 -k4,4 -k2,2n -k3,3n - |
        uniq - |

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $8=="gene" {{
            if ($0~"gene_name") {{
              match($0,/gene_name "([^"]*)".*/,a)
            }}
            else {{
              a[1]=$4
            }} ;
            name[$4]=a[1] ;
            print $1, $2, $3, $4, $3-$2, $6, "gene", $4, name[$4], 0, "basic"   
          }}
          FNR < NR && $7!="transcript" {{
            $5=$3-$2 ; 
            if ($7=="exon") {{ 
              $7="transcript" ; print ;
              $7="exon" ; $8="" ; $9="" ; print ;  
            }} 
            else {{
              $8="" ; $9="" ; print ;
            }}
          }}
        ' {input.bed} - |

        sort -k11,11 -k7,7 -k1,1 -k4,4 -k2,2n -k3,3n - |

        awk -F'\\t' -v OFS='\\t' '
          name[FNR]=$1":"$2"-"$3":"$7 && n!=name[FNR] {{
            print s ;
            t=$10 ; n=name[FNR] ; s=$0
          }}
          name[FNR]=$1":"$2"-"$3":"$7 && n==name[FNR] {{
            t=($10<=t)?$10:t ; $10=t ; s=$0
          }}
          END {{
            print s
          }}
        ' - - |

        sort -k11,11 -k7,7 -k1,1 -k4,4 -k2,2n -k3,3n - |
        uniq - |
        awk -F'\\t' -v OFS='\\t' -f {params.feature_fwd} {output.gene_tab} - |
        sort -k11,11 -k7,7r -k1,1r -k4,4r -k3,3nr -k2,2nr - |
        awk -F'\\t' -v OFS='\\t' -f {params.feature_rev} - |
        sort -k11,11 -k7,7 -k1,1 -k4,4 -k2,2n -k3,3n - > {output.features} &&

#        awk -F'\\t' -v OFS='\\t' '
#          FNR==NR {{ 
#            v[$8]=(v[$8]>$12)?v[$8]:$12 
#          }}
#          FNR < NR {{
#            if (($7!="gene" || $7!="transcript") && (v[$8]>1)) {{
#              print ; $9=($9" var "$12) ; $8=($8"var"$12) ; $7=$7"_var";
#              print
#            }}
#            else {{ 
#             print
#            }}
#          }}' {output.features} {output.features} |
#
#        sort -k7,7 -k1,1 -k4,4 -k2,2n -k3,3n -o {output.features} - &&      

        awk -F'\\t' -v OFS='\\t' '
          $7=="long_intron"{{print}}
        ' {output.features} > {output.long_intron} &&

        awk -F'\\t' -v OFS='\\t' '
          $7=="exon"{{print}}
        ' {output.features} |

        bedtools intersect -s -f 1 -wa -a stdin -b {output.long_intron} |

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{ $7="skip_ex" ; print }}
          FNR< NR {{ print }}
        ' - {output.features} | 

        sort -o {output.features} -k1,1 -k2,2n - &&

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
        
        sort -k7,7 -k4,4 -k2,2n -k3,3n - |
        uniq - |
        sort -o {output.features} -k1,1 -k2,2n -  &&
        cut -f7 {output.features} |
        sort - |
        uniq - > {output.feature_list}
        """ 
