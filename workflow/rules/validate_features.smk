rule star_detected_splice_junctions:
    input:
        sj=lambda wildcards: expand(
                "star/{sample.sample_name}/{sample.unit_name}/{sample.reference}/SJ.out.tab",
                sample=lineage.loc[wildcards.species].loc[wildcards.lineage].itertuples(),
        ),
    output:
        sj="resources/annotations/{species}.{lineage}.star.splice_junctions.bed",
    threads: 1
    resources:
        mem="16G",
        rmem="12G",
    conda:
        "../envs/bedtools.yaml",
    log:
        "logs/awk/{species}_{lineage}_splice_junctions.log",
    shell:
        """
        cat {input.sj} |

        awk -F'\\t' -v OFS='\\t' '
          BEGIN {{
            split("nan|GT:AG|CT:AC|GC:AG|CT:GC|AT:AC|GT:AT",m,"|") ;
          }}
          FNR==NR && $6>=1 {{
            s=$4 ;
            $2=$2-1 ; 
            $4=($4==1)?"+":(($4==2)?"-":"") ;
            name=$1":"$2"-"$3":"$4 ;
            motif=m[($5+1)] ;
            print $1, $2, $3, name, motif, $4, $7, $8, $9
          }}' - |

        sort -k4,4 |

        awk -F'\\t' -v OFS='\\t' '
          name != $4 {{
            save=$0 ; $0=load ; $7=n ; $8=m ; $9=o ;
            if (FNR>1) {{
              print
            }} ;
            $0=save ; name=$4 ; n=$7 ; m=$8 ; o=$9 ; load=$0
          }}
          name==$4 {{
            n+=$7 ; m+=$8 ; o=(o>$9)?o:$9 ;
          }}
          END {{
            $0=load ; $7=n ; $8=m ; $9=o ; print
          }}' - > {output.sj}
        """

rule transcript_bed2fasta:
    input:
        bed = "resources/annotations/{reference}/genome.gtf.bed12",
        fasta = "resources/genomes/{reference}_genome.fasta",
    output:
        fasta = "resources/annotations/{reference}/transcriptome.fasta",
    threads: 1
    resources:
        mem="6G",
        rmem="4G",
    conda:
        "../envs/bedtools.yaml"    
    shell:
        """
        bedtools getfasta -split -nameOnly -fi {input.fasta} -bed {input.bed} > {output.fasta}
        """

rule salmon_lineage_transcriptome_quant:
    input:
        bam="star/{sample}/{unit}/{reference}/Aligned.toTranscriptome.out.bam",
        fasta="resources/annotations/{reference}/transcriptome.fasta",
    output:
        quant="salmon/{sample}/{unit}/{reference}/quant.sf"
    params: 
        libtype="A",
        outdir=lambda wildcards, output: dirname(output.quant),
    threads: 2
    resources:
        mem="8G",
        rmem="6G",
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon quant -t {input.fasta} -l {params.libtype} -a {input.bam} -o {params.outdir} --seqBias --gcBias
        """

rule validate_transcripts:
    input:
        bed="resources/annotations/{reference}/genome.gtf.bed",
        transcripts="resources/annotations/{reference}/genome.gtf.{tag}_transcripts.bed",
        salmon_quant = lambda wildcards: expand(
            "salmon/{sample.sample_name}/{sample.unit_name}/{sample.reference}/quant.sf",
            sample=lineage.loc[get_reference_species(wildcards.reference)].loc[wildcards.lineage][lineage.trs_val.tolist()].itertuples(),
        ),
        gene_tab="resources/annotations/{reference}/genome.gtf.{tag}_gene_info.tab",
        sj=lambda wildcards: "resources/annotations/{species}.{{lineage}}.star.splice_junctions.bed".format(
            species=get_reference_species(wildcards.reference),
        ),
        confident="resources/annotations/{reference}/genome.gtf.{tag}_confident.bed",
    output:
        expressed="resources/annotations/{reference}/{lineage}.gtf.{tag}_expressed.bed",
        rpk_ratio="resources/annotations/{reference}/{lineage}.gtf.{tag}_rpk-ratio.bed",
        principal="resources/annotations/{reference}/{lineage}.gtf.{tag}_principal.bed",
        alternative="resources/annotations/{reference}/{lineage}.gtf.{tag}_alternative.bed",
        main="resources/annotations/{reference}/{lineage}.gtf.{tag}_main.bed",
#        features="resources/annotations/{reference}/{lineage}.gtf.{tag}_validated.bed",
    params:
        min_overhang=config["lineage_feature_validation"]["splicing"]["minimum_overhang"],
        min_int_uniq=config["lineage_feature_validation"]["splicing"]["minimum_unique_splice_reads"],
        min_splice = config["lineage_feature_validation"]["splicing"]["minimum_multimap_splice_reads"],
        min_ret_cov=config["features"]["minimum_retained_intron_coverage"],
        alt_cut=config["lineage_feature_validation"]["genes"]["principal_transcripts_threshold"], 
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
        "logs/awk/{reference}/{lineage}/validate_{tag}_features.log",
    shell:
        """
# Take all salmon quant files and calculate accumulative reads, effective rpk and effective tpm of each transcript
# Print only transcripts with  more than 1 reads into 'expressed' bed
        cat {input.salmon_quant} |
  
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $1!="Name" {{
            rpk[$1]+=($5/$3) ;
            nreads[$1]+=$5 ;
            rpksum+=($5/$3) ; 
          }}
          FNR < NR && $7=="transcript" {{
            if (nreads[$8]>0) {{
              tpm[$8] = rpk[$8]*1000000/rpksum ;
              print $0, nreads[$8], rpk[$8], tpm[$8]  ;
            }} 
          }}' - {input.transcripts} |
          
          sort -k4,4 -k14,14nr > {output.expressed} &&
         
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            rpk[$4] += $14   
          }}
          FNR < NR {{
          if (gene!=$4) {{
            alt=0;
            gene=$4;
          }} ;       
            acum[$4]+=$14/rpk[$4] ;
          if (alt==1) {{
            print $1, $2, $3, $4, $3-$2, $6, $7, $8, $9, $14/rpk[$4], $4, "alternative"  >> "{output.alternative}" ;
          }} else if (alt==0) {{
            print $1, $2, $3, $4, $3-$2, $6, $7, $8, $9, $14/rpk[$4], $4, "principal" >> "{output.principal}" ;
          }} ;
            alt=(acum[$4] >= {params.alt_cut})?1:0 ;
            print $0, $14/rpk[$4], acum[$4];
          }}' {output.expressed} {output.expressed}  > {output.rpk_ratio} &&

# Determine principal transcription/genebody start and end site based on largest span of principal transcripts        
# If a gene does not have detected transcript, use annotated confident transcripts 
# If neither expressed nor confident transcripts exists, use genebody annotation
 
        cat {output.principal} {input.confident} |
        cut -f1-12 |
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $12=="principal" {{
            princ[$4]+=1 ;
            principal[$8] += 1 ;
            p_start[$4]=(p_start[$4]=="" || p_start[$4]>$2)?$2:p_start[$4] ;
            p_end[$4]=(p_end[$4]=="" || p_end[$4]<$3)?$3:p_end[$4] ;
          }}
          FNR==NR && $12 !="principal" {{
            conf[$4] += 1 ;
            confident[$8] += 1 ;
            c_start[$4]=(c_start[$4]=="" || c_start[$4]>$2)?$2:c_start[$4] ;
            c_end[$4]=(c_end[$4]=="" || c_end[$4]<$3)?$3:c_end[$4] ;
          }}      
          FNR < NR && $8=="gene" {{
            if ($0~"gene_name") {{
              match($0,/gene_name "([^"]*)".*/,a)
            }}
            else {{
              a[1]=$4
            }} ;
            name[$4]=a[1] ;
            if (princ[$4] >= 1 ) {{
              print $1, p_start[$4], p_end[$4], $4, p_end[$4] - p_start[$4], $6, "gene", $4, name[$4], 0, $4 ;
            }} else if (conf[$4] >= 1) {{
              print $1, c_start[$4], c_end[$4], $4, c_end[$4] - c_start[$4], $6, "gene", $4, name[$4], 0, $4 ;
            }} else {{ 
              print $1, $2, $3, $4, $3-$2, $6, "gene", $4, name[$4], 0, $4 ;
            }}
          }}' - {input.bed} |
#          FNR < NR && $7!="transcript" {{
#          if ( (principal[$12] >= 1) || (princ[$4] < 1 && confident[$12] >= 1 ) {{ 
#            $10=(principal[$12] >= 1)?1:$10 ;
#            $5=$3-$2 ;
#            if ($7=="exon") {{
#              $7="trscrpt" ; print ;
#              $7="exon" ; $8="" ; $9="" ; print ;
#            }}
#            else {{
#              $8="" ; $9="" ; print ;
#            }}
#          }}
#        ' - {input.transcripts} |

        

        sort -k7,7 -k4,4 -k2,2n -k3,3n |

        awk -F'\\t' -v OFS='\\t' '
          {{
            name=$1":"$2"-"$3":"$7;
            if (n != name) {{
              print s ;
              t=$10 ; n = name ; s=$0
            }}
            else {{
              t = ($10<=t)? $10 : t ; $10=t ; s=$0
            }} ;
          }}
          END {{
            print s
          }}
        ' - > {output.main} &&

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            u[$4]==i

        sort -k7,7 -k4,4 -k2,2n -k3,3n - |
        uniq - |
        awk -F'\\t' -v OFS='\\t' -f {params.feature_fwd} {input.gene_tab} - |
        sort -k7,7r -k4,4r -k3,3nr -k2,2nr - |
        awk -F'\\t' -v OFS='\\t' -f {params.feature_rev} - |
        sort -k7,7 -k4,4 -k2,2n -k3,3n - > {output.features}

        

# Define genebody range using confident (MANE entries or correct biotype and top tsl transcripts) transcript ranges
#        awk -F'\\t' -v OFS='\\t' '
#          FNR==NR && $10=={{ 
#            
#        cat {input.sj} |

        """ 

