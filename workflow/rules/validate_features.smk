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
        
#        valid_features="resources/annotations/{reference}/{lineage}.gtf.{tag}_validated.bed"
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
            print $1, $2, $3, $4, $3-$2, $6, $8, $14/rpk[$4] >> "{output.alternative}" ;
          }} else if (alt==0) {{
            print $1, $2, $3, $4, $3-$2, $6, $8, $14/rpk[$4] >> "{output.principal}" ;
          }} ;
            alt=(acum[$4] >= {params.alt_cut})?1:0 ;
            print $0, $14/rpk[$4], acum[$4];
          }}' {output.expressed} {output.expressed}  > {output.rpk_ratio} &&
          
          awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
                 


# Define genebody range using confident (MANE entries or correct biotype and top tsl transcripts) transcript ranges
#        awk -F'\\t' -v OFS='\\t' '
#          FNR==NR && $10=={{ 
#            
#        cat {input.sj} |

        """ 

