rule star_detected_splice_junctions:
    input:
        sj=lambda wildcards: expand(
                "star/{sample.sample_name}/{sample.unit_name}/{sample.reference}/SJ.out.tab",
                sample=lineage[lineage.splice.tolist()].loc[wildcards.species].loc[wildcards.lineage].itertuples(),
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

rule rmats_lineage_intron_retention_exon_skip:
    input:
        prep=lambda wildcards: expand(
            "rmats/prep/{sample.sample_name}/{sample.unit_name}/{sample.reference}/post",
            sample=lineage[lineage.trs_val.tolist()].loc[wildcards.species].loc[wildcards.lineage].itertuples(),
        ),
    output:
       ret_int="resources/annotations/{species}.{lineage}.rmats.ret_int.bed",
       skip_ex="resources/annotations/{species}.{lineage}.rmats.skip_ex.bed",
       skip_in="resources/annotations/{species}.{lineage}.rmats.skip_in.bed",
       ret_ex="resources/annotations/{species}.{lineage}.rmats.ret_ex.bed",
    threads: 2
    resources:
        mem="8G",
        rmem="6G",
    shell:
        """    
        for i in {input.prep} ; do cat $i/RI.MATS.JCEC.txt ; done |
        awk -F'\\t' -v OFS='\\t' '
          BEGIN {{
            list[0]==""
          }}
          FNR==1 {{ 
            for (i=1 ; i<= NR ; i++) {{
              if ($i=="riExonStart_0Base") {{ rx_start=i }} ;
              if ($i=="riExonStart_0Base") {{ rx_end=i }} ;
              if ($i=="upstreamEE") {{ ri_start=i }} ;
              if ($i=="downstreamES") {{ ri_end=i }} ;
              if ($i=="chr") {{ chr=i }} ;
              if ($i=="strand") {{ strand=i }} ;
              if ($i=="IJC_SAMPLE_1") {{ incl=i }} ;
              if ($i=="SJC_SAMPLE_1") {{ skip=i }} ;
              if ($i=="IncFormLen") {{ inclen=i }} ;
              if ($i=="SkipFormLen") {{ skiplen=i }} ;
            }}
          }}        
          FNR>1 && $1!="ID" {{
            if ($incl > 0 || $skip > 0 ) {{
              match($chr,/chr(.*)/,c) ;
              id=c[1]":"$rx_start"-"$rx_end":"$ri_start"-"$ri_end":"$strand ;
              if (seen[id]==0) {{ 
                list[length(ri_list)]=id;
                seen[id]=1 ;
              }} ;
              inc_rpk[id]+= ($incl/($inclen+1)) ;
              skip_rpk[id]+= ($skip/($skiplen+1)) ;
            }} ;
          }}
          END {{
            for ( i in list ) {{
              if ( i != 0 ) {{
                match(list[i], /^([^:]*):([^-]*)-([^:]*):([^-]*)-([^:]*):(.*)/,x) ;
                print x[1], x[4],x[5],list[i],x[5] - x[4],x[6],inc_rpk[list[i]],skip_rpk[list[i]],inc_rpk[list[i]]/(skip_rpk[list[i]] + inc_rpk[list[i]]) >> "{output.ret_int}" ; 
                print x[1], x[2], x[3], list[i],x[3] - x[2],x[6],inc_rpk[list[i]],skip_rpk[list[i]],inc_rpk[list[i]]/(skip_rpk[list[i]] + inc_rpk[list[i]]) >> "{output.ret_ex}" ;   
              }} ;
            }} ;
          }}' - ;
    
        for i in {input.prep} ; do cat $i/SE.MATS.JCEC.txt ; done |
        awk -F'\\t' -v OFS='\\t' '
          BEGIN {{
            list[0]==""
          }}
          FNR==1 {{
            for (i=1 ; i<= NR ; i++) {{
              if ($i=="ExonStart_0Base") {{ sx_start=i }} ;
              if ($i=="ExonStart_0Base") {{ sx_end=i }} ;
              if ($i=="upstreamEE") {{ si_start=i }} ;
              if ($i=="downstreamES") {{ si_end=i }} ;
              if ($i=="chr") {{ chr=i }} ;
              if ($i=="strand") {{ strand=i }} ;
              if ($i=="IJC_SAMPLE_1") {{ incl=i }} ;
              if ($i=="SJC_SAMPLE_1") {{ skip=i }} ;
              if ($i=="IncFormLen") {{ inclen=i }} ;
              if ($i=="SkipFormLen") {{ skiplen=i }} ;
            }}
          }}
          FNR>1 && $1!="ID" {{
            if ($incl > 0 || $skip > 0 ) {{
              match($chr,/chr(.*)/,c) ;
              id=c[1]":"$sx_start"-"$sx_end":"$si_start"-"$si_end":"$strand ;
              if (seen[id]==0) {{
                list[length(si_list)]=id;
                seen[id]=1 ;
              }} ;
              inc_rpk[id]+= ($incl/($inclen+1)) ;
              skip_rpk[id]+= ($skip/($skiplen+1)) ;
            }} ;
          }}
          END {{
            for ( i in list ) {{
              if ( i != 0 ) {{
                match(list[i], /^([^:]*):([^-]*)-([^:]*):([^-]*)-([^:]*):(.*)/,x) ;
                print x[1], x[4],x[5],list[i],x[5] - x[4],x[6],inc_rpk[list[i]],skip_rpk[list[i]],inc_rpk[list[i]]/(skip_rpk[list[i]] + inc_rpk[list[i]]) >> "{output.skip_in}" ;
                print x[1], x[2], x[3], list[i],x[3] - x[2],x[6],inc_rpk[list[i]],skip_rpk[list[i]],inc_rpk[list[i]]/(skip_rpk[list[i]] + inc_rpk[list[i]]) >> "{output.skip_ex}" ;
              }} ;
            }} ;
          }}' -
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
        long_intron="resources/annotations/{reference}/{lineage}.gtf.{tag}_long_intron.bed",
        sj_int="resources/annotations/{reference}/{lineage}.gtf.{tag}_sj_intron.bed",
        long_exon="resources/annotations/{reference}/{lineage}.gtf.{tag}_long_exon.bed",
        main="resources/annotations/{reference}/{lineage}.gtf.{tag}_main.bed",
        features="resources/annotations/{reference}/{lineage}.gtf.{tag}_validated.bed",
    params:
        min_overhang=config["lineage_feature_validation"]["splicing"]["minimum_overhang"],
        min_uniq=config["lineage_feature_validation"]["splicing"]["minimum_unique_splice_reads"],
        min_mult = config["lineage_feature_validation"]["splicing"]["minimum_multimap_splice_reads"],
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
        cat - {input.bed} |
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
          FNR == NR && $8=="gene" {{
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
          }}
          FNR < NR && $7!="transcript" {{
          if ( (principal[$8] >= 1) || (princ[$4] < 1 && confident[$8] >= 1 )) {{ 
            $10=(principal[$8] >= 1)?1:$10 ;
            $5=$3-$2 ;
            if ($7=="exon") {{
              $7="trscrpt" ; print ;
              $7="exon" ; $8="" ; $9="" ; print ;
            }}
            else {{
              $8="" ; $9="" ; print ;
            }}
          }}
        }}' - {input.transcripts} |

        cut -f1-11 |
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

        sort -k7,7 -k4,4 -k2,2n -k3,3n {output.main} |
        uniq - |
        awk -F'\\t' -v OFS='\\t' -f {params.feature_fwd} {input.gene_tab} - |
        sort -k7,7r -k4,4r -k3,3nr -k2,2nr - |
        awk -F'\\t' -v OFS='\\t' -f {params.feature_rev} - |
        awk -F'\\t' -v OFS='\\t' '
#          $7=="long_intron"{{print >> "{output.long_intron}"}}
          $7=="long_exon" {{ print >> "{output.long_exon}"}}
          $7!~/^long_/ {{ print }}
        ' - |
        sort -o {output.main} -k8,8 - &&

        awk -F'\\t' -v OFS='\\t' '
          BEGIN{{
            id[0]=""
          }}
          FNR==NR {{
            v[$8]=(v[$8]>$14)?v[$8]:$14
          }}
          FNR < NR {{
            if (($7!="gene" || $7!="trscrpt") && (v[$8]>1)) {{
              if (seen[$8]==0) {{ 
                id[length(id)]=$8 ; 
                seen[$8]=1  ;
                entry[$8]=$0 ;
              }} ; 
              $9=($9" var "$14) ; $8=($8"var"$14) ; $7=$7"_var"; print ;
              start[$8]=(start[$8]==0 || start[$8]>$2)?$2:start[$8] ;
              end[$8]=(end[$8]==0 || end[$8]<$3)?$3:end[$8] ;	      
            }}
            else {{
              print
            }}
          }}
          END {{
            for (i in id) {{
              if id[i]!="" {{
                $0=entry[id[i]] ; 
                $2=start[id[i]] ;
                $3=end[id[i]] ;
                $14=1 ;
                print $0 ;
              }}
            }}
          }}' {output.main} {output.main} |

        sort -o {output.main} -k1,1 -k2,2n - &&
     
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $7=="intron" {{
            id=$1":"$2"-"$3":"$6 ;
            seen[id]=1 ;             
          }}
          FNR<NR && $9>{params.min_overhang} && ($7>{params.min_uniq} || $8 > {params.min_mult}) {{
            if (seen[$4]==0) {{
              print $1, $2, $3, "intron", ".", $6
            }}
          }} ' {output.main} {input.sj} > {output.sj_int} &&

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $7=="intron" {{
            print
          }}' {output.main} |
 
        bedtools intersect -a {output.sj_int} -b - -wa -c -F 1 |
 
        awk -F'\\t' -v OFS='\\t' '
          $7 ==1 {{ print >> "{output.alt_int}" }} 
          $7 >=2 {{ print >> "{output.long_intron}" }}
        ' - &&

        awk -F'\\t' -v OFS='\\t' '
          {{print $1,$2,$3, $1":"$2"-"$3":"$6, ".", $6 }}
        ' {output.long_intron} |

        sort | 
        uniq |
        sort -o {output.long_intron} -k1,1 -k2,2n - &&      

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR {{
            n_ex[$1]=$4
          }}
          FNR < NR && ($7=="intron" || $7=="exon") {{
            if (n_ex[$4] >= 2) {{
              print
            }}
          }}' {output.main} |

        bedtools intersect -a {output.long_intron} -b - -wa -wb -F 1 |

        sort -k13,13 -k10,10 -k4,4 -k18,18n |

        awk -F'\\t' -v OFS='\\t' '
          $13=="intron" {{
            if ( id=="" || (id!=$4 && id!="")) {{
              if (id!=$4 && id!="") {{
                  part[id]=first[id] ; 
                  for (i=first[id]+1 ; i <= last[id] ; i++) {{
                    part[id]=part[id]"+"i
                  }} ;
                print chr[id], start[id], end[id], gene[id], end[id]-start[id], strand[id], "big_int", gene[id]":intron:"part[id], gname[id]" "intron" "part[id], 1, gene[id], first[id], last_rev[id], 1,
              }} ;
              id=$4 ; 
              chr[$4]=$1 ;
              gene[$4]=$10 ;
              strand[$4]=$6 ; 
              match($15,/^([^ ]*) /,n) ;
              gname[$4]=n[1] ;
            }} ;
            start[$4]=(start[$4]==0 || start[$4] >= $8)? $8 : start[$4] ; 
            end[$4]=(end[$4]==0 || end[$4] <= $8)? $8 : end[$4] ;
            first[$4]=(first[$4]==0 || first[$4] >= $18)?$18:first[$4] ;
            last[$4]=(last[$4]==0 || last[$4] <= $19)?$19:last[$4] ;
            last_rev[$4]=(last_rev[$4]==0 || last_rev[$4] >= $19)? $19 : last_rev[$4] ;
          }}
          $13=="exon" {{
            print $7,$8,$9,$10, $11, $12, "skip_ex", $14, $15, $16, $17, $18, $19 ,$20
          }}
          END {{
            part[id]=first[id] ;
            for (i=first[id]+1 ; i <= last[id] ; i++) {{
              part[id]=part[id]"+"i
            }} ;
            print chr[id], start[id], end[id], gene[id], end[id]-start[id], strand[id], "big_int", gene[id]":intron:"part[id], gname[id]" "intron" "part[id], 1, gene[id], first[id], last_rev[id], 1,
          }} 
        ' - >> {output.main} && 

        awk -F'\\t'

        """ 

