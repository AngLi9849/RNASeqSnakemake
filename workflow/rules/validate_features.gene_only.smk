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

rule salmon_index:
    input:
        trs = "resources/annotations/{reference}/transcriptome.{tag}_{type}.fasta",
        decoy="resources/salmon/{reference}/transcriptome.{tag}_{type}.decoy.txt"
    output:
        idx=directory("resources/salmon/{reference}/transcriptome.{tag}_{type}.idx"),
    params:
        idx="resources/salmon/{reference}/transcriptome.{tag}_{type}.idx",
    threads: 8
    resources:
        mem="16G",
        rmem="8G",
    conda:
        "../envs/salmon.yaml"    
    shell:
        """
        salmon index -t {input.trs} --decoys {input.decoy} -i {params.idx} -k 31
        """

rule salmon_lineage_transcriptome_quant:
    input:
        idx="resources/salmon/{reference}/transcriptome.{tag}_{type}.idx",
        fq=get_salmon_fq
    output:
        quant="salmon/{sample}/{unit}/{reference}/{tag}_{type}/quant.sf"
    params: 
        libtype="A",
        outdir=lambda wildcards, output: dirname(output.quant),
        fqs=lambda wildcards, input: get_salmon_input(wildcards, input)
    threads: 8
    resources:
        mem="8G",
        rmem="4G",
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon quant -i {input.idx} -l {params.libtype} {params.fqs} -o {params.outdir} --seqBias --gcBias
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
            for (i=1 ; i<= NF ; i++) {{
              if ($i=="riExonStart_0base") {{ rx_start=i }} ;
              if ($i=="riExonEnd") {{ rx_end=i }} ;
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
                list[length(list)]=id;
                seen[id]=1 ;
              }} ;
              include[id] += $incl ;
              skipping[id] += $skip ;
              inc_rpk[id]+= ($incl/($inclen+1)) ;
              skip_rpk[id]+= ($skip/($skiplen+1)) ;
            }} ;
          }}
          END {{
            for ( i in list ) {{
              if ( i != 0 ) {{
                match(list[i], /^([^:]*):([^-]*)-([^:]*):([^-]*)-([^:]*):(.*)/,x) ;
                print x[1], x[4],x[5],list[i],x[5] - x[4],x[6],include[list[i]], skipping[list[i]], inc_rpk[list[i]],skip_rpk[list[i]],inc_rpk[list[i]]/(skip_rpk[list[i]] + inc_rpk[list[i]]) >> "{output.ret_int}" ; 
                print x[1], x[2], x[3], list[i],x[3] - x[2],x[6],include[list[i]], skipping[list[i]], inc_rpk[list[i]],skip_rpk[list[i]],inc_rpk[list[i]]/(skip_rpk[list[i]] + inc_rpk[list[i]]) >> "{output.ret_ex}" ;   
              }} ;
            }} ;
          }}' - ;
    
        for i in {input.prep} ; do cat $i/SE.MATS.JCEC.txt ; done |
        awk -F'\\t' -v OFS='\\t' '
          BEGIN {{
            list[0]==""
          }}
          FNR==1 {{
            for (i=1 ; i<= NF ; i++) {{
              if ($i=="exonStart_0base") {{ sx_start=i }} ;
              if ($i=="exonEnd") {{ sx_end=i }} ;
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
                list[length(list)]=id;
                seen[id]=1 ;
              }} ;
              include[id] += $incl ;
              skipping[id] += $skip ;
              inc_rpk[id]+= ($incl/($inclen)) ;
              skip_rpk[id]+= ($skip/($skiplen)) ;
            }} ;
          }}
          END {{
            for ( i in list ) {{
              if ( i != 0 ) {{
                match(list[i], /^([^:]*):([^-]*)-([^:]*):([^-]*)-([^:]*):(.*)/,x) ;
                print x[1], x[4],x[5],list[i],x[5] - x[4],x[6],include[list[i]], skipping[list[i]], inc_rpk[list[i]],skip_rpk[list[i]],inc_rpk[list[i]]/(skip_rpk[list[i]] + inc_rpk[list[i]]) >> "{output.skip_in}" ;
                print x[1], x[2], x[3], list[i],x[3] - x[2],x[6],include[list[i]], skipping[list[i]], inc_rpk[list[i]],skip_rpk[list[i]],inc_rpk[list[i]]/(skip_rpk[list[i]] + inc_rpk[list[i]]) >> "{output.skip_ex}" ;
              }} ;
            }} ;
          }}' -
        """  

rule validate_main_transcripts:
    input:
        bed="resources/annotations/{reference}/genome.gtf.bed",
        transcripts="resources/annotations/{reference}/genome.gtf.{tag}_transcripts.bed",
        trs_idx="resources/annotations/{reference}/genome.gtf.{tag}_transcripts.indexed.bed",
        trs_sj="resources/annotations/{reference}/genome.gtf.{tag}_trs_sj.bed",
        salmon_all = lambda wildcards: expand(
            "salmon/{sample.sample_name}/{sample.unit_name}/{sample.reference}/{{tag}}_transcripts/quant.sf",
            sample=lineage[lineage.trs_val.tolist()].loc[get_reference_species(wildcards.reference)].loc[wildcards.lineage].itertuples(),
        ),
        salmon_confident = lambda wildcards: expand(
            "salmon/{sample.sample_name}/{sample.unit_name}/{sample.reference}/{{tag}}_confident/quant.sf",
            sample=lineage[lineage.trs_val.tolist()].loc[get_reference_species(wildcards.reference)].loc[wildcards.lineage].itertuples(),
        ),
        gene_tab="resources/annotations/{reference}/genome.gtf.{tag}_gene_info.tab",
        sj=lambda wildcards: "resources/annotations/{species}.{{lineage}}.star.splice_junctions.bed".format(
            species=get_reference_species(wildcards.reference),
        ),
        confident="resources/annotations/{reference}/genome.gtf.{tag}_confident.bed",
    output:
        expressed="resources/annotations/{reference}/{lineage}.gtf.{tag}_expressed.bed",
        form="resources/annotations/{reference}/{lineage}.gtf.{tag}_form.tab",
        features="resources/annotations/{reference}/{lineage}.gtf.{tag}_validated.bed",
    params:
        alt_cut=config["lineage_feature_validation"]["genes"]["principal_transcripts_threshold"], 
        min_reads=config["lineage_feature_validation"]["genes"]["minimum_transcript_reads"],
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
# Take all salmon quant files for confident transcripts and calculate accumulative reads, effective rpk and effective tpm of each transcript
# Print rpk ratio of a gene of each confident transcript form (defined by their first and last splice sites)
 
        cat {input.salmon_confident} |

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $1!="Name" {{
            rpk[$1]+=($5/$3) ;
            nreads[$1]+=$5 ;
            rpksum+=($5/$3) ;
          }}
          FNR < NR {{
            if (nreads[$8]>0) {{
              nreads[$4]+=nreads[$8] ;
              form[$4][$15]= 1 ;
              form_rpk[$4][$15]+=rpk[$8] ;
              form_reads[$4][$15]+=nreads[$8] ;
            }} ;
          }}
          END {{
            for (i in form) {{
              for (j in form[i]) {{
                print i, j, form_rpk[i][j], form_reads[i][j], form_reads[i][j]/nreads[i];
              }}
            }}
          }}
        ' - {input.transcripts} {input.transcripts} |

        sort -k1,1 -k4,4nr > {output.form} |
 

# Take all salmon quant files and calculate accumulative reads, effective rpk and effective tpm of each transcript
# Print only transcripts with  more than 1 reads into 'expressed' bed
        cat {input.salmon_all} |
  
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $1!="Name" {{
            rpk[$1]+=($5/$3) ;
            nreads[$1]+=$5 ;
            rpksum+=($5/$3) ; 
          }}
          FNR < NR {{
            if (nreads[$8]> {params.min_reads}) {{
              tpm[$8] = rpk[$8]*1000000/rpksum ; $9=$15 ; $10=$14 ;
              $12="expressed" ; $13=nreads[$8] ; $14=rpk[$8] ; $15=tpm[$8] ;
              print ; 
            }} ; 
          }}' - {input.transcripts} |
          
        sort -k4,4 -k10,10nr -k14,14nr > {output.expressed} &&

        cat {output.form} {output.expressed} | 
         
        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && NF<8 {{
            read_acum[$1]+=$5 ;
            principal[$1][$2]=(alt[$1]==1)?0:1 ;
            alt[$1]=(read_acum[$1] >= {params.alt_cut})?1:0 ;  
          }} 
          FNR == NR && NF > 4 && $12=="expressed" {{
            if ( (selected[$4]-0)==0) {{
              if (principal[$4][$9]==1) {{
                model[$8]=1 ;
                selected[$4]=1 ;
              }}
            }}
          }}
          FNR<NR && NF<10 {{
            if (selected[$1]==0) {{
              print $5, $6, $7, $1, $7-$6, $8, "gene", $1, $2, 1, $1,  1, 1, 1 ;
            }} ; 
            name[$1]=$2 ;
          }}
          FNR<NR && NF>12 {{
            if (model[$4]==1) {{
              if ($7=="transcript") {{ 
                print $1, $2, $3, $11, $5, $6, "gene", $11, name[$11], 0, $11, 1, 1, 1 ;
              }} else {{
                $5=$3-$2 ; $4=$11 ; $8=$11":"$7$12 ; $9=name[$11]" "$7" "$12 ; $10=0 ; print ;
              }}
            }}
          }}
        ' - {input.gene_tab} {input.trs_idx}  > {output.features}

        
        """ 
