rule fasta2bed:
    input:
        "{prefix}.fasta",
    output:
        "{prefix}.fasta.bed",
    log:
        "logs/bedops/{prefix}_fasta2bed.log",
    resources:
        mem="6G",
        rmem="4G",
    threads: 1
    shell:
        """
        awk '
          $0 ~ ">" {{
            if (NR > 1) {{
              print c;
            }} ;  
            c = 0 ;
            printf substr( $1, 2, length($1) )"\t"0"\t" ; 
            }} 
          $0 !~ ">" {{
            c+=length($0);
          }} 
          END {{ 
            print c; 
          }}
        ' {intput} > {output}
        """
 
rule transcript_bed12_fasta_decoy:
    input:
        fasta = "resources/genomes/{reference}_genome.fasta",
        bed="resources/annotations/{reference}/genome.gtf.bed",
        transcript="resources/annotations/{reference}/genome.gtf.{tag}_{type}.bed",
    output:
        bed12="resources/annotations/{reference}/genome.gtf.{tag}_{type}.bed12",
        fasta = "resources/annotations/{reference}/transcriptome.{tag}_{type}.fasta",
        decoy="resources/salmon/{reference}/transcriptome.{tag}_{type}.decoy.txt",
    log:
        "logs/bed12/resources/annotations/{reference}/genome_{tag}_{type}_bed12.log",
    params:
        temp="resources/annotations/{reference}/genome.gtf.{tag}_{type}.temp.txt",
    resources:
        mem="6G",
        rmem="4G",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '
          $8=="exon" && match($0,/transcript_id "([^"]*)"/,t) {{
            print $1, $2, $3, t[1], $3-$2, $6, $8;
          }}' {input.bed} |
        sort -k7,7 -k1,1 -k4,4 -k2,2n > {params.temp} &&

        cat {params.temp} {input.transcript} | 

        awk -F'\\t' -v OFS='\\t' '
          FNR==NR && $7=="exon" {{
            a=(l[$4]=="")?$2:a ;
            l[$4] = (l[$4]=="")?$5:(l[$4]","$5) ;
            s[$4] = (s[$4]=="")?0:(s[$4]","($2-a)) ;
            n[$4] +=1 ;
          }}
          FNR==NR && $7=="transcript" {{
            seen[$8]==1
          }}
          FNR < NR {{ 
            if (seen[$4]==1) {{
              $7=$2 ; 
              print $0,$3 , 0, n[$4], l[$4], s[$4]
            }}
          }}' - {params.temp} > {output.bed12} &&

        bedtools getfasta -split -nameOnly -fi {input.fasta} -bed {output.bed12} |
        cat - {input.fasta} > {output.fasta} &&

        cat {input.fasta} |
        awk 'match($0,/^>([^ ]*)/,a) {{print a[1]}}' > {output.decoy}
        """


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
        awk -F'\\t' -v OFS='\\t' '{{print $8,$1,$2+1,$3,$6}}' {input} |
        sed '1 i\\GeneID\\tChr\\tStart\\tEnd\\tStrand' > {output}
        """

rule pptx_template:
    output:
        pptx="resources/templates/{width}cm_wide.{height}cm_tall.pptx",
    conda:
        "../envs/ms_office.yaml"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/pptx/{width}cm_w.{height}cm_h_pptx.log"
    script:
        "../scripts/py/pptx.py"

rule docx_template:
    output:
        docx="resources/templates/{font_colour}_{font}.docx",
    conda:
        "../envs/ms_office.yaml"
    resources:
        mem="6G",
        rmem="4G",
    log:
        "logs/pptx/{font_colour}_{font}_docx.log"
    script:
        "../scripts/py/docx.py"

