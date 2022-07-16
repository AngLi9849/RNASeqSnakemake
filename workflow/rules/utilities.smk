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

