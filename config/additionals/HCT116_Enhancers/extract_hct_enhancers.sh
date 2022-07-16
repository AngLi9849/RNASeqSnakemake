#!/bin/bash
zcat ENCFF513PJK.bed.gz > HCT116_chromHMM.bed &&
awk -F'\t' -v OFS='\t' ' \
  FNR==NR { type[$1] = $3 } \
  FNR<NR && type[$4] == "protein_coding" {print} \
' ./ensembl_homo_sapiens_genome.gtf.basic_gene_info.tab ./ensembl_homo_sapiens_genome.custom.annotated_basic.PROMPT-TSS.bed > HCT116_protein_coding_PROMP-TSS.bed
