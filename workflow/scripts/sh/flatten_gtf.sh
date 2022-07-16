awk -F'\\t' -v OFS='\\t' '$8=="exon" && match($0,/transcript_id "([^"]*)".*exon_number "([^"]*)".*transcript_support_level "([^"]).*".*/,a) {{ if (a[3]==1) { print $1, $2, $3, $4, $5, $6, "transcript", a[1],"exon",a[2] }} $8!="exon"&&$8!="CDS"&&match($0,/transcript_id "([^"]*)".*transcript_support_level "([^"]).*".*/,a) { if (a[2]==1) { print $1, $2, $3, $4, $5, $6, "transcript", a[1], $8, 1 }}' resources/annotations/ensembl_homo_sapiens_genome.gtf.bed |  sort -k1,1 -k8,8 -k9,9 -k2,2n -k3,3n

awk -F'\\t' -v OFS='\\t' '$8=="exon"&&match($0,/transcript_support_level "([^"]).*".*/,a){{if(a[1]==1){{print $1,$2,$3,$4,$5,$6}}}}' {input.bed} - |

sort -k1,1 -k4,4 -k2,2n -k3,3n - |

uniq - |

awk -F'\\t' -v OFS='\\t' '
  {{
    if (gene!=$4) {{
      gene=$4 ; a=$2 ; b=$3 ; n=1 ; v=1 ; print $0, "exon", n, v
    }}
    else if (gene==$4) {{
      gene==$4 ; a=$2 ;
      if (a>(b+2)) {{
        n+=1 ; b=$3 ; v=1 ; print $0, "exon", n, v
      }}
      else if (a<=(b+2)&&a>b) {{
        next
      }}
      else if (a<=b) {{
        n+=0 ; b=b ; v+=1 ; print $0, "exon", n, v
      }} ;
    }}
  }}' - |

sort -k1,1r -k4,4r -k3,3nr -k2,2nr - |

uniq - |
        
awk -F'\\t' -v OFS='\\t' '
  $6=="-"{{
    if (gene!=$4) {{
      gene=$4 ; a=$3 ; b=$2 ;  n=1 ; v=1 ; k=$8 ; print $1, $2, $3, $4, $5, $6, $7, n, $8, v
    }}
    else if (gene==$4) {{
      gene==$4 ; a=$3 ;
      if (a<(b-2)) {{
        n+=1 ; b=$2 ; v=1 ; k=$8 ; print $1,$2,$3,$4,$5,$6, $7, n, $8, v
      }}
      else if (a>=(b-2)) {{
        n+=0 ; b=b  ; k=k  ;
        if (k>$8) {{
          next
        }}
        else if (k==$8) {{
          v+=1 ; print $1,$2,$3,$4,$5,$6, $7, n, $8, v
        }}
              }}
            }}
          }}
          $6=="+" {{if (gene!=$4) {{
            gene=$4 ; c=$3 ; d=$2 ; w=$8 ; m=1 ; $10=$9 ; $8=w ; $9=m ; print
            }}
            else if (gene==$4) {{
            gene==$4 ; c=$3
              if (c<(d-2)) {{
               m+=1 ;  w=$8 ; d=$2 ; $10=$9 ; $8=w ; $9=m ; print
              }}
              else if (c>=(d-2)) {{
                m+=0 ; w=w ;
                if (w>$8) {{
                  next
                }}
                else if (w==$8) {{
                  w=$8 ; $10=$9 ; $8=w ; $9=m ; print
                }}
              }}
            }}
          }}' - |
          sort -k1,1 -k2,2n - |
          awk -F'\\t' -v OFS='\\t' 'FNR==NR{{name[$1]=$2; type[$1]=$3 ; num[$1]=$4>1?"Multiexonic":"Monoexonic"}} FNR<NR{{printf "%s\\t%s\\t%s\\t%s\\t%s %s %s\\t%s\\t%s %s%s%s%s\\n", $0,num[$4],type[$4],name[$4],name[$4],"exon",$8,"variant",name[$4],"Ex",$8,"var",$10}}' {input.gene_type} - > {output.exon} &&
          sort -k1,1 -k4,4 -k2,2n -k3,3n {output.exon} |
          awk -F'\\t' -v OFS='\\t' '{{
            if (gene!=$4) {{
              gene=$4 ; a=$2 ; b=$3
             }}
           else if (gene==$4) {{
             gene==$4 ; a=$2 ;
              if (a>b) {{
                n=($6=="+")?($8-1):$8 ; m=($6=="+")?$9:($9-1) ; printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s %s %s\\t%s\\t%s %s%s%s\\n", $1, b+1, a-1, $4, $5, $6, "intron", n, m, $10, $11, $12, $13, $13, "intron", n, "variant", $13, "In", n, "var1" ; b=$3 ;
              }}
              else if (a<=b) {{
                b=$3 ; next
              }}
            }}
          }}' - |
        sort -k1,1 -k2,2n > {output.intron}
        awk  -F'\\t' -v OFS='\\t' 'FNR==NR{{name[$1]=$2; type[$1]=$3 ; num[$1]=$4>1?"Multiexonic":"Monoexonic"}} FNR<NR&&$8=="gene"{{print $1, $2, $3, $4, $5, $6, "gene", 1, 1, 1, num[$4],type[$4],name[$4],name[$4],"body",name[$4]}}' {input.gene_type} {input.bed} |
        awk -F'\\t' -v OFS='\\t' 'FNR==NR{{ $7="gene" ; $15="exon" ; $16=$14 ; $14=$13 ; print ; if ( $8==1 ) {{ $15="first_exon" ; print }} else if ( $9==1 ) {{ $15="last_exon" }} }} FNR<NR{{print}}' {output.exon} - |
        awk -F'\\t' -v OFS='\\t' 'FNR==NR{{ $7="gene" ; $15="intron" ; $16=$14 ; $14=$13 ; print }} FNR<NR{{print}}' {output.intron} - |
        sort -k1,1 -k2,2n - > {output.genes}
        """

