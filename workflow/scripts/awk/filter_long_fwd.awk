FNR==NR { name[$1]=$2 }

FNR<NR && ($7 == "gene" || $7=="trscrpt") { print $0, 1, 1, 1} 

FNR < NR && ($7 != "gene" && $7 != "trscrpt") {
  if (gene!=$4) {
    gene=$4 ; a=$2 ; b=$3 ; n=1 ; v=1 ;  print $0, n, 1, v
  }
  else if (gene==$4) {
    a=$2 ;
    if (a>b) {
      n+=1 ; b=$3 ; v=1 ; print $0, n, 1, v
    }
    else if (a==b && ($7=="exon" || $7=="intron") ) {
      $7="long_"$7 ; $8=$4":"$7 ; $9=name[$4]" "$7 ; print $0, 1, 1, 1
    }
    else if (a<b && $3 >= b) {
      n+=0 ; v+=1 ; print $0, n, 1, v
    } 
    else if (a<b && $3 < b) {
      n+=0 ; v+=1 ; b=$3 ; print $0, n, 1, v
    } ;
  }
}

