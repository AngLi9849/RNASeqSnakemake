($7 == "gene" || $7 == "trscrpt") { print }

$6=="-" && ($7 != "gene" && $7 != "trscrpt") {
  if ($7 ~ "long_") {
    print
  }
  else if (gene!=$4) {
    gene=$4 ; a=$3 ; b=$2 ;  n=1 ; v=1 ; k=$12 ; $13=k ; $12=n ; $14=v ; $8=$8""n ; $9=$9" "n ; print ; 
    if (k>1 && ($7=="exon" || $7=="intron")) { 
      $7="first_"$7 ; print ;
    } 
  }
  else if (gene==$4) {
    a=$3 ;  
    if (a<b) {
      n+=1 ; b=$2 ; v=1 ; k=$12 ; $13=k ; $12=n ; $14=v ; $8=$8""n ; $9=$9" "n ; print
    }
    else if (a>=b) {
      n+=0 ;
      if (k>$12 && ($7=="exon" || $7=="intron") ) { 
        $7="long_"$7 ; $8=$4":"$7 ; $9=$8" "$7 ; $12=1 ; $13=1 ; $14=1 ; print
      }
      else if (k==$12) {
        v+=1 ;  $13=k ; $12=n ; $14=v ; $8=$8""n ; $9=$9" "n ; print
      }
    } ;
    if (k==1 && n>1 && ($7=="exon" || $7=="intron")) {
      $7="last_"$7 ; print
    } ;
  }
}

$6=="+" && ($7 != "gene" && $7 != "trscrpt") {
  if ($7 ~ "long_") {
    print
  }
  else if (gene!=$4) {
    gene=$4 ; c=$3 ; d=$2 ; w=$12 ; m=1 ; $13=m ; $8=$8""w ; $9=$9" "w ; print
    if (w>1 && ($7=="exon" || $7=="intron")) {
      $7="last_"$7 ; print ;
    }
  }
  else if (gene==$4) {
    c=$3 ;
   if (c<d) {
      m+=1 ;  w=$12 ; d=$2 ; $13=m ; $8=$8""w ; $9=$9" "w ; print
    }
    else if (c>=d) {
      m+=0 ; w=w ;
      if (w>$12 && ($7=="exon" || $7=="intron") ) {
        $7="long_"$7 ; $8=$4":"$7 ; $9=$8" "$7 ; $12=1 ; $13=1 ; $14=1 ; print
      }
      else if (w==$12) {
        w=$12 ; $13=m ; $8=$8""w ; $9=$9" "w ; print
      }
    } ; 
    if (w==1 && m>1 && ($7=="exon" || $7=="intron")) {
      $7="first_"$7 ; print
    } ;
  }
}

